#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""DESeq2-FromScratch 命令行入口。

用法
----
    python deg.py --counts data/counts.csv --design data/design.csv \
                  --out data/my_result.csv

输出列与 DESeq2 的 results() 完全对齐：
    gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

另外可以加 --dump-intermediates <目录> 把每一步的中间量导出，用于和官方
DESeq2 的 sizeFactors / dispGeneEst / dispFit / dispersions 逐段对拍。

重要约定
--------
* log2FoldChange 的方向：design.csv 里 group = 1 的组 相对于 group = 0 的组。
  即 log2FC > 0 表示在 group = 1 中表达更高。
* 全 0 基因、以及无法拟合的基因，输出空值（与 DESeq2 输出 NA 一致），
  而不是被悄悄丢掉，保证两条路径的基因集合完全一致，才能逐行对拍。
* 本实现不做 DESeq2 的 independent filtering 与 Cook's 距离离群替换，
  所以对拍时应当使用官方的 deseq2_result_noFilter.csv（见 REPORT.md）。
"""

import argparse
import math
import os
import platform
import sys
import time

# 允许直接以脚本方式运行（python deg.py），也能作为模块被导入
if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from src import __version__
    from src.dispersion import estimate_dispersions
    from src.glm_nb import fit_nb_glm, wald_test
    from src.io_utils import (align_design, read_counts_csv, read_design_csv,
                              write_generic_csv, write_results_csv)
    from src.linalg import build_design_matrix
    from src.multitest import bh_adjust
    from src.normalize import base_mean, normalized_counts, size_factors
else:
    from . import __version__
    from .dispersion import estimate_dispersions
    from .glm_nb import fit_nb_glm, wald_test
    from .io_utils import (align_design, read_counts_csv, read_design_csv,
                           write_generic_csv, write_results_csv)
    from .linalg import build_design_matrix
    from .multitest import bh_adjust
    from .normalize import base_mean, normalized_counts, size_factors


def run_pipeline(counts, groups, limit=None, progress=None, quiet=False):
    """跑完 M1 -> M2 -> M3 -> M4，返回结果行列表与中间量。

    参数
    ----
    counts : counts[基因][样本]
    groups : 各样本的 0/1 分组
    limit  : 只跑前 N 个基因（用于快速冒烟测试）
    """
    n_genes_all = len(counts)
    if limit is not None and limit > 0:
        counts = counts[:limit]
        gene_slice = list(range(limit))
    else:
        gene_slice = list(range(n_genes_all))
    n_genes = len(counts)

    # ---------- M1：归一化 ----------
    factors = size_factors(counts)

    # ---------- M2：离散度 ----------
    disp = estimate_dispersions(counts, factors, groups,
                                quiet=quiet, progress=progress)

    # ---------- M3：GLM + Wald ----------
    X = build_design_matrix(groups)
    log_sf = [math.log(f) for f in factors]
    norm = normalized_counts(counts, factors)
    base_means = base_mean(norm)

    rows = []
    pvalues = []
    lfc_list = []
    se_list = []
    stat_list = []
    beta_list = []
    beta_se_list = []
    iterations = []
    for i in range(n_genes):
        if disp["allZero"][i]:
            rows.append({"baseMean": 0.0, "log2FoldChange": None,
                         "lfcSE": None, "stat": None, "pvalue": None,
                         "padj": None})
            pvalues.append(None)
            lfc_list.append(None)
            se_list.append(None)
            stat_list.append(None)
            beta_list.append(None)
            beta_se_list.append(None)
            iterations.append(0)
            continue

        alpha = disp["dispersion"][i]
        if alpha != alpha or alpha <= 0.0:
            rows.append({"baseMean": base_means[i], "log2FoldChange": None,
                         "lfcSE": None, "stat": None, "pvalue": None,
                         "padj": None})
            pvalues.append(None)
            lfc_list.append(None)
            se_list.append(None)
            stat_list.append(None)
            beta_list.append(None)
            beta_se_list.append(None)
            iterations.append(0)
            continue

        fit = fit_nb_glm(counts[i], X, log_sf, alpha)
        lfc, se, stat, p = wald_test(fit["beta"], fit["cov"], coef_index=1)
        rows.append({"baseMean": base_means[i], "log2FoldChange": lfc,
                     "lfcSE": se, "stat": stat, "pvalue": p, "padj": None})
        pvalues.append(p)
        lfc_list.append(lfc)
        se_list.append(se)
        stat_list.append(stat)
        beta_list.append(fit["beta"][1])
        beta_se_list.append(
            math.sqrt(fit["cov"][1][1]) if fit["cov"][1][1] == fit["cov"][1][1]
            and fit["cov"][1][1] > 0 else None)
        iterations.append(fit["iterations"])

    # ---------- M4：BH 校正 ----------
    padj = bh_adjust(pvalues)
    for i in range(n_genes):
        rows[i]["padj"] = padj[i]

    intermediates = {
        "factors": factors,
        "dispersion": disp,
        "base_means": base_means,
        "pvalues": pvalues,
        "beta": beta_list,
        "beta_se": beta_se_list,
        "iterations": iterations,
        "gene_slice": gene_slice,
    }
    return rows, intermediates


def dump_intermediates(outdir, gene_ids, intermediates):
    """把 M1/M2/M3 的中间量写成 CSV，用于与官方 DESeq2 逐段对拍。"""
    os.makedirs(outdir, exist_ok=True)
    slice_ = intermediates["gene_slice"]
    factors = intermediates["factors"]
    disp = intermediates["dispersion"]

    # M1：size factor（按样本）
    write_generic_csv(os.path.join(outdir, "sizeFactors_mine.csv"),
                      ["sample_index", "sizeFactor"],
                      [[j + 1, factors[j]] for j in range(len(factors))])

    # M2：三种离散度
    header = ["gene_id", "baseMean", "baseVar", "allZero",
              "dispGeneEst", "dispFit", "dispMAP", "dispersion", "dispOutlier"]
    rows = []
    for k, gi in enumerate(slice_):
        rows.append([
            gene_ids[gi],
            disp["baseMean"][k],
            disp["baseVar"][k],
            "TRUE" if disp["allZero"][k] else "FALSE",
            disp["dispGeneEst"][k],
            disp["dispFit"][k],
            disp["dispMAP"][k],
            disp["dispersion"][k],
            "TRUE" if disp["dispOutlier"][k] else "FALSE",
        ])
    write_generic_csv(os.path.join(outdir, "dispersion_mine.csv"), header, rows)

    # M3：GLM 的原始系数与迭代次数
    header = ["gene_id", "beta_ln", "betaSE_ln", "iterations", "pvalue"]
    rows = []
    for k, gi in enumerate(slice_):
        rows.append([gene_ids[gi], intermediates["beta"][k],
                     intermediates["beta_se"][k], intermediates["iterations"][k],
                     intermediates["pvalues"][k]])
    write_generic_csv(os.path.join(outdir, "glm_mine.csv"), header, rows)

    # M2 的先验方差（标量）
    write_generic_csv(os.path.join(outdir, "dispersion_prior_mine.csv"),
                      ["key", "value"],
                      [["fitType", disp["fitType"]],
                       ["dispPriorVar", disp["dispPriorVar"]],
                       ["dispPriorVarMethod",
                        disp.get("dispPriorVarMethod", "trigamma")],
                       ["varLogDispEsts", disp["varLogDispEsts"]]])


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="用纯 Python 标准库从零实现 RNA-seq 差异表达分析"
                    "（不调用 DESeq2）")
    parser.add_argument("--counts", required=True, help="计数矩阵 CSV")
    parser.add_argument("--design", required=True, help="分组表 CSV")
    parser.add_argument("--out", required=True, help="结果表输出路径")
    parser.add_argument("--dump-intermediates", default=None,
                        help="把 M1/M2/M3 的中间量导出到这个目录")
    parser.add_argument("--limit", type=int, default=None,
                        help="只分析前 N 个基因（快速冒烟测试用）")
    parser.add_argument("--quiet", action="store_true", help="不打印进度")
    parser.add_argument("--version", action="version",
                        version="DESeq2-FromScratch %s" % __version__)
    args = parser.parse_args(argv)

    t_start = time.time()
    log = (lambda *a: None) if args.quiet else (lambda *a: print(*a, flush=True))

    log("=" * 68)
    log("DESeq2-FromScratch v%s" % __version__)
    log("=" * 68)
    log("Python      : %s" % sys.version.split()[0])
    log("平台        : %s %s" % (platform.system(), platform.release()))
    log("计数矩阵    : %s" % args.counts)
    log("分组表      : %s" % args.design)
    log("开始时间    : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))

    gene_ids, sample_ids, counts = read_counts_csv(args.counts)
    design = read_design_csv(args.design)
    groups, conditions, cells = align_design(sample_ids, design)

    log("-" * 68)
    log("基因数      : %d" % len(gene_ids))
    log("样本数      : %d" % len(sample_ids))
    log("样本顺序    : %s" % ", ".join(sample_ids))
    log("分组(1/0)   : %s" % ", ".join(str(g) for g in groups))
    log("方向约定    : log2FC > 0 表示在 group=1 中更高")
    if len(set(groups)) != 2:
        log("!! 分组不是两组，退出")
        return 2

    def progress(done, total):
        log("  [M2] 离散度估计 %d / %d (%.0f%%)"
            % (done, total, 100.0 * done / max(total, 1)))

    rows, intermediates = run_pipeline(counts, groups, limit=args.limit,
                                       progress=progress, quiet=args.quiet)

    out_ids = gene_ids[:args.limit] if args.limit else gene_ids
    write_results_csv(args.out, out_ids, rows)
    log("-" * 68)
    log("M1 size factor : %s" % ", ".join("%.4f" % f
                                          for f in intermediates["factors"]))
    d = intermediates["dispersion"]
    log("M2 趋势拟合    : %s" % d["fitType"])
    log("M2 先验方差    : %.4f" % d["dispPriorVar"])
    n_sig = sum(1 for r in rows if r["pvalue"] is not None)
    n_padj = sum(1 for r in rows if r["padj"] is not None)
    log("M3 可检验基因  : %d / %d（其余为全 0 或无法拟合）" % (n_sig, len(rows)))
    log("M4 校正后有效  : %d" % n_padj)
    log("显著基因数     : padj<0.05 且 |log2FC|>=1 的有 %d 个"
        % sum(1 for r in rows
              if r["padj"] is not None and r["padj"] < 0.05
              and r["log2FoldChange"] is not None
              and abs(r["log2FoldChange"]) >= 1.0))
    log("结果已写出     : %s" % args.out)

    if args.dump_intermediates:
        dump_intermediates(args.dump_intermediates, gene_ids, intermediates)
        log("中间量已导出   : %s" % args.dump_intermediates)

    log("总耗时         : %.1f 秒" % (time.time() - t_start))
    log("=" * 68)
    return 0


if __name__ == "__main__":
    sys.exit(main())
