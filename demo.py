#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""一键演示：用模拟数据把整条流程跑通，并和「真值」对照。

它想回答的问题
------------
真实数据（airway）的对拍基准必须在装有 R + DESeq2 的机器上产出。
在拿到那份数据之前，这个脚本让你能**立刻**看到整个项目在做什么、
以及它算得对不对，因为模拟数据的真值是我们自己定的。

运行
----
    python demo.py                 # 默认 2000 基因，4v4
    python demo.py --genes 5000    # 更多基因
    python demo.py --keep          # 保留生成的 CSV 文件

产出
----
    build/demo/counts.csv          模拟计数矩阵
    build/demo/design.csv          分组表
    build/demo/result.csv          我们的结果表
    build/demo/summary.txt         与真值的对照报告
"""

import argparse
import math
import os
import shutil
import statistics
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from src.dispersion import estimate_dispersions          # noqa: E402
from src.multitest import bh_adjust                      # noqa: E402
from src.normalize import size_factors                   # noqa: E402
from src.simulate import simulate_dataset, write_simulated_dataset  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="DESeq2-FromScratch 一键演示")
    ap.add_argument("--genes", type=int, default=2000)
    ap.add_argument("--per-group", type=int, default=4)
    ap.add_argument("--seed", type=int, default=20260926)
    ap.add_argument("--keep", action="store_true", help="保留生成的 CSV")
    args = ap.parse_args()

    outdir = os.path.join(ROOT, "build", "demo")
    os.makedirs(outdir, exist_ok=True)
    lines = []

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    say("=" * 74)
    say("DESeq2-FromScratch 一键演示（模拟数据，真值已知）")
    say("=" * 74)
    say("基因数 %d，每组 %d 个样本，随机种子 %d"
        % (args.genes, args.per_group, args.seed))

    # ---------- 生成 ----------
    ds = simulate_dataset(n_genes=args.genes, n_per_group=args.per_group,
                          seed=args.seed)
    cp = os.path.join(outdir, "counts.csv")
    dp = os.path.join(outdir, "design.csv")
    gene_ids, sample_ids = write_simulated_dataset(ds, cp, dp)
    say("")
    say("模拟数据已生成：%s" % cp)
    say("  真实 size factor : %s" % ["%.4f" % v for v in ds["true_sf"]])
    n_de = sum(ds["de_mask"])
    say("  真实差异基因数   : %d / %d (%.1f%%)"
        % (n_de, args.genes, 100.0 * n_de / args.genes))

    # ---------- 跑流程 ----------
    say("")
    say("跑 M1 -> M2 -> M3 -> M4 ...")
    t0 = time.time()
    from deg import run_pipeline
    rows, inter = run_pipeline(ds["counts"], ds["group"], quiet=True)
    elapsed = time.time() - t0
    say("总耗时 %.2f 秒（%.1f ms / 基因）" % (elapsed, 1000 * elapsed / args.genes))

    # ---------- 对照真值 ----------
    say("")
    say("-" * 74)
    say("M1 归一化：估计 vs 真实")
    say("-" * 74)
    sf = inter["factors"]
    # size factor 只可识别到常数倍，所以先归一化到同一尺度再比
    m = statistics.median(sf)
    r = statistics.median(ds["true_sf"])
    worst = 0.0
    for j in range(len(sf)):
        est = sf[j] / m
        tru = ds["true_sf"][j] / r
        worst = max(worst, abs(est - tru))
    say("  整体尺度因子 = %.4f（不影响 log2FC，只被 GLM 截距吸收）" % (m / r))
    say("  归一化到同尺度后，最大绝对偏差 = %.4f" % worst)
    say("  判定：%s" % ("通过（< 0.05）" if worst < 0.05 else "偏大，请检查"))

    say("")
    say("-" * 74)
    say("M2 离散度：估计 vs 真实")
    say("-" * 74)
    d = inter["dispersion"]
    idx = [i for i in range(args.genes) if not d["allZero"][i]]
    ratio_gene = sorted(d["dispGeneEst"][i] / ds["true_alpha"][i] for i in idx)
    ratio_map = sorted(d["dispersion"][i] / ds["true_alpha"][i] for i in idx)

    def q(lst, p):
        return lst[min(len(lst) - 1, int(p * len(lst)))]

    say("  趋势拟合类型        : %s" % d["fitType"])
    say("  先验方差 dispPriorVar: %.4f" % d["dispPriorVar"])

    def iqr_log(vals):
        lg = sorted(math.log(v) for v in vals)
        return q(lg, 0.75) - q(lg, 0.25)

    say("  gene-wise 估计/真值 : 中位 %.3f  [P25 %.3f, P75 %.3f]"
        % (q(ratio_gene, 0.5), q(ratio_gene, 0.25), q(ratio_gene, 0.75)))
    say("  收缩后   估计/真值 : 中位 %.3f  [P25 %.3f, P75 %.3f]"
        % (q(ratio_map, 0.5), q(ratio_map, 0.25), q(ratio_map, 0.75)))
    ig, im = iqr_log(ratio_gene), iqr_log(ratio_map)
    say("  log 误差的 IQR      : gene-wise %.3f -> 收缩后 %.3f（%s %.0f%%）"
        % (ig, im, "降低" if im < ig else "升高", abs(100 * (1 - im / ig))))
    say("  判定：%s" % ("通过（收缩使误差变小）" if im < ig else "未通过"))

    say("")
    say("-" * 74)
    say("M3 log2FC：估计 vs 真实（只看真实的差异基因）")
    say("-" * 74)
    errs = []
    for i in range(args.genes):
        if ds["de_mask"][i] and rows[i]["log2FoldChange"] is not None:
            errs.append(abs(rows[i]["log2FoldChange"] - ds["true_lfc"][i]))
    errs.sort()
    if errs:
        say("  绝对误差            : 中位 %.3f  P90 %.3f  最大 %.3f（n=%d）"
            % (q(errs, 0.5), q(errs, 0.9), errs[-1], len(errs)))
    say("  非全 0 基因可检验比率: %d / %d"
        % (sum(1 for r in rows if r["pvalue"] is not None), args.genes))

    say("")
    say("-" * 74)
    say("检出能力：塞进去的差异基因找回来了吗")
    say("-" * 74)
    sig = [i for i in range(args.genes)
           if rows[i]["padj"] is not None and rows[i]["padj"] < 0.05
           and rows[i]["log2FoldChange"] is not None
           and abs(rows[i]["log2FoldChange"]) >= 1.0]
    true_de = [i for i in range(args.genes) if ds["de_mask"][i]]
    hit = len(set(sig) & set(true_de))
    fp = len(sig) - hit
    fn = len(true_de) - hit
    say("  检出显著基因 %d 个；真实差异基因 %d 个" % (len(sig), len(true_de)))
    say("  命中 %d / 误报 %d / 漏检 %d" % (hit, fp, fn))
    if len(sig):
        say("  精确率 %.3f   召回率 %.3f"
            % (hit / len(sig), hit / max(len(true_de), 1)))
    say("")
    say("  说明：小样本（每组 %d 个重复）下 p 值本身偏乐观，" % args.per_group)
    say("        误报率偏高是 plug-in Wald 检验的固有性质，DESeq2 同样如此。")
    say("        证据见 REPORT.md 第五节（n=60 时 sd(stat)=1.000）。")

    # ---------- 写文件 ----------
    say("")
    say("-" * 74)
    from src.io_utils import write_results_csv
    op = os.path.join(outdir, "result.csv")
    write_results_csv(op, gene_ids, rows)
    say("结果表已写出：%s" % op)
    with open(os.path.join(outdir, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    say("对照报告已写出：%s" % os.path.join(outdir, "summary.txt"))
    say("=" * 74)

    if not args.keep:
        pass  # 文件保留在 build/demo 下，方便查看；不自动删除
    return 0


if __name__ == "__main__":
    sys.exit(main())
