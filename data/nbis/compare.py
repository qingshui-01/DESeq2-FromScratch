"""对拍第二数据集（NBISweden/workshop-RNAseq，DSS day00 vs day07，3v3）。

官方结果表 data/nbis/dge_results.csv 的性质（已实测核实）：
  * 完整 DESeq2 结果列：baseMean / log2FoldChange / lfcSE / stat / pvalue / padj
  * padj 是「独立过滤后」的 BH：baseMean <= 21.3 的基因被过滤，padj 记为 1；
    保留的 5642 个基因的 padj 与本项目 bh_adjust(their_pvalues) 逐位一致
  * 因此对拍口径：
      - log2FC / lfcSE / stat / pvalue：直接逐基因比
      - padj：在官方保留集上比（我方结果先按同一过滤阈值过滤再做 BH）
      - 显著集合：padj<0.05 的 Jaccard
"""
import csv
import math
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from src.multitest import bh_adjust

HERE = os.path.dirname(os.path.abspath(__file__))


def load_results(path):
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        r = csv.reader(f)
        header = [h.strip().strip('"') for h in next(r)]
        id_col = 0  # 第一列是基因 ID（我们的叫 gene_id，官方的是空列名）
        for row in r:
            if not row:
                continue
            d = dict(zip(header[1:], [v.strip().strip('"') for v in row[1:]]))
            out[row[0].strip().strip('"')] = d
    return out


def fnum(v):
    if v is None:
        return None
    try:
        x = float(v)
    except ValueError:
        return None
    return None if math.isnan(x) else x


def corr(xs, ys):
    return statistics.correlation(xs, ys)


def main():
    mine = load_results(os.path.join(HERE, "my_result.csv"))
    ref = load_results(os.path.join(HERE, "dge_results.csv"))
    common = sorted(set(mine) & set(ref))
    print("基因交集: %d（我 %d / 官方 %d）" % (len(common), len(mine), len(ref)))

    cols = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue"]
    print("\n=== 逐列对拍（全部可比基因）===")
    for c in cols:
        xs, ys = [], []
        for g in common:
            a, b = fnum(mine[g].get(c)), fnum(ref[g].get(c))
            if a is not None and b is not None:
                xs.append(a)
                ys.append(b)
        r = corr(xs, ys)
        diffs = sorted(abs(a - b) for a, b in zip(xs, ys))
        med = diffs[len(diffs) // 2]
        p99 = diffs[int(len(diffs) * 0.99)]
        print("  %-14s n=%5d  r=%.6f  |diff|中位=%.3g  |diff|P99=%.3g"
              % (c, len(xs), r, med, p99))

    # log2FC 符号一致率
    xs, ys = [], []
    for g in common:
        a, b = fnum(mine[g].get("log2FoldChange")), fnum(ref[g].get("log2FoldChange"))
        if a is not None and b is not None:
            xs.append(a)
            ys.append(b)
    same = sum(1 for a, b in zip(xs, ys) if (a > 0) == (b > 0))
    print("\nlog2FC 符号一致率: %.2f%%" % (100.0 * same / len(xs)))

    # p 值 -log10 尺度
    lp_m, lp_r = [], []
    for g in common:
        a, b = fnum(mine[g].get("pvalue")), fnum(ref[g].get("pvalue"))
        if a is not None and b is not None and a > 0 and b > 0:
            lp_m.append(-math.log10(a))
            lp_r.append(-math.log10(b))
    print("p 值相关系数（-log10 尺度）: r = %.6f (n=%d)"
          % (corr(lp_m, lp_r), len(lp_m)))

    # padj：在官方保留集（padj<1 的基因，等价于独立过滤后的集合）上比
    kept = [g for g in common
            if (fnum(ref[g].get("padj")) or 1.0) < 1.0 - 1e-12
            and fnum(ref[g].get("pvalue")) is not None]
    ref_pv = {g: fnum(ref[g]["pvalue"]) for g in kept}
    # 我方按同一批基因做 BH
    my_p = [fnum(mine[g].get("pvalue")) for g in kept]
    my_adj = bh_adjust(my_p)
    ref_adj = [fnum(ref[g]["padj"]) for g in kept]
    lm = [-math.log10(max(a, 1e-300)) for a in my_adj]
    lr = [-math.log10(max(a, 1e-300)) for a in ref_adj]
    print("\n=== padj（官方过滤集 %d 个基因，-log10 尺度）===" % len(kept))
    print("  padj r = %.6f" % corr(lm, lr))

    # 显著集合 Jaccard（各自口径：官方 padj<0.05；我方在过滤集上的 BH<0.05）
    sig_ref = {g for g in kept if fnum(ref[g]["padj"]) < 0.05}
    sig_me = {g for g, a in zip(kept, my_adj) if a < 0.05}
    j = len(sig_ref & sig_me) / len(sig_ref | sig_me)
    print("  显著基因: 官方 %d / 我方 %d，Jaccard = %.4f"
          % (len(sig_ref), len(sig_me), j))

    # 无条件 padj 口径（我方全量 BH vs 官方过滤 BH）显著集差异说明
    my_all_adj = bh_adjust([fnum(mine[g].get("pvalue")) for g in common])
    sig_me_all = {g for g, a in zip(common, my_all_adj)
                  if a is not None and a < 0.05}
    print("  （参考：我方不做独立过滤时显著基因 %d 个）" % len(sig_me_all))


if __name__ == "__main__":
    main()
