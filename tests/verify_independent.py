# -*- coding: utf-8 -*-
"""独立复算：不经过 src/ 与 tests/compare_baseline.py 的任何函数，
只用标准库（csv / math / statistics）另算一遍关键结论。

动机
------------
`tests/compare_baseline.py` 既是「产出报告数字的工具」，又是「验证报告数字的
工具」：这两件事由同一段代码做，是自证。本脚本换一条完全独立的代码路径
（`statistics.correlation`、独立写的 size factor 公式）重算一遍。

要验的 5 件事：

  A. M1 size factor：按定义独立重算，与官方 sizeFactors.csv 比
  B. log2FC 相关系数：用 statistics.correlation 重算
  C. |Δlog2FC| 中位数、Jaccard、**两个反向控制**
  D. my_result.csv 内部自洽：列间关系、padj 单调性、与 counts 的行对齐
  E. 与中间量对账：log2FC == beta_ln/ln2、dispGeneEst 上界

WARNING: 这里的期望值是**写死的**：它们就是 REPORT.md 第八节发布的数字。
数字一变就说明「报告与实际不符」，必须同步改报告或查原因。换数据集时需要
同步更新这些期望值。

用法：python tests/verify_independent.py      （不依赖官方基准之外的东西）
"""

import csv
import math
import os
import random
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

FAILED = []


def check(name, ok, detail):
    print("  [%s] %-46s %s" % ("PASS" if ok else "FAIL", name, detail))
    if not ok:
        FAILED.append(name)


def read_rows(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    return rows[0], rows[1:]


def to_f(s):
    s = (s or "").strip()
    if s == "" or s.upper() in ("NA", "NAN", "NULL"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def main():
    print("=" * 76)
    print("独立复算（只用标准库 csv / math / statistics）")
    print("=" * 76)

    need = ["data/counts.csv", "data/my_result.csv",
            "data/intermediates/dispersion_mine.csv",
            "data/intermediates/glm_mine.csv",
            "data/reference/sizeFactors.csv",
            "data/reference/deseq2_result_noFilter.csv"]
    missing = [p for p in need if not os.path.exists(os.path.join(ROOT, p))]
    if missing:
        print("缺少文件：%s" % missing)
        print("（先跑 python deg.py ... 与笔记本的基准脚本）")
        return 2

    # ------------------------------------------------------------ A. M1
    hdr, rows = read_rows(os.path.join(ROOT, "data", "counts.csv"))
    samples = hdr[1:]
    genes = [r[0] for r in rows]
    counts = [[int(v) for v in r[1:]] for r in rows]
    n_gene, n_smp = len(counts), len(samples)
    print("\nA. M1 size factor 独立重算（%d 基因 x %d 样本）" % (n_gene, n_smp))

    loggeo = [None] * n_gene
    for i, row in enumerate(counts):
        if all(v > 0 for v in row):
            loggeo[i] = sum(math.log(v) for v in row) / n_smp
    mine_sf = []
    for j in range(n_smp):
        ratios = [math.log(counts[i][j]) - loggeo[i]
                  for i in range(n_gene) if loggeo[i] is not None]
        mine_sf.append(math.exp(statistics.median(ratios)))

    with open(os.path.join(ROOT, "data", "reference", "sizeFactors.csv"),
              "r", encoding="utf-8-sig", newline="") as fh:
        ref_sf = [float(x["sizeFactor"]) for x in csv.DictReader(fh)]
    check("A1 官方 sizeFactors 行数 = 样本数", len(ref_sf) == n_smp,
          "%d vs %d" % (len(ref_sf), n_smp))
    m = sum(mine_sf) / n_smp
    rr = sum(ref_sf) / n_smp
    rel = [abs(mine_sf[i] / m - ref_sf[i] / rr) / (ref_sf[i] / rr)
           for i in range(n_smp)]
    check("A2 独立重算 vs 官方：最大相对偏差 < 1e-12", max(rel) < 1e-12,
          "max = %.3e" % max(rel))

    from src.normalize import size_factors as src_sf
    impl = src_sf(counts)
    mi = sum(impl) / n_smp
    rel2 = [abs(impl[i] / mi - mine_sf[i] / m) / (mine_sf[i] / m)
            for i in range(n_smp)]
    check("A3 独立重算 vs src/normalize.py：< 1e-12", max(rel2) < 1e-12,
          "max = %.3e" % max(rel2))

    # ------------------------------------------------------ B/C. 结果表
    print("\nB/C. 结果表独立重算")
    mine, ref = {}, {}
    for path, box in ((os.path.join(ROOT, "data", "my_result.csv"), mine),
                      (os.path.join(ROOT, "data", "reference",
                                    "deseq2_result_noFilter.csv"), ref)):
        with open(path, "r", encoding="utf-8-sig", newline="") as fh:
            for row in csv.DictReader(fh):
                box[row["gene_id"]] = row

    xs, ys = [], []
    for g, rv in ref.items():
        a = to_f(mine.get(g, {}).get("log2FoldChange"))
        b = to_f(rv.get("log2FoldChange"))
        if a is not None and b is not None:
            xs.append(a)
            ys.append(b)
    r_stat = statistics.correlation(xs, ys)
    check("B1 log2FC 可比基因数 = 33469", len(xs) == 33469, "%d" % len(xs))
    check("B2 statistics.correlation = 0.9997 ± 0.0002",
          abs(r_stat - 0.9997) < 2e-4, "r = %.6f" % r_stat)
    same = sum(1 for a, b in zip(xs, ys) if (a > 0) == (b > 0))
    check("B3 符号一致率 = 99.72% ± 0.02%",
          abs(100.0 * same / len(xs) - 99.72) < 0.02,
          "%.2f%%" % (100.0 * same / len(xs)))
    d = sorted(abs(a - b) for a, b in zip(xs, ys))
    check("C1 |Δlog2FC| 中位 = 0.00055 ± 0.0002",
          abs(statistics.median(d) - 0.00055) < 2e-4,
          "median = %.6f" % statistics.median(d))
    check("C2 |Δlog2FC| > 0.5 的基因数 = 0",
          sum(1 for v in d if v > 0.5) == 0,
          "count = %d" % sum(1 for v in d if v > 0.5))

    def sigset(box):
        out = set()
        for g, v in box.items():
            p, l = to_f(v.get("padj")), to_f(v.get("log2FoldChange"))
            if p is not None and p < 0.05 and l is not None and abs(l) >= 1.0:
                out.add(g)
        return out

    sa, sb = sigset(mine), sigset(ref)
    j = len(sa & sb) / len(sa | sb)
    check("C3 显著集 Jaccard = 0.9017 ± 0.001", abs(j - 0.9017) < 1e-3,
          "J = %.6f（我 %d / 官方 %d）" % (j, len(sa), len(sb)))

    # 反向控制：打乱配对后 r 必须塌掉。若这个相关系数是用「两份各自排好序的
    # 列表」算出来的（经典错误写法），打乱后依然会很高，这里必须塌。
    rng = random.Random(20260926)
    shuf = ys[:]
    rng.shuffle(shuf)
    check("C4 反向控制：打乱配对后 |r| < 0.05",
          abs(statistics.correlation(xs, shuf)) < 0.05,
          "r_shuffled = %.6f" % statistics.correlation(xs, shuf))
    r_neg = statistics.correlation(xs, [-v for v in ys])
    check("C5 反向控制：方向翻转后 r = -0.9997",
          abs(r_neg + r_stat) < 1e-9, "r = %.6f" % r_neg)

    # ---------------------------------------------------- D. 内部自洽
    print("\nD. my_result.csv 内部自洽")
    check("D1 行数 = counts 行数", len(mine) == len(counts),
          "%d vs %d" % (len(mine), len(counts)))
    check("D2 基因顺序与 counts 完全一致", list(mine.keys()) == genes,
          "前 3 个：%s" % list(mine.keys())[:3])
    allzero = [g for i, g in enumerate(genes) if sum(counts[i]) == 0]
    az_ok = all(to_f(mine[g]["log2FoldChange"]) is None
                and to_f(mine[g]["padj"]) is None for g in allzero)
    check("D3 全 0 基因（%d 个）各项结果均为空" % len(allzero), az_ok,
          "全部为空" if az_ok else "有非空值")

    ln2 = math.log(2.0)
    bad_stat = bad_p = 0
    pairs = []
    for g in genes:
        row = mine[g]
        lfc, se = to_f(row["log2FoldChange"]), to_f(row["lfcSE"])
        st, pv, pj = to_f(row["stat"]), to_f(row["pvalue"]), to_f(row["padj"])
        if lfc is None:
            continue
        if abs(st - lfc / se) > 1e-9 * max(1.0, abs(st)):
            bad_stat += 1
        if abs(pv - math.erfc(abs(st) / math.sqrt(2.0))) > 1e-12 * max(1e-12, pv):
            bad_p += 1
        if pj is not None:
            pairs.append((pv, pj))
    check("D4 stat == log2FoldChange / lfcSE", bad_stat == 0,
          "违规 %d 个" % bad_stat)
    check("D5 pvalue == erfc(|stat|/sqrt2)", bad_p == 0, "违规 %d 个" % bad_p)
    # WARNING: 踩过坑：BH 的单调性是「按 p 值升序后 padj 不降」。第一版按 counts 的
    # 基因顺序逐个比，得到 2 万多个"非单调"：那是验证脚本自己写错了顺序。
    pairs.sort()
    viol_mono = sum(1 for k in range(1, len(pairs))
                    if pairs[k][1] < pairs[k - 1][1] - 1e-15)
    viol_range = sum(1 for pv, pj in pairs
                     if pj < pv - 1e-12 or pj > 1.0 + 1e-12)
    check("D6 padj ∈ [pvalue, 1] 且按 p 升序单调不降",
          viol_range == 0 and viol_mono == 0,
          "越界 %d / 非单调 %d（%d 个非空）"
          % (viol_range, viol_mono, len(pairs)))

    # -------------------------------------------------  E. 中间量对账
    print("\nE. 结果表与中间量对账")
    bad = 0
    with open(os.path.join(ROOT, "data", "intermediates", "glm_mine.csv"),
              "r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            b = to_f(row["beta_ln"])
            lfc = to_f(mine[row["gene_id"]]["log2FoldChange"])
            if b is None or lfc is None:
                continue
            if abs(lfc - b / ln2) > 1e-9 * max(1.0, abs(lfc)):
                bad += 1
    check("E1 log2FoldChange == beta_ln / ln2", bad == 0, "违规 %d 个" % bad)

    mx = 0.0
    with open(os.path.join(ROOT, "data", "intermediates",
                           "dispersion_mine.csv"),
              "r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            v = to_f(row["dispGeneEst"])
            if v is not None:
                mx = max(mx, v)
    check("E2 dispGeneEst 最大值 = 10.0（上界生效）", abs(mx - 10.0) < 1e-9,
          "max = %.6f" % mx)

    print("\n" + "=" * 76)
    if FAILED:
        print("结论：%d 项未通过 -> %s" % (len(FAILED), FAILED))
        return 1
    print("结论：全部通过（独立复算 × 对拍脚本 × 官方基准，三方一致）")
    print("=" * 76)
    return 0


if __name__ == "__main__":
    sys.exit(main())
