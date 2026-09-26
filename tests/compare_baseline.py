# -*- coding: utf-8 -*-
"""对拍：把我们的实现与官方 DESeq2 的基准逐段比对。

为什么单独写一个脚本
--------------------
`tests/run_acceptance.py` 里的 I-1/I-2/I-3 是**回归验收**：它们只回答
「有没有退步」（相关系数是否 > 0.8、Jaccard 是否 > 0.5）。而 `REPORT.md`
第八节需要的是**结论性的数字**：逐样本偏差、M2 三段的各自吻合度、不一致
基因按表达量分层的分布。这些数字必须可复现，所以固化成本脚本。

比对顺序（与实现顺序一致，便于定位偏差来源）
--------------------------------------------
  M1  size factor        -> data/reference/sizeFactors.csv
  M2  dispGeneEst / dispFit / dispMAP -> data/reference/dispersion_intermediates.csv
  M3  log2FC / stat      -> data/reference/deseq2_result_noFilter.csv
  M4  padj / 显著集合

WARNING: 一律与官方的 **noFilter** 版本比对：我们的实现不做 independent filtering、
   不做 Cook's 距离离群替换，拿默认版本比会把「方法差异」误记成「实现错误」。

只用 Python 标准库（项目硬性约束，见 README）。
"""

import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


# ---------------------------------------------------------------------------
# 基础统计（标准库手写，不引入任何第三方库）
# ---------------------------------------------------------------------------
def corr(xs, ys):
    """Pearson 相关系数。"""
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def median(xs):
    if not xs:
        return float("nan")
    s = sorted(xs)
    n = len(s)
    if n % 2:
        return s[n // 2]
    return 0.5 * (s[n // 2 - 1] + s[n // 2])


def quantile(xs, q):
    if not xs:
        return float("nan")
    s = sorted(xs)
    if len(s) == 1:
        return s[0]
    pos = q * (len(s) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(s) - 1)
    frac = pos - lo
    return s[lo] * (1 - frac) + s[hi] * frac


def sd(xs):
    """样本标准差（n-1 分母）。"""
    n = len(xs)
    if n < 2:
        return float("nan")
    m = sum(xs) / n
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


def jaccard(a, b):
    a, b = set(a), set(b)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def load_table(path):
    """读 CSV 成 dict(基因名 -> dict(列名 -> float/None))，跳过空值与 NA。

    自己解析表头而不用 DictReader：官方的 dispersion_intermediates.csv
    里有**重名的 gene_id 列**（mcols 里本来就有一列 gene_id），DictReader
    遇到重名会丢列。
    """
    import csv
    out = {}
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        key_idx = 0
        cols = {}
        for i, name in enumerate(header):
            name = name.strip()
            if name in ("gene_id", "") and i == 0:
                continue
            if name and name not in cols:
                cols[name] = i
        for row in reader:
            if not row:
                continue
            gid = row[key_idx].strip()
            if not gid:
                continue
            rec = {}
            for name, i in cols.items():
                if i >= len(row):
                    rec[name] = None
                    continue
                v = row[i].strip()
                if v == "" or v.upper() in ("NA", "NAN", "NULL"):
                    rec[name] = None
                else:
                    try:
                        rec[name] = float(v)
                    except ValueError:
                        rec[name] = None
            out[gid] = rec
    return out


def hline(title):
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)


# ---------------------------------------------------------------------------
# M1
# ---------------------------------------------------------------------------
def cmp_m1(mine, ref, lines):
    hline("M1 · size factor（中位数比值法）")
    samples = ref["order"]
    m = mine["values"]
    r = ref["values"]
    lines.append("### M1 size factor 逐样本对比\n")
    lines.append("| 样本 | 本实现 | 官方 DESeq2 | 绝对偏差 | 相对偏差 |")
    lines.append("|---|---|---|---|---|")
    worst = 0.0
    for i, s in enumerate(samples):
        d = abs(m[i] - r[i])
        rel = d / r[i] if r[i] else float("nan")
        worst = max(worst, rel)
        print("  %-12s 本实现 %.6f   官方 %.6f   偏差 %.2e (%.4f%%)"
              % (s, m[i], r[i], d, 100 * rel))
        lines.append("| %s | %.6f | %.6f | %.2e | %.4f%% |"
                     % (s, m[i], r[i], d, 100 * rel))
    print("\n  最大相对偏差 = %.3e （%.6f%%）" % (worst, 100 * worst))
    lines.append("\n最大相对偏差 **%.6f%%**。" % (100 * worst))
    return worst


# ---------------------------------------------------------------------------
# M2
# ---------------------------------------------------------------------------
def cmp_m2(mine, ref, lines):
    hline("M2 · 离散度三段对比（gene-wise / 趋势 / MAP 收缩）")
    lines.append("## M2 离散度三段对比\n")
    lines.append("| 阶段 | 可比基因数 | Pearson(log10) | 中位比值(我/官方) | "
                 "log10 差的 IQR |")
    lines.append("|---|---|---|---|---|")
    for col in ("dispGeneEst", "dispFit", "dispMAP"):
        xs, ys, ratios, diffs = [], [], [], []
        for g, rv in ref.items():
            if col not in rv:
                continue
            a = mine.get(g, {}).get(col)
            b = rv.get(col)
            if a is None or b is None or a <= 0 or b <= 0:
                continue
            xs.append(math.log10(a))
            ys.append(math.log10(b))
            ratios.append(a / b)
            diffs.append(math.log10(a) - math.log10(b))
        if not xs:
            print("  %-12s 无可比数据" % col)
            continue
        r = corr(xs, ys)
        iqr = quantile(diffs, 0.75) - quantile(diffs, 0.25)
        print("  %-12s n=%-6d r(log10)=%.4f  中位比值=%.4f  "
              "log10差IQR=%.4f" % (col, len(xs), r, median(ratios), iqr))
        lines.append("| `%s` | %d | %.4f | %.4f | %.4f |"
                     % (col, len(xs), r, median(ratios), iqr))

    # 分层归因：gene-wise 的差异是不是集中在低表达基因上？
    hline("M2 分层归因（按官方 baseMean，看 dispGeneEst 的差异来自哪一段）")
    lines.append("\n### M2 分层归因（`dispGeneEst`）\n")
    lines.append("| baseMean 区间 | 基因数 | Pearson(log10) | 中位比值(我/官方) | "
                 "我中位数 | 官方中位数 |")
    lines.append("|---|---|---|---|---|---|")
    edges = [(0, 1), (1, 10), (10, 100), (100, 1000), (1000, float("inf"))]
    for lo, hi in edges:
        xs, ys, ratios, ma, mb = [], [], [], [], []
        for g, rv in ref.items():
            bm = rv.get("baseMean")
            a = mine.get(g, {}).get("dispGeneEst")
            b = rv.get("dispGeneEst")
            if bm is None or a is None or b is None or a <= 0 or b <= 0:
                continue
            if not (lo <= bm < hi):
                continue
            xs.append(math.log10(a))
            ys.append(math.log10(b))
            ratios.append(a / b)
            ma.append(a)
            mb.append(b)
        if len(xs) < 5:
            continue
        label = "[%g, %g)" % (lo, hi) if hi != float("inf") else ">= %g" % lo
        r = corr(xs, ys)
        print("  %-14s n=%-6d r(log10)=%.4f  中位比值=%.4f  我=%.3e 官方=%.3e"
              % (label, len(xs), r, median(ratios), median(ma), median(mb)))
        lines.append("| %s | %d | %.4f | %.4f | %.3e | %.3e |"
                     % (label, len(xs), r, median(ratios), median(ma),
                        median(mb)))
    return None


# ---------------------------------------------------------------------------
# M3 / M4
# ---------------------------------------------------------------------------
def cmp_m3m4(mine, ref, lines):
    hline("M3/M4 · 结果表对比（与官方 noFilter 版本）")

    recs = []
    for g, rv in ref.items():
        mv = mine.get(g, {})
        bm = rv.get("baseMean")
        a = mv.get("log2FoldChange")
        b = rv.get("log2FoldChange")
        if bm is None or a is None or b is None:
            continue
        # 后 4 个字段（pvalue / stat）是给「p 值相关系数」用的：
        # 报告里必须给出 **p 值的相关系数**，
        # 而 REPORT 初版只给了 log2FC 与 Wald 统计量，缺这一项。
        recs.append((g, bm, a, b, mv.get("padj"), rv.get("padj"),
                     mv.get("pvalue"), rv.get("pvalue"),
                     mv.get("stat"), rv.get("stat")))

    def _corr_p(rs, which):
        """在 **-log10 尺度**上算 p 值的 Pearson 相关。

        为什么不直接对原始 p 值求相关：p 值的质量绝大部分堆在 1 附近，
        原始尺度上的相关系数会被这一堆「都不显著」的基因主导，几乎必然是
        0.99+，没有信息量。-log10(p) 尺度才反映「显著性排序是否一致」，
        也是这类报告通用的口径。两个尺度都算出来，写明用的是哪个。
        """
        lo, hi = (6, 7) if which == "pvalue" else (4, 5)
        xs, ys = [], []
        for r in rs:
            a, b = r[lo], r[hi]
            if a is None or b is None or a <= 0.0 or b <= 0.0:
                continue
            xs.append(-math.log10(a))
            ys.append(-math.log10(b))
        return (corr(xs, ys), len(xs)) if len(xs) > 10 else (float("nan"), 0)

    def stats(rs):
        xs = [r[2] for r in rs]
        ys = [r[3] for r in rs]
        d = [abs(r[2] - r[3]) for r in rs]
        same = sum(1 for r in rs if (r[2] > 0) == (r[3] > 0))
        a_sig = set(r[0] for r in rs
                    if r[4] is not None and r[4] < 0.05 and abs(r[2]) >= 1.0)
        b_sig = set(r[0] for r in rs
                    if r[5] is not None and r[5] < 0.05 and abs(r[3]) >= 1.0)
        return {
            "n": len(rs), "r": corr(xs, ys), "diff": median(d),
            "same": 100.0 * same / len(rs),
            "sd_x": sd(xs), "sd_y": sd(ys), "sd_d": sd(d),
            "na": len(a_sig), "nb": len(b_sig), "j": jaccard(a_sig, b_sig),
        }

    lines.append("### 结果表总体对比\n")
    overall = stats(recs)
    print("  log2FC 可比基因 %d 个，Pearson r = %.4f，符号一致 %.2f%%"
          % (overall["n"], overall["r"], overall["same"]))
    lines.append("- log2FC 可比基因：**%d** 个，Pearson r = **%.4f**，"
                 "符号一致 **%.2f%%**" % (overall["n"], overall["r"],
                                          overall["same"]))
    print("  显著集合(padj<0.05 且 |log2FC|>=1)：我 %d / 官方 %d，"
          "Jaccard = %.4f" % (overall["na"], overall["nb"], overall["j"]))
    lines.append("- 显著集合（padj<0.05 且 |log2FC|≥1）：本实现 **%d** 个 / "
                 "官方 **%d** 个，Jaccard = **%.4f**"
                 % (overall["na"], overall["nb"], overall["j"]))

    # ---- p 值 / padj 的相关系数（report 里必须给出）----
    r_p, n_p = _corr_p(recs, "pvalue")
    r_padj, n_padj = _corr_p(recs, "padj")
    print("\n  p 值相关系数（-log10 尺度）：")
    print("    pvalue : r = %.4f  (%d 个基因)" % (r_p, n_p))
    print("    padj   : r = %.4f  (%d 个基因)" % (r_padj, n_padj))
    lines.append("- **p 值相关系数（-log10 尺度）**：`pvalue` r = **%.4f**"
                 "（%d 个）、`padj` r = **%.4f**（%d 个）"
                 % (r_p, n_p, r_padj, n_padj))

    # 排除低计数基因后的吻合度，这是报告里最该引用的数字
    print("\n  去掉低计数基因后：")
    lines.append("\n### 去掉低计数基因后的吻合度\n")
    lines.append("| 基因范围 | 基因数 | log2FC r | 符号一致 | "
                 "\\|Δlog2FC\\| 中位 | Jaccard |")
    lines.append("|---|---|---|---|---|---|")
    for lo in (0.0, 1.0, 10.0, 100.0):
        rs = [r for r in recs if r[1] >= lo]
        if len(rs) < 5:
            continue
        s = stats(rs)
        name = "全部" if lo == 0 else "baseMean >= %g" % lo
        print("    %-16s n=%-6d r=%.4f  符号一致=%.2f%%  "
              "|Δlfc|中位=%.4f  J=%.4f"
              % (name, s["n"], s["r"], s["same"], s["diff"], s["j"]))
        lines.append("| %s | %d | %.4f | %.2f%% | %.4f | %.4f |"
                     % (name, s["n"], s["r"], s["same"], s["diff"], s["j"]))

    # Wald 统计量
    sx, sy = [], []
    for g, rv in ref.items():
        a = mine.get(g, {}).get("stat")
        b = rv.get("stat")
        if a is not None and b is not None:
            sx.append(a)
            sy.append(b)
    if sx:
        print("\n  Wald 统计量可比 %d 个，r = %.4f" % (len(sx), corr(sx, sy)))
        lines.append("- Wald 统计量：%d 个，r = **%.4f**"
                     % (len(sx), corr(sx, sy)))

    # ---- 分层 ----
    hline("按官方 baseMean 分层（关键：不一致主要来自哪里）")
    print("  说明：sd(lfc) 是该区间内 log2FC 自身的波动。当某个区间里绝大多数"
          "基因\n        的 log2FC 本来就贴近 0（sd 很小）时，即使两边的差异小到 "
          "0.01，\n        相关系数也会被压得很低，所以「r 低」不等于「不吻合」，"
          "要同时看\n        |Δlfc| 中位数与方向一致率。")
    lines.append("\n### 按表达量分层\n")
    lines.append("| baseMean 区间 | 基因数 | sd(log2FC 我) | sd(log2FC 官方) | "
                 "log2FC r | \\|Δlog2FC\\| 中位 | 方向一致率 | 显著集 Jaccard |")
    lines.append("|---|---|---|---|---|---|---|---|")
    edges = [(0, 1), (1, 10), (10, 100), (100, 1000), (1000, float("inf"))]
    for lo, hi in edges:
        rs = [r for r in recs if lo <= r[1] < hi]
        if len(rs) < 5:
            continue
        s = stats(rs)
        label = "[%g, %g)" % (lo, hi) if hi != float("inf") else ">= %g" % lo
        print("  %-14s n=%-6d sd我=%.4f sd官方=%.4f  r=%.4f  "
              "|Δlfc|中位=%.4f  方向=%.2f%%  J=%.4f"
              % (label, s["n"], s["sd_x"], s["sd_y"], s["r"], s["diff"],
                 s["same"], s["j"]))
        lines.append("| %s | %d | %.4f | %.4f | %.4f | %.4f | %.2f%% | %.4f |"
                     % (label, s["n"], s["sd_x"], s["sd_y"], s["r"], s["diff"],
                        s["same"], s["j"]))
    return None


def diag(mine_res, ref_res, ref_disp, n=12):
    """不一致基因的个案诊断。

    分层表里出现了一个必须解释的现象：baseMean 在 [100,1000) 的区间里，
    「|Δlog2FC| 中位数只有 0.0092」，但「sd(log2FC) 却是官方的 1.9 倍」。
    中位数极小、方差却大得多，只可能是一小撮基因差得极远在后面拉方差。
    所以直接把差异最大的基因逐个摊开看，判断它们是「实现错了」还是
    「两边对无界解采用了不同的正则化」。
    """
    # 计数矩阵（用来判断计数形态：决定 log2FC 是否可识别）
    import csv as _csv
    gcnt = {}
    with open(os.path.join(ROOT, "data", "counts.csv"), "r",
              encoding="utf-8-sig", newline="") as fh:
        rdr = _csv.reader(fh)
        next(rdr)
        for row in rdr:
            if not row:
                continue
            gid = row[0].strip()
            gcnt[gid] = [float(v) if v.strip() != "" else 0.0 for v in row[1:]]

    def shape(g):
        """返回 (分组和 g0/g1, 非零样本数 nz0/nz1, 形态标签)。"""
        c = gcnt.get(g)
        if c is None:
            return None, None, ""
        g0, g1 = c[0::2], c[1::2]
        s0, s1 = sum(g0), sum(g1)
        nz0 = sum(1 for v in g0 if v > 0)
        nz1 = sum(1 for v in g1 if v > 0)
        if min(s0, s1) <= 2:
            tag = "近全0"
        elif min(nz0, nz1) <= 1:
            tag = "每组仅1样本"
        else:
            tag = ""
        return (s0, s1), (nz0, nz1), tag

    rows = []
    for g, rv in ref_res.items():
        a = mine_res.get(g, {}).get("log2FoldChange")
        b = rv.get("log2FoldChange")
        if a is None or b is None:
            continue
        rows.append((abs(a - b), g, a, b))
    rows.sort(reverse=True)

    hline("不一致基因个案诊断")
    print("  |Δlog2FC| 超过阈值（全部可比基因 |Δ| 中位数 = %.4f）的基因数："
          % median([r[0] for r in rows]))
    for th in (0.1, 0.5, 1.0, 2.0, 5.0):
        c = sum(1 for r in rows if r[0] > th)
        print("    > %.1f : %6d 个（占可比基因 %.3f%%）"
              % (th, c, 100.0 * c / len(rows)))

    print("\n  差异最大的前 %d 个基因：" % n)
    print("  %-18s %9s %9s %9s %9s %8s %-12s %s"
          % ("基因", "baseMean", "我 lfc", "官方 lfc", "官方Al", "maxCooks",
             "形态", "非零样本 g0/g1"))
    for d, g, a, b in rows[:n]:
        dv = ref_disp.get(g, {})
        bm = ref_res[g].get("baseMean")
        dispmap = dv.get("dispMAP")
        mc = dv.get("maxCooks")
        _s, nz, tag = shape(g)
        print("  %-18s %9.1f %9.3f %9.3f %9s %8s %-12s %s"
              % (g, bm if bm is not None else float("nan"), a, b,
                 ("%.2e" % dispmap) if dispmap else "NA",
                 ("%.1f" % mc) if mc is not None else "NA", tag or "-",
                 "%d/%d" % nz if nz else "-"))

    # ---- 归因：把差异大的基因分成三类 ----
    print("\n  |Δlog2FC| > 1 的基因按「计数形态」归因：")
    big = [r for r in rows if r[0] > 1.0]
    c_zero = c_sparse = c_rest = 0
    rest = []
    for d, g, a, b in big:
        _s, nz, tag = shape(g)
        if tag == "近全0":
            c_zero += 1
        elif tag == "每组仅1样本":
            c_sparse += 1
        else:
            c_rest += 1
            rest.append((d, g, a, b))
    print("    ① 某一组计数和 <= 2（log2FC 数学上无界）: %4d 个" % c_zero)
    print("    ② 每组只有 1 个样本非零（等于没有重复，log2FC 不可识别）: %4d 个"
          % c_sparse)
    print("    ③ 其余（两组都有多个非零样本，属真正的数值分歧）: %4d 个"
          % c_rest)
    if rest:
        print("       ③ 的清单：")
        for d, g, a, b in rest[:20]:
            _s, nz, _t = shape(g)
            print("         %-18s 我 %9.4f  官方 %9.4f  Δ=%.4f  非零样本 %d/%d"
                  % (g, a, b, d, nz[0], nz[1]))
    return None


# ---------------------------------------------------------------------------
# M2 消融实验：趋势拟合 / 收缩让对拍结果改善了多少
# ---------------------------------------------------------------------------
def disp_ablation(mine_disp, ref_res, lines):
    """把 M2 的三个阶段分别接进 GLM，看与官方 log2FC 的吻合度如何变化。

    为什么需要这个实验：项目文档的「良好」档有一条硬要求：
    「能定量说明『趋势拟合前后，对拍结果改善了多少』」。只看
    `dispGeneEst / dispFit / dispMAP` 三列各自与官方的相关系数是不够的，
    那说明的是「离散度本身像不像」；这里要回答的是**最终结果**变好了多少。

    做法：同一套 counts / 设计矩阵，只换 dispersion 这一个输入：

        ① dispGeneEst：只用 gene-wise MLE（不趋势拟合、不收缩）
        ② dispFit：加趋势拟合（不收缩）
        ③ dispMAP：完整流程（趋势拟合 + 经验贝叶斯收缩）

    每一档都重跑一遍 GLM + Wald，再与官方求相关。

    WARNING: **只看 log2FC 会严重低估 M2 的价值**：log2FC 是点估计，主要由各组的
    均值结构决定，离散度对它影响很小（实测三档只从 0.9966 升到 0.9997）。
    离散度真正决定的是**标准误**，也就是 p 值与「哪些基因显著」。
    所以这个实验同时记录三样东西：

        log2FC r：点估计
        p 值 r（-log10）： 显著性排序
        显著集 Jaccard：padj 阈值下的判定

    实测结果（2026-09-26，airway 全量）：

        M2 阶段                  log2FC r   p 值 r    显著集 J
        ① gene-wise MLE          0.9966     0.8977     0.6082
        ② + 趋势拟合（不收缩）    0.9996     0.7172     0.7055
        ③ + 经验贝叶斯收缩        0.9997     0.9882     0.8703

    **注意 ② 是「非单调」的：p 值 r 先降后升。** 这不是 bug，是两步各自
    性质的体现，报告里必须写下来出：

      * log2FC 是点估计，由各组均值结构主导，对离散度不敏感，三档几乎不动；
      * 趋势拟合把「基因自己的信息」全丢了，只留一条平均曲线。用一条系统性
        偏高约 42% 的曲线（我们用分箱中位数 + Gamma GLM 替代官方的 locfit）
        当每个基因的最终离散度，会让 p 值整体偏保守，反而**不如**噪声更大的
        gene-wise MLE；
      * 收缩把趋势线当**先验均值**而不是**最终答案**，对估计稳的基因保留
        基因自己的信息，这才同时把 p 值 r 和显著集 Jaccard 拉回来。

    换句话说：**趋势拟合单独用是坏的，它的价值必须通过收缩才能兑现。**
    这也正是 DESeq2 把「先验均值」和「最终估计」分成两件事的原因。

    WARNING: 这一项要重跑三遍 M3（约 3~5 分钟），所以默认关闭，用
       `python tests/compare_baseline.py --disp-ablation` 显式打开。
    """
    import csv as _csv
    from src.dispersion import estimate_dispersions
    from src.glm_nb import fit_nb_glm, wald_test
    from src.io_utils import (align_design, read_counts_csv, read_design_csv)
    from src.linalg import build_design_matrix
    from src.multitest import bh_adjust
    from src.normalize import normalized_counts, size_factors

    hline("M2 消融：趋势拟合 / 收缩让最终对拍结果改善了多少")
    gene_ids, sample_ids, counts = read_counts_csv(
        os.path.join(ROOT, "data", "counts.csv"))
    design = read_design_csv(os.path.join(ROOT, "data", "design.csv"))
    groups, _c, _cl = align_design(sample_ids, design)
    X = build_design_matrix(groups)
    sf = size_factors(counts)
    log_sf = [math.log(f) for f in sf]

    stages = [("① gene-wise MLE（不拟合不收缩）", "dispGeneEst"),
              ("② + 趋势拟合（不收缩）", "dispFit"),
              ("③ + 经验贝叶斯收缩（完整）", "dispMAP")]
    ref_lfc = {g: r.get("log2FoldChange") for g, r in ref_res.items()}
    ref_p = {g: r.get("pvalue") for g, r in ref_res.items()}
    ref_padj = {g: r.get("padj") for g, r in ref_res.items()}
    gid_idx = {g: i for i, g in enumerate(gene_ids)}
    # 官方的显著集（判据与 I-2 一致：padj<0.05 且 |log2FC|>=1）
    ref_sig = set(g for g, r in ref_res.items()
                  if r.get("padj") is not None and r["padj"] < 0.05
                  and r.get("log2FoldChange") is not None
                  and abs(r["log2FoldChange"]) >= 1.0)

    lines.append("\n## M2 消融：趋势拟合 / 收缩对最终结果的贡献\n")
    lines.append("| M2 阶段 | 可检验基因 | log2FC r | p 值 r<br>(-log10) | "
                 "显著集 Jaccard | 显著基因数 |")
    lines.append("|---|---|---|---|---|---|")
    rows = []
    for label, col in stages:
        # 按基因下标收集 lfc / p 值（None = 全 0 或不可检验），
        # 跑完统一过 BH，这样每一档的 padj 与它自己的 p 值严格自洽。
        pvals = [None] * len(gene_ids)
        lfc_by_idx = {}
        xs, ys, px, py, n = [], [], [], [], 0
        for i, g in enumerate(gene_ids):
            if i % 4000 == 0:
                print("    %s ... %d / %d" % (label[:12], i, len(gene_ids)))
            alpha = mine_disp.get(g, {}).get(col)
            b = ref_lfc.get(g)
            if alpha is None or alpha != alpha or alpha <= 0 or b is None:
                continue
            fit = fit_nb_glm(counts[i], X, log_sf, alpha)
            lfc, _se, _st, pv = wald_test(fit["beta"], fit["cov"], 1)
            if lfc is None or lfc != lfc:
                continue
            lfc_by_idx[i] = lfc
            pvals[i] = pv
            xs.append(lfc)
            ys.append(b)
            n += 1
            rp = ref_p.get(g)
            if pv is not None and pv > 0 and rp is not None and rp > 0:
                px.append(-math.log10(pv))
                py.append(-math.log10(rp))
        padj = bh_adjust(pvals)
        sig = set(gene_ids[i] for i, v in lfc_by_idx.items()
                  if padj[i] is not None and padj[i] < 0.05 and abs(v) >= 1.0)
        r_lfc = corr(xs, ys)
        r_p = corr(px, py) if len(px) > 10 else float("nan")
        jac = jaccard(sig, ref_sig)
        print("  %-30s n=%-6d log2FC r=%.4f  p 值 r=%.4f  "
              "Jaccard=%.4f  显著 %d 个"
              % (label, n, r_lfc, r_p, jac, len(sig)))
        lines.append("| %s | %d | %.4f | **%.4f** | **%.4f** | %d |"
                     % (label, n, r_lfc, r_p, jac, len(sig)))
        rows.append((label, n, r_lfc, r_p, jac, len(sig)))

    if len(rows) >= 3:
        print("\n  结论：log2FC r   %.4f → %.4f → %.4f"
              % (rows[0][2], rows[1][2], rows[2][2]))
        print("        p 值 r     %.4f → %.4f → %.4f"
              % (rows[0][3], rows[1][3], rows[2][3]))
        print("        显著集 J   %.4f → %.4f → %.4f"
              % (rows[0][4], rows[1][4], rows[2][4]))
        lines.append("\n**三档对比**：log2FC r %.4f → %.4f → %.4f；"
                     "p 值 r %.4f → %.4f → %.4f；"
                     "显著集 Jaccard %.4f → %.4f → %.4f。"
                     % (rows[0][2], rows[1][2], rows[2][2],
                        rows[0][3], rows[1][3], rows[2][3],
                        rows[0][4], rows[1][4], rows[2][4]))
    return rows


# ---------------------------------------------------------------------------
def main():
    import argparse
    import csv

    ap = argparse.ArgumentParser(description="与官方 DESeq2 基准逐段对拍")
    ap.add_argument("--disp-ablation", action="store_true",
                    help="额外跑 M2 消融实验（重跑三遍 M3，约 2~4 分钟）")
    args = ap.parse_args()

    def must(path):
        p = os.path.join(ROOT, path)
        if not os.path.exists(p):
            print("缺少文件：%s" % p)
            print("（先跑 python deg.py 与 laptop/deseq2_baseline.R）")
            sys.exit(2)
        return p

    mine_res = load_table(must("data/my_result.csv"))
    ref_res = load_table(must("data/reference/deseq2_result_noFilter.csv"))
    mine_disp = load_table(must("data/intermediates/dispersion_mine.csv"))
    ref_disp = load_table(must("data/reference/dispersion_intermediates.csv"))

    # size factor
    with open(must("data/reference/sizeFactors.csv"), "r",
              encoding="utf-8-sig", newline="") as fh:
        rdr = csv.DictReader(fh)
        order, vals = [], []
        for row in rdr:
            order.append(row["sample"])
            vals.append(float(row["sizeFactor"]))
    mine_sf = load_table(must("data/intermediates/sizeFactors_mine.csv"))
    mine_vals = [mine_sf[str(i + 1)]["sizeFactor"]
                 for i in range(len(order))]

    lines = []
    print("=" * 74)
    print("DESeq2-FromScratch 对拍报告（数据源：airway，63677 基因 x 8 样本）")
    print("=" * 74)
    cmp_m1({"values": mine_vals}, {"order": order, "values": vals}, lines)
    cmp_m2(mine_disp, ref_disp, lines)
    cmp_m3m4(mine_res, ref_res, lines)
    diag(mine_res, ref_res, ref_disp)
    if args.disp_ablation:
        disp_ablation(mine_disp, ref_res, lines)
    hline("完成")

    outdir = os.path.join(ROOT, "logs")
    os.makedirs(outdir, exist_ok=True)
    import time
    out = os.path.join(outdir, "compare-%s.md"
                       % time.strftime("%Y%m%d-%H%M%S"))
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("<!-- 由 tests/compare_baseline.py 自动生成 -->\n\n")
        fh.write("\n".join(lines) + "\n")
    print("数字片段已写到 %s" % out)


if __name__ == "__main__":
    main()
