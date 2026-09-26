"""M1：归一化（median-of-ratios，即 DESeq2 的 size factor）。

要解决的问题
------------
两个样本的测序深度不同时，计数不能直接比。样本 B 测到的总 read 是 A 的两倍，
那 B 里每个基因的计数天然就是 A 的两倍，与生物学无关。
size factor 就是每个样本的一个缩放系数，计数除以它之后，样本间才同尺度。

算法（按源码核自 DESeq2 官方源码 R/core.R 的 estimateSizeFactorsForMatrix）
------------------------------------------------------------------------
默认 type = "ratio"、locfunc = median：

    1. loggeomeans[g] = mean_j( log(counts[g, j]) )
       注意：只要某个样本里该基因为 0，log(0) = -Inf，行均值就是 -Inf，
       该基因会被第 3 步的 is.finite 判断自动排除。：这就是「几何均值遇到 0 会塌」在官方实现里的处理方式。

    2. 若所有基因的 loggeomeans 都是 -Inf -> 报错：
       "every gene contains at least one zero, cannot compute log geometric means"

    3. 对每个样本 j：
           sf[j] = exp( median( (log(counts[g, j]) - loggeomeans[g])
                                [有限 loggeomeans 且 counts[g, j] > 0] ) )

    4. 官方实现默认**不再**把 sf 归一化到几何均值 1
       （只有调用方自己传入 geoMeans 时才会走那段缩放）。

为什么取中位数而不是均值：绝大多数基因在两个样本间其实没有变化，用中位数
可以抗住少数真正差异表达的基因把整体比例带偏。

第 3 步里「先取对数比、再取中位数、再 exp」和「先取比值的中位数」是等价的，
因为 median 和 exp 都是单调变换。这里按官方写法先取对数。
"""

import math

from .nbmath import median


def log_geomeans(counts, n_samples):
    """第 1 步：每个基因跨样本的对数几何均值；含 0 的基因返回 -inf。"""
    out = []
    append = out.append
    for row in counts:
        s = 0.0
        ok = True
        for v in row:
            if v <= 0:
                ok = False
                break
            s += math.log(v)
        append(s / n_samples if ok else float("-inf"))
    return out


def size_factors(counts):
    """计算每个样本的 size factor。

    参数
    ----
    counts : list of list，counts[基因][样本]，非负整数

    返回
    ----
    list of float，长度 = 样本数
    """
    if not counts:
        raise ValueError("counts 为空")
    n_samples = len(counts[0])
    if n_samples == 0:
        raise ValueError("counts 没有样本列")
    for row in counts:
        if len(row) != n_samples:
            raise ValueError("counts 各行长度不一致")

    lg = log_geomeans(counts, n_samples)
    if all(v == float("-inf") for v in lg):
        raise ValueError(
            "每个基因都至少含一个 0，无法计算对数几何均值"
            "（与 DESeq2 的同名报错一致）")

    n_genes = len(counts)
    factors = []
    for j in range(n_samples):
        ratios = []
        append = ratios.append
        for g in range(n_genes):
            lg_g = lg[g]
            if lg_g == float("-inf"):
                continue
            v = counts[g][j]
            if v > 0:
                append(math.log(v) - lg_g)
        if not ratios:
            raise ValueError(
                "样本 %d 没有任何可用于估计 size factor 的基因"
                "（该样本全部为 0？）" % j)
        factors.append(math.exp(median(ratios)))
    return factors


def normalized_counts(counts, factors):
    """计数除以自身的 size factor，得到归一化计数。

    对应 DESeq2 的 counts(dds, normalized=TRUE)。
    """
    n_samples = len(factors)
    return [[row[j] / factors[j] for j in range(n_samples)] for row in counts]


def base_mean(norm_counts):
    """每个基因在所有样本上的归一化计数均值。

    对应 DESeq2 结果表里的 baseMean 列。
    """
    n_samples = len(norm_counts[0])
    return [sum(row) / n_samples for row in norm_counts]
