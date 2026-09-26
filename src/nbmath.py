"""负二项分布相关的纯标准库数学工具。

所有公式的实现依据（均已核对 DESeq2 官方源码，不是凭记忆）：

1. 负二项分布与离散度的参数化
   DESeq2 使用「均值-离散度」参数化：  Var = mu + alpha * mu^2
   等价于 R 的 dnbinom(mu = mu, size = 1/alpha)。
   见 DESeq2 源码 R/core.R 中的 nbinomLogLike()：
       rowSums(dnbinom(counts, mu = mu, size = 1/disp, log = TRUE))

2. 离散度的取值范围
   minDisp = 1e-8；搜索下界为 log(minDisp/10) = log(1e-9)。
   见 R/core.R estimateDispersionsGeneEst(minDisp = 1e-8) 及
   fitDisp(..., min_log_alpha = log(minDisp/10))。

3. Wald 检验的 p 值使用标准正态：
       WaldPvalue <- 2 * pnorm(abs(WaldStatistic), lower.tail = FALSE)
   见 R/core.R nbinomWaldTest（约第 1513 行）。
   本文件用 math.erfc 实现：
       2 * (1 - Phi(|z|)) = erfc(|z| / sqrt(2))
"""

import math

# 与 DESeq2 一致的最小离散度
MIN_DISP = 1e-8
# 离散度搜索的下界（DESeq2 用 log(minDisp/10)）
MIN_DISP_SEARCH = MIN_DISP / 10.0


# --------------------------------------------------------------------------
# 负二项分布
# --------------------------------------------------------------------------
def nb_logpmf(y, mu, size):
    """单个观测的负二项对数概率密度。

    参数
    ----
    y    : 观测计数（非负整数）
    mu   : 均值（> 0）
    size : R 的 size 参数，等于 1/alpha

    P(y) = C(y + size - 1, y) * p^size * (1-p)^y,  p = size / (size + mu)
    """
    # lgamma(y + size) - lgamma(size) - lgamma(y + 1)
    log_c = math.lgamma(y + size) - math.lgamma(size) - math.lgamma(y + 1.0)
    p = size / (size + mu)
    # size * log(p) + y * log(1-p)
    return log_c + size * math.log(p) + y * math.log1p(-p)


def nb_loglik(counts, mu, alpha):
    """一个基因在整个样本向量上的负二项对数似然之和。

    counts, mu : 等长的序列（长度 = 样本数）
    alpha      : 离散度（> 0），size = 1/alpha

    与 DESeq2 的 nbinomLogLike() 等价（不加权重的情形）。
    """
    size = 1.0 / alpha
    total = 0.0
    lgamma = math.lgamma
    log = math.log
    log1p = math.log1p
    for i in range(len(counts)):
        yi = counts[i]
        mui = mu[i]
        if mui <= 0.0:
            mui = 1e-300
        p = size / (size + mui)
        total += lgamma(yi + size) - lgamma(size) - lgamma(yi + 1.0) \
                 + size * log(p) + yi * log1p(-p)
    return total


# --------------------------------------------------------------------------
# 正态分布
# --------------------------------------------------------------------------
def normal_cdf(x):
    """标准正态分布函数 Phi(x) = 0.5 * erfc(-x / sqrt(2))。

    用于阈值检验（论文 "Composite null hypotheses"）里的单侧积分：
    P(Z <= z)。数值上与 R 的 pnorm(z) 一致。
    """
    return 0.5 * math.erfc(-x / math.sqrt(2.0))


def normal_two_sided_p(z):
    """双侧 p 值 2 * (1 - Phi(|z|))，等价于 R 的 2*pnorm(abs(z), lower.tail=FALSE)。"""
    if z is None:
        return None
    z = abs(z)
    if math.isinf(z):
        return 0.0
    # 2 * 0.5 * erfc(z / sqrt(2)) == erfc(z / sqrt(2))
    return math.erfc(z / math.sqrt(2.0))


# --------------------------------------------------------------------------
# 特殊函数
# --------------------------------------------------------------------------
def trigamma(x):
    """三伽马函数 psi'(x)（digamma 的导数）。

    用递推把 x 抬到 6 以上，再用渐近级数：
        psi'(x) ~ 1/x + 1/(2x^2) + 1/(6x^3) - 1/(30x^5) + 1/(42x^7)
                  - 1/(30x^9) + 5/(66x^11)
    DESeq2 用它计算基因级离散度估计的期望抽样方差：
        expVarLogDisp <- trigamma((m - p)/2)
    """
    if x <= 0:
        raise ValueError("trigamma 要求 x > 0")
    result = 0.0
    while x < 6.0:
        result += 1.0 / (x * x)
        x += 1.0
    inv = 1.0 / x
    inv2 = inv * inv
    return result + inv * (
        1.0
        + 0.5 * inv
        + inv2 / 6.0
        - inv2 * inv2 / 30.0
        + inv2 ** 3 / 42.0
        - inv2 ** 4 / 30.0
        + 5.0 * inv2 ** 5 / 66.0
    )


def _gser(a, x, itmax=300, eps=3e-14):
    """正则化下不完全 gamma 函数 P(a, x)（级数展开）。

    算法取自 Numerical Recipes（gser），用于 chi2_sf 的小 x 分支。
    """
    ap = a
    total = 1.0 / a
    term = total
    for _ in range(itmax):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * eps:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gcf(a, x, itmax=300, eps=3e-14):
    """正则化上不完全 gamma 函数 Q(a, x)（连分式展开）。

    算法取自 Numerical Recipes（gcf），用于 chi2_sf 的大 x 分支。
    """
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, itmax + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def chi2_sf(x, df):
    """卡方分布上尾概率 P(X^2_df > x)，即 R 的 pchisq(x, df, lower.tail=FALSE)。

    用于似然比检验（LRT）：stat = 2*(ll_full - ll_reduced) ~ 卡方(df)。
    实现 = 正则化上不完全 gamma Q(df/2, x/2)。
    已知闭式可用于自检：df=1 时等于 erfc(sqrt(x/2))；df=2 时等于 exp(-x/2)。
    """
    if x <= 0.0:
        return 1.0
    if df <= 0:
        raise ValueError("chi2_sf 要求 df > 0")
    a = df / 2.0
    xx = x / 2.0
    if xx < a + 1.0:
        return 1.0 - _gser(a, xx)
    return _gcf(a, xx)


def _betacf(a, b, x, itmax=300, eps=3e-14):
    """正则化不完全 beta 函数的连分式部分（Numerical Recipes betacf）。"""
    tiny = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betai(a, b, x):
    """正则化不完全 beta 函数 I_x(a, b)，即 R 的 pbeta(x, a, b)。

    用于 F 分布的累积函数：Cook's 距离的 0.99 分位阈值
    （官方 DESeq2 用 qf(0.99, p, m-p) 作为 cooksCutoff）。
    """
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def f_cdf(x, d1, d2):
    """F 分布累积函数 P(F_{d1,d2} <= x)，即 R 的 pf(x, d1, d2)。"""
    if x <= 0.0:
        return 0.0
    return betai(d1 / 2.0, d2 / 2.0, d1 * x / (d1 * x + d2))


def f_quantile(q, d1, d2):
    """F 分布分位数（R 的 qf(q, d1, d2)），用二分法在累积函数上求逆。

    Cook's 距离的判定阈值：qf(0.99, p, m-p)。
    """
    if not (0.0 < q < 1.0):
        raise ValueError("f_quantile 要求 0 < q < 1")
    lo, hi = 0.0, 1.0
    while f_cdf(hi, d1, d2) < q:
        hi *= 2.0
        if hi > 1e300:
            raise ValueError("f_quantile 未收敛")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f_cdf(mid, d1, d2) < q:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12 * max(1.0, hi):
            break
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------
# 稳健统计量
# --------------------------------------------------------------------------
def median(values):
    """中位数（不依赖 statistics，避免 Import 顺序问题，也便于就地优化）。"""
    vals = sorted(values)
    n = len(vals)
    if n == 0:
        raise ValueError("median 需要至少一个元素")
    mid = n // 2
    if n % 2 == 1:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


def mad(values, center=None):
    """中位数绝对偏差，常数与 R 的 mad() 默认值一致（1.4826）。

    实现依据 DESeq2 源码 R/methods.R dispFun.replace()：
        varLogDispEsts <- mad(dispResiduals[aboveMinDisp], na.rm=TRUE)^2
    R 的 mad() 默认 constant = 1.4826，center = median(x)。
    """
    vals = [v for v in values if v is not None and not math.isnan(v)]
    if not vals:
        raise ValueError("mad 需要至少一个有效元素")
    if center is None:
        center = median(vals)
    return 1.4826 * median([abs(v - center) for v in vals])


def sample_variance(values):
    """样本方差（分母 n-1），等价于 R 的 var() / rowVars()。"""
    n = len(values)
    if n < 2:
        return float("nan")
    m = sum(values) / n
    return sum((v - m) ** 2 for v in values) / (n - 1)
