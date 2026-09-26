"""M3 的扩展诊断：Cook's 距离离群点检测。

来源：DESeq2 论文（Love et al. 2014）"Detection of count outliers" 一节。
Cook's 距离衡量「删掉某个样本重新拟合后，系数向量会移动多远」，
用来发现「单个样本主导了某个基因的 LFC 估计」这种情况。

官方 DESeq2 的行为（阈值 = qf(0.99, p, m-p)）：
  * 某组重复数 <= 6：把整个基因标记为离群（p 值置 NA，不参与后续分析）
  * 某组重复数 >= 7：用「截断均值（按 size factor 缩放）」替换离群计数后重估

**本项目只实现诊断部分**（计算每个样本的 Cook 距离 + 阈值判定），
不做自动替换/删除，那会改变结果表口径，与对拍基准（noFilter 口径）不一致。
需要过滤时由使用者在自己的分析里决定如何处理（这在报告里有说明）。

公式（Cook & Weisberg 1982；与官方 DESeq2 相同）：
    D_j = (R_j^2 / p) * h_jj / (1 - h_jj)^2
    R_j = (K_j - mu_j) / sqrt(V(mu_j)),      V(mu) = mu + alpha * mu^2
    H   = W^(1/2) X (X' W X)^(-1) X' W^(1/2),   w_j = mu_j / (1 + alpha * mu_j)
其中 mu 由「无 LFC 先验的负二项 GLM」拟合；alpha 论文要求用稳健矩估计
（α_rob = max((s_rob^2 - mean) / mean^2, 0)），本模块提供 robust_moments_dispersion
供调用方选择，也允许直接传入最终离散度。
"""

import math

from .glm_nb import fit_nb_glm
from .linalg import invert
from .nbmath import f_quantile, mad


def robust_moments_dispersion(norm_counts):
    """稳健矩估计离散度：alpha_rob = max((s_rob^2 - mean) / mean^2, 0)。

    norm_counts : 一个基因的**归一化**计数（原始计数 / size factor）
    s_rob       : 稳健方差，这里用 MAD 的平方（常数 1.4826 与 R 的 mad 一致）

    论文用于 Cook's 距离里的 Pearson 残差与权重（对离群样本稳健）。
    """
    n = len(norm_counts)
    if n < 2:
        return 0.0
    mean = sum(norm_counts) / n
    if mean <= 0.0:
        return 0.0
    s_rob = mad(norm_counts)
    v = (s_rob * s_rob - mean) / (mean * mean)
    return v if v > 0.0 else 0.0


def cooks_threshold(n_samples, n_coef, q=0.99):
    """判定阈值：F(p, m-p) 分布的 q 分位（官方用 qf(0.99, p, m-p)）。

    n_samples : m
    n_coef    : p（含截距）
    退化情形（m <= p）返回 inf（无法判定）。
    """
    if n_samples <= n_coef:
        return float("inf")
    return f_quantile(q, n_coef, n_samples - n_coef)


def cooks_distance(counts, X, log_sf, alpha, beta=None):
    """逐样本 Cook's 距离。返回 D_j 的列表（顺序与样本一致）。

    counts / X / log_sf : 该基因的计数、设计矩阵（m×p）、log size factor
    alpha               : 用于残差与权重的离散度（论文建议稳健矩估计）
    beta                : 可选，复用已拟合的系数；None 时内部重新拟合
    """
    m = len(counts)
    p = len(X[0])
    if beta is None:
        fit = fit_nb_glm(counts, X, log_sf, alpha)
        beta = fit["beta"]
    mu = []
    for j in range(m):
        e = math.exp(min(sum(X[j][a] * beta[a] for a in range(p)) + log_sf[j], 700.0))
        mu.append(e if e > 0.0 else 1e-300)

    # 权重 w_j = mu_j / (1 + alpha * mu_j)，头矩阵对角元 h_jj = w_j x_j' inv(X'WX) x_j
    w = [mu[j] / (1.0 + alpha * mu[j]) for j in range(m)]
    XtWX = [[sum(X[j][a] * w[j] * X[j][b] for j in range(m)) for b in range(p)]
            for a in range(p)]
    try:
        inv = invert(XtWX)
    except ValueError:
        return [float("nan")] * m

    out = []
    for j in range(m):
        # x_j' inv x_j
        tmp = [sum(X[j][a] * inv[a][b] for a in range(p)) for b in range(p)]
        quad = sum(tmp[b] * X[j][b] for b in range(p))
        h = w[j] * quad
        if h >= 1.0 - 1e-12:
            out.append(float("inf"))
            continue
        var_j = mu[j] + alpha * mu[j] * mu[j]
        r = (counts[j] - mu[j]) / math.sqrt(var_j)
        out.append((r * r / p) * h / ((1.0 - h) * (1.0 - h)))
    return out


def flag_cooks_outliers(counts, X, log_sf, alpha, q=0.99):
    """计算 Cook 距离并做阈值判定。

    返回 dict：
        cooks     : 逐样本 Cook 距离
        threshold : qf(q, p, m-p)
        flags     : 逐样本布尔标记（cooks > threshold）
        max_cooks : 该基因的最大 Cook 距离（官方结果表里的 maxCooks 列）
    """
    cook = cooks_distance(counts, X, log_sf, alpha)
    thr = cooks_threshold(len(counts), len(X[0]), q=q)
    flags = [(d == d) and d > thr for d in cook]
    mx = max([d for d in cook if d == d], default=float("nan"))
    return {"cooks": cook, "threshold": thr, "flags": flags, "max_cooks": mx}