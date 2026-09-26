"""M3：负二项 GLM 的拟合（IRLS / Fisher scoring）与 Wald 检验。

模型
----
对每个基因独立拟合：

    counts[j] ~ NegativeBinomial(mean = mu[j], Var = mu[j] + alpha * mu[j]^2)
    log(mu[j]) = log(sizeFactor[j]) + X[j] · beta

X 是设计矩阵。本项目只做两组比较，所以 X 只有两列：

    X = [1, group]      group = 1 表示处理组，0 表示对照组

于是 beta[1] 就是 log(处理组均值 / 对照组均值)，也就是对数倍数变化。
DESeq2 的结果表里 log2FoldChange 是 beta[1]/ln2（见 R/fitNbinomGLMs.R 第 194 行：
    betaMatrix <- log2(exp(1)) * betaRes$beta_mat
），本模块同样在最后一步换算到 log2 尺度。

为什么必须迭代
--------------
GLM 没有解析解。这里用 Fisher scoring（IRLS 的一种）：

    1. eta = log(sf) + X·beta,  mu = exp(eta)（并下限截断到 minmu）
    2. 权重       w[j] = mu[j] / (1 + alpha * mu[j])
       （因为 dmu/deta = mu，Var = mu + alpha*mu^2，故 (dmu/deta)^2/Var = mu/(1+alpha*mu)）
    3. 工作响应   z[j] = eta[j] + (y[j] - mu[j]) / mu[j]
    4. 解加权最小二乘  (X' W X) beta_new = X' W z
    5. 若偏差（deviance）反而变大，则步长减半重试（step halving）

停止准则与 DESeq2 一致（见 R/fitNbinomGLMs.R 第 14–15 行注释）：
    abs(dev - dev_old) / (abs(dev) + 0.1) < betaTol        betaTol = 1e-8

Wald 检验
--------
    统计量  stat = beta / SE(beta)
    p 值    p = 2 * (1 - Phi(|stat|))     （标准正态，双侧）
依据 R/core.R 的 nbinomWaldTest：
    WaldPvalue <- 2 * pnorm(abs(WaldStatistic), lower.tail = FALSE)

阈值检验（论文 "Composite null hypotheses"）
------------------------------------------------
* 找「效应量显著超过 θ」的基因（H0: |beta| <= theta）：
      p = min(1, 2 * (1 - Phi((|beta| - theta) / SE)))
  θ = 0 时与上面的标准双侧 p 值完全一致。
* 找「效应量显著弱于 θ」的基因（H0: |beta| >= theta）：
      p = max( Phi((beta - theta)/SE),  1 - Phi((beta + theta)/SE) )
  注意：该检验要求 beta 是 MLE（未收缩），否则零中心先验会偏袒备择假设。

似然比检验（LRT）
-----------------
    统计量  stat = 2 * (ll_full - ll_reduced)  ~ 卡方(df)，df = 两组设计矩阵列数差
    p 值    p = P(chi2_df > stat)
官方 nbinomLRT 同样对两个模型使用相同的（最终）离散度 alpha。
"""

import math

from .linalg import build_design_matrix, invert, mat_vec, solve
from .nbmath import chi2_sf, normal_cdf, normal_two_sided_p

# 与 DESeq2 fitNbinomGLMs 的默认值一致
BETA_TOL = 1e-8
MAXIT = 100
MINMU = 0.5
# DESeq2 fitBeta（src/DESeq2.cpp）里的 `double large = 30.0`
LARGE_BETA = 30.0
LOG2E = 1.0 / math.log(2.0)


def nb_deviance(counts, mu, alpha):
    """负二项分布的偏差（deviance）。

    D = 2 * sum_j [ y_j * log(y_j / mu_j)
                    - (y_j + 1/alpha) * log((y_j + 1/alpha) / (mu_j + 1/alpha)) ]
    其中 y_j = 0 时第一项取 0（极限）。
    """
    size = 1.0 / alpha
    total = 0.0
    for j in range(len(counts)):
        y = counts[j]
        m = mu[j]
        if y > 0:
            total += y * math.log(y / m)
        total -= (y + size) * math.log((y + size) / (m + size))
    return 2.0 * total


def _initial_beta(counts, X, log_sf):
    """用 log((y + 0.5) / sf) 做普通最小二乘，得到 beta 的初值。

    这是一个标准的 GLM 初值策略：对数尺度上先把量级摆对，
    后面的 Fisher scoring 只需要少量迭代即可收敛。

    实测对比（airway 全量）：用 OLS 初值时 baseMean >= 10 的基因与官方
    log2FC 相关系数 0.997；改成官方那种「从 0 出发」反而掉到 0.72，
    且符号一致率从 97% 掉到 68%（因为从 0 出发时，退化基因会在离真值
    很远的地方提前停下）。所以这里保留 OLS 初值，属于**明示差异**。
    """
    m = len(counts)
    p = len(X[0])
    z = [math.log((counts[j] + 0.5) / math.exp(log_sf[j])) for j in range(m)]
    # 正规方程 (X'X) b = X'z
    XtX = [[sum(X[j][a] * X[j][b] for j in range(m)) for b in range(p)]
           for a in range(p)]
    Xtz = [sum(X[j][a] * z[j] for j in range(m)) for a in range(p)]
    return solve(XtX, Xtz)


def fit_nb_glm(counts, X, log_sf, alpha,
               maxit=MAXIT, beta_tol=BETA_TOL, minmu=MINMU, beta_init=None):
    """对单个基因拟合负二项 GLM。

    参数
    ----
    counts  : 该基因在各样本的原始计数
    X       : 设计矩阵（m × p）
    log_sf  : 各样本 size factor 的自然对数（作为 offset）
    alpha   : 该基因的离散度（固定值，M3 不在这里估）
    minmu   : mu 的下限，与 DESeq2 默认 0.5 一致

    返回
    ----
    dict，含 beta（ln 尺度）、cov（系数的协方差矩阵）、mu、deviance、
          iterations、converged、loglik
    """
    m = len(counts)
    p = len(X[0])

    beta = list(beta_init) if beta_init is not None else _initial_beta(counts, X, log_sf)

    dev_old = None
    converged = False
    iteration = 0
    mu = None

    for iteration in range(1, maxit + 1):
        # ---- 1. 计算 mu = sf * exp(X·beta)，并下限截断 ----
        eta = mat_vec(X, beta)
        mu = []
        for j in range(m):
            e = math.exp(min(eta[j] + log_sf[j], 700.0))
            mu.append(e if e > minmu else minmu)

        # ---- 2. 权重 ----
        weights = [mu[j] / (1.0 + alpha * mu[j]) for j in range(m)]

        # ---- 3. 工作响应 ----
        # WARNING: 对数项必须用**截断之后**的 mu，不能直接用 η。
        # 官方 fitBeta 写的是：
        #     z = arma::log(mu_hat / nfrow) + (yrow - mu_hat) / mu_hat;
        # 其中 mu_hat 已经过 minmu 截断过。
        # 若写成 η + log_sf（= log(未截断的 mu)），一旦某样本的 mu 触底
        # （minmu = 0.5），偏差就不再随 β 变化、但 z 仍随 η 线性下移，
        # IRLS 会沿着这个"偏差不再改善"的方向把 β 一路推下去，实测能推到
        # |β| = 30（log2FC = ±43），而官方结果里 |log2FC| 最大只有 9.61。
        # 用 log(mu) - log_sf 之后，mu 触底时 z 自动停止移动，问题消失。
        z = [(math.log(mu[j]) - log_sf[j]) + (counts[j] - mu[j]) / mu[j]
             for j in range(m)]

        # ---- 4. 解加权最小二乘 ----
        XtWX = [[sum(X[j][a] * weights[j] * X[j][b] for j in range(m))
                 for b in range(p)] for a in range(p)]
        XtWz = [sum(X[j][a] * weights[j] * z[j] for j in range(m))
                for a in range(p)]
        try:
            beta_new = solve(XtWX, XtWz)
        except ValueError:
            # 矩阵奇异（极端数据），保持上一步结果并停止
            break

        # ---- 5. step halving：偏差不许变大 ----
        dev_prev = nb_deviance(counts, mu, alpha)
        step = 1.0
        for _ in range(10):
            beta_try = [beta[a] + step * (beta_new[a] - beta[a]) for a in range(p)]
            eta_try = mat_vec(X, beta_try)
            mu_try = []
            for j in range(m):
                e = math.exp(min(eta_try[j] + log_sf[j], 700.0))
                mu_try.append(e if e > minmu else minmu)
            dev_try = nb_deviance(counts, mu_try, alpha)
            if dev_try <= dev_prev or step < 1e-4:
                break
            step *= 0.5

        # ---- 5b. 系数上界（DESeq2 fitBeta 的 large = 30）----
        # 官方 src/DESeq2.cpp 的 fitBeta 里有：
        #     double large = 30.0;
        #     solve(beta_hat, r, gamma_hat);
        #     if (sum(abs(beta_hat) > large) > 0) { iter(i) = maxit; break; }
        # 也就是说 |β| 一旦超过 30，官方认为这一步不可信、直接停止迭代。
        # 这条保护是必需的：对「某一组几乎全是 0」的基因，负二项似然在 β 的
        # 某个方向上是**无界**的（我们验证过：目标函数会一直下降，
        # 偏差从 220 一路降到 190），没有这个上界，IRLS 会跑到 β = -47
        # （log2FC = -68），而官方结果里 |log2FC| 最大只有 9.61。
        # 我们的处理与官方一致地「不接受这一步」，但更保守一点：
        # 官方 break 时 beta_hat 已是越界值，我们保留上一轮的 β。
        if max(abs(beta_try[a]) for a in range(p)) > LARGE_BETA:
            converged = False
            break

        # ---- 6. 收敛判定 ----
        # WARNING:WARNING: 顺序极重要：delta_beta 必须在「更新 beta」**之前**算。
        # 早先写成先 `beta = beta_try` 再 `delta_beta = |beta_try - beta|`，
        # 结果 delta_beta 恒等于 0，第二个判据永远成立，所有基因都在第 2 次
        # 迭代就退出，等于没有收敛判定。多数基因因为初值已经不错（对数尺度
        # OLS）而看不出差别，但在极端过散布的基因上会一路跑到 beta = -47
        # （log2FC = -68），真实数据上一跑就暴露了。
        delta_beta = max(abs(beta_try[a] - beta[a]) for a in range(p))
        beta = beta_try
        mu = mu_try
        dev = dev_try

        # (a) DESeq2 的官方准则（R/fitNbinomGLMs.R 第 14-15 行）：
        #     abs(dev - dev_old) / (abs(dev) + 0.1) < betaTol
        # (b) 补充准则：系数本身几乎不再变化。
        #     只有 (a) 是不够的，当偏差本身很小（低计数基因），分母里的 +0.1
        #     会让 (a) 实际上要求 1e-10 量级的绝对精度，double 精度下永远达不到，
        #     导致大量基因白白迭代到上限。单元测试里记录了这一点。
        if dev_old is not None:
            rel_dev = abs(dev - dev_old) / (abs(dev) + 0.1)
            if rel_dev < beta_tol or delta_beta < beta_tol:
                converged = True
                break
        dev_old = dev

    # 收敛后的协方差矩阵 = (X' W X)^-1
    weights = [mu[j] / (1.0 + alpha * mu[j]) for j in range(m)]
    XtWX = [[sum(X[j][a] * weights[j] * X[j][b] for j in range(m))
             for b in range(p)] for a in range(p)]
    try:
        cov = invert(XtWX)
    except ValueError:
        cov = [[float("nan")] * p for _ in range(p)]

    loglik = _nb_loglik(counts, mu, alpha)

    return {
        "beta": beta,
        "cov": cov,
        "mu": mu,
        "deviance": nb_deviance(counts, mu, alpha),
        "iterations": iteration,
        "converged": converged,
        "loglik": loglik,
    }


def _nb_loglik(counts, mu, alpha):
    """内联实现，避免每次调用都做函数查找（热路径）。"""
    size = 1.0 / alpha
    total = 0.0
    lgamma = math.lgamma
    log = math.log
    for j in range(len(counts)):
        y = counts[j]
        m = mu[j]
        p = size / (size + m)
        total += lgamma(y + size) - lgamma(size) - lgamma(y + 1.0) \
                 + size * log(p) + y * math.log1p(-p)
    return total


def wald_test(beta, cov, coef_index=1):
    """对某个系数做 Wald 检验，返回 (log2FoldChange, lfcSE, stat, pvalue)。

    log2FC 与 lfcSE 换算到 log2 尺度：
        log2FC = beta / ln2,   lfcSE = sqrt(var) / ln2
    （依据 R/fitNbinomGLMs.R 第 194/198 行的 log2(exp(1)) 因子）
    """
    b = beta[coef_index]
    var = cov[coef_index][coef_index]
    if var is None or var != var or var <= 0.0:
        return (b * LOG2E, None, None, None)
    se_ln = math.sqrt(var)
    stat = b / se_ln
    return (b * LOG2E, se_ln * LOG2E, stat, normal_two_sided_p(stat))


def _coef_se(beta, cov, coef_index):
    """取系数的点估计与标准误（ln 尺度）；协方差无效时返回 (b, None)。"""
    b = beta[coef_index]
    var = cov[coef_index][coef_index]
    if var is None or var != var or var <= 0.0:
        return b, None
    return b, math.sqrt(var)


def wald_test_above_threshold(beta, cov, theta, coef_index=1):
    """阈值检验：H0 为 |beta| <= theta（找效应量「显著超过」阈值的基因）。

    返回 (log2FoldChange, lfcSE, pvalue)。公式（论文 Composite null hypotheses）：
        p = min(1, 2 * (1 - Phi((|beta| - theta) / SE)))
    theta = 0 时退化为标准双侧 Wald p 值（数值上完全相等）。
    """
    b, se_ln = _coef_se(beta, cov, coef_index)
    if se_ln is None:
        return (b * LOG2E, None, None)
    z = (abs(b) - theta) / se_ln
    p = 2.0 * (1.0 - normal_cdf(z))
    if p > 1.0:
        p = 1.0
    return (b * LOG2E, se_ln * LOG2E, p)


def wald_test_below_threshold(beta, cov, theta, coef_index=1):
    """阈值检验：H0 为 |beta| >= theta（找证据表明效应量「弱于」阈值的基因）。

    返回 (log2FoldChange, lfcSE, pvalue)。公式（论文 Composite null hypotheses）：
        p1 = Phi((beta - theta) / SE)
        p2 = 1 - Phi((beta + theta) / SE)
        p = max(p1, p2)

    WARNING: 该检验必须使用 **MLE** 的 beta（未收缩）：零中心的 LFC 先验本身
    偏袒「效应量小」的备择假设，用了收缩估计就不是数据在说话。
    本项目的主流程不启用 LFC 先验，所以直接可用。
    """
    b, se_ln = _coef_se(beta, cov, coef_index)
    if se_ln is None:
        return (b * LOG2E, None, None)
    p1 = normal_cdf((b - theta) / se_ln)
    p2 = 1.0 - normal_cdf((b + theta) / se_ln)
    p = p1 if p1 > p2 else p2
    return (b * LOG2E, se_ln * LOG2E, p)


def lrt_test(counts, X_full, X_reduced, log_sf, alpha,
             maxit=MAXIT, beta_tol=BETA_TOL, minmu=MINMU):
    """似然比检验：比较全模型与简化模型的负二项对数似然。

    参数
    ----
    X_full    : 全设计矩阵（如 [1, group]）
    X_reduced : 简化设计矩阵（如 [1]，即去掉待检验的列）
    alpha     : 该基因的最终离散度（两个模型共用，与官方 nbinomLRT 一致）

    返回
    ----
    dict：stat / df / pvalue / loglik_full / loglik_reduced

    stat = 2 * (ll_full - ll_reduced)，在 H0 下近似服从卡方(df)，
    df = 两组设计矩阵的列数差（本项目两组设计下 df = 1）。
    """
    full = fit_nb_glm(counts, X_full, log_sf, alpha,
                      maxit=maxit, beta_tol=beta_tol, minmu=minmu)
    red = fit_nb_glm(counts, X_reduced, log_sf, alpha,
                     maxit=maxit, beta_tol=beta_tol, minmu=minmu)
    stat = 2.0 * (full["loglik"] - red["loglik"])
    if stat < 0.0:
        # 数值噪声：两个模型各自在收敛容差内停止，理论上 ll_full >= ll_reduced
        stat = 0.0
    df = len(X_full[0]) - len(X_reduced[0])
    p = chi2_sf(stat, df) if df > 0 else None
    return {
        "stat": stat,
        "df": df,
        "pvalue": p,
        "loglik_full": full["loglik"],
        "loglik_reduced": red["loglik"],
    }
