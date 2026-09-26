"""M2：离散度估计与经验贝叶斯收缩。

离散度（dispersion，记作 alpha）描述「这批计数比泊松分布乱多少」：

    Var = mu + alpha * mu^2

alpha 趋近 0 就退化成泊松分布。它之所以是项目的核心难点，是因为**只用每个基因
自己那几个重复去估 alpha，噪声极大**：有的基因碰巧重复间很接近，方差被严重低估，
p 值随后假阳性爆炸。

DESeq2 的解法是把上万个基因放在一起看，分三步（本模块严格对应这三步）：

    2a  gene-wise：每个基因先用矩估计/粗略估计定起点，再做最大似然得到 dispGeneEst
    2b  trend    ：把 dispGeneEst 对 baseMean 拟合一条曲线（默认 y = a + b/mean），
                    得到 dispFit
    2c  MAP      ：把每个基因的估计朝这条曲线收缩，得到最终的 dispersion

公式依据（逐条核对自 DESeq2 官方源码，不是凭记忆）
------------------------------------------------
* 矩估计（momentsDispEstimate，R/core.R:2459）
      xim = mean(1/sizeFactors)
      momentsDisp = (baseVar - xim * baseMean) / baseMean^2
* 粗略估计（roughDispEstimate，R/core.R:2442）
      mu = 线性模型拟合的归一化计数（下限 1）
      roughDisp = rowSums(((y - mu)^2 - mu) / mu^2) / (m - p)
      alpha 起点 = pmin(roughDisp, momentsDisp)
* 离散度方差-均值的趋势（parametricDispersionFit，R/core.R:2186）
      y = asymptDisp + extraPois / mean，用 Gamma(link="identity") GLM 迭代拟合，
      初值 (0.1, 1)，每轮剔除残差不在 (1e-4, 15) 内的基因，最多 10 轮
* 先验方差（R/methods.R:176-184 与 R/core.R:1136-1208）
      dispResiduals = log(dispGeneEst) - log(dispFit)
      varLogDispEsts = mad(dispResiduals[dispGeneEst >= 1e-6])^2
      dispPriorVar = max(varLogDispEsts - trigamma((m - p)/2), 0.25)
      WARNING: 特例：当残差自由度 m-p <= 3 时，trigamma 近似会低估先验方差，
      官方改用「模拟 + KL 匹配」（论文 "Three or less residual degrees of
      freedom"）。本模块用 _prior_var_by_simulation 复现这一分支；
      返回值里的 dispPriorVarMethod 记录实际走的分支。
* 收缩（MAP，R/core.R:944-1131）
      在 log(alpha) 上最大化「负二项对数似然 + Cox-Reid 校正 + 正态先验」：
          argmax_t [ logLik(alpha = e^t) + cr(t)
                     - (t - log(dispFit))^2 / (2 * dispPriorVar) ]
      起始值：dispGeneEst（若它低于 dispFit 一个数量级则改用 dispFit）
      结果截断到 [minDisp, max(10, 样本数)]
      WARNING: 同一个上限也必须作用在**基因级 MLE（dispGeneEst）**上，不能只作用在
      MAP 之后：DESeq2 两边都截。漏掉前者会让极端过散布基因（MLE 可达 142）
      污染趋势拟合，并让下游 GLM 的 log2FC 跑到 ±68。真实数据上踩过，见 P8。
* Cox-Reid 校正（src/DESeq2.cpp fitDisp，useCR=TRUE 是官方默认）
      目标函数里额外加一项 cr = -0.5 * log(det(X' W X))，
      w_j = 1/(1/mu_j + alpha)。它随 alpha 单调上升，会把估计系统性推高；
      缺少它会让 dispGeneEst/dispMAP 系统性偏低（对拍时中位比约 0.5/0.75）。
* mu 的来源（R/core.R:736-763、858、1014）
      两组设计下 modelMatrixGroups 的水平数 == 系数个数，官方走
      linearMu=TRUE 分支：mu 是**归一化计数在设计矩阵上的 OLS 拟合**
      再乘回 size factor（linearModelMuNormalized），下限 minmu=0.5；
      不是用 NB GLM 迭代出来的 mu。gene-wise 与 MAP 两步共用这份 mu。
* noIncrease（R/core.R:828-831，niter=1 时）
      若优化后的目标函数值没有超过起点 |lp_init|/1e6 以上，
      dispGeneEst 回退为起点 alpha_init（防数值噪声推动平坦似然的基因）。
* 离群点（R/core.R:1112-1116）
      dispOutlier = log(dispGeneEst) > log(dispFit) + 2 * sqrt(varLogDispEsts)
      被标为离群点的基因，最终离散度直接用 dispGeneEst（不收缩）
"""

import math
import random

from .linalg import build_design_matrix, determinant, mat_vec, ols_hat_matrix
from .nbmath import (MIN_DISP, MIN_DISP_SEARCH, mad, median,
                     sample_variance, trigamma)
from .normalize import base_mean as compute_base_mean
from .normalize import normalized_counts

# 与 DESeq2 一致的常量
OUTLIER_SD = 2.0
MIN_PRIOR_VAR = 0.25
MINMU = 0.5                              # 离散度估计里 mu 的下限（core.R:663）
USE_FOR_FIT_THRESHOLD = MIN_DISP * 100.0   # dispGeneEst >= 1e-6 才参与拟合
# 离散度搜索区间（对数尺度）
LOG_ALPHA_LO = math.log(1e-9)
LOG_ALPHA_HI = math.log(1e6)


# --------------------------------------------------------------------------
# 一维最大化：在 log(alpha) 上找目标函数的最大值
# --------------------------------------------------------------------------
def _golden_section(objective, lo, hi, iterations=28, rel_tol=1e-7):
    """黄金分割搜索（求最大值）。假设 [lo, hi] 内单峰。"""
    invphi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc = objective(c)
    fd = objective(d)
    for _ in range(iterations):
        if (b - a) <= rel_tol * (abs(a) + abs(b) + 1e-300):
            break
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = objective(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = objective(d)
    return c if fc > fd else d


def _maximize_over_log_alpha(objective, grid_points=25):
    """先用粗网格定位，再黄金分割细化。返回最优 t = log(alpha)。

    比单用黄金分割更稳：目标函数在极端低表达基因上可能出现平台或双峰，
    粗网格能先把搜索限制在正确的峰上。
    """
    best_t = LOG_ALPHA_LO
    best_v = float("-inf")
    step = (LOG_ALPHA_HI - LOG_ALPHA_LO) / grid_points
    prev_t = LOG_ALPHA_LO
    prev_v = objective(LOG_ALPHA_LO)
    if prev_v > best_v:
        best_v, best_t = prev_v, prev_t
    for i in range(1, grid_points + 1):
        t = LOG_ALPHA_LO + step * i
        v = objective(t)
        if v > best_v:
            best_v, best_t = v, t
        if v < prev_v:
            # 已越过峰顶，在 [前一点, 当前点] 内细化
            lo = max(LOG_ALPHA_LO, prev_t - step)
            hi = min(LOG_ALPHA_HI, t + step)
            refined = _golden_section(objective, lo, hi)
            vr = objective(refined)
            return refined if vr > best_v else best_t
        prev_t, prev_v = t, v
    # 一直没下降：峰在右端
    return best_t


# --------------------------------------------------------------------------
# 2a-1 矩估计与粗略估计（用于给最大似然定起点）
# --------------------------------------------------------------------------
def moments_disp_estimate(norm_counts, base_means, base_vars, factors):
    """矩估计：(baseVar - xim * baseMean) / baseMean^2。

    对应 DESeq2 的 momentsDispEstimate()。
    """
    xim = sum(1.0 / f for f in factors) / len(factors)
    out = []
    for i in range(len(norm_counts)):
        bm = base_means[i]
        bv = base_vars[i]
        if bm <= 0.0:
            out.append(float("nan"))
        else:
            out.append((bv - xim * bm) / (bm * bm))
    return out


def rough_disp_estimate(norm_counts, X, hat, n_coef):
    """粗略估计：rowSums(((y - mu)^2 - mu) / mu^2) / (m - p)，下限 0。

    对应 DESeq2 的 roughDispEstimate()（R/core.R:2442）。
    mu 由归一化计数在设计矩阵上的最小二乘拟合给出，并下限截断到 1。
    """
    m = len(X)
    p = n_coef
    out = []
    for row in norm_counts:
        beta = mat_vec(hat, row)
        mu = mat_vec(X, beta)
        s = 0.0
        for j in range(m):
            muj = mu[j]
            if muj < 1.0:
                muj = 1.0
            s += ((row[j] - muj) ** 2 - muj) / (muj * muj)
        e = s / (m - p)
        out.append(e if e > 0.0 else 0.0)
    return out


# --------------------------------------------------------------------------
# Cox-Reid 校正项（DESeq2 src/DESeq2.cpp fitDisp）
# --------------------------------------------------------------------------
def _disp_loglik(counts, mu, alpha):
    """与官方 DESeq2.cpp log_posterior 的 ll_part 逐字一致的形式：

        sum( lgamma(y + 1/α) - lgamma(1/α)
             - y*log(mu + 1/α) - (1/α)*log(1 + mu*α) )

    它省略了 -lgamma(y+1) 与 +y*log(mu) 两个「对 alpha 是常数」的项。
    优化结果不受影响，但**绝对值不同**：官方 noIncrease 判定用的是
    |initial_lp|/1e6 这种与绝对值挂钩的阈值，所以这里必须用官方的形式，
    否则回退判定的灵敏度会差出几个数量级。
    """
    size = 1.0 / alpha
    total = 0.0
    lgamma = math.lgamma
    log = math.log
    for j in range(len(counts)):
        y = counts[j]
        m = mu[j]
        total += (lgamma(y + size) - lgamma(size)
                  - y * log(m + size) - size * math.log1p(m * alpha))
    return total


def _cox_reid_term(mu, alpha, X):
    """cr = -0.5 * log(det(X' W X))，其中 W = diag(mu/(1+alpha*mu))。

    当 det <= 0（数值异常）时返回 0.0，即跳过校正。
    对应官方 useCR=TRUE（默认）的行为。
    """
    m = len(mu)
    p = len(X[0])
    w = [mu[j] / (1.0 + alpha * mu[j]) for j in range(m)]
    XtWX = [[sum(X[j][a] * w[j] * X[j][b] for j in range(m))
             for b in range(p)] for a in range(p)]
    det = determinant(XtWX)
    if det <= 0.0:
        return 0.0
    return -0.5 * math.log(det)


# --------------------------------------------------------------------------
# 2a-2 每个基因的最大似然估计
# --------------------------------------------------------------------------
def gene_wise_dispersion(counts, mu, alpha_start, X=None,
                         use_cr=True, initial_lp=None):
    """给定 mu，在 log(alpha) 上最大化负二项对数似然 + Cox-Reid 校正。

    参数
    ----
    X       : 设计矩阵（用于 CR 校正）；若为 None 则 skip CR
    use_cr  : 是否启用 Cox-Reid 校正（官方默认 TRUE）
    initial_lp: 若提供，niter=1 时执行 noIncrease 回退（lp_final < lp_init+|lp_init|/1e6）

    返回
    ----
    alpha（dispGeneEst），最终边界 pmax(minDisp) 且 pmin(maxDisp) 在 estimate_dispersions 里做
    """
    def objective(t):
        alpha = math.exp(t)
        ll = _disp_loglik(counts, mu, alpha)
        if use_cr and X is not None:
            ll += _cox_reid_term(mu, alpha, X)
        return ll

    if alpha_start is None or alpha_start != alpha_start or alpha_start <= 0.0:
        alpha_start = MIN_DISP_SEARCH
    if alpha_start < MIN_DISP_SEARCH:
        alpha_start = MIN_DISP_SEARCH

    t_best = _maximize_over_log_alpha(objective)
    alpha = math.exp(t_best)
    if alpha < MIN_DISP_SEARCH:
        alpha = MIN_DISP_SEARCH

    # noIncrease 回退（core.R:828-831，niter=1 时）
    if initial_lp is not None:
        final_lp = objective(t_best)
        if final_lp < initial_lp + abs(initial_lp) / 1e6:
            alpha = alpha_start

    return alpha


# --------------------------------------------------------------------------
# 2b 趋势拟合
# --------------------------------------------------------------------------
def parametric_dispersion_fit(means, disps, max_iter=10, tol=1e-8):
    """拟合 y = a + b / mean（Gamma GLM，identity 链接）。

    返回 (a, b)。完全对应 DESeq2 的 parametricDispersionFit()：
    Gamma 分布的方差函数是 mu^2，identity 链接下 dmu/deta = 1，
    所以 IRLS 权重 w = 1/mu^2，工作响应 z = eta + (y - mu)。
    """
    a, b = 0.1, 1.0
    for _ in range(max_iter):
        # 剔除残差异常的基因（与官方一致：residuals = disps / (a + b/means)）
        keep_m, keep_d = [], []
        for i in range(len(means)):
            fit = a + b / means[i]
            if fit <= 0.0:
                continue
            r = disps[i] / fit
            if 1e-4 < r < 15.0:
                keep_m.append(means[i])
                keep_d.append(disps[i])
        if len(keep_d) < 3:
            raise ValueError("parametric dispersion fit failed: 可用基因太少")

        a_new, b_new = _gamma_identity_irls(keep_m, keep_d, a, b, tol)
        if not (a_new > 0.0 and b_new > 0.0):
            raise ValueError("parametric dispersion fit failed: 系数非正")
        # 收敛判据与官方一致：sum(log(new/old)^2) < 1e-6
        if a > 0 and b > 0:
            rel = math.log(a_new / a) ** 2 + math.log(b_new / b) ** 2
            a, b = a_new, b_new
            if rel < 1e-6:
                return a, b
        else:
            a, b = a_new, b_new
    return a, b


def _gamma_identity_irls(means, disps, a0, b0, tol=1e-8, maxit=50):
    """Gamma 分布 + identity 链接的加权最小二乘迭代。"""
    m = len(means)
    coef = [a0, b0]
    dev_old = None
    for _ in range(maxit):
        mu = [coef[0] + coef[1] / means[i] for i in range(m)]
        if any(v <= 0.0 for v in mu):
            raise ValueError("Gamma GLM 发散：mu 出现非正值")
        w = [1.0 / (v * v) for v in mu]
        z = [mu[i] + (disps[i] - mu[i]) for i in range(m)]   # dmu/deta = 1

        sw, swx, swy, swxx, swxy = 0.0, 0.0, 0.0, 0.0, 0.0
        for i in range(m):
            x = 1.0 / means[i]
            wi = w[i]
            sw += wi
            swx += wi * x
            swy += wi * z[i]
            swxx += wi * x * x
            swxy += wi * x * z[i]
        det = sw * swxx - swx * swx
        if abs(det) < 1e-300:
            raise ValueError("Gamma GLM 设计矩阵奇异")
        coef = [(swxx * swy - swx * swxy) / det,
                (sw * swxy - swx * swy) / det]

        mu = [coef[0] + coef[1] / means[i] for i in range(m)]
        if any(v <= 0.0 for v in mu):
            raise ValueError("Gamma GLM 发散：mu 出现非正值")
        dev = 0.0
        for i in range(m):
            dev += 2.0 * ((disps[i] - mu[i]) / mu[i] - math.log(disps[i] / mu[i]))
        if dev_old is not None and abs(dev - dev_old) / (abs(dev) + 0.1) < tol:
            break
        dev_old = dev
    return coef[0], coef[1]


def local_dispersion_fit(means, disps, n_bins=20):
    """均匀分箱 + 中位数的「局部」趋势，作为 parametric 拟合失败时的降级方案。

    说明：DESeq2 在 parametric 失败时改用 locfit 做局部回归。本项目不能引入
    第三方库，因此用一个可解释的替代方案：按 log(baseMean) 分箱取中位数，
    再在 log-log 尺度上线性插值。REPORT.md 里会写明这一处差异。
    """
    pairs = sorted(zip(means, disps), key=lambda p: p[0])
    n = len(pairs)
    if n == 0:
        raise ValueError("local_dispersion_fit: 无数据")
    bin_size = max(1, n // n_bins)
    xs, ys = [], []
    for start in range(0, n, bin_size):
        chunk = pairs[start:start + bin_size]
        if not chunk:
            continue
        xs.append(math.log(median([c[0] for c in chunk])))
        ys.append(math.log(median([c[1] for c in chunk])))

    def disp_function(mean):
        if mean <= 0.0:
            return disps[0]
        t = math.log(mean)
        if t <= xs[0]:
            return math.exp(ys[0])
        if t >= xs[-1]:
            return math.exp(ys[-1])
        # 线性插值
        lo = 0
        for i in range(1, len(xs)):
            if t <= xs[i]:
                lo = i - 1
                break
        x0, x1 = xs[lo], xs[lo + 1]
        y0, y1 = ys[lo], ys[lo + 1]
        w = 0.0 if x1 == x0 else (t - x0) / (x1 - x0)
        return math.exp(y0 + w * (y1 - y0))

    return disp_function


# --------------------------------------------------------------------------
# 2c-0 小残差自由度下的先验方差：模拟估计
#      （论文 "Three or less residual degrees of freedom"）
# --------------------------------------------------------------------------
def _prior_var_by_simulation(residuals, m_minus_p, n_sim=10000, grid_points=200,
                             n_bins=80):
    """当 1 <= m-p <= 3 时，用模拟 + KL 匹配估计先验方差 sigma_d^2。

    论文指出：残差自由度很小时，``s_lr^2 - psi_1((m-p)/2)`` 这个近似会
    低估 sigma_d^2。做法是把观测到的对数残差分布与
    「log(chi2_{m-p}) + N(0, sigma_d^2) - log(m-p)」的模拟分布做匹配，
    取 KL 散度最小的 sigma_d^2。

    实现要点（对应官方 R 实现 R/core.R:1158-1192）：
      * 直方图区间 [-10, 10]、80 个 bin（官方 brks = seq(-20, 20, 0.5)）
      * sigma_d^2 网格 [0, 8] 共 200 点（官方 obsVarGrid）
      * 每格模拟 n_sim 个样本；这里用**公共随机数**（同一批 log-chi2 与
        标准正态），比官方每格独立重抽更平滑，且结果确定可复现
      * KL 散度按官方写法：small 取两边正密度最小值，逐 bin 累加
      * 平滑：官方用 loess(span=0.2)，这里用等权滑动平均近似
      * 下限 0.25（与官方一致）

    WARNING: 已知差异：官方用 R 的 set.seed(2) 固定随机序列，我们用
    Python random.Random(2)：序列不同，具体数值不会与官方逐位一致
    （MC 噪声量级），但方法与结果性质一致、可复现。

    返回 None 表示数据不足，调用方应回退到 trigamma 公式。
    """
    lo, hi = -10.0, 10.0
    width = (hi - lo) / n_bins

    def hist_density(vals):
        counts_ = [0] * n_bins
        for v in vals:
            if lo <= v < hi:
                k = int((v - lo) / width)
                if k >= n_bins:
                    k = n_bins - 1
                counts_[k] += 1
        tot = sum(counts_)
        if tot == 0:
            return None
        return [c / (tot * width) for c in counts_]

    obs = [r for r in residuals if lo <= r < hi]
    if len(obs) < 20:
        return None
    obs_dens = hist_density(obs)
    if obs_dens is None:
        return None

    rng = random.Random(2)
    log_chi2 = [math.log(rng.gammavariate(m_minus_p / 2.0, 2.0))
                for _ in range(n_sim)]
    z = [rng.gauss(0.0, 1.0) for _ in range(n_sim)]
    shift = math.log(m_minus_p)

    kl = []
    for gi in range(grid_points):
        sigma2 = 8.0 * gi / (grid_points - 1)
        s = math.sqrt(sigma2)
        sim = [log_chi2[i] + s * z[i] - shift for i in range(n_sim)]
        dens = hist_density(sim)
        if dens is None:
            kl.append(float("inf"))
            continue
        both = [d for d in obs_dens + dens if d > 0.0]
        small = min(both) if both else 0.0
        total = 0.0
        for k in range(n_bins):
            if obs_dens[k] > 0.0:
                total += obs_dens[k] * (math.log(obs_dens[k] + small)
                                        - math.log(dens[k] + small))
        kl.append(total)

    # 滑动平均平滑（近似官方的 loess, span = 0.2）
    half = max(1, int(0.2 * grid_points) // 2)
    smoothed = []
    for i in range(grid_points):
        a = max(0, i - half)
        b = min(grid_points, i + half + 1)
        chunk = [v for v in kl[a:b] if v == v and v != float("inf")]
        smoothed.append(sum(chunk) / len(chunk) if chunk else float("inf"))

    best_i = min(range(grid_points), key=lambda i: smoothed[i])
    if smoothed[best_i] == float("inf"):
        return None
    best = 8.0 * best_i / (grid_points - 1)
    return best if best > MIN_PRIOR_VAR else MIN_PRIOR_VAR


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def estimate_dispersions(counts, factors, groups, quiet=True, progress=None):
    """M2 全流程：gene-wise -> trend -> MAP 收缩。

    参数
    ----
    counts  : counts[基因][样本]，原始计数
    factors : 每个样本的 size factor（M1 的输出）
    groups  : 每个样本的 0/1 分组（1 = 处理组）
    progress: 可选回调 progress(done, total)

    返回
    ----
    dict，包含 baseMean / baseVar / allZero / dispGeneEst / dispFit /
    dispMAP / dispersion / dispOutlier / dispPriorVar / fitType
    """
    n_genes = len(counts)
    n_samples = len(factors)
    X = build_design_matrix(groups)
    n_coef = len(X[0])
    m_minus_p = n_samples - n_coef

    # ---- 归一化与基本统计量 ----
    norm = normalized_counts(counts, factors)
    base_means = compute_base_mean(norm)
    base_vars = [sample_variance(row) for row in norm]
    all_zero = [base_means[i] <= 0.0 for i in range(n_genes)]

    # DESeq2 对离散度的取值范围是 [minDisp, max(10, 样本数)]（core.R:728-729、
    # 786、849、1101-1102）。这个范围在起点、gene-wise MLE、MAP 三处都生效。
    # WARNING: 上限必须在**基因级估计这一步**就截断，不能只在 MAP 收缩那一步截。
    # 真实数据暴露的后果（2026-09-26，airway）：像 ENSG00000229807 这种
    # 「同一组内一个样本 3929、另一个样本 2」的极端过散布基因，MLE 会跑到
    # 142（官方被截在 10）；这个虚高的值随后有三个坏影响：
    #   ① 进入趋势拟合，把 dispFit 整体抬高约 1.4 倍；
    #   ② 让该基因 GLM 的权重 w = mu/(1+alpha*mu) 退化成常数，
    #      log2FC 一路跑到 -68（官方 -0.405）；
    #   ③ 使 MAP 的起点与先验都失真。
    max_disp = 10.0 if 10.0 > n_samples else float(n_samples)

    # ---- 2a：起点（矩估计 与 粗略估计取小），并截断到 [minDisp, maxDisp] ----
    moments = moments_disp_estimate(norm, base_means, base_vars, factors)
    try:
        hat = ols_hat_matrix(X)
        rough = rough_disp_estimate(norm, X, hat, n_coef)
    except ValueError:
        # 设计矩阵病态时退化为只用矩估计
        rough = [float("inf")] * n_genes

    alpha_start = []
    for i in range(n_genes):
        r = rough[i]
        mo = moments[i]
        if r != r:
            r = float("inf")
        if mo != mo or mo <= 0.0:
            mo = float("inf")
        v = r if r < mo else mo
        if v == float("inf"):
            v = MIN_DISP
        if v < MIN_DISP:
            v = MIN_DISP
        elif v > max_disp:
            v = max_disp
        alpha_start.append(v)

    # ---- 2a：mu 的估计（与 DESeq2 的 linearMu=TRUE 分支一致）----
    # DESeq2 对「分组数 == 系数个数」的设计（本项目正是）走 linearMu：
    # 用归一化计数在设计矩阵上做 OLS，得到归一化 mu，再乘回 size factor
    # 还原到原始计数尺度，下限 minmu=0.5。与 NB GLM 迭代出的 mu 不同。
    try:
        hat = ols_hat_matrix(X)
        mu_linear = []
        for i in range(n_genes):
            if all_zero[i]:
                mu_linear.append(None)
                continue
            beta_hat = mat_vec(hat, norm[i])
            mu_norm = mat_vec(X, beta_hat)
            mu_linear.append([
                max(mu_norm[j] * factors[j], MINMU) for j in range(n_samples)
            ])
    except ValueError:
        # 设计矩阵病态时退化：用归一化均值直接展开
        mu_linear = []
        for i in range(n_genes):
            if all_zero[i]:
                mu_linear.append(None)
                continue
            mu_linear.append([
                max(base_means[i] * factors[j], MINMU) for j in range(n_samples)
            ])

    # ---- 2a：先用 alpha_start 拟合 GLM 得到 mu，再做最大似然 ----
    # WARNING: DESeq2 对离散度的取值范围是 [minDisp, max(10, 样本数)]，这个上限必须
    # 在**基因级估计这一步**就截断，不能只在 MAP 收缩那一步截。
    # 真实数据暴露的后果（2026-09-26，airway）：像 ENSG00000229807 这种
    # 「同一组内一个样本 3929、另一个样本 2」的极端过散布基因，MLE 会跑到
    # 142（官方被截在 10）；这个虚高的值随后有三个坏影响：
    #   ① 进入趋势拟合，把 dispFit 整体抬高约 1.4 倍；
    #   ② 让该基因 GLM 的权重 w = mu/(1+alpha*mu) 退化成常数，
    #      log2FC 一路跑到 -68（官方 -0.405）；
    #   ③ 使 MAP 的起点与先验都失真。
    max_disp = 10.0 if 10.0 > n_samples else float(n_samples)
    disp_gene = [float("nan")] * n_genes
    mu_cache = [None] * n_genes
    for i in range(n_genes):
        if all_zero[i]:
            continue
        # DESeq2 linearMu=TRUE 时，mu 直接来自 OLS，不用 fit_nb_glm
        mu = mu_linear[i]
        mu_cache[i] = mu
        # 初始对数似然（用于 noIncrease 判定；必须用与官方 ll_part
        # 一致的形式，否则 |initial_lp|/1e6 这个阈值的灵敏度对不上）
        init_alpha = max(alpha_start[i], MIN_DISP_SEARCH)
        initial_ll = _disp_loglik(counts[i], mu, init_alpha)
        if X is not None:
            initial_ll += _cox_reid_term(mu, init_alpha, X)
        v = gene_wise_dispersion(counts[i], mu, alpha_start[i],
                                 X=X, use_cr=True, initial_lp=initial_ll)
        # 官方最终边界：pmin(pmax(dispGeneEst, minDisp), maxDisp)
        v = max(v, MIN_DISP)
        v = min(v, max_disp)
        disp_gene[i] = v
        if progress is not None and (i + 1) % 2000 == 0:
            progress(i + 1, n_genes)

    # ---- 2b：趋势拟合 ----
    fit_means, fit_disps = [], []
    for i in range(n_genes):
        if all_zero[i]:
            continue
        d = disp_gene[i]
        if d == d and d > USE_FOR_FIT_THRESHOLD:
            fit_means.append(base_means[i])
            fit_disps.append(d)

    fit_type = "parametric"
    if len(fit_disps) >= 3:
        try:
            a, b = parametric_dispersion_fit(fit_means, fit_disps)
            disp_function = lambda q: a + b / q
        except ValueError:
            fit_type = "local"
            disp_function = local_dispersion_fit(fit_means, fit_disps)
    elif len(fit_disps) > 0:
        fit_type = "mean"
        mean_disp = sum(fit_disps) / len(fit_disps)
        disp_function = lambda q: mean_disp
    else:
        raise ValueError(
            "所有基因的离散度估计都低于 1e-6，趋势拟合无法进行"
            "（与 DESeq2 的同名报错一致）")

    disp_fit = [float("nan")] * n_genes
    for i in range(n_genes):
        if not all_zero[i]:
            v = disp_function(base_means[i])
            if v < MIN_DISP:
                v = MIN_DISP
            disp_fit[i] = v

    # ---- 2c：先验方差 ----
    residuals = []
    for i in range(n_genes):
        if all_zero[i]:
            continue
        d = disp_gene[i]
        if d == d and d >= USE_FOR_FIT_THRESHOLD:
            residuals.append(math.log(d) - math.log(disp_fit[i]))
    if not residuals:
        raise ValueError("没有离散度估计高于 minDisp，无法估计先验方差")

    var_log_disp = mad(residuals) ** 2
    prior_method = "trigamma"
    if 1 <= m_minus_p <= 3:
        # 论文 "Three or less residual degrees of freedom"：小自由度下
        # trigamma 近似会低估 sigma_d^2，改用模拟 + KL 匹配
        sim_var = _prior_var_by_simulation(residuals, m_minus_p)
        if sim_var is not None:
            disp_prior_var = sim_var
            prior_method = "simulation"
        else:
            exp_var = trigamma(m_minus_p / 2.0)
            disp_prior_var = var_log_disp - exp_var
            if disp_prior_var < MIN_PRIOR_VAR:
                disp_prior_var = MIN_PRIOR_VAR
            prior_method = "trigamma-fallback"
    elif m_minus_p > 0:
        exp_var = trigamma(m_minus_p / 2.0)
        disp_prior_var = var_log_disp - exp_var
        if disp_prior_var < MIN_PRIOR_VAR:
            disp_prior_var = MIN_PRIOR_VAR
    else:
        disp_prior_var = var_log_disp
        prior_method = "noSubtraction"

    # ---- 2c：MAP 收缩 ----
    sqrt_var = math.sqrt(var_log_disp)
    disp_map = [float("nan")] * n_genes
    disp_outlier = [False] * n_genes
    dispersion = [float("nan")] * n_genes

    for i in range(n_genes):
        if all_zero[i]:
            continue
        dg = disp_gene[i]
        dfit = disp_fit[i]

        # 离群点：高于趋势线 2 倍标准差的基因保留 gene-wise 值（不收缩）
        if dg == dg and math.log(dg) > math.log(dfit) + OUTLIER_SD * sqrt_var:
            disp_outlier[i] = True
            disp_map[i] = dg
            dispersion[i] = dg
            continue

        # 起始值：低于趋势线一个数量级时改用趋势值
        init = dg if dg > 0.1 * dfit else dfit
        if init <= 0.0 or init != init:
            init = dfit
        log_prior_mean = math.log(dfit)

        def objective(t, log_prior_mean=log_prior_mean, mu=mu_cache[i]):
            alpha = math.exp(t)
            return (_disp_loglik(counts[i], mu, alpha)
                    + _cox_reid_term(mu, alpha, X)
                    - (t - log_prior_mean) ** 2 / (2.0 * disp_prior_var))

        # 先按起始值附近做一次粗网格，再细化
        t_best = _maximize_over_log_alpha(objective, grid_points=20)
        value = math.exp(t_best)
        if value < MIN_DISP:
            value = MIN_DISP
        elif value > max_disp:
            value = max_disp
        disp_map[i] = value
        dispersion[i] = value

    return {
        "baseMean": base_means,
        "baseVar": base_vars,
        "allZero": all_zero,
        "dispGeneEst": disp_gene,
        "dispFit": disp_fit,
        "dispMAP": disp_map,
        "dispersion": dispersion,
        "dispOutlier": disp_outlier,
        "dispPriorVar": disp_prior_var,
        "dispPriorVarMethod": prior_method,
        "varLogDispEsts": var_log_disp,
        "fitType": fit_type,
    }
