"""方案 C：模拟数据生成器（带已知真值）。

为什么必须有它
--------------
对真实数据你并不知道正确答案，所以没法判断 M2/M3 到底是「算对了」还是「碰巧接近」。
只有自己造一批「真值已知」的数据，才能回答下面这些关键问题：

    M1：我把 size factor 恢复出来了吗？（真值就是生成时用的 sf）
    M2：我把 alpha 估回去了吗？（真值就是生成时用的 alpha）
    M3：我把 log2FC 估回去了吗？塞进去的差异基因被检出了吗？

生成模型
--------
    对基因 g、样本 j：
        mu[g][j]   = mean_g * sizeFactor_j * 2 ** (log2FC_g * group_j)
        counts[g][j] ~ NegativeBinomial(mean = mu[g][j], Var = mu + alpha_g * mu^2)

    其中 alpha_g 按 log 空间正态分布从一条「趋势线」上采样，
    模拟真实数据里「低表达基因离散度高、高表达基因离散度低」的规律。

用的是标准库 random，且固定随机种子，所以结果完全可复现。
"""

import math
import random

from .nbmath import MIN_DISP


def trend_dispersion(base_mean, asympt=0.01, extra=5.0):
    """离散度的趋势线：alpha = asympt + extra / baseMean。

    与 DESeq2 的 parametricDispersionFit 用同一个函数族，
    这样模拟数据里就天然含有「离散度随表达量下降」的结构。
    """
    return asympt + extra / max(base_mean, 1e-8)


def _sample_nb(rand, mu, alpha):
    """从负二项分布抽一个计数。

    用 Gamma-Poisson 混合：先抽 lambda ~ Gamma(shape=1/alpha, scale=alpha*mu)，
    再抽 y ~ Poisson(lambda)。这是负二项的精确等价定义。
    """
    if alpha <= 0.0:
        # 退化为泊松
        return _sample_poisson(rand, mu)
    shape = 1.0 / alpha
    # Gamma(shape, scale) 用 Marsaglia-Tsang 方法
    lam = _sample_gamma(rand, shape, alpha * mu)
    if lam < 0.0:
        lam = 0.0
    return _sample_poisson(rand, lam)


def _sample_gamma(rand, shape, scale):
    """Marsaglia & Tsang (2000) 的 Gamma 采样（shape >= 1）。"""
    if shape < 1.0:
        # 用 boost：Gamma(a) = Gamma(a+1) * U^(1/a)
        u = rand.random()
        if u <= 0.0:
            u = 1e-12
        return _sample_gamma(rand, shape + 1.0, scale) * (u ** (1.0 / shape))
    d = shape - 1.0 / 3.0
    c = 1.0 / math.sqrt(9.0 * d)
    while True:
        x = rand.gauss(0.0, 1.0)
        v = (1.0 + c * x) ** 3
        if v <= 0.0:
            continue
        u = rand.random()
        if u == 0.0:
            u = 1e-12
        if math.log(u) < 0.5 * x * x + d - d * v + d * math.log(v):
            return d * v * scale


def _sample_poisson(rand, lam):
    """Knuth 算法；lam 较大时用正态近似（够用且快）。"""
    if lam <= 0.0:
        return 0
    if lam < 30.0:
        limit = math.exp(-lam)
        k = 0
        p = 1.0
        while True:
            k += 1
            p *= rand.random()
            if p <= limit:
                return k - 1
    # 正态近似 + 修正
    v = rand.gauss(lam, math.sqrt(lam))
    if v < 0.0:
        return 0
    return int(v + 0.5)


def simulate_dataset(n_genes=2000, n_per_group=4, seed=20260926,
                     n_de=None, true_lfc=1.5,
                     size_factors=None, base_means=None,
                     asympt=0.01, extra=5.0, disp_scale=0.4,
                     zero_fraction=0.05):
    """生成一份模拟计数数据。

    参数
    ----
    n_genes      : 基因数
    n_per_group  : 每组的样本数（默认 4，与 airway 的 4v4 一致）
    seed         : 随机种子（固定 -> 完全可复现）
    n_de         : 有多少个基因是真的差异表达（默认 n_genes 的 5%）
    true_lfc     : 差异基因的真实 |log2FC|
    size_factors : 各样本真实的 size factor（默认随机取 0.7 ~ 1.4）
    base_means   : 各基因真实的平均表达量（默认对数均匀分布）
    asympt/extra : 离散度趋势线参数
    disp_scale   : 基因级离散度围绕趋势线的波动（对数正态的标准差）
    zero_fraction: 有多少比例的基因被设成「全 0」（用于测试边界处理）

    返回
    ----
    dict，含
        counts        counts[基因][样本]
        group         各样本的 0/1 分组（前 n_per_group 个是组 1）
        true_sf       真实的 size factor
        true_alpha    每个基因真实的离散度
        true_lfc      每个基因真实的 log2FC
        true_base_mean 每个基因真实的基准表达量
        de_mask       哪些基因是真的差异基因
        all_zero_mask 哪些基因被人为设成全 0
    """
    rand = random.Random(seed)
    n_samples = 2 * n_per_group

    if size_factors is None:
        size_factors = [math.exp(rand.gauss(0.0, 0.2)) for _ in range(n_samples)]
    if base_means is None:
        # 对数均匀：模拟真实数据里「少数高表达 + 大量低表达」的分布
        base_means = [math.exp(rand.uniform(math.log(1.0), math.log(5000.0)))
                      for _ in range(n_genes)]
    if n_de is None:
        n_de = max(1, int(n_genes * 0.05))

    de_mask = [False] * n_genes
    de_indices = rand.sample(range(n_genes), min(n_de, n_genes))
    for i in de_indices:
        de_mask[i] = True

    lfc_out = []
    for i in range(n_genes):
        if de_mask[i]:
            # 上/下调各半
            sign = 1.0 if rand.random() < 0.5 else -1.0
            lfc_out.append(sign * true_lfc * rand.uniform(0.8, 1.2))
        else:
            lfc_out.append(0.0)

    true_alpha = []
    for i in range(n_genes):
        base = trend_dispersion(base_means[i], asympt, extra)
        a = base * math.exp(rand.gauss(0.0, disp_scale))
        if a < MIN_DISP:
            a = MIN_DISP
        true_alpha.append(a)

    # 全 0 基因（用于验证边界处理）
    all_zero_mask = [False] * n_genes
    n_zero = int(n_genes * zero_fraction)
    for i in rand.sample(range(n_genes), n_zero):
        all_zero_mask[i] = True

    group = [1] * n_per_group + [0] * n_per_group

    counts = []
    for i in range(n_genes):
        if all_zero_mask[i]:
            counts.append([0] * n_samples)
            continue
        row = []
        for j in range(n_samples):
            mu = base_means[i] * size_factors[j] * \
                (2.0 ** (lfc_out[i] * group[j]))
            if mu < 1e-6:
                mu = 1e-6
            row.append(_sample_nb(rand, mu, true_alpha[i]))
        counts.append(row)

    return {
        "counts": counts,
        "group": group,
        "true_sf": size_factors,
        "true_alpha": true_alpha,
        "true_lfc": lfc_out,
        "true_base_mean": base_means,
        "de_mask": de_mask,
        "all_zero_mask": all_zero_mask,
        "n_per_group": n_per_group,
    }


def write_simulated_dataset(dataset, counts_path, design_path,
                            sample_prefix="SIM"):
    """把模拟数据写成与真实数据完全相同的 CSV 格式。"""
    from .io_utils import write_generic_csv

    n_samples = len(dataset["group"])
    sample_ids = ["%s%02d" % (sample_prefix, j + 1) for j in range(n_samples)]
    gene_ids = ["GENE%05d" % (i + 1) for i in range(len(dataset["counts"]))]

    write_generic_csv(counts_path,
                      [""] + sample_ids,
                      [[gid] + row for gid, row in zip(gene_ids, dataset["counts"])])

    design_rows = []
    for j, sid in enumerate(sample_ids):
        g = dataset["group"][j]
        design_rows.append([sid, "trt" if g == 1 else "untrt",
                            "cell%d" % (1 if g == 1 else 2), g])
    write_generic_csv(design_path, ["sample", "condition", "cell", "group"],
                      design_rows)
    return gene_ids, sample_ids
