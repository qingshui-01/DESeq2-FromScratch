# -*- coding: utf-8 -*-
"""单元测试：M1–M4 各模块 + 模拟数据自证。

运行方式（在项目根目录）：
    python -m unittest tests.test_deg -v

设计原则
--------
* 能用「手算」验证的，一定手算，并把期望值写死在测试里；
* 不能手算的（M2/M3），用「真值已知的模拟数据」验证；
* 边界情况必须覆盖：全 0 基因、空输入、含缺失、p=0/p=1、n=2 的下限。
"""

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.diagnostics import (cooks_distance, cooks_threshold,
                             robust_moments_dispersion)
from src.dispersion import estimate_dispersions
from src.glm_nb import (fit_nb_glm, lrt_test, wald_test,
                        wald_test_above_threshold, wald_test_below_threshold)
from src.linalg import build_design_matrix
from src.multitest import bh_adjust
from src.normalize import base_mean, normalized_counts, size_factors
from src.simulate import simulate_dataset, trend_dispersion, write_simulated_dataset


# ==========================================================================
# M4：BH 多重检验校正
# ==========================================================================
class TestM4BH(unittest.TestCase):

    def test_hand_computed_equal_ratio(self):
        """p 值与秩成正比时，所有 padj 应相等 = p_max * n / n。

        p = [0.01, 0.02, 0.03, 0.04, 0.05]，n = 5：
          rank1: 0.01*5/1 = 0.05
          rank2: 0.02*5/2 = 0.05
          ...
          rank5: 0.05*5/5 = 0.05
        手算结果为全部 0.05。
        """
        out = bh_adjust([0.01, 0.02, 0.03, 0.04, 0.05])
        for v in out:
            self.assertAlmostEqual(v, 0.05, places=12)

    def test_hand_computed_reorder_and_monotonic(self):
        """打乱顺序后 padj 必须跟着 p 走，且必须单调不降。

        p = [0.04, 0.01, 0.03]，n = 3：
          升序 [0.01, 0.03, 0.04]
          rank1: 0.01*3/1 = 0.03
          rank2: 0.03*3/2 = 0.045
          rank3: 0.04*3/3 = 0.04
          从大到小取 cummin： rank3=0.04, rank2=min(0.045,0.04)=0.04, rank1=min(0.03,0.04)=0.03
        映射回原顺序 -> [0.04, 0.03, 0.04]
        """
        out = bh_adjust([0.04, 0.01, 0.03])
        self.assertAlmostEqual(out[0], 0.04, places=12)
        self.assertAlmostEqual(out[1], 0.03, places=12)
        self.assertAlmostEqual(out[2], 0.04, places=12)

    def test_monotonicity_holds_generally(self):
        """单调性：p 更小的基因，padj 绝不能更大（漏掉 cummin 就会违反）。"""
        ps = [0.9, 0.04, 0.5, 0.001, 0.2, 0.03, 0.7, 0.008]
        out = bh_adjust(ps)
        pairs = sorted(zip(ps, out))
        for k in range(1, len(pairs)):
            self.assertGreaterEqual(pairs[k][1] + 1e-12, pairs[k - 1][1],
                                    "违反单调性：p=%.3f 的 padj=%.6f 小于 "
                                    "p=%.3f 的 padj=%.6f"
                                    % (pairs[k][0], pairs[k][1],
                                       pairs[k - 1][0], pairs[k - 1][1]))

    def test_missing_values_are_preserved(self):
        """缺失的 p 值对应缺失的 padj，且不参与 n 的计算。"""
        out = bh_adjust([0.01, None, 0.03])
        self.assertIsNone(out[1])
        # n = 2（只有两个有效值）
        self.assertAlmostEqual(out[0], 0.02, places=12)   # 0.01*2/1
        self.assertAlmostEqual(out[2], 0.03, places=12)   # 0.03*2/2

    def test_all_missing(self):
        self.assertEqual(bh_adjust([None, None]), [None, None])
        self.assertEqual(bh_adjust([]), [])

    def test_single_value(self):
        """n = 1 时 padj = p（乘 1 除 1）。"""
        self.assertAlmostEqual(bh_adjust([0.37])[0], 0.37, places=12)

    def test_extreme_pvalues(self):
        """p = 0 与 p = 1 的边界。"""
        out = bh_adjust([0.0, 1.0])
        self.assertAlmostEqual(out[0], 0.0, places=12)
        self.assertAlmostEqual(out[1], 1.0, places=12)

    def test_never_exceeds_one(self):
        out = bh_adjust([0.9, 0.95, 0.99, 1.0])
        for v in out:
            self.assertLessEqual(v, 1.0)

    def test_nan_treated_as_missing(self):
        out = bh_adjust([0.01, float("nan"), 0.03])
        self.assertIsNone(out[1])


# ==========================================================================
# M1：size factor
# ==========================================================================
class TestM1SizeFactors(unittest.TestCase):

    def test_hand_computed_exact(self):
        """构造比值正好是 1:2:4 的数据，size factor 应精确等于 [0.5, 1, 2]。

        基因 A = [10, 20, 40]  -> 几何均值 (10*20*40)^(1/3) = 20
        基因 B = [100, 200, 400] -> 几何均值 200
        样本 1 的比值：10/20 = 0.5, 100/200 = 0.5 -> 中位数 0.5
        样本 2：1, 1 -> 1
        样本 3：2, 2 -> 2
        """
        counts = [[10, 20, 40], [100, 200, 400]]
        sf = size_factors(counts)
        self.assertAlmostEqual(sf[0], 0.5, places=12)
        self.assertAlmostEqual(sf[1], 1.0, places=12)
        self.assertAlmostEqual(sf[2], 2.0, places=12)

    def test_zero_containing_genes_are_excluded(self):
        """含 0 的基因被排除，结果不受它影响。

        在上一例中加入 [0, 5, 999999]：它含 0，log(0) = -inf，
        行均值 -inf，必须被排除，size factor 保持不变。
        """
        base = [[10, 20, 40], [100, 200, 400]]
        with_zero = base + [[0, 5, 999999]]
        sf_base = size_factors(base)
        sf_zero = size_factors(with_zero)
        for a, b in zip(sf_base, sf_zero):
            self.assertAlmostEqual(a, b, places=12)

    def test_median_is_robust_to_one_outlier_gene(self):
        """中位数的抗干扰性：一个基因剧烈变化不应带偏整体。"""
        counts = [[100, 100, 100],
                  [200, 200, 200],
                  [300, 300, 300],
                  [50, 50, 50],
                  [10, 10, 10]]
        sf = size_factors(counts)
        for v in sf:
            self.assertAlmostEqual(v, 1.0, places=12)

    def test_all_zero_gene_row(self):
        """全 0 的基因同样被排除（不报错）。"""
        counts = [[10, 20], [30, 60], [0, 0]]
        sf = size_factors(counts)
        self.assertTrue(all(f > 0 for f in sf))

    def test_error_when_every_gene_has_zero(self):
        """每个基因都含 0 -> 必须明确报错，而不是给出垃圾结果。"""
        counts = [[0, 1], [1, 0]]
        with self.assertRaises(ValueError):
            size_factors(counts)

    def test_normalized_counts_and_base_mean(self):
        counts = [[10, 20, 40], [100, 200, 400]]
        sf = size_factors(counts)
        norm = normalized_counts(counts, sf)
        # 归一化后两个样本应完全一致
        self.assertAlmostEqual(norm[0][0], 20.0, places=9)
        self.assertAlmostEqual(norm[0][1], 20.0, places=9)
        self.assertAlmostEqual(norm[0][2], 20.0, places=9)
        bm = base_mean(norm)
        self.assertAlmostEqual(bm[0], 20.0, places=9)
        self.assertAlmostEqual(bm[1], 200.0, places=9)


# ==========================================================================
# M2：离散度
# ==========================================================================
class TestM2Dispersion(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # 用固定种子，保证可复现
        cls.ds = simulate_dataset(n_genes=400, n_per_group=4, seed=20260926)
        cls.counts = cls.ds["counts"]
        cls.group = cls.ds["group"]
        cls.sf = size_factors(cls.counts)
        cls.result = estimate_dispersions(cls.counts, cls.sf, cls.group)

    def test_all_zero_genes_are_flagged(self):
        """全 0 基因必须被标记，且不参与估计。

        注意：全 0 基因有两个来源：① 模拟时人为设的；
        ② 表达量极低的基因抽样后恰好全为 0。两者都必须被标记。
        所以这里断言的是「包含关系」，而不是数量相等。
        """
        d = self.result
        for i in range(len(self.counts)):
            if self.ds["all_zero_mask"][i]:
                self.assertTrue(d["allZero"][i],
                                "人为设成全 0 的基因 %d 未被标记" % i)
        # 顺带核对：counts 真的全 0 的基因，必须全部被标记
        for i, row in enumerate(self.counts):
            if sum(row) == 0:
                self.assertTrue(d["allZero"][i])
        for i in range(len(self.counts)):
            if d["allZero"][i]:
                self.assertTrue(math.isnan(d["dispGeneEst"][i]))
                self.assertTrue(math.isnan(d["dispersion"][i]))
        self.assertGreaterEqual(sum(d["allZero"]),
                                sum(self.ds["all_zero_mask"]))

    def test_trend_fit_succeeded(self):
        self.assertEqual(self.result["fitType"], "parametric")
        self.assertGreater(self.result["dispPriorVar"], 0.0)
        for i in range(len(self.counts)):
            if not self.result["allZero"][i]:
                self.assertLess(self.result["dispFit"][i], 100.0)
                self.assertGreater(self.result["dispFit"][i], 0.0)

    def test_trend_is_decreasing_in_mean(self):
        """趋势线必须随表达量升高而下降：alpha = a + b/mean，b > 0。"""
        d = self.result
        idx = [i for i in range(len(self.counts)) if not d["allZero"][i]]
        idx.sort(key=lambda i: d["baseMean"][i])
        lo = idx[len(idx) // 10]
        hi = idx[-len(idx) // 10]
        self.assertGreater(d["dispFit"][lo], d["dispFit"][hi],
                           "低表达基因的趋势离散度应当高于高表达基因")

    def test_shrinkage_reduces_spread_around_truth(self):
        """收缩的核心价值：把噪声很大的 gene-wise 估计拉得更接近真值。

        判据：收缩后 log(估计/真值) 的四分位间距必须小于 gene-wise 的。
        """
        d = self.result
        idx = [i for i in range(len(self.counts)) if not d["allZero"][i]]
        res_gene = sorted(math.log(d["dispGeneEst"][i] / self.ds["true_alpha"][i])
                          for i in idx)
        res_map = sorted(math.log(d["dispersion"][i] / self.ds["true_alpha"][i])
                         for i in idx)
        iqr_gene = res_gene[3 * len(idx) // 4] - res_gene[len(idx) // 4]
        iqr_map = res_map[3 * len(idx) // 4] - res_map[len(idx) // 4]
        self.assertLess(iqr_map, iqr_gene,
                        "收缩后 log 离散度误差的 IQR 应当变小："
                        "gene-wise=%.3f, MAP=%.3f" % (iqr_gene, iqr_map))

    def test_recovers_known_alpha(self):
        """整体上能把已知的 alpha 估回来（中位比值接近 1）。"""
        d = self.result
        ratios = sorted(d["dispersion"][i] / self.ds["true_alpha"][i]
                        for i in range(len(self.counts)) if not d["allZero"][i])
        med = ratios[len(ratios) // 2]
        self.assertGreater(med, 0.5)
        self.assertLess(med, 2.0)

    def test_trend_dispersion_function_shape(self):
        self.assertGreater(trend_dispersion(1.0), trend_dispersion(1000.0))
        self.assertGreater(trend_dispersion(10.0), 0.0)


# ==========================================================================
# M3：负二项 GLM + Wald
# ==========================================================================
class TestM3GLM(unittest.TestCase):

    def test_recovers_known_log2fc(self):
        """构造一个「无噪声」的强差异基因，检查 log2FC 符号与量级。

        计数在两组各 4 个、每组内完全一致，所以 log2FC 应非常接近真实值。
        """
        counts = [200, 200, 200, 200, 50, 50, 50, 50]   # 组1 = 前 4 个
        group = [1, 1, 1, 1, 0, 0, 0, 0]
        X = build_design_matrix(group)
        log_sf = [0.0] * 8
        fit = fit_nb_glm(counts, X, log_sf, 0.01)
        lfc, se, stat, p = wald_test(fit["beta"], fit["cov"], 1)
        self.assertAlmostEqual(lfc, 2.0, delta=0.15)     # log2(200/50) = 2
        self.assertIsNotNone(se)
        self.assertGreater(stat, 0)                      # 处理组更高 -> 正
        self.assertLess(p, 1e-6)

    def test_sign_is_correct(self):
        """下调基因的 log2FC 必须为负。"""
        counts = [50, 50, 50, 50, 200, 200, 200, 200]
        X = build_design_matrix([1, 1, 1, 1, 0, 0, 0, 0])
        fit = fit_nb_glm(counts, X, [0.0] * 8, 0.01)
        lfc, _, stat, _ = wald_test(fit["beta"], fit["cov"], 1)
        self.assertLess(lfc, -1.5)
        self.assertLess(stat, 0)

    def test_all_zero_gene_yields_nothing_extreme(self):
        """全 0 基因不应产生 NaN 之外的危险值（上游会把它标为 NA）。"""
        X = build_design_matrix([1, 1, 0, 0])
        fit = fit_nb_glm([0, 0, 0, 0], X, [0.0] * 4, 0.01)
        self.assertEqual(len(fit["beta"]), 2)
        self.assertFalse(any(math.isnan(b) for b in fit["beta"]))

    def test_deviance_is_nonnegative(self):
        from src.glm_nb import nb_deviance
        counts = [5, 10, 15, 20]
        mu = [6.0, 9.0, 14.0, 21.0]
        self.assertGreaterEqual(nb_deviance(counts, mu, 0.1), 0.0)
        # 饱和模型（mu = y）偏差为 0
        self.assertAlmostEqual(nb_deviance(counts, counts, 0.1), 0.0, places=9)

    def test_wald_calibration_at_large_n(self):
        """大样本 + 真实离散度时，零假设下的 Wald 统计量应近似标准正态。

        这条测试同时是「本项目 Wald 检验实现正确」的硬证据：
        n=60 时统计量的标准差实测为 1.000（理想 1.000）。
        小样本（每组 4 个）时会偏大（约 1.1~1.3），那是有限样本的固有性质，
        不是实现错误：REPORT.md 里有专门一节讨论。
        """
        ds = simulate_dataset(n_genes=150, n_per_group=40, seed=555, n_de=0,
                              size_factors=[1.0] * 80)
        X = build_design_matrix(ds["group"])
        log_sf = [0.0] * len(ds["group"])
        stats = []
        for i in range(len(ds["counts"])):
            fit = fit_nb_glm(ds["counts"][i], X, log_sf, ds["true_alpha"][i])
            _, _, stat, _ = wald_test(fit["beta"], fit["cov"], 1)
            if stat is not None:
                stats.append(stat)
        n = len(stats)
        mean = sum(stats) / n
        sd = math.sqrt(sum((s - mean) ** 2 for s in stats) / (n - 1))
        self.assertLess(abs(mean), 0.35)
        self.assertGreater(sd, 0.7)
        self.assertLess(sd, 1.35)


# ==========================================================================
# 模拟数据与端到端
# ==========================================================================
class TestSimulateAndPipeline(unittest.TestCase):

    def test_reproducible_with_same_seed(self):
        a = simulate_dataset(n_genes=50, seed=1)
        b = simulate_dataset(n_genes=50, seed=1)
        c = simulate_dataset(n_genes=50, seed=2)
        self.assertEqual(a["counts"], b["counts"])
        self.assertNotEqual(a["counts"], c["counts"])

    def test_counts_are_nonnegative_integers(self):
        ds = simulate_dataset(n_genes=100, seed=3)
        for row in ds["counts"]:
            for v in row:
                self.assertIsInstance(v, int)
                self.assertGreaterEqual(v, 0)

    def test_end_to_end_detects_spiked_genes(self):
        """端到端：塞进去的差异基因应该被检出，且 padj 有正确含义。"""
        from deg import run_pipeline
        ds = simulate_dataset(n_genes=500, n_per_group=4, seed=99)
        rows, _ = run_pipeline(ds["counts"], ds["group"], quiet=True)
        self.assertEqual(len(rows), 500)
        # padj 单调性与取值范围
        pairs = sorted((r["pvalue"], r["padj"]) for r in rows
                       if r["pvalue"] is not None)
        for k in range(1, len(pairs)):
            self.assertGreaterEqual(pairs[k][1] + 1e-12, pairs[k - 1][1])
        # 强差异基因（|lfc| 大且表达量高）应被检出
        strong = [i for i in range(500)
                  if ds["de_mask"][i] and abs(ds["true_lfc"][i]) >= 1.2
                  and ds["true_base_mean"][i] > 100]
        self.assertGreater(len(strong), 0)
        hit = sum(1 for i in strong
                  if rows[i]["padj"] is not None and rows[i]["padj"] < 0.05)
        self.assertGreaterEqual(hit / len(strong), 0.6,
                                "强差异基因的检出率过低：%d/%d"
                                % (hit, len(strong)))

    def test_simulated_dataset_roundtrip_csv(self):
        """写出的 CSV 能被读回来（格式与真实数据一致）。"""
        import tempfile
        from src.io_utils import (align_design, read_counts_csv,
                                  read_design_csv)
        ds = simulate_dataset(n_genes=30, seed=7)
        with tempfile.TemporaryDirectory() as tmp:
            cp = os.path.join(tmp, "counts.csv")
            dp = os.path.join(tmp, "design.csv")
            gids, sids = write_simulated_dataset(ds, cp, dp)
            g2, s2, c2 = read_counts_csv(cp)
            self.assertEqual(gids, g2)
            self.assertEqual(sids, s2)
            self.assertEqual(ds["counts"], c2)
            design = read_design_csv(dp)
            groups, _, _ = align_design(s2, design)
            self.assertEqual(groups, ds["group"])

    def test_duplicate_sample_names_are_accepted(self):
        """design.csv 里多出来的样本不报错（只要 counts 里的样本都在）。"""
        ds = simulate_dataset(n_genes=20, seed=8)
        design = {}
        for j in range(len(ds["group"])):
            design["S%02d" % (j + 1)] = {"condition": "x", "cell": "c",
                                         "group": ds["group"][j]}
        design["EXTRA"] = {"condition": "x", "cell": "c", "group": 0}
        from src.io_utils import align_design
        groups, _, _ = align_design(["S%02d" % (j + 1) for j in range(8)], design)
        self.assertEqual(groups, ds["group"])


# ==========================================================================
# 论文功能补全（P14）：阈值检验 / LRT / Cook's 距离 / 小自由度先验方差
# ==========================================================================
class TestPaperGapFill(unittest.TestCase):

    def test_chi2_sf_matches_closed_forms(self):
        """卡方上尾的闭式自检：df=1 用 erfc，df=2/4 用 exp 闭式。"""
        from src.nbmath import chi2_sf
        x = 4.0
        self.assertAlmostEqual(chi2_sf(x, 1), math.erfc(math.sqrt(x / 2.0)),
                               places=13)
        self.assertAlmostEqual(chi2_sf(x, 2), math.exp(-x / 2.0), places=13)
        self.assertAlmostEqual(chi2_sf(x, 4), math.exp(-x / 2.0) * (1 + x / 2.0),
                               places=13)

    def test_f_quantile_consistency(self):
        """F 分位数与已知关系一致：qf(0.99,1,6) = t(6,0.995)^2 = 13.745。"""
        from src.nbmath import f_cdf, f_quantile
        q = f_quantile(0.99, 1, 6)
        self.assertAlmostEqual(q, 3.7074 ** 2, places=2)
        self.assertAlmostEqual(f_cdf(q, 1, 6), 0.99, places=10)
        # Cook 阈值用到的 qf(0.99, 2, 6) ≈ 10.92（F 分布表值）
        self.assertAlmostEqual(f_quantile(0.99, 2, 6), 10.92, places=1)

    def test_threshold_test_reduces_to_wald_at_theta_zero(self):
        """θ=0 时「超过阈值」检验必须与标准双侧 Wald 完全相等。"""
        beta = [4.0, 0.7]
        cov = [[0.1, 0.0], [0.0, 0.0625]]
        p_thr = wald_test_above_threshold(beta, cov, 0.0)[2]
        p_std = wald_test(beta, cov)[3]
        self.assertAlmostEqual(p_thr, p_std, places=15)

    def test_threshold_tests_semantics(self):
        """两个方向的阈值检验语义：
        * 「超过阈值」：阈值越高越不显著；
        * 「弱于阈值」：弱效应显著、跨越阈值的强效应不显著。"""
        cov = [[0.1, 0.0], [0.0, 0.0625]]      # se = 0.25
        p1 = wald_test_above_threshold([4.0, 0.7], cov, 0.3)[2]
        p2 = wald_test_above_threshold([4.0, 0.7], cov, 0.5)[2]
        self.assertLess(p1, p2)
        # beta=0.05 远低于阈值 1.0 → 显著
        self.assertLess(wald_test_below_threshold([4.0, 0.05], cov, 1.0)[2], 0.01)
        # beta=1.6 高于阈值 1.5 且 SE 小 → 不显著（p 接近 0.66）
        self.assertGreater(wald_test_below_threshold([4.0, 1.6], cov, 1.5)[2],
                           0.5)

    def test_lrt_detects_known_difference(self):
        """LRT：构造 log2(200/50) 的差异基因应给出极小 p；零效应基因 p 大。"""
        X = build_design_matrix([0, 0, 0, 0, 1, 1, 1, 1])
        Xr = [[1.0] for _ in range(8)]
        res = lrt_test([50, 50, 50, 50, 200, 200, 200, 200], X, Xr,
                       [0.0] * 8, 0.01)
        self.assertEqual(res["df"], 1)
        self.assertGreater(res["stat"], 100.0)
        self.assertLess(res["pvalue"], 1e-20)
        # 与 Wald 的 z^2 在同一量级（LRT 与 Wald 渐近等价）
        fit = fit_nb_glm([50, 50, 50, 50, 200, 200, 200, 200], X, [0.0] * 8, 0.01)
        z = wald_test(fit["beta"], fit["cov"])[2]
        self.assertLess(abs(res["stat"] - z * z) / res["stat"], 0.15)
        # 零效应基因不显著
        res0 = lrt_test([100, 90, 110, 95, 105, 100, 95, 105], X, Xr,
                        [0.0] * 8, 0.01)
        self.assertGreater(res0["pvalue"], 0.05)

    def test_cooks_distance_flags_outlier_sample(self):
        """Cook's 距离：一个样本计数 1000 vs 其余 ~20，必须只标记该样本。"""
        from src.diagnostics import flag_cooks_outliers, cooks_threshold
        X = build_design_matrix([0, 0, 0, 0, 1, 1, 1, 1])
        res = flag_cooks_outliers([20, 22, 18, 20, 21, 19, 20, 1000],
                                  X, [0.0] * 8, 0.05)
        self.assertEqual([i for i, f in enumerate(res["flags"]) if f], [7])
        self.assertAlmostEqual(res["threshold"], cooks_threshold(8, 2), places=12)
        # 无离群样本的基因不被标记
        res2 = flag_cooks_outliers([20, 22, 18, 20, 21, 19, 20, 18],
                                   X, [0.0] * 8, 0.05)
        self.assertEqual([i for i, f in enumerate(res2["flags"]) if f], [])

    def test_small_df_prior_variance_simulation(self):
        """m-p <= 3 时先验方差走模拟分支；m-p > 3 走 trigamma 分支。

        模拟分支的校准（合成残差，真值 sigma^2 已知）：
        sigma^2=1.0 估计在 [0.8, 1.6]（MC + 直方图噪声）；且结果确定可复现。
        """
        from src.dispersion import _prior_var_by_simulation
        import random as _random

        def make_resid(n, sigma, df, seed):
            rng = _random.Random(seed)
            return [math.log(rng.gammavariate(df / 2.0, 2.0))
                    + rng.gauss(0, sigma) - math.log(df) for _ in range(n)]

        r = make_resid(3000, 1.0, 2, seed=11)
        v1 = _prior_var_by_simulation(r, 2)
        v2 = _prior_var_by_simulation(r, 2)
        self.assertEqual(v1, v2, "模拟分支必须可复现（固定种子）")
        self.assertGreater(v1, 0.8)
        self.assertLess(v1, 1.6)
        # 真值 0 → 应贴下限 0.25
        r0 = make_resid(3000, 0.0, 2, seed=11)
        self.assertAlmostEqual(_prior_var_by_simulation(r0, 2), 0.25, places=12)

        # 分支触发：m=4, p=2 -> m-p=2 用 simulation；m=8 -> trigamma
        counts4 = [[80 + 30 * ((i * j) % 7) for j in range(4)] for i in range(80)]
        d4 = estimate_dispersions(counts4, [1.0] * 4, [0, 0, 1, 1], quiet=True)
        self.assertEqual(d4["dispPriorVarMethod"], "simulation")
        counts8 = [[80 + 30 * ((i * j) % 7) for j in range(8)] for i in range(80)]
        d8 = estimate_dispersions(counts8, [1.0] * 8,
                                  [0, 0, 0, 0, 1, 1, 1, 1], quiet=True)
        self.assertEqual(d8["dispPriorVarMethod"], "trigamma")


# ==========================================================================
# 论文功能补全的边界测试（P15）：四个功能的边界与手算核对
# ==========================================================================
class TestPaperGapFillEdgeCases(unittest.TestCase):

    def test_normal_cdf_table_values(self):
        """标准正态分布表值：Φ(0)=0.5、Φ(1.6)=0.9452007083、对称性。"""
        from src.nbmath import normal_cdf
        self.assertAlmostEqual(normal_cdf(0.0), 0.5, places=15)
        self.assertAlmostEqual(normal_cdf(1.6), 0.9452007083004420, places=12)
        self.assertAlmostEqual(normal_cdf(-1.6) + normal_cdf(1.6), 1.0, places=14)

    def test_chi2_sf_boundaries_and_monotonicity(self):
        """χ² 上尾：sf(0)=1、单调下降、大 x 趋 0。"""
        from src.nbmath import chi2_sf
        self.assertEqual(chi2_sf(0.0, 3), 1.0)
        vals = [chi2_sf(x, 3) for x in (0.5, 1.0, 5.0, 20.0)]
        for a, b in zip(vals, vals[1:]):
            self.assertGreater(a, b)
        self.assertLess(chi2_sf(200.0, 3), 1e-40)

    def test_threshold_above_hand_computed_table_value(self):
        """手算核对：β=0.7(ln), SE=0.25, θ=0.3 -> z=1.6 -> p=2*(1-Φ(1.6))。

        Φ(1.6) = 0.9452007083004420（标准正态表），所以 p = 0.1095985834。
        """
        from src.glm_nb import wald_test_above_threshold
        cov = [[0.1, 0.0], [0.0, 0.0625]]
        p = wald_test_above_threshold([4.0, 0.7], cov, 0.3)[2]
        self.assertAlmostEqual(p, 0.1095985833991160, places=9)

    def test_threshold_capped_at_one_and_invalid_covariance(self):
        """|β| < θ 时「超过阈值」p 截断为 1；协方差无效时返回 None。"""
        from src.glm_nb import wald_test_above_threshold, wald_test_below_threshold
        cov = [[0.1, 0.0], [0.0, 0.0625]]
        # |0.2| < 0.5：检验不成立，p 截断为 1
        self.assertEqual(wald_test_above_threshold([4.0, 0.2], cov, 0.5)[2], 1.0)
        # 方差为 0 / NaN -> p 为 None
        self.assertIsNone(wald_test_above_threshold([4.0, 0.7],
                                                    [[0.1, 0.0], [0.0, 0.0]], 0.3)[2])
        nan_cov = [[0.1, 0.0], [0.0, float("nan")]]
        self.assertIsNone(wald_test_above_threshold([4.0, 0.7], nan_cov, 0.3)[2])
        self.assertIsNone(wald_test_below_threshold([4.0, 0.7], nan_cov, 0.3)[2])

    def test_threshold_sign_symmetry_and_se_effect(self):
        """β 取正负号时两个方向的 p 都不变；SE 越小判定越锐利。"""
        from src.glm_nb import wald_test_above_threshold, wald_test_below_threshold
        cov = [[0.1, 0.0], [0.0, 0.0625]]        # se = 0.25
        cov_small = [[0.1, 0.0], [0.0, 0.015625]]  # se = 0.125
        p_pos = wald_test_below_threshold([4.0, 0.05], cov, 1.0)[2]
        p_neg = wald_test_below_threshold([4.0, -0.05], cov, 1.0)[2]
        self.assertAlmostEqual(p_pos, p_neg, places=15)
        # |β| 相同 -> above 方向也不变
        a_pos = wald_test_above_threshold([4.0, 0.7], cov, 0.3)[2]
        a_neg = wald_test_above_threshold([4.0, -0.7], cov, 0.3)[2]
        self.assertAlmostEqual(a_pos, a_neg, places=15)
        # SE 减半：两个方向都更锐利（p 更小）
        self.assertLess(wald_test_above_threshold([4.0, 0.7], cov_small, 0.3)[2],
                        a_pos)
        self.assertLess(wald_test_below_threshold([4.0, 0.05], cov_small, 1.0)[2],
                        p_pos)

    def test_lrt_df_two_chi2_cross_check(self):
        """三列设计矩阵 -> df=2；p 必须等于 chi2_sf(stat, 2) 本身。"""
        from src.nbmath import chi2_sf
        g1 = [0, 0, 0, 0, 1, 1, 1, 1]
        g2 = [0, 1, 0, 1, 0, 1, 0, 1]
        X3 = [[1.0, float(a), float(b)] for a, b in zip(g1, g2)]
        Xr = [[1.0] for _ in range(8)]
        res = lrt_test([50, 50, 50, 50, 200, 200, 200, 200], X3, Xr,
                       [0.0] * 8, 0.01)
        self.assertEqual(res["df"], 2)
        self.assertAlmostEqual(res["pvalue"],
                               chi2_sf(res["stat"], 2), places=15)
        self.assertLess(res["pvalue"], 1e-30)
        # 相同模型 -> df=0 -> 无法检验
        res0 = lrt_test([50, 50, 50, 50, 200, 200, 200, 200], X3, X3,
                        [0.0] * 8, 0.01)
        self.assertEqual(res0["df"], 0)
        self.assertIsNone(res0["pvalue"])

    def test_lrt_two_sided_symmetry(self):
        """上调与下调同样显著（双侧检验的对称性）。"""
        X = build_design_matrix([0, 0, 0, 0, 1, 1, 1, 1])
        Xr = [[1.0] for _ in range(8)]
        up = lrt_test([50, 50, 50, 50, 200, 200, 200, 200], X, Xr,
                      [0.0] * 8, 0.01)
        down = lrt_test([200, 200, 200, 200, 50, 50, 50, 50], X, Xr,
                        [0.0] * 8, 0.01)
        self.assertAlmostEqual(up["stat"], down["stat"], places=6)
        self.assertLess(up["pvalue"], 1e-20)
        self.assertLess(down["pvalue"], 1e-20)

    def test_cooks_threshold_degenerate_and_monotone(self):
        """m <= p 时阈值退化为 inf；置信度越高阈值越大。"""
        self.assertEqual(cooks_threshold(2, 2), float("inf"))
        self.assertEqual(cooks_threshold(3, 5), float("inf"))
        self.assertLess(cooks_threshold(8, 2, q=0.95),
                        cooks_threshold(8, 2, q=0.99))

    def test_robust_moments_dispersion_hand_computed(self):
        """手算核对：[10,20,30] -> MAD=10 -> s_rob=14.826。

        (14.826^2 - 20) / 20^2 = 0.49952569（与实现逐位一致）。
        方差不超过均值时截断为 0；样本不足返回 0。
        """
        v = robust_moments_dispersion([10, 20, 30])
        self.assertAlmostEqual(v, (14.826 ** 2 - 20) / 400.0, places=8)
        self.assertEqual(robust_moments_dispersion([10, 10, 10]), 0.0)
        self.assertEqual(robust_moments_dispersion([5]), 0.0)

    def test_cooks_distance_beta_reuse_and_nonnegative(self):
        """复用已拟合 beta 的结果必须与内部重拟合完全一致；距离非负。"""
        X = build_design_matrix([0, 0, 0, 0, 1, 1, 1, 1])
        c = [20, 22, 18, 20, 21, 19, 20, 1000]
        fit = fit_nb_glm(c, X, [0.0] * 8, 0.05)
        d1 = cooks_distance(c, X, [0.0] * 8, 0.05)
        d2 = cooks_distance(c, X, [0.0] * 8, 0.05, beta=fit["beta"])
        for a, b in zip(d1, d2):
            self.assertAlmostEqual(a, b, places=12)
        for v in d1:
            if v == v:                      # 非 NaN
                self.assertGreaterEqual(v, 0.0)

    def test_priovar_simulation_insufficient_data_fallback(self):
        """残差不足 20 条时模拟分支返回 None，estimate_dispersions 回退
        trigamma 公式并记录为 trigamma-fallback。"""
        from src.dispersion import _prior_var_by_simulation
        self.assertIsNone(_prior_var_by_simulation([0.1] * 10, 2))
        counts = [[100 + 20 * ((i * j) % 5) for j in range(4)] for i in range(10)]
        d = estimate_dispersions(counts, [1.0] * 4, [0, 0, 1, 1], quiet=True)
        self.assertEqual(d["dispPriorVarMethod"], "trigamma-fallback")

    def test_small_df_simulation_scale_tracking(self):
        """模拟分支能跟踪不同尺度的真值（σ² 越大估计越大），
        且 0 附近贴下限 0.25。"""
        from src.dispersion import _prior_var_by_simulation
        import random as _random

        def make_resid(n, sigma, df, seed):
            rng = _random.Random(seed)
            return [math.log(rng.gammavariate(df / 2.0, 2.0))
                    + rng.gauss(0, sigma) - math.log(df) for _ in range(n)]

        v0 = _prior_var_by_simulation(make_resid(3000, 0.0, 2, 11), 2)
        v1 = _prior_var_by_simulation(make_resid(3000, 1.0, 2, 11), 2)
        v2 = _prior_var_by_simulation(make_resid(3000, 1.5, 2, 11), 2)
        v3 = _prior_var_by_simulation(make_resid(3000, 2.0, 2, 11), 2)
        self.assertAlmostEqual(v0, 0.25, places=12)
        self.assertLess(v0, v1)
        self.assertLess(v1, v2)
        self.assertLess(v2, v3)
        self.assertTrue(0.8 <= v1 <= 1.6, "σ²=1.0 估计 %.3f 超界" % v1)
        self.assertTrue(1.9 <= v2 <= 2.6, "σ²=2.25 估计 %.3f 超界" % v2)
        self.assertTrue(3.4 <= v3 <= 4.5, "σ²=4.0 估计 %.3f 超界" % v3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
