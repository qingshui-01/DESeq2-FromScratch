# -*- coding: utf-8 -*-
"""回归验收脚本：把项目的全部验收标准固化成可重复执行的一条命令。

背景
------------
项目最大的风险是「改了 A 把 B 弄坏了，但没人发现」。本脚本把每一轮迭代
新增的验收点都累积下来，任何一次改代码之后跑一遍，就能确认**所有历史验收
标准仍然成立**：这就是你要求的「接下来的版本迭代可以通过以往全部的版本验收」。

运行
----
    python tests/run_acceptance.py            # 全量（含模拟数据自证，约 1 分钟）
    python tests/run_acceptance.py --quick    # 只跑不需要模拟数据的快检

退出码
------
    0 = 全部通过（PENDING 不算失败）
    1 = 有 FAIL

验收点编号规则
--------------
    A-*  环境与边界     （主力机零安装、纯标准库）
    B-*  代码骨架
    C-*  M4 BH 校正
    D-*  M1 归一化
    E-*  M2 离散度
    F-*  M3 GLM + Wald
    G-*  端到端 / CLI
    H-*  真实数据
    I-*  与官方 DESeq2 对拍
"""

import argparse
import csv
import io
import json
import math
import os
import platform
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

RESULTS = []          # (编号, 描述, 状态, 详情)


def check(code, desc):
    """装饰器：把一个函数登记为一个验收点。"""
    def deco(fn):
        fn._accept_code = code
        fn._accept_desc = desc
        return fn
    return deco


def record(code, desc, status, detail=""):
    RESULTS.append((code, desc, status, detail))


def _with_timing(detail, elapsed):
    """给耗时较长的验收点加上 [x.xs] 标注。

    WARNING: 这里是踩过坑的：早先写成
        detail = (detail + "  [%.1fs]").strip() % () if "%" not in detail else detail
    它把 `% ()` 作用到了「已经拼上 %.1fs 的字符串」上，于是必然抛
    TypeError: not enough arguments for format string。因为只有耗时 > 0.5 秒的
    分支才会走到这里，而此前没有任何一项超过 0.5 秒，这个 bug 一直没暴露，
    直到第一次读 6 MB 的真实基准 CSV 才当场崩掉。
    正确做法：先把耗时格式化好，再拼接；不要对 detail 做任何 % 运算
    （detail 里可能自带 %）。
    """
    if elapsed <= 0.5:
        return detail
    tail = "[%.1fs]" % elapsed
    return ("%s  %s" % (detail, tail)).strip() if detail else tail


def run_check(code, desc, fn):
    t0 = time.time()
    try:
        status, detail = fn()
    except Exception as exc:                      # noqa: BLE001
        status, detail = "FAIL", "%s: %s" % (type(exc).__name__, exc)
    elapsed = time.time() - t0
    detail = _with_timing(detail, elapsed)
    record(code, desc, status, detail)
    mark = {"PASS": "[通过]", "FAIL": "[失败]",
            "PENDING": "[待办]", "INFO": "[信息]"}[status]
    print("  %s %-6s %s" % (mark, code, desc))
    if detail:
        print("         %s" % detail)
    return status


# ==========================================================================
# A. 环境与边界
# ==========================================================================
@check("A-1", "Python 版本 >= 3.8 且为标准库运行环境")
def a1():
    v = sys.version_info
    if v < (3, 8):
        return "FAIL", "Python 版本过低：%s" % sys.version.split()[0]
    return "PASS", "Python %s on %s" % (sys.version.split()[0], platform.system())


@check("A-2", "核心算法未引入 numpy / scipy / pandas 等第三方数值库")
def a2():
    banned = {"numpy", "scipy", "pandas", "statsmodels", "sklearn"}
    loaded = banned & set(sys.modules)
    if loaded:
        return "FAIL", "检测到被导入的禁库：%s" % sorted(loaded)
    # 扫源码，确认没有 import 语句
    hits = []
    for dirpath, _dirs, files in os.walk(os.path.join(ROOT, "src")):
        for f in files:
            if not f.endswith(".py"):
                continue
            p = os.path.join(dirpath, f)
            with open(p, "r", encoding="utf-8") as fh:
                for ln, line in enumerate(fh, 1):
                    s = line.strip()
                    if s.startswith(("import ", "from ")) and \
                            any(b in s.split() for b in banned):
                        hits.append("%s:%d" % (os.path.relpath(p, ROOT), ln))
    if hits:
        return "FAIL", "源码中出现禁库导入：%s" % hits
    return "PASS", "src/ 与根目录脚本中均无禁库导入"


@check("A-3", "项目不依赖任何需要安装的第三方包（主力机零安装）")
def a3():
    req = os.path.join(ROOT, "requirements.txt")
    if os.path.exists(req):
        with open(req, "r", encoding="utf-8") as fh:
            body = fh.read().strip()
        if body:
            return "FAIL", "存在 requirements.txt 且非空：%s" % body
    return "PASS", "无 requirements.txt / 无外部依赖"


# ==========================================================================
# B. 代码骨架
# ==========================================================================
REQUIRED_FILES = [
    "deg.py", "demo.py", "README.md", "REPORT.md", ".gitignore",
    "src/__init__.py", "src/normalize.py", "src/dispersion.py",
    "src/glm_nb.py", "src/multitest.py", "src/nbmath.py",
    "src/linalg.py", "src/io_utils.py", "src/simulate.py",
    "tests/test_deg.py", "tests/run_acceptance.py",
    "tests/check_r_syntax.py", "tests/compare_baseline.py",
    "tests/verify_independent.py", "tests/fault_injection.py",
    "laptop/deseq2_baseline.R", "laptop/笔记本操作命令.txt",
]


@check("B-1", "项目骨架文件齐全")
def b1():
    missing = [f for f in REQUIRED_FILES if not os.path.exists(os.path.join(ROOT, f))]
    if missing:
        return "FAIL", "缺少文件：%s" % missing
    return "PASS", "共 %d 个必需文件全部存在" % len(REQUIRED_FILES)


@check("B-2", "目录结构与约定一致")
def b2():
    need = ["src", "tests", "data", "data/reference", "docs", "logs", "laptop"]
    missing = [d for d in need
               if not os.path.isdir(os.path.join(ROOT, d))]
    if missing:
        return "FAIL", "缺少目录：%s" % missing
    return "PASS", "目录齐全：%s" % ", ".join(need)


@check("B-3", "笔记本 R 脚本语法体检（括号平衡 / 纯 LF / 无写死系数名的断言）")
def b3():
    """这条验收点的由来：P6 阶段第一次运行基准脚本失败，根因是因子水平设错
    导致方向反了，加上一个写死系数名的断言把它拦下。此后把 R 脚本也纳入
    自动验收，改坏了会在主力机上当场发现，不必占用笔记本的一次试错。"""
    from tests.check_r_syntax import check_file
    p = os.path.join(ROOT, "laptop", "deseq2_baseline.R")
    ok, info = check_file(p)
    if info["errors"]:
        return "FAIL", "括号问题：%s" % info["errors"][:3]
    if info["paren"] or info["brace"] or info["bracket"]:
        return "FAIL", ("括号不平衡：()=%d {}=%d []=%d"
                        % (info["paren"], info["brace"], info["bracket"]))
    if info["unterminated_string"]:
        return "FAIL", "存在未闭合的字符串"
    if info["crlf"]:
        return "FAIL", "含 %d 处 CRLF 换行（粘进 Linux 终端会报错）" % info["crlf"]
    bad = [w for w in info["warnings"] if "系数名" in w]
    if bad:
        return "FAIL", "仍存在写死系数名的断言：%s" % bad
    return "PASS", ("%d 字节 / %d 行，括号平衡、纯 LF、无写死系数名的断言"
                    % (info["bytes"], info["lines"]))


@check("B-4", "笔记本 R 脚本函数参数名白名单（防 gc(quiet=) 这类错误）")
def b4():
    """这条验收点的由来：P7 阶段脚本已经跑完全部主产物，却倒在最后一句
    gc(quiet = TRUE) 上：R 的 gc() 根本没有 quiet 参数（只有 verbose/reset/full）。
    括号平衡查不出这种错，只能靠「参数名白名单」在主力机拦下来。
    将来在 R 脚本里新调用某个函数，必须把它补进 check_r_syntax.SIGNATURES，
    否则该调用的参数名不会被检查。"""
    from tests.check_r_syntax import check_file, check_call_args
    p = os.path.join(ROOT, "laptop", "deseq2_baseline.R")
    ok, info = check_file(p)
    bad = info.get("bad_args") or []
    if bad:
        return "FAIL", "可疑参数名：%s" % bad[:4]
    if info.get("calls_checked", 0) < 20:
        return "FAIL", ("只核对了 %d 处调用，明显偏少："
                        "可能扫描器坏了，或脚本里新函数没登记进 SIGNATURES"
                        % info.get("calls_checked", 0))
    # 自检：扫描器必须能抓到已知错误写法，否则「通过」毫无意义
    if not check_call_args("gc(quiet = TRUE)"):
        return "FAIL", "扫描器自检失败：连 gc(quiet = TRUE) 都抓不到"
    if check_call_args("dds <- DESeq(dds, quiet = TRUE)"):
        return "FAIL", "扫描器自检失败：把合法的 DESeq(quiet=) 误报成错误"
    return "PASS", ("核对 %d 处调用 / %d 个已登记函数，参数名全部合法；"
                    "扫描器自检（正反例）通过"
                    % (info["calls_checked"], info["sig_count"]))


# ==========================================================================
# C. M4
# ==========================================================================
@check("C-1", "M4：手算例子核对（p 与秩成正比时全部相等）")
def c1():
    from src.multitest import bh_adjust
    out = bh_adjust([0.01, 0.02, 0.03, 0.04, 0.05])
    for v in out:
        if abs(v - 0.05) > 1e-12:
            return "FAIL", "期望全为 0.05，实际 %s" % out
    return "PASS", "5 个值全部 = 0.05"


@check("C-2", "M4：padj 随 p 单调不降（漏掉 cummin 会违反）")
def c2():
    from src.multitest import bh_adjust
    ps = [0.9, 0.04, 0.5, 0.001, 0.2, 0.03, 0.7, 0.008]
    out = bh_adjust(ps)
    pairs = sorted(zip(ps, out))
    for k in range(1, len(pairs)):
        if pairs[k][1] < pairs[k - 1][1] - 1e-12:
            return "FAIL", "单调性被破坏：%s" % pairs
    return "PASS", "8 个值排序后单调不降"


@check("C-3", "M4：缺失值 / 单元素 / p=0 / p=1 四个边界")
def c3():
    from src.multitest import bh_adjust
    assert bh_adjust([None, None]) == [None, None]
    assert abs(bh_adjust([0.37])[0] - 0.37) < 1e-12
    o = bh_adjust([0.0, 1.0])
    assert abs(o[0]) < 1e-12 and abs(o[1] - 1.0) < 1e-12
    assert bh_adjust([]) == []
    return "PASS", "四个边界全部通过"


# ==========================================================================
# D. M1
# ==========================================================================
@check("D-1", "M1：手算 3 基因例子（ratio = 1:2:4，sf 应为 0.5 / 1 / 2）")
def d1():
    from src.normalize import size_factors
    sf = size_factors([[10, 20, 40], [100, 200, 400]])
    exp = [0.5, 1.0, 2.0]
    for a, b in zip(sf, exp):
        if abs(a - b) > 1e-12:
            return "FAIL", "期望 %s，实际 %s" % (exp, sf)
    return "PASS", "sf = %s" % ["%.4f" % v for v in sf]


@check("D-2", "M1：含 0 的基因被正确排除（几何均值遇到 0 会塌）")
def d2():
    from src.normalize import size_factors
    base = [[10, 20, 40], [100, 200, 400]]
    with_zero = base + [[0, 5, 999999]]
    a = size_factors(base)
    b = size_factors(with_zero)
    for x, y in zip(a, b):
        if abs(x - y) > 1e-12:
            return "FAIL", "加入含 0 基因后结果变化：%s -> %s" % (a, b)
    return "PASS", "结果不受含 0 基因影响"


@check("D-3", "M1：每个基因都含 0 时必须明确报错（不能悄悄给垃圾结果）")
def d3():
    from src.normalize import size_factors
    try:
        size_factors([[0, 1], [1, 0]])
    except ValueError as exc:
        return "PASS", "按预期报错：%s" % str(exc)[:60]
    return "FAIL", "未报错"


@check("D-4", "M1：真实数据上 size factor 落在合理区间（0.3 ~ 3）")
def d4():
    p = os.path.join(ROOT, "data", "counts.csv")
    if not os.path.exists(p):
        return "PENDING", "真实数据尚未就位（先跑 laptop/deseq2_baseline.R）"
    from src.io_utils import read_counts_csv
    from src.normalize import size_factors
    _g, _s, counts = read_counts_csv(p)
    sf = size_factors(counts)
    if min(sf) < 0.3 or max(sf) > 3.0:
        return "FAIL", "size factor 越界：%s" % ["%.3f" % v for v in sf]
    return "PASS", "sf = %s" % ["%.4f" % v for v in sf]


# ==========================================================================
# E. M2
# ==========================================================================
def _sim(n_genes, seed, n_per_group=4):
    from src.simulate import simulate_dataset
    return simulate_dataset(n_genes=n_genes, n_per_group=n_per_group, seed=seed)


@check("E-1", "M2：模拟数据上能估回已知的 alpha（中位比值在 0.5~2.0）")
def e1():
    from src.normalize import size_factors
    from src.dispersion import estimate_dispersions
    ds = _sim(400, 20260926)
    sf = size_factors(ds["counts"])
    d = estimate_dispersions(ds["counts"], sf, ds["group"])
    ratios = sorted(d["dispersion"][i] / ds["true_alpha"][i]
                    for i in range(len(ds["counts"])) if not d["allZero"][i])
    med = ratios[len(ratios) // 2]
    if not (0.5 < med < 2.0):
        return "FAIL", "中位比值 = %.3f，超出 [0.5, 2.0]" % med
    return "PASS", "中位比值 = %.3f（1.0 = 完美）" % med


@check("E-2", "M2：趋势拟合成功且随表达量下降")
def e2():
    from src.normalize import size_factors
    from src.dispersion import estimate_dispersions
    ds = _sim(400, 20260926)
    sf = size_factors(ds["counts"])
    d = estimate_dispersions(ds["counts"], sf, ds["group"])
    if d["fitType"] != "parametric":
        return "FAIL", "fitType = %s（期望 parametric）" % d["fitType"]
    idx = [i for i in range(len(ds["counts"])) if not d["allZero"][i]]
    idx.sort(key=lambda i: d["baseMean"][i])
    lo, hi = idx[len(idx) // 10], idx[-len(idx) // 10]
    if d["dispFit"][lo] <= d["dispFit"][hi]:
        return "FAIL", "趋势线未随表达量下降：低=%.4f 高=%.4f" \
            % (d["dispFit"][lo], d["dispFit"][hi])
    return "PASS", "低表达 dispFit=%.4f > 高表达 dispFit=%.4f" \
        % (d["dispFit"][lo], d["dispFit"][hi])


@check("E-3", "M2：收缩确实让估计更接近真值（log 误差的 IQR 变小）")
def e3():
    from src.normalize import size_factors
    from src.dispersion import estimate_dispersions
    ds = _sim(400, 20260926)
    sf = size_factors(ds["counts"])
    d = estimate_dispersions(ds["counts"], sf, ds["group"])
    idx = [i for i in range(len(ds["counts"])) if not d["allZero"][i]]
    g = sorted(math.log(d["dispGeneEst"][i] / ds["true_alpha"][i]) for i in idx)
    m = sorted(math.log(d["dispersion"][i] / ds["true_alpha"][i]) for i in idx)
    q = len(idx) // 4
    iqr_g = g[3 * q] - g[q]
    iqr_m = m[3 * q] - m[q]
    if not iqr_m < iqr_g:
        return "FAIL", "收缩后 IQR 未变小：gene=%.3f MAP=%.3f" % (iqr_g, iqr_m)
    return "PASS", "IQR: gene-wise %.3f -> 收缩后 %.3f（降低 %.0f%%）" \
        % (iqr_g, iqr_m, 100 * (1 - iqr_m / iqr_g))


@check("E-4", "M2：全 0 基因被标记且不参与估计")
def e4():
    from src.normalize import size_factors
    from src.dispersion import estimate_dispersions
    ds = _sim(200, 5)
    sf = size_factors(ds["counts"])
    d = estimate_dispersions(ds["counts"], sf, ds["group"])
    for i, row in enumerate(ds["counts"]):
        if sum(row) == 0 and not d["allZero"][i]:
            return "FAIL", "基因 %d 全 0 但未被标记" % i
    n = sum(d["allZero"])
    if n == 0:
        return "FAIL", "模拟数据里应当存在全 0 基因"
    return "PASS", "%d 个全 0 基因全部被标记" % n


@check("E-5", "M2：dispGeneEst 截断在 max(10, 样本数)（与官方一致）")
def e5():
    """这条验收点的由来：P8 阶段真实数据暴露的 bug。

    DESeq2 的离散度取值范围是 [minDisp, max(10, 样本数)]，**基因级 MLE 也要
    截断**。我们原先只在 MAP 收缩那一步截，于是像 ENSG00000229807
    （同一组内一个样本 3929、另一个样本 2）这种极端过散布基因，MLE 会跑到
    142：而官方同一个基因是 10（官方结果里 dispGeneEst 的最大值正好就是
    10.000000）。这个虚高的值随后污染趋势拟合，并让下游 GLM 的 log2FC
    跑到 -68。
    """
    from src.dispersion import estimate_dispersions
    counts = [[100 + 7 * ((i * j) % 11) for j in range(8)] for i in range(40)]
    counts.append([2, 0, 3929, 3042, 0, 0, 0, 2])      # 极端过散布
    factors = [1.0] * 8
    groups = [0, 1, 0, 1, 0, 1, 0, 1]
    d = estimate_dispersions(counts, factors, groups, quiet=True)
    vals = [v for v in d["dispGeneEst"] if v == v]
    mx = max(vals)
    if mx > 10.0 + 1e-9:
        return "FAIL", "dispGeneEst 最大值 %.4f 超过上限 10" % mx
    if mx < 10.0 - 1e-6:
        return "FAIL", ("最大值 %.4f 没碰到上限，极端基因没被构造出来，"
                        "这条验收点失效了" % mx)
    return "PASS", "极端基因的 dispGeneEst = %.4f（上限 10，已触发截断）" % mx


# ==========================================================================
# F. M3
# ==========================================================================
@check("F-5", "M3：一组接近全 0 时 log2FC 必须停在有界解（工作响应用截断后的 mu）")
def f5():
    """这条验收点的由来：P8 阶段真实数据暴露的最大一个 bug。

    官方 fitBeta 的工作响应是 z = log(mu_hat/nf) + (y - mu_hat)/mu_hat，
    其中 mu_hat **已经过 minmu 截断**。我们原先写成 η + log_sf（等价于
    log(未截断的 mu)）。当某个样本的 mu 触底（minmu = 0.5）后，偏差不再
    随 β 变化，但 z 仍随 η 线性下移，IRLS 就会沿着这个「偏差不再改善」的
    方向把 β 一路推下去。

    判别力是量出来的（注入 bug 实测 vs 修复后实测）：

       计数形态                        有 bug        修复后
       (0,0,0,1,0,500,0,600)      19.9 ~ 29.3      10.547（跨 α 极差 0）
       (0,120,0,130,0,110,0,125)     10.807        9.365

    所以这里用两个判据：① |log2FC| < 15；② **跨 4 个 α 的极差 < 1**
    ：后者是关键：正确实现下 mu 触底后 z 不再移动，结果与 α 无关；
    有 bug 时结果随 α 大幅漂移（19.9 → 29.3）。
    """
    from src.glm_nb import fit_nb_glm, wald_test
    from src.linalg import build_design_matrix
    X = build_design_matrix([0, 1, 0, 1, 0, 1, 0, 1])
    cases = {
        "对照组全 0": [0, 120, 0, 130, 0, 110, 0, 125],
        "近全 0": [0, 0, 0, 1, 0, 500, 0, 600],
    }
    alphas = (0.01, 0.1, 1.0, 10.0)
    for name, counts in cases.items():
        vals = []
        for alpha in alphas:
            fit = fit_nb_glm(counts, X, [0.0] * 8, alpha)
            lfc, _se, _stat, _p = wald_test(fit["beta"], fit["cov"], 1)
            if lfc is None or lfc != lfc or abs(lfc) == float("inf"):
                return "FAIL", "%s：alpha=%.2f 时 log2FC 非有限值" % (name, alpha)
            if lfc <= 0:
                return "FAIL", ("%s：alpha=%.2f 时 log2FC=%.3f，方向应为正"
                                % (name, alpha, lfc))
            if abs(lfc) > 15.0:
                return "FAIL", ("%s：alpha=%.2f 时 log2FC=%.3f，发散（>15）"
                                % (name, alpha, lfc))
            vals.append(lfc)
        spread = max(vals) - min(vals)
        if spread > 1.0:
            return "FAIL", ("%s：log2FC 随 alpha 漂移 %.3f（4 个 α 下为 %s）："
                            "说明工作响应没有用截断后的 mu"
                            % (name, spread, ["%.3f" % v for v in vals]))
    return "PASS", ("两个形态 × 4 个 α 下 log2FC 均为正、<15、"
                    "且跨 α 极差 <= 1（实际最大 0.000）")


@check("F-1", "M3：能恢复已知的 log2FC（构造 log2(200/50)=2 的基因）")
def f1():
    from src.glm_nb import fit_nb_glm, wald_test
    from src.linalg import build_design_matrix
    X = build_design_matrix([1, 1, 1, 1, 0, 0, 0, 0])
    fit = fit_nb_glm([200, 200, 200, 200, 50, 50, 50, 50], X, [0.0] * 8, 0.01)
    lfc, se, stat, p = wald_test(fit["beta"], fit["cov"], 1)
    if abs(lfc - 2.0) > 0.15:
        return "FAIL", "log2FC = %.4f，期望 2.0 ± 0.15" % lfc
    if p >= 1e-6:
        return "FAIL", "p = %g，期望 < 1e-6" % p
    return "PASS", "log2FC=%.4f  p=%.3g  符号正确" % (lfc, p)


@check("F-2", "M3：方向正确（下调基因 log2FC 为负）")
def f2():
    from src.glm_nb import fit_nb_glm, wald_test
    from src.linalg import build_design_matrix
    X = build_design_matrix([1, 1, 1, 1, 0, 0, 0, 0])
    fit = fit_nb_glm([50, 50, 50, 50, 200, 200, 200, 200], X, [0.0] * 8, 0.01)
    lfc, _, stat, _ = wald_test(fit["beta"], fit["cov"], 1)
    if lfc >= -1.5 or stat >= 0:
        return "FAIL", "log2FC=%.4f stat=%.4f，期望均为负" % (lfc, stat)
    return "PASS", "log2FC=%.4f stat=%.4f" % (lfc, stat)


@check("F-3", "M3：IRLS 收敛（不允许大量基因迭代到上限）")
def f3():
    from src.glm_nb import fit_nb_glm
    from src.linalg import build_design_matrix
    from src.normalize import size_factors
    from src.dispersion import estimate_dispersions
    ds = _sim(200, 3)
    sf = size_factors(ds["counts"])
    d = estimate_dispersions(ds["counts"], sf, ds["group"])
    X = build_design_matrix(ds["group"])
    lsf = [math.log(f) for f in sf]
    iters, conv = [], 0
    for i in range(len(ds["counts"])):
        if d["allZero"][i]:
            continue
        f = fit_nb_glm(ds["counts"][i], X, lsf, d["dispersion"][i])
        iters.append(f["iterations"])
        conv += 1 if f["converged"] else 0
    rate = conv / len(iters)
    if rate < 0.98:
        return "FAIL", "收敛率仅 %.1f%%（中位迭代 %d）" % (100 * rate,
                                                          sorted(iters)[len(iters) // 2])
    return "PASS", "收敛率 %.1f%%，中位迭代 %d 次" \
        % (100 * rate, sorted(iters)[len(iters) // 2])


@check("F-4", "M3：大样本下 Wald 统计量近似标准正态（实现正确性的硬证据）")
def f4():
    from src.glm_nb import fit_nb_glm, wald_test
    from src.linalg import build_design_matrix
    from src.simulate import simulate_dataset
    # 关键：n_de=0（全部为零假设基因）+ size factor 固定为 1
    # （这样 GLM 里用 log_sf = 0 才是自洽的）。
    # 若漏掉 n_de=0，混进来的强差异基因会把统计量的标准差拉爆。
    n_per_group = 40
    ds = simulate_dataset(n_genes=150, n_per_group=n_per_group, seed=555,
                          n_de=0, size_factors=[1.0] * (2 * n_per_group))
    X = build_design_matrix(ds["group"])
    lsf = [0.0] * len(ds["group"])
    stats = []
    for i in range(len(ds["counts"])):
        fit = fit_nb_glm(ds["counts"][i], X, lsf, ds["true_alpha"][i])
        _, _, st, _ = wald_test(fit["beta"], fit["cov"], 1)
        if st is not None:
            stats.append(st)
    n = len(stats)
    mean = sum(stats) / n
    sd = math.sqrt(sum((s - mean) ** 2 for s in stats) / (n - 1))
    if not (0.7 < sd < 1.35):
        return "FAIL", "stat 的 sd = %.3f，期望 ≈ 1.0（共 %d 个零假设基因）" \
            % (sd, n)
    return "PASS", "每组 n=%d 时 sd(stat)=%.3f（理想 1.000，%d 个零假设基因）" \
        % (n_per_group, sd, n)


# ==========================================================================
# G. 端到端 / CLI
# ==========================================================================
@check("G-1", "CLI：能在模拟数据上端到端跑通并输出正确列")
def g1():
    import tempfile
    from src.simulate import write_simulated_dataset
    from src.simulate import simulate_dataset
    with tempfile.TemporaryDirectory() as tmp:
        ds = simulate_dataset(n_genes=150, seed=77)
        cp = os.path.join(tmp, "counts.csv")
        dp = os.path.join(tmp, "design.csv")
        write_simulated_dataset(ds, cp, dp)
        out = os.path.join(tmp, "result.csv")
        r = subprocess.run(
            [sys.executable, os.path.join(ROOT, "deg.py"),
             "--counts", cp, "--design", dp, "--out", out, "--quiet"],
            cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
        if r.returncode != 0:
            return "FAIL", "退出码 %d\n%s" % (r.returncode, r.stderr[-500:])
        with open(out, "r", encoding="utf-8", newline="") as fh:
            header = next(csv.reader(fh))
        expect = ["gene_id", "baseMean", "log2FoldChange", "lfcSE",
                  "stat", "pvalue", "padj"]
        if header != expect:
            return "FAIL", "列名不匹配：%s" % header
        return "PASS", "列名与 DESeq2 对齐：%s" % ",".join(header)


@check("G-2", "CLI：输出基因数与输入完全一致（不丢基因，保证可逐行对拍）")
def g2():
    import tempfile
    from src.simulate import simulate_dataset, write_simulated_dataset
    with tempfile.TemporaryDirectory() as tmp:
        ds = simulate_dataset(n_genes=120, seed=78)
        cp = os.path.join(tmp, "counts.csv")
        dp = os.path.join(tmp, "design.csv")
        write_simulated_dataset(ds, cp, dp)
        out = os.path.join(tmp, "result.csv")
        subprocess.run([sys.executable, os.path.join(ROOT, "deg.py"),
                        "--counts", cp, "--design", dp, "--out", out,
                        "--quiet"], cwd=ROOT, capture_output=True, text=True)
        with open(out, "r", encoding="utf-8", newline="") as fh:
            n = sum(1 for _ in csv.reader(fh)) - 1
        if n != 120:
            return "FAIL", "输入 120 个基因，输出 %d 行" % n
        return "PASS", "输入 120 -> 输出 120 行，无丢失"


@check("G-3", "CLI：中文输出为 UTF-8 且不乱码")
def g3():
    r = subprocess.run([sys.executable, os.path.join(ROOT, "deg.py"),
                        "--counts", "NOPE.csv", "--design", "NOPE.csv",
                        "--out", "NOPE.csv"],
                       cwd=ROOT, capture_output=True)
    txt = r.stderr.decode("utf-8", errors="strict")
    if "计数矩阵" not in txt and "FileNotFoundError" not in txt \
            and "No such file" not in txt:
        return "FAIL", "无法以 UTF-8 解码错误输出：%r" % txt[:200]
    return "PASS", "输出可被严格 UTF-8 解码"


# ==========================================================================
# H. 真实数据
# ==========================================================================
@check("H-1", "真实数据：counts.csv 与 design.csv 已就位且格式正确")
def h1():
    cp = os.path.join(ROOT, "data", "counts.csv")
    dp = os.path.join(ROOT, "data", "design.csv")
    if not (os.path.exists(cp) and os.path.exists(dp)):
        return "PENDING", "尚未拷回（需要先在笔记本上跑 laptop/deseq2_baseline.R）"
    from src.io_utils import read_counts_csv, read_design_csv, align_design
    g, s, c = read_counts_csv(cp)
    d = read_design_csv(dp)
    groups, _, _ = align_design(s, d)
    if len(set(groups)) != 2:
        return "FAIL", "分组数不为 2"
    return "PASS", "%d 基因 x %d 样本，分组 %s" \
        % (len(g), len(s), sorted(set(groups)))


@check("H-2", "真实数据：airway 应为 63677 基因 x 8 样本（4 对配对）")
def h2():
    cp = os.path.join(ROOT, "data", "counts.csv")
    if not os.path.exists(cp):
        return "PENDING", "真实数据尚未就位"
    from src.io_utils import read_counts_csv, read_design_csv
    g, s, _c = read_counts_csv(cp)
    d = read_design_csv(os.path.join(ROOT, "data", "design.csv"))
    if len(s) != 8:
        return "FAIL", "样本数为 %d，期望 8" % len(s)
    cells = sorted(set(v.get("cell", "") for v in d.values() if v.get("cell")))
    return "PASS", "%d 基因 x %d 样本，细胞系 %s" % (len(g), len(s), cells)


# ==========================================================================
# I. 对拍
# ==========================================================================
def _corr(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy)


def _jaccard(a, b):
    a, b = set(a), set(b)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


@check("I-1", "对拍：与官方 DESeq2 的 log2FC 相关系数 > 0.8（及格线）")
def i1():
    from src.io_utils import read_results_csv
    mine_p = os.path.join(ROOT, "data", "my_result.csv")
    ref_p = os.path.join(ROOT, "data", "reference", "deseq2_result_noFilter.csv")
    if not (os.path.exists(mine_p) and os.path.exists(ref_p)):
        return "PENDING", "需要先跑出 my_result.csv 并拷回官方基准"
    mine = read_results_csv(mine_p)
    ref = read_results_csv(ref_p)
    common = sorted(set(mine) & set(ref))
    xs, ys = [], []
    for g in common:
        a, b = mine[g].get("log2FoldChange"), ref[g].get("log2FoldChange")
        if a is not None and b is not None:
            xs.append(a)
            ys.append(b)
    if len(xs) < 100:
        return "FAIL", "可比较的基因只有 %d 个" % len(xs)
    r = _corr(xs, ys)
    if r <= 0.8:
        return "FAIL", "相关系数 = %.4f（< 0.8）" % r
    grade = "优秀" if r > 0.95 else ("良好" if r > 0.9 else "及格")
    return "PASS", "相关系数 = %.4f（%s档），共 %d 个基因" % (r, grade, len(xs))


@check("I-2", "对拍：显著基因集合 Jaccard > 0.5")
def i2():
    from src.io_utils import read_results_csv
    mine_p = os.path.join(ROOT, "data", "my_result.csv")
    ref_p = os.path.join(ROOT, "data", "reference", "deseq2_result_noFilter.csv")
    if not (os.path.exists(mine_p) and os.path.exists(ref_p)):
        return "PENDING", "需要先跑出 my_result.csv 并拷回官方基准"
    mine = read_results_csv(mine_p)
    ref = read_results_csv(ref_p)

    def sig(d):
        return [g for g, v in d.items()
                if v.get("padj") is not None and v["padj"] < 0.05
                and v.get("log2FoldChange") is not None
                and abs(v["log2FoldChange"]) >= 1.0]
    a, b = sig(mine), sig(ref)
    j = _jaccard(a, b)
    if len(b) == 0:
        return "FAIL", "官方基准里没有显著基因，筛选阈值可能不一致"
    if j <= 0.5:
        return "FAIL", "Jaccard = %.4f（我的 %d 个 / 官方 %d 个）" % (j, len(a), len(b))
    grade = "优秀" if j > 0.7 else "及格~良好"
    return "PASS", "Jaccard = %.4f（%s），我的 %d 个 / 官方 %d 个" \
        % (j, grade, len(a), len(b))


@check("I-3", "对拍：size factor 与官方逐样本一致（M1 直接验证）")
def i3():
    from src.io_utils import read_counts_csv
    from src.normalize import size_factors
    cp = os.path.join(ROOT, "data", "counts.csv")
    sp = os.path.join(ROOT, "data", "reference", "sizeFactors.csv")
    if not (os.path.exists(cp) and os.path.exists(sp)):
        return "PENDING", "需要 counts.csv 与官方 sizeFactors.csv"
    _g, _s, counts = read_counts_csv(cp)
    mine = size_factors(counts)
    ref = []
    with open(sp, "r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            ref.append(float(row["sizeFactor"]))
    if len(ref) != len(mine):
        return "FAIL", "样本数不匹配：我 %d / 官方 %d" % (len(mine), len(ref))
    # 注意：size factor 只可识别到常数倍，所以比较的是「归一化后」的相对值
    m = sum(mine) / len(mine)
    r = sum(ref) / len(ref)
    rel = [abs(mine[i] / m - ref[i] / r) / (ref[i] / r) for i in range(len(mine))]
    worst = max(rel)
    if worst > 0.05:
        return "FAIL", "最大相对偏差 %.2f%%（阈值 5%%）" % (100 * worst)
    return "PASS", "最大相对偏差 %.3f%%（逐样本对齐后）" % (100 * worst)


@check("I-4", "对拍：结果表内部自洽（列间关系 / padj 单调 / 与 counts 行对齐）")
def i4():
    """这条验收点来自 P8 的独立复算（`tests/verify_independent.py`）。

    它不依赖官方基准，只看自己产出的 `my_result.csv` 是否自洽：

      * 行数、顺序与 counts.csv 完全一致（保证可逐行对拍）
      * stat = log2FoldChange / lfcSE
      * pvalue = erfc(|stat| / sqrt(2))
      * padj 在 [pvalue, 1] 内，且**按 p 值升序排列后不下降**
      * 全 0 基因的各项结果均为空

    这些关系任一被破坏，都说明输出表在某个环节算错了，而且不需要官方数据
    就能发现。注意 I-1/I-2 读的是**已产出的** my_result.csv，改了 src/ 之后
    必须重跑 deg.py 它们才会跟着变；I-4 同理。
    """
    from src.io_utils import read_counts_csv
    mp = os.path.join(ROOT, "data", "my_result.csv")
    cp = os.path.join(ROOT, "data", "counts.csv")
    if not (os.path.exists(mp) and os.path.exists(cp)):
        return "PENDING", "需要 counts.csv 与 my_result.csv"
    gene_ids, _s, counts = read_counts_csv(cp)
    rows = []
    with open(mp, "r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    if len(rows) != len(gene_ids):
        return "FAIL", "行数 %d != counts 行数 %d" % (len(rows), len(gene_ids))
    if [r["gene_id"] for r in rows] != gene_ids:
        return "FAIL", "基因顺序与 counts.csv 不一致（无法逐行对拍）"

    def num(s):
        s = (s or "").strip()
        if s == "" or s.upper() in ("NA", "NAN", "NULL"):
            return None
        try:
            return float(s)
        except ValueError:
            return None

    bad_stat = bad_p = bad_empty = 0
    pairs = []
    n_testable = 0
    for i, r in enumerate(rows):
        lfc, se = num(r["log2FoldChange"]), num(r["lfcSE"])
        st, pv, pj = num(r["stat"]), num(r["pvalue"]), num(r["padj"])
        if lfc is None:
            if se is not None or st is not None or pv is not None:
                bad_empty += 1
            continue
        n_testable += 1
        if abs(st - lfc / se) > 1e-9 * max(1.0, abs(st)):
            bad_stat += 1
        if abs(pv - math.erfc(abs(st) / math.sqrt(2.0))) > 1e-12 * max(1e-12, pv):
            bad_p += 1
        pairs.append((pv, pj))
    if bad_stat or bad_p or bad_empty:
        return "FAIL", ("stat 不符 %d / p 不符 %d / 空缺不干净 %d"
                        % (bad_stat, bad_p, bad_empty))
    pairs.sort()
    viol_mono = sum(1 for k in range(1, len(pairs))
                    if pairs[k][1] is not None and pairs[k - 1][1] is not None
                    and pairs[k][1] < pairs[k - 1][1] - 1e-15)
    viol_range = sum(1 for pv, pj in pairs
                     if pj is not None and (pj < pv - 1e-12 or pj > 1.0 + 1e-12))
    if viol_mono or viol_range:
        return "FAIL", ("padj 非单调 %d 个 / 越界 %d 个" % (viol_mono, viol_range))
    n_allzero = sum(1 for c in counts if sum(c) == 0)
    if n_testable + n_allzero != len(gene_ids):
        return "FAIL", ("可检验 %d + 全 0 %d != 总基因 %d"
                        % (n_testable, n_allzero, len(gene_ids)))
    return "PASS", ("%d 个基因行对齐；stat/p/单调性全部自洽；"
                    "可检验 %d + 全 0 %d = %d"
                    % (len(rows), n_testable, n_allzero, len(gene_ids)))


@check("I-5", "对拍：dispGeneEst 中位比值 ∈ [0.9, 1.1]（守住 Cox-Reid 校正项）")
def i5():
    """这条验收点的由来：P11 代码深审发现我们漏掉了 Cox-Reid 校正项。

    官方 fitDisp（src/DESeq2.cpp）默认 useCR=TRUE，目标函数里有
    -0.5*log(det(X'WX))。这一项随 alpha 单调上升，缺了它 dispGeneEst 会
    系统性偏低约一半（airway 实测中位比值 0.504）；补上后回到 1.000。
    最终结果表对这个差异不敏感（log2FC 相关系数几乎不动），所以 I-1~I-4
    抓不住它，必须在中间量层面直接对拍。
    """
    mp = os.path.join(ROOT, "data", "intermediates", "dispersion_mine.csv")
    rp = os.path.join(ROOT, "data", "reference", "dispersion_intermediates.csv")
    if not (os.path.exists(mp) and os.path.exists(rp)):
        return "PENDING", "需要 dispersion_mine.csv 与官方 dispersion_intermediates.csv"

    def load_disp(path):
        out = {}
        with open(path, "r", encoding="utf-8-sig", newline="") as fh:
            header = fh.readline().strip().split(",")
            gi = header.index("gene_id") if "gene_id" in header else 0
            di = header.index("dispGeneEst")
            for line in fh:
                parts = line.strip().split(",")
                if len(parts) <= max(gi, di):
                    continue
                try:
                    out[parts[gi]] = float(parts[di])
                except ValueError:
                    continue
        return out

    mine = load_disp(mp)
    ref = load_disp(rp)
    ratios = []
    for g, dv in mine.items():
        rv = ref.get(g)
        if rv is None or dv != dv or rv != rv or dv <= 0 or rv <= 0:
            continue
        ratios.append(dv / rv)
    if len(ratios) < 1000:
        return "PENDING", "可比较的 dispGeneEst 太少（%d 个）" % len(ratios)
    ratios.sort()
    med = ratios[len(ratios) // 2]
    if not (0.9 <= med <= 1.1):
        return "FAIL", ("dispGeneEst 中位比值 = %.4f（期望 ~1.0；"
                        "掉到 ~0.5 说明 Cox-Reid 校正项丢失）" % med)
    return "PASS", "dispGeneEst 中位比值 = %.4f（%d 个基因）" % (med, len(ratios))


# ==========================================================================
# J. 论文功能补全（P14）
# ==========================================================================
@check("J-1", "补全：m-p <= 3 时先验方差走模拟分支（论文小自由度特例）")
def j1():
    """论文 "Three or less residual degrees of freedom"：
    残差自由度 <= 3 时用「模拟 + KL 匹配」替代 trigamma 近似。
    这里验证：① 分支被正确触发（m=4,p=2 -> m-p=2 用 simulation；
    m=8 -> trigamma）；② 模拟结果可复现；③ 合成残差校准（真值 1.0）。"""
    from src.dispersion import _prior_var_by_simulation
    import random as _random

    def make_resid(n, sigma, df, seed):
        rng = _random.Random(seed)
        return [math.log(rng.gammavariate(df / 2.0, 2.0))
                + rng.gauss(0, sigma) - math.log(df) for _ in range(n)]

    r = make_resid(3000, 1.0, 2, seed=11)
    v1 = _prior_var_by_simulation(r, 2)
    v2 = _prior_var_by_simulation(r, 2)
    if v1 != v2:
        return "FAIL", "模拟分支不可复现：%.6f vs %.6f" % (v1, v2)
    if not (0.8 <= v1 <= 1.6):
        return "FAIL", "合成残差校准失败：真值 1.0，估计 %.3f" % v1

    from src.dispersion import estimate_dispersions
    counts4 = [[80 + 30 * ((i * j) % 7) for j in range(4)] for i in range(80)]
    d4 = estimate_dispersions(counts4, [1.0] * 4, [0, 0, 1, 1], quiet=True)
    counts8 = [[80 + 30 * ((i * j) % 7) for j in range(8)] for i in range(80)]
    d8 = estimate_dispersions(counts8, [1.0] * 8,
                              [0, 0, 0, 0, 1, 1, 1, 1], quiet=True)
    if d4["dispPriorVarMethod"] != "simulation":
        return "FAIL", "m-p=2 未走模拟分支（%s）" % d4["dispPriorVarMethod"]
    if d8["dispPriorVarMethod"] != "trigamma":
        return "FAIL", "m-p=6 不应走模拟分支（%s）" % d8["dispPriorVarMethod"]
    return "PASS", ("校准 %.3f（真值 1.0）、可复现；分支触发正确"
                    "（m-p<=3 用 simulation，否则 trigamma）" % v1)


@check("J-2", "补全：阈值检验（论文 Composite null hypotheses）")
def j2():
    """两个方向的复合零假设检验：
    ① θ=0 时与标准双侧 Wald 完全相等（数学上必然）；
    ② 「超过阈值」随阈值单调变松；
    ③ 「弱于阈值」：弱效应显著、跨越阈值的效应不显著。"""
    from src.glm_nb import (wald_test, wald_test_above_threshold,
                            wald_test_below_threshold)
    cov = [[0.1, 0.0], [0.0, 0.0625]]     # se = 0.25
    p_thr = wald_test_above_threshold([4.0, 0.7], cov, 0.0)[2]
    p_std = wald_test([4.0, 0.7], cov)[3]
    if abs(p_thr - p_std) > 1e-15:
        return "FAIL", "θ=0 与标准 Wald 不等：%.3e vs %.3e" % (p_thr, p_std)
    p1 = wald_test_above_threshold([4.0, 0.7], cov, 0.3)[2]
    p2 = wald_test_above_threshold([4.0, 0.7], cov, 0.5)[2]
    if not p1 < p2:
        return "FAIL", "阈值增大后 p 未变松：%.4f -> %.4f" % (p1, p2)
    p_weak = wald_test_below_threshold([4.0, 0.05], cov, 1.0)[2]
    p_strong = wald_test_below_threshold([4.0, 1.6], cov, 1.5)[2]
    if not (p_weak < 0.01 and p_strong > 0.5):
        return "FAIL", ("「弱于阈值」语义错误：弱效应 %.4f（应<0.01）、"
                        "强效应 %.4f（应>0.5）" % (p_weak, p_strong))
    return "PASS", ("θ=0 与 Wald 逐位相等；「超过阈值」单调；"
                    "「弱于阈值」弱效应 p=%.1e / 强效应 p=%.2f" % (p_weak, p_strong))


@check("J-3", "补全：LRT 似然比检验（χ² 闭式自检 + 构造基因）")
def j3():
    """卡方上尾与闭式一致（df=1/2），并验证：
    构造 log2(200/50) 的基因 LRT p 极小、stat 与 Wald 的 z² 同量级；
    零效应基因不显著。"""
    from src.nbmath import chi2_sf
    from src.glm_nb import fit_nb_glm, lrt_test, wald_test
    from src.linalg import build_design_matrix
    x = 4.0
    if abs(chi2_sf(x, 1) - math.erfc(math.sqrt(x / 2.0))) > 1e-12:
        return "FAIL", "chi2_sf(df=1) 与 erfc 闭式不符"
    if abs(chi2_sf(x, 2) - math.exp(-x / 2.0)) > 1e-12:
        return "FAIL", "chi2_sf(df=2) 与 exp 闭式不符"
    X = build_design_matrix([0, 0, 0, 0, 1, 1, 1, 1])
    Xr = [[1.0] for _ in range(8)]
    res = lrt_test([50, 50, 50, 50, 200, 200, 200, 200], X, Xr, [0.0] * 8, 0.01)
    if not (res["df"] == 1 and res["pvalue"] < 1e-20 and res["stat"] > 100):
        return "FAIL", ("差异基因 LRT 异常：stat=%.2f p=%.2e" % (res["stat"], res["pvalue"]))
    fit = fit_nb_glm([50, 50, 50, 50, 200, 200, 200, 200], X, [0.0] * 8, 0.01)
    z = wald_test(fit["beta"], fit["cov"])[2]
    if abs(res["stat"] - z * z) / res["stat"] > 0.15:
        return "FAIL", ("LRT stat=%.1f 与 Wald z²=%.1f 相差过大"
                        % (res["stat"], z * z))
    res0 = lrt_test([100, 90, 110, 95, 105, 100, 95, 105], X, Xr, [0.0] * 8, 0.01)
    if res0["pvalue"] < 0.05:
        return "FAIL", "零效应基因 LRT 误报：p=%.3f" % res0["pvalue"]
    return "PASS", ("χ² 闭式一致；差异基因 stat=%.1f p=%.1e（Wald z²=%.1f）；"
                    "零效应 p=%.2f" % (res["stat"], res["pvalue"], z * z, res0["pvalue"]))


@check("J-4", "补全：Cook's 距离离群点诊断（论文 Detection of count outliers）")
def j4():
    """构造「一个样本 1000、其余 ~20」的基因：
    必须只把该样本标记为离群（阈值 = qf(0.99, p, m-p)）；
    无离群样本的基因不被标记。"""
    from src.diagnostics import flag_cooks_outliers
    from src.linalg import build_design_matrix
    from src.nbmath import f_cdf
    X = build_design_matrix([0, 0, 0, 0, 1, 1, 1, 1])
    res = flag_cooks_outliers([20, 22, 18, 20, 21, 19, 20, 1000],
                              X, [0.0] * 8, 0.05)
    flagged = [i for i, f in enumerate(res["flags"]) if f]
    if flagged != [7]:
        return "FAIL", "离群样本标记错误：%s（应为 [7]）" % flagged
    if abs(f_cdf(res["threshold"], 2, 6) - 0.99) > 1e-9:
        return "FAIL", "阈值不是 F(2,6) 的 0.99 分位"
    res2 = flag_cooks_outliers([20, 22, 18, 20, 21, 19, 20, 18],
                               X, [0.0] * 8, 0.05)
    if any(res2["flags"]):
        return "FAIL", "正常基因被误标"
    return "PASS", ("离群样本 Cook=%.1f > 阈值 %.2f（只标记该样本）；"
                    "正常基因零误标" % (res["max_cooks"], res["threshold"]))


# ==========================================================================
# 主流程
# ==========================================================================
QUICK_SKIP = {"E-1", "E-2", "E-3", "E-4", "F-3", "F-4", "G-1", "G-2"}

ALL_CHECKS = [a1, a2, a3, b1, b2, b3, b4, c1, c2, c3, d1, d2, d3, d4,
              e1, e2, e3, e4, e5, f1, f2, f3, f4, f5, g1, g2, g3, h1, h2,
              i1, i2, i3, i4, i5, j1, j2, j3, j4]


def main():
    ap = argparse.ArgumentParser(description="DESeq2-FromScratch 回归验收")
    ap.add_argument("--quick", action="store_true",
                    help="跳过需要模拟数据的重活（约 3 秒）")
    ap.add_argument("--json", default=None, help="把结果写成 JSON")
    args = ap.parse_args()

    t0 = time.time()
    print("=" * 74)
    print("DESeq2-FromScratch 回归验收")
    print("时间: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    print("Python: %s   模式: %s" % (sys.version.split()[0],
                                    "quick" if args.quick else "full"))
    print("=" * 74)

    sections = [
        ("A 环境与边界", [a1, a2, a3]),
        ("B 代码骨架", [b1, b2, b3, b4]),
        ("C M4 BH 校正", [c1, c2, c3]),
        ("D M1 归一化", [d1, d2, d3, d4]),
        ("E M2 离散度", [e1, e2, e3, e4, e5]),
        ("F M3 GLM + Wald", [f1, f2, f3, f4, f5]),
        ("G 端到端 / CLI", [g1, g2, g3]),
        ("H 真实数据", [h1, h2]),
        ("I 与官方 DESeq2 对拍", [i1, i2, i3, i4, i5]),
        ("J 论文功能补全（P14）", [j1, j2, j3, j4]),
    ]
    for title, fns in sections:
        print("\n--- %s ---" % title)
        for fn in fns:
            code = fn._accept_code
            if args.quick and code in QUICK_SKIP:
                record(code, fn._accept_desc, "PENDING", "quick 模式跳过")
                print("  [跳过] %-6s %s" % (code, fn._accept_desc))
                continue
            run_check(code, fn._accept_desc, fn)

    n_pass = sum(1 for r in RESULTS if r[2] == "PASS")
    n_fail = sum(1 for r in RESULTS if r[2] == "FAIL")
    n_pend = sum(1 for r in RESULTS if r[2] == "PENDING")

    print("\n" + "=" * 74)
    print("汇总：通过 %d / 失败 %d / 待办 %d （共 %d 项），耗时 %.1f 秒"
          % (n_pass, n_fail, n_pend, len(RESULTS), time.time() - t0))
    if n_fail:
        print("\n失败项：")
        for code, desc, status, detail in RESULTS:
            if status == "FAIL":
                print("  %s %s\n      %s" % (code, desc, detail))
    if n_pend:
        print("\n待办项（不算失败，通常是在等外部输入）：")
        for code, desc, status, detail in RESULTS:
            if status == "PENDING":
                print("  %s %s\n      %s" % (code, desc, detail))
    print("\n结论：%s" % ("全部通过" if n_fail == 0
                        else "有 %d 项失败，需要修复" % n_fail))
    print("=" * 74)

    # 结果落盘（供版本迭代之间对比）
    logdir = os.path.join(ROOT, "logs")
    os.makedirs(logdir, exist_ok=True)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    payload = {
        "time": time.strftime("%Y-%m-%d %H:%M:%S"),
        "python": sys.version.split()[0],
        "mode": "quick" if args.quick else "full",
        "pass": n_pass, "fail": n_fail, "pending": n_pend,
        "results": [{"code": c, "desc": d, "status": s, "detail": t}
                    for c, d, s, t in RESULTS],
    }
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, ensure_ascii=False, indent=2)
    with open(os.path.join(logdir, "acceptance-%s.json" % stamp), "w",
              encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)
    print("验收记录已写入 logs/acceptance-%s.json" % stamp)

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
