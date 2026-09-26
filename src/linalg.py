"""极小的线性代数工具（纯标准库）。

这个项目里的矩阵都很小：设计矩阵是 m（样本数，8）× p（系数个数，2），
所以用最朴素的高斯消元就够了，不需要任何第三方库。

约定：矩阵用「行优先的 list of list」表示，A[i][j] 表示第 i 行第 j 列。
"""

import math


def solve(A, b):
    """解线性方程组 A x = b（A 为方阵，b 为向量）。返回 x。

    用带部分主元的高斯-约当消元。A 会被复制，不修改入参。
    """
    n = len(A)
    if n == 0:
        raise ValueError("solve: 矩阵为空")
    if len(b) != n:
        raise ValueError("solve: 维度不匹配")

    M = [list(row) + [float(b[i])] for i, row in enumerate(A)]
    for col in range(n):
        # 选主元
        pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[pivot][col]) < 1e-300:
            raise ValueError("solve: 矩阵奇异，无法求解")
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
        # 归一化主元行
        pv = M[col][col]
        inv_pv = 1.0 / pv
        row_col = M[col]
        for j in range(col, n + 1):
            row_col[j] *= inv_pv
        # 消去其他行
        for r in range(n):
            if r == col:
                continue
            factor = M[r][col]
            if factor == 0.0:
                continue
            row_r = M[r]
            for j in range(col, n + 1):
                row_r[j] -= factor * row_col[j]

    return [M[i][n] for i in range(n)]


def invert(A):
    """求方阵的逆。对奇异矩阵抛 ValueError。"""
    n = len(A)
    if n == 0:
        raise ValueError("invert: 矩阵为空")
    M = [list(row) + [1.0 if i == j else 0.0 for j in range(n)]
         for i, row in enumerate(A)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[pivot][col]) < 1e-300:
            raise ValueError("invert: 矩阵奇异")
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
        pv = M[col][col]
        inv_pv = 1.0 / pv
        row_col = M[col]
        for j in range(col, 2 * n):
            row_col[j] *= inv_pv
        for r in range(n):
            if r == col:
                continue
            factor = M[r][col]
            if factor == 0.0:
                continue
            row_r = M[r]
            for j in range(col, 2 * n):
                row_r[j] -= factor * row_col[j]

    return [[M[i][n + j] for j in range(n)] for i in range(n)]


def determinant(A):
    """方阵行列式（带部分主元的高斯消元）。

    用于 Cox-Reid 校正项里的 det(X' W X)（DESeq2 src/DESeq2.cpp fitDisp）。
    矩阵奇异地返回 0.0。
    """
    n = len(A)
    if n == 0:
        raise ValueError("determinant: 矩阵为空")
    M = [row[:] for row in A]
    det = 1.0
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[pivot][col]) < 1e-300:
            return 0.0
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
            det = -det
        pv = M[col][col]
        det *= pv
        for r in range(col + 1, n):
            factor = M[r][col] / pv
            if factor == 0.0:
                continue
            row_r = M[r]
            row_c = M[col]
            for j in range(col, n):
                row_r[j] -= factor * row_c[j]
    return det


def ols_hat_matrix(X):
    """最小二乘的投影算子 H = (X'X)^-1 X'，形状 p×m。

    用于把一组观测向量投影到设计矩阵张成的空间上（等价于 R 的
    linearModelMu()，DESeq2 源码 R/core.R 用它给 roughDispEstimate 估 mu）。
    """
    p = len(X[0]) if X else 0
    m = len(X)
    # XtX = X' X
    XtX = [[sum(X[i][a] * X[i][b] for i in range(m)) for b in range(p)]
           for a in range(p)]
    # Xt = X'
    Xt = [[X[i][j] for i in range(m)] for j in range(p)]
    XtX_inv = invert(XtX)
    # H = XtX_inv @ Xt
    return [[sum(XtX_inv[a][k] * Xt[k][j] for k in range(p)) for j in range(m)]
            for a in range(p)]


def mat_vec(M, v):
    """矩阵乘向量。"""
    return [sum(row[j] * v[j] for j in range(len(v))) for row in M]


def is_full_rank(X, tol=1e-10):
    """检查设计矩阵是否列满秩（秩 = 列数）。"""
    if not X:
        return False
    m = len(X)
    p = len(X[0])
    if p > m:
        return False
    XtX = [[sum(X[i][a] * X[i][b] for i in range(m)) for b in range(p)]
           for a in range(p)]
    # 用高斯消元算秩
    M = [row[:] for row in XtX]
    rank = 0
    for col in range(p):
        pivot = None
        for r in range(rank, p):
            if abs(M[r][col]) > tol:
                pivot = r
                break
        if pivot is None:
            continue
        M[rank], M[pivot] = M[pivot], M[rank]
        pv = M[rank][col]
        for j in range(col, p):
            M[rank][j] /= pv
        for r in range(p):
            if r != rank and abs(M[r][col]) > 0:
                f = M[r][col]
                for j in range(col, p):
                    M[r][j] -= f * M[rank][j]
        rank += 1
        if rank == p:
            break
    return rank == p


def build_design_matrix(groups):
    """由 0/1 分组向量构建设计矩阵 [截距, 组别]。

    与 R 的 model.matrix(~ group) 等价：第 1 列全 1，第 2 列为 treat 指示。
    log2FoldChange 的含义 = log2(group==1 的均值 / group==0 的均值)。
    """
    return [[1.0, float(g)] for g in groups]
