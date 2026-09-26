"""DESeq2-FromScratch：纯 Python 标准库实现的差异表达分析。

四个核心模块：
    normalize.py：M1 归一化（median-of-ratios size factor）
    dispersion.py：M2 离散度估计（gene-wise MLE + 趋势拟合 + 经验贝叶斯收缩）
    glm_nb.py：M3 负二项 GLM（IRLS）+ Wald 检验
    multitest.py：M4 BH 多重检验校正

约束：核心算法只允许使用 Python 标准库（math / statistics / csv / json / argparse /
unittest / random），禁止 numpy、scipy 等第三方数值库。
"""

__version__ = "1.0.0"
