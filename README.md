# DESeq2-FromScratch

**不调用 DESeq2，用纯 Python 标准库从零实现 RNA-seq 差异表达分析的核心统计流程，
然后与官方 DESeq2 的结果对拍。**

这不是一个「用来替代 DESeq2 的工具」。DESeq2 更快、更全、经过十几年实战检验，
做分析就该用它。本项目的产出是**一份用数字写成的证据**：
「我的实现和官方在 log2FoldChange 上的相关系数是 X，显著基因集合重合率是 Y，
不一致的基因集中在低表达区，原因是 Z。」

---

## 一、它做什么

输入一张计数矩阵（基因 × 样本）和一张分组表，输出每个基因的差异表达统计量：

```powershell
python deg.py --counts data/counts.csv --design data/design.csv --out data/my_result.csv
```

输出列与 DESeq2 `results()` **完全对齐**，因此可以逐列对拍：

| 列 | 含义 |
|---|---|
| `gene_id` | 基因 ID |
| `baseMean` | 该基因在所有样本上的归一化计数均值 |
| `log2FoldChange` | log2（组1均值 / 组0均值） |
| `lfcSE` | log2FC 的标准误 |
| `stat` | Wald 统计量 = log2FC / lfcSE |
| `pvalue` | 双侧 p 值，标准正态 |
| `padj` | BH 校正后的 p 值（FDR） |

---

## 二、四个核心模块

| 模块 | 文件 | 做什么 | 对拍目标（官方中间量） |
|---|---|---|---|
| **M1** | `src/normalize.py` | 归一化：median-of-ratios size factor | `sizeFactors.csv` |
| **M2** | `src/dispersion.py` | 离散度：gene-wise MLE → 趋势拟合 → 经验贝叶斯收缩 | `dispGeneEst` / `dispFit` / `dispersion` |
| **M3** | `src/glm_nb.py` | 负二项 GLM（IRLS/Fisher scoring）+ Wald 检验 | 结果表的 `log2FoldChange` 等列 |
| **M4** | `src/multitest.py` | Benjamini-Hochberg 多重检验校正 | `padj` |

辅助模块：

| 文件 | 作用 |
|---|---|
| `src/nbmath.py` | 负二项对数似然、trigamma、mad、中位数 |
| `src/linalg.py` | 极小规模线性代数（高斯消元、求逆、最小二乘投影） |
| `src/io_utils.py` | CSV 读写与样本对齐 |
| `src/simulate.py` | **方案 C 模拟数据生成器**（真值已知，用于自证 M2/M3） |

**算法里的每一条公式，都是对着 DESeq2 官方源码逐条核过的**（不是凭记忆写的），
推导依据写在每个模块的 docstring 里，例如
[M1 的确切规则](src/normalize.py)、
[M2 的三步与先验方差](src/dispersion.py)。

---

## 三、环境与运行

**零安装。** 只需要 Python 3.8+ 标准库（本项目在 Python 3.12.10 / Windows 下开发验证）。

```powershell
git clone <repo>            # 或直接拷贝目录
cd DESeq2-FromScratch

# 0) 一键演示：用模拟数据把整条流程跑通，并和已知真值对照（几秒）
python demo.py

# 1) 单元测试（50 项）
python -m unittest tests.test_deg -v

# 2) 回归验收：跑完全部历史验收点（38 项）
python tests/run_acceptance.py

# 3) 真实数据（需要先按第四节拿到数据）
python deg.py --counts data/counts.csv --design data/design.csv `
              --out data/my_result.csv --dump-intermediates data/intermediates
```

`demo.py` 的输出会直接告诉你 M1/M2/M3 各算得对不对（对照真值），
适合第一次接触本项目时用来看清全貌。

**硬性约束**：核心算法**只允许使用标准库**，
禁止 numpy / scipy / pandas。`tests/run_acceptance.py` 的 A-2 验收点会扫描源码
强制执行这一条。

---

## 四、输入端：数据从哪来

用的是 Bioconductor 官方的 `airway` 数据集（人类气道平滑肌细胞、地塞米松处理，
GEO: GSE52778），**自带现成计数矩阵**，不需要下载 FASTQ、不需要比对。

`airway` 的真实结构是 **4 个细胞系 × 2 种处理 = 8 个样本**，
即 **4 对配对样本**，不是两个独立组各 4 个。本项目在 `data/design.csv` 里
用 `group` 列把它当作 4 vs 4 处理（**不控制细胞系**），理由与代价写在
[REPORT.md](REPORT.md#三样本方案为什么是-4v4-而不控制细胞系) 里。

由于计数矩阵与官方基准只能在装有 R + DESeq2 的环境里产出，
这一步用一个独立脚本完成（详见 [laptop/deseq2_baseline.R](laptop/deseq2_baseline.R)）：

```bash
# 在装有 R 4.6 + DESeq2 + airway 的机器上
mkdir -p ~/deseq2-baseline
Rscript laptop/deseq2_baseline.R 2>&1 | tee ~/deseq2-baseline/run.log
```

它产出的文件与用途：

| 文件 | 用途 |
|---|---|
| `counts.csv` | 计数矩阵（63677 × 8），**两条路径共用的唯一输入** |
| `design.csv` | 分组表（`sample,condition,cell,group`） |
| `deseq2_result.csv` | 官方结果（默认设置，含 independent filtering） |
| `deseq2_result_noFilter.csv` | 官方结果（关闭过滤）：**对拍用这个**，因为我们的实现不做这些过滤 |
| `deseq2_result_paired.csv` | `~cell+dex` 版本（配对对照，用于说明控制细胞系的代价） |
| `sizeFactors.csv` | 官方 size factor，**M1 的直接对拍目标** |
| `dispersion_intermediates.csv` | 官方离散度的三个中间量，**M2 分段对拍的目标** |

把 `counts.csv`、`design.csv` 放到 `data/`，其余放入 `data/reference/`。

另外还准备了一个**第二数据集**用于跨物种复核（小鼠 DSS 模型，day00 vs day07，3v3）：
`data/nbis/` 下是 counts / design / 第三方发布的官方结果，以及转换与对拍脚本
（`convert.py` / `compare.py`）。它不需要 R 环境，数据文件都从
[NBISweden/workshop-RNAseq](https://github.com/NBISweden/workshop-RNAseq) 直接下载：

```powershell
# 两个文件下载到 data/nbis/ 后，依次运行 convert.py → deg.py → compare.py
# https://raw.githubusercontent.com/NBISweden/workshop-RNAseq/master/data/gene_counts.csv
# https://raw.githubusercontent.com/NBISweden/workshop-RNAseq/master/data/dge_results.csv
```

对拍结果（log2FC r = 0.999999 等）见 [REPORT.md 8.11 节](REPORT.md)。

---

## 五、目录结构

```
DESeq2-FromScratch/
├── README.md                  本文件
├── 功能介绍.md                 给零基础用户的说明书（能干什么 + 怎么用）
├── REPORT.md                  对拍验证报告（核心产出）
├── deg.py                     命令行入口
├── demo.py                    一键演示（模拟数据 + 真值对照）
├── src/                       四个模块的实现 + 辅助
├── tests/
│   ├── test_deg.py            单元测试（50 项）
│   ├── run_acceptance.py      回归验收脚本（38 个验收点）
│   ├── compare_baseline.py    与官方基准逐段对拍（产出 REPORT.md 第八节的数字）
│   ├── verify_independent.py  独立复算（只用标准库另算一遍报告里的数字）
│   └── fault_injection.py     故障注入反向验证（证明验收点不是摆设）
├── data/
│   ├── counts.csv             计数矩阵（需从笔记本拷回）
│   ├── design.csv             分组表
│   ├── reference/             官方 DESeq2 的输出（对拍基准）
│   └── nbis/                  第二数据集（NBIS）：counts + 第三方官方结果 + 对拍脚本
├── laptop/
│   └── deseq2_baseline.R      在 R 环境里产出数据与基准
├── build/demo/                demo.py 的产物（可重新生成）
└── logs/                      每次验收的记录（JSON）
```

---

## 六、我做出来的实现与 DESeq2 的已知差异

差异清单是理解对拍数字的前提。完整讨论见 [REPORT.md](REPORT.md#五我的实现与-deseq2-的差异清单)。

| 差异 | 影响 | 处理方式 |
|---|---|---|
| 不做 **independent filtering** | 官方的 `padj` 会有一批 NA（低表达基因被过滤掉），我们没有 | 对拍时用官方的 `_noFilter` 版本，两边同一口径；过滤的代价在报告里单独量化 |
| 不做 **Cook's 距离离群替换** | 官方会把有离群样本的基因 p 值设为 NA，我们不会 | 同上 |
| 不做 **log2FC 收缩**（`lfcShrink`） | 官方 `results()` 默认也是 MLE，两边一致；但我们也不做 apeglm | 一致，无需处理 |
| 参数化趋势拟合失败时的降级方案不同 | 官方用 locfit 做局部回归；我们用「分箱取中位数 + log-log 插值」 | 保留第三方库约束的代价，报告里声明 |
| 优化器实现不同 | 官方 MLE 用 C++ 线搜索；我们用「粗网格 + 黄金分割」 | 目标函数相同，解应一致；差异体现在收敛精度上 |
| 不做 t 近似 | 官方默认也是标准正态（`2*pnorm`） | 一致 |

---

## 七、验收标准（分级）

脚本化的验收标准见 [tests/run_acceptance.py](tests/run_acceptance.py)。
按三级验收标准（及格 / 良好 / 优秀）：

| 等级 | 标准 | 当前状态 |
|---|---|---|
| 及格 | M1/M3/M4 正确 + 测试通过 + 真实数据能跑 + log2FC 相关系数 > 0.8 + 报告诚实 | 达标：log2FC r = **0.9997**（airway）/ **0.999999**（NBIS 第二数据集） |
| 良好 | 以上 + M2 趋势拟合 + 相关系数 > 0.9 + 能定量说明改善 | 达标：含 M2 完整收缩与消融实验（REPORT 8.9） |
| 优秀 | 以上 + 完整经验贝叶斯收缩 + Jaccard > 0.7 + 逐类分析不一致基因 + 性能对比表 | 全部达标：Jaccard = **0.9017**、分层归因与性能表齐全、含 Git 提交历史（详见 REPORT 第十节） |

---

## 八、许可与声明

本项目是**学习/验证性质的对拍实现**，用于理解 DESeq2 的算法。
做真实的差异表达分析请使用官方 DESeq2：
<https://bioconductor.org/packages/release/bioc/html/DESeq2.html>

算法依据：Love MI, Huber W, Anders S. *Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2.* Genome Biology 15:550 (2014).
DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)
