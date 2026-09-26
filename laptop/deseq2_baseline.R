# ============================================================================
# DESeq2-FromScratch：旧笔记本侧：导出数据 + 产出对拍基准
# ============================================================================
#
# 【这个脚本干什么】
#   在装有 R 4.6 + DESeq2 + airway 的旧笔记本上运行，产出四类东西：
#     1. 计数矩阵与分组表（主力机上那份纯 Python 实现要吃的唯一输入）
#     2. 官方 DESeq2 的差异表达结果（对拍基准）
#     3. M1/M2 的中间量（size factor、离散度），用于**逐段**验证我们的实现
#     4. 运行耗时与内存峰值
#
# 【怎么运行】
#   cd ~
#   Rscript deseq2_baseline.R 2>&1 | tee ~/deseq2-baseline/run.log
#
#   注意：输出目录由脚本自己创建，但 tee 先要目录存在，所以推荐：
#   mkdir -p ~/deseq2-baseline && Rscript ~/deseq2_baseline.R \
#       2>&1 | tee ~/deseq2-baseline/run.log
#
# 【产出文件清单】（全部在 ~/deseq2-baseline/ 下）
#   counts.csv                     计数矩阵（63677 基因 x 8 样本，未预过滤）
#   design.csv                     分组表（sample, condition, cell, group）
#   sample_meta.csv                完整样本元信息（留档）
#   deseq2_result.csv              ~dex 基准，默认设置（含 independent filtering）
#   deseq2_result_noFilter.csv     ~dex 基准，关闭 independentFiltering 与
#                                  Cook's cutoff：**对拍用这个**
#   deseq2_result_paired.csv       ~cell+dex 基准（配对对照，加分材料）
#   sizeFactors.csv                官方 size factor（M1 的直接答案）
#   dispersion_intermediates.csv   dispGeneEst / dispFit / dispersion（M2 三段答案）
#   run.log                        全部终端输出（含耗时、内存、sessionInfo）
#
# 【两个关键约定，改动会毁掉对拍】
#   ① 方向：全部用 contrast = c("dex","trt","untrt")，
#      即 log2FC > 0 表示「处理组 trt 比对照组 untrt 更高」。
#      DESeq2 默认是 untrt vs trt（反的），不写死这一行相关系数会变成 -1。
#   ② 基因集合：不做任何预过滤，全部 63677 个基因都导出。
#      我们的实现同样不丢基因（全 0 基因输出 NA），这样两边才能逐行对齐。
#
# 【预计耗时】（2026-09-26 实测）
#   基准 A（~dex）：15.4 秒。8 个样本的规模下，DESeq2 的耗时由「基因数 × 样本数」
#   共同决定，基因多但样本只有 8 个并不会爆炸，原先估的「十几分钟」高估了
#   约两个数量级。
#   基准 B（~cell+dex）：系数从 2 个变成 5 个，预计 1~2 分钟。
# ============================================================================

outdir <- file.path(path.expand("~"), "deseq2-baseline")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# 记录峰值内存（VmHWM = peak resident set size）
peak_mem_kb <- function() {
  x <- grep("VmHWM", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (length(x) == 0) return(NA_integer_)
  as.integer(sub(".*?([0-9]+).*", "\\1", x[1]))
}

t_script <- Sys.time()
cat("##############################################################\n")
cat("DESeq2-FromScratch 基准生成\n")
cat("开始时间:", format(t_script, "%Y-%m-%d %H:%M:%S"), "\n")
cat("输出目录:", outdir, "\n")
cat("##############################################################\n\n")

suppressMessages({
  library(airway)
  library(DESeq2)
})
data(airway)

cat("--- 环境 ---\n")
cat("R 版本        :", R.version.string, "\n")
cat("DESeq2 版本   :", as.character(packageVersion("DESeq2")), "\n")
cat("airway 版本   :", as.character(packageVersion("airway")), "\n\n")

# ---------------------------------------------------------------------------
# 1. 固定样本顺序并导出数据
# ---------------------------------------------------------------------------
# 排序规则：先按细胞系、再按处理。固定顺序是为了让 counts.csv 与 design.csv
# 的列一一对应，且每一次导出的结果完全一致（可复现）。
# WARNING: 关键：因子水平必须设在 airway 对象「本身」上，不能只设在本地副本 cd 上。
#    否则 DESeqDataSet() 读到的仍是原始的 dex，R 会按字母序把 "trt" 当成参照水平
#    （"t" < "u"），resultsNames 会变成 "dex_untrt_vs_trt"：方向正好反了，
#    与我们 design.csv 里 group=1 表示 trt 的约定冲突。
airway$dex  <- factor(airway$dex, levels = c("untrt", "trt"))
airway$cell <- factor(airway$cell)
cd <- as.data.frame(colData(airway))
ord <- order(cd$cell, cd$dex)
airway <- airway[, ord]
cd <- cd[ord, , drop = FALSE]

cnt <- assay(airway, "counts")

cat("--- 数据自检 ---\n")
cat("对象类型        :", class(airway), "\n")
cat("维度(基因x样本) :", paste(dim(airway), collapse = " x "), "\n")
cat("存储类型        :", typeof(cnt), "\n")
cat("NA 个数         :", sum(is.na(cnt)), "\n")
cat("负值个数        :", sum(cnt < 0), "\n")
cat("是否全为非负整数:", all(cnt >= 0) & all(cnt == floor(cnt)), "\n")
cat("全 0 基因个数   :", sum(rowSums(cnt) == 0), "\n")
cat("样本顺序        :", paste(colnames(cnt), collapse = ", "), "\n")
cat("分组(trt/untrt) :", paste(as.character(cd$dex), collapse = ", "), "\n\n")

stopifnot(sum(is.na(cnt)) == 0)
stopifnot(all(cnt >= 0))

# 计数矩阵：第一列表头留空（与 read.csv 的惯例一致），行名是基因 ID
write.csv(cnt, file.path(outdir, "counts.csv"), quote = FALSE)

design <- data.frame(
  sample    = colnames(cnt),
  condition = as.character(cd$dex),
  cell      = as.character(cd$cell),
  group     = as.integer(cd$dex) - 1L,     # untrt -> 0, trt -> 1
  stringsAsFactors = FALSE
)
write.csv(design, file.path(outdir, "design.csv"),
          row.names = FALSE, quote = FALSE)
write.csv(cd, file.path(outdir, "sample_meta.csv"), quote = FALSE)

cat("已导出 counts.csv / design.csv / sample_meta.csv\n")
cat("design.csv 内容:\n")
print(design)
cat("\n")

# ---------------------------------------------------------------------------
# 2. 基准 A：~dex（**对拍主基准**）
# ---------------------------------------------------------------------------
cat("##############################################################\n")
cat("基准 A：design = ~ dex （全部 8 个样本，4 vs 4）\n")
cat("##############################################################\n")

dds <- DESeqDataSet(airway, design = ~ dex)

t0 <- Sys.time()
dds <- DESeq(dds, quiet = FALSE)
sec_a <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\n>>> DESeq(~dex) 全量耗时: %.1f 秒 (%.1f 分钟)\n\n",
            sec_a, sec_a / 60))
cat(">>> 结果表可用系数名 (resultsNames):\n")
print(resultsNames(dds))
cat("\n")

# 统一方向：trt 相对于 untrt（log2FC > 0 表示处理组更高）
# 说明：这里刻意不再用「系数名」做断言，因子水平设置一变，系数名就会变
# （"dex_trt_vs_untrt" 或 "dex_untrt_vs_trt"），断死名字只会误拦。
# 方向是否正确，改由下面的「方向自检」用数据本身来验证，更可靠。

res_default <- results(dds, contrast = c("dex", "trt", "untrt"))
res_nofilter <- results(dds, contrast = c("dex", "trt", "untrt"),
                        independentFiltering = FALSE,
                        cooksCutoff = FALSE,
                        alpha = 0.05)

# ---- 方向自检（重要，不要删）----
# 用「trt 组归一化计数均值 / untrt 组归一化计数均值」手工算一个对数比值，
# 再和 DESeq2 的 log2FoldChange 求相关：
#     接近 +1 -> 方向正确
#     接近 -1 -> 方向反了（这是最致命、又最难肉眼发现的错误）
# 上一次运行就是因为方向问题（因子水平没设在 airway 对象上）触发了中断，
# 这一步以后能自动抓住同类问题。
.norm_cnt <- counts(dds, normalized = TRUE)
.trt_cols <- dds$dex == "trt"
.untrt_cols <- dds$dex == "untrt"
.manual_lfc <- log2((rowMeans(.norm_cnt[, .trt_cols, drop = FALSE]) + 0.5) /
                    (rowMeans(.norm_cnt[, .untrt_cols, drop = FALSE]) + 0.5))
.ok_dir <- is.finite(.manual_lfc) & is.finite(res_default$log2FoldChange)
r_dir <- cor(.manual_lfc[.ok_dir], res_default$log2FoldChange[.ok_dir])
cat(sprintf(">>> 方向自检：手工组均值比 vs DESeq2 log2FC 相关系数 = %.4f\n", r_dir))
cat("    （接近 +1 表示方向正确；接近 -1 表示方向反了）\n\n")
if (!is.finite(r_dir) || r_dir < 0) {
  stop("方向自检失败：log2FC 与「trt/untrt 组均值比」方向不一致，请检查 dex 的因子水平")
}
rm(.norm_cnt, .trt_cols, .untrt_cols, .manual_lfc, .ok_dir)

cat("--- 默认设置的结果概览 (summary) ---\n")
print(summary(res_default))
cat("\n--- 关闭 independentFiltering / Cook's cutoff 的概览 ---\n")
print(summary(res_nofilter))
cat("\n")

as_df <- function(r) {
  d <- as.data.frame(r)
  data.frame(gene_id = rownames(d), d, check.names = FALSE)
}

write.csv(as_df(res_default), file.path(outdir, "deseq2_result.csv"),
          quote = FALSE)
write.csv(as_df(res_nofilter), file.path(outdir, "deseq2_result_noFilter.csv"),
          quote = FALSE)

# --- M1 的直接答案：size factor ---
sf <- sizeFactors(dds)
write.csv(data.frame(sample = names(sf), sizeFactor = as.numeric(sf)),
          file.path(outdir, "sizeFactors.csv"), row.names = FALSE, quote = FALSE)
cat("--- size factors（这是 M1 的直接对拍目标）---\n")
print(round(sf, 6))
cat("\n")

# --- M2 的三段答案：gene-wise / 趋势 / 收缩 ---
md <- as.data.frame(mcols(dds))
md <- data.frame(gene_id = rownames(md), md, check.names = FALSE)
write.csv(md, file.path(outdir, "dispersion_intermediates.csv"), quote = FALSE)
cat("--- dispersion_intermediates.csv 的列 ---\n")
print(colnames(md))
cat("\n各离散度列的摘要:\n")
for (col in c("dispGeneEst", "dispFit", "dispersion", "dispMAP")) {
  if (col %in% colnames(md)) {
    cat(sprintf("  %-12s 中位数=%.6f  最小值=%.6f  最大值=%.6f\n",
                col, median(md[[col]], na.rm = TRUE),
                min(md[[col]], na.rm = TRUE), max(md[[col]], na.rm = TRUE)))
  }
}
cat("\n")

# 显著基因数（多组阈值，便于和我们自己的结果对照）
cat("--- 显著基因数（基于 res_default）---\n")
for (padj_cut in c(0.05, 0.01)) {
  for (lfc_cut in c(0, 1, 1.5)) {
    n <- sum(!is.na(res_default$padj) & res_default$padj < padj_cut &
             abs(res_default$log2FoldChange) >= lfc_cut)
    cat(sprintf("  padj<%.2f 且 |log2FC|>=%.1f : %d 个\n",
                padj_cut, lfc_cut, n))
  }
}
cat("\n")

rm(dds)
# 注意：R 的 gc() 没有 quiet 参数，只有 verbose / reset / full。
# 写成 gc(quiet = TRUE) 会直接报「参数没有用(quiet = TRUE)」并停止执行。
invisible(gc())

# ---------------------------------------------------------------------------
# 3. 基准 B：~cell + dex（配对对照，加分材料，只用于报告对比）
# ---------------------------------------------------------------------------
cat("##############################################################\n")
cat("基准 B：design = ~ cell + dex （统计上更正确，但我们的实现复现不了）\n")
cat("##############################################################\n")

# 这一节只是「加分材料」，而基准 A 的全部对拍产物在它之前就已经落盘。
# 所以整段用 tryCatch 兜住：万一这里再出问题，也绝不能让它把后面的汇总
# 一起带走（前两次失败都是「主产物写完了、脚本却在后面倒下」，代价太大）。
sec_b <- NA_real_
r_lfc <- NA_real_
.tb_err <- NULL

tryCatch({

  dds2 <- DESeqDataSet(airway, design = ~ cell + dex)
  t0 <- Sys.time()
  dds2 <- DESeq(dds2, quiet = TRUE)
  sec_b <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf(">>> DESeq(~cell+dex) 全量耗时: %.1f 秒 (%.1f 分钟)\n\n",
              sec_b, sec_b / 60))

  res_paired <- results(dds2, contrast = c("dex", "trt", "untrt"))
  write.csv(as_df(res_paired), file.path(outdir, "deseq2_result_paired.csv"),
            quote = FALSE)

  cat("--- 配对版显著基因数 ---\n")
  for (lfc_cut in c(0, 1, 1.5)) {
    n <- sum(!is.na(res_paired$padj) & res_paired$padj < 0.05 &
             abs(res_paired$log2FoldChange) >= lfc_cut)
    cat(sprintf("  padj<0.05 且 |log2FC|>=%.1f : %d 个\n", lfc_cut, n))
  }
  cat("\n")

  # 两份基准的一致程度（这是报告里的现成素材）
  # WARNING: 先转成普通 data.frame 再取列，不要对 DESeqResults 直接做
  #    res[行, "列"] 这种二维取子集：DESeq2 改写过它的 `[` 方法，
  #    返回类型不保证是数值向量，cor() 有可能直接报错。
  #    改成按 gene_id 用 merge 对齐，全程只碰普通向量，行为完全确定。
  df_a <- as_df(res_default)[, c("gene_id", "log2FoldChange")]
  df_b <- as_df(res_paired)[, c("gene_id", "log2FoldChange")]
  mm <- merge(df_a, df_b, by = "gene_id", suffixes = c("_a", "_b"))
  ok_m <- is.finite(mm$log2FoldChange_a) & is.finite(mm$log2FoldChange_b)
  r_lfc <- cor(mm$log2FoldChange_a[ok_m], mm$log2FoldChange_b[ok_m],
               method = "pearson")
  cat(sprintf("两份基准的 log2FC 相关系数 = %.4f （%d 个基因参与计算）\n",
              r_lfc, sum(ok_m)))
  cat("（这个数字本身就有意义：它说明「控制细胞系」对结果的影响有多大）\n")
  # 这个相关系数同时充当「配对版方向自检」：只要有一份结果的方向反了，它就会变成负数
  if (!is.finite(r_lfc) || r_lfc < 0) {
    cat("!! 警告：两份基准的 log2FC 相关系数为负或不可计算，其中一份方向可能反了\n")
  }
  cat("\n")

  rm(dds2, res_paired, df_a, df_b, mm, ok_m)
  invisible(gc())

}, error = function(e) {
  .tb_err <<- conditionMessage(e)
  cat("\n!! 基准 B（~cell+dex 配对版）失败，已跳过。\n")
  cat("   错误原文：", .tb_err, "\n")
  cat("   基准 A 的产物不受影响，可以照常拿回去做对拍。\n\n")
})
invisible(gc())

# ---------------------------------------------------------------------------
# 4. 收尾
# ---------------------------------------------------------------------------
total_sec <- as.numeric(difftime(Sys.time(), t_script, units = "secs"))
cat("##############################################################\n")
cat("汇总\n")
cat("##############################################################\n")
cat(sprintf("DESeq(~dex) 耗时      : %.1f 秒\n", sec_a))
if (is.finite(sec_b)) {
  cat(sprintf("DESeq(~cell+dex) 耗时 : %.1f 秒\n", sec_b))
  cat(sprintf("两份基准 log2FC 相关  : %.4f\n", r_lfc))
} else {
  cat("DESeq(~cell+dex) 耗时 : 已跳过（基准 B 未成功，错误原文见上方）\n")
}
cat(sprintf("整个脚本耗时          : %.1f 秒 (%.1f 分钟)\n",
            total_sec, total_sec / 60))
cat(sprintf("进程峰值内存 (VmHWM)  : %.1f MB\n", peak_mem_kb() / 1024))
cat("\n产出文件:\n")
for (f in sort(list.files(outdir, full.names = FALSE))) {
  cat(sprintf("  %-32s %10.1f KB\n", f,
              file.info(file.path(outdir, f))$size / 1024))
}
cat("\n会话信息 (sessionInfo):\n")
print(sessionInfo())
cat("\n>>> 完成。请把整个 ~/deseq2-baseline 目录拷回主力机的 data/reference/ 下。\n")
