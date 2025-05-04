# 1. 加载需要的包
library(tximport)
library(readr)
library(DESeq2)

# 2. 设置文件路径与样本名
rsem_dir <- "/home/future/rnaseq_batch_250424/08_rsem"
samples <- c("SRR23558706", "SRR23558705", "SRR23558703", "SRR23558702")  # 根据你实际样本添加更多

# 3. 构建文件路径命名列表
files <- file.path(rsem_dir, paste0(samples, ".genes.results"))
names(files) <- samples

# 4. 构建样本信息表（colData）
sampleTable <- data.frame(
  sample = samples,
  condition = c("control", "control", "treated", "treated"),  # 请根据你的实验设计修改
  batch = c("batch1", "batch2", "batch1", "batch2"), # 请根据你的实验设计修改 （注意共线性问题）
  row.names = samples
)
sampleTable$condition <- factor(sampleTable$condition)
sampleTable$batch <- factor(sampleTable$batch)

# 检查共线性问题
table(sampleTable$condition, sampleTable$batch)

# 5. 使用 tximport 导入 RSEM 结果
txi <- tximport(files,
                type = "rsem",
                txIn = FALSE,       # RSEM是gene-level文件
                txOut = FALSE
)

# 检查并过滤掉 length == 0 的基因
keep <- rowSums(txi$length > 0) == ncol(txi$length)
txi_filtered <- list(
  counts = txi$counts[keep, ],
  abundance = txi$abundance[keep, ],
  length = txi$length[keep, ],
  countsFromAbundance = txi$countsFromAbundance  # 添加这个字段
)

# 6. 构建 DESeq2 对象
dds <- DESeqDataSetFromTximport(txi_filtered,
                                colData = sampleTable,
                                design = ~ batch + condition)

# 7. 查看对象
dds

# 8. 差异分析核心步骤
dds <- DESeq(dds)

# 9. 列出可提取的比较
resultsNames(dds)

# 10. 提取差异表达分析结果
res <- results(dds, name="condition_treated_vs_control")
#res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm") # 对 log2FoldChange 结果进行收缩

head(res)[1:5,]

### 查看结果
summary(res)
# 调整后的 p 值小于 0.1
sum(res$padj < 0.1, na.rm=TRUE)

# FDR 5%：padj 小于 0.05 的基因，就会被标记为显著差异表达
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

### 输出结果
# 1.保存完整结果
res_df <- as.data.frame(res)
# 根据padj升序排序（padj小的排在前面）
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file="DESeq2_all_results.csv", row.names=TRUE)

# 2.保存筛选后的显著差异表达基因
# 筛选 padj < 0.05 且 log2FC大于1或小于-1的基因
sig_res <- res05[
  which(res05$padj < 0.05 & abs(res05$log2FoldChange) > 1),
]
sig_res_df <- as.data.frame(sig_res)
sig_res_df <- sig_res_df[order(sig_res_df$padj), ]
write.csv(sig_res_df, file="DESeq2_significant_DEGs.csv", row.names=TRUE)

# 3.区分上调和下
# 上调基因
up_genes <- sig_res_df[sig_res_df$log2FoldChange > 1, ]
# 下调基因
down_genes <- sig_res_df[sig_res_df$log2FoldChange < -1, ]
write.csv(up_genes, file="DESeq2_upregulated_genes.csv", row.names=TRUE)
write.csv(down_genes, file="DESeq2_downregulated_genes.csv", row.names=TRUE)
