### PCA
# 用 RSEM 得到的 TPM 矩阵做 PCA 分析

tpm_matrix <- read.csv("11_matrix/rsem_tpm_geneName.csv", header = T, row.names = 1)
head(tpm_matrix)
dim(tpm_matrix)
# 构建样本信息表
sample_metadata <- data.frame(
  sample_id = c("SRR23558702", "SRR23558703", "SRR23558705", "SRR23558706"),
  group = c("hypoxic", "hypoxic", "aerobic", "aerobic"),
  replicate = c("rep2", "rep1", "rep2", "rep1")
)
rownames(sample_metadata) <- sample_metadata$sample_id

# ----prcomp() 函数-------------------------------------------------------------
# step1:log 转换
tpm_log <- log2(tpm_matrix + 1)
# step2:z-score(中心化+标准化)
tpm_scaled <- t(scale(t(tpm_log), center = TRUE, scale = TRUE))
# step3:去除 NA
tpm_scaled_clean <- tpm_scaled[complete.cases(tpm_scaled), ]
# step4:运行 PCA
pca_res <- prcomp(t(tpm_scaled_clean), center = FALSE, scale. = FALSE)

# 结果
summary(pca_res) # 主成分贡献度
head(pca_res$x) # 各样本在 PC 空间的坐标

# 可视化
# 1.快速简易散点图
plot(pca_res$x[,1], pca_res$x[,2],
     xlab = "PC1", ylab = "PC2",
     main = "PCA of RNA-seq (TPM)")
# 2.ggplot2散点图
library(ggplot2)
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text(vjust = -0.5) +
  theme_minimal() +
  ggtitle("PCA plot of RNA-seq samples")
# 3.factoextra包可视化
library(factoextra)
fviz_pca_ind(pca_res,
             geom.ind = "point",
             col.ind = sample_metadata$group,
             addEllipses = TRUE,
             legend.title = "Group")
fviz_screeplot(pca_res, addlabels = TRUE) # scree plot


# ----PCAtools包----------------------------------------------------------------
library(PCAtools)
# step1:log 转换
tpm_log <- log2(tpm_matrix + 1)
# step2:去掉方差为 0 的基因
gene_variances <- apply(tpm_log, 1, var)
sum(gene_variances == 0)
tpm_log_clean <- tpm_log[gene_variances != 0, ]
# 运行PCA(去掉前 10% 低方差基因（可调），并进行 z-score 标准化)
p <- pca(tpm_log_clean, metadata = sample_metadata, removeVar = 0.1, scale = TRUE)

#可视化
screeplot(p) # 主成分方差解释比例
biplot(p, colby = 'group', shape = 'replicate', legendPosition = 'right') # 二维 PCA 图
