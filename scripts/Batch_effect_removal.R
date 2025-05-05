### Batch effect removal

tpm_matrix <- read.csv("11_matrix/rsem_tpm_geneName.csv", header = T, row.names = 1)
head(tpm_matrix)
dim(tpm_matrix)
# 构建样本信息表
sample_metadata <- data.frame(
  sample_id = c("SRR23558702", "SRR23558703", "SRR23558705", "SRR23558706"),
  group = c("hypoxic", "hypoxic", "aerobic", "aerobic"),
  batch = c("rep2", "rep1", "rep2", "rep1")
)
rownames(sample_metadata) <- sample_metadata$sample_id
sample_metadata$batch <- factor(sample_metadata$batch)

# ----sva::ComBat()-------------------------------------------------------------
# step1:数据准备
tpm_log <- log2(tpm_matrix + 1)
# 检查并去掉方差低的基因（保留方差最高的 90%）
vars <- apply(tpm_log, 1, var)
threshold <- quantile(vars, 0.1)
tpm_log_filtered <- tpm_log[vars > threshold, ]
# step2:批次效应校正前 PCA
library(PCAtools)
# 注意：PCAtools 的 pca() 输入是样本 × 基因，需要转置
p_before <- pca(tpm_log_filtered, metadata = sample_metadata, removeVar = 0.1, scale = TRUE)
biplot(p_before, colby = 'group', shape = 'batch', legendPosition = 'right') + ggtitle('PCA Before Batch Correction')
# step3:使用 ComBat 去除批次效应
library(sva)
# 构建设计矩阵
mod <- model.matrix(~ group, data = sample_metadata)
# ComBat 批次校正
combat_data <- ComBat(dat = as.matrix(tpm_log_filtered),
                      batch = sample_metadata$batch,
                      mod = mod,
                      par.prior = TRUE,
                      prior.plots = FALSE)
# 校正后 PCA
p_after <- pca(combat_data, metadata = sample_metadata, removeVar = 0.1, scale = TRUE)
biplot(p_after, colby = 'group', shape = 'batch', legendPosition = 'right') + ggtitle('PCA After Batch Correction')
