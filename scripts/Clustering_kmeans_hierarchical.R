### Clustering
# 表达矩阵行为基因，列为样本
# 对样本进行聚类

tpm_matrix <- read.csv("11_matrix/rsem_tpm_geneName.csv", header = T, row.names = 1)
head(tpm_matrix)
dim(tpm_matrix)

# ----K-means-------------------------------------------------------------------
# step1:log 转换
tpm_log <- log2(tpm_matrix + 1)
# step2:去除表达量全为 0 和低方差基因
tpm_log <- tpm_log[rowSums(tpm_log) > 0, ]
library(matrixStats)
gene_vars <- rowVars(as.matrix(tpm_log))
threshold <- quantile(gene_vars, 0.1)
tpm_filtered <- tpm_log[gene_vars > threshold, ]# 去掉低方差基因（例如前10%）
# step3:z-score
tpm_scaled <- t(scale(t(tpm_filtered)))

# step4:转置——行是样本，列是基因(对样本进行聚类)!!!
tpm_for_clustering <- t(tpm_scaled)  # 行是样本，列是基因

# step5:确定最佳 k 值(Elbow Method)
set.seed(123)
wss <- sapply(1:10, function(k){
  kmeans(tpm_for_clustering, centers = k, nstart = 25)$tot.withinss
})
# 看拐点（elbow） → 当增加 k 后，wss 不再大幅下降，就是合理的 k
plot(1:10, wss, type = 'b', pch = 19, frame = FALSE,
     xlab = "Number of Clusters K",
     ylab = "Total Within-Clusters Sum of Squares")
# step4:k-means 聚类
set.seed(123)
kmeans_result <- kmeans(tpm_for_clustering, centers = 3, nstart = 25) #确定合适的k

# 可视化
library(factoextra)
fviz_cluster(kmeans_result, data = tpm_for_clustering,
             palette = "jco", # 颜色方案
             ellipse.type = "euclid",
             star.plot = TRUE,
             repel = TRUE)

# ----hierarchical clustering---------------------------------------------------
# step1:log 转换
tpm_log <- log2(tpm_matrix + 1)
# step2:取方差最高的10%的基因差
gene_variances <- apply(tpm_log, 1, var)
num_top <- ceiling(0.1 * nrow(tpm_log))
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:num_top]
tpm_top <- tpm_log[top_genes, ]
# step3:转置--对样本聚类
tpm_top_t <- t(tpm_top)
# step4:计算样本间距离矩阵（默认欧氏距离，你也可以用 correlation 等）
dist_matrix <- dist(tpm_top_t, method = "euclidean")
# step5:分层聚类（"ward.D2" → 方差最小化（推荐用于表达数据），可以换成complete、average、single）
hc <- hclust(dist_matrix, method = "ward.D2")
# 可视化基础树状图
plot(hc, labels = colnames(tpm_matrix), main = "Hierarchical Clustering Dendrogram",
     xlab = "", sub = "", cex = 0.8)
# step6:切割树分组
clusters <- cutree(hc, k = 2)
print(clusters)
# 可视化着色树状图
library(dendextend)
dend <- as.dendrogram(hc)
dend_colored <- color_branches(dend, k = 2)
plot(dend_colored, main = "Colored Hierarchical Clustering Dendrogram")
