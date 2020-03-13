setwd("/Users/yanfangfang/Downloads/Empty_Droplet/data")
# Load R libs ####
library(Seurat)
library(Matrix)
library(keras)
library(dplyr)

# Load count table ####
counts <- readMM(file = "Stoeckius_RNA.mtx")
barcodes <- read.csv("Stoeckius_barcodes.txt")[,2]
genes <- read.csv("Stoeckius_genes.txt")[,2]
rownames(counts) <- barcodes
colnames(counts) <- genes

# keep barcodes that have HTO counts ####
jointbcs <- read.csv("joint_barcodes.txt",header=F)[,1]
counts <- counts[rownames(counts) %in% jointbcs,]

# Calculate total reads ####
n_counts <- Matrix::rowSums(counts)
names(n_counts) <- rownames(counts)

# Remove cells with less than 10 reads ####
counts <- counts[which(n_counts > 10), ]
n_counts <- n_counts[which(n_counts > 10)]

# Visualize distribution of total counts ####
plot(n_counts, log = "y", ylab = "Total counts", xlab = "Sorted barcodes")
abline(h = 200, col ="red")

# Keep 5k most frequent genes ####
counts_5k <- counts[, names(tail(sort(colSums(counts)), 5000))]

# Load classification.global ####
global <- read.delim("classification.global.txt",sep="\t",header=T)
global <- global[rownames(global) %in% rownames(counts_5k),1,drop=F]

# Split barcodes into quantiles by total counts ####
above_200 <- which(n_counts >= 200)
below_200 <- which(n_counts < 200)
split_into_quantiles <- function(tmp){
  qs <- seq(0, 1, length = 20)
  split(names(tmp), cut(tmp, breaks=c(quantile(tmp, probs = qs))))
}
top_splits <- split_into_quantiles(n_counts[above_200])
top_splits <- top_splits[rev(1:length(top_splits))]
bottom_splits <- split_into_quantiles(n_counts[below_200])

# Training set and Testing set ####
test <- c(unlist(top_splits[c(18,19)]),unlist(bottom_splits[c(18,19)]))
train  <- setdiff(rownames(counts_5k), test)
x_train <- data.matrix(counts_5k[train,])
x_test <- data.matrix(counts_5k[test,])
y_train <- rep(0, nrow(x_train))
names(y_train) <- rownames(x_train)
y_train[intersect(names(which(n_counts > 200)), rownames(x_train))] <- 1
y_train <- to_categorical(y_train)
y_test <- global[test,1]
y_test[which(y_test=='Doublet')] <- "Singlet"
y_test <- droplevels(y_test)

# Create NN ####
create_model <- function(){
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 256, activation = 'relu', input_shape = c(ncol(counts_5k))) %>% 
    layer_dropout(rate = 0.4) %>% 
    layer_dense(units = 128, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 2, activation = 'softmax')
  model %>% compile(
    loss = 'binary_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  model
} 

# Non-curriculum learning ####
model_non_curri <- create_model()
history <- model_non_curri %>% fit(
  x_train, y_train,
  batch_size = 64,
  epochs = 17,
  #steps_per_epoch = 10,
  verbose = 2,
  validation_split = 0.1
)
y_test_pred <- model_non_curri %>% predict(x_test)
acc <- unclass(table(y_test_pred[,2] > 0.5, y_test))
non_curri_testing_acc <- sum(diag(acc))/sum(acc)

# Curriculum learning ####
model_curri <- create_model()
accs <- lapply(1:17, function(x){
  print(x)
  training_accs <- unlist(lapply(1:10, function(k){
    top_barcodes <- sample(top_splits[[x]], 64)
    bottom_barcodes <- sample(bottom_splits[[x]], 64)
    batch_indices <- c(match(top_barcodes, rownames(x_train)), match(bottom_barcodes, rownames(x_train)))
    y <- c(rep(1, 64), rep(0, 64))
    r <- sample(1:length(batch_indices))
    batch_indices <- batch_indices[r]
    y <- y[r]
    stats <- train_on_batch(model_curri, x_train[batch_indices,], 
                            to_categorical(y), class_weight = NULL, sample_weight = NULL)
    stats[[2]]
  }))
  
  training_acc <- mean(training_accs)
  
  test_top <- unlist(top_splits[c(18, 19)])
  x_test <- data.matrix(counts_5k[test_top,])
  preds_top <- model_curri %>% predict(x_test)
  
  test_bottom <- unlist(bottom_splits[c(18, 19)])
  x_test <- data.matrix(counts_5k[test_bottom,])
  preds_bottom <- model_curri %>% predict(x_test)
  
  outcome <- global[c(test_top,test_bottom),1]
  outcome[which(outcome=='Doublet')] <- "Singlet"
  outcome <- droplevels(outcome)
  pred_final <- c(preds_top[,2], preds_bottom[,2])
  
  testing_acc <- unclass(table(pred_final > 0.5, outcome))
  testing_acc <- sum(diag(testing_acc))/sum(testing_acc)
  return(c(training_acc, testing_acc))
})

# Plot training/testing accuracies ####
accs <- do.call(rbind, accs)
plot(accs[,1], ylim = c(0,1), col = "blue", pch = 16, 
     xlab="quantile",ylab="accuracy",
     main="Traning and Testing accuracy
     for each quantile of Curriculum learning")
points(accs[,2], col = "red", pch = 16)
legend("bottomright", c('training', "testing"), col = c("blue", "red"), pch = 16, bty = "n")

# Save model ####
#model %>% save_model_hdf5("model_0309.h5")
#model <- load_model_hdf5("model_0309.h5")

# Boxplot of accuracies ####
test <- c(unlist(top_splits[c(18, 19)]),unlist(bottom_splits[c(18, 19)]))
labels_tests <- global[test,1,drop=F]
x_test <- data.matrix(counts_5k[rownames(labels_tests),])
preds <- model %>% predict(x_test)
labels_tests$preds <- preds[,1]
boxplot(preds ~ classification.global, labels_tests,outline=F,
        xlab="True labels",ylab="Predicted probability")
fisher.test(table(labels_tests$preds > 0.5, labels_tests$classification.global))

# New calculated HTO tsne plot of barcodes in the cutoff area ####
htos <- read.csv(file = "GSM2895283_Hashtag-HTO-count_transpose.csv",sep = ",",header = TRUE, row.names = 1)
htos <- as.matrix(t(htos)[-c(9,10,11), test])
rownames(htos) <- unlist(lapply(rownames(htos), function(x) strsplit(x, ".", fixed = T )[[1]][1]))
rownames(htos) <- gsub("Batch", "HTO_", rownames(htos))
tmp <- CreateSeuratObject(counts = htos)
tmp <- NormalizeData(tmp, normalization.method = "CLR")
dist.mtx <- as.matrix(dist(t(GetAssayData(object = tmp))))
tmp <- RunTSNE(tmp, distance.matrix = dist.mtx, perplexity = 100)
DimPlot(tmp)
embed <- tmp@reductions$tsne@cell.embeddings
embed <- data.frame(embed,preds)
ggplot(embed,aes(tSNE_1,tSNE_2,color=Predicted_single))+geom_point()+
  scale_color_gradient(low="grey", high="aquamarine3")+ggtitle("barcodes in the cutoff area")
ggplot(embed,aes(tSNE_1,tSNE_2,color=Predicted_label))+geom_point()+ggtitle("barcodes in the cutoff area")

# Recovered 687 cells ####
recover <- rownames(subset(preds,rownames(preds) %in% unlist(bottom_splits[c(18, 19)]) & Predicted_label=="single"))
recover_counts <- t(counts[recover,])
# remove mitochondrial genes and ribosomal genes
mito.genes <- grep(pattern = "^MT-", rownames(recover_counts), ignore.case=TRUE, value = TRUE)
RPS.genes <- grep(pattern = "^RPS", rownames(recover_counts), ignore.case=TRUE, value = TRUE)
RPL.genes <- grep(pattern = "^RPL", rownames(recover_counts), ignore.case=TRUE, value = TRUE)
recover_counts <- recover_counts[!(rownames(recover_counts) %in% c(mito.genes,RPS.genes,RPL.genes)),]
recover_res <- CreateSeuratObject(counts = recover_counts)
recover_res <- NormalizeData(recover_res)
genes.use <- read.delim("genes.use.txt",header=F)[,1]
recover_res <- ScaleData(recover_res, features = genes.use)
recover_res <- RunPCA(recover_res, features = genes.use)
recover_res <- FindNeighbors(recover_res, dims = 1:10)
recover_res <- FindClusters(recover_res, resolution = 1)
recover_res <- RunTSNE(recover_res, dims = 1:10)
DimPlot(recover_res,label=TRUE)
markers <- FindAllMarkers(recover_res, only.pos = TRUE)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
genes <- c("IL7R","CD14", 'LYZ', "MS4A1","CD79A","CD79B","CD8A",
           "GNLY", "NKG7", "FCER1A", "CST3", "PPBP")
FeaturePlot(recover_res,features=genes)
new.cluster.ids <- c("CD14 Mono","Memory CD4","NK","B cells")
names(new.cluster.ids) <- levels(recover_res)
recover_res <- RenameIdents(recover_res, new.cluster.ids)
DimPlot(recover_res,label=TRUE)
#saveRDS(recover_res, file = "recovered.rds")

# Lost ####
lost <- rownames(subset(preds,rownames(preds) %in% unlist(top_splits[c(18, 19)]) & Predicted_label=="non-single"))
lost_counts <- t(counts[lost,])
lost_res <- CreateSeuratObject(counts = lost_counts)
lost_res <- NormalizeData(lost_res)
lost_res <- ScaleData(lost_res, features = genes.use)
lost_res <- RunPCA(lost_res, features = genes.use)
lost_res <- FindNeighbors(lost_res, dims = 1:10)
lost_res <- FindClusters(lost_res, resolution = 1)
lost_res <- RunTSNE(lost_res, dims = 1:10,check_duplicates = FALSE)
lost_res <- readRDS("loss.rds")
DimPlot(lost_res,label=TRUE)
genes <- c("GNLY","NKG7","LYZ","S100A9","MS4A1","CD79A")
FeaturePlot(lost_res,features=genes)
#saveRDS(lost_res, file = "loss.rds")


