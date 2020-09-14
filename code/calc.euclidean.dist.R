# Load R libs ####
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

# Load data #### 
# The "Stoeckius.result.RData" file contains
# 1. Stoeckius.counts (38k barcodes, unfiltered)
# 2. singlet (Seurat object, 13k singlets in the original paper)
# 3. neg.res.5ite.30splits (negative learning results)
load(file="../data/Stoeckius.result.RData")

# recovered ####
keep.vector <- neg.res.5ite.30splits$nn.keep
names(keep.vector) <- rownames(Stoeckius.counts)
n_counts <- Matrix::rowSums(Stoeckius.counts)
recover.bcs <- intersect(names(n_counts[n_counts<200]),
                         names(keep.vector[keep.vector=="TRUE"]))
recover <- CreateSeuratObject(counts = t(Stoeckius.counts[recover.bcs,]))
recover <- NormalizeData(recover)
genes.use <- VariableFeatures(singlet)
recover <- ScaleData(recover, features = genes.use)
recover <- RunPCA(recover,features=genes.use)
recover <- FindNeighbors(recover, dims = 1:10)
recover <- FindClusters(recover, resolution = 0.2)
recover <- RunTSNE(recover, dims = 1:10, check_duplicates = FALSE)
new.cluster <- c("recovered CD14 Mono","recovered CD4","recovered NK","recovered B","recovered platelet")
names(new.cluster) <- levels(recover)
recover <- RenameIdents(recover,new.cluster)
DimPlot(recover, label=TRUE)+labs(title="EmptyNN recovered 885 cells")+NoLegend()

# Project recovered bcs into space defined by singlets ####
sub.genes.use <- rownames(recover[["RNA"]]@scale.data)
prj.embed <- t(as.matrix(recover[["RNA"]]@scale.data[sub.genes.use,])) %*% 
        Loadings(singlet[["pca"]])[sub.genes.use,]

# combined seurat object 
combined <- merge(singlet, y = recover)
combined <- ScaleData(combined,features=genes.use)
combined <- RunPCA(combined,features=genes.use)
combined@reductions$pca@cell.embeddings <- rbind(singlet@reductions$pca@cell.embeddings,prj.embed)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:10)
combined <- FindClusters(combined, resolution = 0.35)
combined <- RunTSNE(combined, reduction = "pca", 
                    dims = 1:10,check_duplicates = FALSE)
combined$source <- c(rep("original",ncol(singlet)),rep("recovered",ncol(recover)))
combined$source2 <- c(as.character(Idents(singlet)),as.character(Idents(recover)))

# Calculate euclidean distance ####
require(graphics)
dist.mtx <- as.matrix(dist(combined@reductions$pca@cell.embeddings[,1:10],
                           method = "euclidean",upper=TRUE))

# Define function for downsampling count matrix ####
Down_Sample_Matrix <-function (expr_mat) {
  min_lib_size <- min(colSums(expr_mat))
  down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(unlist(lapply(x, function(y) {
      rbinom(1, y, prob)
    })))
  }
  down_sampled_mat <- apply(expr_mat, 2, down_sample)
  return(down_sampled_mat)
}

# Distance comparison ####
# original with recover
row <- names(combined$source2[combined$source2=="recovered B"])
column <- names(combined$source[combined$source=='original'])
df <- dist.mtx[row,column]
dist1 <- data.frame("dist"=colMeans(df))
idx <- match(rownames(dist1),names(combined$source2))
dist1$label <- combined$source2[idx]
ggplot(dist1, aes(x=label, y=dist)) + geom_boxplot(outlier.shape=NA)

# original with original
column <- names(combined$source2[combined$source2=='B cells'])
df <- dist.mtx[column,column]
dist2 <- as.vector(df)
dist2 <- data.frame(dist2,"B cells")
colnames(dist2) <- c("dist",'label')

# recover with recover
column <- names(combined$source2[combined$source2=='recovered B'])
df <- dist.mtx[column,column]
dist3 <- as.vector(df)
dist3 <- data.frame(dist3,"B cells")
colnames(dist3) <- c("dist",'label')

# boxplot
tmp <- rbind(dist1,dist2,dist3)
tmp$method <- c(rep("recovered with original",nrow(dist1)),
                rep("original with original",nrow(dist2)),
                rep("recover with recover",nrow(dist3)))
ggplot(tmp, aes(x=label, y=dist,fill=method)) + 
        geom_boxplot(outlier.shape=NA)+
        xlab("original clusters")+ylab("Euclidean distance")


# Extract distances for each cell type ####
tmp <- dist.mtx[which(combined@meta.data$source == 'original' &
                        combined@meta.data$source2 == 'B cells'),]

meta <- combined@meta.data
asplit <- split(rownames(meta), meta$source2)
cells <- unlist(lapply(asplit, function(x) sample(x, 10)))

dist_tmp <- dist.mtx[which(combined@meta.data$source == 'original' &
                             combined@meta.data$source2 == 'B cells'), cells]

aframe <- data.frame(combined@meta.data[cells, c('source', 'source2')], t(dist_tmp))

tmp <- melt(aframe, measure.vars = rownames(dist_tmp))
ggplot(tmp, aes(source2, value, fill = source2)) +
  facet_wrap(~ source) +
  geom_boxplot() + theme_bw()

# ####
# Run Seurat on downsampled count object ####
counts <- data.matrix(combined@assays$RNA@counts)
counts <- counts[rowSums(counts) > 0, ]
counts_ds <- Down_Sample_Matrix(counts)

seurat_ds <- CreateSeuratObject(counts = counts_ds)
seurat_ds <- NormalizeData(seurat_ds)
seurat_ds@meta.data <- combined@meta.data
genes.use <- VariableFeatures(singlet)
seurat_ds <- ScaleData(seurat_ds, features = genes.use)
seurat_ds <- RunPCA(seurat_ds,features=genes.use)
seurat_ds <- FindNeighbors(seurat_ds, dims = 1:10)
seurat_ds <- FindClusters(seurat_ds, resolution = 0.2)
seurat_ds <- RunTSNE(seurat_ds, dims = 1:10, check_duplicates = FALSE)

DimPlot(recover, label=TRUE)+labs(title="EmptyNN recovered 885 cells")+NoLegend()



