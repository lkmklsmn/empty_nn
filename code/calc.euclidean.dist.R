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
load(file="Stoeckius.result.RData")

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
combined@reductions$pca@cell.embeddings <- rbind(subset@reductions$pca@cell.embeddings,prj.embed)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:10)
combined <- FindClusters(combined, resolution = 0.35)
combined <- RunTSNE(combined, reduction = "pca", 
                    dims = 1:10,check_duplicates = FALSE)
combined$source <- c(rep("original",ncol(singlet)),rep("recovered",ncol(recover)))
combined$source2 <- c(as.character(Idents(singlet)),as.character(Idents(recover)))

# Calculate euclidean distance ####
require(graphics)
dist.mtx <- as.matrix(dist(combined@reductions$tsne@cell.embeddings,
                           method = "euclidean",upper=TRUE))

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
dist2 <- data.frame(dist2,"CD14 Mono")
colnames(dist2) <- c("dist",'label')

# recover with recover
column <- names(combined$source2[combined$source2=='recovered B'])
df <- dist.mtx[column,column]
dist3 <- as.vector(df)
dist3 <- data.frame(dist3,"CD14 Mono")
colnames(dist3) <- c("dist",'label')

# boxplot
tmp <- rbind(dist1,dist2,dist3)
tmp$method <- c(rep("recovered with original",nrow(dist1)),
                rep("original with original",nrow(dist2)),
                rep("recover with recover",nrow(dist3)))
ggplot(tmp, aes(x=label, y=dist,fill=method)) + 
        geom_boxplot(outlier.shape=NA)+
        xlab("original clusters")+ylab("Euclidean distance")


