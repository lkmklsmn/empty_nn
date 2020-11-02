# Load R libs ####
library(EmptyNN)
library(Seurat)
library(Matrix)
library(keras)
library(DropletUtils)

# Load data ####
counts <- Read10X_h5("./../data/pbmc_10k_raw.h5",use.names = TRUE, unique.features = TRUE)

# EmptyNN ####
nn.res <- emptynn(t(counts),threshold=100,k=5,iteration=5,verbose=TRUE)
nn.keep <- nn.res$nn.keep

# CellRanger 2.0 ####
ranger.keep <- defaultDrops(counts, expected=10000)

# EmptyDrops ####
e.out <- emptyDrops(counts)
e.keep <- e.out$FDR <= 0.001

# CellBender ####
bender <- Read10X_h5("./../data/pbmc.10k.CellBender_filtered.h5")
bender.keep <- colnames(counts) %in% colnames(bender)
names(bender.keep) <- colnames(counts)

# Comparison ####
keep <- nn.keep | ranger.keep | (e.keep & !is.na(e.keep)) | bender.keep
retained <- runSeurat(counts=counts[,keep],resolution=0.1)
retained$emptynn <- colnames(retained) %in% colnames(counts)[nn.keep]
retained$cellranger <- colnames(retained) %in% colnames(counts)[ranger.keep]
retained$emptydrops <- colnames(retained) %in% colnames(counts)[e.keep & !is.na(e.keep)]
retained$cellbender <- colnames(retained) %in% colnames(counts)[bender.keep]
DimPlot(retained, group.by="emptynn",reduction = 'tsne',
        cols=c('grey','steelblue3'))+ggtitle("EmptyNN")+NoLegend()
