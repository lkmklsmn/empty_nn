# Load R libs ####
library(Seurat)
library(Matrix)
library(keras)
library(DropletUtils)

# Download data from website ####
website <- "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.h5"
download.file(website,destfile="10x/pbmc_10k_raw.h5")
counts <- Read10X_h5("10x/pbmc_10k_raw.h5",use.names = TRUE, unique.features = TRUE)

# Empty_NN ####
source("train.R")
res <- Empty_NN(t(counts),batch_size = 64,threshold = 500)

# Recovered barcodes with total counts below threshold(500) ####
threshold=500
preds <- intersect(names(res$preds[res$preds>0.5]),
                   names(res$n_counts[res$n_counts<threshold]))
# Find highly variable genes
tmp <- CreateSeuratObject(counts = counts[,res$nn.keep])
tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 1000)
gene.use <- VariableFeatures(tmp)

# Run Seurat ####
recover <- runSeurat(counts=counts[-grep("^MT", rownames(counts)),preds],
                     gene.use=gene.use,resolution=0.5)
DimPlot(recover,label=TRUE)


