# Load R libs ####
library(Seurat)
library(Matrix)
library(keras)
library(dplyr)
library(DropletUtils)

# Download dataset from website ####
website <- "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_raw_gene_bc_matrices_h5.h5"
download.file(website,destfile="10x/neurons_900_raw_gene_bc_matrices_h5.h5")
counts <- Read10X_h5("10x/neurons_900_raw_gene_bc_matrices_h5.h5", 
                     use.names = TRUE, unique.features = TRUE)

# Empty_NN ####
source("train.R")
# Input user-defined threshold
res <- Empty_NN(t(counts),batch_size = 32,threshold=3000)
# Input expected no. of cells
res <- Empty_NN(t(counts),batch_size = 32,expected = 900)

# Find highly variable genes ####
tmp <- CreateSeuratObject(counts = counts[,res$nn.keep])
tmp <- NormalizeData(tmp)
tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 1000)
gene.use <- VariableFeatures(tmp)

# Recovered barcodes with total counts below threshold(200) ####
preds <- intersect(names(res$preds[res$preds>0.5]),
                  names(res$n_counts[res$n_counts<threshold]))
recover <- runSeurat(counts=counts[,preds],gene.use=gene.use,resolution=0.2)
DimPlot(recover,label=TRUE)

# interneuron marker genes ####
genes <- c("Gad1", "Gad2", "Slc6a1")
FeaturePlot(recover,features=genes)

