path = "/Users/yanfangfang/Downloads/Empty_Droplet/"
# Load R libs ####
library(Seurat)
library(Matrix)
library(keras)
library(DropletUtils)
library(dplyr)

# Download data from website ####
website <- "http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_raw_feature_bc_matrix.h5"
download.file(website,destfile="10x/pbmc_10k_raw.h5")
counts <- Read10X_h5("10x/pbmc_10k_raw.h5",use.names = TRUE, unique.features = TRUE)

# Negative learning ####
source(paste0(path,"code/0629/Negative_learning.R"))
neg.res <- neg_learning(t(counts),neg.threshold=100,n_splits=5,iteration=5,verbose=TRUE)

# Curriculum Learning Empty_NN ####
source(paste0(path,"code/0629/Curri_NN.R"))
## input can be "expected no. of cells (eg.expected=10000)" or "user defined threshold (eg.threshold=3000)
curri.res <- Empty_NN(t(counts),batch_size = 64,expected = 10000)

# EmptyDrops ####
e.out <- emptyDrops(counts)
e.keep <- e.out$FDR <= 0.001

# save(neg.res,curri.res,e.keep,file="pbmc_result.RData")

# Benchmark analysis ####
# Input is a vector showing cells or not
benchmark <- function(counts,keep.vector,model){
        counts.keep <- counts[,keep.vector & !is.na(keep.vector)]
        # recovered barcodes with low total counts
        n_counts <- Matrix::colSums(counts)
        names(n_counts) <- colnames(counts)
        recover_bcs <- intersect(colnames(counts.keep),names(n_counts[n_counts<threshold]))
        print(paste0("recovered ",length(recover_bcs)," barcodes"))
        if (length(recover_bcs)==0){
                stop("no recovered barcodes")
        }
        # Run Seurat
        # find gene.use
        tmp <- CreateSeuratObject(counts = counts.keep)
        tmp <- NormalizeData(tmp)
        tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = 1000)
        gene.use <- VariableFeatures(tmp)
        recover <- runSeurat(counts=counts[-grep("^MT", rownames(counts)),recover_bcs],
                             gene.use=gene.use,resolution=0.2)
        DimPlot(object = recover, reduction = "umap",label=TRUE)+
                labs(title = paste0("Recovered by",model))
        return(recover)
}

# Recovered barcodes by each method ####
recover <- benchmark(counts,neg.res$nn.keep,model="negatove learning")
recover <- benchmark(counts,curri.res$nn.keep,model="curriculum learning")
recover <- benchmark(counts,e.keep,model="EmptyDrops")





