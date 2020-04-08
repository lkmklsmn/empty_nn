path = "/Users/yanfangfang/Downloads/Empty_Droplet/data/50,000/"
# Load R libs ####
library(Seurat)
library(Matrix)
library(keras)
library(ggplot2)
library(patchwork)

# Load count table ####
counts <- readMM(file = paste0(path,"Stoeckius_RNA.mtx"))
barcodes <- read.csv(paste0(path,"Stoeckius_barcodes.txt"))[,2]
genes <- read.csv(paste0(path,"Stoeckius_genes.txt"))[,2]
rownames(counts) <- barcodes
colnames(counts) <- genes
# keep barcodes that have HTO counts
jointbcs <- read.csv(paste0(path,"joint_barcodes.txt",header=F))[,1]
counts <- counts[rownames(counts) %in% jointbcs,]

# Load classification.global ####
global <- read.delim(paste0(path,"classification.global.txt",sep="\t",header=T))
global$label <- ifelse(global$classification.global=="Negative","negative","singlet")
global <- global[rownames(global) %in% rownames(counts),]

# Empty_NN ####
source(paste0(path,"train.R"))
res <- Empty_NN(counts,batch_size = 64,threshold=200,k=4)

# Boxplot showing prediction accuracies ####
bcs <- intersect(names(res$n_counts[res$n_counts >100 & res$n_counts <400]),names(res$preds))
nn.res <- global[bcs,2,drop=F]
nn.res$preds <- res$preds[bcs]
t <- table(nn.res$preds>0.5,nn.res$label)
acc <- sum(diag(t))/sum(t)
ggplot(nn.res, aes(preds>0.5, ..count..)) + geom_bar(aes(fill = label), position = "dodge") +
        ylab("no. of cells")+ylim(0,3600)+
        labs(title = "Empty_NN prediction",
        subtitle = paste0("total_counts:from 100 to 400, accuracy: ",round(acc,3)))












