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
# Calculate total reads ####
n_counts <- Matrix::rowSums(counts)
names(n_counts) <- rownames(counts)

# Load classification.global ####
global <- read.delim(paste0(path,"classification.global.txt",sep="\t",header=T))
global$label <- ifelse(global$classification.global=="Negative","negative","singlet")
global <- global[rownames(global) %in% rownames(counts),]

# Empty_NN ####
source(paste0(path,"train.R"))
res <- Empty_NN(counts,batch_size = 64)

# Empty_NN Results ####
preds <- res[[2]]
bcs <- intersect(names(n_counts[n_counts >100 & n_counts <400]),names(preds))
nn.res <- global[bcs,2,drop=F]
nn.res$preds <- preds[bcs]
t <- table(nn.res$preds>0.5,nn.res$label)
sum(diag(t))/sum(t)
p1 <- ggplot(nn.res, aes(preds>0.5, ..count..)) + geom_bar(aes(fill = label), position = "dodge") +
        ylab("no. of cells")+ylim(0,3600)+labs(title = "Empty_NN prediction",
                                  subtitle = "total_counts:from 100 to 400, accuracy: 0.83")

# Compare with EmptyDrops ####
library(DropletUtils)
tmp <- t(counts)
br.out <- barcodeRanks(tmp)
knee <- br.out$knee #  15516 
inflection <- br.out$inflection #  6667
e.out <- emptyDrops(tmp)
e.keep <- e.out$FDR <= 0.001
summary(e.keep)
#Mode   FALSE    TRUE    NA's 
#logical   11288   11619   16935 
e.res <- data.frame(br.out,e.keep)
rownames(e.res) <- rownames(counts)
e.res$label <- global[match(rownames(e.res),rownames(global)),2]
sub.e.res <- subset(e.res, total > 100 & total < 400)
t <- table(sub.e.res$e.keep,sub.e.res$label)
sum(diag(t))/sum(t)
p2 <- ggplot(sub.e.res, aes(e.keep, ..count..)) + geom_bar(aes(fill = label), position = "dodge") +
        ylab("no. of cells")+ylim(0,3600)+labs(title = "EmptyDrops prediction",
                                  subtitle = "total_counts:from 100 to 400, accuracy: 0.55")
p1|p2

# CellRanger 2.0 ####
c.keep <- defaultDrops(tmp, expected=16000)
summary(c.keep)
# cutoff: 3903










