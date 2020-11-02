# Load R libs ####
library(EmptyNN)
library(Seurat)
library(keras)
library(DropletUtils)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)

# Load data ####
# The file can be downloaded from
# https://drive.google.com/file/d/12y0fW_Y9OdhBLns_2gpjo2Xq25c4qnGY/view?usp=sharing
# It contains two objects: Stoeckius.counts, label
load("./../data/Stoeckius/Stoeckius.RData")

# EmptyNN ####
nn.res <- emptynn(Stoeckius.counts,threshold=50,k=30,iteration=5)
nn.keep <- nn.res$nn.keep

# CellRanger 2.0 ####
n_counts <- rowSums(Stoeckius.counts)
ranger.keep <- n_counts>=200

# EmptyDrops ####
e.out <- emptyDrops(t(Stoeckius.counts))
e.keep <- e.out$FDR <= 0.001

# CellBender ####
bender <- Read10X_h5("./../data/Stoeckius.CellBender_filtered.h5")
bender.keep <- rownames(Stoeckius.counts) %in% colnames(bender)
names(bender.keep) <- rownames(Stoeckius.counts)

# Comparison ####
# Overall ROC curve ####
library('pROC')
bcs <- rownames(nn.res$prediction)
rocobj1 <- roc(label[bcs,2], neg.res.5ite.30splits$prediction[bcs,'mean.crossval'])
rocobj2 <- roc(label[names(ranger.keep),2], as.numeric(ranger.keep))
bcs <- rownames(e.out[!is.na(e.out$FDR),])
rocobj3 <- roc(label[bcs,2], e.out[bcs,'FDR'])
rocobj4 <- roc(label[names(bender.keep),2], as.numeric(bender.keep))
ggroc(list("EmptyNN, 0.9473"=rocobj1,
           "CellRanger 2.0, 0.8688"=rocobj2,
           "EmptyDrops, 0.7967"=rocobj3,
           "CellBender, 0.888"=rocobj4),legacy.axes = TRUE) +
        labs(x = "FPR", y = "TPR")+
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype="dashed")

# tSNE plot ####
keep <- nn.keep | ranger.keep | (e.keep & !is.na(e.keep)) | bender.keep
retained <- runSeurat(counts=counts[,keep],genes.use=genes.use,resolution=0.2)
retained$emptynn <- colnames(retained) %in% colnames(counts)[nn.keep]
retained$cellranger <- colnames(retained) %in% colnames(counts)[ranger.keep]
retained$emptydrops <- colnames(retained) %in% colnames(counts)[e.keep & !is.na(e.keep)]
retained$cellbender <- colnames(retained) %in% colnames(counts)[bender.keep]
DimPlot(retained, group.by="emptynn",reduction = 'tsne',
        cols=c('grey','steelblue3'))+ggtitle("EmptyNN")+NoLegend()
