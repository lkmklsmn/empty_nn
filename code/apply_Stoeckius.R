path = "/Users/yanfangfang/Downloads/Empty_Droplet/"
# Load R libs ####
library(Seurat)
library(Matrix)
library(keras)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(DropletUtils)

# Load count and label ####
# The file can be downloaded from 
# https://drive.google.com/file/d/12y0fW_Y9OdhBLns_2gpjo2Xq25c4qnGY/view?usp=sharing
# It contains two objects: Stoeckius.counts, label
load(paste0(path, "data/Stoeckius/Stoeckius.RData"))

# Negative Learning ####
source(paste0(path,"code/0629/Negative_learning.R"))
neg.res <- neg_learning(Stoeckius.counts,neg.threshold=50,n_splits=5,iteration=5,
                        verbose=TRUE,training_verbose = 0)

# Curriculum Learning EmptyNN ####
source(paste0(path,"code/0629/Curri_NN.R"))
curri.res <- Empty_NN(Stoeckius.counts,batch_size = 64,threshold=200,k=4)

# EmptyDrops ####
e.out <- emptyDrops(t(Stoeckius.counts))
e.keep <- e.out$FDR <= 0.001

# save(neg.res,curri.res,e.keep,file="stoeckius_result.RData")
load(file="stoeckius_result.RData")

# Benchmark Analysis ####
# accuracy in each quantile
acc_quantile <- function(counts,keep.vector,model){
        n_counts <- Matrix::rowSums(counts)
        res <- data.frame("predictions"=keep.vector,"counts"=n_counts)
        idx <- match(rownames(res),rownames(label))
        res$label <- label[idx,2]
        res <- res[which(res$counts>10),]
        res$quantile <- ntile(-res$counts, 20)
        res$quantile <- paste0("q_",res$quantile)
        level_order <- c(paste0("q_",seq(1,20)))
        df <- data.frame("quantile"=level_order,"balanced_accuracy"=0)
        for (i in seq(1,nrow(df))){
                tmp <- res[res$quantile==df[i,1],]
                TP <- nrow(tmp[(tmp$predictions=="TRUE" & tmp$label=="singlet"),])
                TN <- nrow(tmp[(tmp$predictions=="FALSE" & tmp$label=="negative"),])
                FP <- nrow(tmp[(tmp$predictions=="TRUE" & tmp$label=="negative"),])
                FN <- nrow(tmp[(tmp$predictions=="FALSE" & tmp$label=="singlet"),])
                if (TP+FN==0){
                        sensitivity = 0
                } else {sensitivity <- TP / (TP + FN)}
                if (TN + FP==0){
                        specificity = 0
                } else {specificity <- TN / (TN + FP)}
                balanced_accuracy <- (sensitivity + specificity) / 2
                df[i,2] <- balanced_accuracy
        }
        p <- ggplot(data=df, aes(x=factor(quantile,level = level_order), y=balanced_accuracy, fill=quantile)) +
                geom_bar(stat="identity",show.legend = FALSE)+ylim(0,1)+
                xlab("quantile")+ggtitle(paste0("Stoeckius et al ", model))
        return(p)
}
neg.plot <- acc_quantile(Stoeckius.counts,neg.res$nn.keep,model="negative learning")
curri.plot <- acc_quantile(Stoeckius.counts,curri.res$nn.keep,model="curriculum learning")
e.keep[which(is.na(e.keep))] <- "FALSE"
e.plot <- acc_quantile(Stoeckius.counts,e.keep,model="EmptyDrops")
acc_in_quantile <- neg.plot|curri.plot|e.plot

