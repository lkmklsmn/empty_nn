# EmptyNN
**EmptyNN** is a novel **cell-calling algorithm** based on **Positive-unlabeled (PU) learning** to **remove cell-free droplets and recover lost cells** in droplet-based single cell RNA sequencing data.

<p align="center">
<img src="Figure 1.png">
</p>

**Workflow. EmptyNN leverages positive unlabeled learning to classify cell-free and cell-containing droplets.** (a) Cells and barcodes are combined in oil droplets. Some droplets may lack a cell but contain ambient RNA. The EmptyNN classifier distinguishes cell-free droplets from cell-containing droplets. (b) Schematic describes the workflow of EmptyNN. Blue curve represents the distribution of total counts (y-axis) across sorted barcodes (x-axis). The grey bar represents sets of barcodes with very low total counts, set *P*. The yellow bar represents barcodes with higher total counts consisting of cell-containing (black dots) and cell-free droplets (open dots), namely set *U*. EmptyNN trains a classifier, where barcodes from *P* are labeled as cell-free droplets. A random sample of barcodes from *U* is labeled as cell-containing droplets. The classifier is applied to the remaining barcodes and the scores are recorded. During each k training fold, each barcode in *U* receives (k-1) scores. The above process repeats *N* iterations. Based on the average prediction, each barcode in *U* will be classified as a cell-free or cell-containing droplet.

## Reproducibility
To reproduce the analysis and figures presented in our manuscript please see the [*Reproducibility*](https://github.com/lkmklsmn/emptynn/blob/master/Reproducibility) folder.

## Tutorial
Check out our jupyter notebook (in R environment) tutorial at [*EmptyNN - Cell Hashing Dataset Tutorial*](https://github.com/lkmklsmn/empty_nn/blob/master/Reproducibility/EmptyNN%20-%20Cell%20Hashing%20Dataset%20Tutorial.ipynb).

## Installation
EmptyNN is an R package. The required packages include **keras** and **Matrix**.

### option 1
```
$ git clone http://github.com/lkmklsmn/emptynn
$ cd emptynn
## enter R and install packages
$ R
> install.packages("EmptyNN_1.0.tar.gz", repos = NULL, type = "source")
```
### option 2
```
> install.packages("devtools")
> library(devtools)
> install_github("lkmklsmn/emptynn")
```
## Input
1. Raw count expression matrix in mtx or h5 format

## Output
1. A boolean vector showing cell-free or cell-containing droplets
2. Predicted probabilities for each barcode set *U*

## Download example datasets
```
$ cd emptynn
$ Rscript "./code/download.data.R"
```
## Usage
```
library(EmptyNN)
# Load data
counts <- Read10X_h5("./data/example_data.h5", use.names = TRUE, unique.features = TRUE)
# Run emptynn()
nn.res <- emptynn(counts,threshold=100,k_folds=10,iteration=10,verbose=TRUE)
# Downstream analysis
retained <- runSeurat(counts=counts[,nn.res$nn.keep],resolution=0.2)
DimPlot(retained,reduction = 'tsne')+ggtitle("EmptyNN")+NoLegend()
```
