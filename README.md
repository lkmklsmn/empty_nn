# EmptyNN
EmptyNN is a novel **cell-calling algorithm** based on **Positive-unlabeled (PU)** learning which **removes cell-free droplets and recovers lost cells** in droplet-based single cell RNA sequencing data. For more info please see our [EmptyNN preprint](https://www.biorxiv.org/content/10.1101/2021.01.15.426387v1).

<p align="center">
<img src="Figure 1.png">
</p>

**Workflow**. EmptyNN leverages positive unlabeled learning to classify cell-free and cell-containing droplets. (a) Cells and barcodes are combined in oil droplets. Some droplets may lack a cell but contain ambient RNA. The EmptyNN classifier distinguishes cell-free from cell-containing droplets. (b) Schematic describes the workflow of EmptyNN. The black curve represents the distribution of total counts (y-axis) across sorted barcodes (x-axis). The blue bar represents sets of barcodes with very low total counts, set *P*. The grey bar represents barcodes with higher total counts consisting of cell-containing and cell-free droplets, set *U*. EmptyNN trains a classifier, where barcodes from *P* are labeled as cell-free droplets (blue) and a fraction of barcodes from *U* is labeled as cell-containing droplets (pink). The classifier is applied to the remaining barcodes in *U* and the predictions are recorded. During each *k* fold, each barcode in *U* is predicted *k*-1 times. The above process is repeated for *N* iterations (default: 10). The average prediction probability of each barcode in *U* defines each barcode as a cell-free or cell-containing droplet. 

## Reproducibility
To reproduce the analysis and figures presented in our manuscript please see the [*Reproducibility*](https://github.com/lkmklsmn/empty_nn/tree/master/Reproducibility) folder.

## Tutorial
Check out our jupyter notebook (in R environment) tutorial at [*EmptyNN - Cell Hashing Dataset Tutorial*](https://github.com/lkmklsmn/empty_nn/blob/master/tutorial/EmptyNN%20-%20Cell%20Hashing%20Dataset%20Tutorial.ipynb).

## Installation
EmptyNN is implemented in R and depends on the **keras** and **Matrix** R packages.

#### Option 1
```
$ git clone http://github.com/lkmklsmn/empty_nn
$ cd empty_nn

## enter R and install packages
$ R

> install.packages("EmptyNN_1.0.tar.gz", repos = NULL, type = "source")
```
#### Option 2
```
> install.packages("devtools")
> library(devtools)
> install_github("lkmklsmn/empty_nn")
```

## Usage

#### Input
1. Raw unfiltered count matrix in mtx or h5 format

#### Output
1. A boolean vector showing the predictions for cell-free or cell-containing droplets
2. Probabilities for each barcode in set *U*

#### Download example datasets
```
$ cd empty_nn
$ sh ./code/download_data.sh
```

```
library(EmptyNN)

# Load data
counts <- Read10X_h5("./data/example_data.h5", use.names = TRUE, unique.features = TRUE)

# Run emptynn()
nn.res <- emptynn(counts, threshold = 100, k_folds = 10, iteration = 10, verbose = TRUE)

# Downstream analysis
retained <- runSeurat(counts = counts[, nn.res$nn.keep], resolution = 0.2)
DimPlot(retained,reduction = 'tsne') + ggtitle("EmptyNN") + NoLegend()
```
