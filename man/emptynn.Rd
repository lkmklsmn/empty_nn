\name{emptynn}
\alias{emptynn}
\title{emptynn}
\description{
a PU-learning-based method to remove empty droplets in single cell RNA-seq data
}
\usage{emptynn(counts,threshold=100,k=5,iteration=5)}
\arguments{
  \item{counts}{
raw count matrix, with rows representing barocdes and columns representing genes
}
  \item{threshold}{
a number, default is 100.
Barcodes with total counts below the threshold is P set (empty droplets).
Barcodes with total counts above the threshold is U set.
}
  \item{k}{
number of k-folds, default is 5.
}
  \item{iteration}{
number of iterations, default is 5.
}
}
\details{
}
\value{
this function returns a list of 2 elements.
The 1st elelment is a boolean vector with the same length of total barcodes, showing cells or not.
The 2nd element is the predicted probabilities for each barcode in U set.
}
\references{
}
\author{
Fangfang Yan, Lukas M. Simon
}
\seealso{
github link for emptynn
https://github.com/lkmklsmn/empty_nn
}
\examples{
library(EmptyNN)
library(Matrix)
## generate a count matrix with 10 genes and 2000 barcodes
count <- matrix(rbinom(10*2000, 10, 0.2), 10, 2000)
res <- emptynn(counts,threshold=100,k=5,iteration=5)
filtered <- count[res$nn.keep]
## filtered is the matrix with empty droplets removed.
}
\keyword{ scRNA-seq; empty droplets }
