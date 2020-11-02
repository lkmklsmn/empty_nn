\name{runSeurat}
\alias{runSeurat}
\title{runSeurat}
\description{
run Seurat analysis pipeline}
\usage{
runSeurat(counts,gene.use,resolution)
}
\arguments{
  \item{counts}{
count matrix
}
  \item{genes.use}{
highly variable genes to use for downstream analysis, eg PCA and tSNE
}
  \item{resolution}{
a number from 0 to 1.
high resolution will result in more clusters.
}
}
\details{
}
\value{
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
res <- runSeurat(counts,resolution=0.2)
}

