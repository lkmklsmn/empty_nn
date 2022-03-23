\name{runSeurat}
\alias{runSeurat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
standard seurat preprocessing step
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runSeurat(counts, gene.use = FALSE, resolution)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{gene.use}{
%%     ~~Describe \code{gene.use} here~~
}
  \item{resolution}{
%%     ~~Describe \code{resolution} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (counts, gene.use = FALSE, resolution) 
{
    require(Seurat)
    tmp <- CreateSeuratObject(counts = counts)
    tmp <- tmp[-grep("^RPS", rownames(tmp)), ]
    tmp <- tmp[-grep("^RPL", rownames(tmp)), ]
    tmp <- NormalizeData(tmp)
    tmp <- FindVariableFeatures(tmp)
    if (!gene.use) {
        gene.use <- VariableFeatures(tmp)
    }
    else {
        gene.use <- gene.use[gene.use \%in\% rownames(tmp)]
    }
    tmp <- ScaleData(tmp, features = gene.use)
    tmp <- RunPCA(tmp, features = gene.use)
    tmp <- FindNeighbors(tmp, dims = 1:10)
    tmp <- FindClusters(tmp, resolution = resolution)
    tmp <- RunTSNE(tmp, dims = 1:10, check_duplicates = FALSE)
    return(tmp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
