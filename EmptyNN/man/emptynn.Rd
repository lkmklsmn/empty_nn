\name{emptynn}
\alias{emptynn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EmptyNN is a novel cell-calling algorithm which removes cell-free droplets and recovers lost cells in droplet-based single cell RNA sequencing data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
emptynn(counts, threshold = 100, k = 5, iteration = 5, batch_size = 16, epoch = 10, verbose = TRUE, training_verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{iteration}{
%%     ~~Describe \code{iteration} here~~
}
  \item{batch_size}{
%%     ~~Describe \code{batch_size} here~~
}
  \item{epoch}{
%%     ~~Describe \code{epoch} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{training_verbose}{
%%     ~~Describe \code{training_verbose} here~~
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
function (counts, threshold = 100, k = 5, iteration = 5, batch_size = 16, 
    epoch = 10, verbose = TRUE, training_verbose = 0) 
{
    require(keras)
    require(Matrix)
    if (nrow(counts) < ncol(counts)) {
        stop(paste0("Please transpose counts matrix before running EmptyNN\n", 
            "  rows are barcodes, columns are genes"))
    }
    n_counts <- Matrix::rowSums(counts)
    names(n_counts) <- rownames(counts)
    n_counts <- n_counts[n_counts >= 10]
    negative <- names(n_counts[which(n_counts <= threshold)])
    if (verbose) {
        print(paste0("there are ", length(negative), " in P set"))
    }
    if (length(negative) >= 10000) {
        negative <- sample(negative, 10000)
    }
    else {
        negative <- negative
    }
    unlabel <- names(n_counts[which(n_counts > threshold)])
    if (verbose) {
        print(paste0("there are ", length(unlabel), " in U set"))
        print(paste0("Samples in U set were split into ", k, 
            " folds"))
    }
    gene.use <- names(tail(sort(Matrix::colSums(counts[negative, 
        ])), 2000))
    counts_2k <- counts[names(n_counts), gene.use]
    if (verbose) {
        print("data normalization")
    }
    norm.counts.2k <- sweep(counts_2k, 1, n_counts, "/")
    cv_p <- data.frame(matrix(0, ncol = iteration, nrow = length(unlabel)))
    rownames(cv_p) <- unlabel
    colnames(cv_p) <- paste0("iteration", seq(1, ncol(cv_p)))
    if (verbose) {
        print("start training")
    }
    for (i in seq(iteration)) {
        if (verbose) {
            print(paste0("iteration ", i))
        }
        df <- data.frame(total_counts = n_counts[unlabel], label = 1)
        require(caret)
        df$folds <- createFolds(factor(df$label), k = k, list = FALSE)
        for (j in seq(1, k)) {
            if (verbose) {
                print(paste0("training fold ", j))
            }
            train <- c(negative, rownames(df[df$folds == j, ]))
            test <- rownames(df[df$folds != j, ])
            x_train <- data.matrix(norm.counts.2k[train, ])
            x_test <- data.matrix(norm.counts.2k[test, ])
            y_train <- c(rep(0, length(negative)), rep(1, nrow(df[df$folds == 
                j, ])))
            random_indices <- sample(1:length(y_train))
            x_train <- x_train[random_indices, ]
            y_train <- y_train[random_indices]
            y_train <- to_categorical(y_train)
            model_neg <- neg_create_model(input = dim(x_train)[2])
            model_neg \%>\% fit(x_train, y_train, batch_size = batch_size, 
                epochs = epoch, verbose = training_verbose, validation_split = 0.2)
            y_test_pred <- model_neg \%>\% predict(x_test)
            cv_p[test, i] <- y_test_pred[, 2]
        }
    }
    cv_p$mean.crossval <- apply(cv_p, 1, mean)
    nn_bcs <- rownames(cv_p[cv_p$mean.crossval > 0.5, ])
    nn.keep <- rownames(counts) \%in\% nn_bcs
    return(list(nn.keep = nn.keep, prediction = cv_p))
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
