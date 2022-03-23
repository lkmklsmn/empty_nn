\name{neg_create_model}
\alias{neg_create_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create deep learning using positive unlabeled learning
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
neg_create_model(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
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
function (input) 
{
    model <- keras_model_sequential()
    model \%>\% layer_dense(units = 256, activation = "relu", input_shape = input) \%>\% 
        layer_dropout(rate = 0.4) \%>\% layer_dense(units = 128, 
        activation = "relu") \%>\% layer_dropout(rate = 0.3) \%>\% 
        layer_dense(units = 2, activation = "softmax")
    model \%>\% compile(loss = "binary_crossentropy", optimizer = optimizer_rmsprop(), 
        metrics = c("accuracy"))
    model
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
