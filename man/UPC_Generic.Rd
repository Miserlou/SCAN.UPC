\name{UPC_Generic}
\alias{UPC_Generic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generic function to apply Universal exPression Codes (UPC) transformation
}
\description{
This function can be used to derive UPC values to any type of gene-expression data. It requires the user to specify expression values for many genes (or transcripts, exons, or probes). And optionally, the user can specify the length and/or GC content (proportion of G or C bases) for the corresponding genomic region (e.g., gene). If these values are specified, the UPC algorithm corrects for biases resulting from length or GC content.
}
\usage{
UPC_Generic(expressionValues, lengths = NULL, gcContent = NULL, modelType = "nn",
  convThreshold = 0.001, higherValuesIndicateHigherExpression = TRUE, verbose = TRUE)
}
\arguments{
  \item{expressionValues}{
     A numeric vector that indicates expression levels. This is the only required parameter.
}
  \item{lengths}{
     An integer vector that indicates the length of the genomic region that corresponds to each expression value. The order of these values must be identical to the order of the expression values. This parameter is optional.
}
  \item{gcContent}{
     A numeric vector that indicates the proportion of bases in the genomic region that are a G or C nucleotide. The order of these values must be identical to the order of the expression values. This parameter is optional.
}
  \item{modelType}{
     Various models can be used for the mixture model to differentiate between active and inactive probes. The default is the normal-normal model (``nn''), which uses the normal distribution. Other available options are log-normal (``ln''), negative-binomial (``nb''), and normal-normal Bayes (``nn_bayes'').
}
  \item{convThreshold}{
     Convergence threshold that determines at what point the mixture-model parameters have stabilized. The default value should be suitable in most cases. However, if the model fails to converge (or converges too quickly), it may be useful to adjust this value. Optional.
}
  \item{higherValuesIndicateHigherExpression}{
     In most cases, higher expression values indicate relatively high expression. However, if higher values indicate relatively low expresssion, this parameter can be used to indicate such. Accordingly, UPC values closer to one will indicate higher expression, and UPC values closer to zero will indicate the opposite. Optional.
}
  \item{verbose}{
     Whether to output more detailed status information as processing occurs. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An numeric vector with UPC values whose order corresponds with the order of the input values.
}
\references{
  Piccolo SR, Withers MR, Francis OE, Bild AH and Johnson WE. Multi-platform single-sample estimates of transcriptional activation. Proceedings of the National Academy of Sciences of the United States of America, 2013, 110:44 17778-17783.
}
\author{
  Stephen R. Piccolo
}
\note{
  The \code{modelType} parameter indicates which type of mixture model to use for UPC transformation. The "nn_bayes" model type is an experimental new approach intended for experiments where a subset of genes are expressed at extreme levels.
}
% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%\dontrun{
%upcValues = UPC_Generic(exprValues)
%}
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
