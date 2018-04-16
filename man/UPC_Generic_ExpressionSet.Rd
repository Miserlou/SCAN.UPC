\name{UPC_Generic_ExpressionSet}
\alias{UPC_Generic_ExpressionSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generic function to apply Universal exPression Codes (UPC) transformation to ExpressionSet objects
}
\description{
  This function can be used to derive UPC values for any type of gene-expression data. The key prerequisite is that the data be placed into an ExpressionSet object. And optionally, the user can specify the length and/or GC content (proportion of G or C bases) for the corresponding genomic region (e.g., gene). If these values are specified, the UPC algorithm will correct for biases resulting from length or GC content.
}
\usage{
UPC_Generic_ExpressionSet(expressionSet, sequenceFeatureName = NA, modelType = "nn",
  convThreshold = 0.001, higherValuesIndicateHigherExpression = TRUE, verbose = TRUE)
}
\arguments{
  \item{expressionSet}{
     An \code{ExpressionSet} object that contains data to be UPC transformed. Required.
}
  \item{sequenceFeatureName}{
     Optionally, the DNA sequence corresponding to each feature (e.g., microarray probe) can be specified in the \code{ExpressionSet} object's metadata. If specified, the UPC transformation will take the G/C content and length of each feature into account. Optional.
}
  \item{modelType}{
     Various models can be used for the mixture model to differentiate between active and inactive probes. The default is the normal-normal model (``nn''), which uses the normal distribution. Other available options are log-normal (``ln''), negative-binomial (``nb''), and normal-normal Bayes (``nn_bayes'').
}
  \item{convThreshold}{
     Convergence threshold that determines at what point the mixture-model parameters have stabilized. The default value should be suitable in most cases. However, if the model fails to converge (or converges too quickly), it may be useful to adjust this value. (This parameter is optional.)
}
  \item{higherValuesIndicateHigherExpression}{
     In most cases, higher expression values indicate relatively high expression. However, if higher values indicate relatively low expresssion, this parameter can be used to indicate such. Accordingly, UPC values closer to one will indicate higher expression, and UPC values closer to zero will indicate the opposite. (This parameter is optional.)
}
  \item{verbose}{
     Whether to output more detailed status information as processing occurs. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An \code{ExpressionSet} object with UPC values.
}
\references{
  Piccolo SR, Withers MR, Francis OE, Bild AH and Johnson WE. Multi-platform single-sample estimates of transcriptional activation. Proceedings of the National Academy of Sciences of the United States of America, 2013, 110:44 17778-17783.
}
\author{
  Stephen R. Piccolo
}
\note{
  The example below illustrates how an \code{ExpressionSet} object can be retrieved from a given Gene Expression Omnibus experiment and then UPC transformed, taking each probe's sequence into account.
}
% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
upcData = UPC_Generic_ExpressionSet(getGEO("GSE56234", GSEMatrix=TRUE)[[1]], sequenceFeatureName="SEQUENCE")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
