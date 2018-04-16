\name{UPC_RNASeq_Single}
\alias{UPC_RNASeq_Single}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Apply Universal exPression Codes (UPC) transformation to a single RNA-Seq sample
}
\description{
This function is used to derive UPC values for a single RNA-Seq sample. It requires an input vector that specifies a read count for each genomic region (e.g., gene). Optionally, this function can correct for the GC content and length of each genomic region.
}
\usage{
UPC_RNASeq_Single(expressionValues, featureNames, lengths = NULL,
  gcContent = NULL, modelType = "nn", convThreshold = 0.01,
  ignoreZeroes = FALSE, verbose = TRUE)
}
\arguments{
  \item{expressionValues}{
     A vector of RNA-Seq count values. Required.
}
  \item{featureNames}{
     A vector of unique names that correspond to the count values. Required.
}
  \item{lengths}{
     A vector indicating the length (in genomic bases) of the genomic region that corresponds to the count values.
}
  \item{gcContent}{
     A vector indicating the number of G/C bases in the genomic region that corresponds to the count values.
}
  \item{modelType}{
     Various models can be used for the mixture model to differentiate between active and inactive probes. The default is the normal-normal model (``nn''), which uses the normal distribution. Other available options are log-normal (``ln''), negative-binomial (``nb''), and normal-normal Bayes (``nn_bayes'').
}
  \item{convThreshold}{
     Convergence threshold that determines at what point the mixture-model parameters have stabilized. The default value should be suitable in most cases. However, if the model fails to converge (or converges too quickly), it may be useful to adjust this value. (This parameter is optional.)
}
  \item{ignoreZeroes}{
     Whether to ignore read counts equal to zero when performing UPC calculations. Default is FALSE.
}
  \item{verbose}{
     Whether to output more detailed status information as files are normalized. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A vector that contains a UPC value for each probeset/gene/transcript.
}
\references{
  Piccolo SR, Withers MR, Francis OE, Bild AH and Johnson WE. Multi-platform single-sample estimates of transcriptional activation. Proceedings of the National Academy of Sciences of the United States of America, 2013, 110:44 17778-17783.
}
\author{
  Stephen R. Piccolo
}
\note{
  RNA-Seq data by nature have a lot of zero read counts. Samples with an excessive number of zeroes may lead to error messages because genes cannot be allocated properly to bins. The user can specify \code{ignoreZeroes=TRUE} to avoid this error. In practice, we have seen that the resulting UPC values are similar with either approach.

  The \code{modelType} parameter indicates which type of mixture model to use for UPC transformation. The "nn_bayes" model type is an experimental new approach intended for experiments where a subset of genes are expressed at extreme levels.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%\dontrun{
%result = UPC_RNASeq("ReadCounts.txt", "Annotation.txt")
%}
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
