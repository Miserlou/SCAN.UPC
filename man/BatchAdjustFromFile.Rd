\name{BatchAdjustFromFile}
\alias{BatchAdjustFromFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Helper function for applying the ComBat batch-adjusting methodology to an ExpressionSet object
}
\description{
  This is a convenience function for adjusting an ExpressionSet object for batch effects.
}
\usage{
BatchAdjustFromFile(expressionSet, batchFilePath)
}
\arguments{
  \item{expressionSet}{
    An \code{ExpressionSet} object that contains data to be adjusted for batch effects. This object should contain at least one \code{phenoData} column that indicates the batch number for each sample. Covariates may also be specified in other \code{phenoData} columns. Required.
}
  \item{batchFilePath}{
    Tab-separated text file that indicates batch (and optionally, covariate information) for each sample. Required.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  An ExpressionSet object indicating expression levels after batch adjustment.
}
%\references{
%  TBD
%}
\author{
  Stephen R. Piccolo
}
\note{
  This function provides a convenient way to batch adjust an \code{ExpressionSet} object. It invokes the \code{ComBat} function within the \code{sva} package. Please see that package for additional details about how batch adjusting is performed.

  The batch file must contain a row for each of the samples in the \code{ExpressionSet} object, and one of the columns must use the name ``batch'' to indicate that the values represent batch number. The remaining column(s), if present, will be treated as covariates.  If multiple covariates are specified, a model matrix will be constructed that accounts for all covariates.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%\dontrun{
%pkgName = InstallBrainArrayPackage(celFilePath, "17.0.1", "hs", "entrezg")
%}
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
