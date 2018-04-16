\name{BatchAdjust}
\alias{BatchAdjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Helper function for applying the ComBat batch-adjusting methodology to an ExpressionSet object
}
\description{
  This is a convenience function for adjusting an ExpressionSet object for batch effects.
}
\usage{
BatchAdjust(expressionSet, batchVariableName, covariateVariableNames=c())
}
\arguments{
  \item{expressionSet}{
    An \code{ExpressionSet} object that contains data to be adjusted for batch effects. This object should contain at least one \code{phenoData} column that indicates the batch number for each sample. Covariates may also be specified in other \code{phenoData} columns. Required.
}
  \item{batchVariableName}{
    Name of \code{phenoData} column that indicates the batch number for each sample. Required.
}
  \item{covariateVariableNames}{
    Name(s) of \code{phenoData} column(s) that indicate the covariate information for each sample. Optional.
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
  This function provides a convenient way to batch adjust an \code{ExpressionSet} object. It invokes the \code{ComBat} function within the \code{sva} package. Please see that package for additional details about how batch adjusting is performed. If multiple covariates are specified, a model matrix will be constructed that accounts for all covariates.
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
