\name{UPC_TwoColor}
\alias{UPC_TwoColor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Universal exPression Codes (UPC) for two-channel microarrays
}
\description{
  This function is used to normalize two-channel expression microarrays (from Agilent) using the Universal exPression Codes (UPC) approach. In raw form, such microarray data come in the form of tab-separate data files.
}
\usage{
UPC_TwoColor(inFilePattern, outFilePath = NA, modelType="nn",
  convThreshold=0.01, batchFilePath = NA, verbose = TRUE)
}
\arguments{
  \item{inFilePattern}{
     Absolute or relative path to the input file to be processed. To process multiple files, wildcard characters can be  used (e.g., "*.txt"). Alternatively, a Gene Expression Omnibus identifier (e.g., GSE39655 or GSM1072833) can be specified.(This is the only required parameter.)
}
  \item{outFilePath}{
     Absolute or relative path where the output file will be saved. (This parameter is optional.)
}
  \item{modelType}{
     Various models can be used for the mixture model to differentiate between active and inactive probes. The default is the normal-normal model (``nn''), which uses the normal distribution. Other available options are log-normal (``ln''), negative-binomial (``nb''), and normal-normal Bayes (``nn_bayes'').
}
  \item{convThreshold}{
     Convergence threshold that determines at what point the mixture-model parameters have stabilized. The default value should be suitable in most cases. However, if the model fails to converge, it may be useful to adjust this value. Optional.
}
  \item{batchFilePath}{
     Absolute or relative path to a tab-separated text file that indicates batch (and optionally, covariate information) for each sample. Optional.
}
  \item{verbose}{
     Whether to output more detailed status information as files are processed. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list is returned, containing two elements: a matrix containing UPC values and a vector of probe names that correspond to each row of the matrix. The matrix will contain two columns---one corresponding to each channel---for each sample. When the array design uses duplicate probe names (this is common for control probes), the vector of probe names will also contain duplicates.
}
\references{
  Piccolo SR, Withers MR, Francis OE, Bild AH and Johnson WE. Multi-platform single-sample estimates of transcriptional activation. Proceedings of the National Academy of Sciences of the United States of America, 2013, 110:44 17778-17783.
}
\author{
  Stephen R. Piccolo
}
\note{
  If a Gene Expression Omnibus (GEO) identifier is specified for the \code{inFilePattern} parameter, an attempt will be made to download the sample(s) directly from GEO. If a study identifier (e.g., GSE39655) is specified, all CEL files from that study will be downloaded. If a sample identifier (e.g., GSM1072833) is specified, only that sample will be downloaded.

  The \code{batchFilePath} parameter provides a convenient way to adjust the data for batch effects. It invokes the \code{ComBat} function within the \code{sva} package. Please see that package for additional details about how batch adjusting is performed. Batch adjusting is performed before UPC transformation occurs.

  The \code{modelType} parameter indicates which type of mixture model to use for UPC transformation. The "nn_bayes" model type is an experimental new approach intended for experiments where a subset of genes are expressed at extreme levels.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
# Normalize a file from GEO and save output to a file
result = UPC_TwoColor("GSM1072833", "output_file.txt")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
