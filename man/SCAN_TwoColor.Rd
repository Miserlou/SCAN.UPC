\name{SCAN_TwoColor}
\alias{SCAN_TwoColor}
\title{
  Single-Channel Array Normalization (SCAN) for Agilent two-color expression microarrays
}
\description{
  This function is used to normalize Agilent two-color expression microarrays via the SCAN method.
}
\usage{
SCAN_TwoColor(inFilePattern, outFilePath = NA, batchFilePath = NA, verbose = TRUE)
}
\arguments{
  \item{inFilePattern}{
     Absolute or relative path to the input file to be processed. To process multiple files, wildcard characters can be  used (e.g., "*.txt"). Alternatively, a Gene Expression Omnibus identifier (e.g., GSE39655 or GSM1072833) can be specified. This is the only required parameter.
}
  \item{outFilePath}{
     Absolute or relative path where the output file will be saved. This is optional.
}
  \item{batchFilePath}{
     Absolute or relative path to a tab-separated text file that indicates batch (and optionally, covariate information) for each sample. Optional.
}
  \item{verbose}{
     Whether to output more detailed status information as files are normalized. Default is TRUE.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list is returned, containing two elements: a matrix containing normalized data values and a vector of probe names that correspond to each row of the matrix. The matrix will contain two columns---one corresponding to each channel---for each sample. When the array design contains duplicate probe names (this is common for control probes), the vector of probe names will also contain duplicates.
}
\references{
  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
}
\author{
  Stephen R. Piccolo
}
\note{
  If a Gene Expression Omnibus (GEO) identifier is specified for the \code{inFilePattern} parameter, an attempt will be made to download the sample(s) directly from GEO. If a study identifier (e.g., GSE39655) is specified, all CEL files from that study will be downloaded. If a sample identifier (e.g., GSM1072833) is specified, only that sample will be downloaded.

  The \code{batchFilePath} parameter provides a convenient way to adjust the data for batch effects. It invokes the \code{ComBat} function within the \code{sva} package. Please see that package for additional details about how batch adjusting is performed. Batch adjusting is performed after values have been SCAN normalized.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
# Normalize a file from GEO and save output to a file
result = SCAN_TwoColor("GSM1072833", "output_file.txt")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
