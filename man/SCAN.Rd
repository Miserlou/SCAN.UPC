\name{SCAN}
\alias{SCAN}
\alias{SCANfast}
\alias{UPC}
\alias{UPCfast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Single-Channel Array Normalization (SCAN) and Universal exPression Codes (UPC) for Affymetrix microarrays
}
\description{
  This function is used to normalize Affymetrix .CEL files via the SCAN and UPC methods.
}
\usage{
SCAN(celFilePattern, outFilePath = NA, convThreshold = 0.01, annotationPackageName = NA,
     probeSummaryPackage = NA, probeLevelOutDirPath = NA, exonArrayTarget=NA,
     batchFilePath=NA, verbose = TRUE)
SCANfast(celFilePattern, outFilePath = NA, convThreshold = 0.50, annotationPackageName = NA,
     probeSummaryPackage = NA, probeLevelOutDirPath = NA, exonArrayTarget=NA,
     batchFilePath=NA, verbose = TRUE)
UPC(celFilePattern, outFilePath = NA, convThreshold = 0.01, annotationPackageName = NA,
     probeSummaryPackage = NA, probeLevelOutDirPath = NA, exonArrayTarget = NA,
     modelType="nn", batchFilePath=NA, verbose = TRUE)
UPCfast(celFilePattern, outFilePath = NA, convThreshold = 0.50, annotationPackageName = NA,
     probeSummaryPackage = NA, probeLevelOutDirPath = NA, exonArrayTarget = NA,
     modelType="nn", batchFilePath=NA, verbose = TRUE)
}
\arguments{
  \item{celFilePattern}{
     Absolute or relative path to the input file to be processed. This is the only required parameter. To process multiple files, wildcard characters can be  used (e.g., "*.CEL"). Alternatively, a Gene Expression Omnibus identifier (e.g., GSE22309 or GSM555237) can be specified. 
}
  \item{outFilePath}{
     Absolute or relative path where the output file will be saved. This is optional.
}
  \item{convThreshold}{
     Convergence threshold that determines at what point the mixture-model parameters have stabilized. The default value should be suitable in most cases. However, if the model fails to converge, it may be useful to adjust this value. (This parameter is optional.)
}
  \item{annotationPackageName}{
     The name of an annotation package that specifies the layout and sequences of the probes. This is optional. By default, the correct annotation package should be identified in most cases. However, with this option allows the user to specify the package explicitly if needed.
}
  \item{probeSummaryPackage}{
     An R package that specifies alternative probe/gene mappings. This is optional. See note below for more details.
}
  \item{probeLevelOutDirPath}{
     Absolute or relative path to a directory where probe-level normalized values can be saved. This is optional. By default, the probe-level values will be discarded after they have been summarized. However, if the user has a need to repeatedly process the same file (perhaps to try various probe/gene mappings), this option can be useful because \code{SCAN} will retrieve previously normalized values if a probe-level file exists, rather than renormalize the raw data. The user should be aware that probe-level files may consume a considerable amount of disk space.
}
  \item{exonArrayTarget}{
     The type of probes to be used. This parameter is optional and should only be specified when Affymetrix Exon 1.0 ST arrays are being processed. This parameter allows the user to specify the subset of probes that should be used and how the probes should be grouped. Available options are NA, "core", "extended", "full", or "probeset". When "probeset" is used, all probes will be used, and the probes will be grouped according to the Affymetrix probeset definitions. When "core", "extended", or "full" are used, the probes that Affymetrix has defined to fall within each classification will be used, and probes will be grouped by Entrez Gene IDs (as defined in the corresponding annotation package). It is recommended to specify "probeset" when the \emph{probeSummaryPackage} parameter is being used so that all probes will be considered.
}
  \item{modelType}{
    The type of mixture model to be used. This value can be either "nn" (default) or "nn_bayes." 
}
  \item{batchFilePath}{
    Absolute or relative path to a tab-separated text file that indicates batch (and optionally, covariate information) for each sample. Optional.
}
  \item{verbose}{
     Whether to output more detailed status information as files are normalized. Default is TRUE.
}
}
%\details{
%}
\value{
  An \code{ExpressionSet} object that contains a row for each probeset/gene/transcript and a column for each input file. SCAN values will be on a log2 scale, centered at zero. UPC values will range between zero and one (lower values indicate that the gene is inactive and higher values indicating that the gene is active).
}
\references{
  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
  Piccolo SR, Withers MR, Francis OE, Bild AH and Johnson WE. Multi-platform single-sample estimates of transcriptional activation. Proceedings of the National Academy of Sciences of the United States of America, 2013, 110(44):11778-17783.
}
\author{
  Stephen R. Piccolo
}
\note{
If a Gene Expression Omnibus (GEO) identifier is specified for the \code{celFilePattern} parameter, an attempt will be made to download the sample(s) directly from GEO. If a study identifier (e.g., GSE22309) is specified, all CEL files from that study will be downloaded. If a sample identifier (e.g., GSM555237) is specified, only that sample will be downloaded.

By default, \code{SCAN} and \code{UPC} use the default mappings between probes and genes that have been provided by the manufacturer. However, these mappings may be outdated or may include problematic probes (for example, those that cross hybridize). The default mappings also may produce multiple summary values per gene. Alternative mappings, such as those provided by the BrainArray resource (see \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp}), allow \code{SCAN} and \code{UPC} to produce a single value per gene and to use updated gene definitions. Users can specify alternative mappings using the \code{probeSummaryPackage} parameter. If specified, this package must conform to the standards of the \code{AnnotationDbi} package. The BrainArray packages can be downloaded from \url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp}. When using BrainArray, be sure to download the R source package for probe-level mappings (see vignette for more information).

Because the SCAN/UPC algorithm accounts for nucleotide-level genomic composition across thousands of probes, it may take several minutes to normalize a sample, depending on the computer's processor speed and the type of microarray. To enable users to normalize samples in a shorter period of time, we have provided alternative functions called \code{SCANfast} and \code{UPCfast}. In this approach, a smaller number of probes is used for normalization, and a less stringent convergence threshold is used by default. We have found that microarrays processed with \code{SCANfast} (using default parameters) require 75\% less processing time (on average) but produce output values that correlate strongly (r = 0.998) with values produced by the \code{SCAN} function for the same arrays.

It is also possible to execute these functions in parallel. This approach uses the \code{foreach} package behind the scenes. If you have registered a parallel backend (for example, via the \code{doParallel} package), multiple CEL files can be processed in parallel. Otherwise, the files will be processed sequentially.

The \code{batchFilePath} parameter provides a convenient way to adjust the data for batch effects. It invokes the \code{ComBat} function within the \code{sva} package. Please see that package for additional details about how batch adjusting is performed. Batch adjusting will be performed after values have been SCAN normalized and summarized at the gene/probeset level. This is also true when \code{UPC} and \code{UPCfast} are being used---the data will be SCAN normalized and summarized, then batch adjusting will be performed, and lastly UPC transformation will occur. This process is different from when \code{UPC} or \code{UPCfast} are invoked without batch information; in this scenario, no SCAN normalization will occur.

The \code{modelType} parameter indicates which type of mixture model to use for UPC transformation. The "nn" model type has been used since the default implementation of this package. The "nn_bayes" model type is an experimental new approach intended for experiments where a subset of genes are expressed at extreme levels. When this model type is used, values will first be SCAN normalized and summarized at the gene/probeset level.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
# SCAN normalize a CEL file from GEO
normalized = SCAN("GSM555237")

# UPC normalize a CEL file from GEO
normalized = UPC("GSM555237")

# Normalize a CEL file and save output to a file
normalized = SCAN("GSM555237", "output_file.txt")

# Normalize a CEL file and summarize at the gene level using BrainArray
# mappings for Entrez Gene. First it is necessary to install the package
# and obtain the package name. For demonstration purposes, this file
# will be downloaded manually from GEO.
tmpDir = tempdir()
getGEOSuppFiles("GSM555237", makeDirectory=FALSE, baseDir=tmpDir)
celFilePath = file.path(tmpDir, "GSM555237.CEL.gz")
pkgName = InstallBrainArrayPackage(celFilePath, "17.0.1", "hs", "entrezg")
normalized = SCAN(celFilePath, probeSummaryPackage=pkgName)

# Normalize multiple files in parallel on multiple cores within a given
# computer. It is also possible using the doParallel package to spread
# the workload across multiple computers on a cluster.
library(doParallel)
registerDoParallel(cores=2)
result = SCAN("GSE22309")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
