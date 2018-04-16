\name{UPC_RNASeq}
\alias{UPC_RNASeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Universal exPression Codes (UPC) for RNA-Seq data
}
\description{
This function is used to derive UPC values for RNA-Seq data. It requires at least one input file that specifies a read count for each genomic region (e.g., gene). This file should list a unique identifier for each region in the first column and corresponding read counts (not RPKM/FPKM values) in the second column.

This function also can correct for the GC content and length of each genomic region. Users who wish to enable this correction must provide a separate annotation file. This tab-separated file should contain a row for each genomic region. The first column should contain a unique identifier that corresponds to identifiers from the read-count input file. The second column should indicate the length of the genomic region. And the third column should specify the number of G or C bases in the region. The \code{\link{ParseMetaFromGtfFile}} function can be used to generate annotation files.
}
\usage{
UPC_RNASeq(inFilePattern, annotationFilePath = NA, outFilePath = NA, modelType = "nn",
  convThreshold = 0.01, ignoreZeroes = FALSE, numDataHeaderRows=0,
  numAnnotationHeaderRows=0, batchFilePath=NA, verbose = TRUE)
}
\arguments{
  \item{inFilePattern}{
     Absolute or relative path to the input file(s) to be processed. The input file(s) can contain one or more columns, where each column would contain data for a given sample. To process multiple files, wildcard characters can be used (e.g., "*.txt"). Required.
}
  \item{annotationFilePath}{
     Absolute or relative path where the annotation file is located. This parameter is optional.
}
  \item{outFilePath}{
     Absolute or relative path where the output file will be saved. This is optional.
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
  \item{numDataHeaderRows}{
    The number of header rows present in the input data file(s). If a header is present, the column names will be used as sample IDs.
}
  \item{numAnnotationHeaderRows}{
    The number of header rows present in the annotation data file (if one has been specified).
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
  An \code{ExpressionSet} object that contains a row for each probeset/gene/transcript and a column for each input file.
}
\references{
  Piccolo SR, Withers MR, Francis OE, Bild AH and Johnson WE. Multi-platform single-sample estimates of transcriptional activation. Proceedings of the National Academy of Sciences of the United States of America, 2013, 110:44 17778-17783.
}
\author{
  Stephen R. Piccolo
}
\note{
  RNA-Seq data by nature have a lot of zero read counts. Samples with an excessive number of zeroes may lead to error messages because genes cannot be allocated properly to bins. The user can specify \code{ignoreZeroes=TRUE} to avoid this error. In practice, we have seen that the resulting UPC values are similar with either approach.

  The \code{batchFilePath} parameter provides a convenient way to adjust the data for batch effects. It invokes the \code{ComBat} function within the \code{sva} package. Please see that package for additional details about how batch adjusting is performed. Batch adjusting is performed before the values are UPC transformed.

  The \code{modelType} parameter indicates which type of mixture model to use for UPC transformation. The "nn_bayes" model type is an experimental new approach intended for experiments where a subset of genes are expressed at extreme levels.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
result = UPC_RNASeq("ReadCounts.txt", "Annotation.txt")
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
