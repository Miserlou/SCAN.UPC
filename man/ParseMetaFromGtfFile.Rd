\name{ParseMetaFromGtfFile}
\alias{ParseMetaFromGtfFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Helper function to parse length and GC content information from a GTF file.
}
\description{
  When applying the \code{\link{UPC_RNASeq}} function, it is possible to correct for the length and GC content of genomic features. To accomplish this, an annotation file indicating these values for each feature must be provided. This helper function enables users to generate an annotation file, using a GTF file and genome FASTA file as references.
}
\usage{
ParseMetaFromGtfFile(gtfFilePath, fastaFilePattern, outFilePath,
                     featureTypes=c("protein_coding"), attributeType="gene_id")
}
\arguments{
  \item{gtfFilePath}{
    Path to the GTF file that will be parsed.
}
  \item{fastaFilePattern}{
    File pattern that indicates where FASTA file(s) for the associated reference genome can be found.
}
  \item{outFilePath}{
   Path where the output file will be stored.
}
  \item{featureTypes}{
   One or more feature types (for example, "protein_coding," "unprocessed_pseudogene") that should be extracted from the GTF file. The default is "protein_coding."
}
  \item{attributeType}{
  The type of attribute ("gene_id", "transcript_id") to be parsed. Values will be grouped according to these attributes.
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%  An \code{ExpressionSet} object that contains a row for each probeset/gene/transcript and a column for each input file.
%}
%\references{
%  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
%}
\author{
  Stephen R. Piccolo
}
%\note{
%}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
ParseMetaFromGtfFile("GRCh37_XY.gtf", "GRCh37.fa", "GRCh37_Annotation.txt")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
