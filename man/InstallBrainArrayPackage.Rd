\name{InstallBrainArrayPackage}
\alias{InstallBrainArrayPackage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Helper function for installing BrainArray packages
}
\description{
  When processing Affymetrix microarrays, users can specify alternative probe/gene mappings via the probeSummaryPackage parameter. Users can download such packages directly from the BrainArray web site and install them manually. Or they can use this helper function to download and install them in a single step.
}
\usage{
InstallBrainArrayPackage(celFilePath, version, organism, annotationSource)
}
\arguments{
  \item{celFilePath}{
    Path to an example CEL file. The Affymetrix version name will be extracted from this file.
}
  \item{version}{
    Version of BrainArray (for example, "17.1.0") you want to download. 
}
  \item{organism}{
    Abbreviated name (for example, "hs" = homo sapiens, "mm" = mus musculus, "gg" = gallus gallus) of the organism for which the microarrays are designed.
}
  \item{annotationSource}{
    Abbreviated name of the annotation source (for example, "entrezg" = Entrez Gene, "ensg" = Ensembl Gene) you want to use. 
}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A character object that indicates the name of the installed package.
}
%\references{
%  Piccolo SR, Sun Y, Campbell JD, Lenburg ME, Bild AH, and Johnson WE. A single-sample microarray normalization method to facilitate personalized-medicine workflows. \emph{Genomics}, 2012, 100:6, pp. 337-344.
%}
\author{
  Stephen R. Piccolo
}
\note{
  Information about BrainArray versions, organism names, and annotation sources can be obtained via the BrainArray web site (\url{http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/genomic_curated_CDF.asp}).
}

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
\dontrun{
pkgName = InstallBrainArrayPackage(celFilePath, "17.0.1", "hs", "entrezg")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
