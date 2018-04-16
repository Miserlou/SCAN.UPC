##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

UPC_RNASeq = function(inFilePattern, annotationFilePath=NA, outFilePath=NA, modelType="nn", convThreshold=0.01, ignoreZeroes=FALSE, numDataHeaderRows=0, numAnnotationHeaderRows=0, batchFilePath=NA, verbose=TRUE)
{
  data = readFilesIntoMatrix(inFilePattern, numHeaderRows=numDataHeaderRows)

  if (!is.na(batchFilePath))
  {
    # Add a tiny bit of noise to the actual values so ComBat doesn't throw an error
    set.seed(1)
    tinyRandomNoise = matrix(runif(length(data)) / 1000000, ncol=ncol(data))

    # Create an ExpressionSet
    expressionSet = ExpressionSet(log2(data + tinyRandomNoise))
    sampleNames(expressionSet) = colnames(data)

    expressionSet = BatchAdjustFromFile(expressionSet, batchFilePath)
    data = 2^exprs(expressionSet)
  }

  annotationData = getAnnotationData(annotationFilePath, numAnnotationHeaderRows)
  lengths = NULL
  gcContent = NULL

  if (!any(is.na(annotationData)))
  {
    commonRowNames = intersect(rownames(data), rownames(annotationData))

    if (length(commonRowNames) < (nrow(data) * 0.5))
      stop("Less than half of the annotations overlap with the data values.")

    data = data[commonRowNames, , drop=FALSE]
    annotationData = annotationData[commonRowNames,]
    lengths = as.numeric(annotationData[,1])
    gcContent = as.numeric(annotationData[,2]) / lengths
  }

  upcData = apply(data, 2, function(x) {
    UPC_RNASeq_Single(x, rownames(data), lengths=lengths, gcContent=gcContent, modelType=modelType, convThreshold=convThreshold, ignoreZeroes=ignoreZeroes, verbose=TRUE)
  })

  if (!is.na(outFilePath))
  {
    if (verbose)
      message(paste("Saving results to", outFilePath))

    write.table(upcData, outFilePath, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  }

  expressionSet = ExpressionSet(upcData)
  featureNames(expressionSet) = rownames(upcData)
  sampleNames = colnames(upcData)

  return(expressionSet)
}

UPC_RNASeq_Single = function(expressionValues, featureNames, lengths=NULL, gcContent=NULL, modelType="nn", convThreshold=0.01, ignoreZeroes=FALSE, verbose=TRUE)
{
  data = cbind(expressionValues, lengths, gcContent)
  rownames(data) = featureNames

  zeroData = NULL
  if (ignoreZeroes)
  {
    zeroData = data[which(expressionValues == 0),, drop=FALSE]
    data = data[which(expressionValues != 0),, drop=FALSE]
  }

  naData = NULL
  if (any(is.na(as.numeric(data[,1]))))
  {
    naData = data[which(is.na(as.numeric(data[,1]))),, drop=FALSE]
    data = data[which(!is.na(as.numeric(expressionValues))),, drop=FALSE]
  }

  lengths = NULL
  if (ncol(data) > 1)
    lengths = as.numeric(data[,2])

  gcContent = NULL
  if (ncol(data) > 2)
    gcContent = as.numeric(data[,3])

  upc = round(UPC_Generic(as.numeric(data[,1]), lengths=lengths, gcContent=gcContent, modelType=modelType, convThreshold=convThreshold), 6)

  names(upc) = rownames(data)
  upc = c(upc, zeroData[,1], naData[,1])
  upc = upc[featureNames]

  return(upc)
}

getAnnotationData = function(annotationFilePath, numAnnotationHeaderRows)
{
  if (is.na(annotationFilePath))
  {
    annotationData = NA
  } else {
    if (!(file.exists(annotationFilePath)))
      stop(paste("No annotation file exists at ", annotationFilePath, ".", sep=""))

    message("Importing annotations...")
    annotationData = as.matrix(read.table(annotationFilePath, sep="\t", header=FALSE, stringsAsFactors=FALSE, row.names=1, quote="\"", check.names=FALSE, skip=numAnnotationHeaderRows))
  }

  return(annotationData)
}
