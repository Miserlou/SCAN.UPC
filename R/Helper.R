##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

shouldDownloadFromGEO = function(inFilePattern)
{
  substr(inFilePattern, 1, 3) %in% c("GSE", "GSM") & length(grep("\\.", inFilePattern)) == 0 & length(grep("\\*", inFilePattern)) == 0
}

downloadFromGEO = function(inFilePattern, expectedFilePrefixPattern="*", expectedFileSuffixPattern="*")
{
  tmpDir=tempdir()
  message(paste("Downloading ", inFilePattern, " directly from GEO to ", tmpDir, ".", sep=""))

  downloadDir=paste(tmpDir, "/", inFilePattern, sep="")
  unlink(downloadDir, recursive=TRUE)

  getGEOSuppFiles(inFilePattern, makeDirectory=FALSE, baseDir=tmpDir)

  if (substr(inFilePattern, 1, 3) == "GSE")
  {
    tarFilePath = file.path(tmpDir, paste(inFilePattern, "_RAW.tar", sep=""))

    if (!file.exists(tarFilePath))
      stop(paste("No raw data files could be downloaded from GEO for ", inFilePattern, sep=""))

    dir.create(downloadDir, recursive=TRUE)
    untar(tarFilePath, exdir=downloadDir)
    inFilePattern = file.path(downloadDir, expectedFilePrefixPattern, sep="")
  }

  if (substr(inFilePattern, 1, 3) == "GSM")
  {
    downloadedFiles = list.files(path=tmpDir, full.names=TRUE, pattern=glob2rx(expectedFileSuffixPattern), ignore.case=TRUE)

    if (length(downloadedFiles) == 0)
      stop(paste("No raw data files could be downloaded from GEO for ", inFilePattern, sep=""))

    inFilePattern = file.path(tmpDir, basename(downloadedFiles))
  }

  inFilePattern
}

downloadBeadChipFromGEO = function(inFilePattern)
{
  tmpDir=tempdir()
  downloadDir=paste(tmpDir, "/", inFilePattern, sep="")

  message(paste("Downloading ", inFilePattern, " directly from GEO to ", downloadDir, ".", sep=""))

  unlink(downloadDir, recursive=TRUE)
  dir.create(downloadDir, recursive=TRUE)

  fileDataFrame = getGEOSuppFiles(inFilePattern, makeDirectory=FALSE, baseDir=downloadDir)

  probeDataFilePath = rownames(fileDataFrame)[grepl("_non-normalized\\.txt\\.gz$", rownames(fileDataFrame))]
  tarFilePath = rownames(fileDataFrame)[grepl("_RAW\\.tar$", rownames(fileDataFrame))]

  if (length(tarFilePath) > 0)
  {
    untar(tarFilePath, exdir=downloadDir)
    tarFileListFilePath = rownames(fileDataFrame)[grepl("filelist\\.txt$", rownames(fileDataFrame))]
    tarFileList = read.table(tarFileListFilePath, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=NULL, check.names=FALSE, comment.char="")

  #wgdaslFilePath = rownames(fileDataFrame)[grepl("WGDASL\\.txt\\.gz$", rownames(fileDataFrame))]
#print(wgdaslFilePath)
#  bgxFilePath = paste(downloadDir, fileList[which(fileList$Type=="BGX"),]$Name, sep="/")
#  idatFilePaths = paste(downloadDir, fileList[which(fileList$Type=="IDAT"),]$Name, sep="/")
  }

##########
# See what tarFileList has when it is a WGDASL file rather than a BGX
##########
print(fileDataFrame)
print(tarFileList)
print(probeDataFilePath)
print(tarFilePath)

stop()

#  return(list(bgxFilePath=bgxFilePath, idatFilePaths=idatFilePaths))
}

InstallBrainArrayPackage = function(celFilePath, version, organism, annotationSource)
{
  platform = cleancdfname(read.celfile.header(celFilePath, info="full")$cdfName)
  platform = sub("cdf", "", platform)
  platform = sub("stv1", "st", platform)
  platform = sub("stv2", "st", platform)

  packageName=paste(platform, organism, annotationSource, "probe", sep="")

  #if (!(packageName %in% rownames(installed.packages())))
  tmpDir = tempdir()

  packageFileName=paste(packageName, "_", version, ".tar.gz", sep="")
  tempPackageFilePath = paste(tmpDir, packageFileName, sep="")
  #packageUrl = paste("http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/", version, "/", annotationSource, ".download/", packageFileName, sep="")
  #http://mbni.org/customcdf/19.0.0/entrezg.download/hgu133ahsentrezgprobe_19.0.0.tar.gz
  packageUrl = paste("http://mbni.org/customcdf/", version, "/", annotationSource, ".download/", packageFileName, sep="")

  download.file(packageUrl, tempPackageFilePath)
  install.packages(tempPackageFilePath, repos=NULL, type="source")

  return(packageName)
}

ParseMetaFromGtfFile = function(gtfFilePath, fastaFilePattern, outFilePath, featureTypes=c("protein_coding"), attributeType="gene_id")
{
  fastaFilePaths = list.files(path=dirname(fastaFilePattern), pattern=glob2rx(basename(fastaFilePattern)), full.names=TRUE)

  if (length(fastaFilePaths) == 0)
    stop(paste("No FASTA file could be found at", fastaFilePattern))

  fastaStrings = readDNAStringSet(fastaFilePaths)

  # This is in case the FASTA description has more than the chromosome name (delimiter may be a space or |)
  names(fastaStrings) = gsub("\\|", " ", names(fastaStrings))
  names(fastaStrings) = sapply(names(fastaStrings), function(x) { strsplit(x, " ")[[1]][1] })

  message("Saving GTF data to temporary files")
  tmpDir = tempdir()
  chromosomes = NULL

  f = file(gtfFilePath, "r")
  stop = FALSE
  lineCount = 0

  while(!stop)
  {
    lineCount = lineCount + 1
    if (lineCount %% 10000 == 0)
      message(paste("Done parsing ", as.character(lineCount), " lines from ", gtfFilePath, sep=""))

    line = readLines(f, n = 1)

    if(length(line) == 0)
    {
      stop = TRUE
      next
    }

    lineItems = strsplit(line, "\t")[[1]]
    chromosome = lineItems[1]
    featureType = lineItems[2]
    subFeatureType = lineItems[3]

    if (chromosome %in% names(fastaStrings) && subFeatureType == "exon" && length(intersect(featureType, featureTypes)) > 0)
    {
      chromosomes = base::unique(c(chromosomes, chromosome))

      startPos = as.integer(lineItems[4])
      stopPos = as.integer(lineItems[5])
      attributes = strsplit(lineItems[9], ";")[[1]]
      attributes = sub("^ ", "", attributes)
      featureID = attributes[which(grepl(paste("^", attributeType, sep=""), attributes, ignore.case=TRUE))]
      featureID = sub(paste(attributeType, " ", sep=""), "", featureID)
      featureID = gsub("\"", "", featureID)

      outRow = c(featureID, startPos, stopPos)
      write.table(matrix(outRow, nrow=1), paste(tmpDir, "/", chromosome, sep=""), sep="\t", row.names=F, col.names=F, quote=F, append=T)
    }
  }

  close(f)

  chromosomes = sort(chromosomes)
  for (chromosome in chromosomes)
  {
    message(paste("Summarizing feature values for chromosome ", chromosome, sep=""))
    ProcessGtfSubset(paste(tmpDir, "/", chromosome, sep=""), fastaStrings, chromosome, outFilePath, chromosome!=chromosomes[1])
  }
}

ProcessGtfSubset = function(tmpFilePath, fastaStrings, chromosome, outFilePath, appendToOut)
{
  data = read.table(tmpFilePath, sep="\t", stringsAsFactors=F, header=F, row.names=NULL, check.names=F)

  featureCoordsList = list()

  f = file(tmpFilePath, "r")
  stop = FALSE

  while (!stop)
  {
    line = readLines(f, n = 1)

    if(length(line) == 0)
    {
      stop = TRUE
      next
    }

    lineItems = strsplit(line, "\t")[[1]]
    featureID = lineItems[1]
    startPos = as.integer(lineItems[2])
    stopPos = as.integer(lineItems[3])
    thisRange = IRanges(startPos, stopPos)

    if (is.null(featureCoordsList[[featureID]]))
    {
      featureCoordsList[[featureID]] = thisRange
    } else {
      featureCoordsList[[featureID]] = IRanges::union(featureCoordsList[[featureID]], thisRange)
    }
  }

  close(f)

  outData = NULL
  for (featureID in names(featureCoordsList))
  {
    featureData = as.matrix(as.data.frame(featureCoordsList[[featureID]]))
    featureData = cbind(featureData, NA)

    for (i in 1:nrow(featureData))
    {
      start = featureData[i,1]
      stop = featureData[i,2]

      subSequence = subseq(fastaStrings[chromosome], start, stop)
      gcCount = sum(letterFrequency(subSequence, c("C", "G")))

      featureData[i,4] = gcCount
    }

    if (nrow(featureData) == 1)
    {
      outRow = featureData[1,3:4]
    } else {
      outRow = apply(featureData[,3:4], 2, sum)
    }

    outData = rbind(outData, c(featureID, outRow))
  }

  write.table(outData, outFilePath, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=appendToOut)
}

BatchAdjust = function(expressionSet, batchVariableName, covariateVariableNames=c())
{
  for (variableName in c(batchVariableName, covariateVariableNames))
    if (!(variableName %in% varLabels(expressionSet)))
      stop(paste("No ", variableName, " variable exists in the phenotype data, so no batch adjustment can be performed.", sep=""))

  mod = NULL

  if (length(covariateVariableNames) > 0)
  {
    modelMatrixCommand = paste("model.matrix(~", paste(covariateVariableNames, collapse=" + "), ", data=pData(expressionSet))", sep="")

    message(paste("Covariate(s) have been specified for batch adjustment. The data will be adjusted using the following model matrix: ", modelMatrixCommand, sep=""))
    eval(parse(text=paste("mod = ", modelMatrixCommand, sep="")))
  } else {
    message("No batch covariates have been specified.")
  }

  message("Getting ready to perform batch adjustment.")
  exprs(expressionSet) = ComBat(exprs(expressionSet), batch=pData(expressionSet)[,batchVariableName], mod=mod)
  message("Done performing batch adjustment.")

  return(expressionSet)
}

BatchAdjustFromFile = function(expressionSet, batchFilePath)
{
  if (class(expressionSet) != "ExpressionSet")
    stop("The object is not an expressionSet.")

  if (!file.exists(batchFilePath))
  {
    warning(paste("No file exists at ", batchFilePath, ".", sep=""))
    return(expressionSet)
  }

  covariateData = read.table(batchFilePath, sep="\t", header=TRUE, row.names=NULL, check.names=FALSE)
  rownames(covariateData) = covariateData[,1]

  if (length(intersect(covariateData[,1], sampleNames(expressionSet))) != ncol(expressionSet))
    stop(paste("For one or more of the sampleIDs in the ExpressionSet object, there was no covariate information in ", batchFilePath, ".", sep=""))

  covariateData = covariateData[sampleNames(expressionSet), ]

  covariateNames = colnames(covariateData)[-1]
  newPhenotypeData = cbind(pData(expressionSet), covariateData[,-1])
  colnames(newPhenotypeData) = c(varLabels(expressionSet), covariateNames)

  pData(expressionSet) = newPhenotypeData

  return(BatchAdjust(expressionSet, batchVariableName=covariateNames[1], covariateVariableNames=covariateNames[-1]))
}

readFilesIntoMatrix = function(inFilePattern, numHeaderRows=0)
{
  inFilePaths = list.files(path=dirname(inFilePattern), pattern=glob2rx(basename(inFilePattern)), full.names=TRUE)

  if (length(inFilePaths) == 0)
    stop("No data files that match the pattern ", inFilePattern, " could be located.")

  header = FALSE
  if (numHeaderRows == 1)
  {
    header = TRUE
    numHeaderRows = numHeaderRows - 1
  }

  matrixData = NULL

  for (inFilePath in inFilePaths)
  {
    message(paste("Reading ", inFilePath, " into a matrix.", sep=""))
    inFileData = as.matrix(read.table(inFilePath, sep="\t", header=header, stringsAsFactors=FALSE, row.names=1, quote="\"", check.names=FALSE, skip=numHeaderRows))

    if (ncol(inFileData) == 1)
    {
      newColumnNames = basename(inFilePath)
    } else {
      if (header)
      {
        newColumnNames = colnames(inFileData)
      } else {
        newColumnNames = paste(basename(inFilePath), "_Sample_", 1:ncol(inFileData), sep="")
      }
    }

    if (is.null(matrixData))
    {
      previousColumnNames = NULL
      matrixData = inFileData
    } else {
      previousColumnNames = colnames(matrixData)

      matrixData = suppressWarnings(merge(matrixData, inFileData, by=0, sort=FALSE, all=TRUE))
      rownames(matrixData) = matrixData[,1]
      matrixData = as.matrix(matrixData[,-1])
    }

    colnames(matrixData) = c(previousColumnNames, newColumnNames)
  }

  return(matrixData)
}
