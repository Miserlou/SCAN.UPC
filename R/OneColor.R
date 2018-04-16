##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

SCAN = function(celFilePattern, outFilePath=NA, convThreshold=0.01, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, batchFilePath=NA, verbose=TRUE)
{
  return(SCANprivate(celFilePattern=celFilePattern, outFilePath=outFilePath, intervalN=50000, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, exonArrayTarget=exonArrayTarget, batchFilePath=batchFilePath, verbose=verbose))
}

SCANfast = function(celFilePattern, outFilePath=NA, convThreshold=0.50, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, batchFilePath=NA, verbose=TRUE)
{
  return(SCANprivate(celFilePattern=celFilePattern, outFilePath=outFilePath, intervalN=10000, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, exonArrayTarget=exonArrayTarget, batchFilePath=batchFilePath, verbose=verbose))
}

SCANprivate = function(celFilePattern, outFilePath=NA, intervalN=50000, convThreshold=0.50, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, batchFilePath=NA, verbose=TRUE)
{
  if (is.na(batchFilePath))
    return(processCelFiles(celFilePattern=celFilePattern, outFilePath=outFilePath, intervalN=intervalN, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, UPC=FALSE, exonArrayTarget=exonArrayTarget, verbose=verbose))

  expressionSet = processCelFiles(celFilePattern=celFilePattern, outFilePath=NA, intervalN=intervalN,convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, UPC=FALSE, exonArrayTarget=exonArrayTarget, verbose=verbose)
  expressionSet = BatchAdjustFromFile(expressionSet, batchFilePath)

  if (!is.na(outFilePath))
    write.table(round(exprs(expressionSet), 8), outFilePath, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  return(expressionSet)
}

UPC = function(celFilePattern, outFilePath=NA, convThreshold=0.01, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, modelType="nn", batchFilePath=NA, verbose=TRUE)
{
  return(UPCprivate(celFilePattern=celFilePattern, outFilePath=outFilePath, intervalN=50000, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, exonArrayTarget=exonArrayTarget, modelType=modelType, batchFilePath=batchFilePath, verbose=verbose))
}

UPCfast = function(celFilePattern, outFilePath=NA, convThreshold=0.50, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, modelType="nn", batchFilePath=NA, verbose=TRUE)
{
  return(UPCprivate(celFilePattern=celFilePattern, outFilePath=outFilePath, intervalN=10000, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, exonArrayTarget=exonArrayTarget, modelType=modelType, batchFilePath=batchFilePath, verbose=verbose))
}

UPCprivate = function(celFilePattern, outFilePath=NA, intervalN=50000, convThreshold=0.01, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, exonArrayTarget=NA, modelType="nn", batchFilePath=NA, verbose=TRUE)
{
  if (!(modelType %in% c("nn", "nn2", "nn_bayes")))
    stop("The modelType parameter must be nn, nn2, or nn_bayes.")

  # Use the traditional UPC approach for Affymetrix arrays?
  if (is.na(batchFilePath) && modelType != "nn_bayes" && modelType != "nn2")
    return(processCelFiles(celFilePattern=celFilePattern, outFilePath=outFilePath, intervalN=intervalN, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, UPC=TRUE, exonArrayTarget=exonArrayTarget, verbose=verbose))

  # nn2 means that we want to do SCAN => UPC_Generic normal-normal UPC rather than the microarray version of normal-normal UPC
  if (modelType == "nn2")
    modelType = "nn"

  # SCAN normalize the data, don't save to file
  expressionSet = processCelFiles(celFilePattern=celFilePattern, outFilePath=NA, intervalN=intervalN, convThreshold=convThreshold, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, probeLevelOutDirPath=probeLevelOutDirPath, UPC=FALSE, exonArrayTarget=exonArrayTarget, verbose=verbose)

  # Adjust for batch
  if (!is.na(batchFilePath))
    expressionSet = BatchAdjustFromFile(expressionSet, batchFilePath)

  # UPC transform using generic model
  expressionSet = UPC_Generic_ExpressionSet(expressionSet, modelType=modelType, convThreshold=convThreshold, verbose=verbose)

  # Save to a file if needed
  if (!is.na(outFilePath))
    write.table(round(exprs(expressionSet), 8), outFilePath, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

  return(expressionSet)
}

processCelFiles = function(celFilePattern, outFilePath=NA, intervalN=50000, convThreshold=0.01, annotationPackageName=NA, probeSummaryPackage=NA, probeLevelOutDirPath=NA, UPC=FALSE, exonArrayTarget=NA, verbose=TRUE)
{
  if (convThreshold >= 1)
    stop("The convThreshold value must be lower than 1.0.")

  fromGEO = shouldDownloadFromGEO(celFilePattern)
  if (fromGEO)
    celFilePattern = downloadFromGEO(celFilePattern, expectedFileSuffixPattern="*.CEL.gz")

  fileNamePattern = sub("\\-", "\\\\-", glob2rx(basename(celFilePattern)))
  fileNamePattern = sub("\\+", "\\\\+", fileNamePattern)
  celFilePaths = list.files(path=dirname(celFilePattern), pattern=fileNamePattern, full.names=TRUE, ignore.case=TRUE)
  celFilePaths = unique(celFilePaths)

  if (length(celFilePaths) == 0)
    stop("No CEL files that match the pattern ", celFilePattern, " could be located.")

  if (!is.na(outFilePath))
    createOutputDir(dirPath=dirname(outFilePath), verbose=verbose)

  if (!is.na(probeLevelOutDirPath))
    createOutputDir(dirPath=probeLevelOutDirPath, verbose=verbose)

  celSummarizedList = foreach(celFilePath=celFilePaths) %dopar% {
    processCelFile(celFilePath=celFilePath, annotationPackageName=annotationPackageName, probeSummaryPackage=probeSummaryPackage, UPC=UPC, intervalN=intervalN, convThreshold=convThreshold, probeLevelOutDirPath=probeLevelOutDirPath, exonArrayTarget=exonArrayTarget, verbose=verbose)
  }

  if (fromGEO)
    unlink(celFilePaths)

  summarized = NULL

  for (i in 1:length(celFilePaths))
  {
    celSummarized = celSummarizedList[[i]]
    celFilePath = celFilePaths[i]

    if (is.null(celSummarized))
      next

    if (is.null(summarized))
    {
      summarized = celSummarized
      colnames(summarized) = basename(celFilePath)
    } else {
      previousColNames = colnames(summarized)
      summarized = merge(summarized, celSummarized, sort=FALSE, by=0)
      rownames(summarized) = summarized[,1]
      summarized = summarized[,-1]
      colnames(summarized) = c(previousColNames, basename(celFilePath))
    }
  }

  if (is.null(summarized))
    return(ExpressionSet())

  if (!is.na(outFilePath))
  {
    write.table(round(summarized, 8), outFilePath, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
    message("Saved output to ", outFilePath, sep="")
  }

  expressionSet = ExpressionSet(as.matrix(summarized))
  sampleNames(expressionSet) = colnames(summarized)
  featureNames(expressionSet) = rownames(summarized)

  return(expressionSet)
}

createOutputDir = function(dirPath, verbose=TRUE)
{
  if (!file.exists(dirPath))
  {
    if (!dir.create(dirPath))
      stop("An output file cannot be saved because a directory does not exist at ", dirPath, " and could not be created.")
  }
}

processCelFile = function(celFilePath, annotationPackageName, probeSummaryPackage, UPC, probeLevelOutDirPath, exonArrayTarget, intervalN, nbins=25, binsize=5000, convThreshold=0.01, verbose=TRUE)
{
#  probeLevelOutFilePath = NA
#  if (!is.na(probeLevelOutDirPath))
#    probeLevelOutFilePath = paste(probeLevelOutDirPath, "/", basename(celFilePath), ".txt", sep="")
#
#  if (!is.na(probeLevelOutFilePath) && file.exists(probeLevelOutFilePath))
#  {
#    data = read.table(probeLevelOutFilePath, sep="\t", stringsAsFactors=FALSE, header=FALSE, row.names=NULL)
#    xyCoord = as.character(data[,1])
#    y_norm = as.numeric(data[,2])
#    gam = as.numeric(data[,3])
#  } else
#  {
    if (is.na(annotationPackageName))
    {
      affyExpressionFS <- read.celfiles(celFilePath)
    } else {
      affyExpressionFS <- read.celfiles(celFilePath, pkgname=annotationPackageName)
    }

    annotationPackageName = affyExpressionFS@annotation

    if (is.na(exonArrayTarget) && (grepl("hugene", annotationPackageName)))
      #exonArrayTarget = "probeset"
      exonArrayTarget = "core"

    if (grepl("^pd\\.hta\\.2\\.0", annotationPackageName) | annotationPackageName == "pd.clariom.d.human") {
      data = getDataForHtaArray(celFilePath, affyExpressionFS, exonArrayTarget, verbose)
    } else {
      data = getDataForMostArrayTypes(celFilePath, affyExpressionFS, exonArrayTarget, verbose)
    }

    if (is.null(data))
      return(NULL)

    chunkSize = 100000
    chunkStartIndex = 1
    chunkEndIndex = chunkStartIndex + chunkSize - 1
    mxStartIndex = 1

    mx = matrix(nrow=nrow(data), ncol=80)

    while (chunkStartIndex <= nrow(data))
    {
      if (chunkEndIndex > nrow(data))
        chunkEndIndex = nrow(data)

      mxChunk = buildDesignMatrix(data[chunkStartIndex:chunkEndIndex,2], verbose=verbose)
      mx[mxStartIndex:(mxStartIndex + nrow(mxChunk) - 1),] = mxChunk

      chunkStartIndex = chunkStartIndex + chunkSize
      chunkEndIndex = chunkEndIndex + chunkSize
      mxStartIndex = mxStartIndex + nrow(mxChunk)
    }

    rawIntensities = as.numeric(data[,1])

    if (min(rawIntensities) == 0)
      rawIntensities = rawIntensities + 1

    my = log2(rawIntensities)
    nGroups = nrow(data) / binsize
    samplingProbeIndices = getSampleIndices(total=nrow(data), intervalN=intervalN, verbose=verbose)

    mixResult = EM_vMix(y=my[samplingProbeIndices], X=mx[samplingProbeIndices,], nbins=nbins, convThreshold=convThreshold, verbose=verbose, demo=length(grep("Vignette_Example", basename(celFilePath))) > 0)

    m1 = mx %*% mixResult$b1
    m2 = mx %*% mixResult$b2

    index = order(m1)
    y_norm = rep(0, length(my))
    for (i in 1:nGroups)
    {
      tmp = index[(binsize * i):min(binsize * i + binsize, length(my))]
      tmpSd = sig(y=my[tmp], m=m1[tmp], verbose=verbose)
      y_norm[tmp] = (my[tmp] - m1[tmp]) / as.vector(tmpSd)
    }

    bin = assign_bin(y=m1, nbins=nbins, verbose=verbose)
    gam = vresp(y=my, X=mx, bin=bin, p=mixResult$p, b1=mixResult$b1, s1=mixResult$s1, b2=mixResult$b2, s2=mixResult$s2, verbose=verbose)[,2]

    y_norm = round(y_norm, 8)
    gam = round(gam, 8)

#    if (!is.na(probeLevelOutFilePath))
#    {
#      message("Outputting probe-level values to ", probeLevelOutFilePath)
#      normOutput = cbind(data[,3], y_norm, gam, round(my, 8), round(m1, 8), round(m2, 8))
#      write.table(normOutput, probeLevelOutFilePath, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
#    }
#  }

  if (UPC)
  {
    data[,1] = gam
  } else {
    data[,1] = y_norm
  }

  if (!any(is.na(probeSummaryPackage)) && probeSummaryPackage != "NA")
  {
    if (is.character(probeSummaryPackage))
    {
      message(paste("Attempting to load", probeSummaryPackage))
      library(probeSummaryPackage, character.only=TRUE)
      probeSummaryPackage = eval(as.name(probeSummaryPackage))
    }

    probeSummaryData = as.data.frame(probeSummaryPackage, stringsAsFactors=FALSE)
    probeSummaryData = cbind(probeSummaryData, paste(probeSummaryData$x, probeSummaryData$y, sep="_"))

    data = as.matrix(merge(data, probeSummaryData, by.x=3, by.y=ncol(probeSummaryData), sort=FALSE)[,c(2,3,1,8)])
  }

  # A trim value of 0.3 seems to work better than 0.1 for HuGene 1.0 arrays (possibly others)
  return(as.matrix(tapply(as.numeric(data[,1]), data[,4], FUN=mean, trim=0.1)))
}

getDataForMostArrayTypes = function(celFilePath, affyExpressionFS, exonArrayTarget, verbose)
{
  if (is.na(exonArrayTarget))
  {
    probeInfo = oligo::getProbeInfo(affyExpressionFS, field=c("x", "y"), probeType="pm")
  } else {
    probeInfo = oligo::getProbeInfo(affyExpressionFS, field=c("x", "y"), probeType="pm", target=exonArrayTarget)
  }

  xyCoord = paste(probeInfo$x, probeInfo$y, sep="_")

  if (is.na(exonArrayTarget))
  {
    pint = oligo::pm(affyExpressionFS)
  } else {
    pint = oligo::pm(affyExpressionFS, target=exonArrayTarget)
  }

  if ((sum(pint==0) / length(pint)) > 0.01)
  {
    message(paste(celFilePath, " has a disproportionate number of zero values, so it cannot be processed.", sep=""))
    return(NULL)
  }

  if (is.na(exonArrayTarget))
  {
    pmSeq = pmSequence(affyExpressionFS)
    keepIndices = which(width(pmSeq) == 25)
    data = cbind(pint[keepIndices], as.character(pmSeq)[keepIndices], xyCoord[keepIndices], probeNames(affyExpressionFS)[keepIndices])
  } else {
    pmSeq = pmSequence(affyExpressionFS, target=exonArrayTarget)
    keepIndices = which(width(pmSeq) == 25)
    data = cbind(pint[keepIndices], as.character(pmSeq)[keepIndices], xyCoord[keepIndices], probeNames(affyExpressionFS, target=exonArrayTarget)[keepIndices])
  }

  return(data)
}

getDataForHtaArray = function(celFilePath, affyExpressionFS, exonArrayTarget, verbose)
{
  if (is.na(exonArrayTarget))
    exonArrayTarget = "probeset"

  if (exonArrayTarget != "probeset")
    stop("Currently, 'probeset' is the only allowed setting for the exonArrayTarget parameter for the Affymetrix HTA 2.0 arrays.")

  if (verbose)
    message(paste("Retrieving probe information for ", celFilePath, sep=""))
  probeInfo = oligo::getProbeInfo(affyExpressionFS, field=c("x", "y"), probeType="pm", target="probeset")

  if (verbose)
    message(paste("Retrieving intensity values for ", celFilePath, sep=""))
  pint = oligo::pm(affyExpressionFS, target="probeset")

  if ((sum(pint==0) / length(pint)) > 0.01)
  {
    message(paste(celFilePath, " has a disproportionate number of zero values, so it cannot be processed.", sep=""))
    return(NULL)
  }

  if (verbose)
    message(paste("Retrieving sequence information for ", celFilePath, sep=""))
  pmSeq = pmSequence(affyExpressionFS)

  if (verbose)
    message(paste("Merging data for ", celFilePath, sep=""))
  data = cbind(pint, as.character(pmSeq), paste(probeInfo$x, probeInfo$y, sep="_"), probeInfo$man_fsetid)

  #if (verbose)
  #  message(paste("Filtering data for ", celFilePath, sep=""))
  #data = data[which(nchar(data[,2]) == 25),]

  #pmIndices = pmindex(affyExpressionFS, target=exonArrayTarget)
  #data = data[pmIndices, ]

  return(data)
}

getSampleIndices = function(total, intervalN, verbose=TRUE)
{
  interval = floor(total / intervalN)
  if (interval <= 1)
    interval = 1

  seq(1, total, interval)
}

assign_bin = function(y, nbins, verbose=TRUE, randomSeed=1)
{
  quans = sort(y)[floor(length(y) * 1:nbins / nbins)]
  bins = sapply(y, function(x) { sum(x>quans) }) + 1

  if (length(table(bins)) != nbins)
  {
    if (verbose)
      message("The values were not separated into enough bins, so a tiny amount of noise will be added to make this possible.")

    set.seed(randomSeed)
    noise = rnorm(length(y)) / 10000000
    y = y + noise
    bins = assign_bin(y, nbins, verbose, randomSeed+1)
  }

  bins
}

buildDesignMatrix = function(seqs, verbose=TRUE)
{
  mx = sequenceDesignMatrix(seqs)
  numA = apply(mx[,1:25], 1, sum)
  numC = apply(mx[,26:50], 1, sum)
  numG = apply(mx[,51:75], 1, sum)
  numT = 25 - (numA + numC + numG)

  mx = cbind(numT, mx, as.integer(numA^2), as.integer(numC^2), as.integer(numG^2), as.integer(numT^2))
  mx = apply(mx, 2, as.integer)
}

EM_vMix = function(y, X, nbins, convThreshold=.01, verbose=TRUE, maxIt=100, demo=FALSE)
{
  if (verbose)
    message("Starting EM")

  quan = sort(y)[floor(0.5 * length(y)) - 1]
  gam = cbind(as.integer(y <= quan), as.integer(y > quan))

  p = apply(gam, 2, mean)

  b1 = mybeta(y=y, X=X, gam=gam[,1], verbose=verbose)
  b2 = mybeta(y=y, X=X, gam=gam[,2], verbose=verbose)
  bin = assign_bin(y=y, nbins=nbins, verbose=verbose&TRUE)
  s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
  s2 = vsig(y=y, X=X, b=b2, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)

  theta_old=c(p, b1, s1, b2, s2)

  it = 0
  conv = 1

  while (conv > convThreshold & it < maxIt)
  {
    # Expectation Step:
    gam = vresp(y=y, X=X, bin=bin, p=p, b1=b1, s1=s1, b2=b2, s2=s2, verbose=verbose)

    #M-Step
    p = apply(gam, 2, mean)
    b1 = vbeta(y=y, X=X, bin=bin, gam=gam[,1], s2=s1, prof=TRUE, verbose=verbose)
    bin = assign_bin(y=(X %*% b1), nbins=nbins, verbose=FALSE)
    b2 = vbeta(y=y, X=X, bin=bin, gam=gam[,2], s2=s2, prof=FALSE, verbose=verbose)
    s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
    s2 = vsig(y=y, X=X, b=b2, gam=gam[,2], bin=bin, nbins=nbins, verbose=verbose)

    theta = c(p, b1, s1, b2, s2)
    conv = max(abs(theta - theta_old) / theta_old)
    theta_old = theta
    it = it + 1

    if (verbose)
      message("Attempting to converge...iteration ", it, ", c = ", round(conv, 6))

    if (demo & (it == 3))
      break
  }

  if (verbose)
  {
    if (demo)
    {
      message("Convergence process halted early while in demo mode.")
    } else
    {
      if (it == maxIt)
      {
        message("Reached convergence limit...", it, " iterations. Proportion of background probes: ", round(p[1], 6))
      } else {
        message("Converged in ", it, " iterations. Proportion of background probes: ", round(p[1], 6))
      }
    }
  }

  list(p=p, b1=b1, b2=b2, s1=s1, s2=s2, bin=bin)
}

mybeta = function(y, X, gam, verbose=TRUE)
{
  sqgam = sqrt(gam)
  Xw = sqgam * X
  yw = sqgam * y

  z = t(Xw) %*% Xw
  a = solve(z)

  b = a %*% t(Xw)
  as.numeric(b %*% yw)
}

vbeta = function(y, X, bin, gam, s2, prof, verbose=TRUE)
{
  vars = sqrt(s2[bin])
  sqgam = sqrt(gam)
  vars_sqgam = vars * sqgam

  Xw = 1 / vars * sqgam * X
  yw = 1 / vars * sqgam * y

  tXw = t(Xw)
  tXwXw = tXw %*% Xw
  stXwXw = solve(tXwXw)
  stXwXwtXw = stXwXw %*% tXw
  result = stXwXwtXw %*% yw

  result
}

vresp = function(y, X, bin, p, b1, s1, b2, s2, verbose=TRUE)
{
  vars0 = s1[bin]
  L0 = dn(y=y, m=(X %*% b1), s2=vars0, verbose=verbose)
  vars1 = s2[bin]
  L1 = dn(y=y, m=(X %*% b2), s2=vars1, verbose=verbose)

  gam1 = p[1] * L0 / (p[1] * L0 + p[2] * L1)
  gam2 = 1 - gam1
  cbind(gam1, gam2)
}

dn = function(y, m, s2, verbose=TRUE)
{
  1 / (sqrt(2 * pi * s2)) * exp(-1 / (2 * s2) * (y - m)^2)
}

vsig = function(y, X, b, gam, bin, nbins, verbose=TRUE)
{
  s2 = NULL

  for (i in 1:nbins)
  {
    ystar = y[bin==i]
    Xstar = X[bin==i,]
    gamstar = gam[bin==i] + .01
    resid = as.numeric(ystar - Xstar %*% b)

    s2 = c(s2, ((resid * gamstar) %*% resid) / sum(gamstar))
  }

  s2
}

sig = function(y, m, verbose=TRUE)
{
  resid = y - m
  sqrt((resid %*% resid) / length(y))
}

#EM_vMix_bayes = function(y, X, nbins, convThreshold=.01, verbose=TRUE, maxIt=100, demo=FALSE,empirical=T, lambda=1, exprProp=0.4, k=1)
#{
#  if (verbose)
#    message("Starting EM")
#  
#  quan = sort(y)[floor(0.5 * length(y)) - 1]
#  gam <<- cbind(as.integer(y <= quan), as.integer(y > quan))
#  
#  p = apply(gam, 2, mean)
#  
#  b1 = mybeta(y=y, X=X, gam=gam[,1], verbose=verbose)
#  X1 = X[y <= quan,] #
#  b2 = mybeta(y=y, X=X, gam=gam[,2], verbose=verbose)
#  X2 = X[y > quan, ]  #
#  bin = assign_bin(y=y, nbins=nbins, verbose=verbose)
#  s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
#  s2 = vsig(y=y, X=X, b=b2, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
#  
#  n_total <- length(y)  #
#  beta_a <- n_total*0.1*exprProp; beta_b <- n_total*0.1*(1-exprProp)  #
#  
#  # prior covariance of beta1 and beta2
#  if (empirical==T){
#    cov_beta <- s2[bin] * solve(t(X1)%*%X1) * n_total/k 
#    cov_beta_inv <- solve(cov_beta)
#    #print(cov_beta_inv)
#    
#  } else {
#    cov_beta_inv <- diag(rep(lambda,NCOL(X)))
#  }
#  
#  theta_old=c(p, b1, s1, b2, s2)
#  
#  it = 0
#  conv = 1
# 
#  while (conv > convThreshold & it < maxIt)
#  {
#    # Expectation Step:
#    gam = vresp(y=y, X=X, bin=bin, p=p, b1=b1, s1=s1, b2=b2, s2=s2, verbose=verbose)
#    
#    #M-Step
#    #p = apply(gam, 2, mean)
#    p = (beta_a + colSums(gam[,2])) / (n_total+beta_a+beta_b) #
#    
#    b1 = vbeta_bayes(y=y, X=X, bin=bin, gam=gam[,1], s2=s1, prof=TRUE, verbose=verbose) #
#    bin = assign_bin(y=(X %*% b1), nbins=nbins, verbose=verbose)
#    b2 = vbeta_bayes(y=y, X=X, bin=bin, gam=gam[,2], s2=s2, prof=FALSE, verbose=verbose)  #
#    s1 = vsig(y=y, X=X, b=b1, gam=gam[,1], bin=bin, nbins=nbins, verbose=verbose)
#    s2 = vsig(y=y, X=X, b=b2, gam=gam[,2], bin=bin, nbins=nbins, verbose=verbose)
#    
#    theta = c(p, b1, s1, b2, s2)
#    conv = max(abs(theta - theta_old) / theta_old)
#    theta_old = theta
#    it = it + 1
#    
#    if (verbose)
#      message("Attempting to converge...iteration ", it, ", c = ", round(conv, 6))
#    
#    if (demo & (it == 3))
#      break
#  }
#  
#  if (verbose)
#  {
#    if (demo)
#    {
#      message("Convergence process halted early while in demo mode.")
#    } else
#    {
#      if (it == maxIt)
#      {
#        message("Reached convergence limit...", it, " iterations. Proportion of background probes: ", round(p[1], 6))
#      } else {
#        message("Converged in ", it, " iterations. Proportion of background probes: ", round(p[1], 6))
#      }
#    }
#  }
#  
#  list(p=p, b1=b1, b2=b2, s1=s1, s2=s2, bin=bin)
#}
#
#vbeta_bayes = function(y, X, bin, gam, s2, prof, verbose=TRUE,cov_beta_inv)
#{
#  vars = sqrt(s2[bin]) ##?why sqrt? s2 is a variance vector with nbin number of elements (need to confirm with Steve)
#  sqgam = sqrt(gam)
#  vars_sqgam = vars * sqgam
#  
#  Xw = 1 / vars * sqgam * X # ??
#  yw = 1 / vars * sqgam * y
#  
#  #   tXw = t(Xw)
#  #   tXwXw = tXw %*% Xw
#  #   stXwXw = solve(tXwXw)
#  #   stXwXwtXw = stXwXw %*% tXw
#  #   result = stXwXwtXw %*% yw
#  
#  result = solve(t(Xw)%*%diag(gam)%*%X+s2[bin]*cov_beta_inv) %*% t(Xw)%*%diag(gam)%*%y
#  
#  result
#}
