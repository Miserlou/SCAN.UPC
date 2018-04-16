### R code from vignette source 'SCAN.vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
## do work in temporary directory
pwd <- setwd(tempdir())


###################################################
### code chunk number 2: download-geo-direct (eval = FALSE)
###################################################
## library(GEOquery)
## tmpDir = tempdir()
## library(GEOquery)
## getGEOSuppFiles("GSM555237", makeDirectory=FALSE, baseDir=tmpDir)
## celFilePath = file.path(tmpDir, "GSM555237.CEL.gz")


###################################################
### code chunk number 3: download-normalize (eval = FALSE)
###################################################
## library(SCAN.UPC)
## normalized = SCAN(celFilePath)


###################################################
### code chunk number 4: scan-geo (eval = FALSE)
###################################################
## normalized = SCAN("GSM555237")


###################################################
### code chunk number 5: download-normalize2 (eval = FALSE)
###################################################
## normalized = SCAN(celFilePath, outFilePath="output_file.txt")


###################################################
### code chunk number 6: download-brainarray (eval = FALSE)
###################################################
## install.packages("hgu95ahsentrezgprobe_15.0.0.tar.gz", repos=NULL, type="source")


###################################################
### code chunk number 7: install-brainarray (eval = FALSE)
###################################################
## pkgName = InstallBrainArrayPackage(celFilePath, "15.0.0", "hs", "entrezg")


###################################################
### code chunk number 8: scan-brainarray (eval = FALSE)
###################################################
## normalized = SCAN(celFilePath, probeSummaryPackage=pkgName)


###################################################
### code chunk number 9: scan-twocolor-main (eval = FALSE)
###################################################
## SCAN_TwoColor("GSM1072833", "output_file.txt")


###################################################
### code chunk number 10: upc-microarray (eval = FALSE)
###################################################
## upc1 = UPC("GSM555237")
## 
## upc2 = UPC_TwoColor("GSM1072833")


###################################################
### code chunk number 11: upc-rnaseq (eval = FALSE)
###################################################
## upc3 = UPC_RNASeq("ReadCounts.txt", "Annotation.txt")


###################################################
### code chunk number 12: scan-parallel (eval = FALSE)
###################################################
## library(doParallel)
## registerDoParallel(cores=2)
## result = SCAN("GSE22309")


