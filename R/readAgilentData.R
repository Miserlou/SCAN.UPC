##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

readAgilentData <-
function (data = 0, filename = "NULL", signal="median", useGC = FALSE, filter= FALSE,filter.cols='default')  ####  
{                                                                                     ####  
    if (filename != "NULL") {                                                         ####  Reads in data if filename is specified.
        data <- read.delim(filename, header = TRUE, sep = "\t",                       ####  'fill-TRUE' allows the rows to have unequal length.
            fill = TRUE, skip = 9)                                                    ####  'skip=9' ignores the first 9 rows of the dataset.
    }                                                                                 ####  
    if (signal == "median") {                                                         ####
    	cols <- match(c("ProbeName","SystematicName","rMedianSignal",
    		 "gMedianSignal", "Sequence"),                                            ####  Determines the column numbers of the desired data.
        	names(data))                                                              ####  
        }                                                                             ####  
    else if (signal == "processed") {                                                 ####  
    	cols <- match(c("ProbeName","SystematicName","rProcessedSignal",
    		 "gProcessedSignal", "Sequence"),                                         ####  Determines the column numbers of the desired data.
        	names(data))                                                              ####  
    	}                                                                             ####  
    else if (signal == "mean") {
    	cols <- match(c("ProbeName","SystematicName","rMeanSignal",
    		 "gMeanSignal", "Sequence"),                
        	names(data))    	
    	}
    else if (is.vector(signal == TRUE)) {
    	cols <- match(c("ProbeName","SystematicName",signal,"Sequence"),names(data))
    	if (length(cols)==0) stop("Signal columns not found",call=FALSE)
    	}
#    else {
#    	cols <- match(c("rMedianSignal", "gMedianSignal", "Sequence"),                ####  Determines the column numbers of the desired data.
#        	names(data))     	
#    	}
    nobs <- length(data[, 1])                                                         ####  Determines length for loop.
    gcCount <- matrix(NA, nobs, 1)                                                    ####  Initializes the gcCount matrix.
    if (useGC != FALSE) {                                                             ####  
        gc <- function(seq) {                                                         ####  
            a <- unlist(strsplit(seq, NULL))                                          ####  Decomposes the string into an atomic vector.
            count <- sum(a == "C" | a == "G")                                         ####  Counts the number of times the sequence contains a 'G' or 'C'.
            return(count = count)                                                     ####  
        }                                                                             ####  
        seq <- as.character(data[, cols[5]])                                          ####  Assigns the 'Sequence' column to the object seq.
#		GC <- sapply(seq, gc)                                                         ####  Applies the function gc to all sequences in the data.
#		for (i in 1:nobs) {                                                           ####  
#			gcCount[i, 1] <- as.numeric(GC[[i]][1])                                   ####  Populates the matrix gcCount with the values returned by the function gc.
#		}
		GC <- apply(as.matrix(seq),1,gc)
		gcCount <- GC                                                                 ####  
    }                                                                                 ####  
	if (filter != FALSE) {
		if (all(filter.cols == "default")) {
			filter.cols <- c("gIsFeatNonUnifOL", "rIsFeatNonUnifOL","gIsPosAndSignif", "rIsPosAndSignif")
			}			
		filter.values.matrix <- matrix(c("gIsFeatNonUnifOL", "rIsFeatNonUnifOL",
			"gIsPosAndSignif", "rIsPosAndSignif"),nrow=2,byrow=TRUE)
		filter.values <- match(c(filter.cols), filter.values.matrix)%%2
		filter.cols <- match(c(as.vector(filter.cols)),names(data))
		filter.values <- as.vector(filter.values)
	missing <- which(apply(t(t(data[,filter.cols])==filter.values),1,sum)>0)
	}
	probeID <- cbind(as.character(data[, cols[1]]),as.character(data[, cols[2]]))
    data <- cbind(log2(data[, cols[3]]), log2(data[, cols[4]]))                       ####  Binds the 'rMedianSignal' and 'gMedianSignal' columns in a matrix.
    if (useGC != FALSE) {data <- cbind(data,gcCount)}
    if (filter != FALSE) {
    	results <- list(probeID=probeID,data=data,missing=missing)
    	return(results)
    	}
    else {
    	results <- list(probeID=probeID,data=data)
    	return(results)
    	}
                                                                          ####  Attaches gcCount and returns the matrix at function end.
}                                                                                     ####  Note: Returns log base 2 of the data.

