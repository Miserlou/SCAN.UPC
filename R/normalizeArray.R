##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

normalizeArray <-
function (data = 0, gcCount = NULL, useGC = FALSE, minGC = 3, signal="median",    ####  
    filter= FALSE ,filter.cols='default', filename = "NULL", type = "c", scale = "n")                                ####  
{                                                                                             ####  
    if (filename != "NULL") {                                                                 ####  
    	temp <- readAgilentData(data=data, filename=filename, signal=signal, useGC=useGC, filter=filter, filter.cols=filter.cols)
        data <- temp$data                                                                     ####  reads data if 'filename' is specified.
        if (useGC != FALSE) {
        gcCount <- data[, 3]                                                                  ####  
        }
        probeID <- temp$probeID
        missing <- temp$missing
    	}                                                                                         ####  
    norm <- NULL                                                                              ####  Clears 'norm' object.
    if (type == "c") {                                                                        ####  
        norm <- channelNormalize(data, gcCount, useGC, minGC)                                 ####  Calls 'channelNormalize' if ' type="c" '.
    	}                                                                                         ####  
    if (type == "q") {                                                                        ####  
        norm <- quantileNormalize(data)                                                       ####  Calls 'quantileNormalize' if ' type="q" '.
    	}                                                                                         ####  
    if (type == "i") {                                                                        ####  
        norm <- iglNormalize(data)                                                            ####  Calls 'iglNormalize' if ' type="i" '.
    	}                                                                                         ####  
    if (scale != "n") {                                                                       ####  
        norm <- madNormalize(norm)                                                            ####  Calls 'madNormalize' if ' scale="n" '.
    	}                                                                                         ####  
    if (filename != "NULL") {
    	results <- list(probeID=probeID,norm=norm,missing=missing)
   		}
    else {
    	results <- list(norm=norm)
    	}
	return(results)                                                                           ####  Returns normalized data.
}                                                                                             ####  

