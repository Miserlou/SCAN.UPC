##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

iglNormalize <-
function (data = 0)                        ####  
{                                                          ####  
    norm <- matrix(NA, length(data[, 1]), 2)               ####  
    r <- data[, 1]                                         ####  'r' = Red Signal
    g <- data[, 2]                                         ####  'g' = Green Signal
    a <- (r + g)/2                                         ####  'a' = mean of corresponding Red and Green signals
    norm[, 1] <- r - loess(r ~ a, span = 0.2)$fitted       ####  norm = deviations from loess fitted values
    norm[, 2] <- g - loess(g ~ a, span = 0.2)$fitted       ####  
    return(norm)                                           ####  
}                                                          ####  

