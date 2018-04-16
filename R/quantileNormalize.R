##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

quantileNormalize <-
function (data = 0)                              ####  
{                                                                     ####  
    xSort <- cbind(sort(data[, 1]), sort(data[, 2]))                  ####  Sorts the values of the red and green signals separately
    xpSort <- apply(xSort, 1, mean)                                   ####  Returns a matrix with the mean of the rows of xSort
    norm <- matrix(NA, length(data[, 1]), 2)                          ####  
    norm[order(data[, 1]), 1] <- xpSort                               ####  Replaces the original value with the mean of it and its 
    norm[order(data[, 2]), 2] <- xpSort                               ####     ordered counterpart.
    return(norm)                                                      ####  
}                                                                     ####  

