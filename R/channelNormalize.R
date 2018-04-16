##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

channelNormalize <-
function (data = 0, gcCount = NULL, useGC = FALSE, minGC = 0)   ####  
{                                                                                   ####  
    norm <- matrix(NA, nrow(data), 2)                                               ####  Declares 'norm' matrix.
    if (useGC == FALSE) {                                                           ####  
        mu <- apply(data[, 1:2], 2, mean)                                           ####  Finds the mean of the red and green signals.
        sigma <- var(data[, 1:2])                                                   ####  Finds the variance covariance matrix of the red and green signals.
        c <- eigen(sigma)$vectors                                                   ####  Finds Eigen vectors of the variance covariance marix.
        d <- diag(eigen(sigma)$values)                                              ####  Finds Eigen values of the variance covariance marix.
        sigma.sqrt <- c %*% sqrt(d) %*% solve(c)                                    ####  Finds the standard deviation matrix
        normalize <- solve(sigma.sqrt) %*% t(cbind(data[, 1] -                      ####  Normalizes the data for the entire data set.
            mu[1], data[, 2] - mu[2]))                                              ####  
        norm[, 1:2] <- t(normalize)                                                 ####  
    }                                                                               ####  
    else {                                                                          ####  
        gcCount <- as.numeric(gcCount)                                              ####  Changes the values of cgCount to numeric variables.
        gcGroups <- NULL                                                            ####  Resets 'gcGroups' object
        for (i in sort(unique(gcCount))) {                                          ####  Runs a loop for every value of gcCount
            if (sum(gcCount == i) >= minGC) {                                       ####  If there are more than 'minGC' sequences with 'i' G's and C's 
                gcGroups <- c(gcGroups, i)                                          ####     then 'i' gets added to gcGroups.
            }                                                                       ####  
        }                                                                           ####  
        for (i in sort(unique(gcCount))) {                                          ####  
            if (sum(gcCount == i) < minGC) {                                        ####  Combines the groups with less than 'minGC' with the 
                gcCount[gcCount == i] <- gcGroups[order(abs(gcGroups -              ####     nearest group with more than 'minGC'
                  i))[1]]                                                           ####  
            }                                                                       ####  
        }                                                                           ####  
        for (i in gcGroups) {                                                       ####  
            y <- data[gcCount == i, 1:2]                                            ####  Normalizes the data for each value of gcGroups separately
            mu <- apply(y, 2, mean)                                                 ####      the same way as above.
            sigma <- var(y)                                                         ####  
            c <- eigen(sigma)$vectors                                               ####  
            d <- diag(eigen(sigma)$values)                                          ####  
            sigma.sqrt <- c %*% sqrt(d) %*% solve(c)                                ####  
            normalize <- solve(sigma.sqrt) %*% t(cbind(y[, 1] -                     ####  
                mu[1], y[, 2] - mu[2]))                                             ####  
            norm[gcCount == i, ] <- t(normalize)                                    ####  
        }                                                                           ####  
    }                                                                               ####  
    return(norm)                                                                    ####  Returns normalized data
}                                                                                   ####  

