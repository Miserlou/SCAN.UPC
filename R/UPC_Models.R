##########################################################################
# Copyright (c) 2012 W. Evan Johnson, Boston University School of Medicine
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##########################################################################

UPC_Generic_ExpressionSet = function(expressionSet, sequenceFeatureName=NA, modelType="nn", convThreshold=0.001, higherValuesIndicateHigherExpression=TRUE, verbose=TRUE)
{
  data = exprs(expressionSet)
  meta = pData(featureData(expressionSet))

  gcContent = NULL
  lengths = NULL

  if (!is.na(sequenceFeatureName) & (sequenceFeatureName %in% colnames(meta)))
  {
    message(paste("Extracting sequence information from meta column", sequenceFeatureName))
    sequences = as.character(meta[,sequenceFeatureName])

    if (any(sequences == ""))
      message(paste(sum(sequences == ""), " features will be excluded because they contain an empty ", sequenceFeatureName, " value.", sep=""))

    keep = which(sequences != "")
    data = data[keep,]
    meta = meta[keep,]
    sequences = sequences[keep]

    sequencesSet = DNAStringSet(sequences)
    gcCount = apply(letterFrequency(sequencesSet, letters=c("C", "G")), 1, sum)
    lengths = width(sequencesSet)
    gcContent = gcCount / lengths

    if (length(unique(lengths)) == 1)
      lengths = NULL # No point in correcting for length if all are the same
  }

  upcData = apply(data, 2, function(x) {
    UPC_Generic(x, lengths=lengths, gcContent=gcContent, modelType=modelType, convThreshold=convThreshold, higherValuesIndicateHigherExpression=higherValuesIndicateHigherExpression, verbose=verbose)
  })

  rownames(upcData) = rownames(data)
  exprs(expressionSet) = upcData

  return(expressionSet)
}

UPC_Generic = function(expressionValues, lengths=NULL, gcContent=NULL, modelType="nn", convThreshold=0.001, higherValuesIndicateHigherExpression=TRUE, verbose=TRUE)
{
  if (!(modelType%in%c("nn", "nn_bayes", "ln", "nb")))
    stop(paste("The specified modelType value (", modelType, ") is invalid. Must be nn, nn_bayes, ln, or nb.", sep=""))

  if (is.null(lengths))
    if (verbose)
      message("No annotation information was present for length, so no correction will be made for this.")

  if (is.null(gcContent))
    if (verbose)
      message("No annotation information was present for GC content, so no correction will be made for this.")

  if (!is.null(lengths))
  {
    if (any(is.na(lengths)))
    {
      warning("Could not adjust for length because at least one NA value was present.")
    } else {
      if (length(unique(lengths))==1)
      {
        warning("Could not adjust for length because all values were the same.")
      }
    }

    if (length(lengths) != length(expressionValues))
      stop(paste("The size of expressionValues (", length(expressionValues), ") is not identical to the size of lengths (", length(lengths), ")", sep=""))

    if (all(lengths > 0))
      lengths = log2(lengths)
  }

  if (!is.null(gcContent))
  {
    if (any(is.na(gcContent)))
    {
      warning("Could not adjust for GC content because at least one NA value was present.")
    } else {
      if (length(unique(gcContent))==1)
      {
        warning("Could not adjust for GC content because all values were the same.")
      }
    }

    if (length(gcContent) != length(expressionValues))
      stop(paste("The size of expressionValues (", length(expressionValues), ") is not identical to the size of gcContent (", length(gcContent), ")", sep=""))
  }

  # Check to see if the data have not already been log transformed. If they haven't, then log transform the data.
  if ((max(expressionValues) - min(expressionValues)) > 50)
    expressionValues = log2(expressionValues - min(expressionValues) + 2)

  if (modelType == "nn")
  {
    message("Calculating UPC values using the normal-normal model.")
    upcs = UPC_nn(expressionValues, l=lengths, gc=gcContent, conv=convThreshold)
  }

  if (modelType == "nn_bayes")
  {
    message("Calculating UPC values using the normal-normal Bayesian model.")
    upcs = UPC_nn_bayes(expressionValues, l=lengths, gc=gcContent, conv=convThreshold)
  }

  if (modelType == "ln")
  {
    message("Calculating UPC values using the log-normal model.")
    upcs = UPC_ln(expressionValues, l=lengths, gc=gcContent, conv=convThreshold)
  }

  if (modelType == "nb")
  {
    message("Calculating UPC values using the negative-binomial model.")

    tryNegBinom = function()
    {
      return(UPC_nb(expressionValues, l=lengths, gc=gcContent, conv=convThreshold))
    }
    retryNegBinom = function(e)
    {
      message(e)
      message("\nRetrying...")
      return(UPC_nb(expressionValues, l=lengths, gc=gcContent, conv=convThreshold))
    }
  
    upcs = tryCatch(tryNegBinom(), error=retryNegBinom)
  }

  if (!higherValuesIndicateHigherExpression)
    upcs = 1 - upcs

  return(upcs)
}

UPC_nn = function(y, l=NULL, gc=NULL, conv=0.001, q=1) {
  groups = rep(1, length(y))
  probs = numeric(length(y))

  for (j in unique(groups)){
    use = (groups==j)

    if (is.null(l) && is.null(gc)) {
      X = matrix(1, nrow=length(y[use]), ncol=1)
    }
    if (is.null(l) && !is.null(gc)) {
      X = 1.*cbind(matrix(1, nrow=length(y), ncol=1), gc)[use,]
    }
    if (!is.null(l) && is.null(gc)) {
      X = 1.*cbind(matrix(1, nrow=length(y), ncol=1), l, l^2)[use,]
    }
    if (!is.null(l) && !is.null(gc)) {
      X = 1.*cbind(matrix(1, nrow=length(y), ncol=1), l, l^2, gc)[use,]
    }

    y1=y[use]

    pi_h=.5

    x_1=cbind(X[y1>median(y1),])
    beta_1=solve(t(x_1)%*%x_1)%*%t(x_1)%*%y1[y1>median(y1)]
    sigma2_1=var(y1[y1>median(y1)])

    x_2=cbind(X[y1<=median(y1),])
    beta_2=solve(t(x_2)%*%x_2)%*%t(x_2)%*%y1[y1<=median(y1)]
    sigma2_2=var(y1[y1<=median(y1)])

    theta=c(pi_h,beta_1,beta_2)
    thetaold=theta+1000
    i<- 0

    gamma=1.*(y1>median(y1))
    x=cbind(1,gamma,X[,-1])
    beta=solve(t(x)%*%(x))%*%t(x)%*%y1

    while(max(abs((theta-thetaold)/thetaold))>conv)
    {
      thetaold=theta 
     
      #E-step
      gamma <- (pi_h*dnorm(y1,cbind(1,1,X[,-1])%*%beta,sqrt(sigma2_1)))/(pi_h*dnorm(y1,cbind(1,1,X[,-1])%*%beta,sqrt(sigma2_1))+(1-pi_h)*dnorm(y1,cbind(1,0,X[,-1])%*%beta,sqrt(sigma2_2)))

      # Occasionally a NaN gamma value will be returned. We set these to zero by default.
      gamma[is.nan(gamma)] = 0

      x=cbind(1,gamma,X[,-1])
     
      #M-step
      beta=solve(t(x)%*%(x))%*%t(x)%*%y1
      sigma2_1  <- sigma2_2 <- sum((y1-x%*%beta)^2)/length(y1)
     
      pi_h <- mean(gamma)

      theta=c(pi_h,beta_1,beta_2)
    
      i <- i+1
      message(paste("Iteration", i))
      message(max(abs((theta-thetaold)/thetaold)))

      # Abort if most values are to one extreme
      if (sum(x[,2] > 0.99) > nrow(x) * 0.99)
      {
        warning("Most values were close to 1.0, so stopping. Perhaps try with a higher convergence threshold.")
        break
      }
      if (sum(x[,2] < 0.01) > nrow(x) * 0.99)
      {
        warning("Most values were close to 0.0, so stopping. Perhaps try with a higher convergence threshold.")
        break
      }

      if (i == 100)
      {
        message("Stopped at 100 iterations")
        break
      }
    }

    probs[use]=gamma
  }

  y.sim=rnorm(length(y),x%*%beta,sqrt(sigma2_1))
  #par(mfrow=c(2,1))
  #hist(y,breaks=250, xlim=c(0,30))
  #hist(y.sim,breaks=250,col=2, xlim=c(0,30))

  return(probs)
}

UPC_nn_bayes = function(y, l=NULL, gc=NULL, conv=0.001, q=1, empirical=T, lambda=1, exprProp=0.4, k=1) {
  groups = rep(1, length(y))
  probs = numeric(length(y))

  for (j in unique(groups)){
    use = (groups==j)

    if (is.null(l) && is.null(gc)) {
      X = matrix(1, nrow=length(y[use]), ncol=1)
    }
    if (is.null(l) && !is.null(gc)) {
      X = 1.*cbind(matrix(1, nrow=length(y), ncol=1), gc)[use,]
    }
    if (!is.null(l) && is.null(gc)) {
      X = 1.*cbind(matrix(1, nrow=length(y), ncol=1), l, l^2)[use,]
    }
    if (!is.null(l) && !is.null(gc)) {
      X = 1.*cbind(matrix(1, nrow=length(y), ncol=1), l, l^2, gc)[use,]
    }

    y1=y[use]
    
    n_total <- length(y1)  
    beta_a <- n_total*0.1*exprProp; beta_b <- n_total*0.1*(1-exprProp)  
    pi_h= beta_a/(beta_a+beta_b) 
    
    x_1=cbind(X[y1>median(y1),])
    beta_1=solve(t(x_1)%*%x_1)%*%t(x_1)%*%y1[y1>median(y1)]
    sigma2_1=var(y1[y1>median(y1)])

    x_2=cbind(X[y1<=median(y1),])
    beta_2=solve(t(x_2)%*%x_2)%*%t(x_2)%*%y1[y1<=median(y1)]
    sigma2_2=var(y1[y1<=median(y1)])
    
    # prior covariance of beta1 and beta2
    if (empirical==T){
      cov_beta <- sigma2_1 * solve(t(x_1)%*%x_1) * n_total/k 
      cov_beta_inv <- solve(cov_beta)      
    } else {
      cov_beta_inv <- diag(rep(lambda,NCOL(X)))
    }
    
    theta=c(pi_h,beta_1,beta_2)
    thetaold=theta+1000
    i<- 0

    gamma=1.*(y1>median(y1))
    
    while(max(abs(theta/thetaold-1))>conv)
    {
      thetaold=theta 
      ptm <- proc.time()
      #E-step
      b_div_a <- (1-pi_h)/pi_h * sqrt(sigma2_2/sigma2_1) * exp(-1/2*((y1-as.numeric(X%*%beta_1))^2/sigma2_1 - (y1-as.numeric(X%*%beta_2))^2/sigma2_2))
      gamma <- as.numeric(1/(1+b_div_a))
      
      # Occasionally a NaN gamma value will be returned. We set these to zero by default.
      gamma[is.nan(gamma)] = 0

      x=cbind(1,gamma,X[,-1])
     
      #M-step
    
      beta_1 <- solve(t(X*(1-gamma))%*%X+as.numeric(sigma2_1)*cov_beta_inv) %*% t(X*(1-gamma))%*%y1
      beta_2 <- solve(t(X*gamma)%*%X+as.numeric(sigma2_2)*cov_beta_inv) %*% t(X*gamma)%*%y1
      sigma2_1 <- (t((y1-X%*%beta_1) * (1-gamma)) %*% (y1-X%*%beta_1)) / sum((1-gamma))
      sigma2_2 <- (t((y1-X%*%beta_2) * gamma) %*% (y1-X%*%beta_2)) / sum(gamma)
      
      pi_h <- (beta_a + sum(gamma)) / (n_total+beta_a+beta_b)
      
      theta=c(pi_h, beta_1, beta_2)
      
      i <- i+1
      message(paste("Iteration", i))
      message(max(abs((theta-thetaold)/thetaold)))
      
      # Abort if most values are to one extreme
      if (sum(x[,2] > 0.99) > nrow(x) * 0.99)
      {
        warning("Most values were close to 1.0, so stopping. Perhaps try with a higher convergence threshold.")
        break
      }
      if (sum(x[,2] < 0.01) > nrow(x) * 0.99)
      {
        warning("Most values were close to 0.0, so stopping. Perhaps try with a higher convergence threshold.")
        break
      }

      if (i == 100)
      {
        message("Stopped at 100 iterations")
        break
      }
    }

    probs[use]=gamma
  }

  return(probs)
}

UPC_ln = function(genecounts, l=NULL, gc=NULL, conv=0.001)
{
  return(RS_BC(genecounts, l, gc, conv=conv)$probs)
}

RS_BC=function(y, l=NULL, gc=NULL, start=NULL, conv){
  if (is.null(l) && is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1),1.*(y>median(y)))
  }
  if (is.null(l) && !is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1),1.*(y>median(y)),gc,gc^2,gc^3)
  }
  if (!is.null(l) && is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1),1.*(y>median(y)),l,l^2,l^3)
  }
  if (!is.null(l) && !is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1),1.*(y>median(y)),l,l^2,l^3,gc,gc^2,gc^3)
  }

  paramsold=Inf
  if(is.null(start)){
    p=.5
    theta=mlln(y,x)
   
  } else {
    p=start[1]
    theta1=start[2:3]
    theta2=start[4:5]
  }

  params=c(p,theta$bhat,theta$shat)
  it=0

  while (max(abs((paramsold-params)/params)) > conv)
  {
    #E-STEP
    x_1=cbind(1,1,x[,-c(1:2)]);
    x_2=cbind(1,0,x[,-c(1:2)]);
    gam=p*dlnorm(y,x_1%*%theta$bhat,theta$shat)/(p*dlnorm(y,x_1%*%theta$bhat,theta$shat)+(1-p)*dlnorm(y,x_2%*%theta$bhat,theta$shat))
    x=cbind(1,gam,x[,-c(1:2)])
    
    #M-STEP
    p=mean(gam)
    theta=mlln(y,x)

    paramsold=params
    params=c(p,theta$bhat,theta$shat)
    it=it+1

    message(paste("Iteration ", it, ": c = ", max(abs((paramsold-params)/params)), sep=""))

    if (it == 100)
    {
      message("Stopped after 100 iterations")
      break
    }
  }

  message(paste("Total iterations:", it))

  y.sim=rlnorm(length(y),x%*%theta$bhat,theta$shat)
  #par(mfrow=c(2,1))
  #hist(y,breaks=250, xlim=c(0,30))
  #hist(y.sim,breaks=250,col=2, xlim=c(0,30))
  x_1=cbind(1,1,x[,-c(1:2)]);
  x_2=cbind(1,0,x[,-c(1:2)]);
  gam=p*dlnorm(y,x_1%*%theta$bhat,theta$shat)/(p*dlnorm(y,x_1%*%theta$bhat,theta$shat)+(1-p)*dlnorm(y,x_2%*%theta$bhat,theta$shat))

  return(list(p=p,theta=theta,probs=gam))
}

mlln=function(y,x,gam=NULL){
 if(is.null(gam))
   gam=rep(1,length(y))

 bhat=solve(t(x)%*%(gam*x))%*%t(x)%*%(gam*log(y))
 shat=sqrt(sum(gam*(log(y)-x%*%bhat)^2)/sum(gam))

 return(list(bhat=bhat,shat=shat))
}

UPC_nb = function(y, l=NULL, gc=NULL, start=NULL, conv=0.001, sam=10000){
  if (is.null(l) && is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1))
  }
  if (is.null(l) && !is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1),gc,gc^2)
  }
  if (!is.null(l) && is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1),l,l^2)
  }
  if (!is.null(l) && !is.null(gc)) {
    x = cbind(matrix(rep(1,length(y)),ncol=1),l,l^2,gc,gc^2)
  }

  y1=round(y)

  paramsold=Inf
  use=1:length(y1)

  if (sam < length(y1)){use=sample(use,sam)}

  fit=suppressWarnings(glm.nb(y1[use]~-1+x[use,]))
  keep1=y>exp(x%*%fit$coefficients)
  p=mean(keep1)
  fit1=suppressWarnings(glm.nb(y1[use]~-1+x[use,],weights=keep1[use]))
  fit2=suppressWarnings(glm.nb(y1[use]~-1+x[use,],weights=1-keep1[use]))

  params=c(p,fit1$theta,fit1$coefficients,fit2$theta,fit2$coefficients)
  it=0

  while (max(abs((paramsold-params)/params)) > conv)
  {
    #E-STEP
    gam=p*dnbinom(y1,size=fit1$theta,mu=exp(x%*%fit1$coefficients))/(p*dnbinom(y1,size=fit1$theta,mu=exp(x%*%fit1$coefficients))+(1-p)*dnbinom(y1,size=fit2$theta,mu=exp(x%*%fit2$coefficients)))

    #M-STEP
    p=mean(gam[use])
    fit1=suppressWarnings(glm.nb(y1[use]~-1+x[use,],weights=gam[use]))
  	fit2=suppressWarnings(glm.nb(y1[use]~-1+x[use,],weights=1-gam[use]))

	paramsold=params
    params=c(p,fit1$theta,fit1$coefficients,fit2$theta,fit2$coefficients)
  	it=it+1

    message(paste("Iteration ", it, ": c = ", max(abs((paramsold-params)/params)),sep=""))

    if (it == 100)
    {
      message("Stopped after 100 iterations")
      break
    }
  }

  message(paste("Total iterations:", it))

  gam=p*dnbinom(y1,size=fit1$theta,mu=exp(x%*%fit1$coefficients))/(p*dnbinom(y1,size=fit1$theta,mu=exp(x%*%fit1$coefficients))+(1-p)*dnbinom(y1,size=fit2$theta,mu=exp(x%*%fit2$coefficients)))

  return(gam)
}
