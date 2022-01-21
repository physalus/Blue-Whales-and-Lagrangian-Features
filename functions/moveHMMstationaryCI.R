stationaryCI <- function(m, alpha=0.95)
{
  
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  beta <- m$mle$beta
  
  if(nrow(beta)==1)
    stop("No covariate effect to plot (nrow(beta)==1).")
   # for differentiation in delta method below
  get_stat <- function(beta,covs,nbStates,i) {
    gamma <- moveHMM:::trMatrix_rcpp(nbStates,beta,covs)[,,1]
    solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))[i]
  }
  
  rawCovs <- m$rawCovs
  gridLength <- 200
  # for confidence intervals
  quantSup <- qnorm(1-(1-alpha)/2)
  
  # loop over covariates
  for(cov in 1:ncol(rawCovs)) {
    inf <- min(rawCovs[,cov],na.rm=TRUE)
    sup <- max(rawCovs[,cov],na.rm=TRUE)
    
    # mean values of each covariate
    meanCovs <- colMeans(rawCovs)
    # set all covariates to their mean, except for "cov"
    # (which takes a grid of values from inf to sup)
    tempCovs <- data.frame(rep(meanCovs[1],gridLength))
    if(length(meanCovs)>1)
      for(i in 2:length(meanCovs))
        tempCovs <- cbind(tempCovs,rep(meanCovs[i],gridLength))
    
    tempCovs[,cov] <- seq(inf,sup,length=gridLength)
    colnames(tempCovs) <- colnames(rawCovs)
    
    desMat <- model.matrix(m$conditions$formula,data=tempCovs)
    
    probs <- stationary(m, covs=desMat)

    # covariance matrix of estimates
    Sigma <- ginv(m$mod$hessian)
    
    # indices corresponding to regression coefficients in m$mod$estimate
    i1 <- length(m$mle$stepPar) + length(m$mle$anglePar) - (!m$conditions$estAngleMean)*nbStates + 1
    i2 <- i1 + length(m$mle$beta) - 1
    gamInd <- i1:i2
    
    lci <- matrix(NA,gridLength,nbStates)
    uci <- matrix(NA,gridLength,nbStates)

    for(state in 1:nbStates) {
      dN <- t(apply(desMat, 1, function(x)
        numDeriv::grad(get_stat,beta,covs=matrix(x,nrow=1),nbStates=nbStates,i=state)))
      
      se <- t(apply(dN, 1, function(x)
        suppressWarnings(sqrt(x%*%Sigma[gamInd,gamInd]%*%x))))
      
      # transform estimates and standard errors to R, to derive CI on working scale,
      # then back-transform to [0,1]
      lci[,state] <- plogis(qlogis(probs[,state]) - quantSup*se/(probs[,state]-probs[,state]^2))
      uci[,state] <- plogis(qlogis(probs[,state]) + quantSup*se/(probs[,state]-probs[,state]^2))
    }
  }
  ciList <- list("xValues" = tempCovs,"lowerCI" = lci, "upperCI" = uci)
  return(ciList)
}
