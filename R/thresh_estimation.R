#' Find a high threshold of R|W which is (approximately) the tau-quantile of this variable
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param method character string specifying either "empirical" for using a binning method, or "KDE" for the method based on kernel density estimation as outlines in Campbell and Wadsworth (2025)
#' @param bww bandwidth for W kernel when method="KDE"
#' @param bin.mesh numerical value affecting number of bins for estimation when method="empirical"
#' @param overlap numerical value affecting overlap of bins for estimation when method="empirical"
#' @param up numerical value giving upper limit for root finding of r0w when method="KDE" (default 30)

#' @return list containing elements r0w (estimated threshold for each given w), r, w, tau and method
#' @export

fit.thresh = function(r,w,tau=0.95,method=c("empirical","KDE"),bww=0.05,bin.mesh=NULL,overlap=NULL,up=30){
  if(any(w<0)){stop("values of w must be in the unit simplex")}
  if(is.vector(w)){
    w<-cbind(w,1-w)
  } else{
    sw<-apply(w,1,sum)
    if(all(sw<1)){
      w<-cbind(w,1-sw)
    }
  }
  
  if(dim(w)[2]==2 & method=="KDE"){
    w = w[,1]
    r0w = radial.thresh.KDE.2d(r=r,w=w,tau=tau,bww=bww,up=up)
  } else if(dim(w)[2]>2 & method=="KDE") {
    r0w = radial.thresh.KDE(r=r,w=w,tau=tau,bww=bww,up=up)
  } else if(dim(w)[2]==2 & method=="empirical"){
    # r0w = emp.thresh(r=r,w=w,tau=tau)
    emp.thresh.output = emp.thresh.2d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]==3 & method=="empirical"){
    emp.thresh.output = emp.thresh.3d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]==4 & method=="empirical"){
    emp.thresh.output = emp.thresh.4d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]==5 & method=="empirical"){
    emp.thresh.output = emp.thresh.5d(r=r,w=w,tau=tau,bin.mesh=bin.mesh,overlap=overlap)
  } else if(dim(w)[2]>5 & method=="empirical"){
    stop("Empirical threshold estimation is not yet implemented for d>5. Use the KDE method instead.")
  }
  if(method=="KDE"){
    return(list(r0w=r0w,
                w=w,
                r=r,
                tau=tau,
                method=method,
                bww=bww))
  } else if(method=="empirical"){
    return(list(r0w=emp.thresh.output$r0w,
                w=w,
                r=r,
                tau=tau,
                quant.grid=emp.thresh.output$quant.grid,  # needed for eval.thresh.emp...
                method=method,
                bin.mesh=bin.mesh,
                overlap=overlap))
  }
}

#' Evaluate a high threshold of R|W at a set of new angles W
#' 
#' @param fit.thresh.out the output of the 'fit.thresh' function
#' @param w.eval a vector of length n or a matrix with n rows
#' @param up numerical value giving upper limit for root finding of r0w when using KDE method (default 30)
#' 
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @export
eval.thresh = function(fit.thresh.out,w.eval,up=30){
  
  # if(any(w.eval<0)){stop("values of w must be in the unit simplex")}
  if(is.vector(w.eval)){
    w.eval<-cbind(w.eval,1-w.eval)
  } else{
    sw<-apply(w.eval,1,sum)
    if(all(sw<1)){
      w.eval<-cbind(w.eval,1-sw)
    }
  }
  
  if(dim(w.eval)[2]==2 & fit.thresh.out$method=="KDE"){
    w.eval = w.eval[,1]
    r0w = KDE.thresh.eval.2d(wpts=w.eval,
                             r=fit.thresh.out$r,w=fit.thresh.out$w,bww=fit.thresh.out$bww,up=up)
  } else if(dim(w.eval)[2]>2 & fit.thresh.out$method=="KDE") {
    r0w = KDE.thresh.eval(wpts=w.eval,
                          r=fit.thresh.out$r,w=fit.thresh.out$w,bww=fit.thresh.out$bww,up=up)
  } else if(dim(w.eval)[2]==2 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.2d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]==3 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.3d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]==4 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.4d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]==5 & fit.thresh.out$method=="empirical"){
    r0w = emp.thresh.eval.5d(wpts=w.eval,quant.grid=fit.thresh.out$quant.grid,bin.mesh=fit.thresh.out$bin.mesh,overlap=fit.thresh.out$overlap)
  } else if(dim(w.eval)[2]>5 & fit.thresh.out$method=="empirical"){
    stop("Empirical threshold estimation is not yet implemented for d>5. Use the KDE method instead.")
  }
  
  return(r0w)
  
}

############################################################################

# Empirical threshold estimation

#' Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  stop("General d-dimension empirical threshold estimation not yet implemented. Check back later.")
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with

  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  qu.r <- array(0, dim = rep(bin.mesh, (dim(w)[2]-1)))
  array.idx = expand.grid(replicate(dim(w)[2]-1, c(1:bin.mesh), simplify = FALSE))
  for(rrow in 1:nrow(array.idx)){
    # TODO: SPEED UP THE NEXT 4 LINES OF CODE
    idx = array.idx[rrow,]
    ind = rep(TRUE,nrow(w))
    for(ii in 1:length(idx)){
      ind = ind & (wsl[as.numeric(idx[ii])] < w[,ii]) & (w[,ii] < wsu[as.numeric(idx[ii])])
    }
    r1 <- r[ind]
    
    if (sum(ind) > 2) {
      qu.r[matrix(as.numeric(idx),1)] <- quantile(r1, tau)
    }
    else {
      qu.r[matrix(as.numeric(idx),1)] <- NA
    }
  }
  
  r0w <- rep(0,nrow(w))
  # array.idx = expand.grid(replicate(dim(w)[2]-1, c(1:bin.mesh), simplify = FALSE))
  for(i in 1:length(r0w)){
    indArray <- array(NA, dim = rep(bin.mesh,dim(w)[2]-1))
    for(rrow in 1:nrow(array.idx)){
      # TODO: SPEED UP THE NEXT 4 LINES OF CODE
      idx = array.idx[rrow,]
      ind = TRUE
      for(ii in 1:length(idx)){
        ind = ind & (wsl[as.numeric(idx[ii])] < w[i,ii]) & (w[i,ii] < wsu[as.numeric(idx[ii])])
      }
      indArray[matrix(as.numeric(idx),1)] = ind
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }
  return(r0w)
}

#' d=2 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.2d = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  if(is.null(dim(w))){
    stop("w needs to be a n x 2 matrix")
  }
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w = w[,1]  # make into a vector
  
  qu.r <- NULL
  for (j in 1:bin.mesh) {
    re <- r[wsl[j] < w & w < wsu[j]]
    qu.r[j] <- quantile(re, tau)
  }
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indvec <- rep(NA, bin.mesh)
    for (j in 1:bin.mesh) {
      indvec[j] <- wsl[j] < w[i] & w[i] < wsu[j]
    }
    r0w[i] <- mean(qu.r[indvec], na.rm = T)
  }
  return(list(r0w=r0w,             # the threshold
              quant.grid=qu.r))
}

#' d=2 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.2d = function(wpts,quant.grid,bin.mesh,overlap){
  if(is.null(dim(wpts))){
    stop("wpts needs to be a n x 2 matrix")
  }
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w = wpts[,1]  # make into a vector
  
  n <- length(w)
  r0w <- numeric(n)
  for (i in 1:n) {
    indvec <- rep(NA, bin.mesh)
    for (j in 1:bin.mesh) {
      indvec[j] <- wsl[j] < w[i] & w[i] < wsu[j]
    }
    r0w[i] <- mean(quant.grid[indvec], na.rm = T)
  }
  
  return(r0w)
}

#' d=3 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.3d = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1} 

  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = w[,1]  # make into a vector
  w2 = w[,2]  # make into a vector
  
  qu.r <- matrix(0, bin.mesh, bin.mesh)
  for (j in 1:bin.mesh) {
    for (k in 1:bin.mesh) {
      ind <- 
        wsl[j] < w1 & w1 < wsu[j] & 
        wsl[k] < w2 & w2 < wsu[k]
      r1 <- r[ind]
      if (sum(ind) > 10) {
        qu.r[j, k] <- quantile(r1, tau)
      }
      else {
        # empt=empt+1
        # qu.r[j, k] <- NA
        qu.r[j, k] <- quantile(r1, tau)
      }
    }
  }
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indMat <- matrix(NA, bin.mesh, bin.mesh)
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        indMat[j, k] <- 
          wsl[j] < w1[i] & w1[i] < wsu[j] &
          wsl[k] < w2[i] & w2[i] < wsu[k]
      }
    }
    r0w[i] <- mean(qu.r[indMat], na.rm = T)
  }
  
  return(list(r0w=r0w,
              quant.grid=qu.r))
}

#' d=3 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.3d = function(wpts,quant.grid,bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = wpts[,1]  # make into a vector
  w2 = wpts[,2]  # make into a vector
  
  n <- length(w1)
  r0w <- numeric(n)
  for (i in 1:n) {
    indMat <- matrix(NA, bin.mesh, bin.mesh)
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        indMat[j, k] <- 
          wsl[j] < w1[i] & w1[i] < wsu[j] &
          wsl[k] < w2[i] & w2[i] < wsu[k]
      }
    }
    r0w[i] <- mean(quant.grid[indMat], na.rm = T)
  }
  
  return(r0w)
}

#' d=4 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.4d = function(r, w, tau = 0.95, bin.mesh, overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1} 
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = w[, 1]
  w2 = w[, 2]
  w3 = w[, 3]
  
  qu.r <- array(0, dim = rep(bin.mesh,3))
  for (j in 1:bin.mesh) {
    for (k in 1:bin.mesh) {
      for (l in 1:bin.mesh) {
        ind <- 
          wsl[j] < w1 & w1 < wsu[j] & 
          wsl[k] < w2 & w2 < wsu[k] & 
          wsl[l] < w3 & w3 < wsu[l]
        r1 <- r[ind]
        if (sum(ind) > 10) {
          qu.r[j, k, l] <- quantile(r1, tau)
        }
        else {
          qu.r[j, k, l] <- NA
        }
      }
    }
  }
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,3))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          indArray[j, k, l] <- 
            wsl[j] < w1[i] & w1[i] < wsu[j] & 
            wsl[k] < w2[i] & w2[i] < wsu[k] &
            wsl[l] < w3[i] & w3[i] < wsu[l]
        }
      }
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }
  
  return(list(r0w=r0w,
              quant.grid=qu.r))
}

#' d=4 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.4d = function(wpts,quant.grid, bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = wpts[, 1]
  w2 = wpts[, 2]
  w3 = wpts[, 3]
  
  n <- length(w1)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,3))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          indArray[j, k, l] <- 
            wsl[j] < w1[i] & w1[i] < wsu[j] & 
            wsl[k] < w2[i] & w2[i] < wsu[k] &
            wsl[l] < w3[i] & w3[i] < wsu[l]
        }
      }
    }
    r0w[i] <- mean(quant.grid[indArray], na.rm = T)
  }
  
  return(r0w)
}

#' d=5 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using an empirical approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.5d = function(r, w, tau = 0.95,bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1} 
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = w[, 1]
  w2 = w[, 2]
  w3 = w[, 3]
  w4 = w[, 4]
  
  qu.r <- array(0, dim = rep(bin.mesh,4))
  for (j in 1:bin.mesh) {
    for (k in 1:bin.mesh) {
      for (l in 1:bin.mesh) {
        for(m in 1:bin.mesh) {
          ind <- 
            (wsl[j] < w1) & (w1 < wsu[j]) &
            (wsl[k] < w2) &( w2 < wsu[k]) &
            (wsl[l] < w3) &( w3 < wsu[l]) &
            (wsl[m] < w4) & (w4 < wsu[m])
          # print(c(,sum(ind)))
          # print(ind)
          r1 <- r[ind]
          if (sum(ind) > 2) {
            qu.r[j, k, l, m] <- quantile(r1, tau)
          }
          else {
            qu.r[j, k, l, m] <- NA
          }
        }
      }
    }
  }
  
  n <- length(r)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,4))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          for (m in 1:bin.mesh) {
            indArray[j, k, l, m] <-
              wsl[j] < w1[i] & w1[i] < wsu[j] &
              wsl[k] < w2[i] & w2[i] < wsu[k] &
              wsl[l] < w3[i] & w3[i] < wsu[l] &
              wsl[m] < w4[i] & w4[i] < wsu[m]
          }
        }
      }
    }
    r0w[i] <- mean(qu.r[indArray], na.rm = T)
  }
  
  return(list(r0w=r0w,
              quant.grid=qu.r))
}

#' d=5 setting: Evaluate a high threshold of R|W at a set of new angles W using the empirical method
#' 
#' @param wpts a vector of length n or a matrix with n rows
#' @param quant.grid the array of quantiles on the overlapping grid used for fitting the threshold
#' @param bin.mesh numerical value affecting number of bins for estimation
#' @param overlap numerical value affecting overlap of bins for estimation

#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
emp.thresh.eval.5d = function(wpts,quant.grid,bin.mesh,overlap){
  
  if(is.null(bin.mesh)){bin.mesh = floor(100/(dim(w)[2]))}
  if(is.null(overlap)){overlap = 0.1}  # <-- this needs to be played with
  
  wsl <- seq(0, 1 - overlap, len = bin.mesh)
  wsu <- seq(overlap, 1, len = bin.mesh)
  
  w1 = wpts[, 1]
  w2 = wpts[, 2]
  w3 = wpts[, 3]
  w4 = wpts[, 4]
  
  n <- length(w1)
  r0w <- numeric(n)
  for (i in 1:n) {
    indArray <- array(NA, dim = rep(bin.mesh,4))
    for (j in 1:bin.mesh) {
      for (k in 1:bin.mesh) {
        for (l in 1:bin.mesh) {
          for (m in 1:bin.mesh) {
            indArray[j, k, l, m] <-
              wsl[j] < w1[i] & w1[i] < wsu[j] &
              wsl[k] < w2[i] & w2[i] < wsu[k] &
              wsl[l] < w3[i] & w3[i] < wsu[l] &
              wsl[m] < w4[i] & w4[i] < wsu[m]
          }
        }
      }
    }
    r0w[i] <- mean(quant.grid[indArray], na.rm = T)
  }
  
  return(r0w)
}

############################################################################

# KDE threshold estimation

#' d=2 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using a KDE approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bww the angular bandwidth (default=0.05)
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
radial.thresh.KDE.2d = function(r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){

  r0w = sapply(w, function(ww){
    weightsw<-dnorm(w,mean=ww,sd=bww)

    ccdf<-function(rc){
      mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,up))
    return(ur$root)
  })

  return(r0w)
}

#' d=2 setting: Evaluate a high threshold of R|W at a set of new angles W using a KDE method
#' 
#' @param wpts a vector of length n or a matrix with n rows and 2 columns
#' @param r Values of radial variable used to create the KDE threshold
#' @param w values of angular variable (on unit simplex) used to create the KDE threshold
#' @param tau level at which to calculate threshold
#' @param bww the angular bandwidth (default=0.05)
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
KDE.thresh.eval.2d = function(wpts,r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){
  # r, w              -> vectors
  # bww, bwr          -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)
  # n.mesh            -> mesh for wpts
  # ker.pdf, ker.cdf  -> kernel pdf and cdf functions
  
  
  r.tau.wpts = sapply(wpts,function(wpts.i){
    weightsw<- dnorm(w,mean=wpts.i,sd=bww)
    pos.weights = weightsw>0
    weightsw = weightsw[pos.weights]
    ccdf<-function(rc){
      ker.vals = pnorm(rc,mean=r,sd=bwr)[pos.weights]
      num = weightsw*ker.vals
      denom = weightsw
      sum(num,na.rm=T)/sum(denom,na.rm=T)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    
    is_error <- FALSE
    tryCatch({
      ur<-uniroot(dummy,interval = c(0,up))$root
    },error=function(e){
      is_error <<- TRUE
    })
    if(is_error) {
      ur=NA
    }
    return(ur)
  })
  
  return(r.tau.wpts)
}

#' d>2 setting: Find a high threshold of R|W which is (approximately) the tau-quantile of this variable using a KDE approach
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param bww the angular bandwidth (default=0.05)
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
radial.thresh.KDE = function(r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){ #ker.pdf=mv.Gaussian.ker.pdf, ker.cdf=Gaussian.ker.cdf){
  require(mvtnorm)
  # r, w      -> vector and matrix, angles defined by the L1 norm
  # bww, bwr  -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)

  num.cols = dim(w)[2]

  r0w = apply(w, 1, function(ww){

    weightsw = dmvnorm(w[,-num.cols],mean=ww[-num.cols],sigma=(bww^2)*diag(num.cols-1))
    
    ccdf<-function(rc){
      mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,up))
    return(ur$root)
  })
  return(r0w)
}

#' d>2 setting: Evaluate a high threshold of R|W at a set of new angles W using a KDE method
#' 
#' @param wpts a vector of length n or a matrix with n rows and 2 columns
#' @param r Values of radial variable used to create the KDE threshold
#' @param w values of angular variable (on unit simplex) used to create the KDE threshold
#' @param tau level at which to calculate threshold
#' @param bww the angular bandwidth (default=0.05)
#' @param bwr the radial bandwidth (default=0.05)
#' @param up upper limit for root finding of r0w (default 30)
#'
#' @return vector containing elements r0w (estimated threshold for each given w.eval)
#' @noRd
KDE.thresh.eval = function(wpts,r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){
  require(mvtnorm)

  n.dims = dim(w)[2]

  r.tau.wpts = apply(wpts,1,function(wpts.i){
    if(sum(wpts.i)>1){
      return(NA)
    } else{
      weightsw = dmvnorm(w[,-n.dims],mean=wpts.i[-n.dims],sigma=(bww^2)*diag(n.dims-1))
      
      ccdf<-function(rc){
        mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
      }
      dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
      is_error <- FALSE
      tryCatch({
        ur<-uniroot(dummy,interval = c(0,up))$root
      },error=function(e){
        is_error <<- TRUE
      })
      if(is_error) {
        ur=NA
      }
      return(ur)
    }
  })
  return(r.tau.wpts)
}

