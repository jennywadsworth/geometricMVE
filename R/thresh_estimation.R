#' Find a high threshold of R|W which is (approximately) the tau-quantile of this variable
#' 
#' @param r Values of radial variable
#' @param w values of angular variable (on unit simplex)
#' @param tau level at which to calculate threshold
#' @param method character string specifying either "empirical" for using a binning mathod, or "KDE" for the method based on kernel density estimation as outlines in Campbell and Wadsworth (2025)
#'
#' @return list containing elements r0w (estimated threshold for each given w), r, w, tau and method
#' @export

fit.thresh = function(r,w,tau=0.95,method=c("empirical","KDE")){
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
    r0w = radial.quants.L1.KDE.2d(r=r,w=w,tau=tau,bww=0.05)
  } else if(dim(w)[2]>2 & method=="KDE") {
    r0w = radial.quants.L1.KDE(r=r,w=w,tau=tau,bww=0.05)
  } else if(method=="empirical"){
    r0w = emp.thresh(r=r,w=w,tau=tau)
  }
  return(list(r0w=r0w,
              w=w,
              r=r,
              tau=tau,
              method=method))
}

eval.thresh = function(fit.thresh.out,w.eval){
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
    r0w = KDE.quant.eval.2d(wpts=w.eval,
                            r=fit.thresh.out$r,w=fit.thresh.out$w)
  } else if(dim(w.eval)[2]>2 & fit.thresh.out$method=="KDE") {
    r0w = KDE.quant.eval(wpts=w.eval,
                         r=fit.thresh.out$r,w=fit.thresh.out$w)
  } else if(fit.thresh.out$method=="empirical"){
    stop("Not yet implemented for empirical threshold.")
  }
  
  return(r0w)

}

emp.thresh = function(r, w, tau = 0.95, bin.mesh = NULL, overlap = NULL){
  
  bin.mesh = floor(100/(dim(w)[2]))

  overlap = 2/(bin.mesh)  # <-- this needs to be played with

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

# kernel density estimation of the cdf when angles are defined by the L1-norm, then invert it
radial.quants.L1.KDE.2d = function(r,w,tau=0.95,bww=0.05,bwr=0.05){
  # r, w              -> vectors
  # bww, bwr          -> bandwidths affects smoothness / how close you can get to "pointy" r_0(w)

  r0w = sapply(w, function(ww){
    weightsw<-dnorm(w,mean=ww,sd=bww)

    ccdf<-function(rc){
      mean(weightsw*pnorm(rc,mean=r,sd=bwr))/mean(weightsw)
    }
    dummy<-function(rc){ccdf(rc) - tau}  # want the root of this
    ur<-uniroot(dummy,interval = c(0,30))
    return(ur$root)
  })

  return(r0w)
}

## General d-dimensional version of QR
radial.quants.L1.KDE = function(r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){ #ker.pdf=mv.Gaussian.ker.pdf, ker.cdf=Gaussian.ker.cdf){
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

#########################

# A function for evaluating r_{\tau}(w) at w

KDE.quant.eval.2d = function(wpts,r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){
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

KDE.quant.eval = function(wpts,r,w,tau=0.95,bww=0.05,bwr=0.05,up=30){
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

