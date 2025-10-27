
#' Simulate new data (using resampled W|R'>k) from the model, conditional on R'>k (R>k*r_0(W))

#' @param fit list output from fit.geometric.par or fit.geometric.pwlin
#' @param nsim number of new datapoints to simulate
#' @param k Value of k (defaults to 1). If k>1 then importance weights are calculated to simulate from the distribution of X|R>kr_0(W) (equivalently X|R'>k)
#' @param nwg number of points to use in expanding wgrid for numerical rescaling of additive mixture in >2 dimensions
#'
#' @return nsim by 2 matrix of values (in exponential margins) simulated from the model (given R>r_0(W), or R>kr_0(W) for k>1)
#'
#' @export
#'

sim.geometric<-function(fit,nsim,k=1,nwg=70)
{
  if(k<1){stop("k should be greater than or equal to 1")}
  # if(fit$class=="parametric")
  # {
    w<-fit$data$w
    r<-fit$data$r
    r0w<-fit$data$r0w

  if(any(w<0)){stop("Invalid values of W")}
  if(is.vector(w)){
    w<-cbind(w,1-w)
  } else{
    sw<-apply(w,1,sum)
    if(all(sw<1)){
      w<-cbind(w,1-sw)
    }
  }
  d<-dim(w)[2]
  if(k>1){
    if(fit$setup$add.gauge)
    {
      if(!fit$fixshape){par1ind<-fit$setup$par1ind-1;par2ind<-fit$setup$par2ind-1
      } else{par1ind<-fit$setup$par1ind;par2ind<-fit$setup$par2ind}
      iw<-iweights(k=k, r0w=r0w, w=w, par=fit$mle,shape=fit$shape,
                   add.gauge=fit$setup$add.gauge,
                   gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,
                   par1ind=par1ind,
                   par2ind=par2ind)
    } else{
      iw<-iweights(k=k, r0w=r0w, w=w, gfun=fit$setup$gfun, par=fit$mle,
                   shape=fit$shape,add.gauge=fit$setup$add.gauge)
    }
    ###
    star.ind<-sample(1:length(w[,1]),size=nsim,replace=T,prob = iw)
    wstar<-w[star.ind,]
    r0w_star<-c(k*r0w)[star.ind]
  }else{
    star.ind<-sample(1:length(w[,1]),size=nsim,replace=T)
    wstar<-w[star.ind,]
    r0w_star<-r0w[star.ind]
  }

  if(!fit$setup$add.gauge | is.null(fit$setup$add.gauge)){
    rate0<-apply(w,1,fit$setup$gfun,par=fit$mle)
    rate<-rate0[star.ind]
    rstar<-qgamma(1-runif(nsim)*pgamma(r0w_star,
                                       shape=fit$shape,rate=rate,
                                       lower.tail=F),
                  shape=fit$shape,rate=rate)
  }
  else{
    wgrid<-create.wgrid(nwg,d=d)
    if(!fit$fixshape){par1ind<-fit$setup$par1ind-1;par2ind<-fit$setup$par2ind-1
    } else{par1ind<-fit$setup$par1ind;par2ind<-fit$setup$par2ind}
    ms<-additivegauge.scaling(gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,
                              weight=fit$mle[length(fit$mle)],par1=fit$mle[par1ind],
                              par2=fit$mle[par2ind],wgrid=wgrid,d=d)
    rate0<-apply(w,1,additivegauge.rescale,
                 par1=fit$mle[par1ind],par2=fit$mle[par2ind],ms=ms,
                 gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,weight=fit$mle[length(fit$mle)])
    rate<-rate0[star.ind]
    rstar<-qgamma(1-runif(nsim)*pgamma(r0w_star,shape=fit$shape,
                                       rate=rate,
                                       lower.tail=F),
                  shape=fit$shape,rate=rate)
  }
  xstar<-rstar*wstar
  return(xstar)
  # }
}

#=========================================================

#==============================================================================
#' Calculate importance weights for simulation of W|R'>k, k>1, (equivalently W|R>k r_0(W), where R'=R/r_0(W)).

#' @param k value of k>1: the importance weights are to simulate from the distribution of W|R'>k
#' @param w Values of W that correspond to the threshold exceedances of R
#' @param r0w The threshold function r_0(w) for each value of W supplied
#' @param gfun fitted gauge function
#' @param par fitted parameter values (first value should be the gamma shape parameter, subsequent values are the gauge function parameters)
#' @param add.gauge logical. If TRUE, then we additively mix (with numerical rescaling) the gauge functions specified by gauge1 and gauge2
#' @param gauge1 first gauge function for additive mixing. Use if add.gauge=TRUE
#' @param gauge2 second gauge function for additive mixing. Use if add.gauge=TRUE
#' @param par1ind vector of indices specifying the position of the parameters for gauge1 in par
#' @param par2ind vector of indices specifying the position of the parameters for gauge2 in par
#' @param nwg number of points to use in expanding wgrid for numerical rescaling of additive mixture in >2 dimensions
#'
#' @return vector of (unnormalized) importance weights to be used in sim.2d
#'


iweights<-function(k, r0w, w, gfun, par, add.gauge=FALSE,gauge1,gauge2,
                   par1ind,par2ind,shape,nwg)
{
  if(k<1){stop("k should exceed 1")}
  if(any(w<0)){stop("Invalid values of W")}
    if(is.vector(w)){
      w<-cbind(w,1-w)
    } else{
      sw<-apply(w,1,sum)
      if(all(sw<1)){
        w<-cbind(w,1-sw)
      }
    }
  d<-dim(w)[2]
  if(add.gauge){
    if(d>2){wgrid<-create.wgrid(nwg,d=d)
    } else{wgrid<-NULL}

    # Check par1ind, par2ind are correct coming from fit$setup

    ms<-additivegauge.scaling(gauge1=gauge1,gauge2=gauge2,
                              weight=par[length(par)],par1=par[par1ind],
                              par2=par[par2ind],d=d,wgrid=wgrid)
    rate<-apply(w,1,additivegauge.rescale,par1=par[par1ind],
                par2=par[par2ind],ms=ms,gauge1=gauge1,gauge2=gauge2,
                weight=par[length(par)])
  }
  else{
    rate<-apply(w,1,gfun,par=par)
  }
  return(pgamma(k*r0w,shape=shape,rate=rate,lower.tail = F)/pgamma(r0w,shape=par[1],rate=rate,lower.tail = F))
}

#==============================================================================
#' Calculate P(R'>k|R'>1) (equivalently P(R>kr_0(W)|R>r_0(w)))

#' @param k value of k>1
#' @param fit list output from fit.geometric.par or fit.geometric.pwlin
#' @param nwg value giving resolution for numerical rescaling of additive gauges (where needed). Larger numbers give greater accuracy, but take longer.
#'
#' @return Estimate of P(R'>k|R'>1)
#'
#' @export

Rexc.prob.k<-function(k, fit, nwg=70)
{
  if(!fit$fixshape){par1ind<-fit$setup$par1ind-1;par2ind<-fit$setup$par2ind-1
  } else{par1ind<-fit$setup$par1ind;par2ind<-fit$setup$par2ind}
  iw<-iweights(k=k, r0w=fit$data$r0w, w=fit$data$w, gfun=fit$setup$gfun,
               par=fit$mle, add.gauge=fit$setup$add.gauge,
               gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,
               par1ind=par1ind,par2ind=par2ind,nwg=nwg)
  return(mean(iw))
}

