
#' Fit truncated gamma model for R|W for bivariate data

#' @param r Values of R|W=w that exceed the threshold r0w
#' @param w Values of W that correspond to the threshold exceedances of R
#' @param r0w The threshold function r_0(w) for each value of W supplied
#' @param model Character string specifying model to use for the gauge function. Options are "log" for logistic-type gauge, "invlog" for inverted logistic, "gauss" for Gaussian, "square" for square. If two options are given, then the fitted gauge will be an additive mixture of these two.
#' @param customgauge Optional argument to specify a custom gauge function. This should have
#' @param fixshape Logical. If TRUE, then the shape parameter of the truncated gamma distribution is fixed to the dimension d, otherwise this is estimated.
#' @param init.val vector of initial gauge function parameter values for optimization (needed when customgauge used, otherwise optional)
#' @param lower.limit vector of lower limits on parameters (when customgauge used; should be same length as init.val)
#' @param upper.limit vector of upper limits on parameters (when customgauge used; should be same length as init.val)
#' @param nwg number of points used in expanding mesh to cover simplex for numerically rescaled additive gauge
#' @param ... Additional arguments to pass to the gauge function (if necessary)
#'
#' @return list containing elements "mle", "nllh", "convergence", and additionally "cov" if argument hessian=T is supplied to optim
#'
#' @export


fit.geometric.par<-function(r, w, r0w, model, customgauge=NULL, fixshape=FALSE,
                            init.val=NULL, lower.limit=NULL,upper.limit=NULL,
                            optimmethod="Nelder-Mead",
                            maxit=1000, hessian=FALSE, nwg=70,...){

  # Handle cases where w is not a matrix with row sums equal to 1
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

  setup<-modeloptions(model,customgauge,init.val,lower.limit,
                      upper.limit,d=d,...)

  if(!is.null(init.val)){setup$init.val<-init.val}
  if(!fixshape){
    setup$init.val=c(d,setup$init.val)
    setup$lower.limit=c(0,setup$lower.limit)
    setup$upper.limit=c(100*d,setup$upper.limit)
    if(!is.null(setup$par1ind)){setup$par1ind<-setup$par1ind+1}
    if(!is.null(setup$par2ind)){setup$par2ind<-setup$par2ind+1}
  }

  if(setup$add.gauge&&d>2){wgrid<-create.wgrid(nwg)
  } else{wgrid<-NULL}

  opt<-optim(rw.gamma.lik, par=setup$init.val,r=r,w=w,r0w=r0w,
             gfun=setup$gfun,
             add.gauge=setup$add.gauge,gauge1=setup$gauge1,gauge2=setup$gauge2,
             par1ind=setup$par1ind,par2ind=setup$par2ind,
             lower.limit=setup$lower.limit,upper.limit=setup$upper.limit,
             d=d,wgrid=wgrid,fixshape=fixshape,method=optimmethod,
             control=list(maxit=maxit),...)
  z<-list()
  if(fixshape){z$fixshape=TRUE;z$mle<-opt$par;z$shape<-d
  }else{z$fixshape=FALSE;z$mle<-opt$par[2:length(opt$par)];z$shape<-opt$par[1]}
  z$nll=opt$value
  z$convergence<-opt$conv
  if(hessian==FALSE){z$cov<-NULL
  } else{z$cov=solve(opt$hessian)}
  z$setup=setup
  z$data=list(r=r,w=w,r0w=r0w)
  z$class="parametric"
  z$dimension=d

  print(list(shape=z$shape,mle=z$mle,nll=z$nll,conv=z$conv))
  return(z)
}


############################################################
############################################################
additivegauge.scaling<-function(gauge1,gauge2,par1,par2,weight,d,wgrid=NULL)
{
  dummy<-function(x)
  {
    weight*gauge1(x,par=par1) + (1-weight)*gauge2(x,par=par2)
  }
  if(d==2)
  {
    wseq1<-seq(0,1,len=501)
    ms<-max(wseq1/apply(cbind(wseq1,1-wseq1),1,dummy))
  }else
  {
    den<-apply(wgrid,1,dummy)
    ms<-apply(wgrid/den,2,max)
  }
  return(ms)
}

# ms is used as input in additivegauge.rescale to get the rescaled gauge function

additivegauge.rescale<-function(x,gauge1,gauge2,par1,par2,weight,ms)
{
  return(weight*gauge1(ms*x,par=par1) + (1-weight)*gauge2(ms*x,par=par2))
}

####################################################################
####################################################################
#' Negative log likelihood function for truncated gamma model for R|W

#' @param psi vector of parameters. First component is the gamma shape parameter
#' @param r Values of R|W=w that exceed the threshold r0w
#' @param w Values of W that correspond to the threshold exceedances of R
#' @param r0w The threshold function r_0(w) for each value of W supplied
#' @param gfun gauge function to be used in the likelihood. Can use inbuilt examples or specify own function of the form gauge<-function(xy,par){}. Leave blank if want to mix gauges additively.
#' @param add.gauge logical. If TRUE, then we additively mix (with numerical rescaling) the gauge functions specified by gauge1 and gauge2
#' @param gauge1 first gauge function for additive mixing. Use if add.gauge=TRUE
#' @param gauge2 second gauge function for additive mixing. Use if add.gauge=TRUE
#' @param par1ind vector of indices specifying the position of the parameters for gauge1 in psi
#' @param par2ind vector of indices specifying the position of the parameters for gauge2 in psi
#' @param pos.par logical. Set to TRUE if all parameters are positive (avoids some warnings in optimization)
#' @param lower.limit optional vector of lower limits on parameters (should be same length as psi)
#' @param upper.limit optional vector of upper limits on parameters (should be same length as psi)
#' @param d dimension
#' @param wgrid matrix containing points that cover the d-1 simplex. Used only if add.gauge=T, to perform numerical rescaling
#' @param fixshape  Logical. If TRUE, then the shape parameter of the truncated gamma distribution is fixed to the dimension d, otherwise this is estimated.

#' @return negative log-likelihood value
#'



rw.gamma.lik<-function(psi,r,w,r0w,gfun,add.gauge=FALSE,gauge1,
                       gauge2,par1ind,par2ind,lower.limit,
                       upper.limit,d,wgrid=NULL,fixshape,...)
{
  if(!is.null(lower.limit)){if(any(psi<lower.limit)){return(10e10)}}
  if(!is.null(upper.limit)){if(any(psi>upper.limit)){return(10e10)}}

  if(fixshape)
  {
    if(add.gauge)
    {
      ms<-additivegauge.scaling(gauge1=gauge1,gauge2=gauge2,
                                weight=psi[length(psi)],par1=psi[par1ind],
                                par2=psi[par2ind],d=d,wgrid=wgrid)
      rate<-apply(w,1,additivegauge.rescale,par1=psi[par1ind],
                  par2=psi[par2ind],ms=ms,gauge1=gauge1,gauge2=gauge2,
                  weight=psi[length(psi)],...)
      ll1<-dgamma(r,shape=d,rate=rate,log=T)
      ll2<-log(pgamma(r0w,shape=d,rate=rate,lower.tail=F))
    }
    else{
      rate<-apply(w,1,gfun,par=psi,...)
      ll1<-dgamma(r,shape=d,rate=rate,log=T)
      ll2<-log(pgamma(r0w,shape=d,rate=rate,lower.tail=F))
    }
  } else{
    if(add.gauge)
    {
      ms<-additivegauge.scaling(gauge1=gauge1,gauge2=gauge2,
                                weight=psi[length(psi)],par1=psi[par1ind],
                                par2=psi[par2ind],d=d,wgrid=wgrid)
      rate<-apply(w,1,additivegauge.rescale,par1=psi[par1ind],
                  par2=psi[par2ind],ms=ms,gauge1=gauge1,gauge2=gauge2,
                  weight=psi[length(psi)],...)
      ll1<-dgamma(r,shape=psi[1],rate=rate,log=T)
      ll2<-log(pgamma(r0w,shape=psi[1],rate=rate,lower.tail=F))
    }
    else{
      rate<-apply(w,1,gfun,par=psi[2:length(psi)],...)
      ll1<-dgamma(r,shape=psi[1],rate=rate,log=T)
      ll2<-log(pgamma(r0w,shape=psi[1],rate=rate,lower.tail=F))
    }
  }
  return(-sum(ll1)+sum(ll2))
}




# wgrid should be a grid of points (matrix with d columns) covering the (d-1) simplex S_(d-1)

create.wgrid<-function(n,d)
{
  wseq<-seq(0,1,len=n)
  wl<-list()
  for(j in 1:(d-1))
  {
    wl[[j]]<-wseq
  }
  wgrid<-expand.grid(wl)
  wgrid<-cbind(wgrid,1-apply(wgrid,1,sum))
  neg.ind<-wgrid[,d]<0
  wgrid<-wgrid[!neg.ind,]
  return(wgrid)
}






