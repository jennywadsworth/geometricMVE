#' Fit piecewise-linear truncated gamma model for R|W, and a piecewise-linear model for W

#' @param r Values of R|W=w that exceed the threshold r0w
#' @param w Values of W that correspond to the threshold exceedances of R
#' @param locs The reference angles (nodes) at which parameters are estimation
#' @param init.val vector of initial parameter values for optimization (when customgauge used)
#' @param fixshape Logical. If TRUE, then the shape parameter of the truncated gamma distribution is fixed to the dimension d, otherwise this is estimated.
#' @param W.fit Logical. If TRUE, then the angular model is fitted.
#' @param joint.fit Logical. If TRUE, the model for R|W and W share the same piecewise-linear gauge function parameters
#' @param joint.fit the optimization method (default is "BFGS")
#' @param pen.const the L2 penalty strength (default is 1)
#' @param bount.fit Logical. If TRUE, the gauge function is estimated such that its unit level set has coordinatewise maximum 1
#'
#' @return list containing elements "mle", "nllh", "convergence", and additionally "cov" if argument hessian=T is supplied to optim
#'
#' @export
fit.geometric.pwl = function(r,w,r0w,locs=NULL,
                             init.val=NULL,fixshape=FALSE,W.fit=FALSE,joint.fit=FALSE,
                             method="BFGS",pen.const=1,bound.fit=FALSE,hessian=FALSE){
  if(any(w<0)){stop("Invalid values of W")}
  if(is.vector(w)){
    w<-cbind(w,1-w)
  } else{
    sw<-apply(w,1,sum)
    if(all(sw<1)){
      w<-cbind(w,1-sw)
    }
  }
  
  if(dim(w)[2]==2){
    if(is.null(locs)){
      locs = seq(0,1,length.out=11)
    }
    fit = fit.pwlin.2d(r=r,w=w[,1],r0w=r0w,locs=locs,init.val=init.val,fixshape=fixshape,
                       W.fit=W.fit,joint.fit=joint.fit,
                       method=method,pen.const=pen.const,bound.fit=bound.fit,hessian=hessian)
    gfun = function(x,par) gfun.2d(x,par,ref.angles=locs)
  } else if(dim(w)[2]>2) {
    if(is.null(locs)){
      stop("define reference angles.")
    }
    fit = fit.pwlin(r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,fixshape=fixshape,
                       W.fit=W.fit,joint.fit=joint.fit,
                       method=method,pen.const=pen.const,bound.fit=bound.fit,hessian=hessian)
    gfun = function(x,par) gfun.pwl(x,par,ref.angles=locs)
  }
  
  fit$fixshape=fixshape
  if(hessian==FALSE){
    fit$cov<-NULL
  } else{
    fit$cov=solve(fit$hessian)
  }
  fit$setup$gfun = gfun
  fit$setup$add.gauge = FALSE
  fit$data = list(r=r,w=w,r0w=r0w)
  fit$class="pwl"
  fit$dimension=dim(w)[2]
  print(list(shape=fit$shape,mle=fit$mle,nll=fit$nll,conv=fit$conv))
  return(fit)
}


nll.pwlin.2d = function(psi, r, r0w, w.adj.angles, locs, 
                        W.fit, joint.fit, 
                        pos.par=TRUE, pen.const=0, fixshape){

  # if(!fixshape & W.fit & !joint.fit){
  #   stop("No shape parameter estimation if only fitting angular model.")
  # }
  
  if(!fixshape){
    shape=psi[1]
    psi=psi[-1]
    if(shape < 0){
      return(1e+11)
    }
  } else {
    shape=2
  }

  # ensure we have no negative parameters
  if(pos.par) {
    if (any(psi <= 0)){
      return(1e+11)
    }
  }

  # rate <- sapply(w, gfun, par = psi)
  if(W.fit & joint.fit){
    # exceedance radii and angles
    rate = pwlin.g.vals.2d(w.adj.angles,psi,locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    G.vol = G.vol.2d(psi,locs)
    llw = -log(2*G.vol) -2*log(rate)
    NLL = -sum(ll1) + sum(ll2) - sum(llw)
  } else if(!W.fit & !joint.fit){
    # exceedance radii only
    rate = pwlin.g.vals.2d(w.adj.angles,psi,locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    NLL = -sum(ll1) + sum(ll2)
  } else if(W.fit & !joint.fit){
    # exceedance angles only
    psi = c(1,psi)  # account for redundancy
    rate = pwlin.g.vals.2d(w.adj.angles,psi,locs)
    G.vol = G.vol.2d(psi,locs)
    llw = -log(2*G.vol) -2*log(rate)
    NLL = -sum(llw)
  }
  if(pen.const>0){
    # pen = abs(diff(psi))
    # pen[c(1,length(pen))] = 2*pen[c(1,length(pen))]  # these appear twice in the penalty
    locs.cart = cbind(locs,1-locs)
    scaled.locs.cart = psi*locs.cart
    coplanar.vecs = apply(scaled.locs.cart,2,diff)
    norm.vecs = cbind(coplanar.vecs[,-1],-coplanar.vecs[,-2])
    grad.g = norm.vecs / apply(scaled.locs.cart[-nrow(scaled.locs.cart),] * norm.vecs, 1,sum)  # gradient of g in each region (is constant/linear in each face)
    pen = apply(grad.g,2,diff)
    # if(pen.norm=="2"){
    pen = apply(pen,1,function(vec) sum(vec^2))  # mean or sum? Sum makes sense.
    # } else if(pen.norm=="1"){
    #   pen = apply(pen,1,function(vec) sum(abs(vec)))
    # }

    # if(pen.adj){
    #   pen.multiplier = sapply(1:length(pen),function(r.idx){
    #     w.low=locs[r.idx]
    #     w.up=locs[r.idx+2]
    #     wexc.vec = sapply(w.adj.angles,function(lst) lst$w)
    #     num = sum(wexc.vec>=w.low & wexc.vec<=w.up)
    #     if(num==0){
    #       return(1)
    #     } else {
    #       return(num)
    #     }
    #   })
    #   pen = pen/pen.multiplier
    # }

    pen = mean(pen)  # used to be sum, but mean makes sense
    return(NLL+(pen.const*pen))
  } else {
    return(NLL)
  }
}

opt.pwl.2d = function(NLL, r, r0w, w, locs,
                      init.val=NULL, fixshape=TRUE, W.fit=FALSE,
                      joint.fit=FALSE,method="BFGS",
                      pen.const=0,hessian=FALSE,...){
  # w.adj.angles -> output of "which.adj.angles.2d"

  # if(is.null(init.val)){
  #   init.val = rep(1,length(locs))
  # }
  w.adj.angles = which.adj.angles.2d(w,locs)
  if(fixshape & !W.fit & !joint.fit){
    # only fit the conditional radial model
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, 
                 W.fit=F,joint.fit=F,pen.const=pen.const,fixshape=fixshape,
                 control=list(maxit=1e6),method=method,hessian=hessian,...)
    #method="L-BFGS-B",gr=NULL,lower=rep(0.01,length(locs)),upper=rep(10,length(locs))
    opt$W.par = NULL  # didn't model angles
    opt$init.val = init.val
    opt$shape=2
  } else if(fixshape & W.fit & joint.fit){
    # only fit the radial model and angular model togethr
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, 
                 W.fit=T,joint.fit=T,pen.const=pen.const,fixshape=fixshape,
                 control=list(maxit=1e6),method=method,hessian=hessian,...)
    opt$W.par = opt$par
    opt$init.val = init.val
    opt$shape=2
  } else if(W.fit & !joint.fit){
    # only fit the angular model together, taking into account the redundancy
    init.val = init.val[-1]  # remove redundancy, fix the first parameter to 1.
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,  
                 W.fit=T,joint.fit=F,pen.const=pen.const,fixshape=fixshape,
                 control=list(maxit=1e6),method=method,hessian=hessian,...)
    opt$W.par = opt$par
    opt$W.par = c(1,opt$W.par)  # account for redundancy
    opt$par = NULL
    opt$init.val = c(1,init.val)
    opt$shape=NULL
  } else if(!fixshape & !W.fit & !joint.fit){
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs, 
                 W.fit=W.fit,joint.fit=W.fit,pen.const=pen.const,fixshape=fixshape,
                 control=list(maxit=1e6),method=method,hessian=hessian,...)
    #method="L-BFGS-B",gr=NULL,lower=rep(0.01,length(locs)),upper=rep(10,length(locs))
    opt$W.par = NULL  # didn't model angles
    opt$init.val = init.val
    opt$shape = opt$par[1]
    opt$par = opt$par[-1]
  } else if(!fixshape & W.fit & joint.fit){
    # only fit the radial model and angular model togethr
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,  
                 W.fit=T,joint.fit=T,pen.const=pen.const,fixshape=fixshape,
                 control=list(maxit=1e6),method=method,hessian=hessian,...)
    opt$init.val = init.val
    opt$shape = opt$par[1]
    opt$par = opt$par[-1]
    opt$W.par = opt$par
  }
  # return(list(mle = opt$par,mle.W = opt$W.par, shape=opt$shape, nll = opt$value, convergence = opt$conv,
  #             aic = opt$aic, init.val = opt$init.val))
  return(list(mle = opt$par, mle.W = opt$W.par, shape=opt$shape,
              nllh = opt$value, convergence = opt$conv,
              init.val = opt$init.val, hessian=opt$hessian))
}


fit.pwlin.2d = function(r,w,r0w,locs,
                     init.val=NULL,fixshape=TRUE,W.fit=FALSE,joint.fit=FALSE,
                     method="BFGS",pen.const=NULL,bound.fit=FALSE,hessian=FALSE){

  if(is.null(init.val) & fixshape){
    init.val = rep(1,length(locs))
  } else if(is.null(init.val) & !fixshape){
    init.val = rep(1,length(locs)+1)
  }

  if(is.null(pen.const)){
    message("searching for gradient penalty constant...")
    pen.const = get_pen_const_2d(r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,W.fit=W.fit,joint.fit=joint.fit)
    message(paste("Fitting with gradient penalty strength",pen.const))
  }

  gfun=function(x,par) {gfun.2d(x,par,ref.angles=locs)}

  t1 = Sys.time()

  nnodes=length(locs)
  w.adj.angles = which.adj.angles.2d(angles=w, locs)
  opt = opt.pwl.2d(NLL=nll.pwlin.2d, r=r, r0w=r0w,w=w,locs=locs,
                   init.val=init.val,W.fit=W.fit,joint.fit=joint.fit,
                   pen.const=pen.const,fixshape=fixshape,hessian=hessian)

  if(bound.fit){
    if(W.fit & !joint.fit){
      stop("No need to bound angular parameters if not jointly fitting.")
    }

    shape.val = opt$shape
    mle=opt$mle

    lik.custom = function(psi, r, r0w, w.adj.angles,
                          locs, 
                          fixed.pars, fixed.pars.idx,
                          W.fit,joint.fit,pen.const,fixshape,hessian){

      if(fixshape){
        par.full = numeric(length(locs))
        par.full[fixed.pars.idx] = fixed.pars
        par.full[-fixed.pars.idx] = psi
      } else {
        par.full = numeric(length(locs)+1)
        par.full[fixed.pars.idx+1] = fixed.pars
        par.full[-(fixed.pars.idx+1)] = psi
      }

      return(nll.pwlin.2d(psi=par.full,r=r,r0w=r0w,w.adj.angles=w.adj.angles,
                          locs=locs,
                          W.fit=W.fit,joint.fit=joint.fit,pen.const=pen.const,
                          fixshape=fixshape))
    }



    # stop("generalise this to use the 1 or 2 norm")

    fixed.pars.idx = c()

    keep.fitting = any(round(as.numeric(apply(cbind(locs,1-locs)/sapply(locs,gfun,par=mle),2,max)),6)!=1)
    
    max.iters=100
    iters=1
    while(keep.fitting){
      if(0.5 %in% locs){

        idx.x  = locs >= 0.5
        idx.y  = locs <= 0.5
        idx.xy = locs == 0.5
        cm.x = max(locs/sapply(locs,gfun,par=mle))
        cm.y = max((1-locs)/sapply(locs,gfun,par=mle))
        cm.xy = 0.5/sapply(0.5,gfun,par=mle)
        if(cm.x==cm.xy & cm.y==cm.xy){
          which.cm = which(locs==0.5)
          mle[which.cm] = mle[which.cm]/cm.xy  # re-scale
          if(is.null(fixed.pars.idx)){
            fixed.pars.idx = which.cm
          } else {
            fixed.pars.idx = sort(c(fixed.pars.idx,which.cm))
          }

        } else {
          # max isn't on the diagonal
          # find the componentiwse max, and re-scale the corresponding MLEs
          idx.x = locs > 0.5
          idx.y = locs < 0.5
          unitg.x = locs/sapply(locs,gfun,par=mle)
          unitg.x[!idx.x] = 0
          unitg.y = (1-locs)/sapply(locs,gfun,par=mle)
          unitg.y[!idx.y] = 0
          which.cm.x = which.max(unitg.x)
          which.cm.y = which.max(unitg.y)
          mle[which.cm.x] = mle[which.cm.x]/(unitg.x[which.cm.x])
          mle[which.cm.y] = mle[which.cm.y]/(unitg.y[which.cm.y])

          # which mle's are re-scaled
          if(is.null(fixed.pars.idx)){
            fixed.pars.idx = sort(c(which.cm.x,which.cm.y))
          } else {
            fixed.pars.idx = sort(c(fixed.pars.idx,which.cm.x,which.cm.y))
          }
        }
      } else {
        # find the componentiwse max, and re-scale the corresponding MLEs
        idx.x = locs > 0.5
        idx.y = locs < 0.5
        unitg.x = locs/sapply(locs,gfun,par=mle)
        unitg.x[!idx.x] = 0
        unitg.y = (1-locs)/sapply(locs,gfun,par=mle)
        unitg.y[!idx.y] = 0
        which.cm.x = which.max(unitg.x)
        which.cm.y = which.max(unitg.y)
        mle[which.cm.x] = mle[which.cm.x]/(unitg.x[which.cm.x])
        mle[which.cm.y] = mle[which.cm.y]/(unitg.y[which.cm.y])

        # which mle's are re-scaled
        if(is.null(fixed.pars.idx)){
          fixed.pars.idx = sort(c(which.cm.x,which.cm.y))
        } else {
          fixed.pars.idx = unique(sort(c(fixed.pars.idx,which.cm.x,which.cm.y)))
        }
      }
      init.vals=mle[-fixed.pars.idx]
      if(!fixshape){
        init.vals = c(shape.val,init.vals)
      }

      # re-fit the model
      if(!W.fit & !joint.fit){
       opt2 = opt.pwl.2d(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                          W.fit=W.fit,joint.fit=joint.fit,pen.const=pen.const,
                          fixshape=fixshape,
                          method=method,
                          fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,
                          hessian=hessian)
        opt2$W.par = NULL  # didn't model angles
      } else if(W.fit & joint.fit){
       opt2 = opt.pwl.2d(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                          W.fit=W.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                          fixshape=fixshape,
                          fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx,
                          hessian=hessian)
        opt2$W.par = opt2$mle
      } else {
        stop("ERROR in likelihood setup.")
      }

      mle[-fixed.pars.idx] = opt2$mle
      shape.val = opt2$shape
      
      locs.cart = cbind(locs,1-locs)
      keep.fitting = any(round(as.numeric(apply(cbind(locs,1-locs)/sapply(locs,gfun,par=mle),2,max)),6)!=1)
    
      if(!keep.fitting){
        # message(paste("comp. min:",paste(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,min)),collapse=" ")))
        # message(paste("comp. max:",paste(as.numeric(apply(locs.cart/sapply(locs,gfun,par=mle),2,max)),collapse=" ")))
        break
      }

      iters=iters+1
      if(iters>max.iters){
        message("Max iters in bounding algorithm reached.")
        break
      }

    }
    t2 = Sys.time()

    if(joint.fit){
      mle.W = mle
    } else {
      mle.W = NULL
    }

    return(list(mle =mle,
                mle.W = mle.W,
                shape=opt2$shape,  # made this opt2 instead of opt
                nll = opt2$nllh,
                convergence = opt2$conv,
                hessian=opt2$hessian))
  } else {
    t2 = Sys.time()

    return(list(mle =opt$mle,
                mle.W = opt$mle.W,
                shape=opt$shape,
                nll = opt$nllh,
                convergence = opt$conv,
                hessian=opt$hessian))
  }
}

get_pen_const_2d = function(r,w,r0w,locs,init.val,W.fit=FALSE,joint.fit=FALSE){
  if(W.fit & !joint.fit){
    pen.consts = seq(0,50,length.out=50)
  } else {
    pen.consts = seq(0,4,length.out=50)
  }

  k.folds = 4
  grps = sample(1:k.folds,length(r),replace=T)

  pen.consts.scores = rep(NA,length(pen.consts))

  iter=1
  for(pen.const.idx in c(1:length(pen.consts))){

    pen.const = pen.consts[pen.const.idx]
    pen.consts.score = c()

    for(K in 1:k.folds){
      # K = 1
      r.fitting = r[grps!=K]
      w.fitting = w[grps!=K]
      r0w.fitting = r0w[grps!=K]
      r.eval = r[grps==K]
      w.eval = w[grps==K]
      r0w.eval = r0w[grps==K]

      is_error <- FALSE
      tryCatch({
        mod.fit = fit.pwlin.2d(r=r.fitting,r0w=r0w.fitting,w=w.fitting,
                                 locs=locs,pen.const=pen.const,method="BFGS",
                                 init.val=init.val,
                               W.fit=W.fit,joint.fit=joint.fit,bound.fit=FALSE)
      },error=function(e){
        is_error <<- TRUE
      })
      if(is_error) {
        # print(pen.const)
        # pen.consts.scores[pen.const.idx] = 1e10
        next
      }

      w.eval.adj.angles = which.adj.angles.2d(angles=w.eval,locs=locs)

      if(W.fit & !joint.fit){
        psi.vals = mod.fit$mle.W
      } else {
        psi.vals = mod.fit$mle
      }

      pen.consts.score = c(pen.consts.score,
        nll.pwlin.2d(psi=psi.vals, r=r.eval, r0w=r0w.eval,
                                       w.adj.angles=w.eval.adj.angles,
                                       locs=locs,
                                       W.fit=W.fit,joint.fit=joint.fit))
    }
    pen.consts.scores[pen.const.idx] = (1/k.folds)*sum(pen.consts.score,na.rm=T)
  }
  return(pen.consts[which.min(pen.consts.scores)])
}

######################################################################

# general d-dimensional model fitting for when w is on the simplex.

pwl.L2.pen = function(psi,locs,del.tri,ij.couples.list){
  require(geometry)
  tri = del.tri$tri

  locs.scaled = locs*psi

  num.cols = dim(locs)[2]

  grad.g = lapply(1:nrow(tri),function(row.idx){
    idx = tri[row.idx,]
    idx.fix = idx[1]
    idx.rest = idx[-1]
    coplanar.mat = do.call(rbind,lapply(idx.rest, function(ii){
      return(locs.scaled[idx.fix,] - locs.scaled[ii,])
    }))
    norm.vec = suppressWarnings({
      c(1,-1)*sapply(c(1:num.cols),function(i){det(coplanar.mat[,-i])})
    })
    grad.g.val = norm.vec / sum(norm.vec * locs.scaled[idx.fix,])
    return(grad.g.val)
  })
  grad.g.mat = do.call(rbind,grad.g)
  pen = sapply(1:nrow(locs), function(l){
    ij.couples = ij.couples.list[[l]]
    if(all(is.na(ij.couples))){
      return(NA)
    }
    pen.val.l = apply(ij.couples, 1, function(ij){
      diff(grad.g.mat[ij,])
    })
    pen.val.l = apply(t(pen.val.l)^2,1,sum)
    return(mean(pen.val.l))
  })
  return(mean(pen,na.rm=T))
}

nll.pwlin = function(psi, r, r0w, w.adj.angles, fixshape, del.tri, locs,
                     ij.couples.list,W.fit, joint.fit, pen.const=0, pos.par=TRUE){

  num.cols = length(w.adj.angles[[1]]$w)

  if(!fixshape){
    shape=psi[1]
    psi=psi[-1]
    if(shape < 0){
      return(1e+11)
    }
  } else {
    shape=num.cols
  }

  # ensure we have no negative parameters
  if(pos.par) {
    if (any(psi <= 0)){
      return(1e+11)
    }
  }

  if(W.fit & joint.fit){
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    G.vol = G.vol(psi,locs)
    llw = -log(num.cols*G.vol) -num.cols*log(rate)
    NLL = -sum(ll1) + sum(ll2) - sum(llw)
  } else if(!W.fit & !joint.fit){
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    ll1 <- dgamma(r, shape = shape, rate = rate, log = T)
    ll2 <- log(pgamma(r0w, shape = shape, rate = rate, lower.tail = F))
    NLL = -sum(ll1) + sum(ll2)
  } else if(W.fit & !joint.fit){
    psi = c(1,psi)  # account for redundancy
    rate <- pwlin.g.vals(w.adj.angles=w.adj.angles,par=psi,par.locs=locs)
    G.vol = G.vol(psi,locs)
    llw = -log(num.cols*G.vol) -num.cols*log(rate)
    NLL = -sum(llw)
  }
  if(pen.const>0){
    L2.pen.val = pwl.L2.pen(psi=psi,locs=locs,del.tri=del.tri,ij.couples.list=ij.couples.list)
    NLL = NLL+(pen.const*L2.pen.val)
  }
  return(NLL)
}

opt.pwl = function(NLL, r, r0w, w, locs, w.adj.angles, del.tri, ij.couples.list,
                     init.val=NULL, fixshape=TRUE, W.fit=FALSE,
                     joint.fit=FALSE,method="BFGS",hessian=FALSE,
                     pen.const=0,...){

  # w.adj.angles = which.adj.angles(w,locs)

  num.cols = length(w.adj.angles[[1]]$w)
  # del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)
  # ij.couples.list = ij.couples(locs)

  if(fixshape & !W.fit & !joint.fit){
    # fit the conditional radial model only
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 ij.couples.list=ij.couples.list,
                 W.fit=F,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixshape=fixshape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,hessian=hessian...)
    opt$W.par = NULL  # didn't model angles
    opt$init.val = init.val
    opt$shape=num.cols
  } else if(fixshape & W.fit & joint.fit){
    # fit the radial model and angular model together
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 ij.couples.list=ij.couples.list,
                 W.fit=T,joint.fit=T,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixshape=fixshape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,hessian=hessian,...)
    opt$W.par = opt$par
    opt$init.val = init.val
    opt$shape=num.cols
  } else if(W.fit & !joint.fit){
    # fit the angular model only
    init.val = init.val[-1]  # remove redundancy, fit the first parameter to 1.
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 ij.couples.list=ij.couples.list,
                 W.fit=T,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixshape=fixshape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,hessian=hessian,...)
    opt$W.par = opt$par
    opt$W.par = c(1,opt$W.par)  # account for redundancy
    opt$par = NULL
    opt$init.val = c(1,init.val)
    opt$shape=NULL
  } else if(!fixshape & !W.fit & !joint.fit){
    # fit the conditional radial model only
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 ij.couples.list=ij.couples.list,
                 W.fit=F,joint.fit=F,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixshape=fixshape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,hessian=hessian,...)
    opt$init.val = init.val
    opt$shape = opt$par[1]
    opt$par = opt$par[-1]
  } else if(!fixshape & W.fit & joint.fit){
    # fit the radial model and angular model together
    opt <- optim(NLL, par = init.val, r = r,
                 r0w = r0w, w.adj.angles=w.adj.angles, locs=locs,
                 ij.couples.list=ij.couples.list,
                 W.fit=T,joint.fit=T,
                 pen.const=pen.const,
                 del.tri=del.tri,
                 fixshape=fixshape,
                 control=list(maxit=1e6,reltol=1e-5),method=method,hessian=hessian,...)
    opt$init.val = init.val
    opt$shape = opt$par[1]
    opt$par = opt$par[-1]
    opt$W.par = opt$par
  }

  opt$init.val = init.val#c(3,init.val)
  return(list(mle = opt$par, mle.W = opt$W.par, shape=opt$shape,
              nllh = opt$value, convergence = opt$conv,
              init.val = opt$init.val,
              hessian=opt$hessian))
}

fit.pwlin = function(r,w,r0w,locs,
                     init.val=NULL,fixshape=TRUE,W.fit=FALSE,joint.fit=FALSE,
                     method="BFGS",pen.const=NULL,bound.fit=FALSE,hessian=FALSE){

  if(W.fit & !joint.fit & bound.fit){
    stop("No need to bound angular parameters if not jointly fitting.")
  }

  t1 = Sys.time()

  if(is.null(init.val) & fixshape){
    init.val = rep(1,nrow(locs))
  } else if(is.null(init.val) & !fixshape){
    init.val = rep(1,nrow(locs)+1)
  }

  if(is.null(pen.const)){
    message("searching for gradient penalty constant...")
    pen.const = get_pen_const(r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,W.fit=W.fit,joint.fit=joint.fit)
    message(paste("Fitting with gradient penalty strength",pen.const))
  }

  num.cols = dim(w)[2]

  w.adj.angles = which.adj.angles(angles=w,locs=locs)
  del.tri = PWLExtremes::delaunayn(p=locs[,-num.cols], output.options=TRUE)
  ij.couples.list = ij.couples(locs)

  # fixed.pars.idx=NULL
  opt <- opt.pwl(NLL=nll.pwlin,r=r,w=w,r0w=r0w,locs=locs,init.val=init.val,
                 w.adj.angles=w.adj.angles, del.tri=del.tri, ij.couples.list=ij.couples.list,
                 fixshape=fixshape,W.fit=W.fit,joint.fit=joint.fit,
                 pen.const=pen.const,method=method,hessian=hessian)

  if(bound.fit){

    shape.val = opt$shape
    mle = opt$mle

    lik.custom = function(psi, r, r0w, w.adj.angles, locs = locs, ij.couples.list=ij.couples.list,
                          fixed.pars,fixed.pars.idx, fixshape, W.fit, joint.fit, pos.par = TRUE,
                          pen.const=0,del.tri=del.tri){

    if(fixshape){
      par.full = numeric(nrow(locs))
      par.full[fixed.pars.idx] = fixed.pars
      par.full[-fixed.pars.idx] = psi
    } else {
      par.full = numeric(nrow(locs)+1)
      par.full[fixed.pars.idx+1] = fixed.pars
      par.full[-(fixed.pars.idx+1)] = psi
    }
 
      return(nll.pwlin(psi=par.full,r=r,r0w=r0w,w.adj.angles=w.adj.angles,
                       locs=locs, ij.couples.list=ij.couples.list, W.fit=W.fit,
                       joint.fit=joint.fit, pos.par=TRUE,
                       pen.const=pen.const,del.tri=del.tri,fixshape=fixshape))

    }
    fixed.pars.idx = c()
    is.bounded=F   # kind of like a dummy variable, will terminate the loop inside
    while(!is.bounded){
      # print(i)

      # evaluate at nodes
      unitg = locs/gfun.pwl(x=locs,par=mle,ref.angles=locs)#apply(locs,1,gfun.4d,par=mle,ref.angles=locs)  # need to round here

      # what is max along nodes
      max.vals = apply(unitg,2, max)
      which.max.vals = apply(unitg,2, which.max)  # row values

      max.df = data.frame(coord=c(1:num.cols),max=max.vals,loc=which.max.vals)
      max.lst = lapply(1:nrow(locs),function(loc){
        cond = max.df$loc==loc
        if(all(cond==F)){
          return(NULL)
        } else {
          return(list(loc=loc,
                      val=max(max.df$max[cond])))
        }
      })
      max.lst = max.lst[lapply(max.lst,length)>0]

      if(all(round(max.vals,num.cols)==1)){
        break  # algorithm complete, all comp. max. vals are 1
      }
      for(lst in max.lst){
        fixed.pars.idx = unique(sort(c(fixed.pars.idx,lst$loc)))
        mle[lst$loc] = mle[lst$loc]/lst$val
      }

      init.vals=mle[-fixed.pars.idx]
      if(!fixshape){
        init.vals = c(shape.val,init.vals)
      }

      if(nrow(locs)==length(fixed.pars.idx)){
        break  # we fixed all the parameters
      }

      if(!W.fit & !joint.fit){
        # fit the conditional radial model only
        init.vals = init.val[-fixed.pars.idx]
        opt2 = opt.pwl(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                      w.adj.angles=w.adj.angles, del.tri=del.tri, ij.couples.list=ij.couples.list,
                      W.fit=W.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,
                      hessian=hessian,
                      fixshape=fixshape,
                      fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx)
        opt2$W.par = NULL  # didn't model angles
      } else if(W.fit & joint.fit){
        # fit the radial model and angular model together
        init.vals = init.val[-fixed.pars.idx]
        opt2 = opt.pwl(NLL=lik.custom,r=r,w=w,r0w=r0w,locs=locs,init.val=init.vals,
                       w.adj.angles=w.adj.angles, del.tri=del.tri, ij.couples.list=ij.couples.list,
                       W.fit=W.fit,joint.fit=joint.fit,pen.const=pen.const,method=method,hessian=hessian,
                       fixshape=fixshape,
                       fixed.pars=mle[fixed.pars.idx], fixed.pars.idx=fixed.pars.idx)
        opt2$W.par = opt2$mle
      } else {
        stop("ERROR in likelihood setup")
      }
      mle[-fixed.pars.idx] = opt2$mle
      shape.val = opt2$shape
      # i=i+1
    }

    t2 = Sys.time()

    if(joint.fit){
      mle.W = mle
    } else {
      mle.W = NULL
    }
    return(list(mle = mle, mle.W = mle.W, shape=shape.val,
                nll = opt2$nllh, convergence = opt2$conv,
                hessian=opt2$hessian))
  } else {
    t2 = Sys.time()
    list(mle = opt$mle, mle.W = opt$mle.W, shape=opt$shape,
         nll = opt$nllh, convergence = opt$convergence, hessian=opt$hessian)
  }
}

get_pen_const = function(r,w,r0w,locs,init.val,W.fit=FALSE,joint.fit=FALSE){
  if(W.fit & !joint.fit){
    pen.consts = seq(0,50,length.out=15)
  } else {
    pen.consts = seq(0,4,length.out=15)
  }

  k.folds = 4
  grps = sample(1:k.folds,length(r),replace=T)

  pen.consts.scores = rep(NA,length(pen.consts))

  iter=1
  for(pen.const.idx in c(1:length(pen.consts))){
    # print(pen.consts.scores)


    pen.const = pen.consts[pen.const.idx]
    # pen.const.W = pen.consts.W[pen.const.idx]
    # lik.scores.R13 = rep(NA,k.folds)
    pen.consts.score = c()
    for(K in 1:k.folds){
      # K = 1
      r.fitting = r[grps!=K]
      w.fitting = w[grps!=K,]
      r0w.fitting = r0w[grps!=K]
      r.eval = r[grps==K]
      w.eval = w[grps==K,]
      r0w.eval = r0w[grps==K]

      is_error <- FALSE
      tryCatch({
        mod.fit = fit.pwlin(r=r.fitting,r0w=r0w.fitting,w=w.fitting,
                               locs=locs,pen.const=pen.const,method="BFGS",
                               init.val=init.val,
                               W.fit=W.fit,joint.fit=joint.fit,bound.fit=FALSE)
      },error=function(e){
        is_error <<- TRUE
      })
      if(is_error) {
        # print(pen.const)
        # pen.consts.scores[pen.const.idx] = 1e10
        next
      }

      w.eval.adj.angles = which.adj.angles(angles=w.eval,locs=locs)

      if(W.fit & !joint.fit){
        psi.vals = mod.fit$mle.W
      } else {
        psi.vals = mod.fit$mle
      }

      pen.consts.score = c(pen.consts.score,
        nll.pwlin(psi=psi.vals, r=r.eval, r0w=r0w.eval,
                     w.adj.angles=w.eval.adj.angles,
                     locs=locs,
                     W.fit=W.fit,joint.fit=joint.fit))
    }
    pen.consts.scores[pen.const.idx] = (1/k.folds)*sum(pen.consts.score,na.rm=T)
  }
  return(pen.consts[which.min(pen.consts.scores)])
}

