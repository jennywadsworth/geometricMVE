

#' PP diagnostic for fitted truncated gamma model

#' @param fit list output from fit.geometric.par or fit.geometric.pwlin
#' @return PP plot with 95\% tolerance intervals
#'
#' @export

ppdiag<-function(fit)
{
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
  if(!fit$setup$add.gauge){
    rate<-apply(w,1,fit$setup$gfun,par=fit$mle)
    num<-pgamma(r,shape=fit$shape,
                rate=rate,lower.tail = F)
    den<-pgamma(r0w,shape=fit$shape,
                rate=rate,lower.tail=F)
  } else{
    if(!fit$fixshape){par1ind<-fit$setup$par1ind-1;par2ind<-fit$setup$par2ind-1
    } else{par1ind<-fit$setup$par1ind;par2ind<-fit$setup$par2ind}
    ms<-additivegauge.scaling(gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,
                              weight=fit$mle[length(fit$mle)],
                              par1=fit$mle[par1ind],par2=fit$mle[par2ind],d=d)
    rate<-apply(w,1,additivegauge.rescale,par1=fit$mle[par1ind],
                par2=fit$mle[par2ind],ms=ms,gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,
                weight=fit$mle[length(fit$mle)])
    num<-pgamma(r,shape=fit$shape,
                rate=rate,lower.tail = F)
    den<-pgamma(r0w,shape=fit$shape,
                rate=rate,lower.tail = F)
  }
  nr<-length(r)
  plot(c(1:nr)/(nr+1),
       sort(1-num/den),xlab="empirical",ylab="model",pch=20,cex=0.8)
  abline(a=0,b=1)

  Ulow<-sapply(1:nr,function(i){qbeta(0.025,i,nr+1-i)})
  Uup<-sapply(1:nr,function(i){qbeta(0.975,i,nr+1-i)})

  lines(c(1:nr)/(nr+1),Ulow,lty=2)
  lines(c(1:nr)/(nr+1),Uup,lty=2)
}


#' QQ diagnostic for fitted truncated gamma model - puts PP plot on any desired quantile scale (default is standard exponential)

#' @param fit list output from fit.geometric.par or fit.geometric.pwlin
#' @param quantilefn a quantile function to specify the desired scale. Defaults to qexp.
#' @return QQ plot with 95\% tolerance intervals
#'
#' @export

qqdiag<-function(fit,quantilefn=qexp)
{
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
    if(!fit$setup$add.gauge){
      rate<-apply(w,1,fit$setup$gfun,par=fit$mle)
      num<-pgamma(r,shape=fit$shape,
                  rate=rate,lower.tail = F)
      den<-pgamma(r0w,shape=fit$shape,
                  rate=rate,lower.tail=F)
    } else{
      if(!fit$fixshape){par1ind<-fit$setup$par1ind-1;par2ind<-fit$setup$par2ind-1
      } else{par1ind<-fit$setup$par1ind;par2ind<-fit$setup$par2ind}
      ms<-additivegauge.scaling(gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,
                                weight=fit$mle[length(fit$mle)],
                                par1=fit$mle[par1ind],par2=fit$mle[par2ind],d=d)
      rate<-apply(w,1,additivegauge.rescale,par1=fit$mle[par1ind],
                  par2=fit$mle[par2ind],ms=ms,gauge1=fit$setup$gauge1,gauge2=fit$setup$gauge2,
                  weight=fit$mle[length(fit$mle)])
      num<-pgamma(r,shape=fit$shape,
                  rate=rate,lower.tail = F)
      den<-pgamma(r0w,shape=fit$shape,
                  rate=rate,lower.tail = F)
    }
  nr<-length(r)
  plot(quantilefn(c(1:nr)/(nr+1)),
       quantilefn(sort(1-num/den)),xlab="empirical",ylab="model",pch=20,cex=0.8)
  abline(a=0,b=1)

  Ulow<-sapply(1:nr,function(i){qbeta(0.025,i,nr+1-i)})
  Uup<-sapply(1:nr,function(i){qbeta(0.975,i,nr+1-i)})

  lines(quantilefn(c(1:nr)/(nr+1)),quantilefn(Ulow),lty=2)
  lines(quantilefn(c(1:nr)/(nr+1)),quantilefn(Uup),lty=2)
}


########################################################################
########################################################################

#' Plot fitted threshold for d=2,3

#' @param thresh.fit output from fit.thresh
#' @param add logical - add to existing plot?
#' @return plot of fitted threshold over data sample
#'
#' @export
#' 
plotfittedthresh=function(thresh.fit, add=FALSE){
  resolution=100
  d = dim(thresh.fit$w)[2]
  if(is.null(d)){
    d=2
  }
  if(d==2){
    wpts = seq(0,1,len=resolution)
    wpts = cbind(wpts,1-wpts)
    thresh.r = eval.thresh(thresh.fit,wpts)
    if(add==FALSE){plot(thresh.fit$r * cbind(thresh.fit$w,1-thresh.fit$w),xlab="x",ylab="y")}
    lines(thresh.r * wpts,col="red",lwd=2)
  } else if(d==3){
    require(rgl)
    wpts = expand.grid(replicate(d-1,seq(0,1,len=resolution),simplify=F))
    wpts = cbind(wpts,1-apply(wpts,1,sum))
    thresh.r = eval.thresh(thresh.fit,wpts)
    cond = apply(wpts,1,function(ww) any(ww<0))
    thresh.r[cond] = NA
    
    if(add==FALSE)
    {
    plot3d(thresh.fit$r * thresh.fit$w,
           xlab="x1",ylab="x2",zlab="x3")
      }
    surface3d(wpts[,1] * matrix(thresh.r,100,100),
              wpts[,2] * matrix(thresh.r,100,100),
              wpts[,3] * matrix(thresh.r,100,100),
              col="red",alpha=0.4,add=add)
  } else {
    stop("Only supported for d=2,3")
  }
}


#' Plot fitted limit set (unit level set of gauge function) for d=2,3

#' @param fit list output from fit.geometric.par or fit.geometric.pwlin
#' @param resolution plotting resolution - the higher the number, the finer the resolution of the plot. For d=2, the default should suffice, but larger numbers should not take too much longer to plot. For d=3, higher numbers will normally take longer to plot, but give better figures.
#' @param add logical - add plot to an existing figure?
#' @param col character string or number giving the colour for the plot
#' @param xlab,ylab,zlab character string giving axis labels
#' @param lwd,type line width and type arguments to plot for d=2
#' @param polygon logical - include shaded polygon of limit set for d=2?
#' @param polygoncol character string or number giving the colour for the polygon
#' @param alphaval alpha value for shading of the polygon or 3d level set
#' @param ... additional arguments to pass to plot for d=2
#' 
#' @return plot of fitted gauge function
#'
#' @export
#' 

plotfittedgauge<-function(fit,resolution=70,add=FALSE,col="red",xlab="x",ylab="y",zlab="z",lwd=2,type="l",polygon=TRUE,polygoncol=rgb(0.5,0,0,alpha=0.5),alphaval=0.5,...)
{
  if(fit$dimension>3){stop("Cannot plot fitted gauge for d>3")}
    
  if(fit$dimension==2)
  {
    wpts<-seq(0,1,len=10*resolution+1)
    if (fit$setup$add.gauge) {
      ms <- additivegauge.scaling(gauge1 = fit$setup$gauge1, gauge2 = fit$setup$gauge2,
                                  weight = fit$mle[length(fit$mle)], par1 = fit$mle[fit$setup$par1ind], par2 = fit$mle[fit$setup$par2ind],d=fit$dimension)
      gfun <- function(x) {
        if(!fit$fixshape){par1ind<-fit$setup$par1ind-1;par2ind<-fit$setup$par2ind-1
        } else{par1ind<-fit$setup$par1ind;par2ind<-fit$setup$par2ind}
        additivegauge.rescale(x = x, gauge1 = fit$setup$gauge1, gauge2 = fit$setup$gauge2,
                              weight = fit$mle[length(fit$mle)], par1 = fit$mle[par1ind],
                              par2 = fit$mle[par2ind], ms = ms)
      }
      if (add) {
        lines(wpts/apply(cbind(wpts, 1 - wpts), 1, gfun),
              (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun),
              col = col, typ = type, lwd = lwd,...)
        if (polygon) {
          polygon(c(0, wpts/apply(cbind(wpts, 1 - wpts), 1,
                                  gfun), 0), c(0, (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun), 0), col = polygoncol, border = NULL)
        }
      } else {
        plot(wpts/apply(cbind(wpts, 1 - wpts), 1, gfun),
             (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun),
             typ = type, col = col, xlab = xlab,
             ylab = ylab, lwd = lwd,...)
        if (polygon) {
          polygon(c(0, wpts/apply(cbind(wpts, 1 - wpts), 1,
                                  gfun), 0), c(0, (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun), 0), col = polygoncol, border = NULL)
        }
      }
    } else{
    if (add) {
      lines(wpts/apply(cbind(wpts, 1 - wpts), 1, fit$setup$gfun, par = fit$mle),
            (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1,  fit$setup$gfun,
                             par = fit$mle), col = col, typ = type, lwd = lwd,...)
      if (polygon) {
        polygon(c(0, wpts/apply(cbind(wpts, 1 - wpts), 1,
                                fit$setup$gfun, par =  fit$mle), 0), c(0, (1 - wpts)/apply(cbind(wpts,
                                                                                  1 - wpts), 1,  fit$setup$gfun, par = fit$mle), 0), col = polygoncol, border = NULL)
      }
    }
    else {
      plot(wpts/apply(cbind(wpts, 1 - wpts), 1,  fit$setup$gfun, par =  fit$mle),
           (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, fit$setup$gfun,par = fit$mle), typ = type, col = col, xlab = xlab,
           ylab = ylab, lwd = lwd,...)
      if (polygon) {
        polygon(c(0, wpts/apply(cbind(wpts, 1 - wpts), 1,
                                fit$setup$gfun, par = fit$mle), 0), c(0, (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1,  fit$setup$gfun, par = fit$mle), 0), col = polygoncol, border = NULL)
      }
    }
  }
    } else if(fit$dimension==3)
  {
      require(rgl)
      
      if(fit$class=="pwl")
      {
      resolution=resolution+30
      wpts = expand.grid(replicate(fit$dimension-1,seq(0,1,len=resolution),simplify=F))
      wpts = cbind(wpts,1-apply(wpts,1,sum))
      g.vals = fit$setup$gfun(wpts,par=fit$mle)
      
      if(add){
        surface3d(wpts[,1] / matrix(g.vals,100,100),
                  wpts[,2] / matrix(g.vals,100,100),
                  wpts[,3] / matrix(g.vals,100,100),
                  col=col,alpha=alphaval,add=add)
      } else{
      plot3d(NA,xlab="x1",ylab="x2",zlab="x3")
      surface3d(wpts[,1] / matrix(g.vals,100,100),
                wpts[,2] / matrix(g.vals,100,100),
                wpts[,3] / matrix(g.vals,100,100),
                col=col,alpha=alphaval)}
      } else if(fit$class=="parametric")
      {
        gfun.plot <- function(x, y, z) {
          fit$setup$gfun(c(x, y, z), par = fit$mle)
        }
        gfun.plot <- Vectorize(gfun.plot)
        x1 <- seq(0, 1, len = resolution)
        if (add) {
          misc3d::contour3d(f = gfun.plot, level = c(1), x = x1,
                            y = x1, z = x1, color = col, alpha = alphaval, add = T)
        }
        else {
          require(rgl)
          misc3d::contour3d(f = gfun.plot, level = c(1), x = x1,
                            y = x1, z = x1, color = col, alpha = alphaval)
          rgl::axes3d()
          rgl::title3d(xlab = xlab, ylab = ylab, zlab = zlab)
      }
  }
    }
}


#' Plot any limit set (unit level set of gauge function) for d=2,3

#' @param gfun function specifying the gauge to plot. Should have arguments x and par
#' @param par vector giving the parameters of the specified gauge function
#' @param d dimension (should be 2 or 3)
#' @param resolution plotting resolution - the higher the number, the finer the resolution of the plot. For d=2, the default should suffice, but larger numbers should not take too much longer to plot. For d=3, higher numbers will normally take longer to plot, but give better figures.
#' @param add logical - add plot to an existing figure?
#' @param col character string or number giving the colour for the plot
#' @param xlab,ylab,zlab character string giving axis labels
#' @param lwd,type line width and type arguments to plot for d=2
#' @param polygon logical - include shaded polygon of limit set for d=2?
#' @param polygoncol character string or number giving the colour for the polygon
#' @param alphaval alpha value for shading of the polygon or 3d level set
#' @param ... additional arguments to pass to plot for d=2
#' 
#' @return plot of fitted gauge function
#'
#' @export
#' 
#' 

plotanygauge<-function(gfun,par,d,resolution=70,add=FALSE,col="red",xlab="x",ylab="y",zlab="z",lwd=2,type="l",polygon=TRUE,polygoncol=rgb(0.5,0,0,alpha=0.5),alphaval=0.5,...)
{
  if(d==2)
  {
    wpts<-seq(0,1,len=10*resolution+1)
  if (add) {
    lines(wpts/apply(cbind(wpts, 1 - wpts), 1, gfun, par=par),
          (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun, par=par),
          col = col, typ = type, lwd = lwd,...)
    if (polygon) {
      polygon(c(0, wpts/apply(cbind(wpts, 1 - wpts), 1,
                              gfun, par=par), 0), c(0, (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun, par=par), 0), col = polygoncol, border = NULL)
    }
  } else {
    plot(wpts/apply(cbind(wpts, 1 - wpts), 1, gfun, par=par),
         (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun, par=par),
         typ = type, col = col, xlab = xlab,
         ylab = ylab, lwd = lwd,...)
    if (polygon) {
      polygon(c(0, wpts/apply(cbind(wpts, 1 - wpts), 1,
                              gfun, par=par), 0), c(0, (1 - wpts)/apply(cbind(wpts, 1 - wpts), 1, gfun, par=par), 0), col = polygoncol, border = NULL)
    }
  }
  } else if(d==3)
  {
    gfun.plot <- function(x, y, z) {
      gfun(c(x, y, z), par = par)
    }
    gfun.plot <- Vectorize(gfun.plot)
    x1 <- seq(0, 1, len = resolution)
    if (add) {
      misc3d::contour3d(f = gfun.plot, level = c(1), x = x1,
                        y = x1, z = x1, color = col, alpha = alphaval, add = T)
    }
    else {
      require(rgl)
      misc3d::contour3d(f = gfun.plot, level = c(1), x = x1,
                        y = x1, z = x1, color = col, alpha = alphaval)
      rgl::axes3d()
      rgl::title3d(xlab = xlab, ylab = ylab, zlab = zlab)
    }
  }
}

##############################################################################

#' Create chi(u) plots

#' @param thresh.fit output from fit.thresh
#' @param fit list output from fit.geometric.par or fit.geometric.pwlin
#' @param dataset dataset x in exponential margins

#' @return plots of chi(u)
#'
#' @export

tail.dep.plots = function(thresh.fit,fit,dataset){
  
  model.fit = fit
  # get the combinations of variable indices
  d = model.fit$dimension
  idx.groups = lapply(c(1:d), function(x) combn(d,x))
  idx.groups = lapply(idx.groups, function(entry) lapply(c(1:ncol(entry)), function(i) entry[,i]))
  idx.groups = unlist(idx.groups, recursive = F)
  idx.groups = idx.groups[sapply(idx.groups,function(vec) length(vec)>1)]
  
  n.mesh=100
  wpts = expand.grid(replicate(d-1,seq(0,1,len=n.mesh),simplify=F))
  wpts = cbind(wpts,1-apply(wpts,1,sum))
  cond = wpts[,d]>=0  # only consider angles in the positive orthant.
  wpts = wpts[cond,]
  thresh.r = eval.thresh(thresh.fit,wpts)
  u.vals.lst = lapply(idx.groups,function(idx){
    u.min = pexp(max(sapply(1:length(thresh.r),function(w.idx){
      w.vec = wpts[w.idx,]
      r0w.val = thresh.r[w.idx]
      return(r0w.val*min(w.vec[idx],na.rm=T))
    }),na.rm=T))
    return(seq(u.min,1,length.out=100))
  })
  
  xstar = sim.geometric(model.fit,nsim=50000)
  pr.exc = length(model.fit$data$r) / dim(dataset)[1]
  
  chi.vals = list(list(NULL))
  grp.no=1
  for(which.cols in idx.groups){
    u.vals = u.vals.lst[[grp.no]]
    
    # empirical
    chi.u.emp = sapply(u.vals, function(u) mean(apply(dataset[,which.cols],1,function(xx) all(pexp(xx)>u)))/(1-u))  # assuming the margins are approximately standard exponential
    
    # estimates
    chi.u.est = pr.exc*sapply(u.vals, function(u) mean(apply(xstar[,which.cols],1,function(xx){all(pexp(xx)>u)}))/(1-u))
    
    col.names = as.character(c(1,2,3))
    col.names = col.names[which.cols]
    name = paste(col.names, collapse = "")
    
    chi.vals[[grp.no]] = list(u.vals=u.vals,
                              chi.u.emp=chi.u.emp,
                              chi.u.est=chi.u.est,
                              name=name)
    
    grp.no=grp.no+1
  }
  
  if(d>1){par(mfrow=rep(ceiling(sqrt(length(idx.groups))),2),pty="m")}
  col.names = as.character(c(1:3))
  grp.no = 1
  for(which.cols in idx.groups){
    
    u.vals = chi.vals[[grp.no]]$u.vals
    chi.u.emp = chi.vals[[grp.no]]$chi.u.emp
    chi.u.est = chi.vals[[grp.no]]$chi.u.est
    name = chi.vals[[grp.no]]$name
    ylab = bquote(.(rlang::sym("chi"))[.(name)](u))
    plot(u.vals,chi.u.emp,type="l",lwd=1,
         xlab="u",ylab=ylab,
         ylim=c(0,1.0))
    lines(u.vals,chi.u.est,col="blue",lty=2)
    grp.no=grp.no+1
  }
}
