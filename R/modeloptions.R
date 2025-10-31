#' Details model options for parametric fitting (to be used internally)

#' @param model Character string specifying which gauge model to use. If two names are provided, these are additively mixed with a weight parameter.
#' @param customgauge Optional custom gauge function
#' @param init.val initial values of custom gauge function parameters (must be supplied if customgauge used)
#' @param lower.limit lower bounds of custom gauge function parameters (must be supplied if customgauge used)
#' @param upper.limit upper bounds of custom gauge function parameters (must be supplied if customgauge used)
#' @param d dimension of model fit
#' @param theta12 used if moel = "alogistic" to specify whether there are joint extremes for this group
#' @param theta13 used if moel = "alogistic" to specify whether there are joint extremes for this group
#' @param theta23 used if moel = "alogistic" to specify whether there are joint extremes for this group
#' @param theta123 used if moel = "alogistic" to specify whether there are joint extremes for this group
#' @noRd
modeloptions<-function(model,customgauge,init.val,lower.limit,
                       upper.limit,d,theta12,theta13,theta23,theta123)
{
  if(!is.null(customgauge))
  {
    gfun=customgauge
    if(is.null(lower.limit)||is.null(upper.limit)){stop("Please provide initial value, upper and lower parameter bounds for the custom gauge function.")}
    init.val=init.val
    lower.limit=lower.limit
    upper.limit=upper.limit
    add.gauge=FALSE
    gauge1=gauge2=par1ind=par2ind=NULL
  } else{
    if(length(model)>2){stop("Cannot mix more than two gauges currently. The user may wish to explore use of the customgauge argument.")}
    if(d>2 & (is.element("expgauss",model)|is.element("expinvlog",model)|is.element("expsquare",model)|is.element("square",model))){stop("Selected gauges only implemented for d=2.")}
    if(d>4 & is.element("gauss",model)){stop("Gaussian gauge currently implemented in dimension up to 4. Use customgauge to implement in higher dimensions.")}

    functionvector<-c(gauge_logistic,gauge_invlogistic,gauge_gaussian,gauge_square,gauge_invclayton,gauge_expgauss,gauge_expinvlog,gauge_expsquare,gauge_alogistic)
    labelvector<-c("log","invlog","gauss","square","invclayton","expgauss", "expinvlog","expsquare","alogistic")
    if(is.element("alogistic",model)){aloglength<-theta12+(theta23+theta13+theta123)*(d==3)} else{aloglength<-1}
    parlengthvector<-c(1,1,d*(d-1)/2,1,1,2,2,2,aloglength)
    lowerlimlist<-list(c(0),c(0),rep(0,d*(d-1)/2),c(0),c(0),c(0,0),c(0,0),c(0,0),c(rep(0,aloglength)))
    upperlimlist<-list(c(1),c(1),rep(0.999,d*(d-1)/2),c(1),c(100),c(10,1),c(10,1),c(10,1),c(rep(1,aloglength)))
    initvallist<-list(c(0.5),c(0.5),rep(0.5,d*(d-1)/2),c(0.5),c(0.8),c(0.8,0.5),c(0.8,0.5),c(0.8,0.5),c(rep(0.5,aloglength)))

    if(length(model)==1)
    {
      gfun<-functionvector[[which(labelvector==model)]]
      init.val=initvallist[[which(labelvector==model)]]
      lower.limit=lowerlimlist[[which(labelvector==model)]]
      upper.limit=upperlimlist[[which(labelvector==model)]]
      add.gauge=FALSE
      gauge1=gauge2=par1ind=par2ind=NULL
    }
    else if(length(model)==2){
      if(is.element("log",model)&is.element("invlog",model)&d==2)
      {
        gfun=gauge_addloginvlog
          init.val=c(0.5,0.5,0.5)
          lower.limit=c(0,0,0)
          upper.limit=c(1,1,1)
          add.gauge=FALSE
        gauge1=gauge2=par1ind=par2ind=NULL
      } else if(is.element("log",model)&is.element("gauss",model)&d==2)
      {
        gfun=gauge_addloggauss
        init.val=c(0.5,0.5,0.5)
        lower.limit=c(0,0,0)
        upper.limit=c(1,1,1)
          add.gauge=FALSE
        gauge1=gauge2=par1ind=par2ind=NULL
      } else if(is.element("log",model)&is.element("square",model)&d==2)
      {
        gfun=gauge_addlogsquare
        init.val=c(0.5,0.5,0.5)
        lower.limit=c(0,0,0)
        upper.limit=c(1,1,1)
          add.gauge=FALSE
        gauge1=gauge2=par1ind=par2ind=NULL
      } else{
        add.gauge=TRUE
        gauge1=functionvector[[which(labelvector==model[1])]]
        gauge2=functionvector[[which(labelvector==model[2])]]
        par1ind= 1:length(initvallist[[which(labelvector==model[1])]])
        par2ind= (length(initvallist[[which(labelvector==model[1])]])+1):(length(initvallist[[which(labelvector==model[1])]])+length(initvallist[[which(labelvector==model[2])]]))
        init.val=c(initvallist[[which(labelvector==model[1])]],initvallist[[which(labelvector==model[2])]],0.5)
        lower.limit=c(lowerlimlist[[which(labelvector==model[1])]],lowerlimlist[[which(labelvector==model[2])]],0)
        upper.limit=c(upperlimlist[[which(labelvector==model[1])]],upperlimlist[[which(labelvector==model[2])]],1)
        gfun=NULL
      }
    }
  }
  return(list(gfun=gfun,init.val=init.val,lower.limit=lower.limit,
              upper.limit=upper.limit,add.gauge=add.gauge,
              gauge1=gauge1,gauge2=gauge2,par1ind=par1ind,par2ind=par2ind))
}





