
#############################################################################
################################ d=2 ########################################
#############################################################################

# given 3 parameters (par) located at 3 reference angles (par.locs),
# return the gauge function value at the point xyz
gfun.simple.d2pwlin.L1 = function(w,par,par.locs){
  # xyz      = 
  # par      = vector of length 2
  # par.locs = vector of length 2
  
  r1 = par[1]
  r2 = par[2]
  
  w.sort = sort(par.locs)
  w.low = w.sort[1]
  w.up = w.sort[2]
  
  dx = r2*w.up - r1*w.low
  dy = r2*(1-w.up) - r1*(1-w.low)
  a = dy/dx
  b = r1*(1-w.low - a*w.low)
  if(dx==0){
    return(w/(r1*w.low))
  } 
  else{
    return((1-w-(a*w))/b)
  }
}


pwlin.g.vals.2d = function(w.adj.angles,par,par.locs){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.
  
  sapply(w.adj.angles, function(lst){
    which.angles = lst$idx.locs
    return(gfun.simple.d2pwlin.L1(w=lst$w,par=par[which.angles],par.locs=par.locs[which.angles]))
  })
}

gfun.2d = function(x, par, ref.angles){
  # x-> matrix of 2 columns
  if(is.null(dim(x)) & length(x)==1){
    x = cbind(x,1-x)
  }
  if(is.null(dim(x)) & length(x)>1){
    x = matrix(x,nrow=1)
  }
  
  rad.inp = apply(x,1,sum)
  angle.inp = x[,1] / rad.inp
  
  w.adj.angles = which.adj.angles.2d(angles=angle.inp, locs=ref.angles)
  
  gval = rad.inp*pwlin.g.vals.2d(w.adj.angles=w.adj.angles,par=par,par.locs=ref.angles)
  
  return(gval)
}

#############################################################################
################################ d>2 ########################################
#############################################################################

gfun.simple.pwlin = function(xyz,par,par.locs){
  # par      = in R^3_+
  # par.locs = angles (in the simplex) where the parameters are location
  
  num.cols = dim(par.locs)[2]

  coplanar.mat = do.call(rbind,lapply(c(1:num.cols)[-1],function(i){
    (par[1]*par.locs[1,])-(par[i]*par.locs[i,])
  }))
  norm.vec = suppressWarnings({
    c(1,-1)*sapply(c(1:num.cols),function(i){det(coplanar.mat[,-i])})
  })
  sum(norm.vec * xyz) / sum(norm.vec * par[1] * par.locs[1,])
}

# evaluating the gauge at a set of angles, returns a vector
pwlin.g.vals = function(w.adj.angles,par,par.locs){
  # w.adj.angles  -> output of which.adj.angles()
  # par           -> parameters at reference angles
  # par.locs      -> 3-column matrix of reference angles.
  
  sapply(w.adj.angles, function(lst){
    which.angles = lst$loc.idx
    if(any(is.na(lst$w)) | all(which.angles==0) | any(lst$w<0) | any(is.na(which.angles))){  # if angle is not in the positive orthant OR if angle lies on line
      return(NA)
    } else {
      return(gfun.simple.pwlin(xyz=lst$w,par=par[which.angles],par.locs=par.locs[which.angles,]))
    }
  })
}

# now evaluate on Cartesian points
gfun.pwl = function(x, par, ref.angles){
  # x-> matrix of d columns
  if(is.null(dim(x))){
    x = matrix(x,nrow=1)
  }
  
  num.cols = dim(x)[2]
  
  rad.inp = apply(x,1,sum)
  angle.inp = x / rad.inp
  is.valid = apply(angle.inp,1,function(vec) all(!is.na(vec)))
  w.adj.angles = replicate(nrow(angle.inp), 
                           list(w=rep(NA,num.cols),  
                                loc.idx=rep(NA,num.cols),
                                vertices=matrix(NA,num.cols,num.cols)), FALSE)
  inn = as.matrix(angle.inp[is.valid,],ncol=num.cols,byrow=T)
  if(ncol(inn)==1){
    inn = t(inn)
  }
  w.adj.angles[is.valid] = which.adj.angles(angles=inn, locs=ref.angles)
  
  gval = rad.inp*pwlin.g.vals(w.adj.angles=w.adj.angles,par=par,par.locs=ref.angles)
  return(gval)
}

############################################################################

## Parametric Gauges in Laplace margins for general d-dimensions (not in geometricMVE)
gauge_gaussian_lap = function(vec,Sigma){
  Q = solve(Sigma)
  dd = length(vec)
  sgn = sign(vec)
  vec.mat = matrix(vec,ncol=1)
  sgn.mat = matrix(sgn,ncol=1)
  t(sgn.mat * sqrt(abs(vec.mat))) %*% Q %*% (sgn.mat * sqrt(abs(vec.mat)))
}
gauge_logistic_lap = function(vec,dep.par){
  dd = length(vec)
  is.pos = vec>0
  is.neg = vec<0
  is.zero = vec==0
  if(sum(is.pos)==dd){
    return((1/dep.par)*sum(vec) + (1-(dd/dep.par))*min(vec))
  } else if(sum(is.neg)==dd){
    return((sum((-vec)^(1/dep.par)))^dep.par)
  } else {
    return((1/dep.par)*sum(vec[is.pos]) + (sum((-vec[is.neg])^(1/dep.par)))^dep.par)
  }
}

##########################################################################

gauss.alog.mixgauge = function(xyz,par){
  require(geometricMVE)
  
  # gauge for d=3 simulation study
  # par -> length 5
  if(length(par) != 5){
    stop("parameters not of correct length")
  }
  min(geometricMVE::gauge_gaussian3d(xyz,par=par[1:3]),
      geometricMVE::gauge_rvad_all3d(xyz,par=par[4:5],
                                     theta1=0,theta2=0,theta3=0,theta13=0,theta23=0))
}

gauss.mixgauge = function(xyz,par){
  require(geometricMVE)
  
  # gauge for d=3 simulation study
  # par -> length 5
  if(length(par) != 6){
    stop("parameters not of correct length")
  }
  min(geometricMVE::gauge_gaussian3d(xyz,par=par[1:3]),
      geometricMVE::gauge_gaussian3d(xyz,par=par[4:6]))
}
