#' Logistic gauge function
#' @param vec values at which to evaluate the gauge function
#' @param par parameter value
#' @noRd

#' @export

gauge_logistic<-function(vec,par)
{
  sum(vec)/par+min(vec)*(1-length(vec)/par)
}

#' Asymmetric logistic gauge function
#' @param vec values at which to evaluate the gauge function (of length 2 or 3)
#' @param par parameter value
#' @param theta1 set equal to 0 if do not want extremes of this group of variables
#' @param theta2 set equal to 0 if do not want extremes of this group of variables
#' @param theta3 set equal to 0 if do not want extremes of this group of variables
#' @param theta12 set equal to 0 if do not want extremes of this group of variables
#' @param theta23 set equal to 0 if do not want extremes of this group of variables
#' @param theta13 set equal to 0 if do not want extremes of this group of variables
#' @param theta123 set equal to 0 if do not want extremes of this group of variables

#' @noRd

#' @export

gauge_alogistic<-function(vec, par, theta1 = 1, theta2 = 1, theta3 = 1, theta12 = 1, 
          theta13 = 1, theta23 = 1, theta123 = 1) 
{
  if(length(vec)==2){
    x <- vec[1]; y <- vec[2]
    mxy <- min(x,y)
    k <- 1
    if (theta12 == 1) {
      alpha12 <- par[k]
      k <- k + 1
    }
    else {
      alpha12 <- 1
    }
    V1 <- c(x[theta1 == 1], (x/alpha12 + mxy * (1 - 1/alpha12))[theta12 == 1])
    V2 <- c(y[theta2 == 1], (y/alpha12 + mxy * (1 - 1/alpha12))[theta12 == 1])
    V12 <- ((x + y)/alpha12 + mxy * (1 - 2/alpha12))[theta12 ==1]
    suppressWarnings(term1 <- min(apply(expand.grid(V1, V2, 
                                                    KEEP.OUT.ATTRS = F), 1, function(x) {
                                                      sum(x)
                                                    })))
    term2 <- V12
    return(min(term1, term2)) 
  }
  if(length(vec)==3){
    x <- vec[1]; y <- vec[2]; z <- vec[3]
    mxy <- min(x, y)
    mxz <- min(x, z)
    myz <- min(y, z)
    mxyz <- min(x, y, z)
    k <- 1
    if (theta12 == 1) {
      alpha12 <- par[k]
      k <- k + 1
    }
    else {
      alpha12 <- 1
    }
    if (theta13 == 1) {
      alpha13 <- par[k]
      k <- k + 1
    }
    else {
      alpha13 <- 1
    }
    if (theta23 == 1) {
      alpha23 <- par[k]
      k <- k + 1
    }
    else {
      alpha23 <- 1
    }
    if (theta123 == 1) {
      alpha123 <- par[k]
      k <- k + 1
    }
    else {
      alpha123 <- 1
    }
    V1 <- c(x[theta1 == 1], (x/alpha12 + mxy * (1 - 1/alpha12))[theta12 == 
                                                                  1], (x/alpha13 + mxz * (1 - 1/alpha13))[theta13 == 1], 
            (x/alpha123 + mxyz * (1 - 1/alpha123))[theta123 == 1])
    V2 <- c(y[theta2 == 1], (y/alpha12 + mxy * (1 - 1/alpha12))[theta12 == 
                                                                  1], (y/alpha23 + myz * (1 - 1/alpha23))[theta23 == 1], 
            (y/alpha123 + mxyz * (1 - 1/alpha123))[theta123 == 1])
    V3 <- c(z[theta3 == 1], (z/alpha23 + myz * (1 - 1/alpha23))[theta23 == 
                                                                  1], (z/alpha13 + mxz * (1 - 1/alpha13))[theta13 == 1], 
            (z/alpha123 + mxyz * (1 - 1/alpha123))[theta123 == 1])
    V12 <- c(((x + y)/alpha12 + mxy * (1 - 2/alpha12))[theta12 == 
                                                         1], ((x + y)/alpha123 + mxyz * (1 - 2/alpha123))[theta123 == 
                                                                                                            1])
    V13 <- c(((x + z)/alpha13 + mxz * (1 - 2/alpha13))[theta13 == 
                                                         1], ((x + z)/alpha123 + mxyz * (1 - 2/alpha123))[theta123 == 
                                                                                                            1])
    V23 <- c(((y + z)/alpha23 + myz * (1 - 2/alpha23))[theta23 == 
                                                         1], ((y + z)/alpha123 + mxyz * (1 - 2/alpha123))[theta123 == 
                                                                                                            1])
    V123 <- ((x + y + z)/alpha123 + mxyz * (1 - 3/alpha123))[theta123 == 
                                                               1]
    suppressWarnings(term1 <- min(apply(expand.grid(V1, V2, V3, 
                                                    KEEP.OUT.ATTRS = F), 1, function(x) {
                                                      sum(x)
                                                    })))
    suppressWarnings(term2 <- min(apply(expand.grid(V1, V23, 
                                                    KEEP.OUT.ATTRS = F), 1, function(x) {
                                                      sum(x)
                                                    })))
    suppressWarnings(term3 <- min(apply(expand.grid(V2, V13, 
                                                    KEEP.OUT.ATTRS = F), 1, function(x) {
                                                      sum(x)
                                                    })))
    suppressWarnings(term4 <- min(apply(expand.grid(V3, V12, 
                                                    KEEP.OUT.ATTRS = F), 1, function(x) {
                                                      sum(x)
                                                    })))
    term5 <- V123
    return(min(term1, term2, term3, term4, term5)) 
  }
}

#' @export

gauge_invlogistic <- function(vec, par)
{
  return(sum(vec^(1/par))^(par))
}

#' Inverted Clayton gauge function
#' @param vec values at which to evaluate the gauge function
#' @param par parameter value
#' @noRd

#' @export

gauge_invclayton<-function(vec,par)
{
  return(max(vec)*(1+length(vec)*par)-sum(vec)*par)
}

#' Gaussian gauge function
#' @param vec values at which to evaluate the gauge function (d=2,3,4)
#' @param par parameter value (correlation parameters filling up matrix)
#' @noRd

#' @export

gauge_gaussian<-function(vec,par)
{
  d=length(vec)
  if(d==2){
    S<-matrix(c(1,par,par[1],1),2,2,byrow=T)
  } else if(d==3){
    S<-matrix(c(1,par[1:2],par[1],1,par[3],par[2],par[3],1),3,3,byrow=T)
  } else if(d==4){
    S <- matrix(c(1, par[1:3],
                  par[1], 1, par[4:5],
                  par[2], par[4], 1, par[6],
                  par[3], par[5], par[6], 1), 4, 4, byrow = T)
  } else{stop("Dimension must be between 2 and 4")}
  return(as.numeric(t(sqrt(vec))%*%solve(S)%*%sqrt(vec)))
}

#' Square gauge function
#' @param vec values at which to evaluate the gauge function (d=2)
#' @param par parameter value
#' @noRd

#' @export

gauge_square<-function(vec,par)
{
  return(max((vec[1]-vec[2])/par,(vec[2]-vec[1])/par,(vec[1]+vec[2])/(2-par)))
}

#' Exponential-Gaussian additive gauge function
#' @param vec values at which to evaluate the gauge function (d=2)
#' @param par parameter value
#' @noRd

#' @export

gauge_expgauss <- function(vec, par){
  x <- vec[1];y <- vec[2]
  Gam=par[1];Rho=par[2]
  
  if(Gam <= (1+Rho)/2 | grepl("No real roots.",quadSol(Gam,Rho))){
    out=(x + y - 2 * Rho * sqrt(x * y))/(1 - Rho^2)
  }else{
    x_sol=quadSol(Gam,Rho)
    y_sol=2*Rho*sqrt((1-Rho^2)*(x_sol-x_sol^2))-(x_sol*(1-2*Rho^2)+Rho^2-1)
    wstar=(x_sol)/(x_sol+y_sol)
    m=(y_sol-Gam)/(x_sol-Gam) # the slope of the tangent line
    if(m!=0){ # equivalently, (1+Rho)/2 < gamma < 1
      if(x <= y){
        out=ifelse(x/(x+y) <= wstar,
                   (x + y - 2 * Rho * sqrt(x * y))/(1 - Rho^2),
                   (y-m*x)/(Gam*(1-m)))
      }else{ 
        out=ifelse(x/(x+y) <= 1-wstar,
                   (y-(1/m)*x)/(Gam*(1-(1/m))),
                   (x + y - 2 * Rho * sqrt(x * y))/(1 - Rho^2))
      } 
    }else{ # equivalently, gamma = 1 (slope = 0)
      if(x <= y){
        out=ifelse(x/(x+y) <= wstar,
                   (x + y - 2 * Rho * sqrt(x * y))/(1 - Rho^2),
                   y)
      }else{
        out=ifelse(x/(x+y) <= 1-wstar,
                   x,
                   (x + y - 2 * Rho * sqrt(x * y))/(1 - Rho^2))
      }
    }
    
    if(Gam > 1){out<-out*Gam}
  }
  return(out)
}

#' Exponential-Inverted Logistic additive gauge function
#' @param vec values at which to evaluate the gauge function (d=2)
#' @param par parameter value
#' @noRd

#' @export

gauge_expinvlog <- function(vec, par){
  x <- vec[1];y <- vec[2]
  Gam=par[1];Theta=par[2]
  
  if(Gam <= 2^(-Theta)){
    out=(x^(1/Theta)+y^(1/Theta))^Theta
  }else{
    Sol=find_solutions(gamma = Gam,theta = Theta) # need to call an outer ft
    x_sol=min(Sol)
    y_sol=(1-x_sol^(1/Theta))^Theta
    wstar=(x_sol)/(x_sol+y_sol)
    m=(y_sol-Gam)/(x_sol-Gam) # the slope of the tangent line
    if(x<=y){
      out=ifelse(x/(x+y) <= wstar,
                 (x^(1/Theta)+y^(1/Theta))^Theta,
                 (y-m*x)/(Gam*(1-m)))
    }else{
      out=ifelse(x/(x+y) <= 1-wstar,
                 (y-(1/m)*x)/(Gam*(1-(1/m))),
                 (x^(1/Theta)+y^(1/Theta))^Theta)
    }
    
    if(Gam==1){
      return(max(x,y))
    }
    
    if(Gam > 1){
      out=ifelse(x <= y, y-((Gam-1)/Gam)*x, x-((Gam-1)/Gam)*y)
      out<-out*Gam
    }
  }
  return(out)
}

#' Exponential-square additive gauge function
#' @param vec values at which to evaluate the gauge function (d=2)
#' @param par parameter value
#' @noRd

#' @export

gauge_expsquare <- function(vec, par){
  x <- vec[1];  y <- vec[2]
  
  Gam=par[1];Theta=par[2]
  if(Gam <= 1-Theta/2){
    out=max((x - y)/Theta, (y - x)/Theta, (x + y)/(2 - Theta))
  }else{
    out=ifelse(x<=y,
               max((x-y)/Theta,(y-x)/Theta,((Theta-1)/(Gam*Theta)+1/Theta)*y+(1/(Gam*Theta)-1/Theta)*x),
               max((x-y)/Theta,(y-x)/Theta,((Theta-1)/(Gam*Theta)+1/Theta)*x+(1/(Gam*Theta)-1/Theta)*y))
    if(Gam > 1){out<-out*Gam}
  }
  return(out)
}

#' Logistic-Gaussian additive gauge function
#' @param vec values at which to evaluate the gauge function (d=2)
#' @param par parameter value
#' @noRd

#' @export

gauge_addloggauss <- function(vec, par){
  x <- vec[1]; y <- vec[2]
  
  Rho=par[1];Gam=par[2];Wgt=par[3]
  
  ab=((1-Wgt)/Wgt)*(1/Gam-1)
  C=(1-Rho^2)*ab
  
  ##  Slope 'kappa'
  kp=(Rho^2)/((1-(1-Rho^2)*(1-Wgt)*(1/Gam-1)/Wgt)^2)
  
  if(Rho==0){
    if(Gam > 1-Wgt){
      ms=(1/Gam+Wgt*((1-Rho^2)^(-1)-1/Gam))^(-1)
    }
    if(Gam == 1-Wgt){
      ms=1/(2-Gam)
    }
    if(Gam < 1-Wgt){
      ms=1/(1+Wgt)
    }
  }
  
  if(C >= 1 & 0 < Rho & Rho < 1){ # k'(z) < 0
    ms=1/(2*Wgt/(1+Rho)+1-Wgt)
  }
  
  if((1+Rho)*ab < 1 & 0 < Rho & Rho < 1){ # kappa \in (0,1)
    num=1/(kp+1)
    den=(Wgt/(1-Rho^2))*(1-2*Rho*sqrt(kp)/(kp+1))+(1-Wgt)*(1/Gam+(1-2/Gam)*(kp/(kp+1)))
    ms=num/den 
  }
  if(1<= (1+Rho)*ab & C < 1 & 0 < Rho & Rho < 1){ # kappa >= 1
    ms=1/(2*Wgt/(1+Rho)+1-Wgt)
  }
  
  mixture.GaussLog=Wgt*(x + y - 2 * Rho * sqrt(x * y))/(1 - Rho^2) + (1-Wgt)*(max(x, y)/Gam + min(x, y) * (1 - 1/Gam))
  return(ms*mixture.GaussLog)
}

#' Logistic-inverted logistic additive gauge function
#' @param vec values at which to evaluate the gauge function (d=2)
#' @param par parameter value
#' @noRd

#' @export

gauge_addloginvlog <- function(vec, par){
  x <- vec[1]; y <- vec[2]
  
  Theta=par[1];Gam=par[2];Wgt=par[3]
  
  ##  Slope 'kappa'
  a=(1-Wgt)/Wgt
  b=(1/Gam-1)
  kp=((a*b)^(1/(Theta-1))-1)^(-Theta)
  
  if(Theta==1){
    if(Gam > 1-Wgt){
      ms=1/(Wgt+(1-Wgt)/Gam)
    }
    if(Gam == 1-Wgt){
      ms=1/(2-Gam)
    }
    if(Gam < 1-Wgt){
      ms=1/(1+Wgt)
    }
  }
  
  if(a*b < 2^(Theta-1) & Theta < 1){
    num=1/(kp+1)
    den=Wgt*((1+kp)^(-1/Theta)+(kp/(kp+1))^(1/Theta))^(Theta)+(1-Wgt)*(1/Gam+(1-2/Gam)*kp/(kp+1))
    ms=num/den
  }
  
  if(a*b >= 2^(Theta-1) & Theta < 1){
    ms=1/(Wgt*(2^Theta)+(1-Wgt)) 
  }
  
  mixture.InvLogLog=Wgt*((x^(1/Theta)+y^(1/Theta))^Theta) + (1-Wgt)*(max(x, y)/Gam + min(x, y)*(1-1/Gam))
  return(ms*mixture.InvLogLog)
}

#' Logistic-square additive gauge function
#' @param vec values at which to evaluate the gauge function (d=2)
#' @param par parameter value
#' @noRd

#' @export

gauge_addlogsquare <- function(vec, par){
  x <- vec[1]; y <- vec[2]
  
  Theta=par[1];Gam=par[2];Wgt=par[3]
  ab=((1-Wgt)/Wgt)*(1-Gam)/Gam
  
  if(Theta==1){
    if(Gam >= 1-Wgt){
      ms=1/(Wgt+(1-Wgt)/Gam)  
    }
    if(Gam < 1-Wgt){
      ms=1/(1+Wgt)
    }
  }
  
  if(ab < (2-Theta)^(-1) & Theta < 1){
    ms=1/(Wgt+(1-Wgt)*(1-Theta+Theta/Gam))
  }
  
  if(ab >= (2-Theta)^(-1) & Theta < 1){
    ms=1/(2*Wgt/(2-Theta)+1-Wgt)
  }
  
  mixture.RectLog=Wgt*(max((x - y)/Theta, (y - x)/Theta, (x + y)/(2 - Theta))) + (1-Wgt)*(max(x, y)/Gam + min(x, y) * (1 - 1/Gam))
  return(ms*mixture.RectLog)
}



