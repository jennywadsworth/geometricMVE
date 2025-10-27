quadSol <- function(gam,rho) {
  
  a=(2*gam-1)^2
  b=(2*(rho^2)*gam*(1-2*gam)+(rho^2)*(2*gam-1)^2-(2*gam-1)^2)
  c=(rho*gam)^2
  
  ##  discriminant
  disc <- (b^2) - (4*a*c)
  
  if(disc < 0) {  #no real roots
    return(paste0("No real roots."))
  }
  else if(disc > 0) { 
    x_solp <- (-b + sqrt(disc)) / (2*a)
    x_solm <- (-b - sqrt(disc)) / (2*a)
    x_sol <- min(x_solp,x_solm) # x* <= (1+rho)/2
    return(x_sol)
  }
  else{
    x_sol <- (-b) / (2*a)
    return(x_sol)
  }
}
