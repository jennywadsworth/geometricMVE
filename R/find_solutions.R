#' Find the real solutions for x
#'
#' This function solves the equation(1-x^(1/theta))^theta - gamma = -((x^(1/theta-1))*(1-x^(1/theta))^(theta-1))*(x-gamma)
#' for x, given the specified constraints on the parameters.
#'
#' @param gamma A numeric value in the interval (2^(-theta), 1].
#' @param theta A numeric value in the interval (0, 1].
#' @param tol The desired numerical tolerance for the solver. Default is 1e-9.
#'
#' @return A numeric vector containing the two solutions for x. 
#' 
#' 
#'
#' @noRd
#' find_solutions(gamma = 0.9, theta = 0.5)
find_solutions <- function(gamma, theta, tol = 1e-9){
  
  #if (theta <= 0 || theta > 1) {
  #  stop("Error: theta must be in the interval (0, 1].")
  #}
  
  # Case: gamma = 1
  # The solutions are known to be 0 and 1.
  if (gamma == 1) {
    return(c(0, 1))
  }
  
  # We need to find the root of h(u) = 0, where u = x^(1/theta)
  # h(u) = (1-u)^(1-theta) + u^(1-theta) - 1/gamma
  p <- 1 - theta
  h <- function(u) {
    # The term u^p can be problematic at u=0 if p is not an integer.
    # We can handle u=0 as a special case.
    if (u == 0) {
      return(1 - 1/gamma)
    }
    return((1 - u)^p + u^p - 1/gamma)
  }
  
  # We find the first root, u1, in the interval [0, 0.5].
  # The other root, u2, will be 1 - u1.
  a <- 0.0
  b <- 0.5
  
  # Bisection loop
  while ((b - a) > tol) {
    mid <- (a + b) / 2
    if (h(mid) == 0) { # Found exact root
      a <- mid
      b <- mid
      break
    }
    # Use sign() to avoid potential floating point issues with multiplication
    if (sign(h(a)) != sign(h(mid))) {
      b <- mid
    } else {
      a <- mid
    }
  }
  
  u1 <- (a + b) / 2
  u2 <- 1 - u1
  
  # Convert u back to x
  x1 <- u1^theta
  x2 <- u2^theta
  
  return(sort(c(x1, x2)))
}