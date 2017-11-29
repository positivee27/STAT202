#' Perform ridge regression of Y on X.
#' 
#' @param matrix n x p matrix of explanatory variables
#' @param matrix n dimensional response vector
#' @param vector regularization parameters
#' @return Returns beta, the ridge regression solution.
#' @export

myRidge <- function(X, Y, lambda){
  n = dim(X)[1]
  p = dim(X)[2]
  z = cbind(rep(1, n), X)
  D = diag(rep(sqrt(lambda), p))
  D = cbind(matrix(rep(0, p), nrow = p), D)
  A = rbind(z, D)
  y = matrix(c(Y, rep(0, p)), ncol = 1)
  A = cbind(A, y)
  qr_result = myQR(A)
  R = qr_result$R
  R1 = R[1:(p+1), 1:(p+1)]
  Y1 = R[1:(p+1), (p+2)]
  beta_ridge = solve(R1, Y1)
  
  return(beta_ridge)
}