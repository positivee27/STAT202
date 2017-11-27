#' Find the lasso solution path for various values of 
#' the regularization parameter lambda.
#'
#' @param matrix n x p matrix of explanatory variables
#' @param matrix n dimensional response vector
#' @param vector Vector of regularization parameters
#' @return a matrix containing the lasso solution vector 
#'   beta for each regularization parameter.
#' @export

myLasso <- function(X, Y, lambda_all){
  n = nrow(X)
  p = ncol(X)
  L = length(lambda_all)
  T = 100
  R = Y
  beta = c(rep(0, p+1))
  beta_all = matrix(rep(0, (p+1) * L), ncol = L)
  
  ##Sum of square
  ss = rep(0, p)
  for (j in 1:p){
    ss[j] = sum(X[ ,j]^2)
  }
  
  for (l in 1:L){
    lambda = lambda_all[l]
    for (t in 1:T){
      for (k in 1:p){
        db = sum(R * X[ ,k]) / ss[k]
        b = beta[k+1] + db
        b = sign(b) * max(0, abs(b) - lambda / ss[k])
        db = b - beta[k+1]
        R = R - X[ ,k] * db
        beta[k+1] = b
      }
    }
    beta[1] = sum(Y - X %*% beta[2:(p+1)]) / (p ** 2)
    beta_all[ ,l] = beta
  }
  return(beta_all)
}