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
  X = cbind(rep(1, n), X)
  p = ncol(X)
  L = length(lambda_all)
  beta_all = matrix(rep(0, p*L), nrow = p)
  lambda_all = sort(lambda_all, decreasing = TRUE)
  T = 100
  R = Y ##all beta = 0, residual = Y
  beta = c(rep(0, p)) ##?
  
  ss = rep(0, p)
  
  for (j in 1:p){
    ss[j] = sum(X[ ,j]^2)
  }
  
  for (l in 1:L){
    lambda = lambda_all[l]
    for (t in 1:T){
      for (k in 1:p){
        db = sum(R * X[ ,k]) / ss[k] ##compute residual without excluding beta_k
        b = beta[k] + db ##
        if (k != 1){
          b = sign(b) * max(0, abs(b) - lambda / ss[k]) ##regularization
        }
        db = b - beta[k]
        R = R-X[ ,k] * db
        beta[k] = b
      }
    }
    #beta[1] = sum(R * X[ ,1]) / ss[1]
    beta_all[ ,l] = beta
  }
  return(beta_all)
}