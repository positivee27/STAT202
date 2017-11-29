#' Perform linear regression of Y on X using QR decomposition
#'
#' @param matrix X is an n x p matrix of explanatory variables
#' @param matrix Y is an n dimensional vector of responses
#' @return a list where the first element is a 1 x (p + 1) vector beta_ls
#'         and the second element is the standard error
#' @export

myLM <- function(X, Y){
  
  n = nrow(X)
  p = ncol(X)
  
  Z = cbind(c(rep(1, n)), X, Y)
  qr_result = myQR(Z)
  R = qr_result$R
  R1 = R[1:(p+1), 1:(p+1)]
  Y1 = R[1:(p+1), (p+2)]
  beta_ls = solve(R1, Y1)
  
  Yhat = cbind(c(rep(1, n)), X) %*% beta_ls
  se = sqrt(diag(sum((Y - Yhat) ** 2) / (n - p - 1) * solve(t(X) %*% X)))
  
  output <- list(beta_ls = beta_ls, se = se)
  return(output)
}