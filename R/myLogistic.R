expit <- function(x){
  1 / (1 + exp(-x))
}

#' Perform the logistic regression of Y on X
#'
#' @param matrix X is an n x p matrix of explanatory variables
#' @param matrix Y is an n dimensional vector of binary responses
#' @param int numIter is the number of iterations
#' @return a list where the first element is a 1 x p vector beta and the second
#'         element is the standard error
#' @export

myLogistic <- function(X, Y){
  
  r = nrow(X)
  c = ncol(X)
  beta = matrix(rep(0, c), nrow = c, ncol = 1)
  err = 1
  
  while (err > 0.000001){
    eta = X %*% beta
    pr = 1 / (1 + exp(-eta))
    w = pr * (1 - pr)
    Z = eta + (Y - pr) / w
    sw = sqrt(w)
    mw = matrix(rep(sw, c), ncol = c)
    X_work = mw * X
    Y_work = sw * Z
    beta_new = myLM(X_work, Y_work)$beta_ls[2:(c+1)]
    err = sum(abs(beta_new - beta))
    beta = beta_new
  }
  
  Yhat = X %*% beta
  se = sqrt(diag(solve(t(X_work) %*% X_work)))
  
  output <- list(beta = beta, se = se)
  return(output)
}