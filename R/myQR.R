#' Perform QR decomposition on matrix A
#'
#' @param matrix input n x m matrix A
#' @return a list with Q.transpose and R, where Q is an orthogonal n x n matrix and
#'         R is an upper triangular n x m matrix
#' @export

myQR <- function(A){
  
  n = nrow(A)
  m = ncol(A)
  R = A
  Q = diag(x = 1, nrow = n, ncol = n)
  for (k in 1:(m - 1)){
    X = matrix(rep(0, n), nrow = n)
    X[k:n, 1] = R[k:n, k]
    V = X
  
    V[k,1] = X[k,1] + sign(X[k]) * norm(X, type = "F")
    S = norm(V, type = "F")
    U = V / S
    
    R = R - 2 * U %*% t(U) %*% R
    Q = Q - 2 * U %*% t(U) %*% Q
  }
  return(list("Q" = t(Q), "R" = R))
}