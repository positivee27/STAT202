#' Perform PCA on matrix A using QR
#'
#' @param matrix A is an square matrix
#' @param int numIter is the number of iterations
#' @return a list with D and V, where D is a vector of eigenvalues of A 
#'         and V is the matrix of eigenvectors of A
#' @export

myEigen_QR <- function(A, numIter = 1000){
  
  r = nrow(A)
  c = ncol(A)
  V = matrix(rnorm(r * r), nrow = r)
  R = matrix()
  Q = matrix()
  for (i in 1:numIter){
    qr_result = myQR(V)
    Q = qr_result$Q
    R = qr_result$R
    V = A %*% Q
  }
  return(list("D" = diag(R), "V" = Q))
}
  