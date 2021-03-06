#' Convert symmetric sparse matrix to full matrix.
#'
#' Convert a symmetric sparse matrix sA to full matrix A.
#'
#' @export full
#' @importFrom methods is
#' @param sA Sparse matrix to convert.
#' @seealso \code{\link{fullup}}
#' @examples
#' sA <- pittsburgh$sA
#' A <- full(sA)

full <- function(sA) {
  if(is(sA, 'sparseMatrix')){
    A <- as.matrix(sA)
  } else {
    if (sA[1,1]!=1){
      ndim <- length(unique(sA[,2]))
    }
    else
    {
      ndim <- length(unique(sA[,1]))
    }
    A <- matrix(0, ncol = ndim, nrow = ndim)
    A[cbind(sA[,1],sA[,2])] <- sA[,3]
  }
  return(A)
}
