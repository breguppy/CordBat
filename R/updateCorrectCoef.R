#' @useDynLib CordBat, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
NULL

#' @export
#' @rdname update.CorrectCoef
#' @title Batch‐effect `a`/`b` updater (R wrapper)
#' @inheritParams updateCorrectCoefCpp
#' @param print.detail logical; echo a message?
#' @return List with components `coef.a` and `coef.b`
update.CorrectCoef <- function(X0.glist,
                               X1.glist,
                               Theta.list,
                               a.i, b.i,
                               penal.ksi,
                               penal.gamma,
                               print.detail = TRUE) {
  # call into C++
  res <- updateCorrectCoefCpp(
    X0.glist,
    X1.glist,
    Theta.list,
    as.numeric(a.i),
    as.numeric(b.i),
    penal.ksi,
    penal.gamma
  )
  
  return(list(coef.a = res$coef.a, coef.b = res$coef.b))
}