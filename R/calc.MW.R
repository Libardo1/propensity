
## Calculate the matching weight as a function of X, Z, and beta
## Xmat is a matrix with each column being Xi
## Z is a vector of length ncol(Xmat)
## beta is a vector of coefficients of length nrow(Xmat)
## Return a vector of length ncol(Xmat)

calc.MW <- function( Xmat, Z, beta ) {
    Xmat <- as.matrix(Xmat) ;
    tmp.ps <- calc.ps.Xbeta( Xmat, beta ) ;
    eta.Ze( Z, tmp.ps ) ;
}
