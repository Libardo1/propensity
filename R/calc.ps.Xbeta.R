# They are used for variance estimation purposes
#
## Calculate the propensity score as e(X, beta)
## X is a vector or a matrix with each column being X_i
## beta is the coefficient, of the same length as the nrow(X)
## return a vector of length ncol(X)
calc.ps.Xbeta <-
function( Xmat, beta ) {
    Xmat <- as.matrix(Xmat) ;
    tmp <- as.numeric( rowVec(beta) %*% Xmat ) ;
    tmp <- exp(tmp) ;
    names(tmp) <- NULL ;
    return( tmp/(1+tmp) ) ;
}
