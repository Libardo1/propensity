# They are used for variance estimation purposes
#
## Calculate the derivative of propensity score w.r.t. beta, e_beta(X)
## X is a matrix with each column being X_i
## beta is the coefficient, of the same length as the nrow(X)
## return a matrix of the same dimension as X
calc.ps.deriv <-
function( Xmat, beta ) {
    Xmat <- as.matrix(Xmat) ;
    tmp.ps <- calc.ps.Xbeta( Xmat, beta ) ;
    ans <- t(Xmat)*(tmp.ps*(1-tmp.ps)) ;
    names(ans) <- rownames(ans) <- NULL ;
    return( t(ans) ) ;
}
