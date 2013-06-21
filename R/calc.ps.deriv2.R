# They are used for variance estimation purposes
#
## Calculate the second derivative of propensity score w.r.t beta,
## i.e., e_beta_beta(Xi), where Xi is a vector (not a matrix!)
## beta is the coefficient, of the same length as Xi
## return a square matrix of length(Xi) by length(beta)
calc.ps.deriv2 <-
function( Xi, beta ) {
    Xi <- colVec(Xi) ;
    tmp.ps <- calc.ps.Xbeta( Xi, beta ) ;
    tmp.deriv <- calc.ps.deriv( Xi, beta ) ;
    ans <- Xi %*% rowVec(tmp.deriv) ;
    names(ans) <- rownames(ans) <- NULL ;
    ans <- ( 1-2*tmp.ps )*ans ;
    return( ans ) ;
}
