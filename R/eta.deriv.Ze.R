
## Calculate the derivative of eta(Z, e)
## Z and e may be vectors, but must of the same lengths
## Return a vector of the same length as Z and e
eta.deriv.Ze <-
function( Z, ps ) {
    ans <- rep(0, length(ps)) ;
    s1 <- Z == 1 ; s0 <- !s1 ;
    if ( any(s1) ) { ans[s1] <- eta1.deriv.fun.approx( ps[s1], setup$delta ) }
    if ( any(s0) ) { ans[s0] <- eta0.deriv.fun.approx( ps[s0], setup$delta ) }
    ans ;
}
