eta0.deriv.fun.approx <-
function( ps, delta = 0.05 ) {
    b <- solve.ab( delta )$b ;
    ans <- eta0.deriv.fun( ps ) ;
    s <- 0.5-delta <= ps & ps <= 0.5+delta ;
    ans[s] <- b[2] + 2*b[3]*ps[s] + 3*b[4]*ps[s]^2 ;
    ans ;
}
