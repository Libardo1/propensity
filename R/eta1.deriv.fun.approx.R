eta1.deriv.fun.approx <-
function( ps, delta = 0.05 ) {
    a <- solve.ab( delta )$a ;
    ans <- eta1.deriv.fun( ps ) ;
    s <- 0.5-delta <= ps & ps <= 0.5+delta ;
    ans[s] <- a[2] + 2*a[3]*ps[s] + 3*a[4]*ps[s]^2 ;
    ans ;
}
