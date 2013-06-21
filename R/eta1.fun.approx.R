eta1.fun.approx <-
function( ps, delta = 0.05 ) {
    a <- solve.ab( delta )$a ;
    ans <- eta1.fun( ps ) ;
    s <- 0.5-delta <= ps & ps <= 0.5+delta ;
    ans[s] <- a[1] + a[2]*ps[s] + a[3]*ps[s]^2 + a[4]*ps[s]^3 ;
    ans ;
}
