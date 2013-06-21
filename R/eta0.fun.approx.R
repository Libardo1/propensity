eta0.fun.approx <-
function( ps, delta = 0.05 ) {
    b <- solve.ab( delta )$b ;
    ans <- eta0.fun( ps ) ;
    s <- 0.5-delta <= ps & ps <= 0.5+delta ;
    ans[s] <- b[1] + b[2]*ps[s] + b[3]*ps[s]^2 + b[4]*ps[s]^3 ;
    ans ;
}
