solve.ab <-
function( delta = 0.05 ) {
    if ( delta < 0.001 ) { stop("*** ERROR in solve.ab: delta too small ***") }
    tmp1 <- matrix( c( 1, (0.5-delta), (0.5-delta)^2, (0.5-delta)^3 ,
                       0, 1, 2*(0.5-delta), 3*(0.5-delta)^2 ,
                       1, (0.5+delta), (0.5+delta)^2, (0.5+delta)^3 ,
                       0, 1, 2*(0.5+delta), 3*(0.5+delta)^2 ) ,
                    nrow = 4, ncol = 4, byrow=T ) ;
    tmp2 <- colVec( c( 1, 0, eta1.fun(0.5+delta),
                       eta1.deriv.fun(0.5+delta) ) ) ;
    tmp3 <- colVec( c( eta0.fun(0.5-delta), eta0.deriv.fun(0.5-delta),
                       1, 0 ) ) ;
    a <- as.numeric( solve(tmp1) %*% tmp2 ) ;
    b <- as.numeric( solve(tmp1) %*% tmp3 ) ;
    list(a = a, b = b)  ;
}
