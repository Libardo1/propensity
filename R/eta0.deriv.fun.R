eta0.deriv.fun <-
function( ps ) {
    ans <- rep( 0, length(ps) ) ;
    tmp <- ps <= 1/2 ;
    ans[tmp] <- 1/(( 1 - ps[tmp] )^2) ;
    ans ;
}
