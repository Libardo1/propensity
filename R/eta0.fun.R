eta0.fun <-
function( ps ) {
    ans <- rep(1, length(ps) ) ;
    tmp <- ps <= 1/2 ;
    ans[ tmp ] <- ps[ tmp ]/( 1 - ps[ tmp ] ) ;
    ans ;
}
