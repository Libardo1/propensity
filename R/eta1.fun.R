eta1.fun <-
function( ps ) {
    ans <- rep(1, length(ps) ) ;
    tmp <- ps > 1/2 ;
    ans[ tmp ] <- 1/ps[ tmp ] - 1 ;
    ans ;
}
