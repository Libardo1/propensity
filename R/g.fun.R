g.fun <-
function( x, fun.type ) {
  if ( length(x) != length(fun.type) ) {
    print("*** ERROR in g.fun() ***") ;
    return( rep(NA, length(x)) ) ;
  } else {
    m <- length(x) ;
    ans <- rep(0, m) ;
    for (j in 1:m) {
      if ( fun.type[j] == "log" ) {
        ans[j] <- log( max( 1e-6, x[j] ) ) ;
      } else if ( fun.type[j] == "logit" ) {
        ans[j] <- logit( min( max( 1e-6, x[j]), 1-1e-6 ) )
      } else {  # identity transformation, by default
        ans[j] <- x[j] ;
      }
    }
    return( ans ) ;
  }
}
