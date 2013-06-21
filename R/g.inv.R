g.inv <-
function( x, fun.type ) {
  if ( length(x) != length(fun.type) ) {
    print("*** ERROR in g.inv() ***") ;
    return( rep(NA, length(x)) ) ;
  } else {
    m <- length(x) ;
    ans <- rep(0, m) ;
    for (j in 1:m) {
      if ( fun.type[j] == "log" ) {
        ans[j] <- exp( x[j] ) ;
      } else if ( fun.type[j] == "logit" ) {
        ans[j] <- inv.logit( x[j] ) ;
      } else {  # identity transformation, by default
        ans[j] <- x[j] ;
      }
    }
    return( ans ) ;
  }
}
