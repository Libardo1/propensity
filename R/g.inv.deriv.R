g.inv.deriv <-
function( x, fun.type ) {
  if ( length(x) != length(fun.type) ) {
    print("*** ERROR in g.inv.deriv() ***") ;
    return( rep(NA, length(x)) ) ;
  } else {
    m <- length(x) ;
    ans <- rep(0, m) ;
    for (j in 1:m) {
      if ( fun.type[j] == "log" ) {
        ans[j] <- exp( x[j] ) ;
      } else if ( fun.type[j] == "logit" ) {
        tmp <- inv.logit( x[j] ) ;
        ans[j] <- tmp*(1-tmp) ;
      } else {  # identity transformation, by default
        ans[j] <- 1 ;
      }
    }
    return( ans ) ;
  }
}
