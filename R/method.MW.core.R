#' the core function of MW
#'
#'
#' @export method.MW.core
method.MW.core <- function( dat, out.ps ) {
    dat$ps.hat <- out.ps$ps.hat ;
            # use estimated PS
    numer <- pmin(1-dat$ps.hat, dat$ps.hat) ;
    denom <- (1-dat$Z)*(1-dat$ps.hat) + dat$Z*dat$ps.hat ;
    dat$W <- numer/denom ;
            # estimated matching weights
    tmp <- dat$W*dat$Z ;
    first.term <- sum( tmp*dat$Y )/sum( tmp ) ;
    tmp <- dat$W*(1-dat$Z) ;
    second.term <- sum( tmp*dat$Y )/sum( tmp ) ;
    est <- first.term - second.term ;
            # point estimator
    mu1.hat <- first.term ;
    mu0.hat <- second.term ;
    beta.hat <- as.numeric(coef(out.ps$fm)) ;
            # prepare for calculating std
    if ( setup$calc.std ) {
        # std estimator from sandwich method
        n <- nrow(dat) ;
        Amat <- Bmat <- 0 ;
        for (i in 1:n) {
            Xi <- as.numeric( dat[i, c("X0", names(coef(out.ps$fm))[-1])] ) ;
                    # "X0" is the key word for intercept
            Zi <- dat[i, "Z"] ;
            Yi <- dat[i, "Y"] ;
            Wi <- calc.MW( Xi, Zi, beta.hat ) ;
            ei <- calc.ps.Xbeta( Xi, beta.hat ) ;
            ei.deriv <- calc.ps.deriv( Xi, beta.hat ) ;
            ei.deriv2 <- calc.ps.deriv2( Xi, beta.hat ) ;
            this.phi <- c( Wi*Zi*(Yi-mu1.hat) ,
                           Wi*(1-Zi)*(Yi-mu0.hat) ,
                           (Zi-ei)/ei/(1-ei)*ei.deriv ) ;
            Bmat <- Bmat + outer( this.phi, this.phi ) ;

            eta.deriv.i <- eta.deriv.Ze( Zi, ei ) ;
            first.row <- c( -Wi*Zi, 0,
                                eta.deriv.i*Zi*(Yi-mu1.hat)*ei.deriv ) ;
            second.row <- c( 0, -Wi*(1-Zi),
                                eta.deriv.i*(1-Zi)*(Yi-mu0.hat)*ei.deriv ) ;
            tmp0 <- matrix(0, nrow=length(beta.hat), ncol=2) ;
            tmp1 <- (Zi-ei)/ei/(1-ei) ;
            tmp2 <- tmp1*ei.deriv2 - tmp1*tmp1*(ei.deriv %*%
                                                 rowVec(ei.deriv)) ;
            third.row <- cbind( tmp0, tmp2 ) ;
            this.phi.deriv <- rbind( first.row, second.row, third.row ) ;
            Amat <- Amat + this.phi.deriv ;
        }
        Bmat <- Bmat/n ;
        Amat <- Amat/n ;
        Amat.inv <- solve(Amat) ;
        var.mat <- ( Amat.inv %*% Bmat %*% t(Amat.inv) )/n ;
        tmp1 <- c(1, -1, rep(0, length(beta.hat))) ;
        var.est <- rowVec(tmp1) %*% var.mat %*% colVec(tmp1) ;
        std <- sqrt( as.numeric(var.est) ) ;
        CI.lower <- est - 1.96*std ;
        CI.upper <- est + 1.96*std ;
    } else {
        std <- CI.lower <- CI.upper <- NA ;
    }
    effect.n1 <- sum( dat$W[ dat$Z == 1 ] ) ;
    effect.n0 <- sum( dat$W[ dat$Z == 0 ] ) ;
                # effective sample size
    ans <- list( est = est, std = std,
                 CI.lower = CI.lower, CI.upper = CI.upper,
                 effect.n1 = effect.n1, effect.n0 = effect.n0 ) ;
    ans ;
}
