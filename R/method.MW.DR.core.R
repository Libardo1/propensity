#' the core function of DR MW
#'
#' @export method.MW.DR.core
method.MW.DR.core <- function( dat, out.ps, out.Y ) {
    dat$ps.hat <- out.ps$ps.hat ;
    dat$Y.hat.m1 <- out.Y$Y.hat.m1 ;
    dat$Y.hat.m0 <- out.Y$Y.hat.m0 ;
            # estimated PS and estimated Y for Z == 1 and Z == 0
    numer <- pmin(1-dat$ps.hat, dat$ps.hat) ;
    denom <- (1-dat$Z)*(1-dat$ps.hat) + dat$Z*dat$ps.hat ;
    dat$W <- numer/denom ;
            # estimated matching weights
    mu1.hat <- sum( dat$W*(dat$Y.hat.m1 - dat$Y.hat.m0) )/sum( dat$W ) ;
    tmp <- dat$W*dat$Z ;
    mu2.hat <- sum( tmp*(dat$Y-dat$Y.hat.m1) )/sum( tmp ) ;
    tmp <- dat$W*(1-dat$Z) ;
    mu3.hat <- sum( tmp*(dat$Y-dat$Y.hat.m0) )/sum( tmp ) ;
    est <- mu1.hat + mu2.hat - mu3.hat ;
            # point estimator
    alpha1.hat <- as.numeric(coef(out.Y$fm1)) ;
    alpha0.hat <- as.numeric(coef(out.Y$fm0)) ;
    beta.hat <- as.numeric(coef(out.ps$fm)) ;
    n.alpha1 <- length(alpha1.hat) ;
    n.alpha0 <- length(alpha0.hat) ;
    n.beta <- length(beta.hat) ;
            # prepare for calculating std
    if ( setup$calc.std ) {
        # std estimator from sandwich method
        n <- nrow(dat) ;
        Amat <- Bmat <- 0 ;
        for (i in 1:n) {
            Xi <- as.numeric( dat[i, c("X0", names(coef(out.ps$fm))[-1])] ) ;
            Vi <- as.numeric( dat[i, c("X0", names(coef(out.Y$fm1))[-1])] ) ;
            Zi <- dat[i, "Z"] ;
            Yi <- dat[i, "Y"] ;
            Wi <- calc.MW( Xi, Zi, beta.hat ) ;
            ei <- calc.ps.Xbeta( Xi, beta.hat ) ;
            ei.deriv <- calc.ps.deriv( Xi, beta.hat ) ;
            ei.deriv2 <- calc.ps.deriv2( Xi, beta.hat ) ;

            this.phi.row1 <- Wi*( sum(Vi*(alpha1.hat-alpha0.hat)) - mu1.hat ) ;
            this.phi.row2 <- Wi*Zi*(Yi-sum(Vi*alpha1.hat)-mu2.hat) ;
            this.phi.row3 <- Wi*(1-Zi)*(Yi-sum(Vi*alpha0.hat)-mu3.hat) ;
            this.phi.row4 <- Zi*Vi*( Yi-sum(Vi*alpha1.hat) ) ;
            this.phi.row5 <- (1-Zi)*Vi*( Yi-sum(Vi*alpha0.hat) ) ;
            this.phi.row6 <- (Zi-ei)/ei/(1-ei)*ei.deriv ;
            this.phi <- c( this.phi.row1, this.phi.row2, this.phi.row3,
                           this.phi.row4, this.phi.row5, this.phi.row6 ) ;
            Bmat <- Bmat + outer( this.phi, this.phi ) ;

            eta.deriv.i <- eta.deriv.Ze( Zi, ei ) ;
            quad1 <- diag( c(-Wi, -Wi*Zi, -Wi*(1-Zi)) ) ;
            quad2 <- matrix(0, nrow=n.alpha1+n.alpha0+n.beta, ncol=3) ;
            tmp <- c( Wi*Vi, -Wi*Vi,
                eta.deriv.i*ei.deriv*(sum(Vi*(alpha1.hat-alpha0.hat))-mu1.hat),
                -Wi*Zi*Vi, rep(0, n.alpha0),
                eta.deriv.i*ei.deriv*Zi*(Yi-sum(Vi*alpha1.hat)-mu2.hat),
                rep(0, n.alpha1), -Wi*(1-Zi)*Vi,
                eta.deriv.i*ei.deriv*(1-Zi)*(Yi-sum(Vi*alpha0.hat)-mu3.hat)
                ) ;
            quad3 <- matrix( tmp, byrow = TRUE, nrow = 3,
                             ncol = n.alpha1 + n.alpha0 + n.beta ) ;
            quad4.blk1 <- -Zi*outer(Vi, Vi) ;
            quad4.blk2 <- -(1-Zi)*outer(Vi, Vi) ;
            tmp1 <- (Zi-ei)/ei/(1-ei) ;
            quad4.blk3 <- tmp1*ei.deriv2 -
                tmp1*tmp1*(ei.deriv %*% rowVec(ei.deriv)) ;
            quad4 <- as.matrix(bdiag( quad4.blk1, quad4.blk2, quad4.blk3 )) ;
            this.phi.deriv <- rbind( cbind(quad1, quad3) ,
                                     cbind(quad2, quad4) ) ;
            Amat <- Amat + this.phi.deriv ;
        }
        Bmat <- Bmat/n ;
        Amat <- Amat/n ;
        Amat.inv <- solve(Amat) ;
        var.mat <- ( Amat.inv %*% Bmat %*% t(Amat.inv) )/n ;
        tmp1 <- c(1, 1, -1, rep(0, length = n.alpha1+n.alpha0+n.beta) ) ;
        var.est <- rowVec(tmp1) %*% var.mat %*% colVec(tmp1) ;
        std <- sqrt( as.numeric(var.est) ) ;
        CI.lower <- est - 1.96*std ;
        CI.upper <- est + 1.96*std ;
    } else {
        std <- CI.lower <- CI.upper <- NA ;
    }
    ans <- list( est = est, std = std,
                 CI.lower = CI.lower, CI.upper = CI.upper ) ;
    ans ;
}
