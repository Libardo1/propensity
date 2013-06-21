## Specify the "correct" Y model (this is real data application, so we do not
## know what is correct or incorrect model)
##
## Again requires a setup object in the .GlobalEnv for logistic regression
correct.Y.model <-
function( dat ) {
  fm1 <- lm( as.formula(paste( "Y ~", paste(setup$X4Y.name, collapse=" + ") )),
             data = dat[ dat$Z == 1, ] ) ;
  fm0 <- lm( as.formula(paste( "Y ~", paste(setup$X4Y.name, collapse=" + ") )),
             data = dat[ dat$Z == 0, ] ) ;
  Y.hat.m1 <- as.numeric( predict( fm1, newdata = dat ) ) ;
  Y.hat.m0 <- as.numeric( predict( fm0, newdata = dat ) ) ;
  return( list( Y.hat.m1 = Y.hat.m1 ,
                Y.hat.m0 = Y.hat.m0 ,
                fm1 = fm1 , fm0 = fm0 ) ) ;
}
