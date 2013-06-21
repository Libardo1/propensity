## Specify the "correct" ps model (this is real data application, so we do not
## know what is correct or incorrect model)
#
# This requires a setup object to be in the .GlobalEnv for the logistic regression
# parameters.
correct.ps.model <-
function( dat ) {
  fm <- glm( as.formula(paste( "Z ~", paste(setup$X4PS.name, collapse=" + ") )),
             data = dat, family = binomial(link = "logit") ) ;
  ps.hat <- as.numeric( predict( fm, newdata = dat, type="response" ) ) ;
  ps.hat <- pmin( pmax( 0.000001, ps.hat), 0.999999 ) ;
  # ps.hat can not be exactly 0 or 1
  return( list(ps.hat = ps.hat, fm = fm) ) ;
}
