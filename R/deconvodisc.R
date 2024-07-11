# This is package deconvodisc 

"deconvodisc" <-
function( count, convmat, silent=TRUE) {
## convmat acts on a _range_ of allowed inputs, and maps them...
## ... to a (presumably bigger) _range_ of outputs

  ndomain <- ncol( convmat)
  nrange <- nrow( convmat)
  
  DROP <- !length( dim( count))
  if( DROP){
    count <- matrix( count, length( count), 1L)
  } else {
stopifnot( length( dim( count))==2)
  }
stopifnot( nrow( count) == nrange)

  ndatasets <- ncol( count)

  faketru <- rep.int( 0L, ndomain)
  xtru <- seq_len( ndomain)

  flunge <- gam( faketru ~ s( xtru), family=poisson( link=log), fit=FALSE)

# AFAICS
# $off: first penalized coef (and col in X)
# $S: smoothing matrix
# $rank: rank of S
# Not sure what happens to off etc with Xtuple smoothers
#
# For nullspace entries (aside from defined "fixed effects") then rank < dim( S) and we could use a pivoted Cholesky decomp
# but actually the rowcol of S will be identically zero I _think_. Piv chol would be safer, but more complicated.

  S <- flunge$S[[1]] # assume only 1
  X <- flunge$X
  off <- flunge$off
  ifix <- 1 %upto% (off-1)

  nfix <- length( ifix) # better be 1!
stopifnot( nfix==1)

  nrand <- ncol( X) - nfix
  
  nullcols_S <- which( rowSums( abs( S))==0)
  S_nonnull <- if( length( nullcols_S)) S[ -nullcols_S, -nullcols_S] else S
  chol_Snn <- chol( S_nonnull)
  chol_S <- 0 * S
  chol_S[ 1:ncol( S_nonnull), 1:ncol( S_nonnull)] <- chol_Snn

  irand_reord <- off-1 + c( 1:nrand %except% nullcols_S, nullcols_S)

  # Don't bother with sparse matrices: small and dense is likely...
  # ... and 'Matrix' is a bit of a PITA
  
  Desmat_fix <- X[, ifix, drop=FALSE]
  Desmat_rand <- X[, irand_reord, drop=FALSE]

  
  # Start with intercept (presumably first col) only
  # Single dataset was:
  # beta_start <- c( log( sum( count)), rep( 0, nfix-1))
  # u_start <- rep( 0, nrand)
  beta_start <- matrix( 0, nfix, ndatasets)
  beta_start[1,] <- log( colSums( count))
  u_start <- matrix( 0, nrand, ndatasets)

  # dyn.load( dynlib( 'ddeconvo')) while testing
  obj <- MakeADFun(
    data=returnList( count, Desmat_fix, Desmat_rand, chol_S, M=convmat),
    parameters=list(u=u_start, beta=beta_start, logsdu=1),
    random="u",
    DLL="deconvodisc"
    ,silent=silent
  )
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  repo <- sdreport( obj, ignore.parm.uncertainty=TRUE)
  u_best <- summary( repo, 'random')[,1]
  dim( u_best) <- dim( u_start)
  beta_best <- summary( repo, 'fixed')[,1] %such.that% (names(.)=='beta')
  dim( beta_best) <- dim( beta_start)
  
  cov <- array( 0, c( ndomain, ndomain, ndatasets))
  starto <- 0
  for( i in 1:ndatasets) {
    cov[,,i] <- repo$cov[ starto + 1:ndomain, starto+1:ndomain]
  }
  
  logtruhat <- Desmat_fix %*% beta_best + Desmat_rand %*% u_best
  
  truhat <- exp( logtruhat)
  pr_truhat <- truhat / rep( colSums( truhat), each=nrow( truhat))
  if( DROP) {
    pr_truhat <- c( pr_truhat)
  }
  
returnList( pr_truhat, cov) # don't bother with SDs etc  
}


"orig_deconvodisc" <-
function( count, convmat, silent=TRUE) {
## This is the original version, which moves any nullspace cfts in the smoother...
## into the outer params. Prolly better not to.

## convmat acts on a _range_ of allowed inputs, and maps them...
## ... to a (presumably bigger) _range_ of outputs

  ndomain <- ncol( convmat)
  nrange <- nrow( convmat)
  
  DROP <- !length( dim( count))
  if( DROP){
    count <- matrix( count, length( count), 1L)
  } else {
stopifnot( length( dim( count))==2)
  }
stopifnot( nrow( count) == nrange)

  ndatasets <- ncol( count)

  faketru <- rep.int( 0L, ndomain)
  xtru <- seq_len( ndomain)

  flunge <- gam( faketru ~ s( xtru), family=poisson( link=log), fit=FALSE)

#	AFAICS
#	$off: first penalized coef (and col in X)
#	$S: smoothing matrix
#	$rank: rank of S
# Not sure what happens to off etc with Xtuple smoothers
#
#	For nullspace entries (aside from defined "fixed effects") then rank < dim( S) and we could use a pivoted Cholesky decomp
#	but actually the rowcol of S will be identically zero I _think_. Piv chol would be safer, but more complicated.

  S <- flunge$S[[1]] # assume only 1
  X <- flunge$X
  off <- flunge$off
  nullcols_S <- which( rowSums( abs( S))==0)
  S_nonnull <- if( length( nullcols_S)) S[ -nullcols_S, -nullcols_S] else S
  chol_S <- chol( S_nonnull)

  reord <- c( 1 %upto% (off-1), off-1 + nullcols_S)

  # Don't bother with sparse matrices: small and dense is likely...
  # ... and 'Matrix' is a bit of a PITA
  
  Desmat_fix <- X[, reord]
  Desmat_rand <- X[, -reord]

  nfix <- off-1 + length( nullcols_S)
  nrand <- ncol( X) - nfix
  
  # Start with intercept (presumably first col) only
  # Single dataset was:
  # beta_start <- c( log( sum( count)), rep( 0, nfix-1))
  # u_start <- rep( 0, nrand)
  beta_start <- matrix( 0, nfix, ndatasets)
  beta_start[1,] <- log( colSums( count))
  u_start <- matrix( 0, nrand, ndatasets)

  # dyn.load( dynlib( 'ddeconvo')) while testing
  obj <- MakeADFun(
    data=returnList( count, Desmat_fix, Desmat_rand, chol_S, M=convmat),
    parameters=list(u=u_start, beta=beta_start, logsdu=1),
    random="u",
    DLL="deconvodisc"
    ,silent=silent
  )
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  repo <- sdreport( obj, ignore.parm.uncertainty=TRUE)
  u_best <- summary( repo, 'random')[,1]
  dim( u_best) <- dim( u_start)
  beta_best <- summary( repo, 'fixed')[,1] %such.that% (names(.)=='beta')
  dim( beta_best) <- dim( beta_start)
  
  cov <- array( 0, c( ndomain, ndomain, ndatasets))
  starto <- 0
  for( i in 1:ndatasets) {
    cov[,,i] <- repo$cov[ starto + 1:ndomain, starto+1:ndomain]
  }
  
  logtruhat <- Desmat_fix %*% beta_best + Desmat_rand %*% u_best
  
  truhat <- exp( logtruhat)
  pr_truhat <- truhat / rep( colSums( truhat), each=nrow( truhat))
  if( DROP) {
    pr_truhat <- c( pr_truhat)
  }

returnList( pr_truhat, cov) # don't bother with SDs etc  
}

