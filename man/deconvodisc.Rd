\name{deconvodisc}
\alias{deconvodisc}
\title{Discrete deconvolution}
\description{
Given samples from an unknown integer-valued distribution, which comes from the convolution of an unknown underlying integer-valued distro and an arbitrary-but-known (discrete) convolution matrix, this estimates the "deconvolved" underlying distro by assuming it's smooth. Can work on multiple datasets from multiple underlying distros (but all having the same convolution matrix, and all are fitted with the same estimated smoothness parameter).

Deconvolution used to be thought of as a very hard statistical problem, but the advent of "regularized methods" (to stop the underlying distro having arbitrarily wiggly solutions) has made it pretty do-able. Here package \pkg{mgcv} is used to build the underlying smooth model; a thin layer of C code using the \pkg{TMB} package overlays it, to add the convolution step.

The (co)variance matrix of the fitted distro (separate for each dataset) is also returned. It comes from \code{TMB::sdreport(...,ignore.parm.uncertainty=TRUE)} which ensures that the estimates across datasets appear independent; in particular, the covariances are computed conditional on the estimated smoopar (i.e. wiggliness penalty). The overall intercepts for each dataset also appear in the "outer" parameters, and are neglected covariance-wise; I \emph{think} that's legit for Poisson/multinomial-equivalence models, but I should really check.
\subsection{To do}{
\itemize{
\item Allow arbitrary and more complicated smoothers, eg random-effect on top of general smoother. Easy enough, at least as long as there's only one smoopar.
\item Separate known convolution matrix for each dataset. (But you can always call \code{deconvodisc} several times, once for each convmat, so that's only relevant if you desperately want to share smoopars across datasets.)
}
}
}
\usage{deconvodisc( count, convmat, silent=TRUE)
}
\arguments{
\item{ count}{An integer-valued vector (or matrix), of bin-counts. If it's a matrix, then each column is treated as a separate dataset, from a separate underlying distro.}
\item{ convmat}{A numeric matrix, of dimension (length(count) X p) where \code{p} is the length of the support of the underlying distro. Note that the latter can be shorter than the observed sample.}
\item{ silent}{passed to \code{TMB::MakeADFun}. Maybe set to FALSE if you are having horrible problems.}
}
\value{A list with elements \code{pr_truhat}, \code{cov}. The former is a matrix with each column (one per dataset) being the estimated probabilities for that underlying distro. The latter is an array with the covariance matrix of the probabilities for the i-th dataset in \code{cov[,,i]}.
}
\examples{
## See "deconvo_test.r" and "deconvo_multi_test.r"
# They are in the namespace, tho should really be in the docs folder
\dontrun{
# Perhaps run them yourself manually, or via the debug package
# like this (but you won't see the comments)
# NFI why I set up the example this way!
#
if( require( debug))\{
  mdrun( deconvodisc:::deconvo_multi_test.r)
\}
}
}
\keyword{misc}
