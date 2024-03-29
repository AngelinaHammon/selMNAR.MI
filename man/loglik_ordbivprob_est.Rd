% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mice_heckman_ordered.R
\name{loglik_ordbivprob_est}
\alias{loglik_ordbivprob_est}
\alias{loglik_ordbivprob}
\alias{loglik_ordbivprob_grad}
\title{Internal function}
\usage{
loglik_ordbivprob_est(data, yS, yO, XO, XS, kS, kO, sel, out, ...)

loglik_ordbivprob(param, yS, yO, XO, XS, kS, kO)

loglik_ordbivprob_grad(param, yS, yO, XO, XS, kS, kO)
}
\description{
This function is a wrapper for estimating the ordered probit model with sample selection. Calculates starting values and conducts numerical optimization. The other two functions calculate the log likelihood and analytic gradients.
}
\keyword{internal}
