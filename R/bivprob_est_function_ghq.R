#' Internal function
#'
#' This function is a wrapper for estimating the censored bivariate 
#' probit model with random intercept using GHQ. Calculates starting 
#' values and conducts numerical optimization.
#'
#' @keywords internal
#' 
#' @importFrom Rcpp sourceCpp
#' @importFrom utils head 
#' @importFrom stats as.formula binomial dnorm glm optim pnorm runif rnorm coef
#' @importFrom lme4 glmer
#' @importFrom parallel makeCluster setDefaultCluster stopCluster
#' @importFrom maxLik numericHessian
#' @importFrom compiler cmpfun
#' @importFrom statmod gauss.quad
#' @importFrom Matrix forceSymmetric nearPD
#' @importFrom corpcor make.positive.definite
loglik_bivprob_est_ghq <- function(data, yS, yO, XO, XS, kS, kO, group,
                                   sel, out, QP,...) {
  
  r_start <- suppressWarnings(lme4::glmer(sel,data=data,
                                          family=binomial(link="probit")))
  sigma1 <- r_start@theta^2  
  if(round(sigma1,6)==0){sigma1 <- 0.0001}
  
  # starting values outcome equation: 
  y_start <- suppressWarnings(lme4::glmer(out,
                                          data=data[which(as.numeric(yS)==1),],
                                          family=binomial(link="probit")))
  sigma2 <- y_start@theta^2 
  if(round(sigma2,6)==0){sigma2 <- 0.0001}
  

  start <- c(r_start@beta, y_start@beta, atanh(0.5),
             log(sigma1), log(sigma2), atanh(0.5))

  opt <- optim(par=start, fn=loglik_bivprob_ghq, gr=NULL, 
               method="BFGS",
               yS, yO, XO, XS, kS, kO, group,QP)
 
  hess <- maxLik::numericHessian(loglik_bivprob_ghq,gr=NULL,
                                 t0=opt$par, 
                                 yS=yS, yO=yO, XO=XO, XS=XS, 
                                 kS=kS, kO=kO, group=group,QP=QP)
  
  res <- list(estimates=opt$par,hessian=hess)
  
  return(res)
}

loglik_bivprob_est_ghq <- compiler::cmpfun(loglik_bivprob_est_ghq)
