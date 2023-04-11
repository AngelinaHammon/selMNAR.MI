#' Internal function
#'
#' This function is a wrapper for estimating the ordered probit model with sample selection and random intercept using GHQ. Calculates starting values and conducts numerical optimization.
#'
#' @keywords internal
#' 
#' @importFrom ordinal clmm VarCorr
loglik_ordbivprob_est_ghq <- function(data, yS, yO, XO, XS, kS, kO, group, sel, out, QP,...) {
  
  r_start <- lme4::glmer(sel,data=data,family=binomial(link="probit"))
  sigma1 <- r_start@theta^2 
  if(sigma1==0){sigma1 <- 0.0001}
  
  # starting values outcome equation: 
  y_start <- ordinal::clmm(out,data=data[which(as.numeric(yS)==1),], link="probit")  
  sigma2 <- ordinal::VarCorr(y_start)$'group'[1]
  if(sigma2==0){sigma2 <- 0.0001}
  
  z_start <- y_start$alpha
  z_trans <- z_start
  
  # transform them: 
  for(i in 2:length(z_start)){
    z_trans[i] <- log(z_start[i]-z_start[i-1])
  }
  
  
  start <- c(r_start@beta,y_start$beta, atanh(0.5), z_trans, log(sigma1), log(sigma2) ,atanh(0.5))
  
  opt <- optim(par=start, fn=loglik_ordbivprob_ghq, gr=NULL, method="BFGS", yS, yO, XO, XS, kS, kO, group,QP)
  
  hess <- maxLik::numericHessian(loglik_ordbivprob_ghq,gr=NULL,t0=opt$par, yS=yS, yO=yO, XO=XO, XS=XS, kS=kS, kO=kO, group=group,QP=QP)
  res <- list(estimates=opt$par,hessian=hess)
  
  return(res)
}

loglik_ordbivprob_est_ghq <- compiler::cmpfun(loglik_ordbivprob_est_ghq)
