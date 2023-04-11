#' Internal function
#'
#' This function calculates the log likelihood of the ordered probit model with sample selection and random intercept using GHQ. Calls a C++ function to speed up the loop.
#'
#' @keywords internal
loglik_ordbivprob_ghq <- function(param, yS, yO, XO, XS, kS, kO,group,QP,...) {
  
  # get categories of yO:
  if(is.factor(yO)==T) {
    cat <- levels(yO)
  }else{
    cat <- unique(yO)[order(unique(yO))]
  }
  
  betaS    <- param[1:kS]
  betaO    <- param[(kS+1):(kS+kO)]
  rho      <- tanh(param[kS+kO+1])
  z        <- param[(kS+kO+2):(kS+kO+2+length(cat)-2)]
  sigma1   <- exp(param[(kS+kO+2+length(cat)-2+1)])
  sigma2   <- exp(param[(kS+kO+2+length(cat)-2+2)])
  tau      <- tanh(param[length(param)])

  # retransform:
  for(i in 2:length(z)){
    z[i] <- exp(z[i])+z[i-1]
  }
  
  group_aux <- c(0,cumsum(table(as.integer(factor(group)))))
  
  
  # calculate nodes and weights for given quadrature points:
  rule1 <- statmod::gauss.quad(QP[1], "hermite")
  rule <- list(w_1=rule1$weights, a_1=rule1$nodes)#, w_2=rule2$weights, a_2=rule2$nodes)
  
  w1 <- rule$w_1
  a1 <- rule$a_1
  
  M <- length(table(group)) 
  N <- length(yS)
  
  # linear predictors:
  linpredS <- XS %*% betaS 
  linpredO <- XO %*% betaO

  if(round(tau,2)==1){
    tau <- 0.9999
  }
  
  if(round(tau,2)==-1){
    tau <- -0.9999
  }
  
  if(round(rho,2)==1){
    rho <- 0.9999
  }
  
  if(round(rho,2)==-1){
    rho <- -0.9999
  }
  
  
  if(param[kS+kO+2]>700){
    sigma1 <- exp(700)
  }
  
  if(param[kS+kO+3]>700){
    sigma2 <- exp(700)
  }
  
  if(param[kS+kO+2] < -700){
    sigma1 <- exp(-700)
  }
  
  if(param[kS+kO+3] < -700){
    sigma2 <- exp(-700)
  }
  
  vc_re <- matrix(0,2,2)
  vc_re[1,1] <- sigma1
  vc_re[2,2] <- sigma2
  vc_re[1,2] <- vc_re[2,1] <- sqrt(sigma1)*sqrt(sigma2)*tau
  
  chol_vc_re <- t(chol(vc_re))

  yO <- as.numeric(as.character(yO))
  
  levels <- as.numeric(as.character(unique(yO[order(yO)])))
  
  loglik <- loop1_ord(yS, yO, linpredS, linpredO, rho, tau, sigma1, sigma2, z, levels, w1, a1, group_aux, group, M, chol_vc_re)

  # Total Log-Likelihood:
  return(loglik)
  
}

