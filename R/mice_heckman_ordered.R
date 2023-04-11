########################################################################
########################################################################
##
## MI function for imputing MNAR ordinal-scaled data 
## (usable with R package mice)
##
## @author Angelina Hammon
##
## 29.09.2021
##
########################################################################
########################################################################



# ordered probit model with sample selection:
#' @export
mice.impute.ordinalMNAR <- function(y, ry, x, type, draw=T, excl="",...) { 
  
  nam <- colnames(x) 
  
  pos.exc <- grep(excl,nam)   
  
  ind_s <- length(nam)+1                
  ind_y <- length(nam)-length(pos.exc)
  
  ry <- as.numeric(ry)
  data <- data.frame(x,y,ry)
  data$y <- as.factor(data$y)
  colnames(data) <- c(colnames(x),"y","ry")
  
  wy <- !ry 
  

  ## Calculate one-step heckman model ##
  exc <- pos.exc
  
  sel <- paste(nam, collapse="+")
  sel <- as.formula(paste("ry","~",sel,sep=""))
  out <- paste(nam[-exc], collapse="+")
  out <- as.formula(paste("y","~",out,sep=""))
  
  XS <- data.matrix(cbind(rep(1,length(y)),x))
  XO <- data.matrix(x[,-exc])
  
  # estimating the model:
  heckit <- loglik_ordbivprob_est(data=data,yS=ry, yO=y, XO=XO, XS=XS, kS=ind_s, kO=ind_y,sel = sel, out=out) 
  
  
  ## Get parameter estimates of the heckman model ## 
  phi <- as.matrix(Matrix::nearPD(MASS::ginv(heckit$hessian))$mat)
  q_hat <- heckit$estimates
  
  ## Draw ALL parameter values ##
  if(draw==T){
    q_star <- mvtnorm::rmvnorm(1,q_hat,phi,method = "chol")
    beta_s_star <- q_star[1:ind_s]          
    beta_star <- q_star[(ind_s+1):(ind_s+ind_y)] 
    rho_star <- tanh(q_star[ind_s+ind_y+1])
    z_star <-  q_star[(ind_s+ind_y+2):length(q_star)]
    for(i in 2:length(z_star)){
      z_star[i] <- exp(z_star[i])+z_star[i-1]
    }
  }else{ ## just draw betas ##
    q_star <- mvtnorm::rmvnorm(1,q_hat,phi,method = "chol")
    beta_s_star <- q_star[1:ind_s]          
    beta_star <- q_star[(ind_s+1):(ind_s+ind_y)] 
    rho_star <- tanh(q_hat[ind_s+ind_y+1])
    z_star <-  q_hat[(ind_s+ind_y+2):length(q_hat)]
    for(i in 2:length(z_star)){
      z_star[i] <- exp(z_star[i])+z_star[i-1]
    }
  }
  
  ## Calculate prob based on the new parameter values ##
  lin.predS <- XS %*% beta_s_star
  lin.predO <- XO %*% beta_star 
  
  lin.predS <- as.numeric(lin.predS)[!ry]
  lin.predO <- as.numeric(lin.predO)[!ry]
  
  fy <- as.factor(y)
  nc <- length(levels(fy))

  z_star <- c(-Inf,z_star,+Inf)
  
  p <- sapply(1:nc, function(j){ sapply(1:sum(!ry),function(i){ 
    (VGAM::pbinorm(-lin.predS[i],(z_star[j+1]-lin.predO[i]), cov12 = rho_star) - 
     VGAM::pbinorm(-lin.predS[i],(z_star[j]-lin.predO[i]), cov12 = rho_star))/pnorm(-lin.predS[i]) })
  })
  
  
  ## Draw Ymis based on new probability ##
  
  un <- rep(runif(sum(wy)), each = nc)
  
  draws <- un > apply(p, 1, cumsum)
  idx <- 1 + apply(draws, 2, sum)
  
  vec <- levels(fy)[idx]
  
  return(vec)
}

mice.impute.ordinalMNAR <- compiler::cmpfun(mice.impute.ordinalMNAR)

#' Internal function
#'
#' This function is a wrapper for estimating the ordered probit model with sample selection. Calculates starting values and conducts numerical optimization. The other two functions calculate the log likelihood and analytic gradients.
#'
#' @keywords internal
#' 
#' @importFrom MASS polr
# estimation function: 
loglik_ordbivprob_est <- function(data, yS, yO, XO, XS, kS, kO, sel, out, ...) {
  
  r_start <- glm(sel, data=data,family=binomial(link="probit"))
  
  y_start <- MASS::polr(out, data=data[which(as.numeric(yS)==1),], method="probit")
  
  z_start <- y_start$zeta
  
  z_trans <- z_start

  # transform them: 
  for(i in 2:length(z_start)){
    z_trans[i] <- log(z_start[i]-z_start[i-1])
  }
  
  start <- c(coef(r_start),coef(y_start),atanh(0.5),z_trans)
  
  opt <- stats::optim(par=start, fn=loglik_ordbivprob, gr=NULL, method=c("BFGS"), yS, yO, XO, XS, kS, kO)
  
  hess <- maxLik::numericHessian(loglik_ordbivprob, grad=NULL, t0=opt$par, yS=yS, yO=yO, XO=XO, XS=XS, kS=kS, kO=kO)
  
  res <- list(estimates=opt$par,hessian=hess)
  
  return(res)     
}


#' @rdname loglik_ordbivprob_est
#' @keywords internal
# likelihood function:
loglik_ordbivprob <- function(param, yS, yO, XO, XS, kS, kO) {
  
  betaS    <- param[1:kS]
  betaO    <- param[(kS+1):(kS+kO)]
  rho      <- tanh(param[kS+kO+1])
  z <- param[(kS+kO+2):length(param)]
  
  # retransform:
  for(i in 2:length(z)){
    z[i] <- exp(z[i])+z[i-1]
  }
  
  # get categories of yO:
  if(is.factor(yO)==T) {
    cat <- levels(yO)
  }else{
    cat <- unique(yO)[order(unique(yO))]
  }
  
  
  # linear predictors:
  lin.predS <- XS %*% betaS
  lin.predO <- XO %*% betaO
  
  
  # Probabilities entering the likelihood function:
  p0 <- ifelse(yS==0,pnorm(-lin.predS),0) ## P(yS=0)
  
  z_aux <- c(-Inf,z,Inf)
  
  for(i in 1:(length(z)+1)) {
    assign(paste0("p1",i-1),VGAM::pbinorm(lin.predS,z_aux[i+1]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[i]-lin.predO, cov12 = -rho))
  }  ## P(yS=1 & yO=h)
  
  p0 <- ifelse(p0==0,1.1111e-300,p0)
  
  for(i in 1:(length(z)+1)) {
    assign(paste0("p1",i-1),eval(parse(text=paste("ifelse(",paste0("p1",i-1),"==0,1.1111e-300,",paste0("p1",i-1),")"))))
  }

  # generate expressions for the respective probabilities:
  aux <- rep("ifelse(yS==",length(z)+1)
  
  for(i in 1:(length(z)+1)) {
    if(i==1){
      aux[i] <- paste0(aux[i],"0",",p0,")
    } 
    
    if(i > 1 & i < (length(z)+1)) {
      aux[i] <- paste0(aux[i],"1 & yO==",cat[i-1],",",paste0("p1",i-2),",")
    }
    
    if(i==(length(z)+1)){
      aux[i] <- paste0(aux[i],"1 & yO==",cat[i-1],",",paste0("p1",i-2),",",paste0("p1",i-1),paste0(rep(")",length(cat)),collapse = ""))
    }
  }
  
  eq <- paste(aux,collapse="")
  # evaluate the created expressions to get log-Likelihood function for every unit i:
  lik <- eval(parse(text=eq))
  
  suppressWarnings(loglik <- ifelse(lik < 1e-08, -1e+100, log(lik)))
  
  # total log likelihood:
  return(-sum(loglik))
  
}


#' @rdname loglik_ordbivprob_est
#' @keywords internal
# analytic gradients: 
loglik_ordbivprob_grad <- function(param, yS, yO, XO, XS, kS, kO) {
  
  betaS    <- param[1:kS]
  betaO    <- param[(kS+1):(kS+kO)]
  rho      <- tanh(param[kS+kO+1])
  z <- param[(kS+kO+2):length(param)]
  z_trans <- z
  # retransform:
  for(i in 2:length(z)){
    z[i] <- exp(z[i])+z[i-1]
  }
  
  # linear predictors:
  lin.predS <- XS %*% betaS
  lin.predO <- XO %*% betaO
  
  # get categories of yO:
  if(is.factor(yO)==T) {
    cat <- levels(yO)
  }else{
    cat <- unique(yO)[order(unique(yO))]
  }
  
  ### Gradients ###
  
  grad <- vector("numeric",length=length(param))
  
  r <- sqrt(1-rho^2)
  
  ## betaS:
  dp0 <- apply(XS,2, function(x) (-x) * (dnorm(-(lin.predS)))/(pnorm(-lin.predS))) ## P(yS=0)
  
  z_aux <- c(-Inf,z,Inf)
  for(i in 1:(length(z)+1)) {
    assign(paste0("dp1",i-1), apply(XS,2, function(x) ((x* dnorm(lin.predS) * pnorm(((z_aux[i+1]-lin.predO)+rho*lin.predS)/r)) - (x* dnorm(lin.predS) * pnorm(((z_aux[i]-lin.predO)+rho*lin.predS)/r)))/
                                      (VGAM::pbinorm(lin.predS,z_aux[i+1]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[i]-lin.predO, cov12 = -rho))))
  }  ## P(yS=1 & yO=h)
  
  dp0 <- ifelse(dp0==0,1.1111e-15,dp0)
  
  for(i in 1:(length(z)+1)) {
    assign(paste0("dp1",i-1),eval(parse(text=paste("ifelse(",paste0("dp1",i-1),"==0,1.1111e-15,",paste0("dp1",i-1),")"))))
  }
  
  aux <- rep("(yS==",length(z)+2)
  
  for(i in 1:(length(z)+2)) {
    if(i==1){
      aux[i] <- paste0(aux[i],"0)","*dp0")
    } 
    
    if(i > 1) {
      aux[i] <- paste0(aux[i],"1 & yO==",cat[i-1],")",paste0("*dp1",i-2))
    }
  }
  
  eq <- paste(aux,collapse="+")
  # execute expressions:
  dprod <- eval(parse(text=eq))
  
  dbetaS <- colSums(dprod)      
  grad[1:kS] <- dbetaS
  
  
  ## betaO:
  z_aux <- c(-Inf,z,Inf)
  for(i in 1:(length(z)+1)) {
    assign(paste0("dp1",i-1), apply(XO,2, function(x) (1/(VGAM::pbinorm(lin.predS,z_aux[i+1]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[i]-lin.predO, cov12 = -rho)))*
                                      ((-x* dnorm((z_aux[i+1]-lin.predO)) * pnorm((lin.predS+rho*((z_aux[i+1]-lin.predO)))/r)) - (-x* dnorm((z_aux[i]-lin.predO)) * pnorm((lin.predS+rho*((z_aux[i]-lin.predO)))/r)))))
  }  ## P(yS=1 & yO=h)

  for(i in 1:(length(z)+1)) {
    assign(paste0("dp1",i-1),eval(parse(text=paste("ifelse(",paste0("dp1",i-1),"==0,1.1111e-15,",paste0("dp1",i-1),")"))))
  }
  
  
  aux <- rep("(yS==",length(z)+1)
  
  for(i in 1:(length(z)+1)) {
      aux[i] <- paste0(aux[i],"1 & yO==",cat[i],")",paste0("*dp1",i-1))
  }
  
  eq <- paste(aux,collapse="+")
  # execute expressions:
  dprod <- eval(parse(text=eq))
  
  dbetaO <- colSums(dprod)    
  
  grad[(kS+1):(kS+kO)] <- dbetaO
  
  
  ## rho: 
  z_aux <- c(-Inf,z,Inf)
  for(i in 1:(length(z)+1)) {
    assign(paste0("dr1",i-1), (1/(VGAM::pbinorm(lin.predS,z_aux[i+1]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[i]-lin.predO, cov12 = -rho))) * 
                              (((-1)*(1-rho^2)*VGAM::dbinorm(lin.predS,z_aux[i+1]-lin.predO, mean1 = 0, mean2 = 0, cov12 = -rho)) - ((-1)*(1-rho^2)*VGAM::dbinorm(lin.predS,z_aux[i]-lin.predO, mean1 = 0, mean2 = 0, cov12 = -rho))))
  }  ## P(yS=1 & yO=h)
  
  for(i in 1:(length(z)+1)) {
    assign(paste0("dr1",i-1),eval(parse(text=paste("ifelse(",paste0("dr1",i-1),"==0,1.1111e-15,",paste0("dr1",i-1),")"))))
  }
  
  aux <- rep("(yS==",length(z)+1)
  
  for(i in 1:(length(z)+1)) {
    aux[i] <- paste0(aux[i],"1 & yO==",cat[i],")",paste0("*dr1",i-1))
  }
  
  eq <- paste(aux,collapse="+")
  # execute expressions:
  dprod <- eval(parse(text=eq))
  
  drho <- sum(dprod)
  
  grad[(kS+kO+1)] <- drho

  
  # threshold parameters: 
  
  z_aux <- c(-Inf,z,Inf)
  for(j in 1:length(z)){
    
      if(j==1){
        assign(paste0("dz",j,j), (1/(VGAM::pbinorm(lin.predS,z_aux[j+1]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[j]-lin.predO, cov12 = -rho))) * 
                 (dnorm(z_aux[j+1]-lin.predO) * pnorm((lin.predS+rho*((z_aux[j+1]-lin.predO)))/r)))
        
        assign(paste0("dz",j,j+1), (1/(VGAM::pbinorm(lin.predS,z_aux[j+2]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[j+1]-lin.predO, cov12 = -rho))) * 
                 (-(dnorm(z_aux[j+1]-lin.predO) * pnorm((lin.predS+rho*((z_aux[j+1]-lin.predO)))/r))))
      }
    
    else{
      assign(paste0("dz",j,j), (1/(VGAM::pbinorm(lin.predS,z_aux[j+1]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[j]-lin.predO, cov12 = -rho))) * 
               (dnorm(z_aux[j+1]-lin.predO) * pnorm((lin.predS+rho*((z_aux[j+1]-lin.predO)))/r)) * exp(z_trans[j]))
      
      assign(paste0("dz",j,j+1), (1/(VGAM::pbinorm(lin.predS,z_aux[j+2]-lin.predO, cov12 = -rho) - VGAM::pbinorm(lin.predS,z_aux[j+1]-lin.predO, cov12 = -rho))) * 
               (-(dnorm(z_aux[j+1]-lin.predO) * pnorm((lin.predS+rho*((z_aux[j+1]-lin.predO)))/r))) * exp(z_trans[j]))
    }
    
  }

  
  for(j in 1:length(z)) {
    assign(paste0("dz",j,j),eval(parse(text=paste("ifelse(",paste0("dz",j,j),"==0,1.1111e-15,",paste0("dz",j,j),")"))))
    assign(paste0("dz",j,j+1),eval(parse(text=paste("ifelse(",paste0("dz",j,j+1),"==0,1.1111e-15,",paste0("dz",j,j+1),")"))))
  }
  
  
  
  aux_list <- vector("list",length(z))
  
  for(j in 1:length(z)){
    aux_list[[j]] <- rep("(yS==",length(z))
  }

  for(j in 1:length(z)) {
    aux_list[[j]][1] <- paste0(aux_list[[j]][1],"1 & yO==",cat[j],")",paste0("*dz",j,j)) 
    aux_list[[j]][2] <- paste0(aux_list[[j]][2],"1 & yO==",cat[j+1],")",paste0("*dz",j,j+1)) 
  }
  

  dprod <- sapply(aux_list,function(x) {
    eq <- paste(x,collapse="+")
    dprod <- eval(parse(text=eq))
  })
  
  dz <- colSums(dprod)    
  
  grad[(kS+kO+2):length(param)] <- dz
  
  
  # return gradient vector:
  return(-grad)
  
}