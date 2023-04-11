%% File Name: mice_heckman_binary_ordered_re.Rd
%% File Version: 1.0

\name{mice.impute.2l.ordinalMNAR}
\alias{mice.impute.2l.ordinalMNAR}

\title{
  Multiple imputation of ordered multilevel MNAR data using GHQ
}

\description{
  The function \code{mice.impute.2l.ordinalMNAR} imputes an ordinal-scaled clustered variable that is supposed to follow an MNAR mechanism using an ordered probit model with sample selection and random effects. 
  To estimate the model Maximum Likelihood and quadrature techniques are used. \code{mice.impute.2l.ordinalMNAR} uses standard Gauss-Hermite quadrature. See Hammon (2022) for details.
}

\usage{
  mice.impute.2l.ordinalMNAR(y, ry, x, type, QP=rep(10,2), draw=T, excl="", ...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    Incomplete data vector of length \code{n}
  }
  \item{ry}{
    Vector of missing data pattern (\code{FALSE} -- missing,
                                     \code{TRUE} -- observed)
  }
  \item{x}{
    Matrix (\code{n} x \code{p}) of complete predictors.
  }
  \item{type}{Type of predictor variable. The cluster identifier
    has type \code{-2} and fixed effects predictors 
    type \code{1}.}
  \item{QP}{
    Number of quadrature points that should be used for the two random intercepts (needs to be a vector of length \code{2}).
    Default number of quadrature points is \code{10}.
  }
  \item{draw}{
    Set \code{TRUE} if all parameters should be drawn from approximate normal posterior.
    Set \code{FALSE} if only betas should be drawn from approximate normal posterior.
    Default is to draw new parameter candidates for all parameters. 
  }
  \item{excl}{
    Provide variable name of exclusion criterion as a string. 
  }
  \item{\dots}{
    Further arguments to be passed
  }
}

%\details{
  %%  ~~ If necessary, more details than the description above ~~
    %}

\value{
  A vector of length \code{nmis=sum(!ry)} with imputed values.
}

\references{
  Hammon, A. (2022). \emph{Multiple Imputation of Ordinal Missing Not at Random Data}.
  AStA Adv Stat Anal.
}

\author{
  Angelina Hammon
  }

%\note{
  %%  ~~further notes~~
    %}

%% ~Make other sections like Warning with \section{Warning }{....} ~
  
  

\examples{
  \dontrun{
    #############################################################################
    # EXAMPLE Imputation of ordered multilevel MNAR data
    #############################################################################
    
    #--- simulate missing values
    set.seed(987)
    
    n <- 2500      # total sample size
    
    m <- 20        # number of clusters
    nj <- n/m      # units per cluster
    
    # covariates:
    x1 <- rnorm(n,0,0.3)
    x2 <- rnorm(n,0,0.8)
    x3 <- rnorm(n,0,4)     # exclusion criterion 
    
    # simulate error terms:
    rho <- 0.6     # correlation between error terms
    vc <- diag(2)
    vc[2,1] <- vc[1,2] <- rho
    eps <- mvtnorm::rmvnorm(n,rep(0,2),vc)
    eps_o <- eps[,1]     
    eps_s <- eps[,2]
  
    # simulate random intercepts:
    tau <- 0.5     # correlation between random intercepts
    vc_re <- matrix(1,2,2)
    vc_re[1,1] <- 0.5
    vc_re[2,2] <- 0.9
    vc_re[2,1] <- vc_re[1,2] <- sqrt(vc_re[1,1])*sqrt(vc_re[2,2])*tau
    
    alpha <- mvtnorm::rmvnorm(m,rep(0,2),vc_re) 
    alpha_s <- rep(alpha[,1], each=nj)  
    alpha_o <- rep(alpha[,2], each=nj)  
    
    # ordinal outcome variable to be imputed:
    y_star <- 1*x1+0.5*x2+eps_o+alpha_o
    z <- c(-0.75,0.5)    # threshold parameters
    z1 <- z[1]
    z2 <- z[2]
    y <- as.factor(ifelse(y_star <= z1 ,0,
                         ifelse(y_star > z1 & y_star <= z2 , 1 ,
                                ifelse(y_star > z2 , 2 ,NA))))
    
    # missing-data indicator of y:          
    r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s+alpha_s
    r <- ifelse(r<0,0,1)
    
    group <- rep(1:m, each=nj)    # cluster IDs
          
    #* create some missings
    y[r==0] <- NA
    data <- data.frame(y,r,x1,x2,x3,group)
    
    #--- prepare imputation in mice

    #* predictor matrix
    pred_MNAR <- make.predictorMatrix(data)
    pred_MNAR["y","r"] <- 0
    pred_MNAR["y","group"] <- -2
    
    #* imputation methods
    meth_MNAR <- make.method(data)
    meth_MNAR["y"] <- "2l.ordinalMNAR"

    #** multiple imputation 
    imp_MNAR_ghq <- mice::mice(data,m=5,maxit=1,method=meth_MNAR,pred=pred_MNAR,print=F,
                               excl="x3",seed=1234)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{mice imputation method}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line