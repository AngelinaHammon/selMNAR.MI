%% File Name: mice_heckman_binary.Rd
%% File Version: 1.0

\name{mice.impute.binaryMNAR}
\alias{mice.impute.binaryMNAR}

\title{
  Multiple imputation of binary MNAR data 
}

\description{
  The function \code{mice.impute.binaryMNAR}
  imputes a binary variable that is supposed to follow an MNAR mechanism using a censored bivariate probit model. 
  To estimate the model Maximum Likelihood is used. 
  See Hammon & Zinn (2020) for details.
}

\usage{
  mice.impute.binaryMNAR(y, ry, x, type, draw=T, excl="", ...)
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
  \item{type}{Type of predictor variable where \code{1} indicates that variable is used as predictor.
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
  Hammon, A., & Zinn, S. (2020). \emph{Multiple imputation of binary multilevel missing not at random data}.
  J. R. Stat. Soc. C, 69: 547-564.
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
    # EXAMPLE Imputation of binary single-level MNAR data
    #############################################################################
    
    #--- simulate data and missing values
    set.seed(987)
    
    n <- 2000       # total sample size
    
    # covariates:
    x1 <- rnorm(n,0,0.3)
    x2 <- rnorm(n,0,1)
    x3 <- rnorm(n,0,4)      # exclusion criterion 
    
    # simulate error terms:
    rho <- 0.6      # correlation between error terms 
    vc <- diag(2)
    vc[2,1] <- vc[1,2] <- rho
    eps <- mvtnorm::rmvnorm(n,rep(0,2),vc)
    eps_o <- eps[,1]
    eps_s <- eps[,2]
    
    # binary outcome variable to be imputed:
    y <- 0.25+1*x1+0.5*x2+eps_o
    y <- as.factor(ifelse(y<0,0,1))
    
    # missing-data indicator of y:
    r <- 0.5+1.5*x1-0.25*x2+0.1*x3+eps_s
    r <- as.factor(ifelse(r<0,0,1))

    #* create some missings
    y[r==0] <- NA
    data <- data.frame(y,r,x1,x2,x3)
    
    #--- prepare imputation in mice
    
    #* predictor matrix
    pred_MNAR <- make.predictorMatrix(data)
    pred_MNAR["y","r"] <- 0
    
    #* imputation methods
    meth_MNAR <- make.method(data)
    meth_MNAR["y"] <- "binaryMNAR"
    
    #** multiple imputation 
    imp_MNAR <- mice::mice(data,m=5,maxit=1,method=meth_MNAR,pred=pred_MNAR,print=F,
                            excl="x3",seed=1234)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{mice imputation method}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line