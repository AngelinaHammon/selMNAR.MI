%% File Name: mice_heckman_ordered.Rd
%% File Version: 1.0

\name{mice.impute.ordinalMNAR}
\alias{mice.impute.ordinalMNAR}


\title{
  Multiple imputation of ordered MNAR data 
}

\description{
  The function \code{mice.impute.ordinalMNAR}
  imputes an ordinal-scaled variable that is supposed to follow an MNAR mechanism using an ordered probit model with sample selection. 
  To estimate the model Maximum Likelihood is used. 
  See Hammon (2022) for details.
}

\usage{
  mice.impute.ordinalMNAR(y, ry, x, type, draw=T, excl="", ...)
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
    # EXAMPLE Imputation of ordered single-level MNAR data
    #############################################################################
    
    #--- simulate missing values
    set.seed(987)

    n <- 2000      # total sample size
    
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
    
    # ordinal outcome variable to be imputed:
    y_star <- 1*x1+0.5*x2+eps_o
    z <- c(-0.75,0.5)      # treshold parameters
    z1 <- z[1]
    z2 <- z[2]
    y <- as.factor( ifelse(y_star <= z1 ,0,
                          ifelse(y_star > z1 & y_star <= z2 , 1 ,
                                ifelse(y_star > z2 , 2 ,NA))))
    
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
    meth_MNAR["y"] <- "ordinalMNAR"
    
    #** multiple imputation 
    imp_MNAR <- mice::mice(data,m=5,maxit=1,method=meth_MNAR,pred=pred_MNAR,print=F,
                           excl="x3",seed=1234)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{mice imputation method}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line