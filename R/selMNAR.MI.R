#' selMNAR.MI: A package for multiple imputation of single-level and multilevel binary and ordinal variables that are supposed to be MNAR.
#'
#' The selMNAR.MI package provides imputation methods and functions based on selection models for binary and ordinal-scaled MNAR data 
#' that can be used within the 'mice' package and its FCS algorithm. There are versions for single-level and multilevel data which both 
#' are based on bivariate selection models. The multilevel models are estimated using quadrature techniques. 
#' The model for binary variables is described in Hammon & Zinn (2020) and the method for ordinal variables can be found in Hammon (2022). 
#' 
#'
#' @docType package
#' @name selMNAR.MI
#' @import mice 
#' @useDynLib selMNAR.MI, .registration=TRUE
NULL
#> NULL
