# selMNAR.MI
R package for multiple imputation of single-level and multilevel binary and ordinal variables that are supposed to be MNAR

The package provides imputation methods based on selection models to multiply impute binary and ordinal-scaled MNAR data in the context of Fully Conditional Specification. The models are also able to reflect hierarchical data structures by incorporating random intercept terms. The imputation functions 'mice.impute.binaryMNAR', 'mice.impute.2l.binaryMNAR', 'mice.impute.ordinalMNAR' and 'mice.impute.2l.ordinalMNAR' are implemented so that they can directly be used within the imputation function of the 'mice' package to impute multivariate missing data.
