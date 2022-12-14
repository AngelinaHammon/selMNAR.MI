# selMNAR.MI
R package for multiple imputation of single-level and multilevel binary and ordinal variables that are supposed to be MNAR

The package provides imputation methods based on selection models to multiply impute binary and ordinal-scaled MNAR data in the context of Fully Conditional Specification. The models are also able to reflect hierarchical data structures by incorporating random intercept terms. The imputation functions 'mice.impute.heckman1step', 'mice.impute.2l.heckman1step_ghq', 'mice.impute.2l.heckman1step_aghq','mice.impute.heckman1step_ord','mice.impute.2l.heckman1step_ord_ghq' and 'mice.impute.2l.heckman1step_ord_aghq' are implemented so that they can directly be used within the imputation function of the 'mice' package to impute multivariate missing data.
