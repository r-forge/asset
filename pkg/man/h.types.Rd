% --- Source file: man/h.types.Rd ---
\name{h.types}
\alias{h.types}
\title{Heterogeneous Subtype Analysis}
\description{
Subset-based Case-control/Case-complement analysis of heterogeneous disease subtypes.
}
\usage{
h.types(dat, response.var, snp.vars, adj.vars, types.lab, cntl.lab, subset=NULL, pool=NULL, side=2, logit=FALSE
	, meth = "Score", zmax.args = NULL, meth.pval = c("DLM", "IS", "B"), pval.args = NULL)
}
\arguments{
  \item{dat}{A data frame containing individual level phenotype, genotype and covariate data. No default.}
  \item{response.var}{A character string giving the name of the response variable containing disease status/subtype information. No default.}
  \item{snp.vars}{A character vector giving the names or positions of the SNP variables to be analyzed. No default.}
  \item{adj.vars}{A character vector giving the names of the columns containing covariates to be adjusted. No default. NULL indicates only intercept.}
  \item{types.lab}{A character vector giving the names/identifier of the disease subtypes in response.var. No default.}
  \item{cntl.lab}{A single character string giving the name/identifier of controls in response.var. No default.}
  \item{subset}{Logical vector. Subset of rows of the data frame to exclude. Default is NULL, all rows are used.}
  \item{pool}{TRUE indicates case-complement analysis, FALSE indicated case-control analysis, NULL indicates both. Default is NULL.}
  \item{side}{Either 1 or 2. For two-tailed tests (where absolute values of Z-scores are maximzed), side should be 2. For one-tailed tests 
  side should be 1 (positive tail is assumed). Default is 2, ignored when search is 2.}
  \item{logit}{If TRUE, overal logistic regression is done for each SNP and results returned. Default is FALSE.}
  \item{meth}{ Either "Score" or "Wald". Score test is much faster when large blocks of SNPs are analyzed at a time. Default is "Score."}
  \item{zmax.args}{Logical. Extra arguments such as \code{sub.def} and \code{subargs.list} to be
  passed to \code{\link{z.max}} as a named list if non-default values are desired}.
  \item{meth.pval}{The method of pvalue computation. Currently the options are DLM (Discrete Local
  Maximum), IS(Importance Sampling) and B (Bonferroni).}
  \item{pval.args}{Logical. Extra arguments such as \code{sub.def} and \code{subargs.list} to be
  passed to \code{\link{p.dlm}} or \code{\link{p.tube}} as a named list if non-default values
  are desired}.  
}
\value{
   A list containing 3 component lists named:
    
   (1) "Overall.Logistic" (standard logistic regression with all cases): 
	   This list is non-null when \code{logit} is TRUE and contains 3 vectors named (pval, beta, sd) of length same as snp.vars.
   
   (2) "Subset.Case.Control" (subset search with case-control analysis):
	   This list is non-null when \code{pool} is NULL or FALSE and contains, 3 vectors named (pval, beta, sd) of length same as snp.vars 
	   and a logical matrix named "pheno" with one row for each snp and one column for each phenotype. For a particular SNP and phenotype this
	   matrix has "TRUE" if this phenotype was in the selected subset, by case-control search.
	   
   (3) "Subset.Case.Complement" (subset search with case-complment analysis):
       This list is non-null when \code{pool} is NULL or FALSE and contains, 3 vectors named (pval, beta, sd) of length same as snp.vars 
	   and a logical matrix named "pheno" with one row for each snp and one column for each phenotype. For a particular SNP and phenotype this
	   matrix has "TRUE" if this phenotype was in the selected subset, by case-control search.
}
\details{
The output standard errors are approximate (based on inverting DLM pvalues) and are used for constructing confidence intervals 
in \code{\link{h.summary}} and \code{\link{types.forest}}.
}

