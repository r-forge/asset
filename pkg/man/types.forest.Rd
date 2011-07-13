% --- Source file: man/types.forest.Rd ---
\name{types.forest}
\alias{types.forest}
\title{Forest Plot for subtype analysis}
\description{
Forest Plot for subtype analysis.
}
\usage{
types.forest(dat, snp.var, response.var, adj.vars, types.lab, cntl.lab, subset=NULL, rlist = NULL
	, level = 0.05, p.adj = TRUE, digits = 2)
}
\arguments{
  \item{dat}{A data frame containing individual level phenotype, genotype and covariate data. No default.}
  \item{snp.var}{A character string giving the name of the SNP variable to be plotted. No default.}
  \item{response.var}{A character string giving the name of the response variable containing disease status/subtype information. No default.}
  \item{adj.vars}{A character vector giving the names of the columns containing covariates to be adjusted. No default. NULL indicates only intercept.}
  \item{types.lab}{A character vector giving the names/identifier of the disease subtypes in response.var. No default.}
  \item{cntl.lab}{A single character string giving the name/identifier of controls in \code{response.var}. No default.}
  \item{subset}{Logical vector. Subset of rows of the data frame to exclude. Default is NULL, all rows are used.}
  \item{rlist}{The list of results returned by \code{\link{h.types}} if already available. SNPs other than \code{snp.var} are ignored} 
  \item{level}{Level for confidence intervals. Default is \code{0.05} for \code{95\%} CI}
  \item{p.adj}{Logical. Whether to report Bonferroni adjusted p-values for each individual subtype. Default is TRUE.}  
  \item{digits}{Number of digit to round off, odds ratios in the plot. Also number of significant digits to retain for displaying p-values.}
}
\value{
Forest Plot for individual, subset-search and overall results of \code{snp.var}.
}
