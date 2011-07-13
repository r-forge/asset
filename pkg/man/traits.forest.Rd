% --- Source file: man/traits.forest.Rd ---
\name{traits.forest}
\alias{traits.forest}
\title{Forest Plot for meta-analysis of heterogenerous traits.}
\description{
Forest Plot for meta-analysis of heterogenerous traits.
}
\usage{
traits.forest(snp.var, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor = NULL, rlist = NULL
		, level = 0.05, p.adj = TRUE, digits = 2)
}
\arguments{
  \item{snp.var}{A character string giving the name of the SNP variable to be plotted. No default.}
  \item{traits.lab}{A character vector giving the names/identifier of the \code{k} traits being analyzed. No default.}
  \item{beta.hat}{A matrix of dimension length(\code{snp.vars}) X (\code{k}) (or a vector of length \code{k} when 1 SNP is passed). Each row gives the coefficients obtained
  from analysis of that SNP using each of the (k) traits individually. No default.}
  \item{sigma.hat}{A vector or matrix of same dimension as beta.hat, giving the corresponding standard errors. No default.}
  \item{ncase}{The number of cases in each of the \code{k} studies. This can be same for each SNP, in which case \code{ncase} is a vector of length \code{k}.
  Alternatively if number of non-missing cases analyzed for each SNP is known, ncase can be a length(\code{snp.vars})X (\code{k}) matrix. No default.}
  \item{ncntl}{Same as ncase (above) for controls. No default.}
  \item{cor}{Either a \code{k} X \code{k} matrix of inter-study correlations or a list containing four case/control overlap matrices 
  named \code{N11}, \code{N00}, \code{N10} and \code{N01}. See details in \code{\link{h.traits}}. Default is NULL, so that studies are assumed to be independent.}
  \item{rlist}{The list of results returned by \code{\link{h.traits}} if already available. SNPs other than \code{snp.var} are ignored} 
  \item{level}{Level for confidence intervals. Default is \code{0.05} for \code{95\%} CI}
  \item{p.adj}{Logical. Whether to report Bonferroni adjusted p-values for each individual subtype. Default is TRUE.}  
  \item{digits}{Number of digit to round off, odds ratios in the plot. Also number of significant digits to retain for displaying p-values.}
}
\value{
Forest Plot for individual, subset-search and overall results of \code{snp.var}.
}
