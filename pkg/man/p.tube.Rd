% --- Source file: man/p.tube.Rd ---
\name{p.tube}
\alias{p.tube}
\title{Importance sampling-based estimated p-value.}
\description{
Function to obtain "Importance Sampling"-based estimates of pvalues for z-scores maximized over subsets (of traits or subtypes), 
with possible restrictions and weights. Should not be called directly. See details.
}
\usage{
p.tube(t.vec, k, search, side, ncase, ncntl, pool, rmat, cor.numr, sizes = rep(1, k), nsamp = 50, sub.def=NULL, sub.args=NULL, wt.def=NULL, wt.args=NULL)
}
\arguments{
  \item{t.vec}{Numeric vector of (positive) points for which to calculate p-values, i.e. generall observed Z-max values. No default.}
  \item{k}{Integer (currently less than 30). The number of studies (traits) or subtypes being analyzed. No default.}  
  \item{search}{0, 1 or 2. Search option, with 0 indicating subtype analysis, 1 and 2 denote one-sided and two-sided subset-search. No default.}
  \item{side}{Either 1 or 2. For two-tailed tests (where absolute values of Z-scores are maximzed), side should be 2. For one-tailed tests 
  side should be 1 (positive tail is assumed). No default. Ignored when search is 2.}
  \item{ncase}{The number of cases in each of the \code{k} studies. This can be same for each SNP, in which case \code{ncase} is a vector of length \code{k}.
  Alternatively if number of non-missing cases analyzed for each SNP is known, ncase can be a length(\code{snp.vars})X (\code{k}) matrix. No default.}
  \item{ncntl}{Same as ncase (above) for controls. No default.}
  \item{pool}{TRUE indicates case-complement analysis, FALSE indicated case-control analysis. No default when \code{search is 0}, ignored otherwise.}
  \item{rmat}{A \code{k} X \code{k} matrix of inter-study correlations. Default is NULL, so that studies are assumed to be independent. Ignored when \code{search} is 0.}
  \item{cor.numr}{Logical. Whether to consider correlation in the numerator of the meta-analysis statistic. No default, ignored when \code{search} is 0.}
  \item{sizes}{Sizes of equivalence classes of traits. By default, no two traits or studies are equivalent. This argument is for internal use.}
  \item{nsamp}{Number of importance sampling replicates. Default is 50. See details.}   
  \item{sub.def}{A function to restrict subsets, e.g., order restrictions in subtype analysis. Should accept a subset (a logical vector of size k) as
  its first argument and should return TRUE if the subset satisfies restrictions and FALSE otherwise. Default is NULL implying all (2^k - 1) subsets
  are considered in the maximum.}
  \item{sub.args}{Other arguments to be passed to \code{sub.def} as list. Default is null (i.e. none).}
  \item{wt.def}{A function that gives weight of one subset with respect to another. Should accept two subsets as first two argumets and return
  a single positive weight. Default NULL. Currently this option is not implemented and the argument is ignored.}
  \item{wt.args}{Other arguments to be passed to \code{wt.def} as a list . Default NULL. Currently ignored.}
  }
\value{
A numeric vector of estimated p-values.
}
\details{ 
The function is vectorized to handle blocks of SNPs at a time. Currently weight options are ignored.
This is a helper function that is called internally by \code{\link{h.traits}} and \code{\link{h.types}} 
and should not be called directly. The arguments of this function that have defaults, can be customized using
the argument \code{pval.args} in \code{\link{h.traits}} and \code{\link{h.types}}.
Using higher value of \code{nsamp} such as 500, will give more accurate answers but can become extremely
slow particularly if \code{k} is high (~10). Generally, even 10 replicates can give reasonably accurate answers. Does not
depend on how small the p-value is.
}
