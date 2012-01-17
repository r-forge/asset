% --- Source file: man/h.summary.Rd ---
\name{h.summary}
\alias{h.summary}
\title{summary results from subset-search.}
\description{
The function produces summarise results from subset-based association analysis.
}
\usage{
h.summary(rlist, level = 0.05, digits = 3)
}
\arguments{
  \item{rlist}{List returned by \code{\link{h.traits}} or \code{\link{h.types}} }
  \item{level}{Level for confidence intervals. Default is \code{0.05} for \code{95\%} CI. The confidence intervals are obtained by inverting the 
corresponding multiple-testing adjusted p-values.}
  \item{digits}{Number of significant digits to retain in odds ratios and confidence intervals in the summary table}
}
\value{
A list of data frames, one for each of the methods specified the original call of the functions \code{\link{h.traits}} or \code{\link{h.types}}. Each row of a data frame corresponds to a SNP and 
the values include p-values for overall association (including component-wise p-values for two-sided search), names of phenotypes or disease subtypes included in
the best-subset, summary regression coefficients (e.g. log-odds-ratio for case-control studies) representing strength of association of a SNP with the
identified subset of traits/subtype  and corresponding confidence intervals.
}
\details{ 
Return data frames having as many rows as SNPs in the return list. 
}
\seealso{\code{\link{h.forestPlot}}, \code{\link{h.traits}}, \code{\link{h.types}} }