% --- Source file: man/h.summary.Rd ---
\name{h.summary}
\alias{h.summary}
\title{Summaries for subset-search.}
\description{
Summarise results of subset-search.
}
\usage{
h.summary(rlist, level = 0.05, digits = 3)
}
\arguments{
  \item{rlist}{List returned by \code{\link{h.traits}} or \code{\link{h.types}}}
  \item{level}{Level for confidence intervals. Default is \code{0.05} for \code{95\%} CI}
  \item{digits}{Number of significant digits to retain in odds ratio-s and pvalues for summary table}
}
\value{
A list of data frames, one for each non-null component of the rlist. Each data frame gives summary of 
phenotypes selected, pvalues, odds ratio-s and CI-s.
}
\details{ 
Return data frames have as many rows as SNPs in the return list. 
}
