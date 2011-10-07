% --- Source file: man/types.forest.Rd ---
\name{types.forest}
\alias{types.forest}
\title{Forest Plot for subtype analysis}
\description{
Forest Plot for subtype analysis.
}
\usage{
types.forest(rlist, snp.var, level=0.05, p.adj=TRUE, digits=2)
}
\arguments{
  \item{rlist}{The list of results returned by \code{\link{h.types}}. SNPs other than \code{snp.var} are ignored}
  \item{snp.var}{A character string giving the name of the SNP variable to be plotted. No default.} 
  \item{level}{Level for confidence intervals. Default is \code{0.05} for \code{95\%} CI}
  \item{p.adj}{Logical. Whether to report Bonferroni adjusted p-values for each individual subtype. Default is TRUE.}  
  \item{digits}{Number of digit to round off, odds ratios in the plot. Also number of significant digits to retain for displaying p-values.}
}
\value{
Forest Plot for a SNP showing case-control odds ratios, confidence intervals and p-values for individual disease subtypes, 
an estimate of an overall case-control odds ratio and confidence interval which includes all disease subtypes,  
and an estimate of the odds ratio and confidence intervals associated with the identified best subset of disease subtypes.
}

\examples{
 # Use the example data
 data(ex_types, package="ASSET")

 # Display the first 10 rows of the data and a table of the subtypes
 data[1:10, ]
 table(data[, "TYPE"])
 
 # Define the input arguments to h.types. 
 adj.vars <- c("CENTER_1", "CENTER_2", "CENTER_3")
 types    <- paste("SUBTYPE_", 1:5, sep="")

 # SUBTYPE_0 will denote the controls
 res <- h.types(data, "TYPE", "SNP_1", adj.vars, types, "SUBTYPE_0", logit=TRUE, 
        test.type="Score", meth.pval="DLM")

 
 types.forest(res, "SNP_1")

}
