% --- Source file: ASSET.Rd ---
\name{ASSET}
\alias{ASSET}
\docType{package}
\title{
Association analysis for SubSETs
}
\description{
This package is for subset-based association analysis of heterogeneous but possibly related traits measured
in different studies or the same study, where each individual can have only one phenotype (say a cancer type
or a histological subtype). Thus, this package is not appropriate for multivatiate phenotype analysis where
the different traits may be measured on the same person.
}
\details{
The main functions for analysis of heterogenous traits (from separate studies) and subtype analysis within
a single case-control study are \code{\link{h.traits}} and \code{\link{h.types}} respectively. The helper functions
\code{\link{z.max}} and \code{\link{p.dlm}} are generic functions called internally for specific definitions of
test statistics and correlations and without any subset restrictions. These functions are provided to facilitate customization
to other application settings that involve z-score maximization over subsets. For example ordering restrictions can be imposed
on subsets by passing a function argument \code{sub.def} to \code{z.max} and \code{p.dlm}.
The functions \code{\link{h.summary}}, \code{\link{types.forest}} and \code{\link{traits.forest}} are also provided for 
a quick summary and visual representation of subset specific effects and their confidence intervals.
}
\author{Samsiddhi Bhattacharjee, Nilanjan Chatterjee and William Wheeler <wheelerb@imsweb.com>}
\keyword{package}