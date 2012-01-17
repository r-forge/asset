% --- Source file: data.Rd ---
\name{data}
\alias{data}
\alias{N00}
\alias{N01}
\alias{N10}
\alias{N11}

%\docType{package}
\title{
  Data for the h.traits example
}
\description{
 Data for \code{\link{h.traits}}
}
\details{ The object data contains estimated log odds-ratios and their standard errors
  for 5 SNPs and 6 traits. The matrices N00, N01, N10, and N11 are the case-control 
  overlap matrices.
}
\examples{

 data(ex_trait, package="ASSET")

 # Display the data, and case/control overlap matrices
 data
 N00
 N11
 N01
 N10

}
\keyword{data}
\seealso{h.traits}



