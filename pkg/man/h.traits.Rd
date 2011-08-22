% --- Source file: man/h.traits.Rd ---
\name{h.traits}
\alias{h.traits}
\title{Heterogeneous traits or studies}
\description{
One-sided or Two-sided Subset-based Meta-analysis of Heterogeneous Traits or Studies.
}
\usage{
h.traits(snp.vars, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=NULL, cor.numr=FALSE
, search=NULL, side=2, meta=FALSE, zmax.args=NULL, meth.pval=c("DLM","IS","B"), pval.args=NULL)
}
\arguments{
  \item{snp.vars}{A character vector giving the names or positions of the SNP variables to be analyzed. No default.}
  \item{traits.lab}{A character vector giving the names/identifier of the \code{k} traits being analyzed. No default.}
  \item{beta.hat}{A matrix of dimension length(\code{snp.vars}) X (\code{k}) (or a vector of length \code{k} when 1 SNP is passed). Each row gives the coefficients obtained
  from analysis of that SNP using each of the (k) traits individually. No default.}
  \item{sigma.hat}{A vector or matrix of same dimension as beta.hat, giving the corresponding standard errors. No default.}
  \item{ncase}{The number of cases in each of the \code{k} studies. This can be same for each SNP, in which case \code{ncase} is a vector of length \code{k}.
  Alternatively if number of non-missing cases analyzed for each SNP is known, ncase can be a length(\code{snp.vars})X (\code{k}) matrix. No default.}
  \item{ncntl}{Same as ncase (above) for controls. No default.}
  \item{cor}{Either a \code{k} X \code{k} matrix of inter-study correlations or a list containing four case/control overlap matrices 
  named \code{N11}, \code{N00}, \code{N10} and \code{N01}. See details. Default is NULL, so that studies are assumed to be independent.}
  \item{cor.numr}{Logical. Whether to use the correlation matrix in the numerator of the meta-analysis statistic
  to gain efficiency. Default is FALSE.}
  \item{search}{1, 2 or NULL. Search option with 1 indicating 1-sided subset search, 2 indicating 2-sided subset-search and NULL indicating both. 
  Default is NULL.}
  \item{side}{Either 1 or 2. For two-tailed tests (where absolute values of Z-scores are maximzed), side should be 2. For one-tailed tests 
  side should be 1 (positive tail is assumed). Default is 2, ignored when search is 2.}
  \item{meta}{Logical. Whether to compute and return standard meta-analysis results. Default is FALSE.}
  \item{zmax.args}{Logical. Extra arguments such as \code{sub.def} and \code{subargs.list} to be
  passed to \code{\link{z.max}} as a named list, if non default values are desired}.
  \item{meth.pval}{The method of pvalue computation. Currently the options are DLM (Discrete Local
  Maximum), IS(Importance Sampling) and B (Bonferroni).}
  \item{pval.args}{Logical. Extra arguments such as \code{sub.def} and \code{subargs.list} to be
  passed to \code{\link{p.dlm}} or \code{\link{p.tube}} as a named list, if non-default values
  are desired}.
}
\value{
   A list containing 3 component lists named:
    
   (1) "Meta" (standard meta-analysis with all traits/studies). 
	   This list is non-null when \code{meta} is TRUE and contains 3 vectors named (pval, beta, sd) of length same as snp.vars.
   
   (2) "Subset.1sided" (one-sided subset search):
	   This list is non-null when \code{search} is NULL or 1 and contains, 3 vectors named (pval, beta, sd) of length same as snp.vars 
	   and a logical matrix named "pheno" with one row for each snp and one column for each phenotype. For a particular SNP and phenotype this
	   matrix has "TRUE" if this phenotype was in the selected subset, by 1-sided search.
	   
   (3) Subset.2sided (two-sided subset search)
       This list is non-null when \code{search} is NULL or 2 and contains 7 vectors named (pval, pval.1, pval.2, beta.1, sd.1, beta.2, sd.2) of length
	   same as snp.vars and two matrices named "pheno.1" and "pheno.2" giving logical indicators of a phenotype being among the positively or negatively
	   associated subsets (respectively) as identified by 2-sided subset search. 
}
\details{
The rows and columns of matrices in \code{cor}  are assumed to be in the same order as \code{traits.lab} and columns 
of \code{beta.hat} and \code{sigma.hat}). By definition,
diagonals of the matrices \code{N11}, \code{N00} can not be zero. Diagonals of \code{N10} and \code{N01} are by definition zero.
The matrices \code{N11}, \code{N00}, \code{N10} and \code{N01} are used to calculate a correlation matrix between studies 
that is held
fixed for all SNPs. This is based on the assumption that a fixed proportion of data is missing from all overlap categories,
so that the correlation is unaffected.
Meta-analysis results are based on the inverse variance weighting approach, i.e. Z-scores are weighted according to the sample
sizes in \code{ncase} and \code{ncntl} using the formula:
The output standard errors are approximate (based on inverting pvalues) and are used for constructing confidence
intervals in \code{\link{h.summary}} and \code{\link{traits.forest}}.
}

\references{ Bhattacharjee S, Chatterjee N and others. A new paradigm for meta-analysis of genetic data in the
presence of heterogeneity. Submitted. \cr
}

\seealso{ \code{\link{h.summary}}, \code{\link{traits.forest}} }

\examples{
 # Use the example data
 data(Xdata, package="ASSET")

 # Note the case control distribution
 Xdata$Subtype <- as.character(Xdata$Subtype)
 table(Xdata$Subtype)
 
 # For illustration, let us obtain SNP level summaries for each of the 6 subtypes, treated as 
 # a separate trait (and study) with the controls being split among the 6 studies, to 
 # create independent studies. The first five subtypes are assined 150 controls each, the 
 # remaining 250 controls are assigned to the sixth subtype.
 
 tlab <- paste("T", 1:6, sep="")
 clab <- paste("C", 1:6, sep="")

 cntl.pos <- which(Xdata$Subtype=="Control")
 for(j in 1:5) { Xdata$Subtype[cntl.pos[(j - 1) * 150 + (1:150)]] <- clab[j] }
 Xdata$Subtype[cntl.pos[5 * 150 + (1:250)]] <- clab[6]
 table(Xdata$Subtype)
 
 # Create summary matrices by running separate logistic regressions on each study,
 snps <- c("SNP1", "SNP2", "SNP3")
 adjust <- c("Age", "Sex", "EV1", "EV2")
 beta <- sd <- ncase <- ncntl <- matrix(NA, 3, 6)
 for(j in 1:6)
 {
	res <- h.types(Xdata, response.var="Subtype", snp.vars=snps, adj.vars=adjust
		, types.lab=tlab[j], cntl.lab=clab[j], pool= -1, side=2, logit=TRUE, meth = "Wald")
		
	ncase[, j] <- sapply(snps, function(snp) nrow(na.omit(Xdata[Xdata$Subtype == tlab[j], c("Subtype", snp, adjust)])))
	ncntl[, j] <- sapply(snps, function(snp) nrow(na.omit(Xdata[Xdata$Subtype == clab[j], c("Subtype", snp, adjust)])))
	beta[, j] <- res$Overall.Logistic$beta
	sd[, j] <- res$Overall.Logistic$sd
 }
 exp(beta)
  
 # Now let us call h.traits on these summary data. The correlation matrix does not need to be 
 # specified as studies are independent.

 res <- h.traits(snp.vars=snps, traits.lab=tlab, beta.hat=beta, sigma.hat=sd, ncase=ncase
			, ncntl=ncntl, cor=NULL, cor.numr=FALSE, search=NULL, side=2, meta=TRUE, zmax.args=NULL
			, meth.pval="DLM")
 
 h.summary(res)

}
 

\keyword{ traits } 
