
ls.meta <- function(sub, gene.vars, p.vals, side=2)
{	
	k <- length(sub)
	ngene <- length(gene.vars)
	if(ngene > 1) stop("Expected ngene=1")

	z <- pval <- rep(NA, ngene)
	nsub <- sum(sub)
	
	pval <- pchisq(sum(-2 * log(p.vals[sub])), df = 2 * nsub, lower.tail=FALSE)
	z <- qnorm(1 - pval)
	list(z = z, pval = pval)
}
