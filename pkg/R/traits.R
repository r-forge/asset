
h.traits <- function(snp.vars, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor=NULL
	, cor.numr=FALSE, search=NULL, side=2, meta=FALSE, zmax.args=NULL
	, meth.pval=c("DLM", "IS", "B"), pval.args=NULL)
{
	k <- length(traits.lab)
	meth.pval <- meth.pval[1]

	nsnp <- length(snp.vars)
	
	beta.hat <- matrix(beta.hat, ncol=k, byrow=FALSE)
	sigma.hat <- matrix(sigma.hat, ncol=k, byrow=FALSE)
	
	if(nrow(beta.hat) > nsnp) beta.hat <- matrix(beta.hat[snp.vars, ], nrow = nsnp, ncol=k, byrow=FALSE)
	if(nrow(sigma.hat) > nsnp) sigma.hat <- matrix(sigma.hat[snp.vars, ], nrow = nsnp, ncol=k, byrow=FALSE)
	
	rownames(beta.hat) <- snp.vars
	rownames(sigma.hat) <- snp.vars
	colnames(beta.hat) <- traits.lab
	colnames(sigma.hat) <- traits.lab
	
	if(!((length(ncase) == k) && is.null(dim(ncase))) && !(is.matrix(ncase) && (ncol(ncase) == k) && (nrow(ncase) == nsnp)))
		stop("ncase should be a vector of size k or a matrix of dimension nsnp X k")
	if(is.null(dim(ncase))) dim(ncase) <- c(1, k)
	if(!((length(ncntl) == k) && is.null(dim(ncntl))) && !(is.matrix(ncntl) && (ncol(ncntl) == k) && (nrow(ncntl) == nsnp)))
		stop("ncntl should be a vector of size k or a matrix of dimension nsnp X k")
	if(is.null(dim(ncntl))) dim(ncntl) <- c(1, k)

	if(is.null(cor)) cor <- diag(k)
	if(!is.matrix(cor) && !is.list(cor)) stop("cor should be a matrix or a list")
	if(is.list(cor))
	{
		if(!all(c("N11", "N00", "N10", "N01") %in% names(cor)))
			stop("Expected four matrices in list cor with names N11, N00, N10 and N01")
		if(!all(sapply(cor[c("N11", "N00", "N10", "N01")], function(m) { is.matrix(m) && nrow(m) == k  && ncol(m) == k })))
			stop("N11, N00, N10 and N01 should be square matrices of dimension k")
		cor <- corr.mat.logit(cor$N11, cor$N00, cor$N10, cor$N01)
	}

	if(!is.null(colnames(cor)) && any(colnames(cor) != traits.lab)) stop("cor should have column names in same order as traits.lab")
	if(!is.null(rownames(cor)) && any(rownames(cor) != traits.lab)) stop("cor should have row names in same order as traits.lab")
	
	
	pval <- pval1 <- pval2 <- rep(NA, nsnp)
	names(pval) <- names(pval1) <- names(pval2) <- snp.vars

	beta <- beta1 <- beta2 <- rep(NA, nsnp)
	names(beta) <- names(beta1) <- names(beta2) <- snp.vars

	sd <- sd1 <- sd2 <- rep(NA, nsnp)
	names(sd) <- names(sd1) <- names(sd2) <- snp.vars
	
	pheno <- pheno1 <- pheno2 <- matrix(FALSE, nsnp, k)
	colnames(pheno) <- colnames(pheno1) <- colnames(pheno2) <- traits.lab
	rownames(pheno) <- rownames(pheno1) <- rownames(pheno2) <- snp.vars

	OS <- (is.null(search) || search==1)
	TS <- (is.null(search) || search==2)
	meta.res <- sub1.res <- sub2.res <- NULL
	if(meta)
	{ 
		for(j in 1:nsnp)
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			res <- traits.meta(rep(TRUE, k), snp.vars[j], beta.hat[j, ], sigma.hat[j,]
							   , ncase[jj, ], ncntl[jj, ], rmat = cor, side = side, cor.numr = cor.numr, wt.sigma = TRUE)
			pval[j] <- res$pval
			beta[j] <- res$beta
			sd[j] <- res$sd
		}
		meta.res <- list(pval=pval, beta=beta, sd=sd)
	}

	subset0 <- data.frame(pval=rep(NA, nsnp), pheno=rep("NA", nsnp), stringsAsFactors=FALSE)
	subset00 <- data.frame(pval=rep(NA, nsnp), pheno1=rep("NA", nsnp), pheno2=rep("NA", nsnp), stringsAsFactors=FALSE)
	if(OS)
	{ 	
		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 1))
		else pval.args$search <- 1

		for(j in 1:nsnp)
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			res <- h.traits1(k, beta.hat[j, ], sigma.hat[j, ], ncase[jj, ], ncntl[jj, ], rmat=cor
							 , cor.numr=cor.numr, side = side, zmax.args=zmax.args
							 , meth.pval=meth.pval, pval.args=pval.args)
			pval[j] <- res$pval
			pheno[j, ] <- res$pheno
			beta[j] <- res$beta
			sd[j] <- res$sd
		}
		sub1.res <- list(pval=pval, beta=beta, sd=sd, pheno = pheno)
	}

	if(TS)
	{ 
		if(is.null(pval.args) || is.null(pval.args$search)) pval.args <- c(pval.args, list(search = 2))
		else pval.args$search <- 2

		for(j in 1:nsnp)
		{
			jj <- if(nrow(ncase) == 1) 1 else j
			res <- h.traits2(k, beta.hat[j, ], sigma.hat[j, ], ncase[jj, ], ncntl[jj, ], rmat=cor
							 , cor.numr=cor.numr, side = side, zmax.args=zmax.args
							 , meth.pval=meth.pval, pval.args=pval.args)
			pval[j] <- res$pval
			pval1[j] <- res$pval.1
			pval2[j] <- res$pval.2
			
			pheno1[j, ] <- res$pheno.1
			pheno2[j, ] <- res$pheno.2
			
			beta1[j] <- res$beta.1
			beta2[j] <- res$beta.2
			
			sd1[j] <- res$sd.1
			sd2[j] <- res$sd.2
		}
		sub2.res <- list(pval=pval, pval.1 = pval1, beta.1=beta1, sd.1=sd1, pheno.1 = pheno1
						 , pval.2 = pval2, beta.2=beta2, sd.2 = sd2, pheno.2 = pheno2)
	}
	
	ret <- list(Meta=meta.res, Subset.1sided=sub1.res, Subset.2sided=sub2.res)
	ret
}


h.traits1 <- function(k, beta.hat, sigma.hat, ncase, ncntl, rmat, cor.numr, side=2, zmax.args=NULL, meth.pval="DLM", pval.args=NULL)
{
	if(k == 0) return(list(pval = 1, pheno = NULL, beta = NA, sd = NA))
	
	res <- do.call(z.max, c(list(k, 1, side=side, meta.def=traits.meta, meta.args=list(beta.hat=beta.hat, sigma.hat=sigma.hat
				 , ncase=ncase, ncntl=ncntl, rmat=rmat, cor.numr=cor.numr, wt.sigma=FALSE)), zmax.args))
	
	zopt <- as.double(res$opt.z)
	pheno <- as.logical(res$opt.s)
	rm(res)

	if(meth.pval == "DLM")
	{
		if(is.null(pval.args) || !("cor.def" %in% names(pval.args))) pval.args <- c(pval.args, list(cor.def=NULL))
		if(!("cor.args" %in% names(pval.args))) pval.args <- c(pval.args, list(cor.args = list(ncase=ncase, ncntl=ncntl
																			   , rmat=rmat, cor.numr=cor.numr)))
		pval <- do.call(p.dlm, c(list(t.vec=abs(zopt), k=k, side = side), pval.args))
		
	}
	if(meth.pval == "IS")
	{
		pval <- do.call(p.tube, c(list(t.vec=abs(zopt), k=k, side = side, ncase=ncase
								, ncntl=ncntl, rmat=rmat, cor.numr=cor.numr), pval.args))
	}
	if(meth.pval=="B") pval <- p.bon(abs(zopt), k, search = 1, side = side)

	beta <- sd <- NA
		
	res <- traits.meta(pheno, 1, beta.hat, sigma.hat, ncase, ncntl, rmat=rmat, wt.sigma=TRUE)
		
	beta <- res$beta
	if(side == 2) sd <- abs(beta)/qnorm(pval/2, lower.tail=FALSE, log.p = FALSE)
	else sd <- beta/qnorm(pval, lower.tail=FALSE, log.p = FALSE)

	ret <- list(pval = pval, pheno = pheno, beta = beta, sd = sd)
}

h.traits2 <- function(k, beta.hat, sigma.hat, ncase, ncntl, rmat, cor.numr, side=2, zmax.args=NULL, meth.pval="DLM", pval.args=NULL)
{
	sub1 <- which(beta.hat >= 0)
	sub2 <- which(beta.hat < 0)

	k1 <- length(sub1)
	k2 <- length(sub2)

	pval.args1 <- pval.args
	pval.args2 <- pval.args

	if(!is.null(pval.args) && !is.null(pval.args$sizes))
	{
		sizes <- pval.args$sizes
		csizes <- cumsum(sizes)
		csizes0 <- c(0, cumsum(sizes[-length(sizes)]))
		sizes1 <- sapply(1:length(sizes), function(i) { sum(sub1 > csizes0[i] & sub1 <= csizes[i]) })
		sizes2 <- sapply(1:length(sizes), function(i) { sum(sub2 > csizes0[i] & sub2 <= csizes[i]) })
		sizes1 <- sizes1[sizes1 > 0]
		sizes2 <- sizes2[sizes2 > 0]
		pval.args1$sizes <- sizes1
		pval.args2$sizes <- sizes2
	}
	res1 <- h.traits1(k1, beta.hat=beta.hat[sub1], sigma.hat=sigma.hat[sub1], ncase=ncase[sub1]
					  , ncntl=ncntl[sub1], rmat=rmat[sub1, sub1]
					  , cor.numr=cor.numr, side=side, zmax.args=zmax.args
					  , meth.pval=meth.pval, pval.args=pval.args1)
	
	res2 <- h.traits1(k2, beta.hat=beta.hat[sub2], sigma.hat=sigma.hat[sub2], ncase=ncase[sub2]
					  , ncntl=ncntl[sub2], rmat[sub2, sub2]
					  , cor.numr=cor.numr, side=side, zmax.args=zmax.args
					  , meth.pval=meth.pval, pval.args=pval.args2)
	
	
	zopt <- if(res1$pval <= 0 || res2$pval <= 0) NA else -2 * (log(res1$pval) + log(res2$pval))		
	
	if(k1 > 0 && k2 > 0) 
	{ 
		if(meth.pval != "B") pval <- pchisq(zopt, df = 4, lower.tail = FALSE)
		if(meth.pval == "B") pval <- pchisq(zopt, df = 4, lower.tail = FALSE) * ((2^k1 - 1) * (2^k2 - 1))
	}
	if(k1 == 0 || k2 == 0 )
	{
		if(meth.pval != "B") pval <- pchisq(zopt, df = 2, lower.tail = FALSE)
		if(meth.pval == "B") pval <- pchisq(zopt, df = 2, lower.tail = FALSE) * ((2^k - 1))
	}

	pheno1 <- rep(FALSE, k) ; pheno1[sub1] <- res1$pheno
	pheno2 <- rep(FALSE, k) ; pheno2[sub2] <- res2$pheno
	
	ret <- list(pval = pval, pval.1=res1$pval, pval.2=res2$pval, pheno.1 = pheno1, pheno.2 = pheno2
				, beta.1 = res1$beta, sd.1 = res1$sd, beta.2 = res2$beta, sd.2 = res2$sd)
	ret
}	

traits.meta <- function(sub, snp.vars, beta.hat, sigma.hat, ncase, ncntl, rmat, side=2, cor.numr=FALSE, wt.sigma=FALSE)
{	
	k <- length(sub)
	if(is.null(dim(rmat))) dim(rmat) <- c(k, k)
	nsnp <- length(snp.vars)
	if(nsnp > 1) stop("Expected nsnp=1")

	z <- beta <- sd <- pval <- rep(NA, nsnp)
	nsub <- sum(sub)
	
	if(!wt.sigma)
	{
		zvec <- beta.hat[sub]/sigma.hat[sub]
		rneff <- sqrt((ncase[sub] * ncntl[sub])/(ncase[sub] + ncntl[sub]))
		if(cor.numr) rneff <- mySolve(rmat[sub, sub], rneff)
		numr <- sum(zvec * rneff)
		denr <- sum(outer(rneff, rneff) * matrix(rmat[sub, sub], nsub, nsub))
		z <- ifelse(is.na(denr) | is.nan(denr) | denr == 0, NA, numr/sqrt(denr))
	} else
	{
		Sigma <- diag(sigma.hat[sub], nsub) %*% matrix(rmat[sub, sub], nsub, nsub) %*% diag(sigma.hat[sub], nsub)
		if(cor.numr) wt <- mySolve(Sigma, rep(1, nsub))
		else wt <- (1/sigma.hat[sub]^2)
		
		numr <- sum(beta.hat[sub] * wt)
		denr <- sum(outer(wt, wt) * matrix(Sigma, nsub, nsub))
		beta <- numr/sum(wt)
		sd <- 1/sqrt(denr)
		z <- ifelse(is.na(denr) | is.nan(denr) | denr == 0, NA, numr/sqrt(denr))
	}
	if(side == 2) pval <- 2 * pnorm(abs(numr/sqrt(denr)), lower.tail=FALSE)
	else pval <- pnorm(numr/sqrt(denr), lower.tail=FALSE)
	
	list(z = z, beta = beta, sd = sd, pval = pval)
}

traits.forest <- function(snp.var, traits.lab, beta.hat, sigma.hat, ncase, ncntl, cor = NULL, rlist = NULL
, level = 0.05, p.adj = TRUE, digits = 2)
{
	if(length(snp.var) > 1) stop("Length of snp.var should be 1")	
	if(!(snp.var %in% rownames(beta.hat))) stop("rownames of beta.hat does not have snp.var")
	k <- length(traits.lab)
	
	if(!is.null(rlist))
	{
		if(is.null(rlist$Meta)) stop("Missing Component in rlist: Meta")
		if(is.null(rlist$Subset.1sided)) stop("Missing Component in rlist: Subset.1sided")
		if(is.null(rlist$Subset.2sided)) stop("Missing Component in rlist: Subset.2sided")
	} else
	{
		rlist <- h.traits(snp.vars = snp.var, traits.lab = traits.lab, beta.hat = beta.hat, sigma.hat = sigma.hat
						  , ncase = ncase, ncntl = ncntl, cor = cor, search = NULL, meta = TRUE)
	}
	
	beta <- beta.hat[snp.var, ]
	sd <- sigma.hat[snp.var, ]
	res <- list(beta = beta, sd = sd, pval = 2 * pnorm(abs(beta/sd), lower.tail=FALSE))
	
	h.forest(k, snp.var, traits.lab, rlist, res, side = 2, level=level, p.adj = p.adj, digits=digits)
}

