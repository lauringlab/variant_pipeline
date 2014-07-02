
## ----setup, include=FALSE, cache=FALSE-----------------------------------
require(knitr)
# set global chunk options
opts_chunk$set(fig.path='tmp/deepSNV-', fig.align='center', fig.show='hold', fig.width=4, fig.height=4, out.width='.4\\linewidth', dpi=150)
options(replace.assign=TRUE,width=75)
knit_hooks$set(nice = function(before, options, envir) {
			if (before) par(mar = c(4, 4, .1, .1), mgp=c(2.5,1,0), bty="n")
		})


## ----echo=FALSE, results='asis'------------------------------------------
	print(citation("deepSNV")[1], style="LaTeX")


## ------------------------------------------------------------------------
    library(deepSNV)
	regions <- data.frame(chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 2074, stop=3585)


## ------------------------------------------------------------------------
# HIVmix <- deepSNV(test = "http://www.bsse.ethz.ch/cbg/software/deepSNV/data/test.bam", 
#                   control = "http://www.bsse.ethz.ch/cbg/software/deepSNV/data/control.bam", 
#                   regions=regions, q=10)


## ------------------------------------------------------------------------
	data(HIVmix) # Attach the data instead, as it could fail in routine checks without internet connection.
	show(HIVmix)


## ------------------------------------------------------------------------
	control(HIVmix)[100:110,]
	test(HIVmix)[100:110,]


## ----HIV, nice=TRUE------------------------------------------------------
	plot(HIVmix)


## ------------------------------------------------------------------------
	SNVs <- summary(HIVmix, sig.level=0.05, adjust.method="BH")
	head(SNVs)
	nrow(SNVs)
	min(SNVs$freq.var)


## ------------------------------------------------------------------------
	sum(RF(test(HIVmix), total=T) > 0.01 & RF(test(HIVmix), total=T) < 0.95)


## ------------------------------------------------------------------------
	data(trueSNVs, package="deepSNV")
	table(p.adjust(p.val(HIVmix), method="BH") < 0.05, trueSNVs)


## ----phiX, dev="jpeg", nice=TRUE-----------------------------------------
	## Load data (unnormalized)
	data(phiX, package="deepSNV")
	plot(phiX, cex.min=.5)
	## Normalize data
	phiN <- normalize(phiX, round=TRUE)
	plot(phiN, cex.min=.5)


## ----pval, nice=TRUE-----------------------------------------------------
	p.norm <- p.val(phiN)
    n <- sum(!is.na(p.norm))
    qqplot(p.norm, seq(1/n,1, length.out=n), log="xy", type="S", xlab="P-value", ylab="CDF")
	p.val <- p.val(phiX) 
    points(sort(p.val[!is.na(p.val)]), seq(1/n,1, length.out=n), pch=16, col="grey", type="S", lty=2)
    legend("topleft", c("raw data", "normalized data"), pch=16, col=c("grey", "black"), bty="n", lty=3)
    abline(0,1)


## ----dev="jpeg", nice=TRUE-----------------------------------------------
	data("RCC", package="deepSNV")
	show(RCC)
	plot(RCC, cex.min=.5) 
	RCC.bb = estimateDispersion(RCC, alternative="two.sided")
	plot(RCC.bb, cex.min=.5) 


## ------------------------------------------------------------------------
	RCC.bb@log.lik 
	RCC@log.lik
    RCC.bb@log.lik - RCC@log.lik
	log(4*nrow(test(RCC)))


## ------------------------------------------------------------------------
	summary(RCC, adjust.method="bonferroni")[,1:6]


## ------------------------------------------------------------------------
	tab <- summary(RCC.bb, adjust.method="bonferroni")[,1:6]
	tab


## ----echo=FALSE, results='asis'------------------------------------------
    toLatex(sessionInfo())


