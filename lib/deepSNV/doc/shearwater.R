
## ----setup, include=FALSE, cache=FALSE-----------------------------------
require(knitr)
# set global chunk options
opts_chunk$set(fig.path='tmp/deepSNV-', fig.align='center', fig.show='hold', fig.width=4, fig.height=4, out.width='.4\\linewidth', dpi=150)
options(replace.assign=TRUE,width=75)
knit_hooks$set(nice = function(before, options, envir) {
			if (before) par(mar = c(4, 4, .1, .1), mgp=c(2.5,1,0), bty="n")
		})


## ----echo=FALSE, results='asis'------------------------------------------
	print(citation("deepSNV")[2], style="LaTeX")


## ----fig.width=5, fig.height=5, out.width='.6\\linewidth'----------------
library(deepSNV)
library(RColorBrewer)
n <-  100 ## Coverage
n_samples <- 1000 ## Assume 1000 samples
x <-  0:20 ## Nucleotide counts
X <-  cbind(rep(x, each = length(x)), rep(x, length(x))) ## All combinations forward and reverse
par(bty="n", mgp = c(2,.5,0), mar=c(3,3,2,2)+.1, las=1, tcl=-.33, mfrow=c(2,2))
for(nu in 10^c(-4,-2)){ ## Loop over error rates
	## Create counts array with errors
	counts = aperm(array(c(rep(round(n_samples*n* c(nu,1-nu,nu,1-nu)), each=nrow(X)), cbind(n - X, X)[,c(3,1,4,2)]), 
					dim=c(nrow(X) ,4,2)), c(3,1,2))
	for(rho in c(1e-4, 1e-2)){ ## Loop over dispersion factors
		## Compute Bayes factors
		BF = bbb(counts, rho=rho, model="OR", return="BF")
		## Plot
		image(z=log10(matrix(BF[2,,1], nrow=length(x))), 
				x=x, 
				y=x, 
				breaks=c(-100,-8:0), 
				col=rev(brewer.pal(9,"Reds")), 
				xlab = "Forward allele count",
				ylab="Backward allele count", 
				main = paste("rho =", format(rho, digits=2), "nu = ", format(nu, digits=2)), 
				font.main=1)
		text(X[,1],X[,2],ceiling(log10(matrix(BF[2,,1], nrow=length(x)))), cex=0.5)
	}
}


## ----fig.width=5,fig.height=5, out.width='.6\\linewidth'-----------------
par(bty="n", mgp = c(2,.5,0), mar=c(3,3,2,2)+.1, las=1, tcl=-.33, mfrow=c(2,2))
for(nu in 10^c(-4,-2)){ ## Loop over error rates
	## Create counts array with errors
	counts = aperm(array(c(rep(round(n_samples*n* c(nu,1-nu,nu,1-nu)), each=nrow(X)), cbind(n - X, X)[,c(3,1,4,2)]), 
					dim=c(nrow(X) ,4,2)), c(3,1,2))
	for(rho in c(1e-4, 1e-2)){ ## Loop over dispersion factors
		## Compute Bayes factors, mode = "AND"
		BF = bbb(counts, rho=rho, model="AND", return="BF")
		## Plot
		image(z=log10(matrix(BF[2,,1], nrow=length(x))), 
				x=x, 
				y=x, 
				breaks=c(-100,-8:0), 
				col=rev(brewer.pal(9,"Reds")), 
				xlab = "Forward allele count",
				ylab="Backward allele count", 
				main = paste("rho =", format(rho, digits=2), "nu = ", format(nu, digits=2)), 
				font.main=1)
		text(X[,1],X[,2],ceiling(log10(matrix(BF[2,,1], nrow=length(x)))), cex=0.5)
	}
}


## ----rho, fig.width=4,fig.height=4, out.width="6cm", out.height="6cm"----
rho = 10^seq(-6,-1)
rhoHat <- sapply(rho, function(r){
			sapply(1:100, function(i){
						n = 100
						X = rbetabinom(1000, n, 0.01, rho=r)
						X = cbind(X, n-X)
						Y = array(X, dim=c(1000,1,2))
						deepSNV:::estimateRho(Y, Y/n, Y < 1000)[1,1]})
		})
par(bty="n", mgp = c(2,.5,0), mar=c(3,4,1,1)+.1,  tcl=-.33)
plot(rho, type="l", log="y", xaxt="n", xlab="rho", ylab="rhoHat", xlim=c(0.5,6.5), lty=3)
boxplot(t(rhoHat+ 1e-7) ~ rho, add=TRUE, col="#FFFFFFAA", pch=16, cex=.5, lty=1, staplewex=0)
points(colMeans(rhoHat), pch="*", col="red", cex=2)


## ----eval=FALSE----------------------------------------------------------
## ## Not run..
## ## Load TxDb
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
## seqlevels(txdb) <- sub("chr","",seqlevels(txdb))
## 
## ## Make prior
## regions <- reduce(exons(txdb, vals=list(gene_id='7157'))) ## TP53 exons
## cosmic <- readVcf("CosmicCodingMuts_v63_300113.vcf.gz", "hg19", param=ScanVcfParam(which=regions))
## pi <- makePrior(cosmic, regions, pi.gene = 1)


## ----prior, fig.width=8,fig.height=4, out.width="12cm", out.height="6cm"----
## Load pi
data(pi, package="deepSNV")

## Plot
par(bty="n", mgp = c(2,.5,0), mar=c(3,3,2,2)+.1, tcl=-.33)
plot(pi[,1], type="h", xlab="Position", ylab="Prior", col=brewer.pal(5,"Set1")[1], ylim=c(0,0.075))
for(j in 2:5)
	lines(pi[,j], type="h", col=brewer.pal(5,"Set1")[j])
legend("topleft", col=brewer.pal(5,"Set1"), lty=1, bty="n", c("A","T","C","G","del"))


## ------------------------------------------------------------------------
## Load data from deepSNV example
regions <- GRanges("B.FR.83.HXB2_LAI_IIIB_BRU_K034", IRanges(start = 3120, end=3140))
files <- c(system.file("extdata", "test.bam", package="deepSNV"), system.file("extdata", "control.bam", package="deepSNV"))
counts <- loadAllData(files, regions, q=10)
dim(counts)


## ------------------------------------------------------------------------
## Run (bbb) computes the Bayes factor
bf <- bbb(counts, model = "OR", rho=1e-4)
dim(bf)
vcf <- bf2Vcf(bf, counts, regions, cutoff = 0.5, samples = files, prior = 0.5, mvcf = TRUE) 
show(vcf)


## ----fig.width=4,fig.height=4--------------------------------------------
## Shearwater Bayes factor under AND model
bf <- bbb(counts, model = "AND", rho=1e-4)
## deepSNV P-value with combine.method="fisher" (product)
dpSNV <- deepSNV(test = files[1], control = files[2], regions=regions, q=10, combine.method="fisher")
## Plot
par(bty="n", mgp = c(2,.5,0), mar=c(3,3,2,2)+.1, tcl=-.33)
plot(p.val(dpSNV), bf[1,,]/(1+bf[1,,]), log="xy",
		xlab = "P-value deepSNV",
		ylab = "Posterior odds shearwater"
		)


## ----eval=FALSE----------------------------------------------------------
## ## Not run
## files <- dir("bam", pattern="*.bam$", full.names=TRUE)
## MC_CORES <- getOption("mc.cores", 2L)
## vcfList <- list()
## for(gene in levels(mcols(regions)$Gene)){
## 	rgn <-  regions[mcols(regions)$Gene==gene]
## 	counts <-  loadAllData(files, rgn, mc.cores=MC_CORES)
## 	## Split into
## 	BF <-  mcChunk("bbb", split = 200, counts, mc.cores=MC_CORES)
## 	COSMIC <-  readVcf("CosmicCodingMuts_v63_300113.vcf.gz", "GRCh37", param=ScanVcfParam(which=rgn) )
## 	prior <- makePrior(COSMIC, rgn, pi.mut = 0.5)
## 	vcfList[[gene]] <- bf2Vcf(BF = BF, counts=counts, regions=rgn, samples = files, cutoff = 0.5, prior = prior)
## }
## ## Collapse vcfList
## vcf <- do.call(rbind, vcfList)


## ----echo=FALSE, results='asis'------------------------------------------
    toLatex(sessionInfo())


