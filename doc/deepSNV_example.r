library(deepSNV)
flu.regions <- read.table("/Users/pulintz/Box Sync/Lauring/titration_experiment/flu_regions.csv", sep=',', header=T, na.strings="^^")
flu001_011 = deepSNV(test="output-001-CACTGT.bam", control="output-011-ACGCTT.bam", regions=flu.regions)
flu003_011 = deepSNV(test="output-003-ACAGCT.bam", control="output-011-ACGCTT.bam", regions=flu.regions)
flu005_011 = deepSNV(test="output-005-ACCAGT.bam", control="output-011-ACGCTT.bam", regions=flu.regions)
flu001_011.vcf <- summary(flu001_011, value='VCF')
flu003_011.vcf <- summary(flu003_011, value='VCF')
flu005_011.vcf <- summary(flu005_011, value='VCF')
flu001_011_result_df <- summary(flu001_011, value='data.frame')
flu003_011_result_df <- summary(flu003_011, value='data.frame')
flu005_011_result_df <- summary(flu005_011, value='data.frame')


bed <- read.table("/Users/pulintz/Box Sync/Lauring/Miseq March/Variant Positions/flu.bed", 
                sep="\t", stringsAsFactors=FALSE, header=FALSE, na.strings="^^")
key.bed <- paste(as.character(bed[,1]), 
                 as.character(bed[,2]), 
                 sapply(bed[,4], function(x) { unlist(strsplit(x, "->"))[2] }), sep="|")
flu001.key <- paste(as.character(flu001_011_result_df[,"chr"]),
            as.character(flu001_011_result_df[,"pos"]),
            as.character(flu001_011_result_df[,"var"]), sep="|")
flu003.key <- paste(as.character(flu003_011_result_df[,"chr"]),
            as.character(flu003_011_result_df[,"pos"]),
            as.character(flu003_011_result_df[,"var"]), sep="|")
flu005.key <- paste(as.character(flu005_011_result_df[,"chr"]),
            as.character(flu005_011_result_df[,"pos"]),
            as.character(flu005_011_result_df[,"var"]), sep="|")

flu001_011_result_df.bed <- flu001_011_result_df[flu001.key %in% key.bed,]
flu003_011_result_df.bed <- flu003_011_result_df[flu003.key %in% key.bed,]
flu005_011_result_df.bed <- flu005_011_result_df[flu005.key %in% key.bed,]


for.box <- matrix(NA, ncol=3, nrow=length(flu001_011_result_df.bed[,"freq.var"]))
for.box[,1] <- flu001_011_result_df.bed[,"freq.var"]
for.box[1:length(flu003_011_result_df.bed[,"freq.var"]),2] <- flu003_011_result_df.bed[,"freq.var"]
for.box[1:length(flu005_011_result_df.bed[,"freq.var"]),3] <- flu005_011_result_df.bed[,"freq.var"]
colnames(for.box) <- c("0.20", "0.05", "0.0125")

boxplot(for.box, notch=TRUE)


tiff("boxplot_flu.tiff")
boxplot(for.box, notch=TRUE, axes=FALSE, xlab="True SNV Frequency", ylab="Observed SNV Frequency")
axis(2, at=seq(0,1,by=0.1))
axis(1, at=c(1,2,3), colnames(for.box))
abline(h=0.0125, col="gray", lty=3)
abline(h=0.2, col="gray", lty=3)
abline(h=0.05, col="gray", lty=3)
box()
dev.off()

library(vioplot)

tiff("vioplot_flu.tiff")
vioplot(flu001_011_result_df.bed[,"freq.var"], flu003_011_result_df.bed[,"freq.var"], flu005_011_result_df.bed[,"freq.var"],
        names=c("0.20", "0.05", "0.0125"), col="skyblue")
title(xlab="True SNV Frequency", ylab="Observed SNV Frequency")
abline(h=0.0125, col="gray", lty=3)
abline(h=0.2, col="gray", lty=3)
abline(h=0.05, col="gray", lty=3)
dev.off()


abline(h=0, col="gray", lty=2)
abline(h=0, col="gray", lty=2)





 axes=FALSE)
axis(2, at=seq(0,1,by=0.1))
axis(1, at=c(1,2,3), c("0.20", "0.05", "0.0125"))
box()


grid(10,10)

axis(2, at=seq(0,1,by=0.1))

summary(flu003_011_result_df.bed[,"freq.var"])


flu003_011_result_df.bed[,"freq.var"]
flu005_011_result_df.bed[,"freq.var"]
