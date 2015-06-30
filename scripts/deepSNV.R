# Rscript --vanilla --slave ../../lauring-variant-pipeline/bin/deepSNV.r ../../lauring-variant-pipeline/lib/ ../04-mark_duplicates/001-0_2a.marked.bam ../04-mark_duplicates/011-PR8control_a.marked.bam ../flu_regions.csv
#
library(tools)
#set seed to make distribiutions determinsitic
set.seed(42)

args <- commandArgs(TRUE)
if (length(args) != 4) {
    stop(paste("Usage:", "deepSNV.R" ,"{library_location} {reference.fasta} {test.bam} {control.bam}"), call.=FALSE)
} 

#print(args)
library.location <- args[1]
reference.fasta <- args[2]
test.bam <- args[3]
control.bam <- args[4]



test_file_prefix = basename(file_path_sans_ext(test.bam))
control_file_prefix = basename(file_path_sans_ext(control.bam))
output_file_name = paste0("variants/",test_file_prefix)
sample_name=strsplit(test_file_prefix,".",fixed=T)[[1]][1]

if(test_file_prefix==control_file_prefix){
	print("we don't need to run the control!")
	stop()
}

print(test_file_prefix)
print(control_file_prefix)
cat(paste("loading libraries from",library.location,"\n", sep=""))

#library("tools")
#suppressPackageStartupMessages(library("deepSNV", lib.loc=library.location))
#suppressPackageStartupMessages(library("plyr", lib.loc=library.location))
#suppressPackageStartupMessages(library("reshape2", lib.loc=library.location))
library("plyr")
library("deepSNV")
library("reshape2")

cat(paste("loading regions from [", reference.fasta, "]...\n", sep=""))
segments <- fasta.info(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=1, stop=segments, row.names=NULL)
cat(paste("loaded regions: ", paste(regions.bed$chr, collapse=","),"\n"))


cat("calling variants with deepSNV\n")
cat(paste("\ttest [",test.bam,"]\n\tcontrol [",control.bam,"]...\n", sep=""))
deepsnv.result <- deepSNV(test=test.bam, control=control.bam, regions=regions.bed,q=25,pseudo.count=0.5,model='betabin')

consensus_fa<-consensusSequence(test(deepsnv.result,total=T),vector=F,haploid=T)
#cat(paste("saving to [",output_file_name,".vcf].\n", sep=""))
#flu_result.vcf <- summary(deepsnv.result, value='VCF')
#writeVcf(flu_result.vcf, paste(output_file_name,".vcf", sep=""))

deepsnv_sum<-summary(deepsnv.result,sig.level = 0.1, adjust.method="BH")
deepsnv_sum$Id<-sample_name # set the sample name for csv


deepsnv_sum<-subset(deepsnv_sum,var!="-" & ref !="-") # removes the indels
mutate(deepsnv_sum,mutation=paste0(chr,"_",ref,pos,var))->deepsnv_sum




#print(head(deepsnv_sum))
cat(paste("saving to [",output_file_name,".csv].\n", sep=""))
write.csv(deepsnv_sum, paste(output_file_name,".csv", sep=""))

cat(paste("saving to [",output_file_name,".fa].\n",sep=""))
#save(list=consensus_fa,file=paste(output_file_name,".fa",sep=""))
write(as.character(consensus_fa),file=paste(output_file_name,".fa",sep=""))
#cat(paste("saving to [",output_file_name,".RData]\n",sep=""))
#eval(parse(text=paste0("snv_",test_file_prefix,"<-deepsnv.result")))
#save(list=ls(pattern=test_file_prefix),file=paste(output_file_name,".RData",sep=""))

cat("done.\n")
