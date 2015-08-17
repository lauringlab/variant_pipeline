# Rscript --vanilla --slave ../../lauring-variant-pipeline/bin/deepSNV.r ../../lauring-variant-pipeline/lib/ ../04-mark_duplicates/001-0_2a.marked.bam ../04-mark_duplicates/011-PR8control_a.marked.bam ../flu_regions.csv
#
suppressMessages(library(tools))
#set seed to make distribiutions determinsitic
set.seed(42)

args <- commandArgs(TRUE)
if (length(args) != 4) {
    stop(paste("Usage:", "deepSNV.R" ," {reference.fasta} {test.bam} {control.bam} {c(BH,bonferroni)}",sep=""), call.=FALSE)
}

#print(args)
#library.location <- args[1]
reference.fasta <- args[1]
test.bam <- args[2]
control.bam <- args[3]
method<-args[4]


test_file_prefix = basename(file_path_sans_ext(test.bam))
control_file_prefix = basename(file_path_sans_ext(control.bam))
output_file_name = paste0("deepSNV/",test_file_prefix)
sample_name=strsplit(test_file_prefix,".",fixed=T)[[1]][1]
control_name=strsplit(control_file_prefix,".",fixed=T)[[1]][1]
output_file_control=paste0("deepSNV/",control_file_prefix)

print(paste0("test is :",sample_name ))

print(paste0("control is:",control_name))
if(test_file_prefix==control_file_prefix){
	print("we don't need to run the control!")
	stop()
}

#print(test_file_prefix)
#print(control_file_prefix)
#cat(paste("loading libraries from",library.location,"\n", sep=""))

#library("tools")
#suppressPackageStartupMessages(library("deepSNV", lib.loc=library.location))
#suppressPackageStartupMessages(library("plyr", lib.loc=library.location))
#suppressPackageStartupMessages(library("reshape2", lib.loc=library.location))
suppressMessages(library("plyr"))
suppressMessages(library("deepSNV"))
suppressMessages(library("reshape2"))

cat(paste("loading regions from [", reference.fasta, "]...\n", sep=""))
segments <- fasta.info(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=1, stop=segments, row.names=NULL)
cat(paste("loaded regions: ", paste(regions.bed$chr, collapse=","),"\n"))


cat("calling variants with deepSNV\n")
cat(paste("\ttest [",test.bam,"]\n\tcontrol [",control.bam,"]...\n", sep=""))
deepsnv.result <- deepSNV(test=test.bam, control=control.bam, regions=regions.bed,q=25)
deepsnv.result<-estimateDispersion(deepsnv.result,alternative="two.sided")

consensus_fa<-consensusSequence(test(deepsnv.result,total=T),vector=F,haploid=T)
#cat(paste("saving to [",output_file_name,".vcf].\n", sep=""))
#flu_result.vcf <- summary(deepsnv.result, value='VCF')
#writeVcf(flu_result.vcf, paste(output_file_name,".vcf", sep=""))

deepsnv_sum<-summary(deepsnv.result, adjust.method=method)
deepsnv_sum$Id<-sample_name # set the sample name for csv




deepsnv_sum<-subset(deepsnv_sum,var!="-" & ref !="-") # removes the indels
mutate(deepsnv_sum,mutation=paste0(chr,"_",ref,pos,var))->deepsnv_sum


## Coverage ##

cov<-rowSums(test(deepsnv.result,total=T)[,1:4]) # no deletions

# make coverage data.frame

cov.df=data.frame(coverage=cov,concat.pos=1:length(cov))

#setup for segment name and position
prior.seg.length<-c()
for(k in 1:length(regions.bed$chr)){
  prior.seg.length[k]<-sum(regions.bed$stop[1:k])  # the end positions of each segment relative to one sequence not including the trimming step
}
prior.seg.length<-c(0,prior.seg.length)

ddply(cov.df,~concat.pos,summarize, coverage=coverage,concat.pos=concat.pos,chr= as.character(regions.bed$chr[max(which(prior.seg.length<concat.pos))]),
      chr.pos=concat.pos-prior.seg.length[max(which(prior.seg.length<concat.pos))])->cov.df



cov.df$Id<-sample_name # set the sample name for csv

### Control coverage

cov.con<-rowSums(control(deepsnv.result,total=T)[,1:4]) # no deletions

# make coverage data.frame

cov.con.df=data.frame(coverage=cov.con,concat.pos=1:length(cov.con))

ddply(cov.con.df,~concat.pos,summarize, coverage=coverage,concat.pos=concat.pos,chr= as.character(regions.bed$chr[max(which(prior.seg.length<concat.pos))]),
      chr.pos=concat.pos-prior.seg.length[max(which(prior.seg.length<concat.pos))])->cov.con.df



cov.con.df$Id<-control_name # set the sample name for csv


#print(head(deepsnv_sum))
cat(paste("saving summary to [",output_file_name,".csv].\n", sep=""))
write.csv(deepsnv_sum, paste(output_file_name,".csv", sep=""), row.names=FALSE)

#print(head(deepsnv_sum))
cat(paste("saving coverage to [",output_file_name,".cov.csv].\n", sep=""))
write.csv(cov.df, paste(output_file_name,".cov.csv", sep=""), row.names=FALSE)

cat(paste("saving control coverage to [",output_file_control,".cov.csv].\n", sep=""))
write.csv(cov.con.df, paste(output_file_control,".cov.csv", sep=""), row.names=FALSE)

cat(paste("saving to [",output_file_name,".fa].\n",sep=""))
#save(list=consensus_fa,file=paste(output_file_name,".fa",sep=""))
write(as.character(consensus_fa),file=paste(output_file_name,".fa",sep=""))


#cat(paste("saving to [",output_file_name,".RData]\n",sep=""))
#eval(parse(text=paste0("snv_",test_file_prefix,"<-deepsnv.result")))
#save(list=ls(pattern=test_file_prefix),file=paste(output_file_name,".RData",sep=""))

cat("done.\n")
