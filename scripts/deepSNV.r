# Rscript --vanilla --slave ../../lauring-variant-pipeline/bin/deepSNV.r ../../lauring-variant-pipeline/lib/ ../04-mark_duplicates/001-0_2a.marked.bam ../04-mark_duplicates/011-PR8control_a.marked.bam ../flu_regions.csv
#

#set seed to make distribiutions determinsitic
set.seed(42)

args <- commandArgs(TRUE)
if (length(args) != 4) {
    stop(paste("Usage:", "deepSNV.r" ,"{library_location} {reference.fasta} {test.bam} {control.bam}"), call.=FALSE)
} 

#print(args)
library.location <- args[1]
reference.fasta <- args[2]
test.bam <- args[3]
control.bam <- args[4]

cat(paste("loading libraries from [", library.location, "]...\n", sep=""))
library(tools)
suppressPackageStartupMessages(library("deepSNV", lib.loc=library.location))
suppressPackageStartupMessages(library("Biostrings", lib.loc=library.location))


test_file_prefix = basename(file_path_sans_ext(test.bam))
control_file_prefix = basename(file_path_sans_ext(control.bam))
output_file_name = paste(test_file_prefix,"--",control_file_prefix, sep="")

cat(paste("loading regions from [", reference.fasta, "]...\n", sep=""))
segments <- fasta.info(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=1, stop=segments, row.names=NULL)
cat(paste("loaded regions: ", paste(regions.bed$chr, collapse=","),"\n"))


cat("calling variants with deepSNV\n")
cat(paste("\ttest [",test.bam,"]\n\tcontrol [",control.bam,"]...\n", sep=""))
deepsnv.result <- deepSNV(test=test.bam, control=control.bam, regions=regions.bed)

cat(paste("saving to [",output_file_name,".vcf].\n", sep=""))
flu_result.vcf <- summary(deepsnv.result, value='VCF')
writeVcf(flu_result.vcf, paste(output_file_name,".vcf", sep=""))

cat(paste("saving to [",output_file_name,".csv].\n", sep=""))
flu_result.df <- deepSNV::summary(deepsnv.result, value='data.frame')
write.csv(flu_result.df, paste(output_file_name,".csv", sep=""))

cat("done.\n")
