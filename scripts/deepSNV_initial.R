# The goal of this script is to run the initial analyis on bam files 
library(deepSNV)
args<-commandArgs(TRUE)

# args = c(file.bam,control,file.RData)
#################################Functions########################################################
# Make names from the file list for coverage plots
naming<-function(filenames){
  x<-strsplit(filenames,'.',fixed=T)[[1]][1]
}
id<-function(x){
  x<-strsplit(x,".",fixed=T)[[1]][1]
}
#################################Deep-SNV analysis##########################

#-------------SET UP ------------------
chromosomes=c("PB2","PB1","PA","HA","NP","N_A","M","NS")
empieza = c(1,1,1,1,1,1,1,1)
ulitma = c(2341,2341,2233,1775,1565,1413,1027,890)
##Get ranges for chromosomes and save the vectors for future analysis checking the true variants 
trim<- 0 # the amount to leave off deepSNV analysis from either end of 
empieza=empieza+trim
ultima=ultima-trim

regions<-data.frame(chr=chromosomes, start = empieza, stop = ultima)



deepsnv<-deepSNV(test=args[1],control=args[2],regions=regions,q=25,pseudo.count=0.5)



save(list = ls(all = TRUE), file = args[3])



