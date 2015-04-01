# The goal of this script is to run the initial analyis on bam files and hopefully it 
# will be optimized to present the snv's in the samples provided.  The input will be the
# directory output from the pipeline.
library(deepSNV)
args<-commandArgs(TRUE)
#################################Functions########################################################
# Make names from the file list for coverage plots
naming<-function(filenames){
  x<-strsplit(filenames,'.',fixed=T)[[1]][1]
}
id<-function(x){
  x<-strsplit(x,".",fixed=T)[[1]][1]
}
#################################Coverage Plots####################################################
print(paste(args[1],"/05_coverage", sep=""))
setwd(paste(args[1],"/05_coverage", sep=""))

# Get file names
files <- list.files( pattern="*cov$", all.files=T, full.names=F)
#Read in coverage files
coverage.ls<-lapply(files,read.table,header=T,sep='\t')  # read in coverage files and save data.frames in a list
cov.vars<-lapply(files,naming) # make the names for the list of data.frames
names(coverage.ls)<-cov.vars  # name the data frames in the order they were read in 


#################################Deep-SNV analysis##########################


setwd(paste(args[1],"/04_removed_duplicates",sep=""))


##Get ranges for chromosomes and save the vectors for future analysis checking the true variants 
trim<- 20 # the amount to leave off deepSNV analysis from either end of 

# combine coverage data frames into one data frame with idendifiers
coverage.df<-do.call(rbind,coverage.ls)
coverage.df$Id<-unlist(lapply(rownames(coverage.df),id)) # Use the row names to make the Id column

segment_ends<-aggregate( . ~  Segment, coverage.df[,1:2], max)
ultima<-segment_ends$Position-trim
segment_starts<-aggregate( . ~  Segment, coverage.df[,1:2], min)
empieza<-segment_starts$Position+trim

chromosomes<-levels(coverage.df$Segment) # Segment order comes from the reference sequence oder ( HA M N_A NP NS PA PB1 PB2 )
regions<-data.frame(chr=chromosomes, start = empieza, stop = ultima)

# making a list for the deepsnv objects
files <- list.files( pattern="*bam$", all.files=T, full.names=F)
control_file<- list.files(pattern="*con*.*bam$", all.files=T, full.names = T) # sim added for simuated analysis

print(control_file)
deepSnv.ls<-lapply(files,deepSNV,control=control_file,regions=regions,q=25,pseudo.count=0.5) # DeepSNV objects in a list # pseudo.count used in Gerstung examples

deepSnv.vars<-lapply(files,naming) # make the names for the list of data.frames
names(deepSnv.ls)<-deepSnv.vars  # name the data frames in the order they were read in 


#################---COUNTING TRUE and FALSE positives---##################################################

# Get list of all snv objects

summary.ls<-lapply(deepSnv.ls,summary, sig.level = 0.01, adjust.method="bonferroni")
summary.vars<-lapply(files,naming) # make the names for the list of data.frames
names(summary.ls)<-summary.vars  # name the data frames in the order they were read in 

summary.df<-do.call(rbind,summary.ls)
summary.df$Id<-unlist(lapply(rownames(summary.df),id)) # Use the row names to make the Id column


#Save data
tag<-paste(format(Sys.time(), "%Y_%b_%d"),round(runif(1,1,10000),digits = 0),sep="_") # of the form "2015_Jan_19_9781"  where the last 4 digits are a random number to help differentiate between runs done on the same day.

save(list = ls(all = TRUE), file = paste(args[2],"/",tag,".RData",sep=""))

# copy s script to saved directory

script<-"/scratch/alauring_fluxm/mccrone/variant_pipeline/scripts/deepsnv_initial_anlaysis.R"
file.copy(script,paste0(args[2],"/",tag,".R"))
