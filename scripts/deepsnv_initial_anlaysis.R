#################################Coverage Plots####################################################
setwd("/scratch/alauring_fluxm/mccrone/Run_1053/analysis/05_coverage/")

#Get file names
files <- list.files(pattern=".cov", all.files=T, full.names=F)

#Read in coverage files
samples<-c()
c<-1
for(file in files){
  
  #Get name of sample : split file path at slashes and then split file name at "."
  x<-strsplit(file,'.',fixed=T)[[1]][1]
  x<-paste(x,"_cov",sep='')
  samples[c]<-x
  
  
  #Read each file into it's name 
  eval(parse(text=paste(x,'<-read.table(\"',file,'\",header=T,sep=\"\t\")',sep='')))
  c<-c+1
}  
#list of coverage items
cov<-ls(pattern = "_cov")

control_var<-ls(pattern = "Control" )
print(control_var)
control_cov<-eval(parse(text=paste(control_var))) 
 
###########################------Deep-SNV analysis----------##########################

library(deepSNV)
setwd("/scratch/alauring_fluxm/mccrone/Run_1053/analysis/04_mark_duplicates/")

files <- list.files( pattern="bam$", all.files=T, full.names=F)
control_file<- list.files(pattern="Control.*bam$", all.files=T, full.names = F) # sim added for simuated analysis

samples<-c()



##Get ranges for chromosomes and save the vectors for future analysis checking the true variants 

# Segment order comes from the reference sequence oder ( HA M N_A NP NS PA PB1 PB2 )
chromosomes<-levels(control_cov$Segment)
c<-1   #start counter



empieza<-c()
ultima<-c()
max_i<-c() # the max position of segment i 

trim<-10  # the amount to leave off deepSNV analysis from either end of 

for( i in 1:length(chromosomes)){
  seg<-control_cov[control_cov$Segment==chromosomes[i],]
  empieza[c]<-1+trim # There is no 0 position
  max_i[c]<-max(seg$Position)
  ultima[c]<-max_i[c]-trim # no 1 added here because we want analysis to stop here leaving 20 
  c<-c+1
}


regions<-data.frame(chr=chromosomes, start = empieza, stop = ultima)
for(file in files){
  
  #Get name of sample :  split file name at "."
  x<-strsplit(file,'.',fixed=T)[[1]][1]
  x<-paste(x,"_snv",sep='')
  samples[c]<-x
  
  #Read each file DeepSNV into it's name 
  print(paste(x,'<-deepSNV(test=\"',file,'\",control=\"',control_file,'\",regions=regions,q=25)',sep=''))
  eval(parse(text=paste(x,'<-deepSNV(test=\"',file,'\",control=\"',control_file,'\",regions=regions,q=25)',sep='')))
  c<-c+1
}  }  


# Get list of all snv objects

snvs<-ls(pattern="_snv")


## Make summary 


for( mutant in snvs){
  x<-sub(pattern="_snv",replacement="_sum",mutant)
  
  eval(parse(text=paste(x,'<-summary(',mutant,', sig.level = 0.01 , adjust.method=\"BH\")',sep='')))
  
}

#Save data
save(list = ls(all = TRUE), file = "/scratch/alauring_fluxm/mccrone/Run_1053/analysis/Run_1053.RData")


