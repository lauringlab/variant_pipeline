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
# Fill in the data frame so they are all the same size and can be put in a matrix
fill.in<-function(x,endpoint=size,seg_name=seg){
  if(nrow(x) < endpoint){
    diff <-endpoint-nrow(x)
    m<-data.frame(Segment=rep(seg_name,times=diff),Position=rep(0,times=diff),Coverage=rep(0,times=diff))
  x<-rbind(x,m)
  }
  return(x)
}
coverage_plots<-function(seg,data.ls=coverage.ls){
  data.ls<-lapply(data.ls,function(x,segment=seg){ y<- x[x$Segment==segment,]; return(y)}) #Just get the segments of interest
  
  size<-max(vapply(data.ls,function(x){nrow(x)},FUN.VALUE = 1))
  
  data.ls<-lapply(data.ls,fill.in,endpoint=size,seg_name=seg)
  
  position<-vapply(X = data.ls,function(x){y<-x$Position;return(y)},numeric(size)) # one long vector
  coverage<-vapply(X = data.ls,function(x){y<-x$Coverage;return(y)},numeric(size))
  
  plot(position,coverage,main=paste(seg), type='l')
}

nplot<-function(x,t, line){
  plot(x)
  title(t)
  abline(h=line)
}

id<-function(x){
  x<-strsplit(x,".",fixed=T)[[1]][1]
}

freq_pval<-function(sumfile){
  
  plot(sumfile$freq.var,sumfile$p.val, log="xy", main=sumfile$Id[1])
  TP<-sumfile[sumfile$Mutation==TRUE,]   # only look at the True mutants
  PB1<- data.frame(TP$freq.var[TP$chr=="PB1"],TP$p.val[TP$chr=="PB1"])
  PB2<-data.frame(TP$freq.var[TP$chr=="PB2"],TP$p.val[TP$chr=="PB2"])
  PA <-data.frame(TP$freq.var[TP$chr=="PA"],TP$p.val[TP$chr=="PA"])
  HA<-data.frame(TP$freq.var[TP$chr=="HA"],TP$p.val[TP$chr=="HA"])
  NP<-data.frame(TP$freq.var[TP$chr=="NP"],TP$p.val[TP$chr=="NP"])
  N_A<-data.frame(TP$freq.var[TP$chr=="N_A"],TP$p.val[TP$chr=="N_A"])
  NS<-data.frame(TP$freq.var[TP$chr=="NS"],TP$p.val[TP$chr=="NS"])
  M<-data.frame(TP$freq.var[TP$chr=="M"],TP$p.val[TP$chr=="M"])
  deletions<-data.frame(TP$freq.var[TP$var=="-"],TP$p.val[TP$var=="-"])
  coverage_filtered<-data.frame(sumfile$freq.var[sumfile$Cov_filter==FALSE],sumfile$p.val[sumfile$Cov_filter==FALSE])
    
    points(PB1, pch = 16, col = "chartreuse4")
    points(PB2, pch = 16, col = "red")
    points(PA, pch = 16, col = "blue")
    points(HA, pch = 16, col = "gray")
    points(NP, pch = 16, col = "yellow")
    points(N_A, pch = 16, col = "brown")
    points(NS, pch = 16, col = "orange")
    points(M, pch = 16, col = "green")
    points(deletions, pch = 16, col = "purple")
    points(coverage_filtered,pch = 16, col = "pink")
    #legend("bottomleft",c("PB1",'PB2','PA','HA','NP','NA','NS','M','deletions'),pch=16,col= c("chartreuse4","red","blue","gray","yellow","brown","orange","green",'purple'))

}

false.neg<-function(sumfile,true=True_snv){
  x<-merge(x = True_snv,y=sumfile,by.y=c("chr","pos","ref","var"),by.x=c("Name","Seq.Pos","Ref","Allele.1"),all.x = TRUE)[,c("Name","Seq.Pos","Ref","Allele.1","p.val","freq.var","Coverage","Cov_filter","Detected")]
 return(x)
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

#-------------------------Plotting coverage-----------------------------------
#MAKE Plots
#par(mfrow=c(4,2))
#coverage_plots('PB2')
#coverage_plots('PB1')
#coverage_plots('PA')
#coverage_plots('HA')
#coverage_plots('NP')
#coverage_plots('N_A')
#coverage_plots('M')
#coverage_plots('NS')
#par(mfrow = c(1,1)) #reset Par

#################################Deep-SNV analysis##########################


setwd(paste(args[1],"/04_mark_duplicates",sep=""))


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
control_file<- list.files(pattern="wt.*bam$", all.files=T, full.names = T) # sim added for simuated analysis

print(control_file)
deepSnv.ls<-lapply(files,deepSNV,control=control_file,regions=regions,q=25) # DeepSNV objects in a list

deepSnv.vars<-lapply(files,naming) # make the names for the list of data.frames
names(deepSnv.ls)<-deepSnv.vars  # name the data frames in the order they were read in 


#################---COUNTING TRUE and FALSE positives---##################################################

# Get list of all snv objects

summary.ls<-lapply(deepSnv.ls,summary, sig.level = 0.01, adjust.method="BH")
summary.vars<-lapply(files,naming) # make the names for the list of data.frames
names(summary.ls)<-summary.vars  # name the data frames in the order they were read in 

summary.df<-do.call(rbind,summary.ls)
summary.df$Id<-unlist(lapply(rownames(summary.df),id)) # Use the row names to make the Id column


#Save data
save(list = ls(all = TRUE), file = paste(args[1],"/deepsnv_virus_control.RData",sep=""))


