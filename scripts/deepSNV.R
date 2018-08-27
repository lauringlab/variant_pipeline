# Rscript --vanilla --slave ../../lauring-variant-pipeline/bin/deepSNV.r ../../lauring-variant-pipeline/lib/ ../04-mark_duplicates/001-0_2a.marked.bam ../04-mark_duplicates/011-PR8control_a.marked.bam ../flu_regions.csv
#
###############################################################################
################## Load packages and set up arguments #########################
###############################################################################


#set seed to make distribiutions determinsitic
set.seed(42)

args <- commandArgs(TRUE)
if (length(args) != 11) {
    stop(paste("Usage:", "deepSNV.R" ," {reference.fasta} {test.bam} {control.bam} {c(BH,bonferroni)} {p.val.cut} {c(fisher,average,max)}  {c(two.sided,one.sided,bin)} {stringent_freq} output.csv output.fa package_library" ,sep=""), call.=FALSE)
 }
#print(args)
reference.fasta <- args[1]
test.bam <- args[2]
control.bam <- args[3]
method<-args[4]
p.cut<-as.numeric(args[5]) # the p.value cut off for samples to be included in downstream analysis.
p.com.meth<-args[6] # combine method for strands
disp<-args[7]
stringent_freq<-as.numeric(args[8])
csv<-args[9]
fa<-args[10]
library.location <- args[11]

.libPaths(library.location)
print(.libPaths())
suppressMessages(library("tools"))
suppressMessages(library("plyr"))
suppressMessages(library("deepSNV"))
suppressMessages(library("reshape2"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyr"))


# handle the file names to be used now and later.
test_file_prefix = basename(file_path_sans_ext(test.bam))
control_file_prefix = basename(file_path_sans_ext(control.bam))
output_file_name = paste0("deepSNV/",test_file_prefix)
sample_name=strsplit(test_file_prefix,".",fixed=T)[[1]][1]
control_name=strsplit(control_file_prefix,".",fixed=T)[[1]][1]
output_file_control=paste0("deepSNV/",control_file_prefix)

print(paste0("test is :",sample_name ))

print(paste0("control is:",control_name))

###############################################################################
################################ Functions ####################################
###############################################################################

get_counts <- function(x,counts,counts.con){ 
  # This function takes in a data frame of frequencies, and test and control count matrixes. 
  # It uses the row column in the data frame to add base counts for each entry 
  pos<-x$pos
  mat_pos<-x$row
  base = x$var
  n.tst.fw = counts[mat_pos,base]
  n.tst.bw = counts[mat_pos,tolower(base)]
  cov.tst.fw = sum(counts[mat_pos,c(1:5)])
  cov.tst.bw = sum(counts[mat_pos,c(6:10)])
  
  # now for the control
  n.ctrl.fw = counts.con[mat_pos,base]
  n.ctrl.bw = counts.con[mat_pos,tolower(base)]
  cov.ctrl.fw = sum(counts.con[mat_pos,c(1:5)])
  cov.ctrl.bw = sum(counts.con[mat_pos,c(6:10)])
  
  out = data.frame(chr = x$chr,
                   pos = x$pos,
                   var = base,
                   n.tst.fw = n.tst.fw,
                   n.tst.bw = n.tst.bw,
                   cov.tst.fw = cov.tst.fw,
                   cov.tst.bw = cov.tst.bw,
                   n.ctrl.fw = n.ctrl.fw,
                   n.ctrl.bw = n.ctrl.bw,
                   cov.ctrl.fw = cov.ctrl.fw,
                   cov.ctrl.bw = cov.ctrl.bw)
  
}
## Set up coordinates 
set_cord<-function(deepsnv.result,freqs){ # This funciton takes in a deepsnv result and of frequencies (from RF command). It appends the chr and pos (on chr) to the frequency data
  positions<-data.frame(chr=character(),
                      pos = integer())
  cors<-coordinates(deepsnv.result)
  #print(head(cors))
  #print(tail(cors))
  #for(i in 1:(length(names(cors))/2)){ # the names of the cor are "chr.i" and "pos.i" for each segment tested. This takes them in order and puts them in columns chr and pos
  #  seg = paste('chr',i,sep=".")
  #  pos_col = paste('pos',i,sep=".")
  #  seg.df<-data.frame(chr=cors[,seg],pos = cors[,pos_col])
  #  positions<-rbind(positions,seg.df)
  #}
  #cbind(positions,freqs)->frequencies.df
   cbind(cors,freqs)->frequencies.df
}



###############################################################################
########################### Calling Variants ##################################
###############################################################################

# Get regions from alignment reference
cat(paste("loading regions from [", reference.fasta, "]...\n", sep=""))
fastaData<-readDNAStringSet(reference.fasta)
regions.bed <- data.frame(chr = names(fastaData), start=1, stop=width(fastaData), row.names=NULL)
cat(paste("loaded regions: ", paste(regions.bed$chr, collapse=","),"\n"))


cat("calling variants with deepSNV\n")
cat(paste("\ttest [",test.bam,"]\n\tcontrol [",control.bam,"]...\n", sep=""))
deepsnv.result <- deepSNV(test=test.bam, control=control.bam, regions=regions.bed,combine.method=p.com.meth,q=30) # This calls the variants using the binomial

# adjust model if needed.

if(disp=="two.sided"){
deepsnv.result<-estimateDispersion(deepsnv.result,alternative="two.sided")
}else if (disp=="one.sided"){
deepsnv.result<-estimateDispersion(deepsnv.result) # one sided is the default
}else{
deepsnv.result<-deepsnv.result
}

# Get summary data  - this sig.level is just for this stage to limit the size of the output - we will filter more stringently later in sift stage.
deepsnv_sum<-summary(deepsnv.result,sig.level=p.cut, adjust.method=method)
# filter to stringent freq
deepsnv_sum<-subset(deepsnv_sum,freq.var<stringent_freq)
print("132")
###############################################################################
################### Calling non stringent variants ############################
###############################################################################

  ## this analysis replicates the deepsnv summary method for variants above the stringent cut off
RF(test(deepsnv.result),total=T)->frequencies # Get the frequencies of all the bases at each position

# set up coordinates
frequencies.df<-set_cord(deepsnv.result,frequencies)

names(frequencies.df)[names(frequencies.df)=="-"]<-"indel" # for ease of handling below
frequencies.df$row<-1:nrow(frequencies.df) # this row will be used to find the proper position in the count matrices ect. We haven't subsetted the frequencies.df so there is a row for each position querried.

frequencies.df %>% gather(var,freq.var,A:indel)->frequencies.df # long form
subset(frequencies.df,freq.var>=stringent_freq & var!="indel")->less_stringent # subset to only those that qualify for nonstringency
print("148")
if(nrow(less_stringent)>0){ # if the stringent freq is set to 1 or higher than we don't do this and require deepSNV to call variants
  less_stringent$var[less_stringent$var=="indel"]<- "-" # correct name
  # Get the counts for the test and control samples - for all positions
  test_counts<-test(deepsnv.result)
  ctrl_counts<-control(deepsnv.result)

  less_stringent %>% adply(1,get_counts,test_counts,ctrl_counts) -> less_stringent.format # formatted like deepsnv df -  uses row column to find appropriate counts in test control matrices

  ## Now we get the reference base from the control matrix and process it to add chr and pos columns - this could also be do using consensusSequance - I have checked - they give the same results.
  RF(control(deepsnv.result),total=T)->con.freq
  df.con<-set_cord(deepsnv.result,con.freq)

  names(df.con)[ncol(df.con)]<-"indel"
  df.con %>% gather(ref,freq.var,A:indel)->con.freqs
  subset(con.freqs,freq.var>0.5,select=c(chr,pos,ref))->con.major # only want the major bases - consensus.

  less_stringent.final<-join(less_stringent.format,con.major,by = c('chr','pos'),type = 'left') # join by chr and pos
  # Add NA to columns that deepSNV uses but these don't  - p.val sigma.freq.var ect.
  for(c in names(deepsnv_sum)){
    if(c %in% names(less_stringent.final)==F){
      less_stringent.final[,c]=NA
    }
  }

  less_stringent.final<-subset(less_stringent.final,select=-c(row)) # remove the row column
  deepsnv_sum<-rbind(deepsnv_sum,less_stringent.final) # combine stringent and nonstringent column
}
print("176")
deepsnv_sum<-deepsnv_sum[order(deepsnv_sum$chr,deepsnv_sum$pos),] # order the result
print("178")
if(nrow(deepsnv_sum)>0){
  print("180")
  deepsnv_sum$Id<-sample_name # set the sample name for csv
  deepsnv_sum<-subset(deepsnv_sum,!(var=="-" | ref =="-")) # removes the indels that are below the consensus level
  mutate(deepsnv_sum,mutation=paste0(chr,"_",ref,pos,var))->deepsnv_sum
}else{
  deepsnv_sum$Id<-NULL
  deepsnv_mutation<-NULL
}


###############################################################################
########################## Consensus Sequences ################################
###############################################################################
print("193")
consensus_fa<-consensusSequence(test(deepsnv.result,total=T),vector=F,haploid=T)
control_fa<-consensusSequence(control(deepsnv.result,total=T),vector=F,haploid=T)

print('196')
###############################################################################
############################### Coverage ######################################
###############################################################################

cov<-rowSums(test(deepsnv.result,total=T)[,1:4]) # set to 4 for no deletions sum of all bases calls at each position (row)

# make coverage data.frame
cov.df=data.frame(coverage=cov,concat.pos=1:length(cov))

#setup for segment name and position - 

cov.df<-set_cord(deepsnv.result,cov.df)
names(cov.df)[names(cov.df)=="pos"]<-'chr.pos'

cov.df$Id<-sample_name # set the sample name for csv

### Control coverage

cov.con<-rowSums(control(deepsnv.result,total=T)[,1:4]) # set to 4 for no deletions

# make coverage data.frame

cov.con.df=data.frame(coverage=cov.con,concat.pos=1:length(cov.con))

cov.con.df<-set_cord(deepsnv.result,cov.con.df)
names(cov.con.df)[names(cov.con.df)=="pos"]<-'chr.pos'
cov.con.df$Id<-control_name # set the sample name for csv
print("224")
###############################################################################
########################## Saving output ######################################
###############################################################################

cat(paste("saving summary to [",csv,"].\n", sep=""))
write.csv(deepsnv_sum, csv, row.names=FALSE)

cat(paste("saving coverage to [",output_file_name,".cov.csv].\n", sep=""))
write.csv(cov.df, paste(output_file_name,".cov.csv", sep=""), row.names=FALSE)

cat(paste("saving control coverage to [",output_file_control,".cov.csv].\n", sep=""))
write.csv(cov.con.df, paste(output_file_control,".cov.csv", sep=""), row.names=FALSE)

cat(paste("saving to [",fa,"].\n",sep=""))
write(as.character(consensus_fa),file=fa)

cat("done.\n")
