# Rscript --vanilla --slave ../../lauring-variant-pipeline/bin/deepSNV.r ../../lauring-variant-pipeline/lib/ ../04-mark_duplicates/001-0_2a.marked.bam ../04-mark_duplicates/011-PR8control_a.marked.bam ../flu_regions.csv
#
suppressMessages(library("tools"))
suppressMessages(library("plyr"))
suppressMessages(library("deepSNV"))
suppressMessages(library("reshape2"))
suppressMessages(library("magrittr"))
suppressMessages(library("tidyverse"))

#set seed to make distribiutions determinsitic
set.seed(42)

args <- commandArgs(TRUE)
if (length(args) != 10) {
    stop(paste("Usage:", "deepSNV.R" ," {reference.fasta} {test.bam} {control.bam} {c(BH,bonferroni)} {p.val.cut} {c(fisher,average,max)}  {c(two.sided,one.sided,bin)} {stringent_freq} output.csv output.fa",sep=""), call.=FALSE)
 }

################# functions ###########
get_counts <- function(x,deepx){
  pos<-x$pos
  mat_pos<-x$row
  test(deepx)->counts
  base = x$var
  n.tst.fw = counts[mat_pos,base]
  n.tst.bw = counts[mat_pos,tolower(base)]
  cov.tst.fw = sum(counts[mat_pos,c(1:5)])
  cov.tst.bw = sum(counts[mat_pos,c(6:10)])
  
  # now for the control
  control(deepx)->counts.con
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


#print(args)
#library.location <- args[1]
reference.fasta <- args[1]
test.bam <- args[2]
control.bam <- args[3]
method<-args[4]
p.cut<-args[5] # the p.value cut off for samples to be included in downstream analysis.
p.com.meth<-args[6] # combine method for strands
disp<-args[7]
stringent_freq<-args[8]
csv<-args[9]
fa<-args[10]

test_file_prefix = basename(file_path_sans_ext(test.bam))
control_file_prefix = basename(file_path_sans_ext(control.bam))
output_file_name = paste0("deepSNV/",test_file_prefix)
sample_name=strsplit(test_file_prefix,".",fixed=T)[[1]][1]
control_name=strsplit(control_file_prefix,".",fixed=T)[[1]][1]
output_file_control=paste0("deepSNV/",control_file_prefix)

print(paste0("test is :",sample_name ))

print(paste0("control is:",control_name))



cat(paste("loading regions from [", reference.fasta, "]...\n", sep=""))
segments <- fasta.info(reference.fasta)
regions.bed <- data.frame(chr = gsub("[ ].*","", names(segments)), start=1, stop=segments, row.names=NULL)
cat(paste("loaded regions: ", paste(regions.bed$chr, collapse=","),"\n"))


cat("calling variants with deepSNV\n")
cat(paste("\ttest [",test.bam,"]\n\tcontrol [",control.bam,"]...\n", sep=""))
deepsnv.result <- deepSNV(test=test.bam, control=control.bam, regions=regions.bed,combine.method=p.com.meth,q=30) # This calls the variants using the binomial
if(disp=="two.sided"){
deepsnv.result<-estimateDispersion(deepsnv.result,alternative="two.sided")
}else if (disp=="one.sided"){
deepsnv.result<-estimateDispersion(deepsnv.result) # one sided is the default
}else{
deepsnv.result<-deepsnv.result
}


deepsnv_sum<-summary(deepsnv.result,sig.level=as.numeric(p.cut), adjust.method=method)
# filter to stringent freq
deepsnv_sum<-subset(deepsnv_sum,freq.var<stringent_freq)

## Add bases above stringent freq
RF(test(deepsnv.result),total=T)->frequencies # Get the frequencies of all the bases at each prosition
cbind(coordinates(deepsnv.result),frequencies)->frequencies.df

names(frequencies.df)[names(frequencies.df)=="-"]<-"indel" # for ease of handling below
frequencies.df$row<-1:nrow(frequencies.df)
frequencies.df %>% gather(var,freq.var,A:indel)->frequencies.df # long form
subset(frequencies.df,freq.var>=stringent_freq)->less_stringent

less_stringent %>% adply(1,get_counts,deepsnv.result) -> less_stringent.format # formatted like deepsnv df

## Now we get the reference base from the control matrix ######
RF(control(deepsnv.result),total=T)->con.freq
cbind(coordinates(deepsnv.result),con.freq)->df.con

names(df.con)[ncol(df.con)]<-"indel"
df.con %>% gather(ref,freq.var,A:indel)->con.freqs
subset(con.freqs,freq.var>0.5,select=c(chr,pos,ref))->con.major

less_stringent.final<-join(less_stringent.format,con.major)
for(c in names(deepsnv_sum)){
  if(c %in% names(less_stringent.final)==F){
    less_stringent.final[,c]=NA
  }
}

less_stringent.final<-subset(less_stringent.final,select=-c(row))
deepsnv_sum<-rbind(deepsnv_sum,deepsnv_sum)
deepsnv_sum<-deepsnv_sum[order(deepsnv_sum$chr,deepsnv_sum$pos),]


if(dim(deepsnv_sum)[1]>0){ # if varaints were found
    deepsnv_sum$Id<-sample_name # set the sample name for csv
    deepsnv_sum<-subset(deepsnv_sum,!(var=="-" | ref =="-")) # removes the indels that are below the consensus level
    mutate(deepsnv_sum,mutation=paste0(chr,"_",ref,pos,var))->deepsnv_sum
}else{
   deepsnv_sum<-data.frame("chr"=character(),
			   "pos"=integer(),
			   "ref"=character(),
			   "var" = character(),
			   "p.val" = numeric(),
		           "freq.var" = numeric(),
 			   "sigma2.freq.var"= numeric(),
		           "n.tst.fw" = integer(),
			   "cov.tst.fw" = integer(),
		           "n.tst.bw" = integer(),
		           "cov.tst.bw"= integer(),
	                   "n.ctrl.fw"= integer(),
	                   "cov.ctrl.fw" = integer(),
	                   "n.ctrl.bw" = integer(),
		           "cov.ctrl.bw" = integer(),
		           "raw.p.val" = numeric(),
			   "Id"=character(),
			   "mutation"=character() ) # adding the column names even though no variants were found
} 

##### consensus sequences
consensus_fa<-consensusSequence(test(deepsnv.result,total=T),vector=F,haploid=T)
control_fa<-consensusSequence(control(deepsnv.result,total=T),vector=F,haploid=T)




## Coverage ##

cov<-rowSums(test(deepsnv.result,total=T)[,1:4]) # no deletions
#head(cov)
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

#head(cov)

cov.df$Id<-sample_name # set the sample name for csv

### Control coverage

cov.con<-rowSums(control(deepsnv.result,total=T)[,1:4]) # no deletions

# make coverage data.frame

cov.con.df=data.frame(coverage=cov.con,concat.pos=1:length(cov.con))

ddply(cov.con.df,~concat.pos,summarize, coverage=coverage,concat.pos=concat.pos,chr= as.character(regions.bed$chr[max(which(prior.seg.length<concat.pos))]),
      chr.pos=concat.pos-prior.seg.length[max(which(prior.seg.length<concat.pos))])->cov.con.df

cov.con.df$Id<-control_name # set the sample name for csv

#print(head(deepsnv_sum))
cat(paste("saving summary to [",csv,"].\n", sep=""))
write.csv(deepsnv_sum, csv, row.names=FALSE)

#print(head(deepsnv_sum))
cat(paste("saving coverage to [",output_file_name,".cov.csv].\n", sep=""))
write.csv(cov.df, paste(output_file_name,".cov.csv", sep=""), row.names=FALSE)

cat(paste("saving control coverage to [",output_file_control,".cov.csv].\n", sep=""))
write.csv(cov.con.df, paste(output_file_control,".cov.csv", sep=""), row.names=FALSE)

cat(paste("saving to [",fa,"].\n",sep=""))
#save(list=consensus_fa,file=paste(output_file_name,".fa",sep=""))
write(as.character(consensus_fa),file=fa)

#cat(paste("saving to [",output_file_control,".fa].\n",sep=""))
#save(list=consensus_fa,file=paste(output_file_name,".fa",sep=""))
#write(as.character(control_fa),file=paste(output_file_control,".fa",sep=""))

#cat(paste("saving to [",output_file_name,".RData]\n",sep=""))
#eval(parse(text=paste0("snv_",test_file_prefix,"<-deepsnv.result")))
#save(list=ls(pattern=test_file_prefix),file=paste(output_file_name,".RData",sep=""))

cat("done.\n")
