load 'download.config.groovy'

def branchFile = new File(BRANCHES)
def branchMap = [:]
def sraFastqMap = [:]
branchFile.eachLine { line ->
    def (name, sra, firstFastq, secondFastq) = line.split(/\t/)
    branchMap[name] = [sra, firstFastq, secondFastq]
    sraFastqMap[sra+"_1.fastq"] = firstFastq
    sraFastqMap[sra+"_2.fastq"] = secondFastq
}






def names = []
new File(BRANCHES).eachLine{ line ->
  names<< line.split("\t")[0]
}


download = {
	doc="Download SRR files"
	
	branch.sra = branchMap[branch.name][0]
	branch.sra_file = sra +".sra"
	branch.sra_partial = sra[0..5]
	output.dir = OUT
	SAMPLE_NAME=sra
	produce(SAMPLE_NAME+".sra"){
		exec "wget -r -np -nd -k -P $output.dir ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$sra_partial/$sra/$sra_file"
	}
}


split = {
	doc = "Split sra into fastq"
	branch.sra_file = input.split("/").last()
	sra_name = sra_file.split(".sra")[0]
	output.dir = OUT
	produce(sra_name+"_1.fastq", sra_name+"_2.fastq"){
		exec "fastq-dump  $input.sra --split-files -O $output.dir"
	}
}

rename = {
	doc = "replace the naming scheme with the original"
	output.dir = OUT
	branch.sra_fastq_file = input.split("/").last()
	branch.original = sraFastqMap[sra_fastq_file]
	produce(original){
		exec "mv $input $output"
	}
}	

run{
	names * [download] +
	"%.sra" * [split]+
	"%.fastq" * [rename]
}
