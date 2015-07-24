rename = {
	doc "runs basic shell script to rename remove everything after the _ in the sample name"
	exec "${SCRIPTS}/rename1.sh input.fastq"
}
// fastqc generates two zipped directory files using a custom naming convention
// input: two *.fastq files
// output: two *_fastqc.zip files
fastqc = {
    doc "Run FASTQC to generate QC metrics for the fastq files"
    output.dir = "01_fastqc"
    output_dir = "01_fastqc"
    produce("${output_dir}/*_fastqc.zip") {
        exec "fastqc -o ${output_dir} --noextract -f fastq $input1"
        exec "fastqc -o ${output_dir} --noextract -f fastq $input2"
    }
}



// pydmx generates a directory of demultiplexed fastqs named using the specified barcode file
//     it requires python2.7
// input: two *.fastq files and a bars.csv file
// intermediate output:
//        trimmed: two fastq files derived from original fastq pairs trimmed to shortest sequence in the pair
//        multiplexed: two fastqs with control sequences removed and duplication consolidated along with a summary text file
//   final output:
//        demultiplexed: a pair of fastq files for each barcoded sample (defined by specific barcode file)
pydmx = {
    doc "Runs pydmx to reformat the fastq header and generate consensus sequence"
    out_dir = "02-pydmx/demultiplexed"
    output.dir = "02-pydmx/demultiplexed"
    produce("${out_dir}/*.fastq") {
        exec "python2.7 ${PYDMX} -l $input1.fastq -r $input2.fastq -b ${BARS} -o 02-pydmx "
    }
}


// input: 'r1' and 'r2' fastq files
// output: a *.sam file
//bowtie2_M= {
//    doc "Aligns using Bowtie, generating a SAM file.  Note, this file may be very large."
//    output.dir = "03_align"
//    produce ("03_align/*.sam") {
//        exec "bowtie2 --sensitive -x ${REFERENCE} -1 $input1 -2 $input2 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.fastq")[0] + '.sam'
//    }
//}
bowtie2 = {
    doc "Aligns using Bowtie, generating a SAM file.  Note, this file may be very large."
    output.dir = "03_align"
    if(input.input.size == 2){
        produce ("03_align/*.sam") {
            exec "bowtie2 --sensitive -x ${REFERENCE} -1 $input1 -2 $input2 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + '.sam'
        }
    }
    if(input.input.size == 4){
        produce ("03_align/*.sam") {
            exec "bowtie2 --sensitive -x ${REFERENCE} -1 $input1,$input2 -2 $input3,$input4 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + '.sam'
        }
    }
    if(input.input.size == 6){
        produce ("03_align/*.sam") {
            exec "bowtie2 --sensitive -x ${REFERENCE} -1 $input1,$input2,$input3 -2 $input4,$input5,$input6 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + '.sam'
        }
    }
    if(input.input.size == 8){
        produce ("03_align/*.sam") {
            exec "bowtie2 --sensitive -x ${REFERENCE} -1 $input1,$input2,$input3,$input4 -2 $input5,$input6,$input7,$input8 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + '.sam'
        }
    }
}
samtools_filter = {
    doc "Filter the sam file to include only proper pair alignments"
    output.dir = "03_align"
    filter("proper"){
        exec " samtools view -f 0X0002 -Sh $input.sam > $output.sam"
    }
}

picard_sortsam = {
    doc "Sort SAM file so that its in reference order and convert to BAM."
    tmp_dir    = "./tmp"
    output.dir = "03_align"
    transform("bam") {
        exec """
            java -Xmx4g -Djava.io.tmpdir=$tmp_dir -jar ${LIBRARY_LOCATION}/picard-tools-1.133/picard.jar SortSam
            SO=coordinate
            INPUT=$input.sam
            OUTPUT=$output
            VALIDATION_STRINGENCY=LENIENT
            CREATE_INDEX=true
        """
    }
}


picard_removedups = {
    doc "Remove duplicates"
    tmp_dir    = "./tmp"
    output.dir = "04_removed_duplicates"
    filter("removed") {
        exec """
            java -Xmx2g -Djava.io.tmpdir=$tmp_dir -jar ${LIBRARY_LOCATION}/picard-tools-1.133/picard.jar MarkDuplicates
            INPUT=$input.bam
            OUTPUT=$output
            REMOVE_DUPLICATES=true
            CREATE_INDEX=true
            METRICS_FILE=${output}-picard.out.metrics
            VALIDATION_STRINGENCY=LENIENT
        """
    }
}

picard_markdups = {
    doc "Mark  duplicates"
    tmp_dir    = "./tmp"
    output.dir = "04_marked_duplicates"
    filter("marked") {
        exec """
            java -Xmx1g -Djava.io.tmpdir=$tmp_dir -jar ${LIBRARY_LOCATION}/picard-tools-1.115/MarkDuplicates.jar
            INPUT=$input.bam
            OUTPUT=$output
            REMOVE_DUPLICATES=false
            CREATE_INDEX=true
            METRICS_FILE=${output}-picard.out.metrics
            VALIDATION_STRINGENCY=LENIENT
        """
    }
}

samtools_mpileup = {
	doc "mpileup to create the files needed to messure coverage across the samples"
	output.dir = "05_coverage"
	transform(".pileup"){

	exec "samtools mpileup -d 1000000 $input.bam > $output.pileup" // + new File(input1).name.split("\\.[12]\\.bam")[0] + '.pileup'

	}

}


coverage = {
	doc " this trims large .pileup file to a smaller file that is more easily transfered and presents coverage intuitively"
	output.dir = "05_coverage"
	produce("*.csv"){
			exec "${SCRIPTS}/Trim_to_coverage.py $input.pileup $output.cov.csv"
	}
}



//Get the name of the control file. to use in deepSNV
get_control={
new File('.').eachFileRecurse{
	if(it.name=~/.*$CONTROL.*\..*\.bam$/){
			CONTROL_BAM=it.getPath()
			println "found control $CONTROL_BAM"
		}
	}
}

deepsnv = {
	doc "Runs a basic deepSNV script to call variants in each sample saving the outputs as .Rdata files and csv of the summary output"
	output.dir = "variants"
	def in_bam=file(input).name
	def control = file(CONTROL_BAM).name.replace(".bam","")
	def test = file(input.bam).name.replace(".bam","")
	println "test:" + test
	println "control:" + control
	if( test!=control) {
				produce("variants/*.csv","variants/*.fa"){
					exec "Rscript  ${SCRIPTS}/deepSNV.R ${LIBRARY_LOCATION}/R ${REFERENCE_FA} $input1 $CONTROL_BAM bonferroni"
				}
	} else {
		produce("variants/*.csv"){
			exec "touch ${output}.csv"
		}
	}
}

mapq_conditional = {
        doc "runs a python script that for each variant called calculates the average mapped quality of the variant and adds the whole summary line + this value to the summary.csv file "
        output.dir="mapq"
        //I've been having some trouble getting the right inputs and outputs to line up. The best way to do this might be to make a map that explicitly describes how the sample pair up. bUt that might take a little too long to figure out how to do now.
        def csv = file(input.csv).name.replace(".csv","")
        def bam = file(input.bam).name.replace(".bam","")
        def control = file(CONTROL_BAM).name.replace(".bam","")
	if (csv!=control){
		if ( csv == bam) {
                	println "Found match: " + csv + ".csv and " +bam+".bam"
                	transform("mapq.sum.csv","mapq.reads.csv"){
                        	exec " python ${SCRIPTS}/mapq.py $input.csv $input.bam $output1.csv $output2.csv"
        			}
        	} else {
                	println "csv: " + csv + " doesn't match bam: "+bam
        	}
          }
}


combine = {
	doc "combines the coverage, reads and variant calls into one file that can easily be imported into R for analysis"

	exec "python ${SCRIPTS}/combine.py ./mapq reads.csv all.reads.csv "
	exec "python ${SCRIPTS}/combine.py ./mapq sum.csv all.sum.csv"
	exec "python ${SCRIPTS}/combine_add_name.py ./05_coverage/ .csv all.cov.csv"
}
