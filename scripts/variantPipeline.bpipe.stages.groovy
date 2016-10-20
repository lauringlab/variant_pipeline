// fastqc generates two zipped directory files using a custom naming convention
// input: two *.fastq files
// output: two *_fastqc.zip files
fastqc = {
    doc "Run FASTQC to generate QC metrics for the fastq files"
    output.dir = "01_fastqc"
    output_dir = "01_fastqc"
    if(input.input.size == 2){
    	produce("${output_dir}/*_fastqc.zip") {
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input1"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input2"
    		}
	}
    if(input.input.size == 4){
    	produce("${output_dir}/*_fastqc.zip") {
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input1"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input2"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input3"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input4"
    		}
	}
    if(input.input.size == 6){
    	produce("${output_dir}/*_fastqc.zip") {
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input1"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input2"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input3"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input4"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input5"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input6"
    		}
	}
   if(input.input.size == 8){
    	produce("${output_dir}/*_fastqc.zip") {
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input1"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input2"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input3"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input4"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input5"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input6"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input7"
        	exec "fastqc -o ${output_dir} --noextract -f fastq $input8"
    		}
	}
}

samtools_mapq_filter ={
	doc " remove the reads that fail a mapping quality cut off"
	output.dir = "02_filter"
	filter("filtered"){
		exec "samtools view -Shq 30 $input > $output"
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
//        exec "bowtie2 --seed 42 --sensitive -x ${REFERENCE} -1 $input1 -2 $input2 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.fastq")[0] + 'sam' 2> new File(input1).name.split("\.[12]\.[0-9]\.fastq")[0] + '.log'
//    }
//}
bowtie2 = {
    doc "Aligns using Bowtie, generating a SAM file.  Note, this file may be very large."
    output.dir = "03_align"
    if(input.input.size == 2){
 	def sam_out=file(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0]+ '.sam'
	//println "expected outout " + sam_out
	produce(sam_out) {
            exec "bowtie2 --seed 42 --sensitive -x ${REFERENCE} -1 $input1 -2 $input2 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + 'sam' 2> new File(input1).name.split("\.[12]\.[0-9]\.fastq")[0] + '.log'
        }
    }
    if(input.input.size == 4){
        def sam_out=file(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0]+ '.sam'
        println "expected outout " + sam_out
        produce(sam_out) {
            exec "bowtie2 --seed 42 --sensitive -x ${REFERENCE} -1 $input1,$input2 -2 $input3,$input4 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + 'sam' 2> new File(input1).name.split("\.[12]\.[0-9]\.fastq")[0] + '.log'
        }
    }
    if(input.input.size == 6){
        def sam_out=file(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0]+ '.sam'
        println "expected outout " + sam_out
        produce(sam_out) {
            exec "bowtie2 --seed 42 --sensitive -x ${REFERENCE} -1 $input1,$input2,$input3 -2 $input4,$input5,$input6 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + 'sam' 2> new File(input1).name.split("\.[12]\.[0-9]\.fastq")[0] + '.log'
        }
    }
    if(input.input.size == 8){
       def sam_out=file(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0]+ '.sam'
        println "expected outout " + sam_out
        produce(sam_out) {
            exec "bowtie2 --seed 42 --sensitive -x ${REFERENCE} -1 $input1,$input2,$input3,$input4 -2 $input5,$input6,$input7,$input8 -S ./03_align/" + new File(input1).name.split("\\.[12]\\.[0-9]\\.fastq")[0] + 'sam' 2> new File(input1).name.split("\.[12]\.[0-9]\.fastq")[0] + '.log'
        }
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


//Get the name of the control file. to use in deepSNV
get_control={
new File('04_removed_duplicates').eachFileRecurse{
         if(it.name=~/.*$CONTROL.*bam$/){
			CONTROL_BAM=it.getPath()
			println "found control $CONTROL_BAM"
		}
	}
}
get_control.post_align={
new File(INPUT_DIR).eachFileRecurse{ // I\m not groovy with groovy and may need to play with how this input_dir variable is called
         if(it.name=~/.*$CONTROL.*bam$/){
                        CONTROL_BAM=it.getPath()
                        println "found control $CONTROL_BAM"
                }
        }
}

deepsnv = {
	doc "Runs a basic deepSNV script to call variants in each sample saving the outputs as .Rdata files and csv of the summary output"
	output.dir = "deepSNV"
	def in_bam=file(input).name
	def control = file(CONTROL_BAM).name.replace(".bam","")
	def test = file(input.bam).name.replace(".bam","")
	println "test:" + test
	println "control:" + control
				transform("csv","fasta"){
					exec "Rscript  ${SCRIPTS}/deepSNV.R ${REFERENCE_FA} $input1 $CONTROL_BAM bonferroni ${P_CUT} ${P_COM_METH} ${DISP} $output.csv $output.fasta"
				}

}

mapq_conditional = {
        doc "runs a python script that for each variant called calculates the average mapped quality of the variant and adds the whole summary line + this value to the summary.csv file "
        output.dir="Variants"
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


parse= {
	doc "Take the concatenated consensus fasta file from deepSNV and deconcatenate it using the segmented positions in from the coverage file"
	output.dir = "parsed_fa"
	filter("parsed"){
		exec "python ${SCRIPTS}/parse_consensus.py ${REFERENCE_FA} $input $output "
	}
	forward input.csv
}



classification = { 
        doc "Add amino acid data to variant calls based a reference sequence"
        output.dir = "Final_variants"
        def ref = "./parsed_fa/"+file(CONTROL_BAM).name.replace(".bam",".parsed.fasta")    
	filter("AA"){
                exec "python ${SCRIPTS}/AA_var.py  $ref $input.csv $output ${OPTIONS}"
                }   
}

sift = {
	doc "Variant csv file and filter based on quality scores"
	output.dir = "Filter_var"
	filter("filtered"){
		exec "python ${SCRIPTS}/filter_var.py $input $output ${OPTIONS} "
	}
}



combine = {
	doc "combines the coverage, reads and variant calls into one file that can easily be imported into R for analysis"

	exec "python ${SCRIPTS}/combine.py ./Final_variants AA.csv ../all.variants.csv"
	exec "python ${SCRIPTS}/combine.py ./deepSNV cov.csv ../all.coverage.csv"
}
