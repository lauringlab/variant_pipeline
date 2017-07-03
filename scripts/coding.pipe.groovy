
load './fasta.config.groovy'

coding = {
	doc "Trim the parsed sequences to just the coding regions defined by a separate fasta file which contains the coding region for each segment."
	output.dir = "${MAIN_DIR}/coding_fa"
	filter("coding"){
		exec "python ${SCRIPTS}/trim_to_coding.py ${MUSCLE} $input.fasta ${REFERENCE_FA} -out_fa $output.fasta"
	}

}

concatenate = {
	doc "Concatenate all the segements in the file."
	output.dir = "${MAIN_DIR}/concat_coding"

	filter("concat"){
	exec "python ${SCRIPTS}/concat_all.py $input $output"
	}
}

clean_up = {
	doc "clean up the files that were put in the working directory"
	exec "mv ./fasta.config.groovy ${MAIN_DIR}/fasta.config.groovy"
	//exec "mv ./commandlog.txt ${MAIN_DIR}/consensuscommandlog.txt"
}

run{
	"%.fasta"*[coding] + 
	"%.fasta"*[concatenate] + clean_up

}
