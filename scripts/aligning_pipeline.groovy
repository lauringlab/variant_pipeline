

// All the core pipeline stages in the pipeline
load 'variantPipeline.bpipe.stages.groovy'

done = {
    exec "date"
    exec "echo variantPipeline SUCCESS"
}


//For Hiseq data rename prior to running using the sample sheet and the change_hiseq_names.py found in scripts
//For Miseq data where the sample name is derived from what we give the core use change_miseq_names.py
Bpipe.run  {

// QC on reads
    "%.*.fastq" * [cutadapt] +
    "%.*.trimmed.fastq" * [ fastqc ] + //  needs if statements to handle multiple fastq.  More than 2
// Align each pair of input files separately in parallel
    "%.*.trimmed.fastq" * [ bowtie2 ] +
// Sort, remove duplicates, call variants, and set coverage
    "%.sam" * [ picard_sortsam + picard_removedups  ] +
    done
    }
