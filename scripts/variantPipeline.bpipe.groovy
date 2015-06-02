// Executable bpipe script for flu workflow.

load 'variantPipeline.bpipe.config.groovy'

// All the core pipeline stages in the pipeline
load 'variantPipeline.bpipe.stages.groovy'

done = {
    exec "date"
    exec "echo variantPipeline SUCCESS"
}


run  {

   // "%.fastq" * [ rename ] +
    // Align each pair of input files separately in parallel
    "%.*.fastq" * [ fastqc ]  +
    //pydmx + 
    "%.*.fastq" * [ bowtie2 ] + 
    "%.sam" * [ picard_sortsam + picard_removedups  ] + get_control + 
    "%.bam" * [	samtools_mpileup + trim_pileup,deepsnv + "%.csv"*[mapq_conditional]] + 
    done
}
