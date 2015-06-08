// Executable bpipe script for flu workflow.

load 'variantPipeline.bpipe.config.groovy'

// All the core pipeline stages in the pipeline
load 'variantPipeline.bpipe.stages.groovy'

done = {
    exec "date"
    exec "echo variantPipeline SUCCESS"
}


/////// Notifications ####

//notifications {
//  gtalk {
//    to="mccronejt@gmail.com"
//    username="bpipe.notifications@gmail.com"
//    password="bpipe.notificatio"
//    events="STAGE_FAILED"
//  }
//}




run  {

   // "%.fastq" * [ rename ] +
    // Align each pair of input files separately in parallel
   // "%.*.fastq" * [ fastqc ]  + needs if statements to handle multiple fastq.  More than 2 
    //pydmx + 
    "%.*.fastq" * [ bowtie2 ] + 
    "%.sam" * [ picard_sortsam + picard_removedups  ] + get_control + 
    "%.bam" * [	samtools_mpileup + trim_pileup,deepsnv + "%.csv"*[mapq_conditional]] + 
    combine +
     done
}
