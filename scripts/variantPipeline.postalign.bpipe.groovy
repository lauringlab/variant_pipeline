// Executable bpipe script for flu workflow.

load 'variantPipeline.bpipe.config.groovy'

// All the core pipeline stages in the pipeline
load 'variantPipeline.bpipe.stages.groovy'

done = {
    exec "date"
    exec "echo variantPipeline SUCCESS"
}


run  {
    get_control.post_align +
    "%.bam" * [deepsnv + "%.csv"*[mapq_conditional]] +
    combine +
     done
}
