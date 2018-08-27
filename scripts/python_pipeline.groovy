

// All the core pipeline stages in the pipeline
load 'variantPipeline.bpipe.stages.groovy'

done = {
    exec "date"
    exec "echo variantPipeline SUCCESS"
}


//For Hiseq data rename prior to running using the sample sheet and the change_hiseq_names.py found in scripts
//For Miseq data where the sample name is derived from what we give the core use change_miseq_names.py
Bpipe.run  {
    "%.bam" * [ consensus ] +
    "%.bam" * [ position_stats]+
    done
    }
