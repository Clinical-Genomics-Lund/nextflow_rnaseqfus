nohup nextflow run main.nf -c ./configs/nextflow.trannel.config --csv ./RnaSeqValidationSamples_filt2.csv -with-singularity container/rnaseqfus_active.sif -with-timeline timeline.html  -resume > run.log &
#nextflow run main.nf -c ./configs/nextflow.trannel.config --csv ./example.csv -with-singularity container/star-fusion.v1.8.1.simg  -resume

