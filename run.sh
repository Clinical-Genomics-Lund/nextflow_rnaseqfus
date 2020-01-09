#nextflow run main.nf -c ./configs/nextflow.trannel.config --csv ./example.csv -with-singularity container/rnaseqfus_current.sif -resume
nextflow run main.nf -c ./configs/nextflow.trannel.config --csv ./example.csv -with-singularity container/star-fusion.v1.8.1.simg  -resume

