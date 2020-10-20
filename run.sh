nohup nextflow run main.nf -c ./configs/nextflow.trannel.config --csv /path/to/xxx.csv -with-singularity /path/to/container/rnaseqfus_active.sif -with-timeline timeline.html  -resume > run.log &
