module load Java singularity nextflow
nohup nextflow run main.nf --csv ./example.csv -with-singularity /fs1/resources/containers/rnaseqfus_active.sif -with-timeline timeline2.html -with-trace trace2.txt  -resume > out.log

