#!/usr/bin/dev nextflow
//This pipeline contains Starfusion and fusioncatcher tools. 


//JAFFA= "/opt/conda/envs/CMD-RNASEQFUS/share/jaffa-1.09-1"
jaffa_file = params.file("/opt/conda/envs/CMD-RNASEQFUS/share/jaffa-1.09-1/JAFFA_direct.groovy")
params.fusionCatcher_ref= "/data/bnf/dev/sima/rnaSeq_fus/data/fusioncatcher/human_v95"
params.star_fusion_ref = "/data/bnf/dev/sima/rnaSeq_fus/data/starFusion/ctat_genome_lib_build_dir"
params.outdir = "/data/bnf/dev/sima/rnaSeq_fus/results"
params.name = false
params.singleEnd= false

params.reads = "/data/NextSeq1/190808_NB501697_0150_AHN7T7AFXY/Data/Intensities/BaseCalls/ALL354A185_122-59112_S5_R{1,2}_001.fastq.gz"

Channel
        .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .into {read_files_star_fusion; read_files_fusioncatcher; read_files_jaffa }


star_fusion_ref = Channel
            .fromPath(params.star_fusion_ref)
            .ifEmpty { exit 1, "Star-Fusion reference directory not found!" }

fusionCatcher_ref = Channel
				.fromPath(params.fusionCatcher_ref)
				.ifEmpty { exit 1, "Fusioncatcher reference directory not found!" }


process star_fusion{
    tag "$name"
    cpus 4  
    publishDir "${params.outdir}/tools/StarFusion", mode: 'copy'

    //when:
    //params.star_fusion || (params.star_fusion && params.debug)

    input:
    set val(name), file(reads) from read_files_star_fusion
    //file star_index_star_fusion
    file reference from star_fusion_ref

    output:
    file '*fusion_predictions.tsv' optional true into star_fusion_fusions
    file '*.{tsv,txt}' into star_fusion_output

    script:
    //def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    option = params.singleEnd ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    //def extra_params = params.star_fusion_opt ? "${params.star_fusion_opt}" : ''
    """
    STAR-Fusion \\
        --genome_lib_dir ${reference} \\
        ${option}\\
        --CPU ${task.cpus} \\
        --examine_coding_effect \\
        --output_dir . 
    """
}


process fusioncatcher {
    tag "$name"
    cpus 4  
    publishDir "${params.outdir}/tools/Fusioncatcher", mode: 'copy'

    //when:
    params.fusioncatcher || (params.fusioncatcher && params.debug)

    input:
    set val(name), file(reads) from read_files_fusioncatcher
    file data_dir from fusionCatcher_ref

    output:
    file 'final-list_candidate-fusion-genes.txt' optional true into fusioncatcher_fusions
    file '*.{txt,zip,log}' into fusioncatcher_output

    script:
    option = params.singleEnd ? reads[0] : "${reads[0]},${reads[1]}"
    //def extra_params = params.fusioncatcher_opt ? "${params.fusioncatcher_opt}" : ''
    """
    fusioncatcher \\
        -d ${data_dir} \\
        -i ${option} \\
        --threads ${task.cpus} \\
        -o . \\
    """
}

process jaffa{
    publishDir  "${params.outdir}/tools/jaffa", mode: 'copy'

    input:
    set val(name), file(reads) from  read_files_jaffa
    //file groovy from ch_jaffa_direktgroovy
    output:
    file '*.{fasta,csv}' into jaffa_output
    
    script:
    """
    bpipe run m_  ${reads[0]} ${reads[1]}
    """
}
