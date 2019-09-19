#!/usr/bin/dev nextflow
//This pipeline contains Starfusion and fusioncatcher tools. 



jaffa_file = file("/opt/conda/envs/CMD-RNASEQFUS/share/jaffa-1.09-1/JAFFA_direct.groovy")
params.fusionCatcher_ref= "/data/bnf/dev/sima/rnaSeq_fus/data/fusioncatcher/human_v95"
params.star_fusion_ref = "/data/bnf/dev/sima/rnaSeq_fus/data/starFusion/ctat_genome_lib_build_dir"
params.outdir = "/data/bnf/dev/sima/rnaSeq_fus/results"
params.name = false
params.singleEnd= false
params.genome_fasta = "/data/bnf/dev/sima/rnaSeq_fus/data/hg_files/hg38/hg38.fa"
params.genome_gtf = "/data/bnf/dev/sima/rnaSeq_fus/data/hg_files/hg38/gencode.v31.chr_patch_hapl_scaff.annotation.gtf"
params.genome_star_index_dir = "/data/bnf/ref/b37/star"
params.reads = "/data/NextSeq1/190829_NB501697_0156_AH35YWBGXC/Data/Intensities/BaseCalls/ALL354A187_122-60853_S41_R{1,2}_001.fastq.gz"
params.fus = false
params.stra_fusion = false
params.fusioncatcher = false
params.quant= false



Channel
        .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .into {read_files_star_fusion; read_files_fusioncatcher; read_files_jaffa; read_files_star; read_files_star_align }
	
star_fusion_ref = Channel
            .fromPath(params.star_fusion_ref)
            .ifEmpty { exit 1, "Star-Fusion reference directory not found!" }

fusionCatcher_ref = Channel
			.fromPath(params.fusionCatcher_ref)
			.ifEmpty { exit 1, "Fusioncatcher reference directory not found!" }

genome_fasta = Channel.fromPath(params.genome_fasta)

Channel
	.fromPath(params.genome_gtf)
	.into{gtf_indexing; gtf_alignment}

genome_star_index_dir = Channel.fromPath(params.genome_star_index_dir)



process build_star_index {

	publishDir "${params.outdir}/star_index", mode:'copy'
	cpu = 4

	input :

	file genome_fasta 
	file gtf_indexing 


	output:
	file "genome_index" into star_index
	
	script:
	"""
	mkdir genome_index
	STAR \\
	--runMode genomeGenerate \\
	--runThreadN ${task.cpu} \\
	--sjdbGTFfile $gtf_indexing \\
	--genomeDir genome_index/ \\
	--genomeFastaFiles ${genome_fasta} 
	"""	
}



process star_alignment{
	
	tag "$name"
	publishDir "${params.outdir}/tools/star/$name", mode :'copy'
	cpus = 8
	input:
	set val(name), file (reads) from read_files_star_align
	file (index) from  star_index
	file  gtf_alignment
	
	output:
	set file("*Log.final.out"), file ('*.bam') into star_aligned
	file "*.out" into alignment_logs
	file "*SJ.out.tab"
	file "*Log.out" into star_log
	file "*Aligned.sortedByCoord.out.bam" into bam_index_1, bam_index_2
	
	script: 

	"""
	STAR --genomeDir $index \\
	--sjdbGTFfile $gtf_alignment \\
	--readFilesIn $reads \\
	--runThreadN ${task.cpus} \\
	--twopassMode Basic \\
	--outSAMtype BAM SortedByCoordinate  \\
	--readFilesCommand zcat --outFileNamePrefix '$name'
	"""
}

process SamBamBa {
	tag "$name"
	publishDir "${params.outdir}/tools/star/sambamba", mode: 'copy'

	input:
	file reads_bam from bam_index_1
	
	output:
	file "*Aligned.sortedByCoord.out.bam.bai"  into bai_ch
	

	script:
	"""
	sambamba index --show-progress -t 8 $reads_bam 
	"""
}

/*


process star_fusion{
    tag "$name"
    cpus 4  
    publishDir "${params.outdir}/tools/StarFusion/${name}", mode: 'copy'

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
        --output_dir ./${name} 
    """
}



process fusioncatcher {
    tag "$name"
    cpus 4  
    publishDir "${params.outdir}/tools/Fusioncatcher/$name", mode: 'copy'

    //when: params.fusioncatcher || (params.fusioncatcher && params.debug)

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
    tag "$name"
    publishDir  "${params.outdir}/tools/jaffa/${name}", mode: 'copy'

    input:
    set val(name), file(reads) from  read_files_jaffa
    //file groovy from ch_jaffa_direktgroovy
    output:
    file '*.{fasta,csv}' into jaffa_output
    
    script:
    """
    bpipe run -p  genome=hg38 -p refBase="/opt/conda/envs/CMD-RNASEQFUS/share/jaffa-1.09-1/"  $jaffa_file  ${reads[0]} ${reads[1]}  
    """
}

*/
