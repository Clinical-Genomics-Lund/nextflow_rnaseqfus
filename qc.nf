#!/usr/bin/env nextflow 

params.outdir = "/data/bnf/dev/sima/rnaSeq_fus/results" 
params.reads = "/data/NextSeq1/190829_NB501697_0156_AH35YWBGXC/Data/Intensities/BaseCalls/ALL354A187_122-60853_S41_R{1,2}_001.fastq.gz"
params.config = "/data/bnf/dev/sima/test/fa_screen/FastQ_Screen_Genomes/fastq_screen.conf"
config_ch = Channel.fromPath(params.config)
params.bam = "/data/bnf/dev/sima/rnaSeq_fus/results/tools/star/ALL354A187_122-60853_S41_R/ALL354A187_122-60853_S41_RAligned.sortedByCoord.out.bam"
params.ref_bed ="/data/bnf/dev/sima/rnaSeq_fus/data/RseQC/hg38_RefSeq.bed"

ref_RseQC_ch= Channel.fromPath(params.ref_bed) 

Channel
	.fromPath(params.bam)
	.into{bam_ch1; bam_ch2}



params.qc = true
params.fastqscreen = true
params.qualimap = true
params.genome_gtf = "/data/bnf/dev/sima/rnaSeq_fus/data/hg_files/hg38/gencode.v31.chr_patch_hapl_scaff.annotation.gtf"


Channel
	.fromPath(params.genome_gtf)
	.into{gtf_quali}


Channel 
	.fromFilePairs( params.reads )
	.ifEmpty{ error "Cannot find any reads matching: ${params.reads}" }
	.into {read_files_star_fusion; read_files_fusioncatcher; read_files_jaffa; read_files_star; read_files_star_align; read_files_fastqscreen}


 
process fastqscreen{ 
	cpus = 8
	publishDir "${params.outdir}/tools/fastqscreen/$name" , mode :'copy'
	tag "$name"
	
	input:
	set val(name), file(reads) from read_files_fastqscreen
	file config from config_ch
	output:
	
	file '*.{html,png,txt}' into fastq_screen_ch

	when:
	params.qc && params.fastqscreen 
	 
	script:
	"""
	fastq_screen --conf $config --aligner bowtie2  ${reads[0]}   
	"""
}


process qualimap {
	tag "$name"
	publishDir "${params.outdir}/tools/qualimap", mode :'copy'

	
	input:
	file(bam1) from bam_ch1
	file gtf_qualimap from gtf_quali

	output:
	file '*' into qualimap_ch

	script:
	"""
	export JAVA_OPTS='-Djava.io.tmpdir=/data/tmp'
	qualimap --java-mem-size=12G rnaseq -bam ${bam1} -gtf $gtf_qualimap  -pe -outdir . 
	"""

}


	
process rseqc_genebody_coverage {
	publishDir "${params.outdir}/tools/genebody_cov", mode:'copy'
	input :
	file ref_bed from  ref_RseQC_ch
	file bam2 from bam_ch2
	output:
	file '*.pdf' into gene_bodyCov_ch
	script:
	"""
	geneBody_coverage.py -i $bam2 -r $ref_bed  -o output
	"""
}



	
