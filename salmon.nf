#!/usr/bin/env nextflow

params.reads = "/data/NextSeq1/190829_NB501697_0156_AH35YWBGXC/Data/Intensities/BaseCalls/ALL354A187_122-60853_S41_R{1,2}_001.fastq.gz"


Channel
        .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .into {read_files_salmon}


params.salmon = false
params.name = false
decoys_file ="/data/bnf/dev/sima/rnaSeq_fus/data/salmon/decoys.txt"
params.ref_salmon = "/data/bnf/dev/sima/rnaSeq_fus/data/salmon/gentrome.fa"
params.outdir = "/data/bnf/dev/sima/rnaSeq_fus/results"



ref_file_salmon = Channel 
	.fromPath(params.ref_salmon)
	.ifEmpty { exit 1, "Star-Fusion reference directory not found!" }	

process create_refIndex{
	cpus = 8
	publishDir "${params.outdir}/tools/salmon", mode:'copy'
	input:
	file ref from ref_file_salmon 
        file decoys from decoys_file              

	output:
	file 'ref_index' into ref_index_ch

	script:
	"""
	salmon  index  --threads 8  -t $ref -d $decoys -i ref_index
	"""  
	}
process quant{
	tag "$sample_id"
	publishDir "${params.outdir}/tools/salmon/transcripts_quant/$sample_id", mode:'copy'
	cpus = 8
	input:
	set val(sample_id), file (reads) from read_files_salmon
	file (index) from ref_index_ch	

	output:
	file  'transcripts_quant'  into transcripts_quant_ch
	script:
	"""
	salmon quant --threads $task.cpus -i $index -l A -1 ${reads[0]} -2 ${reads[1]} --validateMappings -o transcripts_quant
	"""

}
