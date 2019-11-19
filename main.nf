#!/usr/bin/dev nextflow



smpl_id = 'ALL354A185_122-59112_S5_R'
jaffa_file = "/opt/conda/envs/CMD-RNASEQFUS/share/jaffa-1.09-2/JAFFA_direct.groovy"

/*
//Input files  
params.outdir = "/data/bnf/dev/sima/rnaSeq_fus/results"
params.reads = "/data/NextSeq1/190808_NB501697_0150_AHN7T7AFXY/Data/Intensities/BaseCalls/ALL354A185_122-59112_S5_R{1,2}_001.fastq.gz"
smpl_id = 'ALL354A185_122-59112_S5_R'

params.genome_fasta = "/data/bnf/dev/sima/rnaSeq_fus/data/hg_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa" //fasta file from ensembl
params.genome_gtf = "/data/bnf/dev/sima/rnaSeq_fus/data/hg_files/gtf/Homo_sapiens.GRCh38.98.gtf" //fasta file from ensembl
params.ref_genome_dir = "/data/bnf/dev/sima/rnaSeq_fus/results/star_refGenome_index/star_ref_index" //is used for star alignment

// fastqscreen genome config file 
params.genome_conf = "/data/bnf/dev/sima/rnaSeq_fus/data/fastqScreen/FastQ_Screen_Genomes/fastq_screen.conf"

params.fusionCatcher_ref= "/data/bnf/dev/sima/rnaSeq_fus/data/fusioncatcher/human_v98"
params.star_fusion_ref = "/data/bnf/dev/sima/rnaSeq_fus/data/starFusion/ctat_genome_lib_build_dir"

// Quantificatin files
params.ref_salmon = "/data/bnf/dev/sima/rnaSeq_fus/data/transcriptome_ref/Homo_sapiens.GRCh38.cdna.all.fa"

//Provider  files
params.ref_bed = "/data/bnf/sw/provider/HPA_1000G_final_38.bed"
params.ref_bedXy= "/data/bnf/sw/provider/xy_38.bed"


//BodyCov
params.ref_rseqc_bed = "/data/bnf/dev/sima/rnaSeq_fus/data/RseQC/Homo_sapiens.GRCh38.79.bed"

//jaffa
jaffa_file = "/opt/conda/envs/CMD-RNASEQFUS/share/jaffa-1.09-2/JAFFA_direct.groovy"



// Set running tool flags 
// QC tools
params.qc = false
params.star_inedx = false
params.star = false
params.fastqscreen = false
params.fastqscreen_genomes = true //This flag shows that if config file for fastqscreen already exists or not. 
params.qualimap = false
params.bodyCov = false
params.provider =false
params.combine = false

// Fusion identification tools 
params.fusion = false
params.star_fusion = true
params.fusioncatcher = false
params.jaffa = false

// Reads quantification tool 
params.quant = false

//Other flags 
params.singleEnd= false

*/

/* Define channels */


Channel
        .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .into {read_files_star_fusion; read_files_fusioncatcher; read_files_jaffa; read_files_star; read_files_star_align; read_files_salmon; read_files_fastqscreen}
	
Channel
	.fromPath(params.genome_fasta)
	.into{genome_fasta_ch; genomeRef_salmon}
		

genome_index = Channel.fromPath(params.ref_genome_dir )
			.ifEmpty { exit 1, "genome reference directory not found!" }


Channel
	.fromPath(params.genome_gtf)
	.into{gtf_star_index; gtf38_qualimap}


Channel.fromPath(params.ref_rseqc_bed).set{ref_RseQC_ch} 

//Provider channels
Channel.fromPath(params.ref_bed).set{bed_ch}
Channel.fromPath(params.ref_bedXy).set{bedXy_ch}


star_fusion_ref = Channel
            .fromPath(params.star_fusion_ref)
            .ifEmpty { exit 1, "Star-Fusion reference directory not found!" }
			

fusionCatcher_ref = Channel
			.fromPath(params.fusionCatcher_ref)
			.ifEmpty { exit 1, "Fusioncatcher reference directory not found!" }


transcriptome_ref= Channel 
			.fromPath(params.ref_salmon)
			.ifEmpty { exit 1, " Reference file/directory not found!" }


/* Part1: QC */

/*
process build_star_index {

	publishDir "${params.outdir}/star_refGenome_index", mode:'copy'
	cpus = 8
	when:
	params.qc || params.star_inedx
	input :
	file (fasta) from  genome_fasta_ch
	file (gtf) from gtf_star_index 


	output:
	file "star_ref_index" into star_index
	
	script:
	"""
	mkdir star_ref_index
	STAR \\
		--runMode genomeGenerate \\
		--runThreadN ${task.cpus} \\
		--sjdbGTFfile ${gtf} \\
		--genomeDir star_ref_index/ \\
		--genomeFastaFiles ${fasta} 
	"""	
}
*/

process star_alignment{
	
	tag "$name"
	publishDir "${params.outdir}/${name}/star", mode :'copy'
	cpus = 8

	when:
		params.qc || params.star

	input:
		set val(name), file (reads) from read_files_star_align
		file (index_files) from  genome_index  //star_index
	
	output:
		set file("Log.final.out"), file ('*.bam') into star_aligned
		file "SJ.out.tab"
		file "Log.out" into star_log
		file "Aligned.sortedByCoord.out.bam" into aligned_bam, star_sort_bam ,star_sort_bam_1 , star_sort_bam_2
		file "Log.final.out" into star_logFinalOut_ch
		

	script: 

	"""
	STAR --genomeDir ${index_files} \\
	--readFilesIn ${reads} \\
	--runThreadN ${task.cpus} \\
	--outSAMtype BAM SortedByCoordinate \\
	--readFilesCommand zcat \\
	--limitBAMsortRAM 10000000000 
	"""
	
}

process SamBamBa {
	tag "$smpl_id"
	publishDir "${params.outdir}/${smpl_id}/star", mode: 'copy'
	when:
		params.qc || params.star
	input:
		file (reads_bam) from aligned_bam
	
	output:
		file "Aligned.sortedByCoord.out.bam.bai"  into star_sort_bai

	script:
	"""
	sambamba index --show-progress -t 8 $reads_bam 
	"""
}


if (params.fastqscreen_genomes) {
    Channel
        .fromPath(params.genome_conf)
        .ifEmpty { exit 1, "Fastqscreen genome config file not found: ${params.genome_conf}" }
        .set {fastq_screen_config_ch}
	} 
else {
	process fastqscreen_getGenome{ 
		
		publishDir "/data/bnf/dev/sima/rnaSeq_fus/data", mode: 'copy'
		
		output :
		file "*" into output_ch
		file "fastq_screen.conf" into fastq_screen_config_ch
		script:
		"""
		fastq_screen --get_genomes
		"""
	}

}
process fastqscreen { 
	//errorStrategy 'ignore'
	cpus = 8
	publishDir "${params.outdir}/${name}/qc/fastqscreen" , mode :'copy'
	tag "$name"
	when:
	params.fastqscreen || params.qc 

	input:
	set val(name), file(reads) from read_files_fastqscreen
	file (config) from fastq_screen_config_ch
	output:
	file '*.{html,png,txt}' into fastq_screen_ch

	script:
	"""
	fastq_screen --conf $config --aligner bowtie2 --force ${reads[0]} ${reads[1]}
	"""
}


process qualimap {
	tag  "$smpl_id"
	publishDir "${params.outdir}/${smpl_id}/qc/qualimap", mode :'copy'
	errorStrategy 'ignore'

	when :
		params.qc || params.qualimap

	input:
		file (bam_f) from star_sort_bam
		file (gtf_qualimap) from gtf38_qualimap

	output:
		file '*' into qualimap_ch

	script:

	"""
	export JAVA_OPTS='-Djava.io.tmpdir=/data/tmp'
	qualimap --java-mem-size=12G rnaseq -bam ${bam_f} -gtf ${gtf_qualimap} -pe -outdir . 
	"""
}
	
process rseqc_genebody_coverage{
	tag "$smpl_id"
	publishDir "${params.outdir}/${smpl_id}/qc/genebody_cov", mode:'copy'
	errorStrategy 'ignore'

	when:
		params.qc || params.bodyCov
	
	input :
	
		file (ref_bed) from ref_RseQC_ch
		file (bam_f) from star_sort_bam_1
		file (bai_f) from star_sort_bai
		//star_sort_bam_1
	
	output:
		file "${smpl_id}*.pdf" into gene_bodyCov_ch
		file "${smpl_id}.geneBodyCoverage.txt" into gene_bodyCov_ch_
	
	script:
	"""
	geneBody_coverage.py -i ${bam_f} -r ${ref_bed} -o ${smpl_id}
	"""
}


process provider{

	tag "$smpl_id"
	publishDir "${params.outdir}/${smpl_id}/qc/provider" , mode:'copy'
	errorStrategy 'ignore'

	when:
		params.qc || params.provider

	input:
		file (bam_f) from star_sort_bam_2
		file (bed_f) from bed_ch
		file (bedXy_f) from bedXy_ch	
	
	output:
		file "*.genotypes" into provider_output_ch

	script:

	prefix = "${smpl_id}"
	"""
	provider.pl  --out ${prefix} --bed ${bed_f} --bam ${bam_f} --bedxy ${bedXy_f}
	"""
}

	
/* Part2 : fusion identification part */

process star_fusion {
	//errorStrategy 'ignore'
    tag "$name"
    cpus 16  
    publishDir "${params.outdir}/${name}/fusion/StarFusion/", mode: 'copy'

    when:
    params.star_fusion || params.fusion 

    input:
	set val(name), file(reads) from read_files_star_fusion
	file (reference) from star_fusion_ref
	//file (junction) from star_junction_ch

    output:
	file '*fusion_predictions.tsv' optional true into star_fusion_agg_ch
    file '*.{tsv,txt}' into star_fusion_output

    script:
    //def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    option = params.singleEnd ? "--left_fq ${reads[0]}" : "--left_fq ${reads[0]} --right_fq ${reads[1]}"
    //def extra_params = params.star_fusion_opt ? "${params.star_fusion_opt}" : ''
    """
    STAR-Fusion \\
    --genome_lib_dir ${reference} \\
    ${option} \\
    --CPU ${task.cpus} \\
    --output_dir . \\
	--verbose_level 2 \\
	--FusionInspector validate \\
	--tmpdir /data/bnf/tmp
    """
	/*
	"""
	STAR-Fusion \\
    --genome_lib_dir ${reference} \\
	--J ${junction} \\
    --CPU ${task.cpus} \\
    --output_dir . \\
	--FusionInspector validate  \\
	--verbose_level 2 
	"""
 */		
}


process fusioncatcher{
    errorStrategy 'ignore'
    tag "$name"
    cpus 8 
    publishDir "${params.outdir}/${name}/fusion/FusionCatcher", mode: 'copy'

    when: 
	params.fusioncatcher || params.fusion

    input:
    set val(name), file(reads) from read_files_fusioncatcher
    file (data_dir) from fusionCatcher_ref

    output:
	file 'final-list_candidate-fusion-genes.txt' optional true into fusioncatcher_fusions
	file 'final-list_candidate-fusion-genes.hg19.txt' into final_list_fusionCatcher_ch, final_list_fusionCatcher_agg_ch
   	file '*.{txt,zip,log}' into fusioncatcher_output

    script:

    option = params.singleEnd ? reads[0] : "${reads[0]},${reads[1]}"
    //def extra_params = params.fusioncatcher_opt ? "${params.fusioncatcher_opt}" : ''
    """
   	fusioncatcher.py  -d ${data_dir} -i ${option}  --threads ${task.cpus} -o . 
    """
}

process filter_aml_fusions {
	errorStrategy 'ignore'
	publishDir "${params.outdir}/${smpl_id}/fusion/FusionCatcher", mode: 'copy'

	when: 
		params.fusioncatcher || params.fusion
	input:
		file (fusioncatcher_files) from fusioncatcher_output

	output:
	file "${smpl_id}.fusioncatcher.xls" into filter_fusion_ch
	script:
	"""
	filter_aml_fusions.pl ${fusioncatcher_files} > ${smpl_id}.fusioncatcher.xls
	"""
}

process jaffa {
    tag "$name"
	errorStrategy 'ignore'
    publishDir  "${params.outdir}/${name}/fusion/jaffa", mode: 'copy'

    when:
    	params.jaffa || params.fusion

    input:
    	set val(name), file(reads) from  read_files_jaffa

    output:
    	file "*.csv" into jaffa_csv_ch
    	file "*.fasta" into jaffa_fasta_ch 
    
    script:

   	"""
   	bpipe run  -p  genome=hg38 -p refBase="/data/bnf/dev/sima/rnaSeq_fus/data/hg_files/hg38/"  ${jaffa_file}  ${reads[0]} ${reads[1]}  
   	"""
}
 /* 
process Fuseq{
tag "$name"
errorStrategy 'ignore'
publishDir  "${params.outdir}/${name}/fusion/Fuseq", mode: 'copy'

script:
"""
export LD_LIBRARY_PATH=/data/bnf/sw/fuseq/1.1.0/linux/lib:$LD_LIBRARY_PATH
export PATH=/data/bnf/sw/fuseq/1.1.0/linux/bin:$PATH
FuSeq -i /data/bnf/ref/fuseq/GRCh37/ -l IU -1 <(gunzip -c  /data/NextSeq2/180815_NB501699_0062_AH73F5AFXY/Data/Intensities/BaseCalls//2845-15_S2_R1_001.fastq.gz ) -2 <(gunzip -c  /data/NextSeq2/180815_NB501699_0062_AH73F5AFXY/Data/Intensities/BaseCalls//2845-15_S2_R2_001.fastq.gz ) -p  4 -g /data/bnf/ref/rsem/GRCh37/Homo_sapiens.GRCh37.75.gtf -o /data/bnf/premap/rnaseq/INV3-2845-15_0.fuseqfolder
Rscript FuSeq.R  in=/data/bnf/premap/rnaseq/INV3-2845-15_0.fuseq/folder txfasta=/data/bnf/ref/rsem/GRCh37/GRCh37.transcripts.fa sqlite=/data/bnf/ref/fuseq/Homo_sapiens.GRCh37.75.sqlite txanno=/data/bnf/ref/fuseq/GRCh37/Homo_sapiens.GRCh37.75.txAnno.RData \\ params=/data/bnf/sw/fuseq/1.1.0/R/params.txt out=/data/bnf/premap/rnaseq/INV3-2845-15_0.fuseq
 */
process cleanTranscriptomeVersion{
        publishDir "/data/bnf/dev/sima/rnaSeq_fus/data/transcriptome_ref/", mode : 'copy'
		errorStrategy 'ignore'
        input :
        file (tr_fasta) from transcriptome_ref

        output:
        file 'transcripts.cleanversion.fa' into transcriptomeRef_cleanversion_salmon
        script:
        """
        Rscript /opt/fuseq/FuSeq_v1.1.2_linux_x86-64/R/excludeTxVersion.R ${tr_fasta} transcripts.cleanversion.fa
        """
}
/***************************************************/
/*  Part3 : Expression quantification    */
/****************************************************/

process create_refIndex {
	errorStrategy 'ignore'
	cpus = 8
	publishDir "/data/bnf/dev/sima/rnaSeq_fus/data/transcriptome_ref/", mode:'copy'

	when:
		params.quant 

	input:
		file (ref_genome_fasta) from genomeRef_salmon
		file (transcriptome_fasta) from transcriptomeRef_cleanversion_salmon 
    	            

	output:
		file 'salmon_index' into salmon_index_ch

	script:
	"""
	grep "^>" <(zcat ${ref_genome_fasta}) | cut -d " " -f 1 > decoys.txt
	sed -i -e 's/>//g' decoys.txt
	cat ${transcriptome_fasta} ${ref_genome_fasta} > gentrome.fa.gz
	salmon index --threads 12 -t gentrome.fa.gz -d decoys.txt  -i salmon_index --gencode
	"""  
	}

process quant{
	errorStrategy 'ignore'
	tag "$name"
	publishDir "${params.outdir}/${name}", mode:'copy'
	cpus = 8

	when:
		params.quant 

	input:
		set val(name), file (reads) from read_files_salmon
		file (index) from salmon_index_ch

	output:
		file 'quant'  into transcripts_quant_ch
		file 'quant/libParams/flenDist.txt' into flendist_ch
		file 'quant/quant.sf' into quant_ch, quant_ch_extract

	script:

	"""
	salmon quant --threads ${task.cpus} -i ${index} -l A -1 ${reads[0]} -2 ${reads[1]} --validateMappings -o quant

	"""
}


process extract_expression {
	errorStrategy 'ignore'
	publishDir "${params.outdir}/${smpl_id}/quant", mode:'copy'

	input:
	file (quants) from quant_ch_extract

	when:
		params.quant

	output:
	file "*" into salmon_expr_ch

	script:
	"""
	extract_expression_fusion.R  ${quants}  ${smpl_id}.salmon.expr
	"""

}



/* Part 4: Post processing */  

// Create fusion report  
process create_fusion_report{
	errorStrategy 'ignore'
	publishDir "${params.outdir}/${smpl_id}/final_results" , mode:'copy'

	input:
		file (quant) from quant_ch

	output:
		file "${smpl_id}.STAR.fusionreport" into report

	script:
	"""
	fusion_classifier_report.R  ${smpl_id}  ${quant}  ${smpl_id}.STAR.fusionreport
	"""

}


/* post alignment */
process postaln_qc_rna {
	publishDir "${params.outdir}/${smpl_id}/final_results" , mode:'copy'
	errorStrategy 'ignore'

	when:
		params.combine 

	input:
		file (star_final) from star_logFinalOut_ch
		file (fusion) from final_list_fusionCatcher_ch
		file (geneCov) from gene_bodyCov_ch_
		file (provIder) from provider_output_ch
		file (flendist) from flendist_ch
	
	output:
		file "${smpl_id}.STAR.rnaseq_QC" into final_QC 

	
	script:

	"""
	postaln_qc_rna.R  --star ${star_final} --fusion ${fusion} --id '${smpl_id}'  --provider ${provIder} --flendist ${flendist} --genebody ${geneCov}> '${smpl_id}.STAR.rnaseq_QC'
	"""
} 


/* Part 5 :  Prepare for and upload to Coyote */

// aggregate fusion files
process aggregate_fusion{
	errorStrategy 'ignore'
	publishDir "${params.outdir}/${smpl_id}/final_results" , mode: 'copy'

	input:
		file (fusionCatcher_file) from final_list_fusionCatcher_agg_ch
		file (starFusion_file) from star_fusion_agg_ch
		file (fusionJaffa_file) from jaffa_csv_ch

	output:
		file "${smpl_id}.agg.vcf" into agg_vcf_ch

	script:

	"""
	aggregate_fusions.pl \\
		--fusioncatcher ${fusionCatcher_file} \\
		--starfusion ${starFusion_file} \\
		--jaffa ${fusionJaffa_file} \\
		--priority fusioncatcher,jaffa,starfusion > ${smpl_id}.agg.vcf
	"""
}


/*
//Register to CMD 
process register_to_CMD{
	publishDir "${params.outdir}/${smpl_id}/final_results" , mode:'copy'
	input:
	file (final_QC_file) from final_QC
	output:
	file '*' into  registering_ch
	script:
	"""
	register_sample.pl --run-folder  /data/NextSeq1/181121_NB501697_0089_AHFGY3AFXY  --sample-id ${smpl_id} --assay rnaseq-fusion --qc  ${final_QC_file}
	"""
	}
*/

/*
process import_to_Coyote{
	input:
	file (fusion_agg) from agg_vcf_ch
	file 

import_fusion_to_coyote.pl --classification /data/bnf/postmap/rnaseq/6192-11.STAR.fusionreport --fusions /data/bnf/6192-11.agg.vcf --id 6192-11-fusions --qc /data/bnf/postmap/rnaseq/6192-11.STAR.rnaseq_QC --group fusion
*/ 


//env PERL5LIB= PERL_LOCAL_LIB_ROOT= cpan

