#!/usr/bin/dev nextflow

OUTDIR = params.outdir+'/'+params.subdir

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.clarity_sample_id, file(row.read1), file(row.read2)) }
    .into {reads_starfusion; reads_fusioncatcher; reads_jaffa; reads_align; reads_salmon; reads_fastqscreen}

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{row -> tuple(row.clarity_sample_id,row.id,row.clarity_pool_id,row.assay)}
	.set{coyote_meta}

/**************************/
/* Part1: Alignment & QC  */
/**************************/
process star_alignment{

	errorStrategy 'ignore'
	tag "${smpl_id}"
	publishDir "$OUTDIR/bam", mode :'copy'
	cpus = 16
	memory 40.GB
	when:
		params.qc || params.star

	input:
		set val(smpl_id) , file(read1), file(read2) from reads_align
		
	
	output:
		set val(smpl_id), file("${smpl_id}.Aligned.sortedByCoord.out.bam") into qualimap_bam, provider_bam, geneBody_bam
		set val(smpl_id), file("${smpl_id}.Log.final.out") into star_logFinalOut_ch
		//set val(smpl_id), file("${smpl_id}.Aligned.sortedByCoord.out.bam"),file("${smpl_id}.Aligned.sortedByCoord.out.bam.bai") into geneBody_bambai
		//sambamba index --show-progress -t 8 ${smpl_id}.Aligned.sortedByCoord.out.bam
	script: 

	"""
	STAR --genomeDir ${params.ref_genome_dir} \\
		--readFilesIn ${read1} ${read2} \\
		--runThreadN ${task.cpus} \\
		--outSAMtype BAM SortedByCoordinate \\
		--readFilesCommand zcat \\
		--limitBAMsortRAM 10000000000 

	mv Aligned.sortedByCoord.out.bam  ${smpl_id}.Aligned.sortedByCoord.out.bam
	mv Log.final.out ${smpl_id}.Log.final.out
	"""	
	}

process fastqscreen{ 

	errorStrategy 'ignore'
	//scratch true
	tag "${smpl_id}"
	cpus = 8
	publishDir "$OUTDIR/qc/${smpl_id}.fastqscreen" , mode :'copy'
	
	when:
		params.fastqscreen 

	input:
		set val(smpl_id), file(read1), file(read2) from reads_fastqscreen
		

	output:
		file '*.{html,png,txt}' into fastq_screen_ch

	script:

	"""
	fastq_screen --conf ${params.genome_conf} --aligner bowtie2 --force ${read1} ${read2}
	"""
	}


process qualimap{

	tag  "${smpl_id}"
	publishDir "$OUTDIR/qc", mode :'copy'
	errorStrategy 'ignore'
	memory 18.GB

	when :
		params.qualimap 

	input:
		set val(smpl_id), file(bam_f) from qualimap_bam
		
	output:
		file '*' into qualimap_ch

	script:

	"""
	export JAVA_OPTS='-Djava.io.tmpdir=${params.tmp_dir}'
	qualimap --java-mem-size=18G rnaseq -bam ${bam_f} -gtf ${params.genome_gtf} -pe -outdir ${smpl_id}.qualimap
	"""
	}
	
process rseqc_genebody_coverage{
	tag "${smpl_id}"
	publishDir "$OUTDIR/qc", mode:'copy'
	errorStrategy 'ignore'
	memory 10.GB
	cpus 1
	
	when:
		params.qc || params.bodyCov
	
	input :
	
		set val(smpl_id),file(bam_f) from geneBody_bam
		
			
	output:
		
		set val(smpl_id), file("${smpl_id}.geneBodyCoverage.txt") into gene_bodyCov_ch
	
		//samtools view  -s 0.3 -b ${bam_f} -o ${smpl_id}.subsampled.bam
	script:
	"""
	
	samtools view  -s 0.3 -b ${bam_f} -o ${smpl_id}.subsample.bam
	sambamba index --show-progress -t 8 ${smpl_id}.subsample.bam
	geneBody_coverage.py -i ${smpl_id}.subsample.bam -r ${params.ref_rseqc_bed} -o ${smpl_id}
	"""
	}

process provider{

	tag "${smpl_id}"
	publishDir "$OUTDIR/qc" , mode:'copy'
	errorStrategy 'ignore'

	when:
		params.qc || params.provider

	input:
		set val(smpl_id), file(bam_f) from provider_bam
		
	
	output:
		set val(smpl_id), file("${smpl_id}.genotypes") into provider_output_ch

	script:

	prefix = "${smpl_id}"
	"""
	provider.pl  --out ${prefix} --bed ${params.ref_bed} --bam ${bam_f} --bedxy ${params.ref_bedXy}
	"""
	}
	
	
	

/* ******************************** */	
/* Part2 : fusion identification  */
/* ******************************** */

process star_fusion{
	//errorStrategy 'ignore'
	//scratch true
	tag "${smpl_id}"
	cpus 18
	memory  70.GB
	publishDir "$OUTDIR/fusion", mode: 'copy'
	
	when:
		params.star_fusion || params.fusion 
	
	input:
		set val(smpl_id) , file(read1), file(read2) from reads_starfusion
		
	output:
		set val(smpl_id), file("${smpl_id}.star-fusion.fusion_predictions.tsv") optional true into star_fusion_agg_ch
    	
	script:

    	//def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
   	option = params.singleEnd ? "--left_fq ${read1}" : "--left_fq ${read1} --right_fq ${read2}"
    	//def extra_params = params.star_fusion_opt ? "${params.star_fusion_opt}" : ''
    	"""
    	/usr/local/src/STAR-Fusion/STAR-Fusion \\
		--genome_lib_dir ${params.star_fusion_ref} \\
		${option} \\
		--CPU ${task.cpus} \\
		--output_dir . \\
		--FusionInspector validate \\
		--verbose_level 2  
		
	mv  star-fusion.fusion_predictions.tsv ${smpl_id}.star-fusion.fusion_predictions.tsv 
    	"""
	}


process fusioncatcher {
	errorStrategy 'ignore'
	tag '${smpl_id}'
	cpus 16 
	publishDir "$OUTDIR/fusion", mode: 'copy'
	
	when: 
		params.fusioncatcher || params.fusion
	
	input:
		set val(smpl_id) , file(read1), file(read2) from reads_fusioncatcher
    	
	output:
		set val(smpl_id), file("${smpl_id}.final-list_candidate-fusion-genes.hg19.txt") into final_list_fusionCatcher_ch, final_list_fusionCatcher_agg_ch
		set val(smpl_id), file("${smpl_id}.fusioncatcher.xls") into filter_fusion_ch
	
	script:
	option = params.singleEnd ? read1 : "${read1},${read2}"
    	//def extra_params = params.fusioncatcher_opt ? "${params.fusioncatcher_opt}" : ''
    	"""
   	fusioncatcher.py  -d ${params.fusionCatcher_ref} -i ${option}  --threads ${task.cpus} -o ./${smpl_id}.fusioncatcher
	filter_aml_fusions.pl ./${smpl_id}.fusioncatcher > ${smpl_id}.fusioncatcher.xls
	mv  ./${smpl_id}.fusioncatcher/final-list_candidate-fusion-genes.hg19.txt ${smpl_id}.final-list_candidate-fusion-genes.hg19.txt
    	"""
	}


process jaffa{
	tag "${smpl_id}"
	errorStrategy 'ignore'
	publishDir  "$OUTDIR/fusion", mode: 'copy'
	memory 75.GB 
	cpus 18
	
	when:
		params.jaffa || params.fusion
	
	input:
		set val(smpl_id) , file(read1), file(read2) from  reads_jaffa
	
	output:
		set val(smpl_id) ,file ("${smpl_id}.jaffa_results.csv") into jaffa_csv_ch
		file "*.fasta" into jaffa_fasta_ch 
	
    	script:

   	"""
   	bpipe run -m 64GB -n ${task.cpus} -p genome=hg38 -p refBase="${params.jaffa_base}" ${params.jaffa_file}  ${read1} ${read2}
	mv  jaffa_results.csv ${smpl_id}.jaffa_results.csv
   	"""
	}




/***************************************************/
/*  Part3 : Expression quantification    */
/****************************************************/

process salmon{

	errorStrategy 'ignore'
	tag "${smpl_id}"
	publishDir "$OUTDIR/quant", mode:'copy'
	cpus = 8

	when:
		params.quant 

	input:
		set val(smpl_id) , file(read1), file(read2) from reads_salmon
		
		
	output:
		set val(smpl_id), file ("${smpl_id}.quant.sf") into quant_ch
		set val(smpl_id), file ("${smpl_id}.flenDist.txt") into flendist_ch
		
	script:

	"""
	salmon quant --threads ${task.cpus} -i ${params.salmon_index_dir} -l A -1 ${read1} -2 ${read2} --validateMappings -o .
	mv ./libParams/flenDist.txt ${smpl_id}.flenDist.txt
	mv quant.sf ${smpl_id}.quant.sf
	
	"""

	}

process create_expr_ref {

	publishDir "${params.refbase}/extract_expr_ref", mode:'copy'
	when:
		params.create_exprRef
	output:
		file("reference_expression.all.tsv")
		file("genes_of_interest.tsv")
	script:
	"""
	extract_expression_fusion_ny.R create-reference
	"""
}

process extract_expression {
	
	errorStrategy 'ignore'
	publishDir "$OUTDIR/finalResults", mode:'copy'

	input:
		set val(smpl_id), file(quants) from quant_ch // args[2] in both
		
	when:
		params.quant

	output:

		set val(smpl_id), file("${smpl_id}.salmon.expr") into salmon_expr_ch  
		set val(smpl_id), file("${smpl_id}.expr.classified") into classification_report //args[5] in classifier

	script:

	"""
	extract_expression_fusion_ny.R  ${params.genesOfIntrest}  ${quants}  ${params.reference_expression_all}  ${smpl_id}.salmon.expr
	fusion_classifier_report_ny.R  ${smpl_id} ${quants} ${params.hem_classifier_salmon} ${params.ensembl_annotation} ${smpl_id}.expr.classified
	
	"""
	}


/***************************************************/
/*  Part 4: Post processing                         */
/****************************************************/
 
process postaln_qc_rna {
	publishDir "$OUTDIR/finalResults" , mode:'copy'
	errorStrategy 'ignore'

	when:
		params.combine 

	input:
		set val(smpl_id), file(star_final), file(fusion), file(geneCov), file(provIder), file(flendist) from star_logFinalOut_ch.join(final_list_fusionCatcher_ch.join(gene_bodyCov_ch.join(provider_output_ch.join(flendist_ch))))
	
	
	output:
		set val(smpl_id), file("${smpl_id}.STAR.rnaseq_QC") into final_QC 

	
	script:

	"""
	postaln_qc_rna.R \\
		--star ${star_final} \\
		--fusion ${fusion} \\
		--id '${smpl_id}' \\
		--provider ${provIder} \\
		--flendist ${flendist} \\
		--genebody ${geneCov}> '${smpl_id}.STAR.rnaseq_QC'
	"""
	} 



/***********************************************/
/* Part 5 :  Prepare for and upload to Coyote  */
/***********************************************/
// aggregate fusion files
process aggregate_fusion{
	errorStrategy 'ignore'
	publishDir "$OUTDIR/finalResults" , mode: 'copy'

	when :
	     params.combine

	input:
		set val(smpl_id), file(fusionCatcher_file), file(starFusion_file), file(fusionJaffa_file) from final_list_fusionCatcher_agg_ch.join(star_fusion_agg_ch.join(jaffa_csv_ch)).view()
		

	output:
		set val(smpl_id), file("${smpl_id}.agg.vcf") into agg_vcf_ch 

	script:

	"""
	aggregate_fusions.pl \\
		--fusioncatcher ${fusionCatcher_file} \\
		--starfusion ${starFusion_file} \\
		--jaffa ${fusionJaffa_file} \\
		--priority fusioncatcher,jaffa,starfusion > ${smpl_id}.agg.vcf
	"""
	}


// import result to coyote
process import_to_coyote {
	publishDir "${params.crondir}/coyote", mode: 'copy'

	input:
		set clarity_id, file(class_report), file(agg_vcf), file(rnaseq_QC), file(salmon_expr), lab_id, pool_id, assay from classification_report.join(agg_vcf_ch.join(final_QC.join(salmon_expr_ch.join(coyote_meta))))
		
	when:
		params.coyote

	output:
		file("${id}.fusion.validation.coyote")
	
	script:
		id= "${lab_id}_validation"
		group= 'fusion_validation_nf'
	
	"""
	echo "import_fusion_to_coyote.pl --classification $OUTDIR/finalResults/${class_report} --fusions $OUTDIR/finalResults/${agg_vcf} --id ${id} --qc $OUTDIR/finalResults/${rnaseq_QC} --group ${group} --expr $OUTDIR/finalResults/${salmon_expr} --clarity-sample-id ${clarity_id} --clarity-pool-id ${pool_id}" > ${id}.fusion.validation.coyote

	"""
	}


