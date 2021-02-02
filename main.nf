#!/usr/bin/dev nextflow

OUTDIR = params.outdir+'/'+params.subdir


csv = file(params.csv)
// Print commit-version of active deployment
file(params.git)
    .readLines()
    .each { println "git commit-hash: "+it }
// Print active container
container = file(params.container).toRealPath()
println("container: "+container)

workflow.onComplete {

	def msg = """\
		Pipeline execution summary
		---------------------------
		Completed at: ${workflow.complete}
		Duration    : ${workflow.duration}
		Success     : ${workflow.success}
		scriptFile  : ${workflow.scriptFile}
		workDir     : ${workflow.workDir}
		exit status : ${workflow.exitStatus}
		errorMessage: ${workflow.errorMessage}
		errorReport :
		"""
		.stripIndent()
	def error = """\
		${workflow.errorReport}
		"""
		.stripIndent()

	base = csv.getBaseName()
	logFile = file("/fs1/results/cron/logs/" + base + ".complete")
	logFile.text = msg
	logFile.append(error)
}

Channel
    .fromPath(params.csv)
    .splitCsv(header:true)
    .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }
    .into {reads_starfusion; reads_sub; reads_jaffa; reads_align; reads_salmon; reads_fastqscreen; reads_meta}

Channel
	.fromPath(params.csv)
	.splitCsv(header:true)
	.map{row -> tuple(row.id,row.clarity_sample_id, row.clarity_pool_id)}
	.into{coyote_meta;cdm_meta}




/***********************************/
/*      Part1: Alignment & QC     */
/**********************************/
if (!params.subsampling) {

   reads_fusioncatcher =Channel
   .fromPath(params.csv)
   .splitCsv(header:true)
   .map{ row-> tuple(row.id, file(row.read1), file(row.read2)) }.view()

}else {
      process subsampling_fastqs {
        //Downsample fastqs to 65000000. To change it shall be done in the subsampling.sh script located in bin directory
	memory 75.GB 
	//when:
	//params.subsampling

	input:
		set val(smpl_id) , read1, read2 from reads_sub
	output:
		set val(smpl_id) , file("${smpl_id}_read1_sub.fastq.gz"), file("${smpl_id}_read2_sub.fastq.gz") into reads_fusioncatcher

	script:
	"""
	bash subsampling.sh $read1 $read2
	mv  read1_sub.fastq.gz  ${smpl_id}_read1_sub.fastq.gz
	mv  read2_sub.fastq.gz  ${smpl_id}_read2_sub.fastq.gz
	"""
	}
}

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
		set val(smpl_id), file("${smpl_id}.Aligned.sortedByCoord.out.bam") into qualimap_bam, provider_bam, geneBody_bam, bamidx
		set val(smpl_id), file("${smpl_id}.Log.final.out") into star_logFinalOut_ch
		
		
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
	
process index_bam {

	publishDir "$OUTDIR/bam", mode :'copy'
	 
	input:
		set val(smpl_id), file(bam) from  bamidx
	output:
		file "${smpl_id}.Aligned.sortedByCoord.out.bam.bai"
		
	script:
	"""
	sambamba index --show-progress -t 8 ${bam}
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
	
	cpus 1
	
	when:
		params.qc || params.bodyCov
	
	input :
	
		set val(smpl_id),file(bam_f) from geneBody_bam
		
			
	output:
		
		set val(smpl_id), file("${smpl_id}.geneBodyCoverage.txt") into gene_bodyCov_ch
	
		
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
/* Part2 : fusion identification    */
/* ******************************** */

process star_fusion{
	//errorStrategy 'ignore'
	scratch true
	tag "${smpl_id}"
	cpus 18
	memory  60.GB
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
	tag "${smpl_id}"
	cpus 48 
	memory 180.GB
	scratch true


	publishDir "$OUTDIR/fusion", mode: 'copy'
	
	when: 
		params.fusioncatcher || params.fusion
	
	input:
		set val(smpl_id) , file(read1), file(read2) from reads_fusioncatcher
    	
	output:
		set val(smpl_id), file("${smpl_id}.final-list_candidate-fusion-genes.txt") into final_list_fusionCatcher_ch, final_list_fusionCatcher_agg_ch
		set val(smpl_id), file("${smpl_id}.fusioncatcher.xls") into filter_fusion_ch
	
	script:
	option = params.singleEnd ? read1 : "${read1},${read2}"
    	//def extra_params = params.fusioncatcher_opt ? "${params.fusioncatcher_opt}" : ''
    	"""
   	fusioncatcher.py -d ${params.fusionCatcher_ref} -i ${option}  --threads ${task.cpus} --limitSjdbInsertNsj 50000000 --limitOutSJcollapsed 2000000 -o ./${smpl_id}.fusioncatcher
	filter_aml_fusions.pl ./${smpl_id}.fusioncatcher > ${smpl_id}.fusioncatcher.xls
	mv  ./${smpl_id}.fusioncatcher/final-list_candidate-fusion-genes.txt ${smpl_id}.final-list_candidate-fusion-genes.txt
    	"""
	}



process jaffa{
	tag "${smpl_id}"
	errorStrategy 'ignore'
	publishDir  "$OUTDIR/fusion", mode: 'copy'
	memory 75.GB 
	cpus 4
	
	when:
		params.jaffa 
	
	input:
		set val(smpl_id) , file(read1), file(read2) from  reads_jaffa
	
	output:
		set val(smpl_id) ,file ("${smpl_id}.jaffa_results.csv") into jaffa_csv_ch
		file "*.fasta" into jaffa_fasta_ch 
	
    	script:

   	"""
   	bpipe run -m 75GB -n ${task.cpus} -p genome=hg38 -p refBase="${params.jaffa_base}" ${params.jaffa_file}  ${read1} ${read2}
	mv  jaffa_results.csv ${smpl_id}.jaffa_results.csv
   	"""
	}



/*****************************************/
/*  Part3 : Expression quantification    */
/*****************************************/

process salmon{

	errorStrategy 'ignore'
	tag "${smpl_id}"
	publishDir "$OUTDIR/quant", mode:'copy'
	cpus  8
	memory 24.GB

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
/*          Part 4: Post processing                */
/***************************************************/
 
process postaln_qc_rna {
	publishDir "$OUTDIR/finalResults" , mode:'copy'
	

	when:
		params.combine 

	input:
		set val(smpl_id), file(star_final), file(fusion), file(geneCov), file(provIder), file(flendist) from star_logFinalOut_ch.join(final_list_fusionCatcher_ch.join(gene_bodyCov_ch.join(provider_output_ch.join(flendist_ch))))
	
	
	output:
		set val(smpl_id), file("${smpl_id}.STAR.rnaseq_QC") into final_QC,finalqc_cmd 

	
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
	
	publishDir "$OUTDIR/finalResults" , mode: 'copy'

	when :
	     params.combine

	input:
		set val(smpl_id), file(fusionCatcher_file), file(starFusion_file) from final_list_fusionCatcher_agg_ch.join(star_fusion_agg_ch)
		//set val(smpl_id), file(fusionCatcher_file), file(starFusion_file), file(fusionJaffa_file) from final_list_fusionCatcher_agg_ch.join(star_fusion_agg_ch.join(jaffa_csv_ch)).view()
		

	output:
		set val(smpl_id), file("${smpl_id}.agg.vcf") into agg_vcf_ch 

	script:

	"""
	aggregate_fusions.pl \\
		--fusioncatcher ${fusionCatcher_file} \\
		--starfusion ${starFusion_file} \\
		--priority fusioncatcher,starfusion > ${smpl_id}.agg.vcf
	"""
	}


// import files to coyote
process import_to_coyote {
	publishDir "${params.crondir}/coyote", mode: 'copy'

	input:
		set smpl_id, file(class_report), file(agg_vcf), file(rnaseq_QC), file(salmon_expr), clarity_id, pool_id from classification_report.join(agg_vcf_ch.join(final_QC.join(salmon_expr_ch.join(coyote_meta))))
		
	when:
		params.coyote

	output:
		file("${id}.coyote")
	
	script:
		id= "${smpl_id}-fusions"
		group= 'fusion'
	
	"""
	echo "import_fusion_to_coyote.pl --classification $OUTDIR/finalResults/${class_report} --fusions $OUTDIR/finalResults/${agg_vcf} --id ${id} --qc $OUTDIR/finalResults/${rnaseq_QC} --group ${group} --expr $OUTDIR/finalResults/${salmon_expr} --clarity-sample-id ${clarity_id} --clarity-pool-id ${pool_id}" > ${id}.coyote

	"""
	}


process  register_to_cdm{
	 publishDir "${params.crondir}/qc", mode: 'copy', overwrite: true
	 cpus 1
	 memory '8 GB'
	 time '1h'
	 input:
		set smpl_id, file(postalignqc), fastq_r1, fatq_r2 ,clarity_id, pool_id from finalqc_cmd.join(reads_meta.join(cdm_meta)) 
	 output:
		set val(smpl_id), file("${smpl_id}.cdm")
	 script:
	 
	 parts = fastq_r1.toString().split('/')
	 parts.println()
	 idx= parts.findIndexOf {it ==~ /......_......_...._........../}
	 rundir= parts[0..idx].join("/")
	 
	 
	"""
	echo "--run-folder ${rundir} --sample-id ${smpl_id} --assay rnaseq-fusion --qc $OUTDIR/finalResults/${postalignqc}" > ${smpl_id}.cdm 
	"""

}

