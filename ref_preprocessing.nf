#!/usr/bin/dev nextflow

Channel .fromPath(params.ref_salmon)
			.ifEmpty { exit 1, " Reference file/directory not found!" }
            .into{transcriptome_ref}

Channel
	.fromPath(params.genome_fasta)
	.into{genome_fasta_ch; genomeRef_salmon}

Channel
	.fromPath(params.genome_gtf)
	.into{gtf_star_index}

/*** create  star genome index ***/
process build_star_index {

	publishDir "${params.ref_dir}/star_refGenome_index", mode:'copy'
	cpus = 8
	when:
    params.star_inedx
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

/**** remove version in transcriptome ref file and then index it *****/
process cleanTranscriptomeVersion {
        publishDir "${params.ref_dir}/transcriptome_ref/", mode : 'copy'
		errorStrategy 'ignore'
        when:
            params.salmon_index  
        input :
        file (tr_fasta) from transcriptome_ref

        output:
        file 'transcripts.cleanversion.fa' into transcriptomeRef_cleanversion_salmon
        script:
        """
        Rscript /opt/fuseq/FuSeq_v1.1.2_linux_x86-64/R/excludeTxVersion.R ${tr_fasta} transcripts.cleanversion.fa
        """
}

process create_refIndex {
	errorStrategy 'ignore'
	cpus = 12
	publishDir "${params.ref_dir}/transcriptome_ref/", mode:'copy'

	when:
		params.salmon_index  

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
	salmon index --threads ${task.cpus} -t gentrome.fa.gz -d decoys.txt  -i salmon_index --gencode
	"""  
	}


/**** Download and create genome config file for fastqScreen **********/
process fastqscreen_getGenome{ 
		
		publishDir "${params.ref_dir}", mode: 'copy'
		when :
        params.fastqscreen_getgenome
		output :
		file "*" into output_ch
		file "fastq_screen.conf" into fastq_screen_config_ch
		script:
		"""
		fastq_screen --get_genomes
		"""
	}
