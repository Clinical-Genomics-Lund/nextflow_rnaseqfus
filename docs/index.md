---
title: 'GMS-Myeloid Pipeline'
output: 
  rmdformats::html_clean:
    toc_depth: 2
    toc_float: true
    theme: united
    highlight: tango
    fig_width: 7
    fig_height: 6
    fig_caption: true
    keep_md: true
    self_contained: true
    thumbnails: true
    lightbox: true
    gallery: false
    url: blue

date: '2023-02-06'
---



## Pipeline Details

**Name:** GMS-Myeloid Pipeline

**Version:** v1.0

**Contact Person:** Ram

**Location:** /home/ramsainanduri/Pipelines/Validation_Documents/mkdocs_test

## Description

This is Mk Docs pipeline. Test follows...

## Minimum Requirements

**Nextflow version:** 19

## List of Software Used

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> ID </th>
   <th style="text-align:left;"> Tool/Software </th>
   <th style="text-align:left;"> Version </th>
   <th style="text-align:left;"> Build </th>
   <th style="text-align:left;"> Type </th>
   <th style="text-align:left;"> Availability </th>
   <th style="text-align:left;"> Container Name </th>
   <th style="text-align:left;"> External Contact Person If Any </th>
   <th style="text-align:left;"> Reference </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> <a href="https://github.com/lh3/bwa" target="_blank">BWA</a> </td>
   <td style="text-align:left;"> v0.7.17-r1188 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/lh3/bwa/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.48550/arXiv.1303.3997" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> <a href="https://cnvkit.readthedocs.io/en/stable/" target="_blank">CNVkit</a> </td>
   <td style="text-align:left;"> v0.9.5 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/etal/cnvkit/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1371/journal.pcbi.1004873" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> <a href="https://cnvkit.readthedocs.io/en/stable/" target="_blank">CNVkit</a> </td>
   <td style="text-align:left;"> v0.9.9 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/etal/cnvkit/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1371/journal.pcbi.1004873" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> <a href="https://github.com/dellytools/delly" target="_blank">delly</a> </td>
   <td style="text-align:left;"> v0.8.7 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/dellytools/delly/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1093/bioinformatics/bts378" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:left;"> <a href="https://github.com/freebayes/freebayes" target="_blank">freebayes</a> </td>
   <td style="text-align:left;"> v1.3.5 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/freebayes/freebayes/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.48550/arXiv.1207.390" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:left;"> <a href="https://gatk.broadinstitute.org/hc/en-us" target="_blank">GATK</a> </td>
   <td style="text-align:left;"> v4.1.9.0 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> gatk_4.1.9.0.sif </td>
   <td style="text-align:left;"> https://github.com/broadinstitute/gatk/issues </td>
   <td style="text-align:left;"> <a href="http://dx.doi.org/10.1101/gr.107524.110" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> <a href="https://github.com/Illumina/manta" target="_blank">Manta</a> </td>
   <td style="text-align:left;"> v1.6.0 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/Illumina/manta/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1093/bioinformatics/btv710" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:left;"> <a href="https://melt.igs.umaryland.edu/index.php" target="_blank">MELT</a> </td>
   <td style="text-align:left;"> v2.1.5 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1101%2Fgr.218032.116" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:left;"> <a href="https://github.com/genome/pindel" target="_blank">pindel</a> </td>
   <td style="text-align:left;"> v0.2.5b9, 20160729 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/genome/pindel/issues (kaiye@xjtu.edu.cn) </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1093%2Fbioinformatics%2Fbtp394" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:left;"> <a href="https://www.r-project.org/" target="_blank">R</a> </td>
   <td style="text-align:left;"> v4.1.0 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://www.r-project.org/ </td>
   <td style="text-align:left;"> <a href="-" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:left;"> <a href="https://github.com/biod/sambamba" target="_blank">Sambamba</a> </td>
   <td style="text-align:left;"> v0.8.0 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/biod/sambamba/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1093/bioinformatics/btv098" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:left;"> <a href="https://github.com/samtools/samtools" target="_blank">Samtools</a> </td>
   <td style="text-align:left;"> v1.12 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/samtools/samtools/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1093/bioinformatics/btp352" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:left;"> <a href="https://support.sentieon.com/manual/" target="_blank">Sentieon</a> </td>
   <td style="text-align:left;"> v202010.01 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Commercial </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> - </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1101/115717" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:left;"> <a href="https://github.com/J35P312/SVDB" target="_blank">svdb</a> </td>
   <td style="text-align:left;"> v2.2.0 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/J35P312/SVDB/issues </td>
   <td style="text-align:left;"> <a href="-" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:left;"> <a href="https://github.com/AstraZeneca-NGS/VarDict" target="_blank">Vardict</a> </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/AstraZeneca-NGS/VarDict/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1093/nar/gkw227" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:left;"> <a href="https://github.com/vcftools/vcftools" target="_blank">VCFtools</a> </td>
   <td style="text-align:left;"> v0.1.16 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> https://github.com/vcftools/vcftools/issues </td>
   <td style="text-align:left;"> <a href="http://dx.doi.org/10.1093/bioinformatics/btr330" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:left;"> <a href="https://github.com/Ensembl/ensembl-vep" target="_blank">VEP</a> </td>
   <td style="text-align:left;"> v95.3 </td>
   <td style="text-align:left;"> HG19/GRCh37 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> container_VEP.sif </td>
   <td style="text-align:left;"> https://github.com/Ensembl/ensembl-vep/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1186/s13059-016-0974-4" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:left;"> <a href="https://github.com/Ensembl/ensembl-vep" target="_blank">VEP</a> </td>
   <td style="text-align:left;"> v103.1 </td>
   <td style="text-align:left;"> HG38/GrCh38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> ensembl-vep_release_103.sif </td>
   <td style="text-align:left;"> https://github.com/Ensembl/ensembl-vep/issues </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1186/s13059-016-0974-4" target="_blank">open</a> </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:left;"> <a href="https://genome.sph.umich.edu/wiki/Vt" target="_blank">Vt</a> </td>
   <td style="text-align:left;"> v0.5 </td>
   <td style="text-align:left;"> HG19/HG38 </td>
   <td style="text-align:left;"> Open-Source </td>
   <td style="text-align:left;"> Container/Image </td>
   <td style="text-align:left;"> SomaticPanelPipeline_2021-06-24.sif </td>
   <td style="text-align:left;"> Adiran (atks@umich.edu) </td>
   <td style="text-align:left;"> <a href="https://doi.org/10.1093/bioinformatics/btv112" target="_blank">open</a> </td>
  </tr>
</tbody>
</table>
