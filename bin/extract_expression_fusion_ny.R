#!/usr/bin/env Rscript


suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("jsonlite"))
suppressPackageStartupMessages(library("EnsDb.Hsapiens.v86"))
suppressPackageStartupMessages(library("reshape"))
library(biomaRt)

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
   
  cat("Possible arguments:
    create-reference      - to create reference
    salmon.file output.file - extract expression from salmon file\n")
  
  q(save = "no")
  
}

genes_of_intrest <- c("MECOM","PBX1","MN1","CD34","CD33","IL1RAP","IL3RA","CD14","MPO","GYPA","MPL","CD3G","CD3D","TRA","CD19","MS4A1","NCAM1")
edb <- EnsDb.Hsapiens.v86
gd<-transcripts(edb,return.type="DataFrame")



if(args[1]=="create-reference"){
  
  cat("Creating genes_of_interest.tsv file ...")  

  ensembl=useMart("ensembl")  # using ensembl database data
  ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)   # from ensembl using homosapien gene data
  genes.with.id=getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters = "hgnc_symbol",values =genes_of_intrest, mart= ensembl)

  #genes_file = paste(ref_path, "genes_of_interest.tsv", sep="")
  write.csv(genes.with.id[,c("ensembl_gene_id", "hgnc_symbol")],file="genes_of_interest.tsv",row.names = F)

  cat("Creating reference...")
   
  files <- list.files(path="/fs1/results/rnaseq_fusion/quant/", pattern="*quant.sf$",full.names = T,include.dirs = F)  
  files<-files[!file.info(files)$isdir]

  allData <- tximport(files=files,type="salmon",txOut = F,tx2gene = gd[,c("tx_name","gene_id")])
  colnames(allData$abundance)<-basename(files)
  
  abundance<-allData$abundance
  
  
  #ref_file = paste(ref_path, "reference_expression.all.tsv", sep="")
  write.csv(abundance[intersect(rownames(abundance),genes.with.id$ensembl_gene_id),],"reference_expression.all.tsv")
  cat("Reference written to reference_expression.all.tsv \n ")
  q(save = "no")
  
  #plotMDS(log2(cpm(abundance)))
  
}


if(file.exists(args[2])){

  genes.with.id<-read.csv(args[1])
  rownames(genes.with.id)<-genes.with.id$ensembl_gene_id
  sample <- tximport(files=args[2],type="salmon",txOut = F,tx2gene = gd[,c("tx_name","gene_id")])
  #sample <- tximport(files="/data/bnf/premap/rnaseq/7460-19_0.salmon",type="salmon",txOut = F,tx2gene = gd[,c("tx_name","gene_id")])
  reference <- read.csv(args[3] ,row.names = 1)
  #"/fs1/resources/ref/hg38/data/salmon/reference_expression.all.tsv"
  available_ens <- unique(intersect(rownames(sample$abundance),rownames(reference)))
  
  sample.goi <- data.frame(sample_expression=sample$abundance[available_ens,],genes.with.id[available_ens,])
  sample.goi<-sample.goi[complete.cases(sample.goi),]
  reference.goi <- data.frame(references=reference[available_ens,],genes.with.id[available_ens,])
  reference.goi <- reference.goi[complete.cases(reference.goi),]
  
  
  #POPULATION PARAMETER CALCULATIONS
  reference.goi$reference_sd <- apply(reference.goi[,grep("quant",colnames(reference.goi))],1,sd) 
  reference.goi$reference_mean <- apply(reference.goi[,grep("quant",colnames(reference.goi))],1,mean) 
  reference.goi$reference_median <- apply(reference.goi[,grep("quant",colnames(reference.goi))],1,median) 
  
  sample.goi <- merge(sample.goi,reference.goi[,c("hgnc_symbol", "reference_sd","reference_mean","reference_median")],by="hgnc_symbol")
  sample.goi$reference_mean_mod <- ifelse(sample.goi$reference_mean<1,1,sample.goi$reference_mean)
  sample.goi$sample_mod <- ifelse(sample.goi$sample<1,1,sample.goi$sample)
  sample.goi$z <- (sample.goi$sample_mod - sample.goi$reference_mean_mod)/sample.goi$reference_sd
  
  rownames(sample.goi)<-NULL
  rownames(reference.goi)<-NULL
  
  colnames(reference.goi) <- gsub('\\.',"_",colnames(reference.goi))
  colnames(reference.goi) <- gsub('references_X',"",colnames(reference.goi))
  
  
  results_for_json.list<-list(sample=sample.goi,reference=reference.goi,expression_version="1.0")
  
  write(toJSON(results_for_json.list,pretty = T,auto_unbox=T,),file = args[4])
  
  cat(paste0("Results written to ", args[4]))
  q(save = "no")
  
  
}
