#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  
  
  cat("Possible arguments:
    create-reference      - to create reference
    salmon.file output.file - extract expression from salmon file\n")
  
  q(save = "no")
  
}


genes_of_intrest <- c("ABL1","ABL2","ACTN4","AFDN","AFF1","AFF3","AFF4","AGGF1","ARHGEF12","ATF7IP","BCL9L","BCR","BTBD18","CASC5","CBFB","CBL","CD74","CENPC","CEP164","CEP85L","CREBBP","CSF1R","CT45A3","DCP1A","DEK","DOCK2","DUX4","EBF1","EEFSEC","ELL","EML1","EPS15","ETV6","FGFR1","FIP1L1","FOXO3","FOXO4","FOXP1","FRYL","GAS7","HLF","IGH","IL3","INPP5D","JAK2","KMT2A","LSM14A","MAML2","MAPRE1","MATR3","MBNL1","MECOM","MEF2D","MKL1","MLLT1","MLLT10","MLLT11","MLLT3","MLLT6","MYH11","NUP153","NUP214","PAX","PBX1","PCM1","PDGFRA","PDGFRB","PICALM","PML","PRPF19","PRRC1","PRRC2B","RANBP2","RARA","RBM15","RCSD1","RUNDC3B","RUNX1","RUNX1T1","SACM1L","SATB1","SEPT11","SEPT5","SFPQ","SLC38A2","SLC9A3R1","SNX2","SSBP2","TBL1XR1","TCF3","TET1","TNIP1","TNRC18","USP2","ZBTB16","ZC3HAV1","ZMIZ1","ZMYND8","ZNF384")


suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("jsonlite"))
suppressPackageStartupMessages(library("EnsDb.Hsapiens.v86"))
suppressPackageStartupMessages(library("reshape"))
# library(biomaRt)
# 
edb <- EnsDb.Hsapiens.v86
gd<-transcripts(edb,return.type="DataFrame")
# 
# 
# ensembl=useMart("ensembl")  # using ensembl database data
# ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)   # from ensembl using homosapien gene data
# 
# genes.with.id=getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),filters = "hgnc_symbol",values =genes_of_intrest, mart= ensembl)
# 
# write.csv(genes.with.id[,c("ensembl_gene_id", "hgnc_symbol")],file="/data/bnf/ref/fuseq/genes_of_interest.tsv",row.names = F)

genes.with.id<-read.csv("/data/bnf/ref/fuseq/genes_of_interest.tsv")
rownames(genes.with.id)<-genes.with.id$ensembl_gene_id



if(args[1]=="create-reference"){
  
  cat("Creating reference...")
  
  
  files <- list.files(path="/data/bnf/proj/rnaseq/expression_levels/cmd_pipe_out/premap/rnaseq/", pattern="*salmon$",full.names = T,include.dirs = F)
  files<-files[!file.info(files)$isdir]
  
  allData <- tximport(files=files,type="salmon",txOut = F,tx2gene = gd[,c("tx_name","gene_id")])
  
  colnames(allData$abundance)<-basename(files)
  
  abundance<-allData$abundance
  
  write.csv(abundance[intersect(rownames(abundance),genes.with.id$ensembl_gene_id),],"/data/bnf/ref/salmon/reference_expression.all.tsv")
  
  cat("Reference written to /data/bnf/ref/salmon/reference_expression.all.tsv \n ")
  q(save = "no")
  
  
  plotMDS(log2(cpm(abundance)))
  
}

if(file.exists(args[1])){
  
  sample <- tximport(files=args[1],type="salmon",txOut = F,tx2gene = gd[,c("tx_name","gene_id")])
  #sample <- tximport(files="/data/bnf/premap/rnaseq/7460-19_0.salmon",type="salmon",txOut = F,tx2gene = gd[,c("tx_name","gene_id")])
  reference <- read.csv("/data/bnf/ref/salmon/reference_expression.all.tsv",row.names = 1)
  
  available_ens <- unique(intersect(rownames(sample$abundance),rownames(reference)))
  
  sample.goi <- data.frame(sample_expression=sample$abundance[available_ens,],genes.with.id[available_ens,])
  sample.goi<-sample.goi[complete.cases(sample.goi),]
  reference.goi <- data.frame(references=reference[available_ens,],genes.with.id[available_ens,])
  reference.goi <- reference.goi[complete.cases(reference.goi),]
  
  
  #POPULATION PARAMETER CALCULATIONS
  reference.goi$reference_sd <- apply(reference.goi[,grep("salmo",colnames(reference.goi))],1,sd) 
  reference.goi$reference_mean <- apply(reference.goi[,grep("salmo",colnames(reference.goi))],1,mean) 
  reference.goi$reference_median <- apply(reference.goi[,grep("salmo",colnames(reference.goi))],1,median) 
  
  sample.goi <- merge(sample.goi,reference.goi[,c("hgnc_symbol", "reference_sd","reference_mean","reference_median")],by="hgnc_symbol")
  sample.goi$reference_mean_mod <- ifelse(sample.goi$reference_mean<1,1,sample.goi$reference_mean)
  sample.goi$sample_mod <- ifelse(sample.goi$sample<1,1,sample.goi$sample)
  sample.goi$z <- (sample.goi$sample_mod - sample.goi$reference_mean_mod)/sample.goi$reference_sd
  
  rownames(sample.goi)<-NULL
  rownames(reference.goi)<-NULL
  
  colnames(reference.goi) <- gsub('\\.',"_",colnames(reference.goi))
  colnames(reference.goi) <- gsub('references_X',"",colnames(reference.goi))
  
  
  results_for_json.list<-list(sample=sample.goi,reference=reference.goi,expression_version="1.0")
  
  write(toJSON(results_for_json.list,pretty = T,auto_unbox=T,),file = args[2])
  
  cat(paste0("Results written to ", args[2]))
  
  q(save = "no")
  
  
}
