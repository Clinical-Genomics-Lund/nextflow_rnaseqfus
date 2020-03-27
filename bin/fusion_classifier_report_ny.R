#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#/usr/bin/Rscript /data/bnf/scripts/fusion_classifier_report.R '8948-19-fusions' '/data/bnf/premap/rnaseq/8948-19-fusions_0.salmon' '/data/bnf/postmap/rnaseq/8948-19-fusions.STAR.fusionreport'
#Paths of input files have been changed by Sima. (line 23, 24)

library("reshape")
library("knitr")
library("jsonlite")



options(digits=2)

options(stringsAsFactors = F)


sample_name <- args[1]
quant_file <- args[2]
out_json <- args[5]


load(args[3]) #"/fs1/resources/ref/hg38/data/fusion_classifier/hem_classifier.salmon.20190510.RData"
load(args[4]) #"/fs1/resources/ref/hg38/data/fusion_classifier/ensembl_annotation.RData"


quant <- read.table(quant_file, header = T, row.names = 1)
quant$ensembl_transcript_id <- rownames(quant)
quant <- merge(quant,result,by="ensembl_transcript_id")
sample <- quant$TPM
symbol <- quant$hgnc_symbol


#symbol <- gsub("-","_",symbol)
names(sample)<-symbol
sample <- as.matrix(sample)

# function to remove the rules from SwitchBox classifier when the genes that involved in the rules are not in the sample
keep_rules <- function(SampleData, classifier){
  # each rule has two genes, both of them should be in the sample data
  # otherwise this rule will be remove to avoid any errors
  sum  <- classifier$TSPs[,1] %in% rownames(SampleData) + classifier$TSPs[,2] %in% rownames(SampleData)
  
  # subset the classifier's rules
  classifier$name  <- paste0("TSPs", sum(sum==2))
  classifier$score <- classifier$score[sum==2]
  classifier$TSPs  <- classifier$TSPs[sum==2,]
  return(classifier)
}



# classify new samples
sample_score <- function(data, classifier){
  
  classifier <-  keep_rules(data,classifier)
  trues <- data[classifier$TSPs[,1],]  > data[classifier$TSPs[,2],]
 # print(paste0("# of true rules:", sum(trues)))
  score <- sum(classifier$score[trues])/sum(classifier$score)
  return(c(score,sum(trues),length(classifier$score)))
}


results <- data.frame(matrix(NA, nrow = length(seq(1)), 
                             ncol = length(seq(1:4))))
colnames(results)<-c("classifier", "score", "true_genes","total_genes")

# scores for main MI subtypes


results<-rbind(results,c("BA+BA-like",sample_score(sample,ba.classifier)))
results<-rbind(results,c("HeH",sample_score(sample,heh.classifier)))
results<-rbind(results,c("DUX4-high",sample_score(sample,dux4.classifier)))
results<-rbind(results,c("ETV6-RUNX1",sample_score(sample,ER.classifier)))
results<-rbind(results,c("MLL",sample_score(sample,MLL.classifier)))
results<-rbind(results,c("Nonspecific",sample_score(sample,normal.classifier)))
results<-rbind(results,c("Other",sample_score(sample,other.classifier)))
results<-rbind(results,c("TCF3-PBX1",sample_score(sample,tcf3.classifier)))



results <- results[-1,]
colnames(results)<-c("Klass","Poäng", "Sanna", "Totalt")
results$Poäng<-as.numeric(results$Poäng)
results<-results[order(results$Poäng,decreasing = T),]

results$Poäng<-as.numeric(results$Poäng)
results$Sanna<-as.numeric(results$Sanna)
results$Totalt<-as.numeric(results$Totalt)

results_for_json <- results
colnames(results_for_json)<-c("class","score","true","total")

results_for_json.list<-list(classifier_results=results_for_json,classifier_version="HEM_1.01")

write(toJSON(results_for_json.list,pretty = T,auto_unbox=T),file = out_json)
