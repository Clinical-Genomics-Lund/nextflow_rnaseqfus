#!/usr/bin/env Rscript
suppressMessages(suppressWarnings(library(jsonlite)))
suppressMessages(suppressWarnings(library(optparse)))

option_list <- list( 
  
  make_option(c("-s", "--star"), type="character", default="", 
              help="STAR final log fil [/data/NextSeq2/180711_NB501699_0056_AH3MYMAFXY/tmp/tmp/rnaseq/4177/Log.final.out]",
              metavar="character"),
  make_option(c("-i", "--id"), type="character", default="", 
              help="Sample name [4014-14]",
              metavar="character"),
  make_option(c("-f", "--fusion"), type="character", default="", 
              help="Fusion catcher results [/data/bnf/proj/rnaseq/cmd-pipe-out/tmp/rnaseq/3376_02/fusioncatcher/final-list_candidate-fusion-genes.hg19.txt]",
              metavar="character"),
  make_option(c("-p", "--provider"), type="character", default="", 
              help="Genotypes from provider [/data/bnf/proj/rnaseq/cmdpipe-fix-test/TESTTTT/postmap/rnaseq/HEJHEJHEJ.STAR.provider.genotypes]",
              metavar="character"),
  make_option(c("-g", "--genebody"), type="character", default="", 
              help="genebody coverage [/data/bnf/postmap/rnaseq/19KFU0009.STAR.genebodycov.geneBodyCoverage.curves] ",
              metavar="character"),
  make_option(c("-l", "--flendist"), type="character", default="", 
              help="Frag lengths from Salmon [/data/bnf/premap/rnaseq/7728_14_0.salomon.salomon/libParams/flenDist.txt] ",
              metavar="character")
)


# load("/data/bnf/proj/rnaseq/example_opt.Rdata")

opt <- parse_args(OptionParser(option_list=option_list))

options(stringsAsFactors = F)

if (length(args)==0) {
  stop("At least one argument must be supplied (STAR log file).\n", call.=FALSE)
}

#
# STAR
#

star_qc <- opt$star
fusion_list <- opt$fusion
sample_id <- opt$id
genebody_cov <- opt$genebody
provider <- opt$provider
flendist <- opt$flendist


if (!file.exists(star_qc)){
  cat(paste0("File not found: ",star_qc))
  q(save = "no")
}

qc.log<-read.csv(star_qc,sep="\t",header=F)

qc.data<-data.frame(uppmatt=c(qc.log[grep("Uniquely mapped reads number",qc.log$V1),"V2"],
                              qc.log[grep("Uniquely mapped reads %",qc.log$V1),"V2"],
                              qc.log[grep("% of reads mapped to multiple loci",qc.log$V1),"V2"],
                              qc.log[grep("Mismatch rate per base, %",qc.log$V1),"V2"],
                              qc.log[grep("Number of splices: Annotated",qc.log$V1),"V2"],
                              qc.log[grep("Number of splices: Non-canonical",qc.log$V1),"V2"])
)

qc.data<-rbind(qc.data,c(as.integer(as.numeric(qc.data$uppmatt[5])/as.numeric(qc.data$uppmatt[6])),">100"))

qc.data$uppmatt<-as.numeric(gsub("%","",qc.data$uppmatt))
#qc.data<-gsub("%", "", qc.data$uppmatt)
qc.data.l <- as.list(qc.data$uppmatt)

names(qc.data.l) =  c("tot_reads","mapped_pct", "multimap_pct",
                      "mismatch_pct", "canon_splice", "non_canon_splice",
                      "splice_ratio")

#
# genebodycov
#


if (file.exists(genebody_cov)){
  qc.gbc<-read.csv(genebody_cov,sep="\t",header=T,row.names = NULL)
  
  qc.gbc.res <- unname(as.vector(qc.gbc[,-1]))
  #names(qc.gbc.res)<-""
  qc.data.l[["genebody_cov"]]<-unlist(qc.gbc.res)
  df <- data.frame(perc=1:100,val=unlist(qc.gbc[,-1]))
  df$val <- 1000*df$val/max(df$val)
  qc.data.l[["genebody_cov_slope"]]<-summary(lm(val~perc,df[20:80,]))$coefficients[2,1]
  
}

# 
# 
# par(mfrow=c(3,2))
# for(i in c("18KFU0016","19KFU0013","18KFU0006","19KFU0009","19KFU0012","19KFU0002" )){
# 
#   qc.gbc<-read.csv(paste0("/data/bnf/postmap/rnaseq/",i,".STAR.genebodycov.geneBodyCoverage.txt"),sep="\t",header=T,row.names = NULL)
#   df <- data.frame(perc=1:100,val=unlist(qc.gbc[,-1]))
#   
#   plot(df$val,main=paste(i, summary(lm(val~perc,df[20:80,]))$coefficients[2,1]))
#   abline(lm(val~perc,df[20:80,]),col="red")
#   
# }

#
# Fusioncatcher
#

if (file.exists(fusion_list)){
  nfusions <- system(paste0("wc -l <" ,fusion_list),intern=T)
  #  q(save = "no")
  qc.data.l[["tot_fusions"]]<-as.numeric(nfusions)
}


#
# provider
#

if (file.exists(provider)){
  provider_res<-read.csv(provider,sep="\t",header=T,row.names = NULL)
  colnames(provider_res)<-c("rsid","genotype")
  provider_res_v <- provider_res$genotype
  names( provider_res_v )<-provider_res$rsid
  qc.data.l[["provider_genotypes"]]<-unbox(data.frame(t(provider_res_v)))
  
  qc.data.l[["provider_called_genotypes"]]<-sum(provider_res[,2]!="",na.rm = T)
}


#
# fraglen
#

if (file.exists(flendist)){
  
 flendist_dat <- read.table(flendist)
  
 qc.data.l[["flendist"]]<-seq_along(flendist_dat)[flendist_dat == max(flendist_dat)]

}

### Print

qc.data.l[["sample_id"]]<-sample_id


toJSON(qc.data.l,pretty=T,digits=4,auto_unbox = T)
