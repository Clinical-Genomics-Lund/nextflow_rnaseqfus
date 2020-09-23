
read1=$1 
read2=$2

gunzip < $read1 > read1.fastq
gunzip < $read2 > read2.fastq
    
n_reads=$(cat read1.fastq | grep  '@'|wc -l) 
#n_subsample=$(echo ${n_reads}/4  |  bc  -l  | awk '{print int($1)}')
echo "Number of reads: ${n_reads}" 

if (("$n_reads" > 65000000)); then 
paste read1.fastq read2.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t");} }' | 
awk -v k=65000000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x- 1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | 
awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "read1_sub.fastq" ;print $2"\n"$4"\n"$6"\n"$8 > "read2_sub.fastq" }' 
gzip read1_sub.fastq
gzip read2_sub.fastq
else 
mv read1.fastq > read1_sub.fastq
mv read2.fastq > read2_sub.fastq 
gzip read1_sub.fastq
gzip read2_sub.fastq   
fi
