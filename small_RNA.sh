:qmkdir 0.data 1.fastqc
ls HMZ0039/*/*L001*|rush -k 'cat {} >0.data/{/%}.fastq' --dry-run
nohup ls HMZ0039/*/*L002*|rush -k 'cat {} >>0.data/{/%}.fastq' --dry-run&
##fastqc QC
nohup ls 0.data/*|rush -k 'fastqc {} -o 1.fastqc' &
##install cutadapt and run it
conda install -c bioconda cutadapt
ls 0.data/*|rush -k 'nohup cutadapt -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" -o 2.cutadapter/{%.}_clean.fq {} --discard-untrimmed -j 30 &'
##seqkit 18-60seq
ls 2.cutadapter/*.fq | rush -k "seqkit fx2tab -l {}|awk '\$5>=18'|awk '\$5<=60'|awk '{print \$3}'|sort|uniq -c|awk '{print \">\"\$2\"_\"\$1\"\n\"\$2}' > 3.len18-60/{%@(\d+)_.*}_len18-60.fa" 
##rm structural RNA
ls 3.len18-60/*.fa|rush -k 'bowtie -p 2 /home/lpeng/Gmatic5/genome/tair10/tair10_structure_RNA --un 4.rm_stru_RNA/{%.@(\d+)_len.*}_rmstruRNA.fa -v 3 -f {} > /dev/null'
##small to genome
ls 4.rm_stru_RNA/*.fa|rush -k 'bowtie ~/Gmatic5/genome/tair10/tair10 -k 20 -m 20 -v 0 --best --strata -f {} > 5.map_to_genome/{%:@(.+)_.*}.bwt'
##bwttobed
ls 5.map_to_genome/*.bwt |rush -k "sed 's/_/\\t/' {}|awk '{print \$4\"\\t\"\$5\"\\t\"\$5+length(\$1)\"\\t\"\$1\"\\t\"\$2\"\\t\"length(\$1)\"\\t\"\$NF+1}' > 6.bed/{%.}.bed"
##create windows 500
bedtools makewindows -g /home/lpeng/Gmatic5/genome/tair10/tair10.chromSize -w 500 -i winnum|awk -v OFS="\t" '{print $1,$2,$3,$1"_"$4 }'>tair10_win500.bed
##bed_with_windows
#create all sRNA in all samples
cd 6.bed
cat *.bed|cut -f1-4,6,7|sort -k1,1 -k2,2n |uniq>all.bed
#map all sRNAs in tair10 windows
bedtools intersect -a all.bed -b ../tair10_win500.bed -f 0.5 -wa -wb >all_win.bed
#create counts table in all samples
csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' all_win.bed 1.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 2.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 3.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 4.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 5.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 6.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 7.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 8.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 9.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 10.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 11.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 12.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 13.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 14.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 15.bed -H -k --fill 0|csvtk join -t -f '1,2,3,4,5,6;1,2,3,4,6,7' - 16.bed -H -k --fill 0 > 7.table/sRNA_counts_win.tsv
cat 7.table/sRNA_counts_win_head.txt 7.table/sRNA_counts_win.tsv > 7.table/sRNA_counts_win_w_header.tsv
###annotate gene.bed
less /home/lpeng/Gmatic5/genome/tair10/tair10.gff |grep 'gene'|cut -f 1,4,5,9 > /home/lpeng/Gmatic5/genome/tair10/gene.bed
cat 7.table/sRNA_HNA_win_w_header_counts_up_5.tsv |sed '1d'|less >7.table/sRNA_HNA_win_w_header_counts_up_5_1.tsv
bedtools intersect -a 7.table/sRNA_HNA_win_w_header_counts_up_5_1.tsv -b /home/lpeng/Gmatic5/genome/tair10/gene.bed -f 1 -wa -wb |less -S|cut -f 1-26,30|less  -S|awk -F ";" '{print $NF"\t"$0}'|sed 's/Name=//'|csvtk join -t -H - /home/lpeng/Gmatic5/genome/tair10/tair10_gene_anno.tsv >7.table/sRNA_HNA_win_w_header_counts_up_5_anno.tsv


###DE_with_annotation
less 8.DEsRNA/sRNA_col0_dcl3_10_DE.tsv |sed '1d'|cut -f 2-|bedtools intersect -a - -b /home/lpeng/Gmatic5/genome/tair10/gene.bed -f 1 -wa -wb|awk -F ";" '{print $NF"\t"$0}'|sed 's/Name=//'|csvtk join -t -H - /home/lpeng/Gmatic5/genome/tair10/tair10_gene_anno.tsv |cut -f 2-> 8.DEsRNA/sRNA_col0_dcl3_10_DE_with_anno.tsv
less 8.DEsRNA/sRNA_col0_dcl234_10_DE.tsv |sed '1d'|cut -f 2-|bedtools intersect -a - -b /home/lpeng/Gmatic5/genome/tair10/gene.bed -f 1 -wa -wb|awk -F ";" '{print $NF"\t"$0}'|sed 's/Name=//'|csvtk join -t -H - /home/lpeng/Gmatic5/genome/tair10/tair10_gene_anno.tsv |cut -f 2-> 8.DEsRNA/sRNA_col0_dcl234_10_DE_with_anno.tsv
less 8.DEsRNA/sRNA_col0_hen1_10_DE.tsv |sed '1d'|cut -f 2-|bedtools intersect -a - -b /home/lpeng/Gmatic5/genome/tair10/gene.bed -f 1 -wa -wb|awk -F ";" '{print $NF"\t"$0}'|sed 's/Name=//'|csvtk join -t -H - /home/lpeng/Gmatic5/genome/tair10/tair10_gene_anno.tsv |cut -f 2-> 8.DEsRNA/sRNA_col0_hen1_10_DE_with_anno.tsv


###win_with_annotation

