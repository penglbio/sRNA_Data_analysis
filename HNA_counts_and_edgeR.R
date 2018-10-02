data <- read.table("/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_counts_win_w_header.tsv",header = T)
library(dplyr)
library(pheatmap)
data1<-data[,11:26]
head(data1)
sRNA <- apply(data1,2,function(x){ifelse(x>=5,1,0)})
head(sRNA_5)
sRNA_5 <- data[apply(sRNA,1,sum)>=1,]
write.table(sRNA_5,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_counts_win_w_header_counts_up_5.tsv",sep="\t",quote = FALSE)
head(sRNA_5)
sRNA_5_sum_sample <- apply(sRNA_5[,11:26],2,sum)
sRNA_5_sum_sample
sRNA_5_RPTM <- apply(sRNA_5[,11:26],2,function(x){x/sRNA_5_sum_sample*10000000})
colnames(sRNA_5_RPTM)
colnames(sRNA_5_RPTM)<-paste0("RPTM_",colnames(sRNA_5[,11:26]))
sRNA_counts_RPTM <- cbind(sRNA_5,sRNA_5_RPTM)
head(sRNA_counts_RPTM)
write.table(sRNA_counts_RPTM,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_counts_RPTM_win_w_header_counts_up_5.tsv",sep="\t",quote = FALSE)
sRNA_5_HNA <- apply(sRNA_5_RPTM,2,function(x){x/sRNA_5[,6]})
sum_w_m<-apply(sRNA_5_HNA,1,sum)
colnames(sRNA_5_HNA)<-paste0("HNA_",colnames(sRNA_5[,11:26]))
sRNA_counts_RPTM_HNA <- cbind(sRNA_counts_RPTM,sRNA_5_HNA)
RPTM <- cbind(sRNA_counts_RPTM[,1:10],sRNA_5_RPTM)
HNA <- cbind(sRNA_counts_RPTM[,1:10],sRNA_5_HNA)
write.table(RPTM,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_RPTM_win_w_header_counts_up_5.tsv",sep="\t")
write.table(HNA,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_HNA_win_w_header_counts_up_5.tsv",sep="\t",row.names = FALSE,quote = FALSE)


###windowns plot figure
HNA <- read.table("/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_HNA_win_w_header_counts_up_5.tsv",header = T)
win_HNA_sum <- group_by(HNA,win_id) %>% summarise(test1=sum(HNA_s1),test2=sum(HNA_s2),test3=sum(HNA_s3),test4=sum(HNA_s4),test5=sum(HNA_s5),
                                                                              test6=sum(HNA_s6),test7=sum(HNA_s7),test8=sum(HNA_s8),test9=sum(HNA_s9),test10=sum(HNA_s10),
                                                                              test11=sum(HNA_s11),test12=sum(HNA_s12),test13=sum(HNA_s13),test14=sum(HNA_s14),
                                                                              test15=sum(HNA_s15),test16=sum(HNA_s16)) 
win_HNA_sum <- as.data.frame(win_HNA_sum)
rownames(win_HNA_sum) <- win_HNA_sum[,1]
win_HNA_sum <- win_HNA_sum[2:17]
colnames(win_HNA_sum) <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12","s13","s14","s15","16")
write.table(win_HNA_sum,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_window_HNA_value_counts_up_5.tsv",sep="\t",row.names =T,quote = FALSE)
widetype_win_HNA_sum_mean <- apply(win_HNA_sum[,1:4],1,mean)
mutant1_win_HNA_sum_mean <- apply(win_HNA_sum[,5:8],1,mean)
mutant2_win_HNA_sum_mean <- apply(win_HNA_sum[,9:12],1,mean)
mutant3_win_HNA_sum_mean <- apply(win_HNA_sum[,13:16],1,mean)
win_HNA_sample_mean <- cbind(widetype_win_HNA_sum_mean,mutant1_win_HNA_sum_mean,mutant2_win_HNA_sum_mean,mutant3_win_HNA_sum_mean)
dim(win_HNA_sample_mean) 
head(win_HNA_sample_mean)
win_HNA_sample_mean_w_m1<- win_HNA_sample_mean[,1:2]
win_HNA_sample_mean_w_m2<- win_HNA_sample_mean[,c(1,3)]
win_HNA_sample_mean_w_m3<- win_HNA_sample_mean[,c(1,4)]
win_HNA_sample_mean_w_m1 <- win_HNA_sample_mean_w_m1[rowSums(win_HNA_sample_mean_w_m1)>=400,]
win_HNA_sample_mean_w_m2 <- win_HNA_sample_mean_w_m2[rowSums(win_HNA_sample_mean_w_m2)>=400,]
win_HNA_sample_mean_w_m3 <- win_HNA_sample_mean_w_m3[rowSums(win_HNA_sample_mean_w_m3)>=400,]

FC_mutant1_vs_widetype <- win_HNA_sample_mean_w_m1[,2]/win_HNA_sample_mean_w_m1[,1]
FC_mutant2_vs_widetype <- win_HNA_sample_mean_w_m2[,2]/win_HNA_sample_mean_w_m2[,1]
FC_mutant3_vs_widetype <- win_HNA_sample_mean_w_m3[,2]/win_HNA_sample_mean_w_m3[,1]


mutant1_FC2_winid <-c(FC_mutant1_vs_widetype[FC_mutant1_vs_widetype >=2],FC_mutant1_vs_widetype[FC_mutant1_vs_widetype<=0.5])
mutant2_FC2_winid <-c(FC_mutant2_vs_widetype[FC_mutant2_vs_widetype>=2],FC_mutant2_vs_widetype[FC_mutant2_vs_widetype<=0.5])
mutant3_FC2_winid <- c(FC_mutant3_vs_widetype[FC_mutant3_vs_widetype>=2],FC_mutant3_vs_widetype[FC_mutant3_vs_widetype<=0.5])

##mutant1 dcl3
head(HNA)
win_mutant1_table <- HNA[HNA$win_id %in% names(mutant1_FC2_winid),][c(5,10,11:18)]
##ratio value log2((mutant+1)/(widetype+1))
win_mutant1_table <- cbind(win_mutant1_table[,1:2], m1_mean=rowMeans(win_mutant1_table[,c(7:10)]), w_mean=rowMeans(win_mutant1_table[,c(3:6)]))
win_mutant1_table
win_mutant1_fig_table <- group_by(win_mutant1_table,win_id,length) %>% summarise(m1_means=sum(m1_mean),w_means=sum(w_mean))
ratio <- log2((win_mutant1_fig_table$m1_means+1)/(win_mutant1_fig_table$w_means+1))
win_mutant1_fig_table<- cbind(win_mutant1_fig_table[,c(1:2)],ratio=ratio)



##mutant2-dcl234
win_mutant2_table <- HNA[HNA$win_id %in% names(mutant2_FC2_winid),][c(5,10,11:15,19:22)]

win_mutant2_table <- cbind(win_mutant2_table[,1:2], m1_mean=rowMeans(win_mutant2_table[,c(7:10)]), w_mean=rowMeans(win_mutant2_table[,c(3:6)]))
win_mutant2_fig_table <- group_by(win_mutant2_table,win_id,length) %>% summarise(m1_means=sum(m1_mean),w_means=sum(w_mean))
ratio <- log2((win_mutant2_fig_table$m1_means+1)/(win_mutant2_fig_table$w_means+1))
win_mutant2_fig_table<- cbind(win_mutant2_fig_table[,c(1:2)],ratio=ratio)

#mutant3-hen1
win_mutant3_table <- HNA[HNA$win_id %in% names(mutant3_FC2_winid),][c(5,10,11:15,23:26)]
win_mutant3_table <- cbind(win_mutant3_table[,1:2], m1_mean=rowMeans(win_mutant3_table[,c(7:10)]), w_mean=rowMeans(win_mutant3_table[,c(3:6)]))
win_mutant3_fig_table <- group_by(win_mutant3_table,win_id,length) %>% summarise(m1_means=sum(m1_mean),w_means=sum(w_mean))
ratio <- log2((win_mutant3_fig_table$m1_means+1)/(win_mutant3_fig_table$w_means+1))
win_mutant3_fig_table<- cbind(win_mutant3_fig_table[,c(1:2)],ratio=ratio)
dim(win_mutant3_fig_table)


write.table(win_mutant1_fig_table,file="/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant1_table_fig.tsv",quote =F,sep="\t",row.names = F,col.names = F)
write.table(win_mutant2_fig_table,file="/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant2_table_fig.tsv",quote =F,sep="\t",row.names = F,col.names = F)
write.table(win_mutant3_fig_table,file="/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant3_table_fig.tsv",quote =F,sep="\t",row.names = F,col.names = F)

system("/var/data/06.Richa_sRNA/smallRNA/dataframe_to_matrix.py /var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant1_table_fig.tsv >/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant1_table_fig_1.tsv")
system("/var/data/06.Richa_sRNA/smallRNA/dataframe_to_matrix.py /var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant2_table_fig.tsv >/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant2_table_fig_1.tsv")
system("/var/data/06.Richa_sRNA/smallRNA/dataframe_to_matrix.py /var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant3_table_fig.tsv >/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant3_table_fig_1.tsv")


fig_m1 <- read.table("/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant1_table_fig_1.tsv")
fig_m1_matrix <- matrix(fig_m1$V3, ncol=43, byrow = T)
colnames(fig_m1_matrix)=c(18:60)
rownames(fig_m1_matrix)=unique(fig_m1$V1)
pdf("9.fig/col0_dcl3_dcl234_hen1.pdf")
pheatmap(fig_m1_matrix,cluster_cols = F,show_rownames=F,main="sRNA ratio")


fig_m2 <- read.table("/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant2_table_fig_1.tsv")
fig_m2_matrix <- matrix(fig_m2$V3, ncol=43, byrow = T)
colnames(fig_m2_matrix)=c(18:60)
rownames(fig_m2_matrix)=unique(fig_m2$V1)
pheatmap(fig_m2_matrix,cluster_cols = F,show_rownames=F,main="sRNA ratio")

fig_m3 <- read.table("/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant3_table_fig_1.tsv")
fig_m3_matrix <- matrix(fig_m3$V3, ncol=43, byrow = T)
colnames(fig_m3_matrix)=c(18:60)
rownames(fig_m3_matrix)=unique(fig_m3$V1)
windown_bed <- read.table("tair10_win500.bed")
colnames(windown_bed) <- c("chrom","start","end","winID")
fig_m1_matrix1 <- merge(windown_bed,fig_m1_matrix,by.x=4,by.y=0,all.y=T) 
fig_m2_matrix1 <- merge(windown_bed,fig_m2_matrix,by.x=4,by.y=0,all.y=T)
fig_m3_matrix1 <- merge(windown_bed,fig_m3_matrix,by.x=4,by.y=0,all.y=T)
write.table(fig_m1_matrix1,file="/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant1_matrix_fig.tsv",quote =F,sep="\t",row.names = T,col.names = T)
write.table(fig_m2_matrix1,file="/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant2_matrix_fig.tsv",quote =F,sep="\t",row.names = T,col.names = T)
write.table(fig_m3_matrix1,file="/var/data/06.Richa_sRNA/smallRNA/8.DEsRNA/win_mutant3_matrix_fig.tsv",quote =F,sep="\t",row.names = T,col.names = T)

pheatmap(fig_m3_matrix,cluster_cols = F,show_rownames=F,main="sRNA ratio")
dev.off()
##DEG analysis
dcl3_col0 <- c("s1","s2","s3","s4","s5","s6","s7","s8")
dcl234_col0 <- c("s1","s2","s3","s4","s9","s10","s11","s12")
hen1_col0 <- c("s1","s2","s3","s4","s13","s14","s15","s16")

count_matrix_col_dcl3<- sRNA_5[,dcl3_col0]
head(counts_matrix_col_dcl3_10)
count_matrix_col_dcl3 <- count_matrix_col_dcl3/sRNA_5$Hit
count_matrix_col_dcl3_10 <- apply(count_matrix_col_dcl3,2,function(x){ifelse(x>=10,1,0)})
counts_matrix_col_dcl3_10 <- count_matrix_col_dcl3[apply(count_matrix_col_dcl3_10,1,sum)>=1,]
count <- counts_matrix_col_dcl3_10
group <- c(1,1,1,1,2,2,2,2)
y <- DGEList(counts = count, group = group)

##filter
minGroupSize <- min(table(group))
keep <- rowSums(cpm(y) > 1) >= minGroupSize
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
##test
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=T)
y <- estimateTagwiseDisp(y)
##DE
p_cutoff <- 0.05
fc <- 2
fc_cutoff <- log2(fc)
et <- exactTest(y, pair=c(1, 2))
de <- decideTestsDGE(et, p.value=p_cutoff, lfc=fc_cutoff, adjust.method='fdr')
FDR <- sprintf('%.3e', p.adjust(et$table$PValue, method='fdr'))
a <- cbind(et$table,FDR,de)
de_up <- a[a$de>0,]
de_down <- a[a$de<0,]
name <- sRNA_5[,1:10]
DE_table_name<- name[rownames(a),]
DE_table_count <- counts_matrix_col_dcl3_10[rownames(a),]
Col0_vs_dcl3_DE_table<-cbind(DE_table_name,DE_table_count,a)
write.table(Col0_vs_dcl3_DE_table ,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_col0_dcl3_10_DE.tsv",sep="\t",quote = FALSE)

##dcl234
count_matrix_col_dcl234<- sRNA_5[,dcl234_col0]
head(counts_matrix_col_dcl234_10)
count_matrix_col_dcl234 <- count_matrix_col_dcl234/sRNA_5$Hit
count_matrix_col_dcl234_10 <- apply(count_matrix_col_dcl234,2,function(x){ifelse(x>=10,1,0)})
counts_matrix_col_dcl234_10 <- count_matrix_col_dcl234[apply(count_matrix_col_dcl234_10,1,sum)>=1,]
count <- counts_matrix_col_dcl234_10
group <- c(1,1,1,1,2,2,2,2)
y <- DGEList(counts = count, group = group)
count
##filter
minGroupSize <- min(table(group))
keep <- rowSums(cpm(y) > 1) >= minGroupSize
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
##test
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=T)
y <- estimateTagwiseDisp(y)
##DE
p_cutoff <- 0.05
fc <- 2
fc_cutoff <- log2(fc)
et <- exactTest(y, pair=c(1, 2))
de <- decideTestsDGE(et, p.value=p_cutoff, lfc=fc_cutoff, adjust.method='fdr')
FDR <- sprintf('%.3e', p.adjust(et$table$PValue, method='fdr'))
a <- cbind(et$table,FDR,de)
de_up <- a[a$de>0,]
de_down <- a[a$de<0,]
name <- sRNA_5[,1:10]
DE_table_name<- name[rownames(a),]
DE_table_count <- counts_matrix_col_dcl234_10[rownames(a),]
Col0_vs_dcl234_DE_table<-cbind(DE_table_name,DE_table_count,a)
write.table(Col0_vs_dcl234_DE_table ,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_col0_dcl234_10_DE.tsv",sep="\t",quote = FALSE)

###hen1
count_matrix_col_hen1<- sRNA_5[,hen1_col0]
count_matrix_col_hen1 <- count_matrix_col_hen1/sRNA_5$Hit
count_matrix_col_hen1_10 <- apply(count_matrix_col_hen1,2,function(x){ifelse(x>=10,1,0)})
counts_matrix_col_hen1_10 <- count_matrix_col_hen1[apply(count_matrix_col_hen1_10,1,sum)>=1,]
count <- counts_matrix_col_hen1_10
group <- c(1,1,1,1,2,2,2,2)
y <- DGEList(counts = count, group = group)
count
##filter
minGroupSize <- min(table(group))
keep <- rowSums(cpm(y) > 1) >= minGroupSize
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
##test
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=T)
y <- estimateTagwiseDisp(y)
##DE
p_cutoff <- 0.05
fc <- 2
fc_cutoff <- log2(fc)
et <- exactTest(y, pair=c(1, 2))
de <- decideTestsDGE(et, p.value=p_cutoff, lfc=fc_cutoff, adjust.method='fdr')
FDR <- sprintf('%.3e', p.adjust(et$table$PValue, method='fdr'))
a <- cbind(et$table,FDR,de)
de_up <- a[a$de>0,]
de_down <- a[a$de<0,]
name <- sRNA_5[,1:10]
DE_table_name<- name[rownames(a),]
DE_table_count <- counts_matrix_col_hen1_10[rownames(a),]
Col0_vs_hen1_DE_table<-cbind(DE_table_name,DE_table_count,a)
write.table(Col0_vs_hen1_DE_table ,file="/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_col0_hen1_10_DE.tsv",sep="\t",quote = FALSE)

