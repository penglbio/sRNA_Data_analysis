data <- read.table("/var/data/Richa_sRNA/smallRNA/7.table/sRNA_counts_win_w_header.tsv",header = T)
library(dplyr)
data1<-data[,11:26]
head(data1)
sRNA <- apply(data1,2,function(x){ifelse(x>=5,1,0)})
?apply
sRNA_5 <- data[apply(sRNA,1,sum)>=1,]
write.table(sRNA_5,file="/var/data/Richa_sRNA/smallRNA/7.table/sRNA_counts_win_w_header_counts_up_5.tsv",sep="\t",quote = FALSE)
head(sRNA_5)
sRNA_5_sum_sample <- apply(sRNA_5[,11:26],2,sum)
sRNA_5_sum_sample
sRNA_5_RPTM <- apply(sRNA_5[,11:26],2,function(x){x/sRNA_5_sum_sample*10000000})
colnames(sRNA_5_RPTM)
colnames(sRNA_5_RPTM)<-paste0("RPTM_",colnames(sRNA_5[,11:26]))
sRNA_counts_RPTM <- cbind(sRNA_5,sRNA_5_RPTM)
head(sRNA_counts_RPTM)
write.table(sRNA_counts_RPTM,file="/var/data/Richa_sRNA/smallRNA/7.table/sRNA_counts_RPTM_win_w_header_counts_up_5.tsv",sep="\t",quote = FALSE)
sRNA_5_HNA <- apply(sRNA_5_RPTM,2,function(x){x/sRNA_5[,6]})
sum_w_m<-apply(sRNA_5_HNA,1,sum)
colnames(sRNA_5_HNA)<-paste0("HNA_",colnames(sRNA_5[,11:26]))
sRNA_counts_RPTM_HNA <- cbind(sRNA_counts_RPTM,sRNA_5_HNA)
RPTM <- cbind(sRNA_counts_RPTM[,1:10],sRNA_5_RPTM)
HNA <- cbind(sRNA_counts_RPTM[,1:10],sRNA_5_HNA)
write.table(RPTM,file="/var/data/Richa_sRNA/smallRNA/7.table/sRNA_RPTM_win_w_header_counts_up_5.tsv",sep="\t")
write.table(HNA,file="/var/data/Richa_sRNA/smallRNA/7.table/sRNA_HNA_win_w_header_counts_up_5.tsv",sep="\t",row.names = FALSE,quote = FALSE)

sRNA_counts_RPTM_HNA_cutoff_400<-sRNA_counts_RPTM_HNA[sum_w_m>=400,]
win_HNA_sum <- group_by(sRNA_counts_RPTM_HNA_cutoff_400,win_id) %>% summarise(test1=sum(HNA_s1),test2=sum(HNA_s2),test3=sum(HNA_s3),test4=sum(HNA_s4),test5=sum(HNA_s5),
                                                               test6=sum(HNA_s6),test7=sum(HNA_s7),test8=sum(HNA_s8),test9=sum(HNA_s9),test10=sum(HNA_s10),
                                                               test11=sum(HNA_s11),test12=sum(HNA_s12),test13=sum(HNA_s13),test14=sum(HNA_s14),
                                                               test15=sum(HNA_s15),test16=sum(HNA_s16)) 
win_HNA_sum <- as.data.frame(win_HNA_sum)
rownames(win_HNA_sum) <- win_HNA_sum[,1]
win_HNA_sum <- win_HNA_sum[2:17]
head(win_HNA_sum)
widetype_win_HNA_sum_mean <- apply(win_HNA_sum[,1:4],1,mean)
mutant1_win_HNA_sum_mean <- apply(win_HNA_sum[,5:8],1,mean)
mutant2_win_HNA_sum_mean <- apply(win_HNA_sum[,9:12],1,mean)
mutant3_win_HNA_sum_mean <- apply(win_HNA_sum[,13:16],1,mean)
win_HNA_sample_mean <- cbind(widetype_win_HNA_sum_mean,mutant1_win_HNA_sum_mean,mutant2_win_HNA_sum_mean,mutant3_win_HNA_sum_mean)
 
FC_mutant_vs_widetype <- apply(win_HNA_sample_mean[,2:4], 2, function(x){win_HNA_sample_mean[,1]/x} )
dim(FC_mutant_vs_widetype[FC_mutant_vs_widetype[,3]<=0.5,])
mutant1_FC2_winid <-c(rownames(FC_mutant_vs_widetype[FC_mutant_vs_widetype[,1]>=2,]),rownames(FC_mutant_vs_widetype[FC_mutant_vs_widetype[,1]<=0.5,]))
mutant2_FC2_winid <-c(rownames(FC_mutant_vs_widetype[FC_mutant_vs_widetype[,2]>=2,]),rownames(FC_mutant_vs_widetype[FC_mutant_vs_widetype[,2]<=0.5,]))
mutant3_FC2_winid <- c(rownames(FC_mutant_vs_widetype[FC_mutant_vs_widetype[,3]>=2,]),rownames(FC_mutant_vs_widetype[FC_mutant_vs_widetype[,3]<=0.5,]))
length(mutant1_FC2_winid)
length(mutant2_FC2_winid)
length(mutant3_FC2_winid)
FC_mutant_vs_widetype[,1]>=2


b <- c(1,2,3,4)
b>2
a<-data_frame(x=c(1,2,3,4),x2=c(2,3,4,5))
a[b>2,]


