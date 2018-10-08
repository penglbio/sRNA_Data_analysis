##length destribution PCA 
getwd()
data <- read.table("/var/data/06.Richa_sRNA/smallRNA/7.table/sRNA_counts_win_w_header.tsv",header = T,stringsAsFactors = F)
sample <- read.table("/var/data/06.Richa_sRNA/smallRNA/7.table/sample.tsv",header = T,stringsAsFactors = F)
outfig <- "/var/data/06.Richa_sRNA/smallRNA/9.fig/"
sample_name <- paste(sample$Sample_ID,sample$replicates,sep="_")
library(dplyr)
library("RColorBrewer")
library(ggplot2)
data1<-data[,11:26]
data1 <- cbind(data1,data[,5])
colnames(data1)[17] <- "length"
head(data2)
data2 <- aggregate(data1[1:16],by=list(data1[,17]),FUN=sum)
rownames(data2) <- data2[,1] 
data2 <- data2[2:17]
data2_ratio <- apply(data2,2,function(x){x/sum(x)*100})
color=brewer.pal(4,"Set3")
plot(data2_ratio[,1],type="l",xaxt="n",main =  "sRNA length destribution",xlab = "length",ylab = "Percentage(%)",bty="n",col=color[1])
axis(1,at=1:43,labels = 18:60)
pdf(paste0(outfig,"sRNA length destribution_rep_separ.pdf"),width = 12,height = 5)
par(xpd=T,mar)
?par
for(i in 1:4){
    barplot(t(data2_ratio[,((i-1)*4+1):(i*4)]),beside = T,col = color,main = paste(unique(sample$Sample_ID)[i],"sRNA length destribution")
            ,xlab="length(nt)",ylab="Percentage(%)",border='NA')
    legend(142,20,sample_name[((i-1)*4+1):(i*4)],bty="n",pch=16,col=color)
}
dev.off()

pdf(paste0(outfig,"sRNA length destribution.pdf"),width = 12,height = 5)
for(i in 1:16){
  barplot(t(data2_ratio[,i]),beside = T,col = "grey",main = paste(sample_name[i],"sRNA length destribution")
          ,xlab="length",ylab="Percentage(%)")
}
dev.off()

###PCA

sRNA_5 <- apply(data1[,1:16],2,function(x){ifelse(x>=5,1,0)})
data <- data[,11:26]
sRNA_5 <- data[apply(sRNA_5,1,sum)>=1,]

sRNA_RTPM_5 <- apply(sRNA_5,2,function(x){x/sum(x)*10000000})
sRNA_RTPM_5_lg <- log2(sRNA_RTPM_5+1)
pca <- prcomp(t(sRNA_RTPM_5_lg), center = T, scale. = T)
x <- pca$x
prop_var <- round(summary(pca)$importance[2,1:2]*100,0)
set2_cols <- brewer.pal(4, 'Set1')
cols <- rep(set2_cols, each=4)
pchs <- rep(16:19,each=4)
pdf(paste0(outfig,"pca.pdf"), height =7, width=7)
# layout(matrix(c(1,2),nrow=1), wid=c(5, 2))
# par(mar=c(5,4,4,0))
par(xpd=TRUE)
plot(x[,1], x[,2],xlim=range(x[,1])*1.1, ylim=range(x[,2])*1.1, col=cols, cex=2,pch=pchs,
     xlab=paste0('PC1 (',prop_var[1],'% of Variance)'),
     ylab=paste0('PC2 (',prop_var[2],'% of Variance)')
)
legend(700,1000,unique(sample$Sample_ID),pch = 16:19,col = set2_cols,bty="n")
dev.off()
getwd()
itg<-read.table("/var/data/06.Richa_sRNA/smallRNA/12.sRNA_QC_data/intergenetic_sRNA.tsv")
gene<-read.table("/var/data/06.Richa_sRNA/smallRNA/12.sRNA_QC_data/gene_sRNA.tsv")
TE<-read.table("/var/data/06.Richa_sRNA/smallRNA/12.sRNA_QC_data/TE_sRNA.tsv")
pseudo<-read.table("/var/data/06.Richa_sRNA/smallRNA/12.sRNA_QC_data/pseudogene_sRNA.tsv")

gene_counts<-colSums(gene[,11:26])
itg_counts<-colSums(itg[,5:20])
TE_counts<-colSums(TE[,5:20])
pseudo_counts<-colSums(pseudo[,11:26])
all_couns<- colSums(data)
gene_counts
itg_counts
TE_counts
pseudo_counts
names(all_couns)<-names(gene_counts)
others<-all_couns-rowSums(cbind(gene_counts,pseudo_counts,TE_counts,itg_counts))
classify_table<-cbind(gene_counts,pseudo_counts,TE_counts,itg_counts,others)
rownames(classify_table)<-1:16
pdf("sRNA_in_genome_region_classify.pdf",height = 5,width = 10)
colors=brewer.pal(5, 'Set1')
par(pdx=T,mar=c(3,4,4,6))
for( i in 1:16){
  pie(classify_table[i,],main=paste0("sample",i,"_classify"),init.angle = 90,col=colors,border = "white",labels = NA)
  legend(1,-0.2,legend = c("gene","pseudogene","TE","intergenic","other"),pch=16,col=colors,bty = "n")
}
dev.off()

?legend

