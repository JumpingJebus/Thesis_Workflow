source("https://bioconductor.org/biocLite.R")
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
##########################################################################
####################### CREATING DDS #####################################
##########################################################################

dir="/home/stephen/Downloads/UNEWER/two/"
dds_unstr= as.matrix(read.csv(file.path(dir,'gene_count_matrix.csv'), row.names=1), header= TRUE)
dds_unstr
sample_unstr= data.frame(row.names=c("SRR3374051","SRR3374052","SRR3374053","SRR3374054",
                                     "SRR3374055","SRR3374056","SRR3374057","SRR3374058",
                                     "SRR3374059","SRR3374060","SRR3374061","SRR3374062",
                                     "SRR3374063","SRR3374064","SRR3374065","SRR3374066"),
                         condition= as.factor(c("un","un","un","un",
                                                "six","six","six","six",
                                                "eya","eya","eya","eya",
                                                "sixeya","sixeya","sixeya","sixeya")),
                         treatment= as.factor(c("chx","chx","dex","dex",
                                          "chx","chx","dex","dex",
                                          "chx","chx","dex","dex",
                                          "chx","chx","dex","dex")))
sample_unstr$condition= relevel(sample_unstr$condition, ref= "un")
sample_unstr$treatment= relevel(sample_unstr$treatment, ref= "chx")
sample_unstr
all(rownames(sample_unstr) %in% colnames(dds_unstr))
dds_unstr= dds_unstr[, rownames(sample_unstr)]
all(rownames(sample_unstr)== colnames(dds_unstr))
dds_unstr
dds_unstr= DESeqDataSetFromMatrix(countData = dds_unstr, colData = sample_unstr, design= ~treatment + condition + condition:treatment)
dds_unstr

dds <- estimateSizeFactors(dds_unstr)
dds
nrow(dds_unstr)
idx <- rowSums( counts(dds, normalized=TRUE) >= 6 ) >= 8
dds <- dds_unstr[idx,]
nrow(dds)
dds_unstr <- DESeq(dds)
rld= rlog(dds_unstr, blind = F)
unpca=rld
peeca= plotPCA(rld[,5:16], intgroup= c("condition","treatment"), returnData= T)
percentVarStr= round(100 * attr(peeca, "percentVar"))

head(peeca)
ggplot(peeca, aes(PC1,PC2, color=group, shape= group)) +
  geom_point(size= 7) +
  xlab(paste0("PC1: ", percentVarStr[1], "% variance")) +
  ylab(paste0("PC2: ", percentVarStr[2], "% variance")) +
  theme(legend.text=element_text(size=20))
#DISTANCE HEATMAP
dists = dist(t(assay(rld[,5:16])))
mat = as.matrix(dists)
colnames(mat) = c("six_chx","six_chx","six_dex","six_dex",
                  "eya_chx","eya_chx","eya_dex","eya_dex",
                  "sixeya_chx","sixeya_chx","sixeya_dex","sixeya_dex")
rownames(mat) = colnames(mat)
hc = hclust(dists)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap.2(mat, symm = TRUE,
          trace = "none", col = rev(hmcol), margin = c(9, 9), dendrogram = "both",cexCol = 1.5,cexRow = 1.5)

res= results(dds_unstr)

##########################################################################
####################### SIX1 ANALYSIS ####################################
##########################################################################


resultsNames(dds_unstr)

res_six= results(dds_unstr, contrast= list(c("treatmentdex.conditionsix")))
resOrdered_six= res_six[order(res_six$padj),]
resOrdered_six
head(resOrdered_six)
sum(res_six$padj < 0.05, na.rm = TRUE)
sumsix= sum(res_six$padj < 0.05, na.rm = TRUE)
Sigsix= (resOrdered_six[1:sumsix,])
Sigsix
write.csv(as.data.frame(Sigsix), file= "/home/stephen/Downloads/UNEWER/two/SigSix.csv")
write.csv(as.data.frame(resOrdered_six), file= "/home/stephen/Downloads/UNEWER/two/six.csv")
plotMA(res_six, ylim = c(-15, 15))
#geminin
plotCounts(dds_unstr, gene= "MSTRG.27792", intgroup= c("condition", "treatment"))
plotCounts(dds_unstr, gene = which.min(res_six$padj), intgroup= c("condition", "treatment"))
plotCounts(dds_unstr, gene = which.max(res_six$log2FoldChange), intgroup= c("condition", "treatment"))


mat= assay(rld[,5:16])
head(mat)
sigheat_six= rownames(Sigsix)[which(Sigsix$padj < 0.05)]
pheatmap(mat[sigheat_six[1:1000],],
         colorRampPalette(c("green", "black", "red"))(100),
         show_rownames = F,cluster_cols = F,
         labels_col= c("six-chx","six-chx","six-dex","six-dex",
                       "eya-chx","eya-chx","eya-dex","eya-dex",
                       "sixeya-chx","sixeya-chx","sixeya-dex","sixeya-dex"),
         scale= "row",
         fontsize_col = 15)
##########################################################################
####################### EYA1 ANALYSIS ####################################
##########################################################################


resultsNames(dds_unstr)
res_eya= results(dds_unstr, contrast= list(c("treatmentdex.conditioneya")))
res_eya
resOrdered_eya= res_eya[order(res_eya$padj),]
resOrdered_eya
summary(res_eya)
sum(res_eya$padj < 0.05, na.rm = TRUE)
sumeya= sum(res_eya$padj < 0.05, na.rm = TRUE)
Sigeya= (resOrdered_eya[1:sumeya,])
Sigeya
write.csv(as.data.frame(Sigeya), file= "/home/stephen/Downloads/UNEWER/two/SigEya.csv")
write.csv(as.data.frame(resOrdered_eya), file= "/home/stephen/Downloads/UNEWER/two/eya.csv")
plotMA(res_eya, ylim = c(-15, 15))

plotCounts(dds_unstr, gene = which.min(res_eya$padj), intgroup= c("condition", "treatment"))
plotCounts(dds_unstr, gene = which.max(res_eya$log2FoldChange), intgroup= c("condition", "treatment"))
plotCounts(dds_unstr, gene = which.min(res_eya$log2FoldChange), intgroup= c("condition", "treatment"))

mat= assay(rld)
mat= mat
sigheat= rownames(Sigeya)[which(Sigeya$padj < 0.05)]
pheatmap(mat[sigheat[1:50],],cluster_cols = F,
         color = colorRampPalette(c("green", "black", "red"))(100),
         show_rownames = F,
         labels_col= c("un-chx","un-chx","un-dex","un-dex",
                       "six-chx","six-chx","six-dex","six-dex",
                        "eya-chx","eya-chx","eya-dex","eya-dex",
                        "sixeya-chx","sixeya-chx","sixeya-dex","sixeya-dex"),
         scale="row",
         fontsize_col= 15)



##########################################################################
####################### SIX1EYA1 ANALYSIS ################################
##########################################################################


#SIXEYA RESULTS GENERATION
res_sixeya= results(dds_unstr, contrast= list(c("treatmentdex.conditionsixeya")))

#ORDER SIXEYA BY PADJ
resOrdered_sixeya= res_sixeya[order(res_sixeya$padj),]
write.csv(as.data.frame(resOrdered_sixeya), file= "/home/stephen/Downloads/UNEWER/two/sixeya.csv")



#SUBSET OUT THE SIGNIFICANT GENES
sum(res_sixeya$padj < 0.05)
sumsixeya= sum(res_sixeya$padj < 0.05,na.rm=TRUE)
sumsixeya
SigSixEya= (resOrdered_sixeya[1:sumsixeya,])
write.csv(as.data.frame(SigSixEya), file= "/home/stephen/Downloads/UNEWER/two/SigSixEya.csv")

#MA PLOT OF SIXEYA
plotMA(res_sixeya, ylim = c(-15, 15))

#PLOT COUNTS FOR SIXEYA
plotCounts(dds_unstr, gene = which.min(res_sixeya$padj), intgroup= c("condition", "treatment"))
plotCounts(dds_unstr, gene = which.max(res_sixeya$log2FoldChange), intgroup= c("condition", "treatment"))
plotCounts(dds_unstr, gene = which.min(res_sixeya$log2FoldChange), intgroup= c("condition", "treatment"))
plotCounts(dds_unstr, gene = "MSTRG.18229", intgroup= c("condition", "treatment"))


#HEATMAP FOR SIXEYA
mat= assay(rld[,5:16])
sigheat_sixeya= rownames(SigSixEya)[which(SigSixEya$padj < 0.05)]
pheatmap(mat[sigheat_sixeya[1:2000],] ,
         color = colorRampPalette(c("green", "black", "red"))(100),
         cluster_cols = F,labels_col= c("six-chx","six-chx","six-dex","six-dex",
                                        "eya-chx","eya-chx","eya-dex","eya-dex",
                                        "sixeya-chx","sixeya-chx","sixeya-dex","sixeya-dex"),
         scale="row",
         fontsize_col= 15,
         show_rownames = F)
