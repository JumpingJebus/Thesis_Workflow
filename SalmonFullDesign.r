library(DESeq2)
library(tximport)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
#set the dir to the location of the files that will be used.
dir = "/home/stephen/Downloads"
#create a table in the variable samplesUn which contains the contents of the sampleUN.txt file
samples = read.table(file.path(dir,"samples.txt"))
#add a column to this table called condition which contains the factors relatng to the treatments
samples$condition= factor(c("UN","UN","UN","UN",
                            "six","six","six","six",
                            "eya","eya","eya","eya",
                            "sixeya","sixeya","sixeya","sixeya"))
samples$treatment= as.factor(c("chx","chx","dex","dex",
                                "chx","chx","dex","dex",
                                "chx","chx","dex","dex",
                                "chx","chx","dex","dex"))
#rename the rownames after the first column in the sampleUn variable which is the run IDs
rownames(samples)= samples$V1
#name the columns
samples
colnames(samples)= c("run","condition","treatment")
#create a variable which contains the path to te quants1.sf files which will be used
dir2= "/home/stephen/Downloads/quants"
files= file.path(dir2,"salmon",samples$run,"quants1.sf")
#name each file path after the run ID
names(files)= samples$run
files
#load in the transcript to gene ID file
tx2gene= read.csv(file.path(dir,"quants","t2gone.csv"))
#run tximport and assign the results to a variable
txi= tximport(files,type="salmon",tx2gene=tx2gene)
samples
#create a deseq dataset from the tximport file. Give it the column data from samplesUn and the deisgn data from the condition column
ddsTxi= DESeqDataSetFromTximport( txi, colData= samples, design= ~condition + treatment + treatment:condition)
#count the number of rows in the dataset
nrow(ddsTxi)
#filter out very low/no read count rows
ddsTxi$condition= relevel(ddsTxi$condition, ref= "UN")
ddsTxi$treatment= relevel(ddsTxi$treatment, ref= "chx")
dds <- estimateSizeFactors(ddsTxi)
dds
nrow(ddsTxi)
idx <- rowSums( counts(dds, normalized=TRUE) >= 6 ) >= 8
dds <- ddsTxi[idx,]
nrow(dds)

ddsTxi <- DESeq(dds)

res= results(ddsTxi)
plotCounts(ddsTxi, gene = which.min(res$padj), intgroup= c("condition","treatment"))


rldSal= rlog(ddsTxi, blind= FALSE)
plotPCA(rldSal[,5:16], intgroup= c("condition","treatment"))
#PCA of all treatments
plotPCA(rldSal[,c(3:4,7:8,11:12,15:16)], intgroup= c("condition","treatment"))
#PCA of all controls
plotPCA(rldSal[,c(1:2,5:6,9:10,13:14)], intgroup= c("condition","treatment"))

peecaSal= plotPCA(rldSal[,5:16], intgroup= c("condition","treatment"), returnData= T)
percentVarSal= round(100 * attr(peecaSal, "percentVar"))

ggplot(peecaSal, aes(PC1,PC2, color=group, shape= group)) +
  geom_point(size= 7) +
  xlab(paste0("PC1: ", percentVarSal[1], "% variance")) +
  ylab(paste0("PC2: ", percentVarSal[2], "% variance")) +
  theme(legend.text=element_text(size=20))

distsSal = dist(t(assay(rldSal[,5:16])))
matSal = as.matrix(distsSal)
colnames(matSal) = c("six_chx","six_chx","six_dex","six_dex",
                  "eya_chx","eya_chx","eya_dex","eya_dex",
                  "sixeya_chx","sixeya_chx","sixeya_dex","sixeya_dex")
rownames(matSal) = colnames(matSal)
hc = hclust(distsSal)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap.2(matSal, symm = TRUE,
          trace = "none", col = rev(hmcol), margin = c(9, 9), dendrogram = "both",cexCol = 1.5,cexRow = 1.5)


#get the results by using results()
res_six1= results(ddsTxi, contrast = list(c("conditionsix.treatmentdex")))
#view the result names
resultsNames(ddsTxi)
#order the result by the most significant p adjusted values
resOrdered_six1= res_six1[order(res_six1$padj),]
#view the summary of the results
summary(res_six1)
#view the summary of the significant results
sum(res_six1$padj <0.05, na.rm= TRUE)
write.csv(as.data.frame(resOrdered_six1), file= "/home/stephen/Salmon_deseq2/salmon_six.csv")

sumsix1= sum(res_six1$padj <0.05, na.rm= TRUE)
Sigsix1= (resOrdered_six1[1:sumsix1,])
write.csv(as.data.frame(Sigsix1), file= "/home/stephen/Salmon_deseq2/SigSix.csv")
tail(resOrdered_six1,3100)
plotMA(res_six1, ylim = c(-15, 15))
plotCounts(ddsTxi, gene = which.min(res_six1$padj), intgroup= c("condition","treatment"))
plotCounts(ddsTxi, gene = "notch1.L", intgroup= c("condition","treatment"))
plotCounts(ddsTxi, gene = "Xelaev18022529m.g", intgroup= c("condition","treatment"))
plotCounts(ddsTxi, gene = which.max(res_six1$log2FoldChange), intgroup= c("condition","treatment"))

matsixSal= assay(rldSal[,5:16])

sigheat_sixSal= rownames(Sigsix1)[which(Sigsix1$padj < 0.05)]
pheatmap(matsixSal[sigheat_sixSal[1:1100],],
         colorRampPalette(c("green", "black", "red"))(100),
         show_rownames = F,cluster_cols = F,
         labels_col= c("six-chx","six-chx","six-dex","six-dex",
                       "eya-chx","eya-chx","eya-dex","eya-dex",
                       "sixeya-chx","sixeya-chx","sixeya-dex","sixeya-dex"),
         scale= "row",
         fontsize_col = 15)



res_eya1= results(ddsTxi, contrast = list(c("conditioneya.treatmentdex")))
resultsNames(ddsTxi)
resOrdered_eya1= res_eya1[order(res_eya1$padj),]
summary(res_eya1)
sum(res_eya1$padj <0.05, na.rm= TRUE)
write.csv(as.data.frame(resOrdered_eya1), file= "/home/stephen/Salmon_deseq2/salmon_eya.csv")

sumeya1= sum(res_eya1$padj < 0.05, na.rm = TRUE)
Sigeya1= (resOrdered_eya1[1:sumeya1,])
write.csv(as.data.frame(Sigeya1), file= "/home/stephen/Salmon_deseq2/SigEya.csv")
plotMA(res_eya1, ylim = c(-15, 15))

plotCounts(ddsTxi, gene = which.min(res_eya1$padj), intgroup= c("condition", "treatment"))
plotCounts(ddsTxi, gene = which.max(res_eya1$log2FoldChange), intgroup= c("condition", "treatment"))
plotCounts(ddsTxi, gene = which.min(res_eya1$log2FoldChange), intgroup= c("condition", "treatment"))
plotCounts(ddsTxi, gene = "wnt5b.S", intgroup= c("condition", "treatment"))

mateyaSal= assay(rldSal[,5:16])

sigheat_eyaSal= rownames(Sigeya1)[which(Sigeya1$padj < 0.05)]
pheatmap(matsixSal[sigheat_eyaSal[1:2100],],
         colorRampPalette(c("green", "black", "red"))(100),
         show_rownames = F,cluster_cols = F,
         labels_col= c("six-chx","six-chx","six-dex","six-dex",
                       "eya-chx","eya-chx","eya-dex","eya-dex",
                       "sixeya-chx","sixeya-chx","sixeya-dex","sixeya-dex"),
         scale= "row",
         fontsize_col = 15)


res_sixeya1= results(ddsTxi, contrast = list(c("conditionsixeya.treatmentdex")))
resultsNames(dds)
resOrdered_sixeya1= res_sixeya1[order(res_sixeya1$padj),]
summary(res_sixeya1)
sum(res_sixeya1$padj <0.05, na.rm= TRUE)
write.csv(as.data.frame(resOrdered_sixeya1), file= "/home/stephen/Documents/salmon_sixeya.csv")

#SUBSET OUT THE SIGNIFICANT GENES
sumsixeya1= sum(res_sixeya1$padj < 0.05,na.rm=TRUE)
sumsixeya1
SigSixeya1= (resOrdered_sixeya1[1:sumsixeya1,])
write.csv(as.data.frame(SigSixeya1), file= "/home/stephen/Salmon_deseq2/SigSixEya.csv")

#MA PLOT OF SIXEYA
plotMA(res_sixeya1, ylim = c(-15, 15))

#PLOT COUNTS FOR SIXEYA
plotCounts(ddsTxi, gene = which.min(res_sixeya1$padj), intgroup= c("condition", "treatment"))
plotCounts(ddsTxi, gene = which.max(res_sixeya1$log2FoldChange), intgroup= c("condition", "treatment"))
plotCounts(ddsTxi, gene = which.min(res_sixeya1$log2FoldChange), intgroup= c("condition", "treatment"))
plotCounts(ddsTxi, gene = "cltc.S", intgroup= c("condition", "treatment"))

matsixeyaSal= assay(rldSal[,5:16])

sigheat_sixeyaSal= rownames(SigSixeya1)[which(SigSixeya1$padj < 0.05)]
pheatmap(matsixeyaSal[sigheat_sixeyaSal[1:2200],],
         colorRampPalette(c("green", "black", "red"))(100),
         show_rownames = F,cluster_cols = F,
         labels_col= c("six-chx","six-chx","six-dex","six-dex",
                       "eya-chx","eya-chx","eya-dex","eya-dex",
                       "sixeya-chx","sixeya-chx","sixeya-dex","sixeya-dex"),
         scale= "row",
         fontsize_col = 15)

