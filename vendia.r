install.packages('VennDiagram')
library(VennDiagram)

#load in files
st_eya= read.csv(file= "/home/stephen/Downloads/UNEWER/two/StringtieGenes/StringTieEya.csv", sep=",",stringsAsFactors = F)
sa_eya= read.csv(file= "/home/stephen/Salmon_deseq2/Sal_eya_fix.csv",sep= ",", stringsAsFactors = F)

a_eya= merge(x= sa_eya, y= st_eya, by= "X")

comb_eya= length(a_eya$X)
grid.newpage()
draw.pairwise.venn(length(sa_eya$X),length(st_eya$X),
                   cross.area=comb_eya, category= c("Salmon Eya1","StringTie Eya1"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(200,150),cat.dist=rep(0.035,2),
                   cex= 2, cat.cex=1.5)

st_six= read.csv(file= "/home/stephen/Downloads/UNEWER/two/StringtieGenes/StringTieSix.csv", sep=",",stringsAsFactors = F)
sa_six= read.csv(file= "/home/stephen/Salmon_deseq2/Sal_six_fix.csv",sep= ",", stringsAsFactors = F)

a_six= merge(x= sa_six, y= st_six, by= "X")

comb_six= length(a_six$X)
grid.newpage()
draw.pairwise.venn(length(sa_six$X),length(st_six$X),
                   cross.area=comb_six, category= c("Salmon Six1","StringTie Six1"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(200,150),cat.dist=rep(0.035,2),
                   cex=2, cat.cex= 1.5)

st_sixeya= read.csv(file= "/home/stephen/Downloads/UNEWER/two/StringtieGenes/StringTieSixEya.csv", sep=",",stringsAsFactors = F)
sa_sixeya= read.csv(file= "/home/stephen/Salmon_deseq2/Sal_fix_SixEya.csv",sep= ",", stringsAsFactors = F)

a_sixeya= merge(x= sa_sixeya, y= st_sixeya, by= "X") 

comb_sixeya= length(a_sixeya$X)
grid.newpage()
draw.pairwise.venn(length(sa_sixeya$X),length(st_sixeya$X),
                   cross.area=comb_sixeya, category= c("Salmon Six1+Eya1","StringTie Six1+Eya1"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(200,150),cat.dist=rep(0.045,2),
                   cex = 2,cat.cex=1.5)


sa_gene_no= 45038
st_gene_no= 62283
perc_six_st= 100/(st_gene_no/length(st_six$X))
perc_six_st= round(perc_six_st,digits=2)
perc_eya_st= 100/(st_gene_no/length(st_eya$X))
perc_eya_st= round(perc_eya_st,digits=2)
perc_sixeya_st= 100/(st_gene_no/length(st_sixeya$X))
perc_sixeya_st= round(perc_sixeya_st,digits=2)

perc_six_sa= 100/(45038/length(sa_six$X))
perc_six_sa= round(perc_six_sa,digits=2)
perc_eya_sa= 100/(sa_gene_no/length(sa_eya$X))
perc_eya_sa= round(perc_eya_sa, digits=2)
perc_sixeya_sa= 100/(sa_gene_no/length(sa_sixeya$X))
perc_sixeya_sa= round(perc_sixeya_sa,digits=2)

o= matrix(c(st_gene_no,sa_gene_no,perc_six_st,perc_six_sa,
          perc_eya_st,perc_eya_sa,
          perc_sixeya_st,perc_sixeya_sa),byrow=T,nrow=4)

colnames(o)= c("StringTie","Salmon")
rownames(o)= c("Total Genes", "Six1","Eya1","Six1+Eya1")
o= as.data.frame(o)



st_e_up= read.table("/home/stephen/eyaupgenes.txt",stringsAsFactors = F)
sa_e_up= read.table("/home/stephen/Salmon_deseq2/Sal_eyaupgenes.txt",stringsAsFactors = F)
eyaup= merge(x= st_e_up, y= sa_e_up, by= "V1")
write.table(eyaup, file="/home/stephen/EyaUpcross.txt",row.names = F, quote = F)

grid.newpage()
draw.pairwise.venn(length(sa_e_up$V1),length(st_e_up$V1),
                   cross.area=length(eyaup$V1), category= c("Salmon Eya1 Upregulated","StringTie Eya1 Upregulated"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(200,150),cat.dist=rep(0.035,2),
                   cex=2, cat.cex= 1.5)

st_e_d= read.table("/home/stephen/eyadowngenes.txt",stringsAsFactors = F)
sa_e_d= read.table("/home/stephen/Salmon_deseq2/Sal_eyadowngenes.txt",stringsAsFactors = F)
eyadown= merge(x= st_e_d, y= sa_e_d, by= "V1")
write.table(eyadown, file="/home/stephen/EyaDowncross.txt",row.names = F, quote = F)

grid.newpage()
draw.pairwise.venn(length(sa_e_d$V1),length(st_e_d$V1),
                   cross.area=length(eyadown$V1), category= c("Salmon Eya1 Downregulated","StringTie Eya1 Downregulated"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(200,160),cat.dist=rep(0.035,2),
                   cex=2, cat.cex= 1.2)


st_s_up= read.table("/home/stephen/sixupgenes.txt",stringsAsFactors = F)
sa_s_up= read.table("/home/stephen/Salmon_deseq2/Sal_sixupgenes.txt",stringsAsFactors = F)
sup= merge(x= st_s_up, y= sa_s_up, by= "V1")
write.table(sup, file="/home/stephen/SixUpcross.txt",row.names = F, quote = F)

grid.newpage()
draw.pairwise.venn(length(sa_s_up$V1),length(st_s_up$V1),
                   cross.area=length(sup$V1), category= c("Salmon Six1 Upregulated","StringTie Six1 Upregulated"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(200,160),cat.dist=rep(0.035,2),
                   cex=2, cat.cex= 1.5)

st_s_d= read.table("/home/stephen/sixdowngenes.txt",stringsAsFactors = F)
sa_s_d= read.table("/home/stephen/Salmon_deseq2/Sal_sixdowngenes.txt",stringsAsFactors = F)
sdown= merge(x= st_s_d, y= sa_s_d, by= "V1")
write.table(sdown, file="/home/stephen/SixDowncross.txt",row.names = F, quote = F)

grid.newpage()
draw.pairwise.venn(length(sa_s_d$V1),length(st_s_d$V1),
                   cross.area=length(sdown$V1), category= c("Salmon Six1 Downregulated","StringTie Six1 Downregulated"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(190,180),cat.dist=rep(0.045,2),
                   cex=2, cat.cex= 1.5)



st_se_up = read.table("/home/stephen/sixeyaupgenes.txt",stringsAsFactors = F)
sa_se_up = read.table("/home/stephen/Salmon_deseq2/Sal_sixeyaupgenes.txt",stringsAsFactors = F)
seup= merge(x= st_se_up, y= sa_se_up, by= "V1")
write.table(seup, file="/home/stephen/SixEyaUpcross.txt",row.names = F, quote = F)

grid.newpage()
draw.pairwise.venn(length(sa_se_up$V1),length(st_se_up$V1),
                   cross.area=length(seup$V1), category= c("Salmon Six1+Eya1 Upregulated","StringTie Six1+Eya1 Upregulatedregulated"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(200,170),cat.dist=rep(0.035,2),
                   cex=2, cat.cex= 1.4)



sa_se_d = read.table("/home/stephen/Salmon_deseq2/Sal_sixeyadowngenes.txt",stringsAsFactors = F)
st_se_d = read.table("/home/stephen/sixeyadowngenes.txt",stringsAsFactors = F)
sedown= merge(x= st_se_d, y= sa_se_d, by= "V1")
write.table(sedown, file="/home/stephen/SixEyaDowncross.txt",row.names = F, quote = F)

grid.newpage()
draw.pairwise.venn(length(sa_se_d$V1),length(st_se_d$V1),
                   cross.area=length(sedown$V1), category= c("Salmon Six1+Eya1 Downregulated","StringTie Six1+Eya Downregulated"),
                   lty= rep("blank",2),
                   fill=c("green","red"),alpha=rep(0.5,2),
                   cat.pos= c(190,160),cat.dist=rep(0.035,2),
                   cex=2, cat.cex= 1.4)
