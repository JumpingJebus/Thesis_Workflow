library(data.table)
t2g = read.table(file.path("/home/stephen/Downloads/LAEVIS9.gtf"))
t2g_gene= t2g[,"V10"]
t2g_transcript = t2g[,"V13"]
df= data.frame(t2g_transcript,t2g_gene)
df= df[!duplicated(df),]
write.table(df, "/home/stephen/Downloads/t2gone1.csv", row.names= FALSE, sep=",")


t2g1 = read.table(file.path("/home/stephen/my.gtf"))
t2g1
t2g1_gene= t2g1[,"V16"]
t2g1_transcript = t2g1[,"V13"]
df1= data.frame(t2g1_transcript,t2g1_gene)
df1= df1[!duplicated(df1),]
write.table(df1, "/home/stephen/Downloads/tx2gene.csv", row.names= FALSE, sep=",")

head(t2g1,70)


t2g3 = read.table("/home/stephen/stringtieWUT.gtf", fill= TRUE, header= FALSE, quote="", encoding="UTF-8")
t2g3
mstrg= t2g3$V10
gene_name= t2g3$V18
length(t2g3$V1)
t2g1_gene= t2g1[,"V16"]
t2g1_transcript = t2g1[,"V13"]
df1= data.frame(t2g1_transcript,t2g1_gene)
df1= df1[!duplicated(df1),]
write.table(df1, "/home/stephen/Downloads/t2goneALT.csv", row.names= FALSE, sep=",")