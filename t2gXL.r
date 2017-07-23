
t2g1 = read.table(file.path("/home/stephen/my.gtf"))
t2g1_gene= t2g1[,"V16"]
t2g1_transcript = t2g1[,"V13"]
df1= data.frame(t2g1_transcript,t2g1_gene)
df1= df1[!duplicated(df1),]
write.table(df1, "/home/stephen/Downloads/tx2gene.csv", row.names= FALSE, sep=",")
