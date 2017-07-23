library(dplyr)

map_file= read.table("/home/stephen/GenePageToJgiLaevisScaffoldMapping_9.1.txt")
Transcript_ID= map_file$V3
GeneID= map_file$V2
t2g= data.frame(Transcript_ID, GeneID)
t2g= t2g[!duplicated(t2g),]
write.table(t2g, "/home/stephen/Downloads/UNEWER/two/tx2NEW.csv", row.names= FALSE, sep=",")
t2g= as.data.frame(read.csv(file="/home/stephen/Downloads/UNEWER/two/tx2NEW.csv", sep= ",",stringsAsFactors= FALSE))



#SIX1
six=as.data.frame(read.csv("/home/stephen/Salmon_deseq2/SigSix.csv" , stringsAsFactors= FALSE))

#adding gene name from new t2g mapping file
for(id in 1:nrow(t2g)){
  six$X[six$X %in% t2g$Transcript_ID[id]] <- t2g$GeneID[id]
}

sixdf = as.data.frame(six)
sixup= sixdf[sixdf$log2FoldChange >0,]
sixupgenes= sixup$X

write.table(sixupgenes, file= "/home/stephen/Salmon_deseq2/Sal_sixupgenes.txt", row.names= F,quote = F)
sixdown= sixdf[sixdf$log2FoldChange <0,]
sixdowngenes= sixdown$X

write.table(sixdowngenes, file= "/home/stephen/Salmon_deseq2/Sal_sixdowngenes.txt",row.names=F,quote = F)
write.table(sixdf, file= "/home/stephen/Salmon_deseq2/Sal_sixgenes.txt", row.names = F)

#Write fixed results to a file
write.csv(six, file=("/home/stephen/Salmon_deseq2/Sal_six_fix.csv"), row.names = FALSE)

#EYA1
eya= as.data.frame(read.csv("/home/stephen/Salmon_deseq2/SigEya.csv" , stringsAsFactors= FALSE))


for(id in 1:nrow(t2g)){
  eya$X[eya$X %in% t2g$Transcript_ID[id]] <- t2g$GeneID[id]
}
write.csv(eya, file=("/home/stephen/Salmon_deseq2/Sal_eya_fix.csv"), row.names = FALSE)

eyadf= as.data.frame(eya)
eyaup= eyadf[eyadf$log2FoldChange >0,]
eyaupgenes= eyaup$X
write.table(eyaupgenes, file= "/home/stephen/Salmon_deseq2/Sal_eyaupgenes.txt", row.names= F, quote = F)

eyadown= eyadf[eyadf$log2FoldChange <0,]
eyadowngenes= eyadown$X

write.table(eyadowngenes, file= "/home/stephen/Salmon_deseq2/Sal_eyadowngenes.txt",row.names=F, quote = F)


#SIXEYA
sixeya= as.data.frame(read.csv("/home/stephen/Salmon_deseq2/SigSixEya.csv" , stringsAsFactors= FALSE))

for(id in 1:nrow(t2g)){
  sixeya$X[sixeya$X %in% t2g$Transcript_ID[id]] <- t2g$GeneID[id]
}
head(t2g)
write.csv(sixeya, file=("/home/stephen/Salmon_deseq2/Sal_fix_SixEya.csv"), row.names = FALSE)

sixeyadf= as.data.frame(sixeya)
sixeyaup= sixeyadf[sixeyadf$log2FoldChange >0,]
sixeyaupgenes= sixeyaup$X
write.table(sixeyaupgenes, file= "/home/stephen/Salmon_deseq2/Sal_sixeyaupgenes.txt", row.names= F, quote = F)
sixeyadown= sixeyadf[sixeyadf$log2FoldChange <0,]
sixeyadowngenes= sixeyadown$X
write.table(sixeyadowngenes, file= "/home/stephen/Salmon_deseq2/Sal_sixeyadowngenes.txt",row.names=F, quote= F)