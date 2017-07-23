library(dplyr)
#read in the concatonated map file
data =as.data.frame(read.csv("/home/stephen/Downloads/UNEWER/two/map.csv", sep= ","))
#filter out any non novel genes
data = data %>% group_by(Gene_Name) %>% filter(!("Novel" %in% Gene_Name))
#write the unique mapping file
write.csv(data, file= "/home/stephen/MapFileNoMSTRG.csv", row.names= FALSE)
data= as.data.frame(read.csv("/home/stephen/MapFileNoMSTRG.csv", stringsAsFactors= FALSE))

#SIX1
#read in the stringtie deseq2 results
six= as.data.frame(read.csv("/home/stephen/Downloads/UNEWER/two/six.csv" , stringsAsFactors= FALSE))
#replace any row which matches the map file with the corresponding gene name in the map file
for(id in 1:nrow(six)){
  six$X[six$X %in% data$GeneID[id]] <- data$Gene_Name[id]
}
#write the new deseq2 results file
write.csv(six, file=("/home/stephen/StringTieSix.csv"), row.names = FALSE)