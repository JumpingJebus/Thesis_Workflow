#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:07:40 2017

@author: stephen
"""
import re
import sys
print("Loading in the GTF file")
gtf = open(sys.argv[1], "r")
#gtf= open("/home/stephen/Downloads/two/un/SRR3374051/testing2.gtf")
gtf = gtf.read()

#this will vary depending on the command used to run stringtie. 
StringtieCommand= '# stringtie /data4/Stephen/xl_genome/tophat_output/SRR3374051/SRR3374051.bam -o /data4/Stephen/xl_genome/stringtieALT3/two/SRR3374051/SRR3374051.gtf -G /data4/Stephen/xl_genome/stringtieALT3/stringtieMerged2.gtf -A /data4/Stephen/xl_genome/stringtieALT3/two/SRR3374051/SRR3374051_gene_abund.tab -B -e -p 6 -a 30 -j 50 -c 20 -g 20# StringTie version 1.3.3b\n'

gtf = gtf[(len(StringtieCommand)):]

gtf = gtf.split('\n')
print("Getting GeneID,GeneName and Transcript_ID information")
code=[]
for i in gtf:
    Name = re.findall(r'ref_gene_name "(.*?)";',i)
    ID= re.findall(r'gene_id "(.*?)";',i)
    
    if len(ID) == 0:
        code.append("Error")
    else:
        code.append(str(ID[0]))
      
    if len(Name)==0:
        if "X" in str(ID) or "L" in str(ID) or "l" in str(ID):
            code.append(ID[0])
        else:
            code.append("Novel")
    else:
        code.append(str(Name[0]))
        

print("Matching GeneID to GeneName and Transcript_ID")
mapping= []
for i in range(0,(int(len(code))),2):
    mapping.append(code[i] + "\t" + code[i+1] + ("\n"))

print("Removing duplicates")
mapping = list(set(mapping))
  
print("Writing mapping information to new file")
out= open(sys.argv[2],"w")
#out= open("/home/stephen/Downloads/two/un/SRR3374051/maps.csv", "w")
out.write("GeneID" + "\t" + "Gene_name" + "\n")
for i in mapping:
    out.write(i)

print("Program finished")


