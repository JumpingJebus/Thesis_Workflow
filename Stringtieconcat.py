#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 12:04:54 2017

@author: stephen
"""
import sys
import os
out= open(sys.argv[1], "a")

print("Reading in map.csv files." + "\n")
for i in range(2,(len(sys.argv))):
    print("file number " + str(i))
    Maps= open(sys.argv[i])
    Maps= Maps.readlines()
    print("Concatinating")
    for j in Maps:
        out.write(j)
print("New concatonated file created" +"\n")
out.close()
print("Loading in new concatonated file" + "\n")
infile = open(sys.argv[1],"r")

lines= infile.readlines()
infile.close()
mapping=[]

for t in lines:
    mapping.append(t)
print("Removing duplicates from concatonated file" + "\n")
mapping= list(set(mapping))

print("Writing unique GeneID:Gene_Name matches to Full Map file"+ "\n")
outfile= open(sys.argv[1],"w")
outfile.write("GeneID" + "\t" + "Gene_Name" + "\t" "Transcript_ID + "\n")
for r in mapping:
    outfile.write(r)
outfile.close()

print("Deleting old map files")
for z in range(2,(len(sys.argv))):
    os.remove(sys.argv[z])
    
print("Program finished")
