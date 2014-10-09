#!/usr/bin/env python

# Stephen N. Floor 
# 7 October 2014 
# floor@berkeley.edu 

# TODO: 

# convert UCSC gene names to refseq?

import sys, os
from snfGene import * 

print " ---------------------------------"
print "| Extract Regions from knownGenes |"
print "|  snf   7 October 2014           |"
print " ---------------------------------\n\n"

if len(sys.argv) != 3:
    sys.exit("ERROR: Please provide the knownGenes file and output basename as inputs")


# output filenames: 
utr5FName = sys.argv[2] + "_5utr.bed"
cdsFName = sys.argv[2] + "_cds.bed"
utr3FName = sys.argv[2] + "_3utr.bed"
exonFName = sys.argv[2] + "_exons.bed"
intronFName = sys.argv[2] + "_introns.bed"
codingExonFName = sys.argv[2] + "_codingexons.bed"
codingIntronFName = sys.argv[2] + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA 
noncodingExonFName = sys.argv[2] + "_noncodingexons.bed" 
noncodingIntronFName = sys.argv[2] + "_noncodingintrons.bed" 

#keep track of where we are 
genesRead = 0

#terminate if output files exist

if os.path.exists(utr5FName) or os.path.exists(cdsFName) or os.path.exists(utr3FName) or os.path.exists(exonFName) or os.path.exists(intronFName) \
        or os.path.exists(codingExonFName) or os.path.exists(codingIntronFName) or os.path.exists(noncodingExonFName) or os.path.exists(noncodingIntronFName):
    sys.exit("ERROR: output basename %s files already exist" % sys.argv[2]) 

#process the file

with open(sys.argv[1]) as knownGenesFile, open(utr5FName, "w") as utr5File, open(cdsFName, "w") as cdsFile, \
        open(utr3FName, "w") as utr3File, open(exonFName, "w") as exonFile, open (intronFName, "w") as intronFile, \
        open(codingExonFName, "w") as codingExonFile, open(codingIntronFName, "w") as codingIntronFile, \
        open(noncodingExonFName, "w") as noncodingExonFile, open(noncodingIntronFName, "w") as noncodingIntronFile:

    for line in knownGenesFile:
        # all of the knowngenes parsing and metadata construction is done inside snfGene.py, especially the createGene method 
        gene = createGene(line)
        genesRead += 1

        if(gene.coding): 
            for entry in gene.bedFormat(region="5utr"):
                utr5File.write(entry + "\n")
            for entry in gene.bedFormat(region="cds"):
                cdsFile.write(entry + "\n")
            for entry in gene.bedFormat(region="3utr"):
                utr3File.write(entry + "\n")

            for entry in gene.bedFormat(region="exons"):
                exonFile.write(entry + "\n")
                codingExonFile.write(entry + "\n")
            
            for entry in gene.bedFormat(region="introns"):
                intronFile.write(entry + "\n")
                codingIntronFile.write(entry + "\n")
            
        else: # noncoding transcripts just have exons and introns 
            for entry in gene.bedFormat(region="exons"):
                exonFile.write(entry + "\n")
                noncodingExonFile.write(entry + "\n")
            
            for entry in gene.bedFormat(region="introns"):
                intronFile.write(entry + "\n")
                noncodingIntronFile.write(entry + "\n")

        if (not genesRead % 5000):
            print "Processed %d genes..." %  genesRead

print "Processed %d genes." %  genesRead
