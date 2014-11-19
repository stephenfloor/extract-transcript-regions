#!/usr/bin/env python

# Stephen N. Floor 
# 7 October 2014 
# floor@berkeley.edu 

# TODO: 

# convert UCSC gene names to refseq?
# add a region type which is mrna that contains the whole spliced transcript, and preserve cdsStart and cdsEnd for these.  - this is the same as codingExons, just need to preserve start/stop 

import sys, os
from UCSCKnownGene import * 

print " ---------------------------------"
print "| Extract Regions from knownGenes |"
print "|  snf   7 October 2014           |"
print " ---------------------------------\n\n"

if len(sys.argv) != 3:
    sys.exit("ERROR: Please provide the knownGenes file and output basename as inputs")


# output filenames: 
utr5FName = sys.argv[2] + "_5utr.bed"
utr5StartFName = sys.argv[2] + "_5utr_start.bed"
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

# parameters that should be passed via the cmd line
useBlocks = True 


#terminate if output files exist

if os.path.exists(utr5FName) or os.path.exists(utr5StartFName) or os.path.exists(cdsFName) or os.path.exists(utr3FName) or os.path.exists(exonFName) or os.path.exists(intronFName) \
        or os.path.exists(codingExonFName) or os.path.exists(codingIntronFName) or os.path.exists(noncodingExonFName) or os.path.exists(noncodingIntronFName):
    sys.exit("ERROR: output basename %s files already exist" % sys.argv[2]) 

#process the file

with open(sys.argv[1]) as knownGenesFile, open(utr5FName, "w") as utr5File, open(utr5StartFName, "w") as utr5StartFile, open(cdsFName, "w") as cdsFile, \
        open(utr3FName, "w") as utr3File, open(exonFName, "w") as exonFile, open (intronFName, "w") as intronFile, \
        open(codingExonFName, "w") as codingExonFile, open(codingIntronFName, "w") as codingIntronFile, \
        open(noncodingExonFName, "w") as noncodingExonFile, open(noncodingIntronFName, "w") as noncodingIntronFile:

    for line in knownGenesFile:
        # all of the knowngenes parsing and metadata construction is done inside UCSCKnownGene.py, especially the createGene method 
        gene = createGene(line)
        genesRead += 1

        if (useBlocks): # output all region primitives on the same line by specifying nBlocks and lists inside the BED output 
            if(gene.coding): 
                #blockBedFormat is one line by definition 
                if (gene.utr5Len > 0): utr5File.write(gene.blockBedFormat(region="5utr") + "\n")
                if (gene.utr5startLen > 0): utr5StartFile.write(gene.blockBedFormat(region="5utr_start") + "\n")
                if (gene.cdsLen > 0): cdsFile.write(gene.blockBedFormat(region="cds") + "\n")
                if (gene.utr3Len > 0): utr3File.write(gene.blockBedFormat(region="3utr") + "\n")
                
                if (gene.exonsLen > 0): 
                    exonFile.write(gene.blockBedFormat(region="exons") + "\n")
                    codingExonFile.write(gene.blockBedFormat(region="exons") + "\n")
                
                if (gene.intronsLen > 0):
                    intronFile.write(gene.blockBedFormat(region="introns") + "\n")
                    codingIntronFile.write(gene.blockBedFormat(region="introns") + "\n")
            
            else: # noncoding transcripts just have exons and introns 
                if (gene.exonsLen > 0):
                    exonFile.write(gene.blockBedFormat(region="exons") + "\n")
                    noncodingExonFile.write(gene.blockBedFormat(region="exons") + "\n")
                
                if (gene.intronsLen > 0):
                    intronFile.write(gene.blockBedFormat(region="introns") + "\n")
                    noncodingIntronFile.write(gene.blockBedFormat(region="introns") + "\n")

        else: # output one line per region primitive instead of combining regions via blocks 
            if(gene.coding): 
                for entry in gene.bedFormat(region="5utr"):
                    utr5File.write(entry + "\n")
                for entry in gene.bedFormat(region="5utr_start"):
                    utr5StartFile.write(entry + "\n")
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

        if (not genesRead % 2500):
            print "Processed %d entries..." %  genesRead

print "Processed %d entries." %  genesRead
