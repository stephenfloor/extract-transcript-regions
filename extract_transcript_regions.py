#!/usr/bin/env python

# Stephen N. Floor
# 7 October 2014
# floor@berkeley.edu

import sys, os, argparse

import GTF

from collections import defaultdict 

from Transcript import *

print " ----------------------------------"
print "| Extract Regions from annotations |"
print "|  snf        Fall 2014            |"
print " ----------------------------------\n\n"


# ------ ARGUMENT PARSING ----------

parser = argparse.ArgumentParser(description="Create transcript regions (5' UTR/CDS/3'UTR etc) from knownGenes or a GTF") 

parser.add_argument("-i", "--input", help="input filename", required=True)
parser.add_argument("-o", "--output", help="output basename", required=True) 
parser.add_argument("--ucsc", help="Read from a UCSC knownGenes formatted file (BED)", action="store_true")
parser.add_argument("--gtf", help="Read from a GTF (only tested with Ensembl GTFs)", action="store_true")

args = parser.parse_args() 

if ( (not (args.ucsc or args.gtf)) or (args.ucsc and args.gtf)):
    sys.exit("FATAL: must set one but not both of --ucsc and --gtf") 

# output filenames:
utr5FName = args.output + "_5utr.bed"
utr5StartFName = args.output + "_5utr_start.bed"
cdsFName = args.output + "_cds.bed"
utr3FName = args.output + "_3utr.bed"
exonFName = args.output + "_exons.bed"
intronFName = args.output + "_introns.bed"
codingExonFName = args.output + "_codingexons.bed"
codingIntronFName = args.output + "_codingintrons.bed" # note that these are introns from coding genes, not necessarily introns that make it to mRNA
noncodingExonFName = args.output + "_noncodingexons.bed"
noncodingIntronFName = args.output + "_noncodingintrons.bed"

#keep track of where we are
genesRead = 0

# parameters that should be passed via the cmd line
useBlocks = True


#terminate if output files exist

if os.path.exists(utr5FName) or os.path.exists(utr5StartFName) or os.path.exists(cdsFName) or os.path.exists(utr3FName) or os.path.exists(exonFName) or os.path.exists(intronFName) \
        or os.path.exists(codingExonFName) or os.path.exists(codingIntronFName) or os.path.exists(noncodingExonFName) or os.path.exists(noncodingIntronFName):
    sys.exit("ERROR: output basename %s files already exist" % args.output)

#process the file

with open(utr5FName, "w") as utr5File, open(utr5StartFName, "w") as utr5StartFile, open(cdsFName, "w") as cdsFile, \
        open(utr3FName, "w") as utr3File, open(exonFName, "w") as exonFile, open (intronFName, "w") as intronFile, \
        open(codingExonFName, "w") as codingExonFile, open(codingIntronFName, "w") as codingIntronFile, \
        open(noncodingExonFName, "w") as noncodingExonFile, open(noncodingIntronFName, "w") as noncodingIntronFile:

    def writeOutput(gene):
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


    if (args.ucsc): 
        with open(args.input, "r") as genesFile: 
        
            for line in genesFile:
                # all of the knowngenes parsing and metadata construction is done inside Transcript.py, especially the createGene method

                gene = createUCSCTranscript(line) 
                genesRead += 1

                writeOutput(gene)

                if (not genesRead % 2500):
                    print "Processed %d entries..." %  genesRead

                
    elif (args.gtf): 
            
            # first parse the entire file into a dictionary of lists

        txDict = defaultdict(list) 

        print "Building GTF dictionary..." 

        # the issue here is that lines for various transcripts may be interleaved, so can either create lots of objects, or a giant dict. opted for giant dict. 
        for line in GTF.lines(args.input): 
            # only want to read in lines corresponding to these features
            if line["feature"] in ["exon", "CDS", "start_codon", "stop_codon"]:
                txDict[line["transcript_id"]].append(line)
                genesRead += 1

                if (not genesRead % 25000):
                    print "\tProcessed %d lines..." %  genesRead

        print "Dictionary built." 

        print "Writing transcript properties."
        genesRead = 0
        
        # now create a Transcript object for each transcript and output it 

        for key in txDict: 

            #print key

            tx = createGTFTranscript(txDict[key])

            #print tx 
            writeOutput(tx)
            genesRead += 1
            
            if (not genesRead % 2500):
                print "\tProcessed %d entries..." %  genesRead


print "Processed %d entries." %  genesRead


