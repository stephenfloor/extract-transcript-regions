# class to create genes out of lines in knowngene file 

# Stephen N. Floor
# 8 October 2014 

class UCSCKnownGene: 
    def __init__(self):
        #properties defined in UCSC knowngenes 
        self.name = '' 
        self.chrom = ''
        self.strand = '' 
        self.txStart = 0
        self.txEnd = 0 
        self.cdsStart = 0
        self.cdsEnd = 0
        self.exonCt = 0
        self.exonStarts = []
        self.exonEnds = []
        self.exonLengths = []
        
        #meta properties to be computed during construction.  these are lists of BED first four field tuples with the exception of Len terms which are the length of the total region for the gene 
        self.utr5 = []
        self.utr5Len = 0
        self.utr5start = []
        self.utr5startLen = 0
        self.cds = []
        self.cdsLen = 0
        self.utr3 = []
        self.utr3Len = 0
        self.exons = []
        self.exonsLen = 0
        self.introns = []
        self.intronsLen = 0
    
        self.coding = False


    def __str__(self):  #currently roughly knownGenes format with a second line containing metadata 
        return "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d" % (self.name, self.chrom, self.strand, self.txStart, self.txEnd, self.cdsStart, self.cdsEnd, self.exonCt, self.exonStarts, self.exonEnds, self.utr5, self.utr5Len, self.cds, self.cdsLen, self.utr3, self.utr3Len, self.exons, self.exonsLen, self.introns, self.intronsLen)
    
#BED format output is goal.  Fields are optional after featureEnd 
# chrom    featureStart   featureEnd   nameOfLine   score(0-1000)   strand   thickStart  thickEnd  itemRGBtuple  blockCount  blockSizes   blockStarts 

#this function returns a list of BED-formatted strings for the feature passed as region with multiple entries per region possible, one for each primitive (exon/intron) 
    def bedFormat(self, region="exons"):
        if (not self.coding and (region == "5utr" or region == "cds" or region == "3utr")):
            print "UCSCKnownGene bedFormat error: noncoding transcripts do not have 5utr/cds/3utr"
            return []

        returnVal = []

        if (region == "5utr"):
            for chunk in self.utr5:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_5utr",self.strand, chunk[1],chunk[2]))

        elif (region == "5utr_start"):
            for chunk in self.utr5_start:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_5utr_start",self.strand, chunk[1],chunk[2]))

        elif (region == "cds"):
            for chunk in self.cds:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_cds",self.strand, chunk[1],chunk[2]))

        elif (region == "3utr"):
            for chunk in self.utr3:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_3utr",self.strand, chunk[1],chunk[2]))

        elif (region == "exons"):
            for chunk in self.exons:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_exon",self.strand, chunk[1],chunk[2]))

        elif (region == "introns"):
            for chunk in self.introns:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%c\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_intron",self.strand, chunk[1],chunk[2]))
        else:
            print "UCSCKnownGene bedFormat error: currently only regions 5utr/cds/3utr/exons/introns are supported"
            

        return returnVal

#BED format output is goal.  Fields are optional after featureEnd 
# chrom    featureStart   featureEnd   nameOfLine   score(0-1000)   strand   thickStart  thickEnd  itemRGBtuple  blockCount  blockSizes   blockStarts 

#blockCount - The number of blocks (exons) in the BED line.
#blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
#blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount. 

#this function returns a BED-formatted string for the feature passed as region with blocks defining the exons as per the BED file format 
    def blockBedFormat(self, region="exons"):
        if (not self.coding and (region == "5utr" or region == "cds" or region == "3utr")):
            print "UCSCKnownGene blockBedFormat error: noncoding transcripts do not have 5utr/cds/3utr"
            return ""

        returnVal = ""
        score = 0 
        rgb = 0

        if (region == "5utr"):
            chromStart = self.utr5[0][1] # start of feature is start of first block
            chromEnd = self.utr5[-1][2] # end of feature is end of last block 
            regionName = self.name + "_5utr"

            blockCount = len(self.utr5)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.utr5])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.utr5])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "5utr_start"):
            chromStart = self.utr5start[0][1] # start of feature is start of first block
            chromEnd = self.utr5start[-1][2] # end of feature is end of last block 
            regionName = self.name + "_5utr_start"

            blockCount = len(self.utr5start)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.utr5start])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.utr5start])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)
            
        elif (region == "cds"):
            chromStart = self.cds[0][1] # start of feature is start of first block
            chromEnd = self.cds[-1][2] # end of feature is end of last block 
            regionName = self.name + "_cds"
            blockCount = len(self.cds)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.cds])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.cds])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "3utr"):
            chromStart = self.utr3[0][1] # start of feature is start of first block
            chromEnd = self.utr3[-1][2] # end of feature is end of last block 
            regionName = self.name + "_3utr"
            blockCount = len(self.utr3)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.utr3])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.utr3])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "exons"):
            chromStart = self.exons[0][1] # start of feature is start of first block
            chromEnd = self.exons[-1][2] # end of feature is end of last block 
            regionName = self.name + "_exon"
            blockCount = len(self.exons)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.exons])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.exons])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)

        elif (region == "introns"):
            chromStart = self.introns[0][1] # start of feature is start of first block
            chromEnd = self.introns[-1][2] # end of feature is end of last block 
            regionName = self.name + "_intron"
            blockCount = len(self.introns)
            blockSizes = ''.join(["%d," % (chunk[2]-chunk[1]) for chunk in self.introns])
            blockStarts = ''.join(["%d," % (chunk[1]-chromStart) for chunk in self.introns])

            #print "blockCount %d blockSizes %s blockStarts %s" % (blockCount, blockSizes,  blockStarts)
        else:
            print "UCSCKnownGene blockBedFormat error: currently only regions 5utr/cds/3utr/exons/introns are supported"
            
        returnVal = "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%s\t%d\t%s\t%s" % (self.chrom, chromStart, chromEnd, regionName, score, self.strand, chromStart, chromEnd, rgb, blockCount, blockSizes, blockStarts)

        return returnVal
        
            

#example line format from knownGenes file (from UCSC) 
    #    uc010nxq.1      chr1    +       11873   14409   12189   13639   3       11873,12594,13402,      12227,12721,14409,      B7ZGX9  uc010nxq.1
    # line format
    #    name            chrom   strand  txStart txEnd   cdsStart cdsEnd exonCt  exonStarts              exonEnds                proteinID  alignID 

def createGene(knownGeneLineString):
    foo = UCSCKnownGene()
    line = knownGeneLineString.split()

    # -- read in knownGene fields -- 

    foo.name = line[0]
    foo.chrom = line[1]
    foo.strand = line[2]
    foo.txStart = int(line[3])
    foo.txEnd = int(line[4])
    foo.cdsStart = int(line[5])
    foo.cdsEnd = int(line[6])
    foo.exonCt = int(line[7])
    
    starts = line[8].split(",")
    ends = line[9].split(",")

    for i in range(foo.exonCt):
        foo.exonStarts.append(int(starts[i]))
        foo.exonEnds.append(int(ends[i]))

    # not doing it this way because of the trailing , in the knowngenes format 
    #for i in line[8].split(","): foo.exonStarts.append(i)
    #for i in line[9].split(","): foo.exonEnds.append(i)

    # -- begin computing metadata -- 

    # -- note: chose clarity of code and conditionals here over most efficient computation (i.e. some clauses may be redundant)

    if (foo.strand == "+"): 
        #print ("DBUG - exonCt %d i %d exonEnds[i] %d cdsStart %d exonStarts[i] %d cdsEnd %d") % \
        #    (foo.exonCt, i, foo.exonEnds[i], foo.cdsStart, foo.exonStarts[i], foo.cdsEnd)
        for i in range (foo.exonCt): 
            if (foo.cdsStart != foo.cdsEnd): # if this is a coding transcript
                foo.coding = True
                # -- first compute 5'utr, CDS, 3'utr regions --
                #case 1 - exon spans 5' UTR/CDS/3' UTR
                if (foo.exonStarts[i] < foo.cdsStart and foo.exonEnds[i] > foo.cdsEnd):
                    foo.utr5.append((foo.chrom, foo.exonStarts[i], foo.cdsStart, foo.name))
                    foo.utr5Len += foo.cdsStart - foo.exonStarts[i]
                    foo.utr5start.append((foo.chrom, foo.exonStarts[i], foo.cdsStart + 4, foo.name)) # add four here to include the ATG and the +4 nucleotide 
                    foo.utr5startLen += foo.cdsStart + 4 - foo.exonStarts[i]
                    foo.cds.append((foo.chrom, foo.cdsStart, foo.cdsEnd, foo.name))
                    foo.cdsLen += foo.cdsEnd - foo.cdsStart
                    foo.utr3.append((foo.chrom, foo.cdsEnd, foo.exonEnds[i], foo.name))
                    foo.utr3Len += foo.exonEnds[i] - foo.cdsEnd
                #case 2 - exon spans 5' UTR/CDS junction
                elif (foo.exonStarts[i] < foo.cdsStart and foo.exonEnds[i] >= foo.cdsStart):
                    foo.utr5.append((foo.chrom, foo.exonStarts[i], foo.cdsStart, foo.name))
                    foo.utr5Len += foo.cdsStart - foo.exonStarts[i]
                    foo.utr5start.append((foo.chrom, foo.exonStarts[i], foo.cdsStart + 4, foo.name)) # add four here to include the ATG and the +4 nucleotide 
                    foo.utr5startLen += foo.cdsStart + 4 - foo.exonStarts[i]
                    foo.cds.append((foo.chrom, foo.cdsStart, foo.exonEnds[i], foo.name))
                    foo.cdsLen += foo.exonEnds[i]- foo.cdsStart
                #case 3 - exon spans CDS/3'UTR junction 
                elif (foo.exonStarts[i] >= foo.cdsStart and foo.exonStarts[i] <= foo.cdsEnd and foo.exonEnds[i] > foo.cdsEnd):
                    foo.cds.append((foo.chrom, foo.exonStarts[i], foo.cdsEnd, foo.name))
                    foo.cdsLen += foo.cdsEnd - foo.exonStarts[i]
                    foo.utr3.append((foo.chrom, foo.cdsEnd, foo.exonEnds[i], foo.name))
                    foo.utr3Len += foo.exonEnds[i] - foo.cdsEnd
                #case 4 - exon is 5' UTR only 
                elif (foo.exonStarts[i] < foo.cdsStart and foo.exonEnds[i] < foo.cdsStart): 
                    foo.utr5.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
                    foo.utr5Len += foo.exonEnds[i] - foo.exonStarts[i]
                    foo.utr5start.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name)) # add four here to include the ATG and the +4 nucleotide 
                    foo.utr5startLen += foo.exonEnds[i] - foo.exonStarts[i]
                #case 5 - exon is CDS only
                elif (foo.exonStarts[i] >= foo.cdsStart and foo.exonEnds[i] <= foo.cdsEnd):
                    foo.cds.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
                    foo.cdsLen += foo.exonEnds[i] - foo.exonStarts[i]
                #case 6 - exon is 3' UTR only 
                elif (foo.exonStarts[i] > foo.cdsEnd and foo.exonEnds[i] > foo.cdsEnd):
                    foo.utr3.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
                    foo.utr3Len += foo.exonEnds[i] - foo.exonStarts[i]
                else: 
                    print "Thar be dragons - UCSCKnownGene createGene + stranded gene region parsing" 
            #else: 
            #    print "noncoding + strand transcript" 
                
            # -- generate combined exonic and intronic regions -- 
            #exons are easy 
            foo.exons.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
            foo.exonsLen += foo.exonEnds[i] - foo.exonStarts[i]
                
            #print "DBUG2: i %d foo.exonCt-1 %d foo.exonEnds %s foo.exonStarts %s" % (i, foo.exonCt-1, foo.exonEnds, foo.exonStarts)
        
            if (i < foo.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                foo.introns.append((foo.chrom, foo.exonEnds[i], foo.exonStarts[i+1], foo.name))
                foo.intronsLen += foo.exonStarts[i+1] - foo.exonEnds[i] 

    elif (foo.strand == "-"):
     #uc001ach.2	    chr1    -	    910578  917473  911551  916546  5	    910578,911878,914260,916516,917444,	    911649,912004,916037,916553,917473,	    Q5SV97  uc001ach.2
     #	name		chrom	strand	txStart txEnd	cdsStart foo.cdsEnd exonCt	exonStarts		exonEnds		proteinID  alignID 
     # for the minus strand everything is the same except the order of encountering regions is reversed
     # i.e. 3' UTR -> CDS -> 5' UTR 

        for i in range (foo.exonCt): 
            #print ("DBUG - exonCt %d i %d foo.exonEnds[i] %d foo.cdsStart %d exonStarts[i] %d foo.cdsEnd %d") % \
            #    (foo.exonCt, i, foo.exonEnds[i], foo.cdsStart, foo.exonStarts[i], foo.cdsEnd)
            
            if (foo.cdsStart != foo.cdsEnd):
                foo.coding = True 
                # -- first compute 5'utr, CDS, 3'utr regions --
                # -- this is the same as for + sense except 5' UTR and 3' UTR are swapped throughout
                #case 1 - exon spans 3' UTR/CDS/5' UTR
                if (foo.exonStarts[i] < foo.cdsStart and foo.exonEnds[i] > foo.cdsEnd):
                    foo.utr3.append((foo.chrom, foo.exonStarts[i], foo.cdsStart, foo.name))
                    foo.utr3Len += foo.cdsStart - foo.exonStarts[i]
                    foo.cds.append((foo.chrom, foo.cdsStart, foo.cdsEnd, foo.name))
                    foo.cdsLen += foo.cdsEnd - foo.cdsStart
                    foo.utr5.append((foo.chrom, foo.cdsEnd, foo.exonEnds[i], foo.name))
                    foo.utr5Len += foo.exonEnds[i] - foo.cdsEnd
                    foo.utr5start.append((foo.chrom, foo.cdsEnd - 4 , foo.exonEnds[i], foo.name))
                    foo.utr5startLen += foo.exonEnds[i] - (foo.cdsEnd - 4)
                #case 2 - exon spans 3' UTR/CDS junction
                elif (foo.exonStarts[i] < foo.cdsStart and foo.exonEnds[i] >= foo.cdsStart):
                    foo.utr3.append((foo.chrom, foo.exonStarts[i], foo.cdsStart, foo.name))
                    foo.utr3Len += foo.cdsStart - foo.exonStarts[i]
                    foo.cds.append((foo.chrom, foo.cdsStart, foo.exonEnds[i], foo.name))
                    foo.cdsLen += foo.exonEnds[i]- foo.cdsStart
                #case 3 - exon spans CDS/5'UTR junction 
                elif (foo.exonStarts[i] >= foo.cdsStart and foo.exonStarts[i] <= foo.cdsEnd and foo.exonEnds[i] > foo.cdsEnd):
                    foo.cds.append((foo.chrom, foo.exonStarts[i], foo.cdsEnd, foo.name))
                    foo.cdsLen += foo.cdsEnd - foo.exonStarts[i]
                    foo.utr5.append((foo.chrom, foo.cdsEnd, foo.exonEnds[i], foo.name))
                    foo.utr5Len += foo.exonEnds[i] - foo.cdsEnd
                    foo.utr5start.append((foo.chrom, foo.cdsEnd - 4 , foo.exonEnds[i], foo.name))
                    foo.utr5startLen += foo.exonEnds[i] - (foo.cdsEnd - 4)
                #case 4 - exon is 3' UTR only 
                elif (foo.exonStarts[i] < foo.cdsStart and foo.exonEnds[i] < foo.cdsStart): 
                    foo.utr3.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
                    foo.utr3Len += foo.exonEnds[i] - foo.exonStarts[i]
                #case 5 - exon is CDS only
                elif (foo.exonStarts[i] >= foo.cdsStart and foo.exonEnds[i] <= foo.cdsEnd):
                    foo.cds.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
                    foo.cdsLen += foo.exonEnds[i] - foo.exonStarts[i]
                #case 6 - exon is 5' UTR only 
                elif (foo.exonStarts[i] > foo.cdsEnd and foo.exonEnds[i] > foo.cdsEnd):
                    foo.utr5.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
                    foo.utr5Len += foo.exonEnds[i] - foo.exonStarts[i]
                    foo.utr5start.append((foo.chrom, foo.exonStarts[i] , foo.exonEnds[i], foo.name))
                    foo.utr5startLen += foo.exonEnds[i] - foo.exonStarts[i]
                else: 
                    print "Thar be dragons - UCSCKnownGene createGene - stranded gene region parsing" 
                    
            #else: 
            #    print "- strand noncoding transcript"
                
            # -- generate combined exonic and intronic regions -- 
            #exons are easy 
            foo.exons.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
            foo.exonsLen += foo.exonEnds[i] - foo.exonStarts[i]
            
            if (i < foo.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                foo.introns.append((foo.chrom, foo.exonEnds[i], foo.exonStarts[i+1], foo.name))
                foo.intronsLen += foo.exonStarts[i+1] - foo.exonEnds[i] 
                
    else:
        print "Thar be dragons - UCSCKnownGene createGene strand does not match + or -"
	
    return foo 

