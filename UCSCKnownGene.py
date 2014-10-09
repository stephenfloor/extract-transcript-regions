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
        self.cds = []
        self.cdsLen = 0
        self.utr3 = []
        self.utr3Len = 0
        self.exonicRegions = []
        self.exonicRegionsLen = 0
        self.intronicRegions = []
        self.intronicRegionsLen = 0
    
        self.coding = False


    def __str__(self):  #currently roughly knownGenes format with a second line containing metadata 
        return "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d" % (self.name, self.chrom, self.strand, self.txStart, self.txEnd, self.cdsStart, self.cdsEnd, self.exonCt, self.exonStarts, self.exonEnds, self.utr5, self.utr5Len, self.cds, self.cdsLen, self.utr3, self.utr3Len, self.exonicRegions, self.exonicRegionsLen, self.intronicRegions, self.intronicRegionsLen)
    
#BED format output is goal.  Fields are optional after featureEnd 
# chrom    featureStart   featureEnd   nameOfLine   score(0-1000)   strand   thickStart  thickEnd  itemRGBtuple  blockCount  blockSizes   blockStarts 

#this function returns a list of BED-formatted strings for the feature passed as region 
    def bedFormat(self, region="exons"):
        if (not self.coding and (region == "5utr" or region == "cds" or region == "3utr")):
            print "UCSCKnownGene bedFormat error: noncoding transcripts do not have 5utr/cds/3utr"
            return []

        returnVal = []

        if (region == "5utr"):
            for chunk in self.utr5:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_5utr",chunk[1],chunk[2]))

        elif (region == "cds"):
            for chunk in self.cds:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_cds",chunk[1],chunk[2]))

        elif (region == "3utr"):
            for chunk in self.utr3:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_3utr",chunk[1],chunk[2]))

        elif (region == "exons"):
            for chunk in self.exonicRegions:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_exon",chunk[1],chunk[2]))

        elif (region == "introns"):
            for chunk in self.intronicRegions:
                returnVal.append("%s\t%d\t%d\t%s\t0\t%d\t%d\t0" % (chunk[0], chunk[1], chunk[2], chunk[3]+"_intron",chunk[1],chunk[2]))
        else:
            print "UCSCKnownGene bedFormat error: currently only regions 5utr/cds/3utr/exons/introns are supported"
            

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
                    foo.cds.append((foo.chrom, foo.cdsStart, foo.cdsEnd, foo.name))
                    foo.cdsLen += foo.cdsEnd - foo.cdsStart
                    foo.utr3.append((foo.chrom, foo.cdsEnd, foo.exonEnds[i], foo.name))
                    foo.utr3Len += foo.exonEnds[i] - foo.cdsEnd
                #case 2 - exon spans 5' UTR/CDS junction
                elif (foo.exonStarts[i] < foo.cdsStart and foo.exonEnds[i] >= foo.cdsStart):
                    foo.utr5.append((foo.chrom, foo.exonStarts[i], foo.cdsStart, foo.name))
                    foo.utr5Len += foo.cdsStart - foo.exonStarts[i]
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
            foo.exonicRegions.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
            foo.exonicRegionsLen += foo.exonEnds[i] - foo.exonStarts[i]
                
            #print "DBUG2: i %d foo.exonCt-1 %d foo.exonEnds %s foo.exonStarts %s" % (i, foo.exonCt-1, foo.exonEnds, foo.exonStarts)
        
            if (i < foo.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                foo.intronicRegions.append((foo.chrom, foo.exonEnds[i], foo.exonStarts[i+1], foo.name))
                foo.intronicRegionsLen += foo.exonStarts[i+1] - foo.exonEnds[i] 

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
                else: 
                    print "Thar be dragons - UCSCKnownGene createGene - stranded gene region parsing" 
                    
            #else: 
            #    print "- strand noncoding transcript"
                
            # -- generate combined exonic and intronic regions -- 
            #exons are easy 
            foo.exonicRegions.append((foo.chrom, foo.exonStarts[i], foo.exonEnds[i], foo.name))
            foo.exonicRegionsLen += foo.exonEnds[i] - foo.exonStarts[i]
            
            if (i < foo.exonCt - 1): # only compute introns for nonterminal exons
                # an intron is the region between the end of the current exon and start of the next 
                foo.intronicRegions.append((foo.chrom, foo.exonEnds[i], foo.exonStarts[i+1], foo.name))
                foo.intronicRegionsLen += foo.exonStarts[i+1] - foo.exonEnds[i] 
                
    else:
        print "Thar be dragons - UCSCKnownGene createGene strand does not match + or -"
	
    return foo 

