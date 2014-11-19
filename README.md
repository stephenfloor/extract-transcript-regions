knowngenes-to-transcript-regions
================================

For biological deep sequencing data.  Decompose a UCSC knownGenes file into transcript regions (i.e. exons, introns, UTRs and CDS).

This program takes a knownGene.txt file for some genome from the UCSC genome browser and decomposes it into the following transcript regions: 

  - exons
  - introns
  - exons from coding transcripts
  - introns from coding transcripts
  - exons from noncoding transcripts
  - introns from noncoding transcripts
  - 5' UTRs for coding transcripts
  - 5' UTRs plus start codon plus +4 nucleotide for coding transcripts
  - CDS for coding transcripts
  - 3' UTRs for coding transcripts

Invocation:

  python knowngenes_to_transcript_regions.py knownGenes.txt output_basename

Input: 

  UCSC knownGene.txt
  
Output: 

  ten region files in .bed format 
  
This is built on the UCSCKnownGene class, which constructs gene objects from a knownGene file and computes metadata about the genes during parsing of the knownGene file. 

Notes: 

  - output is now in blockBedFormat by default (i.e. all exons dumped onto one line using block columns at the end) 
  - to change output format, change the value of useBlocks in knowngenes_to_transcript_regions.py 

