extract-transcript-regions
================================

For bioinformatics on genome annotation sets. Decompose a UCSC knownGenes file or Ensembl-derived GTF into transcript regions (i.e. exons, introns, UTRs and CDS).

This program takes either a knownGene.txt file for some genome from the UCSC genome browser or a GTF for transcripts from Ensembl and decomposes it into the following transcript regions: 

  - exons
  - introns
  - exons from coding transcripts
  - introns from coding transcripts
  - exons from noncoding transcripts
  - introns from noncoding transcripts
  - 5' UTRs for coding transcripts
  - 5' UTRs plus start codon plus 27 nucleotides for coding transcripts (to calculate Kozak context & uORF overlap with start codon) 
  - CDS for coding transcripts
  - 3' UTRs for coding transcripts

Invocation:

  python extract_transcript_regions.py -i knownGenes.txt -o output_basename (--gtf || --ucsc)

For help enter extract_transcript_regions.py -h 

Input: 

  UCSC knownGene.txt or Ensembl GTF 
  
Output: 

  ten region files in .bed format 
  
This is built on the Transcript class, which constructs transcript objects from the input file and computes metadata about the transcripts. 

Notes: 

  - output is now in blockBedFormat by default (i.e. all exons dumped onto one line using block columns at the end) 
  - to change output format, change the value of useBlocks in extract_transcript_regions.py 

