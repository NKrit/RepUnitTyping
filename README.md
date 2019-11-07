# RepUnitTyping

## RepUnitTyping version 1.4
**RepUnitTyping.py** predicts copy numbers of repeat units in VNTR loci from PCR-free Illumina short reads.  
This script is using a backbone of SpoTyping-v2.1, a well-known in-silico spoligotyping tool (Xia E et al. Genome Med. 2016).
Entire modules for VNTR repeat-unit prediction were newly built in this script.

### Citation
[Genotyping of _Mycobacterium tuberculosis_ spreading in Hanoi, Vietnam using conventional and whole genome sequencing methods]
(https://doi.org/10.1016/j.meegid.2019.104107)

**ref/rep_unit.fasta** contains a provisional set of repeat unit sequences and their flanking sequences observed in the 33 VNTR loci of _Mycobacterium tuberculosis_ for the 24-locus MIRU-VNTR and other VNTR analyses.

You may prepare another multi-fasta file containing user-specified variants optimized for your own settings.  

Locus labels (=keys):  
[M2, 0424, ETR-C, M4, M40, M10, M16, 1955, 1982, M20, 2074, 2163a, 2163b, ETR-A, 2347, 2372, 2401, ETR-B, M23, M24, M26, M27, 3155, 3171, M31, 3232, 3336, 3690, 3820, 4052, 4120, 4156, M39]
	
* the 24 MIRU-VNTR loci = [MIRU02(M2), Mtub04(0424), ETR-C, MIRU04(M4), MIRU40(M40), MIRU10(M10), MIRU16(M16), Mtub21(1955), MIRU20(M20), QUB11b(2163b), ETR-A, Mtub29(2347), Mtub30(2401), ETR-B, MIRU23(M23), MIRU24(M24), MIRU26(M26), MIRU27(M27), Mtub34(3171), MIRU31(M31), Mtub39(3690), QUB26(4052), QUB4156(4156), MIRU39(M39)]


## Prerequisites:
    Python2.7 or 3.5(3.6)
    BLAST+ [ncbi-blast(-2.8.1+)]
    (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
    samtools-1.9
    (http://www.htslib.org/download/)
    
    python libraries used inside:
	    import sys
	    import os
	    import re
	    from optparse import OptionParser
	    import subprocess
	    import gzip
	    from collections import OrderedDict
	    
## Installation:
```
cd
git clone https://github.com/NKrit/RepUnitTyping.git
```    
## Input:
    Fastq(.gz) file or paired-end files
    Fasta file of a complete genomic sequence or assembled contigs

## Output:
    Output file specified: predicted number of repeat units in VNTR loci
    (In the non-prediction mode, the presence or absence is shown)
    Output log file: number of hits in BLAST for each repeat unit or flanking sequence
    Output log2 file: summary of search results in the fasta mode

## Usage:
```
python RepUnitTyping.py --help

    Usage: python RepUnitTyping.py [options] FASTQ_1/FASTA FASTQ_2(optional)

    Options:
  	--version             show program's version number and exit
	-h, --help            show this help message and exit
  	-p, --pred            set this if you try prediction of the number of repeat
        	              units based on hits on flanking sequences [Default is
                              off]
  	--seq                 set this if input is a fasta file that contains only a
                              complete genomic sequence or assembled contigs
                              [Default is off]
  	-s SWIFT, --swift=SWIFT
                              swift mode, either "on" or "off" [Default: off]
  	-O OUTDIR, --outdir=OUTDIR
                              output directory [Default: running directory]
  	-o OUTPUT, --output=OUTPUT
          	              basename of output files generated [Default:
          	   	      RepUnitTyping]
  	-f, --filter          stringent filtering of reads (used only for low
                              quality reads) [Default is off]
  	--sorted              set this only when the reads are sorted to a reference
      	                      genome [Default is off]
  	-d, --detail          enable detail mode, keeping intermediate files for
      	                      checking [Default is off]
  	-q Q_FASTA, --query=Q_FASTA
                              query file for repeat units [Default is rep_unit.fasta
                              in "ref" subdirectory]
  	-c CUTOFF, --cutoff=CUTOFF
        	              threshold for the presence of each sequence [Default:
                	      0.15] times the average read depth calculated from Mtb
                              genome size
```
## Examples:
```
cd ./RepUnitTyping
python RepUnitTyping.py -s off ../AL123456.3H37Rv_HS25-l150-f200_R1.fastq.gz ../AL123456.3H37Rv_HS25-l150-f200_R2.fastq.gz -q rep_unit.fasta -O RepUnit_out -o 2019RepUnitTyping -p # prediction/non-swift mode for paired fastq files
python RepUnitTyping.py -s on ../AL123456.3H37Rv_HS25-l150-f200_R1.fastq.gz ../AL123456.3H37Rv_HS25-l150-f200_R2.fastq.gz -q rep_unit.fasta -O RepUnit_out -o 2019RepUnitTyping # non-prediction/swift mode for paired fastq files 
python RepUnitTyping.py --seq ../AL123456.3H37Rv.fasta -q rep_unit.fasta -O RepUnit_out -o 2019RepUnitTyping -p # fasta mode mainly for complete genome sequence data
# rep_unit.fasta should be located in the ref subdirectory.
# Output files created from the last command are present in the RepUnit_out directory as an example. You can delete them when new analyses are made.
```
or you may use a shell script, **rep-unit-typing.sh**, to run RepUnitTyping.py interactively and repeatedly.    
```
cd ./RepUnitTyping
sh rep-unit-typing.sh
```
* For good prediction, PCR-free deep sequencing (depth of coverage > 200) is indispensable.
* When inconsistencies with experimental typing results are suspected, incomplete matches due to unidentified repeat unit variants or flanking sequences should be considered, and an optimal **rep_unit.fasta** file should be reconstructed, extracting unlisted variants from de-novo assembled sequences.
* Initial version v1.4 (2019-05-06)
