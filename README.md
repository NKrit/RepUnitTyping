# RepUnitTyping

## RepUnitTyping
RepUnitTyping.py predicts copy numbers of repeat units in VNTR loci from PCR-free Illumina short reads.

This script is using a backbone of SpoTyping-v2.1, a well-known in-silico spoligotyping tool (Xia E et al. Genome Med. 2016).
Entire modules for VNTR repeat-unit prediction were newly built in this script.

Locus labels (=keys): 
[M2, 0424, ETR-C, M4, M40, M10, M16, 1955, 1982, M20, 2074, 2163a, 2163b, ETR-A, 2347, 2372, 2401, ETR-B, M23, M24, M26, M27, 3155, 3171, M31, 3232, 3336, 3690, 3820, 4052, 4120, 4156, M39]

As reference, MIRU24 = [MIRU02(M2), Mtub04(0424), ETR-C, MIRU04(M4), MIRU40(M40), MIRU10(M10), MIRU16(M16), Mtub21(1955), MIRU20(M20), QUB11b(2163b), ETR-A, Mtub29(2347), Mtub30(2401), ETR-B, MIRU23(M23), MIRU24(M24), MIRU26(M26), MIRU27(M27), Mtub34(3171), MIRU31(M31), Mtub39(3690), QUB26(4052), QUB4156(4156), MIRU39(M39)]

## Prerequisites:
    Python2.7 or 3.5
    BLAST+ [ncbi-blast-2.4.0+ or the latest (-2.8.1+)]
    (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
    python libraries used inside:
	    import sys
	    import os
	    import re
	    from optparse import OptionParser
	    import subprocess
	    import gzip
	    from collections import OrderedDict
	    
## Installation:
    $ cd
    $ git clone https://github.com/NKrit/RepUnitTyping.git
    
## Input:
    Fastq file or pair-end fastq files
    Fasta file of a complete genomic sequence or assembled contigs of an Mtb isolate

## Output:
    In the output file specified: predicted number of repeat units (VNTR loci)
    In the output log file: count of hits from BLAST result for each repeat unit and flanking sequence
    In the output log2 file: candidate repeat variants in fasta mode

## Usage:
$ python RepUnitTyping.py --help

    Usage: python RepUnitTyping.py [options] FASTQ_1/FASTA FASTQ_2(optional)

    Options:
    --version           show program's version number and exit
    -h, --help          show this help message and exit
  
    -p, --pred          set this if you try prediction of the number of repeat
                        units based on hits on flanking sequences [Default is
                        off]
                        
    --seq               set this if input is a fasta file that contains only a
                        complete genomic sequence or assembled contigs from an
                        Mtb isolate [Default is off]
                        
    -s SWIFT, --swift=SWIFT
                        swift mode, either "on" or "off" [Default: on]
                        
    -O OUTDIR, --outdir=OUTDIR
                        output directory [Default: running directory]
                        
    -o OUTPUT, --output=OUTPUT
                        basename of output files generated [Default:
                        RepUnitTyping]
                        
    -f, --filter        stringent filtering of reads (used only for low
                        quality reads)[Default is off]
                        
    --sorted            set this only when the reads are sorted to a reference
                        genome [Default is off]
                        
    -d, --detail        enable detail mode, keeping intermediate files for
                        checking [Default is off]
                        
    -q Q_FASTA, --query=Q_FASTA
                        query file for repeat units [Default is rep_unit.fasta
                        in "ref" subdirectory]
                        
## Examples:
    $ cd ./RepUnitTyping
    $ python RepUnitTyping.py -s off ../AL123456.3H37Rv7.5M_R1.fastq.gz ../AL123456.3H37Rv7.5M_R2.fastq.gz -q rep_unit.fasta -O RepUnit_out -o 2019RepUnitTyping -p # prediction/non-swift mode for paired fastq files
    $ python RepUnitTyping.py -s on ../AL123456.3H37Rv7.5M_R1.fastq.gz ../AL123456.3H37Rv7.5M_R2.fastq.gz -q rep_unit.fasta -O RepUnit_out -o 2019RepUnitTyping # non-prediction/swift mode for paired fastq files 
    $ python RepUnitTyping.py --seq ../AL123456.3H37Rv.fasta -q rep_unit.fasta -O RepUnit_out -o 2019RepUnitTyping -p # fasta mode mainly for complete genome sequence data
    $ # rep_unit.fasta should be located in the ref subdirectory.
    $ # Output files created from the last command are kept in the RepUnit_out directory as an example. Please delete them when new analyses are made.
    
or you may use a shell script, rep-unit-typing.sh, to run RepUnitTyping.py more interactively.    

    $ cd ./RepUnitTyping
    $ sh rep-unit-typing.sh

For good prediction, PCR-free deep sequencing (coverage depth > 200) is required.  
