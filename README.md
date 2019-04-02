# RepUnitTyping

## RepUnitTyping
RepUnitTyping.py predicts repeat unit numbers in VNTR loci from PCR-free Illumina short reads.

This script is using a backbone of SpoTyping-v2.1, a well-known in-silico spoligotyping tool (Xia E et al. Genome Med. 2016)

Locus labels (=keys): 
['M2', '0424', 'ETR-C', 'M4', 'M40', 'M10', 'M16', '1955', '1982', 'M20', '2074', '2163a', '2163b', 'ETR-A', '2347', '2372', '2401', 'ETR-B', 'M23', 'M24', 'M26', 'M27', '3155', '3171', 'M31', '3232', '3336', '3690', '3820', '4052', '4120', '4156', 'M39']

## Prerequisites:
    Python2.7
    BLAST+ (ncbi-blast-2.4.0+ to -2.8.1+)
    
## Input:
    Fastq file or pair-end fastq files
    Fasta file of a complete genomic sequence or assembled contigs of an Mtb isolate

## Output:
    In the output file specified: predicted number of repeat units (VNTR loci)
    In the output log file: count of hits from BLAST result for each repeat unit and flanking sequence.
    In the output log2 file: candidate repeat variants in fasta mode

## Usage
$ python RepUnitTyping.py --help

    Usage: python RepUnitTyping.py [options] FASTQ_1/FASTA FASTQ_2(optional)

    Options:
    --version           show program's version number and exit
    -h, --help          show this help message and exit
  
    -p, --pred          Set this if you try prediction of the number of repeat
                        units based on hits on flanking sequences [Default is
                        off]
                        
    --seq               Set this if input is a fasta file that contains only a
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
                        
    -d, --detail          enable detail mode, keeping intermediate files for
                        checking [Default is off]
                        
    -q Q_FASTA, --query=Q_FASTA
                          query file for repeat units [Default is rep_unit.fasta
                        in "ref" subdirectory]
                        
## Example
    $ cd ./RepUnitTyping
    $ python RepUnitTyping.py -s off ../AL123456.3H37Rv_fasta-l150-f400_R1.fq.gz ../AL123456.3H37Rv_fasta-l150-f400_R2.fq.gz -q rep_unit.fasta -O RepUnit_out -o 20190403RepUnitTyping -p
    $ # rep_unit.fasta should be located in the ref subdirectory.
or you can use a shell script, rep-unit-typing.sh, to run RepUnitTyping.py
