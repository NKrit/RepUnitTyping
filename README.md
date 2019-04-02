# RepUnitTyping
RepUnitTyping

$ python RepUnitTyping3.3.py --help
Usage: python RepUnitTyping3.3.py [options] FASTQ_1/FASTA FASTQ_2(optional)

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -p, --pred            Set this if you try prediction of the number of repeat
                        units based on hits on flanking sequences [Default is
                        off]
  --seq                 Set this if input is a fasta file that contains only a
                        complete genomic sequence or assembled contigs from an
                        Mtb isolate [Default is off]
  -s SWIFT, --swift=SWIFT
                        swift mode, either "on" or "off" [Default: on]
  -O OUTDIR, --outdir=OUTDIR
                        output directory [Default: running directory]
  -o OUTPUT, --output=OUTPUT
                        basename of output files generated [Default:
                        RepUnitTyping]
  -f, --filter          stringent filtering of reads (used only for low
                        quality reads)[Default is off]
  --sorted              set this only when the reads are sorted to a reference
                        genome [Default is off]
  -d, --detail          enable detail mode, keeping intermediate files for
                        checking [Default is off]
  -q Q_FASTA, --query=Q_FASTA
                        query file for repeat units [Default is rep_unit.fasta
                        in "ref" subdirectory]
