##
## [RepUnitTyping]
##

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## http://www.opensource.org/licenses/gpl-3.0.html

## The backbone of this program is SpoTyping-v2.1 developed for in-silico spoligotyping by Xia Eryu (Genome Medicine, 2016)
## Entire modules for VNTR repeat-unit prediction were newly created and merged with this backbone by N Keicho (nkeicho-tky@umin.ac.jp) 
## temporary version 3.8 (20190413) -v1.1

## --------------------------------

import sys
import os
import re
from optparse import OptionParser
import subprocess
import gzip
from collections import OrderedDict

## Global variables
dir = os.path.split(os.path.realpath(__file__))[0] # script directory
setlength = 50*4500000                             # base input cut-off for swift mode

## Option variables
usage = "usage: python %prog [options] FASTQ_1/FASTA FASTQ_2(optional)"
parser = OptionParser(usage=usage, version="v1.1")
parser.add_option("-p","--pred",action="store_true",dest="pred",help="set this if you try prediction of the number of repeat units based on hits on flanking sequences [Default is off]")
parser.add_option("--seq",action="store_true",dest="seq",help="set this if input is a fasta file that contains only a complete genomic sequence or assembled contigs from an Mtb isolate [Default is off]")
parser.add_option("-s","--swift",action="store",type="string",dest="swift",default="off",help="swift mode, either \"on\" or \"off\" [Default: off]")
parser.add_option("-O","--outdir",action="store",type="string",dest="outdir",default=".",help="output directory [Default: running directory]")
parser.add_option("-o","--output",action="store",type="string",dest="output",default="RepUnitTyping",help="basename of output files generated [Default: RepUnitTyping]")
parser.add_option("-f","--filter",action="store_true",dest="filtQ",help="stringent filtering of reads (used only for low quality reads)[Default is off]")
parser.add_option("--sorted",action="store_true",dest="sortS",help="set this only when the reads are sorted to a reference genome [Default is off]")
parser.add_option("-d","--detail",action="store_true",dest="detail",help="enable detail mode, keeping intermediate files for checking [Default is off]")
parser.add_option("-q","--query",action="store",type="string",dest="q_fasta",default="rep_unit.fasta",help="query file for repeat units [Default is rep_unit.fasta in \"ref\" subdirectory]")
(options, args) = parser.parse_args()

filt = options.filtQ
sortS = options.sortS
pred = options.pred                     # prediction of the number of repeat units, when reference file contains the 5'- and 3'-flanking sequences
seq = options.seq                       # input file contains a genomic DNA sequence
swift = options.swift                   # swift mode, default: off
outdir = options.outdir                 # output directory
output = options.output                 # basename of output files
detail = options.detail                 # detail mode
q_fasta = options.q_fasta               # catalogue file for repeat units

ver = parser.get_version()
print("version: "+ver)
## Check input fastq files
narg = len(args)
print("no. args = "+str(narg))
if narg == 0:
    print("usage: python RepUnitTyping.py [options] FASTQ_1/FASTA FASTQ_2(optional)")
    print("     -h, --help for further help message")
    print("     --version for this program's version")
    sys.exit()
elif narg == 1:
    input1 = args[0]    # Input fastq file 1
    if not os.path.isfile(input1):
        print("ERROR: Invalid FASTQ_1/FASTA file!")
        sys.exit()
elif narg == 2:
    if seq:
        print("ERROR: Only one fasta file is allowed when '--seq' is set!")
        sys.exit()
    input1 = args[0]    # Input fastq file 1
    input2 = args[1]    # Input fastq file 2
    if not os.path.isfile(input1):
        print("ERROR: Invalid FASTQ_1 file!")
        sys.exit()
    elif not os.path.isfile(input2):
        print("ERROR: Invalid FASTQ_2 file!")
        sys.exit()


############################################################
## Main class is described here
############################################################
class Main:
    def filterQuality(self, instr):
        '''
        This function trims/filters out low quality reads that satisfies one of the conditions below:
        1. Leading and trailing 'N's would be removed.
        2. Any read with more than 3 'N's in the middle would be removed.
        3. Any read with more than 7 consecutive bases identical would be trimmed/filtered out given
           the length of the flanking regions.
        '''
        instr = instr.upper()
        instr = instr.strip('N')
        if instr.count('N') >= 3:
            return ""
        tmp = re.split('A{7,}|T{7,}|C{7,}|G{7,}',instr)
        if len(tmp) == 1:
            return instr.replace('N', '')
        output = []
        for item in tmp:
            if len(item) >= 25:
                item = item.replace('N','')
                output.append(item)
        return "".join(output)


    def concatenation(self,pair_n,in_file,out_handle):  # Concatenate sequences without length check
        if in_file.endswith(".gz"):
            in_handle = gzip.open(in_file, 'rb')
        else:
            in_handle = open(in_file)
        count = 0
        outlength = 0
        splitoff = 450000000 # 450MB=4.5MB x 100 coverage depth
        times = 1
        for line in in_handle:
            line = line.decode('utf8').strip('\n')
            if outlength > splitoff * times:
                out_handle.write("\n>Concat"+str(pair_n)+"_"+str(times)+"\n")
                times += 1
            if count % 4 == 1:
                if filt:
                    line = self.filterQuality(line)
                out_handle.write(line)
                # outlength does not include labels' size
                outlength += len(line)
            count = (count+1) %4
        in_handle.close()
        return outlength

    def concatenation_check(self,in_file,out_handle,cutoff):  # Concatenate sequences with length check
        if in_file.endswith(".gz"):
            in_handle = gzip.open(in_file, 'rb')
        else:
            in_handle = open(in_file)
        count = 0
        outlength = 0
        for line in in_handle:
            line = line.decode('utf8').strip('\n')
            if outlength > cutoff:
                break
            elif count % 4 == 1:
                if filt:
                    line = self.filterQuality(line)
                out_handle.write(line)
                outlength += len(line)
            count = (count+1) %4
        in_handle.close()
        return outlength


    def scan_file(self,in_file,setlength):
        '''
        This function is designed to deal with situations where sequence reads are sorted
        (extracrted from sorted bam files, for example)
        '''
        totalBase = 0

        if in_file.endswith(".gz"):
            in_handle = gzip.open(in_file, 'rb')
        else:
            in_handle = open(in_file)
        count = 0

        for line in in_handle:
            line = line.decode('utf8').strip('\n')
            if count % 4 == 1:
                totalBase += len(line)
            count = (count+1) % 4
        in_handle.close()

        if totalBase < setlength:
            return 1
        return int(totalBase*1.0/setlength+0.5)


    def concatenation_sort(self,in_file,out_handle,setlength,step,start):
        '''
        This function is designed to deal with situations where sequence reads are sorted
        (extracrted from sorted bam files, for example)
        '''
        print("Sorted reads ..")
        if in_file.endswith(".gz"):
            in_handle = gzip.open(in_file, 'rb')
        else:
            in_handle = open(in_file)

        count = 0
        downsample = 0
        outlength = 0

        for line in in_handle:
            line = line.decode('utf8').strip('\n')
            if outlength > setlength:
                break
            elif count % 4 == 1 and downsample % step == start:
                if filt:
                    line = self.filterQuality(line)
                out_handle.write(line)
                outlength += len(line)
                downsample = (downsample+1) %step
            elif count %4 == 1:
                downsample = (downsample+1) %step
            count = (count+1) %4
        in_handle.close()
        return outlength


    def parse_blast(self,in_file,qlabels,log_handle,out_handle):  # Parse blast output file
        #Using ordered dictionary
        record = OrderedDict()
        record_relax = OrderedDict()
        query_n = len(qlabels)
        print("\nNo. query sequences: "+str(query_n))
        #print("\nReference labels: ")
        #print(qlabels)
        for i in range(query_n):
            record[qlabels[i]] = 0
            record_relax[qlabels[i]] = 0

        # 13th column (non-standard option)
        file_blast = open(in_file)
        for line in file_blast:
            line = line.strip("\n")
            if not re.search('#',line):
                tmp = re.split('\s+',line)
                # tmp[2]=pident (Percentage of identical matches)
                tmp[2] = float(tmp[2])
                # tmp[3]=length (Alignment length)
                tmp[3] = int(tmp[3])
                # tmp[12]=qlen (Query sequence length)
                tmp[12]=int(tmp[12])
                thres = 100*(float(tmp[12])-1)/float(tmp[12]) 
                if tmp[2]==100 and tmp[3]==tmp[12]:
                    record[tmp[0]]+=1
                    record_relax[tmp[0]]+=1
                elif (tmp[2]>=thres and tmp[3]==tmp[12]) or (tmp[2]==100 and tmp[3]==tmp[12]-1):
                    record_relax[tmp[0]]+=1
        file_blast.close()

        storage=[]
        for i in range(query_n):
            signal =0
            # min_strict Default: 0.1 * average read depth
            if (record[qlabels[i]] >= min_strict):
                signal = 1
            storage.append(signal)
            log_handle.write("%s\t%d\t%d\t%d\n" % (qlabels[i],record[qlabels[i]],record_relax[qlabels[i]],signal))
        print("\nStorage: ")
        print(storage)
        print("\nmin_strict: "+str(min_strict))

        out_handle.write("%s\n" % '\t'.join([str(item) for item in storage]))
        out_handle.write("%s\t%s\n" % ("ID", '\t'.join([qlabel for qlabel in qlabels])))


    def parse_blast_pred(self,in_file,qlabels,log_handle,out_handle):  # Parse blast output file for prediction of the  number of repeat units
        #Using ordered dictionary
        print("[NOTES]")
        print("Flanking sequence labels = xxx_FL_01, xxx_FL_02 ..")
        print("Repeat sequence labels = xxx_01, xxx_02, xxx_03 ...")
        record = OrderedDict()
        record_relax = OrderedDict()
        record_95 = OrderedDict()
        query_n = len(qlabels)
        print("\nNo. query sequences: "+str(query_n))
        #print("\nReference labels: ")
        #print(qlabels)
        for i in range(query_n):
            record[qlabels[i]] = 0
            record_relax[qlabels[i]] = 0
            record_95[qlabels[i]] = 0

        # 13th column (non-standard option)
        file_blast = open(in_file)
        for line in file_blast:
            line = line.strip("\n")
            if not re.search('#',line):
                tmp = re.split('\s+',line)
                #print(tmp)
                #print(record[tmp[0]])
                # tmp[2]=pident (Percentage of identical matches)
                tmp[2] = float(tmp[2])
                # tmp[3]=length (Alignment length)
                tmp[3] = int(tmp[3])
                # tmp[12]=qlen (Query sequence length)
                tmp[12]=int(tmp[12])
                thres = 100*(float(tmp[12])-1)/float(tmp[12])
                #print("\nthres: ")
                #print(thres)
                if tmp[2]==100 and tmp[3]==tmp[12]:
                    record[tmp[0]]+=1
                    record_relax[tmp[0]]+=1
                    record_95[tmp[0]]+=1
                elif (tmp[2]>=thres and tmp[3]==tmp[12]) or (tmp[2]==100 and tmp[3]>=tmp[12]-1):
                    record_relax[tmp[0]]+=1
                    record_95[tmp[0]]+=1
                #elif (tmp[2]>=95 and tmp[3]==tmp[12]) or (tmp[2]==100 and tmp[3]>=int(0.95*tmp[12])):
                elif (tmp[2]>=95) and (tmp[3]>=int(0.95*tmp[12])):
                    record_95[tmp[0]]+=1

        file_blast.close()

        FlankingHits = OrderedDict()
        FlankingSqs = OrderedDict()
        RepeatHits = OrderedDict()
        FLnums = OrderedDict()
        storage=[]
        for i in range(query_n):
            signal =0
            # min_strict Default: 0.1 * average read depth
            if (record[qlabels[i]] >= min_strict):
                signal = 1
            log_handle.write("%s\t%d\t%d\t%d\t%d\n" % (qlabels[i],record[qlabels[i]],record_relax[qlabels[i]],record_95[qlabels[i]],signal))

        qlabels_sep=[]
        for i in range(query_n):
            qlabels_sep.append(qlabels[i].split("_"))
            FlankingHits[qlabels_sep[i][0]] = 0.0
            FlankingSqs[qlabels_sep[i][0]] = 0.0
            RepeatHits[qlabels_sep[i][0]] = []
            FLnums[qlabels_sep[i][0]] = 0
        print("\nNo. query loci: "+str(len(FlankingHits)))
        print("\nLocus labels (=keys): ")
        print("["+', '.join(FlankingHits.keys())+"]")

        for i in range(query_n):
            #print(qlabels_sep[i])
            if (qlabels_sep[i][1]=="FL"):
            # min_strict Default: 0.1 * average read depth
                if (record[qlabels[i]] >= min_strict):
                    FlankingHits[qlabels_sep[i][0]] += float(record[qlabels[i]])
                    FlankingSqs[qlabels_sep[i][0]] += float(record[qlabels[i]])**2
                    FLnums[qlabels_sep[i][0]] += 1
            else:
                RepeatHits[qlabels_sep[i][0]].append(record[qlabels[i]])

        log_handle.write("## %s" % "[PREDICTION]")
        log_handle.write("\n## %s\t%s\t%s\t%s\n" % ("RepLocus","NoRepeats","NoFL","CV>0.2_unstableFL"))
        #print(RepeatHits)
        for qgroup in FlankingHits:
            FLnum=FLnums[qgroup]
            if FLnum>0:
                FLmean = FlankingHits[qgroup]/FLnum
                #signal_unit = [round(x/FLmean,0) for x in RepeatHits[qgroup]]
                signal_unit = sum([(round(x/FLmean,0)) for x in RepeatHits[qgroup]])
                #print(signal_unit)
                if FLnum>1:
                    # SD = ([sum of squares divided by N] minus [square of mean])^0.5
                    coeff_var = (((FlankingSqs[qgroup]/FLnum)-(FLmean**2))**0.5)/FLmean
                    storage.append([int(signal_unit),round(FLmean,1),round(coeff_var,1)])
                    log_handle.write("%s\t%d\t%d\t%3.1f\n" % (qgroup,signal_unit,FLnum,coeff_var))
                else:
                    storage.append([int(signal_unit),round(FLmean,1),"NA"])
                    log_handle.write("%s\t%d\t%d\t%s\n" % (qgroup,signal_unit,FLnum,"NA"))
            else:
                storage.append(["NA","NA","NA"])
                log_handle.write("%s\t%s\t%d\t%s\n" % (qgroup,"NA",FLnum,"NA"))

        #print("\nSum of Flanking Hits: ")
        #print("["+', '.join([str(x) for x in FlankingHits.values()])+"]")
        print("\nRepeat Hits: ")
        print("["+', '.join([str(x) for x in RepeatHits.values()])+"]")
        print("\nMean of Flanking Hits: ")
        print([item[1] for item in storage])
        print("\nCountable Number of Flanking Sequences: ")
        print("["+', '.join([str(x) for x in FLnums.values()])+"]")
        print("\nCoefficient of Variation (>0.2 Unstable Flanking Hits): ")
        print([item[2] for item in storage])
        print("\nPredicted Number of Repeats: ")
        print([item[0] for item in storage])

        out_handle.write("%s\n" % '\t'.join([str(item[0]) for item in storage]))
        out_handle.write("%s\t%s\n" % ("#CV", '\t'.join([str(item[2]) for item in storage])))
        out_handle.write("%s\t%s\n" % ("#ID", '\t'.join([qgroup for qgroup in RepeatHits])))


    def parse_blast_log2(self,input_file, in_file,log2_handle):
        #Using ordered dictionary
        knownqs = []
        unknownqs = []
        # 13th column (non-standard option)
        file_blast = open(in_file)
        for line in file_blast:
            line = line.strip("\n")
            if not re.search('#',line):
                tmp = re.split('\s+',line)
                # tmp[2]=pident (Percentage of identical matches)
                tmp[2] = float(tmp[2])
                # tmp[3]=length (Alignment length)
                tmp[3] = int(tmp[3])
                # tmp[12]=qlen (Query sequence length)
                tmp[12]=int(tmp[12])
                # two mismatches allowed
                thres = 100*(float(tmp[12])-2)/float(tmp[12]) 
                if (tmp[2]==100 and tmp[3]==tmp[12]):
                    knownqs.append([tmp[0],input_file,tmp[1],tmp[8],tmp[9],tmp[2],round(100.0,1),"Listed"])
                elif (tmp[2]>=thres and tmp[3]>=tmp[12]-2):
                    tmp[0]= str(tmp[0].split("_")[0])
                    matched_len=round(float(tmp[3])*100/float(tmp[12]),1)
                    unknownqs.append([tmp[0],input_file,tmp[1],tmp[8],tmp[9],tmp[2],matched_len,"Unlisted"])
        for unknownq in unknownqs:
            for knownq in knownqs: 
                if knownq[2]==unknownq[2]:
                    if (knownq[3]==unknownq[3] or knownq[4]==unknownq[4]):
                        unknownq[7]="Duplicated"
        print("\nListed hits: "),
        print(sum(knownq[7].count("Listed") for knownq in knownqs))
        print("Dupicated hits: "),
        print(sum(unknownq[7].count("Duplicated") for unknownq in unknownqs))
        print("Possibly unlisted hits: "),
        print(sum(unknownq[7].count("Unlisted") for unknownq in unknownqs))

        file_blast.close()

        for knownq in knownqs:
            log2_handle.write("%s\n" % '\t'.join([str(elem) for elem in knownq]))
            #print(knownq)
        for unknownq in unknownqs:
            log2_handle.write("%s\n" % '\t'.join([str(elem) for elem in unknownq]))
            #print(unknownq)

############################################################
############################################################
## Code starts here
############################################################
if __name__ == "__main__":
    t = Main()

    ## Check name of existing tmp files
    tmpnum = 0
    while os.path.isfile("%s/%s.RepUnitTyping.tmp.%d" % (outdir,output,tmpnum)):
        tmpnum+=1
    tmpfile = "%s/%s.RepUnitTyping.tmp.%d" % (outdir,output,tmpnum)

    ##########################################################
    ## Create a fasta file with the reads concatenated.
    ##########################################################
    if not seq:
        file_tmp = open(tmpfile, 'w')
        file_tmp.write(">Concat\n")

        if swift == 'on':
            ## Deal with sorted
            if sortS:
                step = t.scan_file(input1,setlength)
                if step == 1:
                    out_first = t.concatenation_check(input1,file_tmp,setlength)
                else:
                    out_first = t.concatenation_sort(input1,file_tmp,setlength,step,0)
                    if (narg == 1) and out_first < setlength:
                        t.concatenation_sort(input1,file_tmp,setlength-out_first,step,1)
            ## Not sorted
            else:
                out_first = t.concatenation_check(input1,file_tmp,setlength)

            remaining = setlength - out_first

            if (narg == 2) and (remaining > 0):
                out_second = t.concatenation_check(input2,file_tmp,remaining)
                if out_second < remaining:
                    min_strict = max(1,int((out_first + out_second) * 0.1 / 4500000))
            elif (narg == 1) and (remaining > 0):
                    min_strict = max(1,int(out_first * 0.1 / 4500000))
            else:
                    min_strict = max(1,int(setlength*0.1 / 4500000))
        else:
            out_first = t.concatenation(1,input1,file_tmp)
            if narg == 2:
                out_first += t.concatenation(2,input2,file_tmp)
            min_strict = max(1,int(out_first * 0.1 / 4500000))

        print("\n*** Necessary conditions for accurate prediction of copy numbers")
        print(" 1. PCR-free library preparation")
        print(" 2. coverage depth > 200\n")
        print("out_first: "+str(out_first)+"B in total")
        sub_label = subprocess.check_output(["grep",">",tmpfile])
        print("split labels by 450 MB: \n"+str(sub_label.decode('utf8')))

        file_tmp.write("\n")
        file_tmp.close()

    ##########################################################
    ## Blast the repeat units against the concatenated fasta file.
    ##########################################################
    blastDB = tmpfile
    if seq:
        blastDB = input1
        min_strict = 1

    # max file size = 1GB (default)
    maxfilesz = "1GB"
    threads = "4"
    tmpH = open("%s.blast.out" % tmpfile, 'w')
    subprocess.call(["makeblastdb", "-in", blastDB, "-out", blastDB, "-dbtype", "nucl", "-max_file_sz", maxfilesz])
    print("makeblastdb -in "+blastDB+" -out "+blastDB+" -dbtype nucl -max_file_sz "+maxfilesz)
    outfmt_opt = "7 std qlen"
    subprocess.call(["blastn", "-query", "%s/ref/%s" % (dir, q_fasta), "-db", blastDB, "-task", "blastn", "-dust", "no", "-outfmt", outfmt_opt, "-max_target_seqs", "1000000", "-num_threads", threads], stdout=tmpH)
    tmpH.close()

    ##########################################################
    ## Parsing blast output & write to the output directory
    ##########################################################
    log = open("%s/%s.log" % (outdir,output),'a')
    out_file = open("%s/%s" % (outdir,output),'a')
    #log2name= outdir + '/' + output + '.log2'
    log2 = open("%s/%s.log2" % (outdir,output),'a')

    if narg == 2:
        log.write("## %s %s\n" % (input1,input2))
        log2.write("## %s %s\n" % (input1,input2))
        out_file.write("%s&%s\t" % (os.path.basename(input1),os.path.basename(input2)))
    else:
        log.write("## %s\n" % input1)
        log2.write("## %s\n" % input1)
        out_file.write("%s\t" % os.path.basename(input1))

    log.write("## query=%s\n" % q_fasta)
    log2.write("## query file = %s\n" % q_fasta)
    log.write("## ver=%s | minimum error-free hits=%d | pred=%s | seq=%s | swift=%s | filter=%s| sorted=%s\n" %(ver, min_strict, pred, seq, swift, filt, sortS))

    # query labels without '>' are saved with subprocess.check_out command
    ref_label = subprocess.check_output(["grep",">", "%s/ref/%s" % (dir, q_fasta)])
    ref_label = ref_label.decode('utf8').replace(" ","").replace(">","")
    ref_label = ref_label.split()

    if pred:
        log.write("## RepUnit\tCompleteHits\tOneErrorHits\t95PerCentHits\tSignal01\n")
        t.parse_blast_pred("%s.blast.out" % tmpfile,ref_label,log,out_file)
    else:
        log.write("## RepUnit\tCompleteHits\tOneErrorHits\tSignal01\n")
        t.parse_blast("%s.blast.out" % tmpfile,ref_label,log,out_file)

    if seq:
        log2.write("## Query\tInput_file\tSubject\tSStart\tSEnd\t%identity\tAlign_Query_R\tNote\n")
        t.parse_blast_log2(input1,"%s.blast.out" % tmpfile,log2)
    else:
        log2.write("## information about sequence similarity is obtained only in --seq mode \n")

    log2.close()
    out_file.close()
    log.close()

    print("\n")

    ##########################################################
    ## Cleaning up
    ##########################################################
    if not detail:
        post = ['','.blast.out','.reference']
        for i in post:
            if os.path.isfile("%s%s" % (tmpfile,i)):
                os.remove("%s%s" % (tmpfile,i))

    post2 = ['.nsq','.nhr','.nin']
    for i in post2:
        if os.path.isfile("%s%s" % (blastDB,i)):
            os.remove("%s%s" % (blastDB,i))
