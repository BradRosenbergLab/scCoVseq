#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Phillip Cohen
"""

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys, gzip, fileinput

###Check if R1 fastq has been given
###Code from Jorg's scTCRseq Pipeline Scripts
if not len([i for i in sys.argv if 'R1_input_fastq=' in i])==1:raise Exception('no R1 fastq file')
R1_input_fastq=[i for i in sys.argv if 'R1_input_fastq=' in i][0].replace('R1_input_fastq=','')
if not R1_input_fastq[-9:]=='.fastq.gz':raise Exception('R1 fastq must be a .fastq.gz file')

###Write new R1 fastq file
R1_title = str(R1_input_fastq)
R1_title = R1_title[:-23] + "reformat_" + R1_title[-23:]
R1_reformat_fastq = gzip.open('%s' %R1_title, "wt")

###Open R1 fastq file
R1_input_fastq = gzip.open('%s' %R1_input_fastq, "rt")

###Write R2 fastq file
R2_title = R1_title.replace("R1", "R2", 1)
R2_fastq = gzip.open('%s' %R2_title, "wt")

###Store TSO Sequence
tso_sequence = "TTTCTTATATGGG"
###Calculate the length of the TSO sequence
tso_length = len(tso_sequence)
###Add the length of the TSO to the 26 nt for UMI + cell barcode
length_to_cut_from_R2 = tso_length + 26

while 1:
    line_1 = R1_input_fastq.readline()
    if not line_1:
        R1_input_fastq.close()
        R2_fastq.close()
        break
    line_2 = R1_input_fastq.readline()
    line_3 = R1_input_fastq.readline()
    line_4 = R1_input_fastq.readline()
    ###Write line 1 to R1
    R1_reformat_fastq.write(line_1)
    ###Reformat line 1 in R1 for R2
    R2_line_1 = line_1.replace(" 1:N:0", " 2:N:0")
    R2_fastq.write(R2_line_1)
    ###Split line 2 into R1 and R2
    R1_line_2 = line_2[:26] + "\n"
    R2_line_2 = line_2[length_to_cut_from_R2:]
    ###Take the reverse complement of R2 line 2
    R2_line_2 = Seq(R2_line_2, IUPAC.unambiguous_dna)
    R2_line_2 = R2_line_2.reverse_complement()
    R2_line_2 = str(R2_line_2) + "\n"
    R2_line_2 = R2_line_2.lstrip()
    ###Write line 2 to R1 and R2
    R1_reformat_fastq.write(R1_line_2)
    R2_fastq.write(R2_line_2)
    ###Write line 3 to R1 and R2
    R1_reformat_fastq.write(line_3)
    R2_fastq.write(line_3)
    ###Split line 4 into R1 and R2
    R1_line_4 = line_4[:26] + "\n"
    R2_line_4 = line_4[:length_to_cut_from_R2 - 1:-1] + "\n"
    R2_line_4 = R2_line_4.lstrip()
    ###Write line 4 to R1 and R2
    R1_reformat_fastq.write(R1_line_4)
    R2_fastq.write(R2_line_4)

R1_reformat_fastq.close()
R1_input_fastq.close()

