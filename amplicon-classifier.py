#/usr/bin/env python

import pydna
from Bio import SeqIO
from Bio import Emboss


def fasta_seq(fasta):
    """
    Open fasta format file and returns the sequence.
    """
    f = SeqIO.parse(open(fasta, 'rU'), 'fasta').next()
    return f.seq



#######################################################
gb = pydna.Genbank("andreas.sjodin@gmail.com") # Tell Genbank who you are!

gene = gb.nucleotide("X06997") # Kluyveromyces lactis LAC12 gene for lactose permease.

#######################################################


genome = pydna.read('FSC237.fasta')

#genome = fasta_seq('FSC237.fasta')
primer_f,primer_r = pydna.parse(''' >B4_400_1-F
                                    AGCAGTGCCTGTTGTACC

                                    >B4_400_1-R
                                    AGTTTCTCAACATGGAAT
                                    ''', ds=False)

pcr_prod = pydna.pcr(primer_f,primer_r, genome)



#stssearch  -seqall genome.fasta -infile primers.txt -stdout --auto
primersearch  -seqall genome.fasta -infile primers.txt -mismatchpercent 1 -stdout --auto -mismatchpercent 10
primersearch  -seqall ../../rawdata/FSC770.1.fna -infile primers.txt -mismatchpercent 1 -stdout --auto -mismatchpercent 20 -scircular1



