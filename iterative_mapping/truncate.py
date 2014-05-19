#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

#truncate read
fasta_file = sys.argv[1]
shift_in = sys.argv[2]
length = sys.argv[3]
result_file = sys.argv[4]

shift = int(shift_in) #number of nt to be truncated from the 3' end of each read
length = int(length)  #The maxium length requirement to be not truncated 
    
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta');
h = file(result_file,'w')
for seq in fasta_sequences:
        nuc = seq.id;
        sequence = seq.seq.tostring();
        if (len(sequence)-shift)>=length:
                h.write('>'+nuc)
                h.write('\n')
                h.write(sequence[0:(len(sequence)-shift)])
                h.write('\n')




h.close()




