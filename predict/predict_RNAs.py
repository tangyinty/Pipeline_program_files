#RNA structure prediction & Output and illustrate reactivities

import sys
from parse_dis_pac import *
from read_file import *
from Bio import SeqIO
import os
from rtts_plot import *


id_file = sys.argv[1]
seq_file = sys.argv[2]


flag = 0
if len(sys.argv)>=4: #input reactivity file if provided
    react_file = sys.argv[3]
    react = parse_dist(react_file)
    react = react[1]
    flag = 1
    h = file("transcript_reactivities.txt", 'w')


ids = read_t_file(id_file)
sequences = SeqIO.parse(seq_file, 'fasta')


seqs = {}
for seq in sequences:
    seqs[seq.id] = seq.seq.tostring()

if len(ids)>10: #setup a limit of the number of sequence to be predicted
    print("Number of sequences exceeds limitation!")
    sys.exit(0)
    

#predict RNA structures
for i in range(len(ids)):
    id_s = ids[i][0]
    #Put RNA sequence and reactivities into files
    if id_s in seqs:
        f = file("temp.txt", 'w')
        f.write('>'+id_s)
        f.write('\n')
        f.write(seqs[id_s])
        f.close()
        if flag == 0:
            os.system("fold "+"temp.txt"+" "+id_s+".ct")
        if flag == 1:
            if id_s in react:
                f = file("constraint.txt",'w')
                make_plot(react[id_s],id_s) #make a plot of the distribution of the reactivites of the input RNA
                h.write(id_s)
                h.write('\n')
                for j in range(0, (len(react[id_s]))):
                    if react[id_s][j]!='NA':
                        f.write(str(j+1))
                        f.write('\t')
                        f.write(str(react[id_s][j]))
                        f.write('\n')
                        h.write(str(react[id_s][j])) #Output the reactivities
                        h.write('\t')
                f.close()
                h.write('\n')
                h.write('\n')
                os.system("fold "+"temp.txt"+" -sh"+" "+"constraint.txt"+" "+id_s+".ct")
            else:
                print(id_s+" not in the data of react!")
        os.system("draw "+id_s+".ct "+id_s+".ps")
    else:
        print(id_s+" not in the data of sequences!")

#Remove the unnecessary files
os.system("rm -f temp.txt")
if flag == 1:
    os.system("rm -f constraint.txt")

h.close()
    



