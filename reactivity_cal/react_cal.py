#Calculate reactivity on each nucleotide

import sys
from Bio import SeqIO
import math
from parse_dis_react import *
from react_norm_function import *
import os


dist_file1 = sys.argv[1] #plus library
dist_file2 = sys.argv[2] #minus library
seq_file = sys.argv[3]
nt_spec = sys.argv[4]

if len(sys.argv)>5: #if normalization is needed and threshold is provided
    flag_in = sys.argv[5]
    threshold = sys.argv[6]
else:
    flag_in = 0
    threshold = 7


distri_p = parse_dist(dist_file1) #parse RTSC file into a dictionary
distri_m = parse_dist(dist_file2)
threshold = float(threshold) #threshold to cap the reactivity

h = file("react.txt",'w') #intermediate output file
flag_in = int(flag_in)

seqs = SeqIO.parse(open(seq_file),'fasta');
nt_s = set()
for i in range(len(nt_spec)):
    nt_s.add(nt_spec[i])

flag = 0
trans = []
distri_p = distri_p[1]
distri_m = distri_m[1]


transcripts = {}
for seq in seqs:
    n = seq.id
    trans.append(n)
    transcripts[n] = seq.seq.tostring()
    

for i in range(0, len(trans)):
    h.write(trans[i])
    h.write('\n')
    if (trans[i].find('AT1G29930')==-1) and (trans[i].find('At1g29930')==-1):        
        for j in range(len(distri_p[trans[i]])):
            distri_p[trans[i]][j] = math.log((int(distri_p[trans[i]][j])+1),math.e) #take log to the nubmer of reads mapped to each nucleotide in both libraries
        for j in range(len(distri_m[trans[i]])):
            distri_m[trans[i]][j] = math.log((int(distri_m[trans[i]][j])+1),math.e)       
        s_p = sum(distri_p[trans[i]])
        s_m = sum(distri_m[trans[i]])
        length = len(distri_p[trans[i]])
        if s_p!= 0 and s_m!= 0:
            r = []
            for j in range(0, len(distri_p[trans[i]])):
                f_p = (float(distri_p[trans[i]][j]))/float(s_p)*length #normalize the log number of reads mapped to each position
                f_m = (float(distri_m[trans[i]][j]))/float(s_m)*length
                raw_react = f_p-f_m #Subtract the nubmer of normalized reads on each nt in minus library from that in plus library
                r.append(max(0, raw_react)) #only keep positive value for reactivity
    else:
        for j in range(len(distri_p[trans[i]])):
            distri_p[trans[i]][j] = int(distri_p[trans[i]][j])
        for j in range(len(distri_m[trans[i]])):
            distri_m[trans[i]][j] = int(distri_m[trans[i]][j])       
        s_p = sum(distri_p[trans[i]])
        s_m = sum(distri_m[trans[i]])
        if s_p!= 0 and s_m!= 0:
            r = []
            for j in range(0, len(distri_p[trans[i]])):
                f_p = float(distri_p[trans[i]][j])/float(s_p)
                f_m = float(distri_m[trans[i]][j])/float(s_m)
                r.append((max(0,(f_p-f_m)))*100)
                
    if s_p!= 0 and s_m!= 0:    
        for k in range(1,(len(r)-1)):
            if transcripts[trans[i]][k-1] in nt_s:
                h.write(str(r[k]))
                h.write('\t')
            else:
                h.write('NA')
                h.write('\t')
        k = k+1
        if transcripts[trans[i]][k-1] in nt_s:
            h.write(str(r[k]))
            h.write('\n')
        else:
            h.write('NA')
            h.write('\n')
            

h.close()

if flag_in == 0:
    react_norm("react.txt","reactivity.txt", threshold) #further normalize the reactivities calculated (2%8%) with the capped value provided
    os.system("rm -f react.txt")
else:
    os.system("cp react.txt reactivity.txt") #output the unnormalized reactivities
    os.system("rm -f react.txt")
    
     
            
    
    
        





















        





