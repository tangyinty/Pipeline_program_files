#Get RTSC file from a bam file

import sys
from read_file import *
from Bio import SeqIO
import os

fasta_file = sys.argv[1]
map_file = sys.argv[2]

os.system("samtools view -F 0xfff "+map_file+"|cut -f 3,4 > map_info.txt")  #get mapped reads and mapping infomation (transcript and position mapped)

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta'); #parse reference file
length_seq = {};
for seq in fasta_sequences:
        nuc = seq.id;
        length_seq[nuc] = len(seq.seq.tostring());



mapping = {}
transcripts = []

f = open("map_info.txt");
for aline in f.readlines():   #put mapping information into a dictionary
    tline = aline.strip();
    tl = tline.split('\t');
    if tl[0].strip() not in transcripts:
        transcripts.append(tl[0].strip());
        mapping[tl[0].strip()] = [];

    mapping[tl[0].strip()].append(tl[1].strip());

distribution = {}; #Save the RTSC information in a dictionary first
coverage = {};
for transcript in length_seq:
    distribution[transcript] = [];
    for i in range(0, length_seq[transcript]): #reads distribution initialization
        distribution[transcript].append(0);
    sum_count = float(0);
    if transcript in mapping:
        for j in range(0, len(mapping[transcript])):
            index = mapping[transcript][j];
            #count = reads[mapping[transcript][j][0]];
            sum_count = sum_count + 1;
            distribution[transcript][int(index)-1] = distribution[transcript][int(index)-1] + 1; #The count of reads plus 1 if a read is mapped to the position
        coverage[transcript] = float(sum_count)/float(length_seq[transcript]);  #calculate reads coverage
    else:
        coverage[transcript] = 0

        
        
    

h = file("reads_distribution.txt", 'w')
for transcript in length_seq: #output the RTSC file as a text file
    h.write(transcript);
    h.write('\n')
    for i in range(0, length_seq[transcript]):
        h.write(str(distribution[transcript][i]))
        h.write('\t')
    h.write('\n')
    h.write('\n')



    

f.close();
h.close()




