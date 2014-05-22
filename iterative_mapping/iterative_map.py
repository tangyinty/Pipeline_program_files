#!/usr/bin/env python
#Iterative map the reads to reference library using bowtie

import sys
import os
from read_file import *
from read_s_file import *


seq_file = sys.argv[1]
ref_file = sys.argv[2]
shift = sys.argv[3]
length = sys.argv[4]
type_in = sys.argv[5]

if len(sys.argv)>6: #import bowtie flags
    s = ""
    for i in range(6,len(sys.argv)):
        s = s+sys.argv[i]
        s = s+" "
else:
    s = "-v 3 -a --best --strata "

if type_in == 'FASTQ':
    os.system("fastq_to_fasta -n -i "+seq_file+" -o seq0.fa") #convert to fasta file for easier further treatment
else:
    if type_in == 'FASTA':
        os.system("cp "+seq_file+" seq0.fa") #backup the sequence file
    else:
        print("Unknow sequence type!")
        sys.exit(0)
        
os.system("bowtie-build -f "+ref_file+" ref") #build up reference

k = 0
while(True):
    os.system("bowtie "+s+"-f ref"+" seq"+str(k)+".fa -S > map"+str(k)+".sam")
    os.system("samtools view -Sb -F 0xfff map"+str(k)+".sam > "+"mapped"+str(k)+".bam") #get mapped reads
    os.system("samtools view -Sb -f 0x4 map"+str(k)+".sam > "+"umapped"+str(k)+".bam") #get unmapped reads
    os.system("samtools view -Sb -f 0x10 map"+str(k)+".sam > "+"rmapped"+str(k)+".bam") #get reversed mapped reads
    os.system("samtools merge -f unmapped"+str(k)+".bam "+"umapped"+str(k)+".bam "+"rmapped"+str(k)+".bam") #combine unmapped and reversed mapped reads
    os.system("samtools view -h -o unmapped"+str(k)+".sam "+"unmapped"+str(k)+".bam") #convert to sam format for further treatment
    if k>0:
        os.system("samtools view -h -o mapped"+str(k)+".sam "+"mapped"+str(k)+".bam") #convert to sam format for further treatment
        os.system("cut -f 1 unmapped"+str(k)+".sam >"+"unmapped"+str(k)+".txt") #Get read id from unmapped reads
        os.system("cut -f 1 mapped"+str(k)+".sam >"+"mapped"+str(k)+".txt") 
        os.system("python remove_map.py "+"unmapped"+str(k)+".txt "+"mapped"+str(k)+".txt"+" runmapped"+str(k)+".txt") # Remove the reads that were mapped in the unmapped read set (The case the a read can be mapped to cDNA library both forward and backward)
        os.system("rm mapped"+str(k)+".sam") #remove the reads that won't be used 
        os.system("rm mapped"+str(k)+".txt")
    else:
        os.system("cut -f 1 unmapped"+str(k)+".sam >"+"runmapped"+str(k)+".txt")
    
    os.system("rm unmapped"+str(k)+".bam")
    os.system("rm umapped"+str(k)+".bam")
    os.system("rm rmapped"+str(k)+".bam")
    os.system("python seq_track.py runmapped"+str(k)+".txt seq"+str(k)+".fa "+" unmap_seq"+str(k)+".fa") #get unmapped read sequence
    os.system("python truncate.py unmap_seq"+str(k)+".fa "+shift+" "+length+" seq"+str(k+1)+".fa") #truncate unmapped read
    os.system("rm seq"+str(k)+".fa") #Remove sequences being mapped
    os.system("rm map"+str(k)+".sam") #Remove mapping file
    os.system("rm unmap_seq"+str(k)+".fa") #Remove unmapped sequnce
    os.system("rm unmapped"+str(k)+".txt")
    os.system("rm runmapped"+str(k)+".txt")
    os.system("rm unmapped"+str(k)+".sam")
    
    os.system("wc -l seq"+str(k+1)+".fa > count"+str(k+1)+".txt")
    c = read_sp_file("count"+str(k+1)+".txt")
    if c[0][0] == '0': #If no reads is in the sequence file, stop
        os.system("rm count"+str(k+1)+".txt")
        os.system("rm seq"+str(k+1)+".fa")
        break
    os.system("rm count"+str(k+1)+".txt")
    k = k+1

ss = ""
for i in range(0,k+1):
    ss = ss+" mapped"+str(i)+".bam"


os.system("samtools merge -f mapped_all.bam"+ss) #merge all mapped reads in bam format
os.system("rm ref*")
    

