#Normalize the reactivities calculated (2%8%) with the capped value provided

import sys
from Bio import SeqIO
import math
from parse_dis_react import *

def cap(a,value): #cap the value
    if a>=value:
        return value
    else:
        return a

def react_norm(react_file, result_file, capped_value):
    react1 = parse_dist(react_file)
    react = react1[1]
    h = file(result_file, 'w')

    capped = int(capped_value)

    all_react = []


    for t in react: #put all reactivity into a list
        if react[t]!='null':
            for i in range(len(react[t])):
                if react[t][i]!='NA':
                    all_react.append(float(react[t][i]))

    all_react.sort(reverse = True) #Sort the reactivity

    eight = all_react[int(len(all_react)*0.02):int(len(all_react)*0.1)] #get the reactivities of 8%
    meight = sum(eight)/len(eight) #get the average reactivity of 8%
    print(meight)

    for t in react: #Divide the reactivity on each nt by the average reactivity of 8%
        h.write(t)
        h.write('\n')
        if react[t]!='null':
            if (t.find('AT1G29930')==-1) and (t.find('At1g29930')==-1):
                for i in range((len(react[t])-1)):
                    if react[t][i]!='NA':
                        h.write(str(cap((float(react[t][i])/meight),capped)))
                    else:
                        h.write('NA')
                    h.write('\t')
                if react[t][i+1]!='NA':
                    h.write(str(cap((float(react[t][i+1])/meight),capped)))
                else:
                    h.write('NA')
                h.write('\n')
            else:
                for i in range((len(react[t])-1)):
                    if react[t][i]!='NA':
                        h.write(str(float(react[t][i])*2.6))
                    else:
                        h.write('NA')
                    h.write('\t')
                if react[t][i+1]!='NA':
                    h.write(str(float(react[t][i])*2.6))
                else:
                    h.write('NA')
                h.write('\n')
                
                

    h.close()
        





















        





