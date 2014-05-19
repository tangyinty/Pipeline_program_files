#parse RTSC file or reactivity file into a dictionary

import sys

def parse_dist(in_file):
    result = []
    distribution = {}
    name = []
    f = open(in_file)
    flag = 0
    for aline in f.readlines():
        line = aline.strip()
        dis = line.strip()
        dist = dis.split('\t') #split the line and the reactivites or reads are in a list
        if len(dist) > 0:
            if len(dist) == 1:
                if dist[0].strip().find('coverage')==-1: #if the line is the transcript id line 
                    if flag == 0:
                        name.append(line) #add the name in the name list
                        flag = 1
                        t_name = line
                    else:
                        distribution[t_name] = 'null' #if no read is mapped to the previous trasncript in the library (reactivity cannot be calculated) 
                        name.append(line)
                        flag = 1
                        t_name = line
            else:
                distri = []
                for i in range(0, len(dist)):
                    distri.append(dist[i].strip())
                distribution[t_name] = distri #add the list of reactivities/reads into a dictionary
                flag = 0
    result.append(name)
    result.append(distribution) #Output the dictionary
    f.close()
    return result
                
                







        





