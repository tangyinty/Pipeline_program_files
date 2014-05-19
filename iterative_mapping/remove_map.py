#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from read_file import *

# Remove the reads that were mapped in the unmapped read set (The case the a read can be mapped to cDNA library both forward and backward)
unmap_file = sys.argv[1]
map_file = sys.argv[2]
result_file = sys.argv[3]


unmap = read_t_file(unmap_file)
mapped = read_t_file(map_file)
h = file(result_file, 'w')

maps = set()
for i in range(len(mapped)):
    maps.add(mapped[i][0]) #put all mapped reads in a set


for i in range(len(unmap)):
    name = unmap[i][0]
    if name not in maps: #remove the mapped reads
        h.write(name)
        h.write('\n')


h.close()
