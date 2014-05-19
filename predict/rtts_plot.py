#!/usr/bin/env python
#Make a plot of reactivity distribution

import sys
import numpy as np
import matplotlib
from pylab import *

#Convert the reactivities (Make NA to 0)
def convert_react(a):
    r = []
    for i in range(len(a)):
        if a[i]!='NA':
            r.append(float(a[i]))
        else:
            r.append(float(0))
    return r
        

#Make a plot of the distribution
def make_plot(ar,id_s):
    N = len(ar)
    a = convert_react(ar)
    w = 1
    ind = np.arange(N)
    fig = figure()
    fig, ax = subplots()
    ax.bar(ind, a, width = w, color = 'r',edgecolor = 'r')
    ax.set_ylabel('DMS Reactivity')
    ax.set_xlabel('Nucleotide')
    ax.set_title(id_s+" reactivity distribution")
    savefig(id_s+'.tif')



    
    
    


