#!/usr/bin/env python

import numpy as np
import param, math
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LinearLocator
from matplotlib.ticker import ScalarFormatter


coords = []
with open('1_INIT_gm', 'r') as f:
    line1 = f.readline().split()
    step = line1[2]
    numclust = line1[5]
    f.readline()
    data = f.readlines()
    for line in data:   
        coords.append([ int(line.split()[0]), 
                      float(line.split()[1]), 
                      float(line.split()[2]), 
                      float(line.split()[3]),
                      float(line.split()[4])])

dist = []
for i in range(len(coords)):
    for j in range(i, len(coords)):
        dist.append(((coords[i][2] - coords[j][2])**2 + 
                     (coords[i][2] - coords[j][3])**2)**0.5)

dist.sort()
print (dist)
