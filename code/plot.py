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

print('step=',step)
print('numclust=',numclust)
# plot setting

minorLocator = AutoMinorLocator()
mlx  = MultipleLocator(param.xstep_max)
mly  = MultipleLocator(param.ystep_max)

title_font = {'fontname':'Times New Roman', 'size':'18', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Times New Roman', 'size':'18'}
mpl.rc('font',family='Times New Roman')

ax = plt.subplot() # Defines ax variable by creating an empty plot

# Set the tick labels font
for label in (ax.get_yticklabels() + ax.get_xticklabels()):
    label.set_fontname('Times New Roman')
    label.set_fontsize(16)

#ax.yaxis.set_minor_locator(mly)
#ax.xaxis.set_minor_locator(mlx)

xcoords = []
ycoords = []
Rcoords = []
types = []
for i in range(len(coords)):
    types.append(coords[i][0])
    Rcoords.append(30*coords[i][1]**2)
    xcoords.append(coords[i][2])
    ycoords.append(coords[i][3])

frame1 = plt.gca()
frame1.axes.yaxis.set_ticklabels([])
plt.yticks([])
frame1.axes.xaxis.set_ticklabels([])
plt.xticks([])
plt.ylim(param.miny,param.maxy)
plt.xlim(param.minx,param.maxx)
for i,type in enumerate(types):
    x = xcoords[i]
    y = ycoords[i]
    plt.scatter(xcoords, ycoords, color='blue', s=Rcoords, marker='o')
    plt.text(x, y, type, fontsize=(12+type/3), color='yellow', horizontalalignment='center', 
             verticalalignment='center')
#plt.axvline(x=xdups*param.primcell_a, color='black', linestyle='-')
#plt.axhline(y=ydups*param.primcell_b, color='black', linestyle='-')
#plt.grid(which='minor', axis='both', linestyle='--')
plt.show()
