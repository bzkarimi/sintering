#!/usr/bin/env python

import numpy as np
import constants

def distance(x1, y1, x2, y2):

  r = ((x1 - x2)**2.0 + (y1 - y2)**2.0)**(0.5)

  return r

with open('par_local.xyz.0','r') as f1: 
  num_lines    = sum(1 for line in open('par_local.xyz.0'))
  final_data   = []
  E_data       = []
  
  for i in range(2*constants.NClust_tot): 

    line = f1.readline()
    if ( (np.shape(line.split())) == (2,) ):    
      natoms = int(line.split()[0])
      E_data.append(float(line.split()[1]))

    else:
      raw_data         = []
      unsorted_circles = []

      for j in range(0, natoms):
        raw_data.append(line.split())

        if (j != natoms-1):
          line = f1.readline()

      for j in range(0, natoms):
        for k in range(j+1, natoms):
          r       = distance(float(raw_data[j][1]), float(raw_data[j][2]), float(raw_data[k][1]), float(raw_data[k][2]))
          centerx = ( float(raw_data[j][1]) + float(raw_data[k][1]) )/2.0
          centery = ( float(raw_data[j][2]) + float(raw_data[k][2]) )/2.0
          unsorted_circles.append([r,centerx,centery,natoms])

      final_data.append(sorted(unsorted_circles, key=lambda circle: circle[0])[-1])    #sort circles by R

for i in range(len(final_data)):
  final_data[i].append(E_data[i])
  final_data[i][0] += constants.Ratom

print final_data

with open('DATA','w') as f2:
  f2.write('%-3s  %14s  %14s  %14s  %14s\n' % ('Pt','R','X','Y', 'E'))

  for i in range(len(final_data)):
    f2.write('%3i  %16.8f  %16.8f  %16.8f  %16.8f\n' % ( int(final_data[i][3]), final_data[i][0], final_data[i][1], final_data[i][2], float(final_data[i][4])  ))
