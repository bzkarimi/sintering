#!/usr/bin/env python

import numpy as np


def radius(num_atom):
    a = 2.0975
    b = -0.1963
    r = a * np.log(num_atom) + b
    return r

def E_cluster(num_atom):
    a = -5.28370344 
    b = 2.37917432
    c = 0.65557645
    E_TiO2 = -1723.828749
    E_Pt8_TiO2 = -1761.229018
    Etot = (num_atom * (a + b / num_atom**c) + E_TiO2 - E_Pt8_TiO2) / 27.2114
    return Etot

with open('DATA','a') as f:
    for i in range(101,1001):
        r = radius(i)
        E = E_cluster(i)
        f.write('%-9i  %-5.8f  %14.8f  %18.8f  %14.8f\n' % (i, r, 0.00000000, 0.00000000, E))
