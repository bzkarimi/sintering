#!/usr/bin/env python

'''
================    ===============================
    Variable                  Description
================    ===============================
    MMAX             num of Metropolis moves
    wstep            write metropolis every wstep
    Tot              total num of data points
    maxx             maximum X
    maxy             maximum Y
    minx             minimum X
    miny             minimum Y
    Nptmax           max Pt atom in cluster
    Nznmax           max Zn atom in cluster
    Maxin            num of atoms undergoing Metropolis moves
                     = Nptmax + Nznmax 
    N_mesh           total number of PES points
----------------     -------------------------------
    beta300          Boltzmann factor at 300 K  (ev)
    beta700          Boltzmann factor at 700 K  (ev)
    beta1000         Boltzmann factor at 1000 K (ev)
    gamma300         surface tension of liquid metal at 300 K 
    gamma700         surface tension of liquid metal at 700 K
    gamma1000        surface tension of liquid metal at 1000 K
                     Lida98 & L. Aguilera MS

    ozzy300          Scaling factor of Ost.-Fre. Eq. 
    ozzy700          ( 8 * Pi * gamma * (rave^2) ) / 3
    ozzy1000         for different T = 300,700,1000 K  
  
    lamd300          de Broglie lambda for Zn atom 
    lamd700          at given T for evap-model
    lamd1000         300,700,1000 K
    rave     ????????             average coverage radius/ratio for Pt-Zn 1:3
    rbon            average height of monomer, r_bond
    consph   ????????????    
    xstep    ???????????
    ystep    ??????????
================    ===============================
'''
 
MMAX   =  100000
Tot    =  51200

#maxx   =  26.43   
#minx   = -26.43

#maxy   =  24.78 
#miny   = -24.19 


unitcell_a = 13.216 
unitcell_b = 12.243
unitcell_c = 26.362

primcell_a = 6.607835831
primcell_b = 3.060749557

maxx   =  primcell_a * 6   
minx   =  0.0

maxy   =  primcell_b * 12 
miny   =  0.0

Nptmax =  4
Nznmax =  4 
Maxin  =  8
NClust =  4
NClust_tot = 33
Ratom  =  1.35
wstep  =  1 
N_mesh = 11*11

xstep_max  =  primcell_a / 10.0
ystep_max  =  primcell_b / 10.0

hartree_to_ev = 27.2114
k_J       = 1.38064852e-23
k_ev      = 8.6173303e-5
k_hartree = k_ev/hartree_to_ev

kT300   = k_hartree*300.0
beta300 = 1/kT300

#beta300       = 38.663660 
#beta700       = 16.576243
#beta1000      = 11.604516

gamma300  = 0.052973 
gamma700  = 0.048729
gamma1000 = 0.045545   
ozzy300   = 0.83293830
ozzy700   = 0.76620238
ozzy1000  = 0.71615043
lamd300   = 0.2209 
lamd700   = 0.1446
lamd10000 = 0.1210

rave      = 1.35 
rbon      = 3.25 
consph    = 4.188790

