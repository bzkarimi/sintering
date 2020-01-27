#!/usr/bin/env python

'''
Sintering simulation of metalic clusters on a surface based on NVT ensemble metropolis moves.

Borna Zandkarimi 2020

'''

import numpy as np
import timeit, math, copy
import constants
from scipy.special import expit          # for handling very small exp


#|================|
#|   FUNCTIONS    |
#|================|

def distance(x1, y1, x2, y2):
    r = ((x1 - x2)**2.0 + (y1 - y2)**2.0)**(0.5)
    return r

def PES_finder(X_new, Y_new, PES_copy, N_mesh):
    finder = False
    for ii in range(N_mesh): 
        dx =  np.abs( constants.xstep_max*( math.ceil(X_new/constants.xstep_max) % N_meshx ) - PES_copy[ii][1] )
        dy =  np.abs( constants.ystep_max*( math.ceil(Y_new/constants.ystep_max) % N_meshy ) - PES_copy[ii][2] )
        if ( dx  < 0.1 and dy  < 0.1 ):    
            # PES found 
            E_PES  = PES_copy[ii][3] 
            finder = True
            break    
    if (finder == False):
        raise ValueError('Could not find PES for single atom. Bad move! Check you initial setup')      
    return E_PES

def Cluster_finder(Clusters, OUTPUT_data, indx, which='both'):
     prob1     = []
     prob2     = []
     prob1_tmp = []
     prob2_tmp = []
     for j in range(Ncluster_tot):     
         if (Clusters[j][0] == OUTPUT_data[indx][0] + 1 and (which == 'add' or which == 'both') ) :
             if (prob1_tmp == []):
                 E_min_add = Clusters[j][4]
                 prob1_tmp.append(1.0)     
             else:
                 # choose a larger cluster
                 prob1_tmp.append( expit(-beta * (Clusters[j][4] - E_min_add) ))     

         if (Clusters[j][0] == OUTPUT_data[indx][0] - 1 and (which == 'remove' or which == 'both') ) :
             if (prob2_tmp == []):
                 E_min_rmv = Clusters[j][4]
                 prob2_tmp.append(1.0)     
             else:
                 # choose a smaller cluster
                 prob2_tmp.append( expit(-beta * (Clusters[j][4] - E_min_rmv) ))    

     if ( which == 'add' or which == 'both' ): 
         prob1 = [icount / sum(prob1_tmp) for icount in prob1_tmp]
         # weighted random chosen number
         clust_idx_add = np.random.choice(np.arange(len(prob1)), 1, p=prob1, replace=False)     
         Enew_ad   = - np.log(prob1_tmp[clust_idx_add[0]]) / beta + E_min_add

     if ( which == 'remove' or which == 'both' ): 
         prob2 = [icount / sum(prob2_tmp) for icount in prob2_tmp]
         # weighted random chosen number
         clust_idx_rmv = np.random.choice(np.arange(len(prob2)), 1, p=prob2, replace=False)     
         Enew_remv = - np.log(prob2_tmp[clust_idx_rmv[0]]) / beta + E_min_rmv

           
     remove = False
     add    = False
     for k in range(Ncluster_tot):          
       if ( (Clusters[k][0] == OUTPUT_data[indx][0] - 1) and (remove==False) and (which == 'remove' or which == 'both') ):
           Rnew_remv = Clusters[k+clust_idx_rmv[0]][1]
           remove    = True            
       if ( (Clusters[k][0] == OUTPUT_data[i][0] + 1) and (add==False) and (which == 'add' or which == 'both') ):
           Rnew_ad  = Clusters[k+clust_idx_add[0]][1]
           add      = True 
 
     # return new energy and radius

     if (which == 'both'):
         return Enew_ad, Enew_remv, Rnew_ad, Rnew_remv
     elif (which == 'add'):
         return Enew_ad, Rnew_ad
     elif (which == 'remove'):
         return Enew_remv, Rnew_remv
 
#|================|
#|   READ INPUT   |
#|================|

PES           = []                        # potential energy surface element, x, y, z, E
INIT_data     = []                        # initial cluster R, x, y, E
Clusters      = []                        # all possible cluster R and E
beta          = constants.beta300         # Choose Temperature
hartree_to_ev = constants.hartree_to_ev   
Ncluster      = constants.NClust          # Number of clusters on the surface during simulation which can be changed
Ncluster_tot  = constants.NClust_tot      # Unique different clusters (is constant in the input file)
Metro_Max     = constants.MMAX            # num of Metropolis steps
write_step    = constants.wstep           # writing metropolis each wstep
N_mesh        = constants.N_mesh          # total number of PES points
N_meshx       = N_meshy = (N_mesh)**(0.5) # number of mesh points in x and y direction

with open('INIT','r') as f1:
    f1.readline()
    data = f1.readlines()
    for line in data:  
        INIT_data.append([ int(line.split()[0]), 
                           float(line.split()[1]), 
                           float(line.split()[2]), 
                           float(line.split()[3]),
                           float(line.split()[4]) ])

with open('DATA','r') as f2:
    f2.readline()
    data = f2.readlines()
    for line in data:  
        Clusters.append([ int(line.split()[0]), 
                          float(line.split()[1]), 
                          float(line.split()[2]), 
                          float(line.split()[3]),
                          float(line.split()[4]) ])

with open('PES','r') as f3:
    data = f3.readlines()
    for line in data:   
        PES.append([ str(line.split()[0]), 
                     float(line.split()[1]), 
                     float(line.split()[2]), 
                     float(line.split()[3]) ])

#|===================|
#|  METROPOLIS LOOP  |
#|===================|

np.random.seed()       # seed for random number generator
start_time_MC   = timeit.default_timer()

OUTPUT_data  = copy.deepcopy(INIT_data)
PES_copy     = copy.deepcopy(PES)

with open('metropolis','w') as f4:

  for step in range(Metro_Max+1):

    #|================|
    #|  write output  |
    #|================|

    if ( (step % write_step) == 0 ): 

      f4.write('%-15s  %12i\n' % ('step',step))
      f4.write('%4s  %14s  %14s  %14s  %14s\n' % 
               ('Pt', 'R', 'X','Y','E'))

      for i in range(Ncluster):

        f4.write('%3i  %16.8f  %16.8f  %16.8f  %16.8f\n' % 
                 (OUTPUT_data[i][0], OUTPUT_data[i][1], OUTPUT_data[i][2], OUTPUT_data[i][3], OUTPUT_data[i][4]))

      f4.write('\n')

    #|======================|
    #|  choosing a cluster  |
    #|======================|

    indx = int(np.random.rand() * Ncluster) 
    #indx = 1

    num_atm_temp  = OUTPUT_data[indx][0] 
    R_temp        = OUTPUT_data[indx][1]
    X_temp        = OUTPUT_data[indx][2]
    Y_temp        = OUTPUT_data[indx][3]
    E_temp        = OUTPUT_data[indx][4]

    #|======================|
    #|  choosing step size  |
    #|======================|

    if ( num_atm_temp == 1 ): #For monomer case.  

      px = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < px < 4   
      py = int( 4.0*(2.0*np.random.rand() - 1.0) )   # -4 < py < 4   
                    
      print ('px,py=',px, py)
    else:  #If part of cluster we want bigger step, otherwise never breaks. We ensure that 
         
      px = math.ceil( (9.0 * R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) #the step-size is cluster size dependent by generalizing 
      py = math.ceil( (9.0 * R_temp)  * ( 2.0*np.random.rand() - 1.0 ) ) #the problem of packing circles in square, circle, & chain.
    

    #print 'px,py',px, py
    #print 'xtemp=',X_temp
    #print 'ytemp=',Y_temp
    #print 'xstep=', constants.xstep_max
    #print 'ystep=', constants.ystep_max


    X_new = X_temp + px * constants.xstep_max
    Y_new = Y_temp + py * constants.ystep_max

    #print 'X,Y', X_new, Y_new
 
    if ( X_new == X_temp and Y_new == Y_temp ):
      continue
    #|==============================|
    #| Periodic Boundry Conditions  |
    #|==============================|
    if ( X_new > constants.maxx ):   
      X_new = constants.minx - constants.maxx + X_new
    elif ( X_new < constants.minx ):   
      X_new = constants.maxx - constants.minx + X_new
    if ( Y_new > constants.maxy ):  
      Y_new = constants.miny - constants.maxy + Y_new
    elif ( Y_new < constants.miny ): 
      Y_new = constants.maxy - constants.miny + Y_new
    #|==============================|
    #|   Check for New Clusters     |
    #|==============================|

    inside_cluster = False
    for i in range (Ncluster):
   
      #print '==============================='
      #print 'distance=', distance(X_new, Y_new, OUTPUT_data[i][2]), OUTPUT_data[i][3] )
      #print 'R + Ratom', OUTPUT_data[i][1] + constants.Ratom
      #print 'Xnew,Ynew', X_new, Y_new

      if ( (distance(X_new, Y_new, OUTPUT_data[i][2], OUTPUT_data[i][3] ) <  (OUTPUT_data[i][1] + constants.Ratom)) and (i !=indx) ):   # if it is inside a cluster

        #print 'YES inside a new cluster!'

        Eold = E_temp + OUTPUT_data[i][4]       # total Eold

        if ( OUTPUT_data[indx][0] == 1):
      
          Enew_rmv = 0.0

          Enew_add, Rnew_add = Cluster_finder(Clusters, OUTPUT_data, indx, 'add')

        elif ( OUTPUT_data[indx][0] == 2):

          # CLUSTER FINDER FOR EADD
          Enew_add, Rnew_add = Cluster_finder(Clusters, OUTPUT_data, indx, 'add')

          Enew_rmv = PES_finder(X_new, Y_new, PES_copy, N_mesh)
          Rnew_rmv = constants.Ratom 


        #|======================|
        #| Choose a new cluster |
        #|======================|

        else:

          Enew_add, Enew_rmv, Rnew_add, Rnew_rmv = Cluster_finder(Clusters, OUTPUT_data, indx)

        Enew = Enew_add + Enew_rmv

        #print 'Enew_add=', Enew_add
        #print 'Enew_rmv=', Enew_rmv

        #print 'E1_OLD=', E_temp 
        #print 'E2_OLD=', OUTPUT_data[i][4]       

        #print 'delta E=', Enew-Eold
        #print 'prob=', np.exp( -beta*(Enew-Eold)*hartree_to_ev )
        #print 'Rnewrmv', Rnew_rmv     
        #print 'Rnewadd', Rnew_add 

        if ( (Enew-Eold < 0.0) or ( np.exp( -beta*(Enew-Eold) ) > np.random.rand()) ):  # accept the move
          
          #|========================|
          #|   Modify old clusters  |
          #|========================|

          OUTPUT_data[i][0] = OUTPUT_data[i][0] + 1           # Natom
          OUTPUT_data[i][1] = Rnew_add                           # R
          OUTPUT_data[i][4] = Enew_add                           # E
  
          if ( OUTPUT_data[indx][0] == 1 ):        #Natom = 1

            # DELETE A CLUSTER (single atom)
            del (OUTPUT_data[indx])
            Ncluster -= 1

          else:

            OUTPUT_data[indx][0] = OUTPUT_data[indx][0] - 1     # Natom
            OUTPUT_data[indx][1] = Rnew_rmv                        # R
            OUTPUT_data[indx][4] = Enew_rmv                        # E

        else:
          print('MOVE TO A NEW CLUSTER IS NOT FAVORABLE!')

        inside_cluster = True
        break  
 
    if( inside_cluster == False ):  # we get a new single atom meaing that a new cluster forms
      
      print 'inside cluster = false'
      if ( OUTPUT_data[indx][0] == 1):
        Enew_rmv = 0.0

      elif ( OUTPUT_data[indx][0] == 2):
        Enew_rmv = PES_finder(X_temp, Y_temp, PES_copy, N_mesh)
        Rnew_rmv = constants.Ratom

        #|========================|
        #|  Choose a new cluster  |
        #|========================|
      
      else:
        Enew_rmv, Rnew_rmv = Cluster_finder(Clusters, OUTPUT_data, indx, 'remove')

      E_atom = PES_finder(X_new, Y_new, PES_copy, N_mesh)
      Enew   = E_atom + Enew_rmv 
      Eold   = OUTPUT_data[indx][4]

      if ( (Enew-Eold < 0.0) or ( np.exp( -beta*(Enew-Eold) ) > np.random.rand()) ):  # accept the move

        print('SINGLE ATOM MOVE OUTSIDE OF A CLLUSTER ACCEPTED!')
        #print 'prob look here=', np.exp( -beta*(Enew-Eold))
        #print 'Eold=', Eold
        #print 'Enew=', Enew
 
        #break
        #|========================|
        #|   Modify old clusters  |
        #|========================|

        if ( OUTPUT_data[indx][0] == 1 ):     # Note for a single atom move on surface

          del (OUTPUT_data[indx])
          Ncluster -= 1

        else:

          OUTPUT_data[indx][0] = OUTPUT_data[indx][0] - 1     # Natom
          OUTPUT_data[indx][1] = Rnew_rmv                        # R
          OUTPUT_data[indx][4] = Enew_rmv                        # E

        OUTPUT_data.append([1, constants.Ratom, X_new, Y_new, E_atom])             
        Ncluster += 1
       
      else:
        print('SINGLE ATOM MOVE OUTSIDE OF A CLLUSTER NOT ACCEPTED!')
      
print('final number of clusters=', Ncluster)
print('MC time in seconds =', timeit.default_timer() - start_time_MC)
