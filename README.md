# sintering

Sintering via Ostwald Ripening written in *Python 3*

Copyright (C) 2020 Borna Zandkarimi

##**General comments regarding the code**

This code is a 2D simulation of sintering of nanometer clusters deposited on a surface. Specifically, it has the information needed to simulate Ptn clusters deposited on TiO2(110) surface. It uses the ensemble-average modeling of fluxional clusters (n = 2-8 in this case), i.e. for each cluster size all thermally-accessible isomers, rather than only global minimum structure, are considered at reaction temperature. For clusters larger than that, it is assumed that there is only one isomer per cluster size. Note that this code can be used for all clusters and surfaces. You just need to edit PES and DATA files accordingly. You need to have the PES of a single atom of the cluster deposited on the surface in PES file, and energies corresponding to the clusters of interest in DATA file.

**----------------------------------------------------------------------------------------------------------------------**

**How to run the code**

**Step 0**: Before running the code, you need to set 4 input files: **INIT**, **DATA**, **PES**, and _**param.py**_. Note that **INIT** which contains the initial cluster coordinates, radius, and energies can be generated manually or automatically using _**param.py**_. **DATA** contains cluster size, cluster radius, and cluster energy. **PES** contains the single atom potential energy surface on the support of interest. _**param.py**_ contains the parameters needed to run the simulation. For more detailed description of each file refer to the last section of this file. 

**Step 1**: First, you need to generate the initial clusters along with their coordinates. This can be done automatically using _**autoinit.py**_ file. You just need to set 3 parameters **num_clust**, **num_single_atom**, and **largest_cluster** in *param.py* and then run _**autoinit.py**_. The Description of each parameter can be found in _**param.py**_.

**Step 2**: After running _**param.py**_, a figure of the initial cell with cluster sizes will be shown. If you are not satisfied with the distribution you can change the parameters and rerun _**param.py**_. Also, INIT file will be generated. You can manually change the clusters if you wish.

**Note**: _**param.py**_ will automatically suggest a decent supercell size based on the primitive cell parameters and generated clusters. Make sure to modify the **maxx** and **maxy** parameters in _**param.py**_ accordingly.

**Step 3**: After having **INIT** set up, set all the desired simulation parameters in _**param.py**_. You can find the description of each parameter in _**param.py**_ as a comment line. 

**Step 4**: Now you are ready to run _**sintering.py**_. After the simulation is done, you will get **metropolis** file, which contains information about the cluster distribution in each step and **LOG** file which contains some additional information such as the time of simulation.  

**----------------------------------------------------------------------------------------------------------------------**

**Description of each file (alphabetically)**

_**autoinit.py**_:

**Cluster-energies**:

**DATA**: It's the database which contains all the energies of different cluster sizes. Note that we assume that once the number of atoms in the cluster reaches 9 the larger cluster will grow and the we have only one cluster structure per cluster.

**INIT**: Contains the initial positions of each cluster on the surface along with their radius and energy. Note that the energy should be in hartree.
The file structure is:

NumofAtomsintheCluster      R(A)       X(A)     Y(A)       E(hartree, relative to the most stable cluster)

**LOG**: contains some messages/warnings/errors.

***metropolis**:

_**param.py**_: contains the value of input variables such as total number of Metropolis moves, writing step,

**PES**: potential energy surface (2D) of a single atom on a support.

AtomSymbol    X      Y      E(hartree, relative to the most stable cluster which is Pt8 in this case)

_**sintering.py**_: sintering via Ostwald Ripening code based on metropolis moves.
