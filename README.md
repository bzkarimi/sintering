# sintering

Sintering via Ostwald Ripening written in *Python 3*

Copyright (C) 2020 Borna Zandkarimi

**General comments regarding the code**

This code is a 2D simulation of sintering of nanometer clusters deposited on a surface. Specifically, it has the information needed to simulate Ptn clusters deposited on TiO2(110) surface. It uses the ensemble-average modeling of fluxional clusters (n = 2-8 in this case), i.e. for each cluster size all thermally-accessible isomers, rather than only global minimum structure, are considered at reaction temperature. For clusters larger than that, it is assumed that there is only one isomer per cluster size. Note that this code can be used for all clusters and surfaces. You just need to edit PES and DATA files accordingly. You need to have the PES of a single atom of the cluster deposited on the surface in PES file, and energies corresponding to the clusters of interest in DATA file.

**----------------------------------------------------------------------------------------------------------------------**

**How to Run the code**

**Step 0**: Before running the code, you need to set 4 input files: INIT, DATA, PES, and *param.py*. Note that INIT can be generated manually or automatically using *param.py*.

**Step 1**: First, you need to generate the initial clusters along with their coordinates. This can be done automatically using *autoinit.py* file. You just need to set 3 parameters **num_clust**, **num_single_atom**, and **largest_cluster** in *param.py* and then run *autoinit.py*. The Description of each parameter can be found in *param.py*.

**Step 2**: After running *param.py*, a figure of the initial cell with cluster sizes will be shown. If you are not satisfied with the distribution you can change the parameters and rerun *param.py*. Also, INIT file will be generated. You can manually change the clusters if you wish.

**Note**: *param.py* will automatically suggest a decent supercell size based on the primitive cell parameters and generated clusters. Make sure to modify the **maxx** and **maxy** parameters in *param.py* accordingly.

**Description of each file**

1.INIT: Contains the initial positions of each cluster on the surface along with their radius and energy. Note that the energy should be in hartree.
The file structure is:

NumofAtomsintheCluster      R(A)       X(A)     Y(A)       E(hartree, relative to the most stable cluster)

2.DATA: It's the database which contains all the energies of different cluster sizes. Note that we assume that once the number of atoms in the cluster reaches 9 the larger cluster will grow and the we have only one cluster structure per cluster.

3.PES: potential energy surface (2D) of a single atom on a support.

AtomSymbol    X      Y      E(hartree, relative to the most stable cluster which is Pt8 in this case)

4.param.py: contains the value of input variables such as total number of Metropolis moves, writing step,

5.LOG: contains some messages/warnings/errors.

6.sintering.py: sintering via Ostwald Ripening code based on metropolis moves.

7.autoinit.py:

8.Cluster-energies
