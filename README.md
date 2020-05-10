# Sintering

Sintering via Ostwald Ripening written in *Python 3*

Copyright &copy; 2020 Borna Zandkarimi

## **General comments regarding the code**:

This code is a 2D simulation of sintering of nanometer clusters deposited on the surface. Specifically, it has the information needed to simulate Pt clusters deposited on rutile TiO<sub>2</sub>(110) surface. It uses the ensemble-average modeling of fluxional clusters (n = 2-8 in this case), i.e. for each cluster size all thermally-accessible isomers, rather than only global minimum structure, are considered at reaction temperature. For clusters larger than that, it is assumed that there is only one isomer per cluster size. Note that this code can be used for all clusters and surfaces. You just need to edit **PES** and **DATA** files accordingly. You need to have the PES of a single atom of the cluster deposited on the surface in **PES** file, and energies corresponding to the clusters of interest in **DATA** file.

## **How to run the code**:

* **Step 0**: Before running the code, you need to set 4 input files: **INIT**, **DATA**, **PES**, and _**param.py**_. Note that **INIT** which contains the initial cluster coordinates, radius, and energies can be generated manually or automatically using _**autoinit.py**_. **DATA** contains cluster size, cluster radius, and cluster energy. **PES** contains the single atom potential energy surface on the support of interest. _**param.py**_ contains the parameters needed to run the simulation. For more detailed description of each file refer to the last section of this file.

* **Step 1**: First, you need to generate the initial clusters along with their coordinates. This can be done automatically using _**autoinit.py**_ file. You just need to set 3 parameters **num_clust**, **num_single_atom**, and **largest_cluster** in _**param.py**_ and then run _**autoinit.py**_. The Description of each parameter can be found in _**param.py**_.

* **Step 2**: After running _**param.py**_, a figure of the initial cell with cluster sizes will be shown. If you are not satisfied with the distribution you can change the parameters and rerun _**autoinit.py**_. Also, after **INIT** is generated you can manually change the clusters' information if you wish.

  **Note**: _**autoinit.py**_ will automatically suggest a decent supercell size based on the primitive cell parameters and generated clusters in the output terminal. Make sure to use that to modify the **maxx** and **maxy** parameters in _**param.py**_ accordingly.

* **Step 3**: After having **INIT** set up, set all the desired simulation parameters in _**param.py**_. You can find the description of each parameter as comment line inside _**param.py**_. 

* **Step 4**: Now you are ready to run _**sintering.py**_. After the simulation is done, you will get **metropolis** file, which contains information about the cluster size and energy distribution in each step and **LOG** file which contains some additional information such as the time of simulation.  

## **Description of each file (alphabetically)**:

* _**autoinit.py**_: Automatically generates the initial cluster size, radius, energy and coordinates based on the 3 parameters **num_clust**, **num_single_atom**, and **largest_cluster** in _**param.py**_.

* **Cluster-energies**: Contains the energies of Ptn clusters (n = 2-8) in hartree. This file is not necessary for simulation.

* **DATA**: This file is actually the database which contains all the energies and radii of different cluster sizes up to n = 1000. Note that we assume once the number of atoms in the cluster reaches 9 we have only one cluster structure per cluster. The file structure is:

  Num-of-Atoms-in-the-Cluster &nbsp; &nbsp; R(&angst;) &nbsp; &nbsp; X(&angst;) &nbsp; &nbsp; Y(&angst;) &nbsp; &nbsp; E(hartree, relative to the most stable cluster)

* **INIT**: Contains the initial positions of each cluster on the surface along with their corresponding radius and energy. **Note**: that the energy should be in hartree. The file structure is:

  Num-of-Atoms-in-the-Cluster &nbsp; &nbsp; R(&angst;) &nbsp; &nbsp; X(&angst;) &nbsp; &nbsp; Y(&angst;) &nbsp; &nbsp; E(hartree, relative to the most stable cluster)

* **LOG**: Contains information about the simulation time and it might also contain some messages regarding metropolis move rejections.

* **metropolis**: This file is basically the trajectory file which contains the information about the system in every step. (or every **wstep** which is set in _**param.py**_)

* _**param.py**_: Contains the value of input parameters for generating initial clusters and for simulation such as total number of Metropolis moves, writing step, temperature, etc. Detailed description of each parameter can be found as comment line inside the file.

* **PES**: Potential energy surface (2D) of a single atom on the support obtained from PAW-DFT. In this case it is Pt/TiO<sub>2</sub>(110). The file structure is:

  Atom-Symbol &nbsp; &nbsp; X(&angst;) &nbsp; &nbsp; Y(&angst;) &nbsp; &nbsp; E(hartree, relative to the most stable cluster)

* _**sintering.py**_: The main code of sintering via Ostwald Ripening based on metropolis moves and NVT ensemble.
