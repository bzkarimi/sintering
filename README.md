# sintering
Ostwald sintering code
Copyright (C) 2020 Borna Zandkarimi

To run the code you need 3 input files: INIT, DATA, and PES.

1. INIT: Contains the initial positions of each cluster on the surface along with their radius and energy. Note that the energy should be in hartree.


ElementSymbol   R               X               Y               E
num of atoms in the cluster

2. DATA: It's the database which contains all the energies of different cluster sizes.

3. PES: potential energy surface of a single atom on a support
