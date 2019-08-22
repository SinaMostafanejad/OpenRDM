This folder contains the data files:

1- eref.txt:     The reference electronic energy for H2 molecule
                 calculated using SVWN3/STO-3G level of theory.

2- eclass.txt:   The classical energy which is sum of nuclear repulsion,
                 nuclear attraction, kinetic energy and Hartree Coulomb
                 interaction energies.

3- cmat.txt:     The AO->MO transformation matrix C.
                 The dimensions <#rows,#cols> is shown in the first 
                 line of the file.

4- grids.txt:    The Cartesian grid points with the order w,x,y,z,
                 and phi(r), the value of the orbitals on those points.
                 The number of grid points is in the first line.
