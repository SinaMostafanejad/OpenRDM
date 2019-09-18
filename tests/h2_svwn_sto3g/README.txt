This folder contains the data files:

1- eref.txt:     The reference electronic energy for H2 molecule
                 calculated using SVWN3/STO-3G level of theory.

2- eclass.txt:   The classical energy which is sum of nuclear repulsion,
                 nuclear attraction, kinetic energy and Hartree Coulomb
                 interaction energies.

3- cmat.txt:     The AO->MO transformation matrix C.
                 The dimensions <#rows,#cols> is shown in the first 
                 line of the file.

4- grids.txt:    The Cartesian grid points with the order w,x,y,z
                 with the number of grids on the first line.

5- orbitals.txt  The super phi(p,i) matrix hosts the value of the
                 i-th basis function calculated at each point p.
                 The first line shows the number of grid points and
                 number of basis functions, respectively. 

6- opdm.txt:     The 1-electron reduced-density matrix (1RDM) in the
                 AO basis. The dimensions <#rows,#cols> is shown in 
                 the first line of the file.
