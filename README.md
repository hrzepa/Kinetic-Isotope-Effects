# Kinetic-Isotope-Effects
A program to solve the  Bigeleisen-Mayer equations for kinetic isotope effects.

Calculation of harmonic isotope effects based on the normal vibrational wavenumbers of Reactant and an isotopomer, and  transition state and the equivalent isotopomer.

Original version: H S Rzepa 1975, Austin Texas.
Literature references: doi: 10.1021/ja00493a008  and doi: 10.1021/ja00486a013 (1978).
Formatting and I/O changes: H S Rzepa, 1980, Imperial College.
Comments:  H S Rzepa, July, 2015, Imperial College.

Compilation:   gfortran kinisot.f  -o kinisot.exe (Tested 3 July, 2015 on OS X 10.10.4).  Compiler from http://hpc.sourceforge.net/   Version: gcc version 5.1.0 (GCC) 

Inputs:  isotope.dat.  Output: isotope.out

Format of Input file
1:  48 288 298 300 310 320 330 340 350 360 370     where N=48 is the number of atoms  followed by  10 temperatures. If  number of atoms is zero, program stops.
2: Title for normal isotope reactant
3: 3N-6 values for the normal mode wavenumbers for reactant.
4: Title for isotopomer of reactant
5: 3N-6 values for the normal mode wavenumbers for isotopomeric reactant
6: Title for normal isotope transition state
7: 3N-6 values for the normal mode wavenumbers for transition state, with imaginary mode  listed first (as a -ve number)
8: 3N-6 values for the normal mode wavenumbers for isotopomeric transition state, with imaginary mode  listed first (as a -ve number)
9: Repeat card 1 for new isotopomers.



