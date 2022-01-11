# Lepage_analysis

Numerov library contains many classes in order to solve the Schrodinger equation, for what concern the most complicated functions we furnished both a numba version and a classic python version, in these users without numba are still able to use to most important function of the library.

xtest.py are test small programm testing functionalities of Numerov library.

radial_well.py is a program finding the energy for various radius and energy level for a Coulomb plus radial_well potential

gif.py contain the program that create gif of evolution of the eigenfunction with increasing radius of radial well potential

.gif file are the gif file

pictures is a folder containing many graphs

test_results is a folder containg results of text.py file

results is the folder containg results obtained for the project

delta_coulomb.py contain the test of the potential -1/r-delta(r), as delta function we used a smeared delta function, what we have done is to use different value of a in the gaussain representation of the delta and see the 1S resulting energy eigenvalue.
