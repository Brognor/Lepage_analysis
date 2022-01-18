# Lepage_analysis

Numerov library contains many classes in order to solve the Schrodinger equation, for what concern the most complicated functions we furnished both a numba version and a classic python version, in these users without numba are still able to use to most important function of the library.

xtest.py are test small programm testing functionalities of Numerov library.

radial_well.py is a program finding the energy for various radius and energy level for a Coulomb plus radial_well potential

gif.py contain the program that create gif of evolution of the eigenfunction with increasing radius of radial well potential


pictures is a folder containing many graph, in particular:
	pictures of wavefunctions for diffrent potentials
	potentials graphs
	final resullt of relative energy errors graph
	gif 

test_results is a folder containg results of text.py file

results is the folder containg results obtained for the project

delta_coulomb.py contain the test of the potential -1/r-delta(r), as delta function we used a smeared delta function, what we have done is to use different value of a in the gaussain representation of the delta and see the 1S resulting energy eigenvalue.

delta.py plot many representations of our delta function with different a

delta_1st_order fit find results for a dirac delta plus coulomb potential approxiamted at first order, result of this are in results

c_delta_with_a find results for a regulated coulomb plus a smeared delta potential fitting the best value of c, results are found in results

effective_fit perform the search for best values of c and d that minimize square error

effective_low perfeorm the search for best value of c and d that makes low energy matching
<<<<<<< HEAD
=======

>>>>>>> 79fbb0684134a63b752a7cd7f61709852494256f

operator_expectation compare energy, squared eigenfunction in the origin and expected momentum to the fourth for data and effective case

