#program to test the class auxiliary functions of Numerov library
import Numerov as nu

f=open("test_results/test_auxiliary.txt","w")
func = nu.auxiliary_functions.auxiliary_function_g(1,0,nu.potentials.Coulomb(2.))
f.write(f"The auxiliary function g gfor mass=1, E=0 and the Coulomb potential in r=2 is {func}\n")

func = nu.auxiliary_functions.auxiliary_function_g(1,0,nu.potentials.square_well(2,left=1,right=5,value=12))
f.write(f"The auxiliary function g for mass=1, E=0 and the square well potential in r=2 is {func}\n")


#it is possible to have directly the function vectorized along the grid

grid, step=nu.spatial_objects.uniform_grid(1,10,10,retstep=True)
f.write(f"the grid is {grid}\n")
f.write(f"the step is {step}\n")

func = nu.auxiliary_functions.auxiliary_function_g(1,0,nu.potentials.Coulomb(grid))
f.write(f"The auxiliary function g for mass=1, E=0 and the Coulomb potential along the grid is {func}\n")


func = nu.auxiliary_functions.auxiliary_function_f(func,step)
f.write(f"The auxiliary function f for mass=1, E=0 and the Coulomb potential along the grid is {func}\n")


func = nu.auxiliary_functions.auxiliary_function_g(1,0,nu.potentials.square_well(grid,left=2,right=7,value=4))
f.write(f"The auxiliary function  g for mass=1, E=0 and the square_wellpotential along the grid is {func}\n")

func = nu.auxiliary_functions.auxiliary_function_f(func,step)
f.write(f"The auxiliary function  f for mass=1, E=0 and the square_well potential along the grid is {func}\n")

#let's consider what happen with the logarithm grid
f.write(f"LOGARITHMIC GRID CASE\n")
grid,step=nu.spatial_objects.uniform_grid(0.1,1,10,retstep=True)
f.write(f"now the grid is {grid}\n")
grid=nu.spatial_objects.logarithmic_grid(grid)

f.write(f"the log grid is {grid}\n")
f.write(f"the step is {step}\n")

func = nu.auxiliary_functions.auxiliary_function_g(1,0,nu.potentials.Coulomb(grid))
f.write(f"The auxiliary function g for mass=1, E=0 and the Coulomb potential along the grid is {func}\n")

func = nu.auxiliary_functions.auxiliary_function_f(func,step)
f.write(f"The auxiliary function f for mass=1, E=0 and the Coulomb potential along the grid is {func}\n")

func = nu.auxiliary_functions.auxiliary_function_g(1,0,nu.potentials.square_well(grid,left=2,right=7,value=4))
f.write(f"The auxiliary function  gfor mass=1, E=0 and the sqaure_wellpotential along the grid is {func}\n")

func = nu.auxiliary_functions.auxiliary_function_f(func,step)
f.write(f"The auxiliary function  f for mass=1, E=0 and the square_well potential along the grid is {func}\n")

f.close()


