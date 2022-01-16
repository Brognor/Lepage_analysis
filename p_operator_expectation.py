import Numerov as nu
import numpy as np
import scipy.integrate as int
import matplotlib.pyplot as plt
f =open("results/Operator_analysis.txt","w")
#We compute momentum and eigenfunction in the origin for data
energy_real = []
energy_effective = []
momentum_real = []
momentum_effective = []
origin_real = []
origin_effective = []
for n in range (1,7):
    grid, step = nu.spatial_objects.uniform_grid(0,25*n,10**6,retstep=True)
    coulomb = nu.potentials.Coulomb(grid)
    well = nu.potentials.radial_well(grid,2,-1.2)
    potential = coulomb+well
    energy, eigenfunction = nu.Numerov_algorithm.fast_radial_numerov(grid, step, 1, n ,-2,0,potential, accuracy=10**-7)
    
    eigenfunction = nu.normalize(grid,eigenfunction)
    origin = nu.function_in_origin(grid,eigenfunction) #find value of eigenfunction in the origin
    p_eigenfunc = np.gradient(eigenfunction, step)
    p2_eigenfunc = np.gradient(p_eigenfunc, step)
    p3_eigenfunc = np.gradient(p2_eigenfunc, step)
    p4_eigenfunc = np.gradient(p3_eigenfunc, step) # this is  p^4*eigenfunction
    p_4_mean = int.trapezoid(eigenfunction*p4_eigenfunc,grid)*4*np.pi
    energy_real.append(energy)
    origin_real.append(origin)
    momentum_real.append(p_4_mean)


    #We compute momentum and eigenfunction for effective potential with a=1 with c=-46.27, d=0.28

for n in range (1,7):

    grid, step = nu.spatial_objects.uniform_grid(0,25*n,10**6,retstep=True)
    potential = nu.potentials.effective_potential(grid,1,-46.72,0.28)
    energy, eigenfunction = nu.Numerov_algorithm.fast_radial_numerov(grid, step, 1, n ,-2,0,potential, accuracy=10**-7) 
    eigenfunction = nu.normalize(grid,eigenfunction)
    origin = nu.function_in_origin(grid,eigenfunction) #find value of eigenfunction in the origin
    p_eigenfunc = np.gradient(eigenfunction, step)
    p2_eigenfunc = np.gradient(p_eigenfunc, step)
    p3_eigenfunc = np.gradient(p2_eigenfunc, step)
    p4_eigenfunc = np.gradient(p3_eigenfunc, step) # this is  p^4*eigenfunction
    p_4_mean = int.trapezoid(eigenfunction*p4_eigenfunc,grid)*4*np.pi
    energy_effective.append(energy)
    origin_effective.append(origin)
    momentum_effective.append(p_4_mean)


f.write("ENERGY\n\n")
f.write("Level \t\t Data \t\t Effective\n")
for n in range (1,7):

   
    f.write(f" {n}S \t\t {round(energy_real[n-1],6)} \t {round(energy_effective[n-1],6)}\n")

f.write("\nMOMENTUM TO THE FOURTH\n\n")
f.write("Level \t\t Data \t\t Effective\n")
for n in range (1,7):

    
    f.write(f" {n}S\t \t {round(momentum_real[n-1],6)}\t{round(momentum_effective[n-1],6)}\n")



f.write("\nSQUARED EIGENFUNCTION IN THE ORIGIN\n\n")
f.write("Level \t\t Data \t\t Effective\n")

for n in range (1,7):

    
    f.write(f" {n}S \t\t{round(origin_real[n-1],6)}  \t{round(origin_effective[n-1],6)}\n")

f.close()



