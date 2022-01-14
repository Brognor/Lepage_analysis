import Numerov as nu       
import numpy as np   

f=open("results/data.txt","w")




value =-1.2
radius=2
data = []

f.write(f"Data are generated with a square well potential with value = {value} and radius = {radius}\n\n")

for n in range(1,11):

    grid, step = nu.spatial_objects.uniform_grid(0,50*n,10**7,retstep=True)
    coulomb=nu.potentials.Coulomb(grid)
    well=nu.potentials.radial_well(grid,radius,value)
    potential=coulomb+well
    energy=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-5,0,potential,dy=10**-150,accuracy=10**-10)[0] #data generation
    data.append(energy)



f.write(f"Data are:{data}\n\n")
