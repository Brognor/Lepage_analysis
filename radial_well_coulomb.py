import Numerov as  nu
import matplotlib.pyplot as plt
from tqdm import tqdm
f=open("results/radial_well_coulomb_energy.txt","w")
grid,step=nu.spatial_objects.uniform_grid(0,100,10**6,retstep=True)
value =-1
radius=[0.1,0.2,0.3,0.5,0.75,1.0,1.5,2.0]
coulomb=nu.potentials.Coulomb(grid)
for radius in tqdm(radius):
    well_potential=nu.potentials.radial_well(grid, radius=radius, value=value)
    potential=coulomb+well_potential
    f.write(f"Here radial well potential has  VALUE={value} and RADIUS=={radius}\n\n")
    for n in range(1,6):
        energy,eigenfunc=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-3,0,potential)
        f.write(f"Energy of the {n}S level is {energy}, only Coulomb would be {-1/(2*n**2)}\n")
        eigenfunc=nu.normalize(grid,eigenfunc)

        if (radius==0.1 or radius==0.5 or radius==1.0 or radius==1.5 or radius==2.0):
            coulomb_eigenfunc=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-1,0,coulomb)[1]
            coulomb_eigenfunc=nu.normalize(grid,coulomb_eigenfunc)
            plt.plot(grid,well_potential,label=' well potential')
            plt.plot(grid, eigenfunc, label='coulomb plus well eigenfunction')
            plt.plot(grid, coulomb_eigenfunc, label='coulomb eigenfunction')
            plt.xlim(0,20*n)
            plt.legend()
            plt.title(f"eigenfunction with n={n} and radius of short interaction = {radius}")
            plt.savefig(f"pictures/radial_well/eigenfunc_n={n}_radius={radius}.png")
            plt.close()
        
    f.write("\n#############################################################\n")
