import Numerov as  nu
import matplotlib.pyplot as plt
from tqdm import tqdm
import imageio
f=open("test_results/radial_well_coulomb_energy","w")
grid,step=nu.spatial_objects.uniform_grid(0,100,10**4,retstep=True)
value =-1
coulomb=coulomb=nu.potentials.Coulomb(grid)
radius=[0.05*i for i in range(0,2000)]
filenames = []
numerfig=1000
for radius in tqdm(radius):
    well_potential=nu.potentials.radial_well(grid, radius=radius, value=value)
    potential=coulomb+well_potential

    energy,eigenfunc=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,5,-2,0,potential)
    eigenfunc=nu.normalize(grid,eigenfunc)

    coulomb_eigenfunc=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,5,-1,0,coulomb)[1]
    coulomb_eigenfunc=nu.normalize(grid,coulomb_eigenfunc)
    
    plt.figure()
    filename=f'figures/gif/{numerfig}.png'
    filenames.append(filename)
    plt.plot(grid,well_potential,label=' well potential')
    plt.plot(grid, eigenfunc, label='coulomb plus well eigenfunction')
    plt.plot(grid, coulomb_eigenfunc, label='coulomb eigenfunction')
    plt.xlim(0,100)
    plt.legend()
    plt.title(f"eigenfunction with radius of short interaction = {round(radius,2)}")
    plt.savefig(filename)
    plt.close()
    numerfig+=1

images= []

for filename in filenames:
    images.append(imageio.imread(filename))

imageio.mimsave('figures/gif/radialweel_change5.gif', images,duration=0.1)

filenames = []
radius=[0.05*i for i in range(0,2000)]

numerfig=5000
for radius in tqdm(radius):
    well_potential=nu.potentials.radial_well(grid, radius=radius, value=value)
    potential=coulomb+well_potential

    energy,eigenfunc=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,5,-2,0,potential)
    eigenfunc=nu.normalize(grid,eigenfunc)

    coulomb_eigenfunc=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,5,-1,0,coulomb)[1]
    coulomb_eigenfunc=nu.normalize(grid,coulomb_eigenfunc)
    
    plt.figure()
    filename=f'figures/gif/{numerfig}.png'
    filenames.append(filename)
    plt.plot(grid,well_potential,label=' well potential')
    plt.plot(grid, eigenfunc, label='coulomb plus well eigenfunction')
    plt.plot(grid, coulomb_eigenfunc, label='coulomb eigenfunction')
    plt.xlim(0,100)
    plt.legend()
  
    plt.savefig(filename)
    plt.close()
    numerfig+=1

images= []
for filename in filenames:
    images.append(imageio.imread(filename))

imageio.mimsave('figures/gif/radial_well_notitle5.gif', images,duration=0.01)