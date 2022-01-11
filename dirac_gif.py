import Numerov as  nu
import matplotlib.pyplot as plt
from tqdm import tqdm
import imageio
import os
grid, step = nu.spatial_objects.uniform_grid(-60,3,10**5,retstep=True)
grid=nu.spatial_objects.logarithmic_grid(grid)
coulomb = nu.potentials.Coulomb(grid)

a=[10./i for i in range(1,100)]
filenames = []
numerfig=1
for a in tqdm(a):
    dirac=-nu.potentials.dirac_delta_smeared(grid, a)
    potential=coulomb+dirac

    energy,eigenfunc=nu.Numerov_algorithm.log_fast_radial_numerov(grid,step,1,1,-100,0,potential)
    eigenfunc=nu.normalize(grid,eigenfunc)

    coulomb_eigenfunc=nu.Numerov_algorithm.log_fast_radial_numerov(grid,step,1,1,-1,0,coulomb)[1]
    coulomb_eigenfunc=nu.normalize(grid,coulomb_eigenfunc)
    
    plt.figure()
    filename=f'pictures/gif/delta/{numerfig}.png'
    filenames.append(filename)
    plt.plot(grid,dirac,label='dirac potential')
    plt.plot(grid, eigenfunc, label='coulomb plus dirac eigenfunction')
    plt.plot(grid, coulomb_eigenfunc, label='coulomb eigenfunction')
    plt.xlim(0,4)
    plt.ylim(-10,3)
    plt.legend()
    plt.title(f"eigenfunction with c=-1 and a= {round(a,3)}")
    plt.savefig(filename)
    plt.close()
    numerfig+=1

images= []

for filename in filenames:
    images.append(imageio.imread(filename))


imageio.mimsave('pictures/gif/delta/dirac_change_n=1_c=-1.gif', images,duration=0.1)

for filename in set(filenames):
    os.remove(filename)


filenames = []
a=[10./i for i in range(1,100)]

numerfig=1000
for a in tqdm(a):
    dirac=-nu.potentials.dirac_delta_smeared(grid, a)
    potential=coulomb+dirac

    energy,eigenfunc=nu.Numerov_algorithm.log_fast_radial_numerov(grid,step,1,1,-100,0,potential)
    eigenfunc=nu.normalize(grid,eigenfunc)

    coulomb_eigenfunc=nu.Numerov_algorithm.log_fast_radial_numerov(grid,step,1,1,-1,0,coulomb)[1]
    coulomb_eigenfunc=nu.normalize(grid,coulomb_eigenfunc)
    
    plt.figure()
    filename=f'pictures/gif/delta/{numerfig}.png'
    filenames.append(filename)
    plt.plot(grid,dirac,label='dirac potential')
    plt.plot(grid, eigenfunc, label='coulomb plus dirac eigenfunction')
    plt.plot(grid, coulomb_eigenfunc, label='coulomb eigenfunction')
    plt.xlim(0,4)
    plt.ylim(-10,3)
    plt.legend()
    plt.savefig(filename)
    plt.close()
    numerfig+=1

images= []

for filename in filenames:
    images.append(imageio.imread(filename))

imageio.mimsave('pictures/gif/delta/dirac_notitle_n=1_c=-1.gif', images,duration=0.01)

for filename in set(filenames):
    os.remove(filename)