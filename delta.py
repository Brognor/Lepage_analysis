import Numerov as nu
import matplotlib.pyplot as plt

grid, step = nu.spatial_objects.uniform_grid(-60,3,3*10**6,retstep=True)
grid=nu.spatial_objects.logarithmic_grid(grid)
plt.rcParams['text.usetex'] = True

for a in range(1,17,1):
    a=1/a
    
    delta=nu.potentials.smeared_dirac_delta(grid,a)
    plt.plot(grid,delta,label=f'a={a}')

plt.xlim(0,0.4)
plt.title(r'$\delta(r)$ representations in function of a')
plt.legend()
plt.savefig("pictures/delta/delta.png")