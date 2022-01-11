import Numerov as nu
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
f=open("results/dirac_coulomb_result.txt","a")
f.write("We use succession of gaussians as delta approximation chossing as potential -1/r-1.0*delta(r)\n")

grid, step = nu.spatial_objects.uniform_grid(-60,3,3*10**6,retstep=True)
grid=nu.spatial_objects.logarithmic_grid(grid)
coulomb = nu.potentials.Coulomb(grid)

for a in range(1,17,1):
    a=1/a
    
    delta=nu.potentials.dirac_delta_smeared(grid,a)
    potential=coulomb-delta
    energy,eigenfunction= nu.Numerov_algorithm.log_fast_radial_numerov(grid,step,1,1,-50,0,potential,l=0)
    f.write(f"With a={a} in the gaussian succession, we obtain for the first energy eigenvalue with a potential coulomb-delta: {energy}\n")
    eigenfunction=nu.normalize(grid,eigenfunction)
    plt.plot(grid,eigenfunction,label=f'a={a}')
plt.title(r'1S eigenfunction for a potential $-\frac{1}{r}-\delta(r)$')
plt.xlim(0,4.5)
plt.ylim(0,3)
plt.legend(fontsize=9)
plt.savefig("pictures/delta/1S_eigenfunctions")





    
    

