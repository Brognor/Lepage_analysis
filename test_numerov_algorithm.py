import Numerov as nu
import time
import matplotlib.pyplot as plt
grid,step=nu.spatial_objects.uniform_grid(0,100,10**6,retstep=True)

f=open("test_results/numerov_test.txt","w")
coulomb=nu.potentials.Coulomb(grid)


start=time.time()

nu.Numerov_algorithm.radial_numerov(grid, step, 1, 1, -2, 0, coulomb)
end=time.time()
f.write(f"For the radial algorithm, radius from 0 to 100, step of 10^-4, the time of execution is {end-start} seconds\n")

start=time.time()
nu.Numerov_algorithm.fast_radial_numerov(grid, step, 1, 1, -2, 0, coulomb)
end=time.time()
f.write(f"For the fast radial algorithm, radius from 0 to 100, step of 10^-4 the time of execution is {end-start} seconds\n")

coulomb=nu.potentials.Coulomb(grid)

start=time.time()
for n in range(1,6):
    energy=nu.Numerov_algorithm.fast_radial_numerov(grid, step, 1, n, -1., 0, coulomb)[0]
    f.write(f"The energy eigenvalue for the level S with n={n} is {energy}, the theory provide {-1/(2*n**2)}\n")
end=time.time()

f.write(f"For the fast radial algorithm, radius from 0 to 100, step of 10^-4, the time of execution in order to find all these values is  {end-start} seconds\n")


func=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,3,-1,0,coulomb)[1]

func=nu.normalize(grid,func)

plt.plot(grid,func)
plt.title("Normalized 3S function for Coulomb potential")
plt.show()

grid,step=nu.spatial_objects.uniform_grid(-50,6,10**6,retstep=True)
grid=nu.spatial_objects.logarithmic_grid(grid)
coulomb=nu.potentials.Coulomb(grid)
nu.Numerov_algorithm.log_radial_numerov(grid, step, 1, 1, -1., 0, coulomb,l=2)
f.write(f"For the  radial algorithm using radial grid we find for n=1, l=2 energy {energy}\n")
nu.Numerov_algorithm.log_fast_radial_numerov(grid, step, 1, 1, -1., 0, coulomb,l=2)
f.write(f"For the fast radial algorithm using radial grid we find for n=1, l=2 energy {energy}\n")


f.close()

