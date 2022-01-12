import Numerov as nu       
import numpy as np                                           
from numba import njit

f=open("results/delta_with_a.txt","w")

grid, step = nu.spatial_objects.uniform_grid(-50,6,10**6,retstep=True)
grid=nu.spatial_objects.logarithmic_grid(grid)

value =-1.2
radius=2
#value of c and d to be found
c_delta=0
c_eff=0
d_eff=0

#function to compute the sum of relative errors
@njit
def error_sum(data,energy):
    error=0
    for data, energy in data, energy:
        error+=abs(data-energy)/abs(data)
    return error

#initialite fixed potentials 

coulomb=nu.potentials.Coulomb(grid)
well=nu.potentials.radial_well(grid,radius,value)
potential=coulomb+well

#initialize data array

data = []
coulomb_energy = []
dirac_approx_data = []
smeared_delta_data = []
effective_potential_data = []


for n in range(1,11):
    energy=nu.Numerov_algorithm.log_fast_radial_numerov(grid,step,1,n,-2,0,potential,dy=10**-150,accuracy=10**-10)[0] #data generation
    data.append(abs(energy))
    energy, eigenfunc=nu.Numerov_algorithm.log_fast_radial_numerov(grid,step,1,n,-2,0,coulomb,dy=10**-150,accuracy=10**-10)[0] #coulomb energy
    coulomb_energy.append(abs(energy))



c_max=0
c_min=-50000
errors=[]    
c_values=[]
for c in range(c_min,c_max):
    c=c/10000
    c_values.append(c)
    for n in range(0,10):
        dirac_approx_data.append(coulomb_energy[n]+c*1/(np.sqrt(np.pi)*(n+1)**3))

    errors.append(error_sum(data,dirac_approx_data))
    dirac_approx_data=[]       


pos=np.argmin(errors)
c=c_values[pos]
f.write(f"Minimixe square error  gives  c= {c}\n")

for n in range(0,10):
    dirac_approx_data.append(coulomb_energy[n]+c*1/(np.sqrt(np.pi)*(n+1)**3))

f.write(f"Our approximation gives energies {dirac_approx_data}\n\n")
f.write(f"Data are:{data}\n\n")
n=10
c = (data[n-1]+1/(2*n**2))*np.sqrt(np.pi)*n**3
f.write(f"Fixing c from  level {n}S will give: c= {c}\n")

dirac_approx_data = []

for n in range(0,10):
    dirac_approx_data.append(coulomb_energy[n]+c*1/(np.sqrt(np.pi)*(n+1)**3))
f.write(f"Energy eigenvalues in this case are {dirac_approx_data}\n")

f.close()
    

    

