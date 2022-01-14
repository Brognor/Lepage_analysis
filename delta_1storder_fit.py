import Numerov as nu       
import numpy as np                                           



f=open("results/delta_perturbation.txt","w")

grid, step = nu.spatial_objects.uniform_grid(0,300,10**6,retstep=True)



value =-1.2
radius=2

f.write(f"Data are generated with a square well potential with value = {value} and radius = {radius}\n\n")




#function to compute the sum of relative errors

def error_sum(data,energy):
    error=0
    for datas, energies in zip(data, energy):
        error+=(datas-energies)**2
    return error

#initialite fixed potentials 

coulomb=nu.potentials.Coulomb(grid)


#initialize data array


data = [-1.568949087013607, -0.19602247419243213, -0.07505287685489748, -0.039067748657544143, -0.023876894192653708, 
-0.016084664603113197, -0.011565622189664282, -0.008714093928574584, -0.006800400515203364, -0.0054541801364393905]

coulomb_energy = []
dirac_approx_data = []



for n in range(1,11):
    energy=nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-5,0,coulomb,dy=10**-150,accuracy=10**-10)[0] #coulomb energy
    coulomb_energy.append(energy)
f.write(f"Data are:{data}\n\n")
print(coulomb_energy)
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

n=10
c = (data[n-1]+1/(2*n**2))*np.sqrt(np.pi)*n**3
f.write(f"Fixing c from  level {n}S will give: c= {c}\n")

dirac_approx_data = []

for n in range(0,10):
    dirac_approx_data.append(coulomb_energy[n]+c*1/(np.sqrt(np.pi)*(n+1)**3))
f.write(f"Energy eigenvalues in this case are {dirac_approx_data}\n")

f.close()