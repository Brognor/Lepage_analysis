import Numerov as nu       
import numpy as np                                           
from tqdm import tqdm

f=open("results/effective.txt","a")



value =-1.2
radius=2
#value of c and d to be found


a = [0.1]


#function to compute the sum of relative errors

def error_sum(data,energy):
    error=0
    for datas, energies in zip(data, energy):
        error+=(datas-energies)**2
    return error





#initialize data array

data = [-1.568949087013607, -0.19602247419243213, -0.07505287685489748, -0.039067748657544143, -0.023876894192653708, 
-0.016084664603113197, -0.011565622189664282, -0.008714093928574584, -0.006800400515203364, -0.0054541801364393905]


effective_data = []














for a in a:

    f.write(f"For a={a}\n\n")

    effective_data=[]

    
    errors=[]  
    difference_low_energy = []  
    values=[]


    

    for c in tqdm(range(-11050,-11008)):
        
        c=c/100.
        
        
        for d in range(-2400,-2300):
            d=d/100
            values.append([c,d])





            for n in range(1,11):

                grid, step = nu.spatial_objects.uniform_grid(0,n*50,10**4,retstep=True)
                potential=nu.potentials.effective_potential(grid,a,c,d)#writing down the full effective potential at order a**4

                energy_effective = nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-100,0,potential,accuracy=10**-10,dy=10**-150)[0]
                effective_data.append(energy_effective)

            errors.append(error_sum(data,effective_data)) #compute error for fit
            effective_data=[]  
   


    pos=np.argmin(errors)
    c, d= values[pos]
    
   
    

    for n in range(1,11):
        grid, step = nu.spatial_objects.uniform_grid(0,n*50,10**5,retstep=True)

        potential=nu.potentials.effective_potential(grid,a,c,d)
        energy_effective = nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-100,0,potential,accuracy=10**-10,dy=10**-150)[0]
        effective_data.append(energy_effective)

    
    f.write(f"Minimixe square error  gives  c= {c} and d={d}\n")
    f.write(f"Our approximation gives energies {effective_data}\n\n")
    f.write(f"Error sum {error_sum(data,effective_data)}\n\n")


f.close()
    

    

