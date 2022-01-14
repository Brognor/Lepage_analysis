import Numerov as nu       
import numpy as np                                           
from tqdm import tqdm

f=open("results/delta_with_a.txt","w")
f.write("We now find best value for c for a regulated Coulomb+smeared delta potential for different a\n\n")


value =-1.2
radius=2
#value of c and d to be found


a = [1,0.5,0.1]


#function to compute the sum of relative errors

def error_sum(data,energy):
    error=0
    for datas, energies in zip(data, energy):
        error+=(datas-energies)**2
    return error





#initialize data array

data = [-1.568949087013607, -0.19602247419243213, -0.07505287685489748, -0.039067748657544143, -0.023876894192653708, 
-0.016084664603113197, -0.011565622189664282, -0.008714093928574584, -0.006800400515203364, -0.0054541801364393905]


smeared_delta_data = []




f.write(f"Data are:{data}\n\n")









for a in a:

    f.write(f"For a={a}\n\n")
    c_max=round(-30-3/a)*100
    c_min=round(-70-3/a)*100
    smeared_delta_data=[]

    
    errors=[]  
    difference_low_energy = []  
    c_values=[]


    

    for c in tqdm(range(c_min,c_max)):
        
        c=c/100
        c_values.append(c)




        for n in range(1,11):
            grid, step = nu.spatial_objects.uniform_grid(0,n*50,10**5,retstep=True)
            reg_coulomb=nu.potentials.Coulomb_regulated(grid,a)
            delta=nu.potentials.smeared_dirac_delta(grid,a)
            potential=reg_coulomb+c*delta*a**2 #writing down the potential with regulated coulomb + smeared delta
            
            energy_delta = nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-100,0,potential,accuracy=10**-10,dy=10**-150)[0]
            smeared_delta_data.append(energy_delta)

        errors.append(error_sum(data,smeared_delta_data)) #compute error for fit
        difference_low_energy.append((abs(data[9]-smeared_delta_data[9]))) #compute difference for the 10S level
        smeared_delta_data=[]  
   


    pos=np.argmin(errors)
    c=c_values[pos]
   
    

    for n in range(1,11):
        grid, step = nu.spatial_objects.uniform_grid(0,n*50,10**5,retstep=True)
        reg_coulomb=nu.potentials.Coulomb_regulated(grid,a)
        delta=nu.potentials.smeared_dirac_delta(grid,a)
        potential=reg_coulomb+c*delta*a**2
        energy_delta = nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-100,0,potential,accuracy=10**-10,dy=10**-150)[0]
        smeared_delta_data.append(energy_delta)

    
    f.write(f"Minimixe square error  gives  c= {c}\n")
    f.write(f"Our approximation gives energies {smeared_delta_data}\n\n")
    
    pos=np.argmin(difference_low_energy)
    c=c_values[pos]
    

    smeared_delta_data = []

    for n in range(1,11):
        grid, step = nu.spatial_objects.uniform_grid(0,n*50,10**5,retstep=True)
        reg_coulomb=nu.potentials.Coulomb_regulated(grid,a)
        delta=nu.potentials.smeared_dirac_delta(grid,a)
        potential=reg_coulomb+c*delta*a**2
        energy_delta = nu.Numerov_algorithm.fast_radial_numerov(grid,step,1,n,-100,0,potential,accuracy=10**-10,dy=10**-150)[0]
        smeared_delta_data.append(energy_delta) 
    
    n=10
    
    f.write(f"Fixing c from  level {n}S will give: c= {c}\n")
    f.write(f"Energy eigenvalues in this case are {smeared_delta_data}\n\n")

f.close()
    

    

