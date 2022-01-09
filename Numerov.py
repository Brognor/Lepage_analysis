#Numerov library contain objects and many types of algorithms
#In all the definitions we have used natural units with h_bar=c=1

from numba import njit
import numpy as np
import copy
import scipy.integrate as int



def normalize(grid, func): #return the normalized eigenfunction
    if len(grid)!=len(func):
        raise ValueError("function and grid are not matching")

    func_2=np.multiply(func,func)
    norm=int.trapezoid(func_2,grid)
    return func/norm





class spatial_objects:  #functions that gives back grid-like elements



    
    def uniform_grid(x_min,x_max,N,retstep=False):
        if retstep:
            grid=np.linspace(x_min,x_max,num=N)
            return grid,grid[1]-grid[0]
        else:
            grid=np.linspace(x_min,x_max,num=N)
            return grid

    
    """
    This function return a 1xN  grid dividing uniformly the intervall between x_min and x_max in N point, x_min and x_max included

    If retstep is True the function return a 1x(N,1) array
    the first element is the same grid as before
    the second element is the space step
    """
    
    uniform_grid_step = np.arange

    "alias for np.arange(x_min, x_max, step) and return a grid that goes from start until x_max excluded, with step equal step"

    @njit
    def interval_step(uniform_grid):
        return abs(uniform_grid[1]-uniform_grid[0])
    
    "return the step in the case of a uniform grid, to avoid higher computation time we put no control over the uniformity of the grid"
    

    @njit 
    def logarithmic_grid(uniform_grid,Z=1,r_0=1):
        if 0 in uniform_grid:
            raise ValueError("0 in uniform_grid give -inf in log grid")
        else:
            return r_0/Z*np.exp(uniform_grid)
        

class auxiliary_functions: #in this class we define functions that enter in iteration of Numerov algorithm for Schrodinger equation.

    @njit
    def auxiliary_function_g(m,E,potential): 
        return 2*m*(E-potential)
    
    g_func=auxiliary_function_g #shorter alias
    
    @njit
    def auxiliary_function_f(g,step):
        step_2=step**2
        intermediate=np.multiply(g,step_2/12)
        return np.add(intermediate,1)

    f_func=auxiliary_function_f #shorter alias

    @njit
    def auxiliary_function_glog(m,E,potential,grid,l=0): 
        return 2*m*grid**2*(E-potential)-(l+0.5)**2
    
    g_funclog=auxiliary_function_glog #shorter alias
    
    @njit
    def auxiliary_function_flog(g,step):
        step_2=step**2
        return 1+g*step_2/12

    f_funclog=auxiliary_function_flog #shorter alias


class potentials: #class containig some known potentials
    
    @njit
    def Coulomb(r,alpha=1):
        r=np.reciprocal(r)
        return np.multiply(r,-alpha)



    
    def square_well(r,left=-1,right=1,value=1): #to make function work for both scalar and vector we used an if and control over type, we put also an error control looking if left<right
        if left>=right:
            raise ValueError("Left limit of the interval is greater or equal than right limit")
        
        else:
            R=copy.copy(r) #asarray modify the element that receveived so we need to create a copy of the grid
            a= np.asarray([R]) if np.isscalar(R) else np.asarray(R) #due to this line we cannot use njit here
            
            if len(a)==1:
                if left<=a<=right:
                    return value
                else:
                    return 0
            else:
                for i in range(len(a)):
                    if left<=a[i]<=right:
                        a[i]=value
                    else: 
                        a[i]=0
                return a

    
    
    def radial_well(r,radius=1,value=1): #to make function work for both scalar and vector we used an if and control over type, we put also an error control looking if radius>0
        if radius<0:
            raise ValueError("Negative radius as input")
        
        else:
            R=copy.copy(r) #asarray modify the element that receveived so we need to create a copy of the grid
            a= np.asarray([R]) if np.isscalar(R) else np.asarray(R) #due to this line we cannot use njit here
            
            if len(a)==1:
                if a<=radius: #a should be positive working in radial coordinates, we do not waste computation time adding a control.
                    return value
                else:
                    return 0
            else:
                for i in range(len(a)):
                    if a[i]<=radius:
                        a[i]=value
                    else: 
                        a[i]=0
                return a

    @njit
    def angular_potential(r,m=1,l=0): #define the part of the effective potential with l
        r=np.reciprocal(r**2)
        return np.multiply(r,l*(l+1)/(2*m))
    



class Numerov_algorithm: #here we define algorithm that perform the Numerov algorithm in the context of the Schrodinger equation
    """There are many ways in which one can find the correct energy, here we consider the method that find energy imposing continuity of first derivative."""
    
    def radial_numerov(grid,step,m,n_energy,E_min,E_max,potential,accuracy=10**(-5),dy=10**(-120)): 
        """grid is the grid on which apply the algorithm,  
        step is the grid step  
        m is the mass
        n_energy contain the number of the energy level, 
        dy indicates the starting value of the eigenfunction in the first point different from 0 (ad libitum due to normalization)
        E_min is the minimal value for the energy, E_max is the maximal value for the energy
        potential must be an array containing the value of the potential over the grid, we do not compute it inside the algorithm because it is energy independent, so if we put 
        the algorithm in a for cycle we avoid unuseful computation, of course the grid must be the same and we must put attention over this fact.
        accuracy is the uncertainity on the energy determined"""


        if E_max<E_min :
            raise ValueError("Maximum energy smaller than minimum energy")
        #here we do not consider divergencies in potential because the we allow only potenntial with divergency in r=0 and we put psi(0)=0 by theory


        else:

            energy=0

            node_number = n_energy#number of node, we use this to control if the tested energy is to high or to low
            
            lenght=len(grid)

            while(E_max-E_min>accuracy):

                index_inv=0
                energy=(E_max+E_min)/2.
                print("trying energy eigenvalue ",energy)
                node_counter=1


                #now we find the classical inversion point
                for i in range(lenght):
                    if potential[i]>=energy:
                        index_inv=i
                        
                        break
                    
                if index_inv==0:
                    raise ValueError("No inversion point found")
                
                eigenfunction_left=np.zeros(index_inv+1)
                eigenfunction_left[1]=dy#the lower the value, the lower is the risk to have an overfloor.
                eigenfunction_right=np.zeros(lenght-index_inv)
                eigenfunction_right[lenght-index_inv-2]=dy



                g=auxiliary_functions.g_func(m,energy,potential)
                f=auxiliary_functions.f_func(g,step)

                #due to avoid divergencies in r=0 we need to compute by hand the eigenvalue_left[2]
                eigenfunction_left[2]=(12-10*f[1])*eigenfunction_left[1]/f[2]
                

                for i in range (2,index_inv):
                    eigenfunction_left[i+1]=((12-10*f[i])*eigenfunction_left[i]-f[i-1]*eigenfunction_left[i-1])/f[i+1]
                    #we obtained the eigenfunction for the tested energy unitl inversion point

                #we know nodes must be before the inversion point, so we control if we have the right number of nodes before integrating from x_max backwards
                 #sign of each element
                sign_val=np.sign(eigenfunction_left) #we create an array of sign of the eigenfunction in order to find nodes
                
                for i in range(1,len(sign_val)-1):     #in this cycle we compute the number of node               
                    if sign_val[i+1]*sign_val[i]<=0:                        
                        node_counter+=1 
                
                
                if node_counter<node_number: #we need now to update our interval of interesting energies, before we use the node number
                    print("too little nodes")
                    E_min=energy
                    continue
                elif node_counter>node_number:
                    print("too many nodes")
                    E_max=energy
                    continue

                elif node_counter==node_number: #if the node number is right we proceed with backward integration
                    
                    

                    for i in range (lenght-index_inv-2,0,-1):
                        eigenfunction_right[i-1]=((12-10*f[index_inv+i])*eigenfunction_right[i]-f[index_inv+i+1]*eigenfunction_right[i+1])/f[index_inv+i-1]                         
                        #we obtained the eigenfunction for the tested energy from the last point backward to the inversion point 

                    #we need now to have continuity of the eigenfunction so we rescale the right eigenfunction
                    

                    eigenfunction_right=eigenfunction_right*(eigenfunction_left[index_inv]/eigenfunction_right[0])

                    #Now we compute the discontinuity of the first derivative in the inversion point
                    
                    discontinuity=(eigenfunction_right[1]+eigenfunction_left[index_inv-1]-(14.-12.*f[index_inv])*eigenfunction_left[index_inv])/step*(-1)**(n_energy-1)
                    


                    if discontinuity>0:
                        E_max=energy
                    else:
                        E_min=energy
                   

            eigenfunction_left=np.delete(eigenfunction_left,index_inv)
            eigenfunction=np.append(eigenfunction_left,eigenfunction_right) #put together the two part of eigenfunction
            return energy, eigenfunction

    ##########################################################################################################################################################
    ##########################################################################################################################################################

    @njit
    def fast_radial_numerov(grid,step,m,n_energy,E_min,E_max,potential,accuracy=10**(-6),dy=10**(-120)): 
        """the same function as before but defined using numba (some minor changes must be done to use it)"""
        #file=open("error.txt","w")


        if E_max<E_min :
            raise ValueError("Maximum energy smaller than minimum energy")
        #here we do not consider divergencies in potential because the we allow only potenntial with divergency in r=0 and we put psi(0)=0 by theory


        else:

            energy=0


            node_number = n_energy#number of node, we use this to control if the tested energy is to high or to low
            
            lenght=len(grid)

            while(E_max-E_min>accuracy):
                index_inv=0


                energy=(E_max+E_min)/2.
                node_counter=1
                print("trying energy eigenvalue ",energy)


                #now we find the classical inversion point
                for i in range(lenght):
                    if potential[i]>=energy:
                        index_inv=i
                        break
                if index_inv==0:
                    raise ValueError("No inversion point found")

                
                        

                
                eigenfunction_left=np.zeros(index_inv+1)
                eigenfunction_left[1]=dy#the lower the value, the lower is the risk to have an overfloor.
                eigenfunction_right=np.zeros(lenght-index_inv)
                eigenfunction_right[lenght-index_inv-2]=dy



                g=2*m*(np.add(energy,-potential))
                f=np.add(np.multiply(g,1/12*step**2),1)

                #due to avoid divergencies in r=0 we need to compute by hand the eigenvalue_left[2]
                eigenfunction_left[2]=(12-10*f[1])*eigenfunction_left[1]/f[2]
                

                for i in range (2,index_inv):
                    eigenfunction_left[i+1]=((12.-f[i]*10.)*eigenfunction_left[i]-f[i-1]*eigenfunction_left[i-1])/f[i+1]

                #we obtained the eigenfunction for the tested energy unitl inversion point
                #we know nodes must be before the inversion point, so we control if we have the right number of nodes before integrating from x_max backwards

                 #sign of each element
                sign_val=np.sign(eigenfunction_left) #we create an array of sign of the eigenfunction in order to find nodes
                
                for i in range(1,len(sign_val)-1):     #in this cycle we compute the number of node               
                    if sign_val[i+1]*sign_val[i]<=0:                        
                        node_counter+=1 
                        
                
                
                if node_counter<node_number: #we need now to update our interval of interesting energies, before we use the node number
                    E_min=energy
                    print("too little nodes")
                    continue
                    
                elif node_counter>node_number:
                    E_max=energy
                    print("too many nodes")
                    continue

                elif node_counter==node_number: #if the node number is right we proceed with backward integration
                    
                    

                    for i in range (lenght-index_inv-2,0,-1):
                        eigenfunction_right[i-1]=((12.-f[index_inv+i]*10.)*eigenfunction_right[i]-f[index_inv+i+1]*eigenfunction_right[i+1])/f[index_inv+i-1]                         
                        #we obtained the eigenfunction for the tested energy from the last point backward to the inversion point 

                    #we need now to have continuity of the eigenfunction so we rescale the right eigenfunction
                    
                    
                    eigenfunction_right=np.multiply(eigenfunction_right,(eigenfunction_left[index_inv]/eigenfunction_right[0]))
                    


                    #Now we compute the discontinuity of the first derivative in the inversion point


                    discontinuity=(eigenfunction_right[1]+eigenfunction_left[index_inv-1]-(14.-12.*f[index_inv])*eigenfunction_right[0])/step*(-1)**(n_energy-1)
                    
                    
                    


                    if discontinuity>0.:
                        E_max=energy
                    else:
                        E_min=energy
                   

            eigenfunction_left=np.delete(eigenfunction_left,index_inv)
            eigenfunction=np.append(eigenfunction_left,eigenfunction_right) #put together the two part of eigenfunction
            
            
            return energy, eigenfunction



    def log_radial_numerov(grid,step,m,n_energy,E_min,E_max,potential,l=0,accuracy=10**(-5),dy=10**(-120)): 
        """grid is the grid on which apply the algorithm, a good idea is to make it go from negative to positive values, we will transofrm it
        step is the grid step  
        m is the mass
        l is the angular momentum number,
        n_energy contain the number of the energy level, 
        dy indicates the starting value of the eigenfunction in the first point different from 0 (ad libitum due to normalization)
        E_min is the minimal value for the energy, E_max is the maximal value for the energy
        potential must be an array containing the value of the potential over the grid, we do not compute it inside the algorithm because it is energy independent, so if we put 
        the algorithm in a for cycle we avoid unuseful computation, of course the grid must be the same and we must put attention over this fact.
        accuracy is the uncertainity on the energy determined"""


        if E_max<E_min :
            raise ValueError("Maximum energy smaller than minimum energy")
        #here we do not consider divergencies in potential because the we allow only potenntial with divergency in r=0 and we put psi(0)=0 by theory


        else:

            energy=0

            node_number = n_energy#number of node, we use this to control if the tested energy is to high or to low
            
            lenght=len(grid)
            


            while(E_max-E_min>accuracy):

                index_inv=0
                energy=(E_max+E_min)/2.
                print("trying energy eigenvalue ",energy)


                #now we find the classical inversion point
                for i in range(lenght):
                    if potential[i]>=energy:
                        index_inv=i
                        
                        break
                    
                if index_inv==0:
                    raise ValueError("No inversion point found")
                
                eigenfunction_left=np.zeros(index_inv+1)
                eigenfunction_left[1]=dy#the lower the value, the lower is the risk to have an overfloor.
                eigenfunction_right=np.zeros(lenght-index_inv)
                eigenfunction_right[lenght-index_inv-2]=dy



                g=auxiliary_functions.g_funclog(m,energy,potential,grid,l=l)
                f=auxiliary_functions.f_funclog(g,step)

                #due to avoid divergencies in r=0 we need to compute by hand the eigenvalue_left[2]
                eigenfunction_left[2]=(12-10*f[1])*eigenfunction_left[1]/f[2]
                

                for i in range (2,index_inv):
                    eigenfunction_left[i+1]=((12-10*f[i])*eigenfunction_left[i]-f[i-1]*eigenfunction_left[i-1])/f[i+1]
                    #we obtained the eigenfunction for the tested energy unitl inversion point

                #we know nodes must be before the inversion point, so we control if we have the right number of nodes before integrating from x_max backwards

                sign_array=np.sign(eigenfunction_left) #sign of each element
                
                signchange = (np.absolute(np.roll(sign_array, 1) - sign_array))/2 #the subtraction can give 0 if no change or 2,-2 in there is a change in sign, so we need to take the abs and divde by 2!= 0).astype(int)
                node_counter=np.sum(signchange) #actual number of node in our eigenfunction
                
                
                if node_counter<node_number: #we need now to update our interval of interesting energies, before we use the node number
                    print("too little nodes")
                    E_min=energy
                    continue
                elif node_counter>node_number:
                    print("too many nodes")
                    E_max=energy
                    continue

                elif node_counter==node_number: #if the node number is right we proceed with backward integration
                    
                    

                    for i in range (lenght-index_inv-2,0,-1):
                        eigenfunction_right[i-1]=((12-10*f[index_inv+i])*eigenfunction_right[i]-f[index_inv+i+1]*eigenfunction_right[i+1])/f[index_inv+i-1]                         
                        #we obtained the eigenfunction for the tested energy from the last point backward to the inversion point 

                    #we need now to have continuity of the eigenfunction so we rescale the right eigenfunction
                    

                    eigenfunction_right=eigenfunction_right*(eigenfunction_left[index_inv]/eigenfunction_right[0])

                    #Now we compute the discontinuity of the first derivative in the inversion point
                    
                    discontinuity=(eigenfunction_right[1]+eigenfunction_left[index_inv-1]-(14.-12.*f[index_inv])*eigenfunction_left[index_inv])/step*(-1)**(n_energy-1)
                    


                    if discontinuity>0:
                        E_max=energy
                    else:
                        E_min=energy
                   

            eigenfunction_left=np.delete(eigenfunction_left,index_inv)
            eigenfunction=np.append(eigenfunction_left,eigenfunction_right) #put together the two part of eigenfunction
            return energy, eigenfunction


    @njit
    def log_fast_radial_numerov(grid,step,m,n_energy,E_min,E_max,potential,l=0,accuracy=10**(-6),dy=10**(-120)): 
        """the same function as before but defined using numba (some minor changes must be done to use it)"""
        #file=open("error.txt","w")


        if E_max<E_min :
            raise ValueError("Maximum energy smaller than minimum energy")
        #here we do not consider divergencies in potential because the we allow only potenntial with divergency in r=0 and we put psi(0)=0 by theory


        else:

            energy=0


            node_number = n_energy#number of node, we use this to control if the tested energy is to high or to low
            
            lenght=len(grid)

            while(E_max-E_min>accuracy):
                index_inv=0


                energy=(E_max+E_min)/2.
                node_counter=1
                print("trying energy eigenvalue ",energy)


                #now we find the classical inversion point
                for i in range(lenght):
                    if potential[i]>=energy:
                        index_inv=i
                        break
                if index_inv==0:
                    raise ValueError("No inversion point found")

                
                        

                
                eigenfunction_left=np.zeros(index_inv+1)
                eigenfunction_left[1]=dy#the lower the value, the lower is the risk to have an overfloor.
                eigenfunction_right=np.zeros(lenght-index_inv)
                eigenfunction_right[lenght-index_inv-2]=dy



                g=2*m*grid**2*(energy-potential)-(l+0.5)**2
                f=np.add(np.multiply(g,1/12*step**2),1)

                #due to avoid divergencies in r=0 we need to compute by hand the eigenvalue_left[2]
                eigenfunction_left[2]=(12-10*f[1])*eigenfunction_left[1]/f[2]
                

                for i in range (2,index_inv):
                    eigenfunction_left[i+1]=((12.-f[i]*10.)*eigenfunction_left[i]-f[i-1]*eigenfunction_left[i-1])/f[i+1]

                #we obtained the eigenfunction for the tested energy unitl inversion point
                #we know nodes must be before the inversion point, so we control if we have the right number of nodes before integrating from x_max backwards

                 #sign of each element
                sign_val=np.sign(eigenfunction_left) #we create an array of sign of the eigenfunction in order to find nodes
                
                for i in range(1,len(sign_val)-1):     #in this cycle we compute the number of node               
                    if sign_val[i+1]*sign_val[i]<=0:                        
                        node_counter+=1 
                        
                
                
                if node_counter<node_number: #we need now to update our interval of interesting energies, before we use the node number
                    E_min=energy
                    print("too little nodes")
                    continue
                    
                elif node_counter>node_number:
                    E_max=energy
                    print("too many nodes")
                    continue

                elif node_counter==node_number: #if the node number is right we proceed with backward integration
                    
                    

                    for i in range (lenght-index_inv-2,0,-1):
                        eigenfunction_right[i-1]=((12.-f[index_inv+i]*10.)*eigenfunction_right[i]-f[index_inv+i+1]*eigenfunction_right[i+1])/f[index_inv+i-1]                         
                        #we obtained the eigenfunction for the tested energy from the last point backward to the inversion point 

                    #we need now to have continuity of the eigenfunction so we rescale the right eigenfunction
                    
                    
                    eigenfunction_right=np.multiply(eigenfunction_right,(eigenfunction_left[index_inv]/eigenfunction_right[0]))
                    


                    #Now we compute the discontinuity of the first derivative in the inversion point


                    discontinuity=(eigenfunction_right[1]+eigenfunction_left[index_inv-1]-(14.-12.*f[index_inv])*eigenfunction_right[0])/step*(-1)**(n_energy-1)
                    
                    
                    


                    if discontinuity>0.:
                        E_max=energy
                    else:
                        E_min=energy
                   

            eigenfunction_left=np.delete(eigenfunction_left,index_inv)
            eigenfunction=np.append(eigenfunction_left,eigenfunction_right) #put together the two part of eigenfunction
            
            
            return energy, eigenfunction





                


