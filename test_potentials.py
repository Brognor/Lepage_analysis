#Here we test the objects and functions of class potentials inside Numerov library
import Numerov as nu

f=open("test_results/potentials_test.txt","w")


coulomb=nu.potentials.Coulomb(2.)
f.write(f"The coulomb potential at r=2 is {coulomb}\n")


well=nu.potentials.square_well(1.5,1,2,3)
f.write(f"The square_well potential at r=1.5 is {well}\n")

well=nu.potentials.radial_well(1.5,5,3)
f.write(f"Radial well at 1.5 is {well}\n")


angular=nu.potentials.angular_potential(1)
f.write(f"angular part of the potential at r=1 with l=0 and m=1 grid is {angular}\n")

#We can even obtain a vector-like object for the potential

grid=nu.spatial_objects.uniform_grid(1,10,10)
f.write(f"The grid is {grid}\n")

coulomb=nu.potentials.Coulomb(grid)
f.write(f"Coulomb potential along the grid is {coulomb}\n")

well=nu.potentials.square_well(grid,5,6,3)
f.write(f"Square well along grid is {well}\n")


well=nu.potentials.radial_well(grid,5,3)
f.write(f"Radial well along grid is {well}\n")


angular=nu.potentials.angular_potential(grid)
f.write(f"angular part of the potential with l=0 and m=1 along grid is {angular}\n")

angular=nu.potentials.angular_potential(grid,l=1)
f.write(f"angular part of the potential with l=1 and m=1 along grid is {angular}\n")

f.close()