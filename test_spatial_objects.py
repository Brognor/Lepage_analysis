#program to test functions defined in Numerov.spatial_objects

import Numerov as nu
f=open("test_results/spatial_object_text.txt","w")
grid,step=nu.spatial_objects.uniform_grid(0.,10.,11,retstep=True) 
f.write(f"grid = {grid}\n")
f.write(f"step = {step}\n")

grid=nu.spatial_objects.uniform_grid_step(0,10.1,1)
f.write(f"grid = {grid}\n")

step=nu.spatial_objects.interval_step(grid)
f.write(f"step = {step}\n")

grid=nu.spatial_objects.uniform_grid_step(1,10.1,1) #we have changed the grid because 0 in the grid raise an error as we have choose in library Numerov

log_grid=nu.spatial_objects.logarithmic_grid(grid, 1,1)
f.write(f"log_grid is {log_grid}\n")

