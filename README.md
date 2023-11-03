# Viral_Uptake
These codes generate initial positions of the viral uptake simulations, run dynamics, and analyze the results to quantify the virus wrapping, cell surface wrapping, and their time-dependent properties.


Initial.f90 code generates the initial conditions of the viral uptake simulations. This requires initial files sphere_fibonacci_grid_n1000.xyz, Sphere_Data.dat, and Sphere_Edge.dat, which contain the initial coordinates of the spherical viral surface, the number of nodes, and the total number of edges, and edge definitions, respectively. 

Trig_March.f90 runs the dynamics for the viral uptake and dumps the trajectory data.

Quantify_Endocytosis_New.f90 quantifies the trajectory data by analyzing the final shape of the cell surface while it is wrapped around the virus. It calculates the shape index cell surface area used in uptake, the virus wrapping, how much virus surface is covered with the membrane, and the velocity of viral uptake.
