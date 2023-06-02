# Quasi1DEuler
## Project 2 - SURE Summer 2023

In this project I transferred the Quasi-1D Euler Code for a converging-diverging nozzle from Python to C++. In the process I attempted to employ more OOP techniques.

### General Questions
1. I wasn't always sure of the distinction or advantion of using a class and method as opposed to a function, this may be a point of ignorance for me.
2. Why not always pass by reference? Is it not better for optimization as then the value isn't copied?
3. Proper usage of lambdas? I used when a function was only being used within a limited scope
4. Is it better to use template arguments as I have here or allocate on the heap to get variable array sizes?


### Process
0. calculate initial condition

Each Iteration
1. calculate dt
2. calculate eigenval at boundary
3. calculate numerical flux & multiply by area 
4. Calculate residual and store density residual
5. Update the properties based on the residual value
6. Update the boundary conditions (using the old neighbouring values <- store somewhere)
7. Calculate the new P and speed of sound (try to do in same loop as the property update to speed up)

Break when the count reaches lim or R(rho) -> 0

9. Calculate Mach number 
10. Export P, M, R(rho) to a file for graphing in Python

