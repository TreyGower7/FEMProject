import numpy as np

"""
This project aims to use 1-dimensional Finite Element Methods to solve a heat transfer PDE 
with initial and Dirichlet Boundary Conditions
"""

__author__ = "Trey Gower"


def user_in():
    """ Function to gather user inputs for general purpose Finite Element code
     
     Args: N/A

     Returns: N = total # of global nodes, xl = left boundary, xr = right boundary, 
     fx = external forcing function, Boundary Conditions
    """
    N = input("Enter N the total # of global nodes: \n")
    xl = input("Enter xl the left boundary: \n")
    xr = input("Enter xr the right boundary: \n")
    fx = input("Enter fx the forcing function: \n")
    BCs = input("Enter BCs the boundary conditions: \n")

    return [{"N":N,"xl":xl,"xr":xr,"fx":fx, "BCs": BCs}]

def user_in_heateq(x0):
    """ Function to gather user inputs for the heat equation problem given
     
     Args: x0 (value of x at time t = 0)

     Returns: ux0 = initial condition, t0 = initial time, tf = final time, dt = time step
    """
    ux0 = np.sin(np.pi*x0)
    t0 = input("Enter initial time: \n")
    tf = input("Enter finial time: \n")
    dt = input("Enter time step: \n")

    return [{"ux0":ux0,"t0":t0,"tf":tf,"dt":dt}]

def uni_grid(N, xl,xr,):
    """ Function to create a uniform grid
     
     Args: N = total # of global nodes, xl = left boundary, xr = right boundary

     Returns: iee = the unifrom grid, x (vector) = [xl, xvalues ,xr]
    """
    x = np.array(np.shape(N)) #Define x
    iee = np.array(np.shape(N,2)) # Nx2 matrix since in 1d there are only ever 2 nodes
    Ne = N-1 #total # of 1D elements
    h = (xr-xl)/(Ne)
    x[N] = 0
    iee[Ne][2] = 0

def parent_grid(h):
    """ Function to map our uniform grid to a parent grid [-1, 1] via parent functions

     Args: h = the spacing between each element

     Returns: function value
    """
    phi_1 = (1-z)/2 #first parent function
    phi_2 = (1+z)/2 #second parent function
    dxdz = h/2
    dzdx = 2/h

    #Main loop to create the uniform grid and connectivity map
    for i in range(1,N-1,1):
        x[i] = xl + (i-1)*h
        iee[i][0] = i
        iee[i][1] = i+1
    
    x[N] = xr

def main():
    """ Main entry point of the app """
    print("hello world")


if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()