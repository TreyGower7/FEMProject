import numpy as np
"""
Module Docstring
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

    return {"N":int(N),"xl":int(xl),"xr":int(xr),"fx":fx, "BCs": BCs}

def uni_grid(N, xl,xr,):
    """ Function to create a uniform grid
     
     Args: N = total # of global nodes, xl = left boundary, xr = right boundary

     Returns: iee = the unifrom grid, x (vector) = [xl, xvalues ,xr]
    """
    x = np.zeros((N,1)) #Define x initializing to zero
    iee = np.zeros((N,2)) # Nx2 matrix since in 1d there are only ever 2 nodes
    Ne = N-1 #total # of 1D elements
    h = (xr-xl)/(Ne)
    #Main loop to create the uniform grid and connectivity map
    for i in range(0,N,1):
        x[i] = xl + (i)*h
        iee[i][0] = i
        iee[i][1] = i+1
    
    x[N-1] = xr

    return {'iee': iee, 'x': x}
def te():
    print(h)
def main():
    """ Main entry point of the app """
    uin = user_in()
    iee = uni_grid(uin['N'],uin['xl'], uin['xr'])
    global h
    h = (uin['xr']-uin['xl'])/uin['N']

    print(uin['N'])
    print('')
    print(iee['iee'])
    print('')
    print(iee['x'])
    
    te()



if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()