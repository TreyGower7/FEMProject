import numpy as np

"""
This project aims to use 1-dimensional Finite Element Methods to solve a heat transfer PDE 
given initial and Dirichlet Boundary Conditions
"""

__author__ = "Trey Gower"

def fxt(x,t):
    return (((np.pi)**2-1)*np.exp(-t)*np.sin(np.pi*x))

def user_in():
    """ Function to gather user inputs for general purpose Finite Element code
     
     Args: N/A

     Returns: N = total # of global nodes, xl = left boundary, xr = right boundary
    """
    N = input("Enter N the total # of global nodes: \n")
    xl = input("Enter xl the left boundary: \n")
    xr = input("Enter xr the right boundary: \n")

    return {"N":N,"xl":xl,"xr":xr}

def user_in_heateq(x):
    """ Function to gather user inputs for the heat equation problem given
     
     Args: x0 (value of x at time t = 0)

     Returns: ux0 = initial condition, t0 = initial time, tf = final time, dt = time step
    """
    #using the dirichlet boundary condition to initialize u_n 
    u_n = np.array([np.sin(np.pi*x_i) for x_i in x]) 

    t0 = input("Enter initial time: \n")
    tf = input("Enter finial time: \n")
    dt = input("Enter how many time steps: \n")
    dt = 1/int(dt)
    return {"u_n":u_n,"t0":int(t0),"tf":int(tf),"dt":dt}

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

def basis(z):
    """ Function to map our uniform grid to a parent grid [-1, 1] via parent functions

     Args: h = the spacing between each element

     Returns: mapped values from x-space to xi-space
    """
    phi_1 = (1-z)/2 #first node parent function
    phi_2 = (1+z)/2 #second node parent function
    dxdz = h/2 #h is global value in main function
    dzdx = 2/h
    return {'phi_1': int(phi_1), 'phi_2': int(phi_2),'dxdz': int(dxdz), 'dzdx': int(dzdx)}

def element_mats(N,iee):
    Ne = N-1 #total # of 1D elements
    Kloc = np.zeros((2,2))
    Kglob = np.zeros((N,N))
    Mloc = np.zeros((2,2))
    Mglob = np.zeros((N,N))

    for i in range(Ne):
        for l in range(1):
            for m in range(1):
                Kloc[l][m] = 
                Kloc[l][m] = 


        #Global assembly 
        for l in range(1):
            global_node1 = iee[i][l]
            for m in range(1):
                global_node2 = iee[i][m]
                Kglob[global_node1][global_node2] += Kloc[l][m]
                Mglob[global_node1][global_node2] += Mloc[l][m]
    return Kglob, Mglob

def Assembly_f(uin, iee, uin_h):
    """ Function to assemble the global FE Mesh

     Args: N = total # of global nodes, iee = the uniform grid

     Returns: the global forcing function vector 
    """
    Nt = int((uin_h['tf']-uin_h['t0'])/(uin_h['dt']))
    ctime = uin_h['t0']
    N = uin['N']
    Ne = N-1 #total # of 1D elements
    floc = np.zeros((2,1))
    fglob = np.zeros((N,1))
    for n in range(Nt):
        ctime = uin_h['t0'] + n*uin_h['dt']
        for i in range(Ne):
            #Local calculations for each element 
            for l in range(1):
                floc[l] = solve_quad(uin['xl'],uin['xr'],ctime)

            #Global assembly 
            for l in range(1):
                global_node1 = int(iee[i][l])
                fglob[global_node1] += floc[l]
    return fglob
            

def solve_quad(a,b,t):
    """ Function to solve an integral with guassian quadrature and change variable to parent space

     Args: 

     Returns: Numerically approximated integral value
    """
    n = 2
    m = (b-a)/2
    c = (a+b)/2
    xi = np.zeros((n,1))
    guassT = np.array([-.5774,.5774]) #guassian quadrature points
    f = np.zeros((n,1)) #function values for mapped values
    guassW = np.array([1,1]) #weights
    result = 0
    #computing parent mapped x's and summing
    for i in range(n):
        b = basis(xi[i])
        xi[i] = m*guassT[i]+c
        f[i] = fxt(xi[i],t)
        result += guassW[i]*f[i]*b(['phi_' + str(i)])
    return b(['dxdz'])*result*m

def main():
    """ Main entry point of the app """
    

    uin = user_in()

    print("Making Unifrom Grid\n")
    print("--------------------\n")
    grid = uni_grid(uin['N'], uin['xl'], uin['xr'])

    #appending values for specific heat equation problem inputing x0
    uin_h = user_in_heateq(grid['x'])

    global h #make h global since it doesnt change and its the easiest solution
    h = (uin['xr']-uin['xl'])/uin['N']

    print("Creating Mass and Stiffness Matrices\n")
    print("------------------------------\n")
    M, K = element_mats(uin['N'],grid['iee'])

    print("Creating forcing vector via guass quadrature\n")
    print("------------------------------\n")
    f = Assembly_f(uin, grid['iee'], uin_h)

if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
