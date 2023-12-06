import numpy as np
import matplotlib.pyplot as plt
"""
This project aims to use 1-dimensional Finite Element Methods to solve a heat transfer PDE 
given initial and Dirichlet Boundary Conditions
"""

__author__ = "Trey Gower"

def fxt(x,t):
    return (((np.pi)**2-1)*np.exp(-t)*np.sin(np.pi*x))

def analytical_soln(x, t):
    return (np.exp(-t) * np.sin(np.pi * x))

def user_in():
    """ Function to gather user inputs for general purpose Finite Element code
     
     Args: N/A

     Returns: N = total # of global nodes, xl = left boundary, xr = right boundary
    """
    N = input("Enter N the total # of global nodes: \n")
    xl = input("Enter xl the left boundary: \n")
    xr = input("Enter xr the right boundary: \n")
    emeth = input("Enter (FE) for Foward or (BE) for Backward Euler Method: \n")

    return {"N":int(N),"xl":int(xl),"xr":int(xr), 'emeth': emeth}

def user_in_heateq(N,x):
    """ Function to gather user inputs for the heat equation problem given
     
     Args: x0 (value of x at time t = 0)

     Returns: ux0 = initial condition, t0 = initial time, tf = final time, dt = time steps
    """

    t0 = input("Enter initial time: \n")
    tf = input("Enter finial time: \n")
    ts = input("Enter number of time steps: \n")
    dt = 1/int(ts)
    #using the dirichlet boundary condition to initialize u_n 
    u_n = np.array([np.sin(np.pi*x_i) for x_i in x]) 

    return {"u_n":u_n,"t0":int(t0),"tf":int(tf),"dt":dt}

def uni_grid(N, xl,xr,):
    """ Function to create a uniform grid
     
     Args: N = total # of global nodes, xl = left boundary, xr = right boundary

     Returns: iee = the unifrom grid, x (vector) = [xl, xvalues ,xr]
    """
    x = np.zeros((N,1)) #Define x initializing to zero
    iee = np.zeros((N,2)) # Nx2 matrix since in 1d there are only ever 2 nodes
    Ne = N-1 #total # of 1D elements
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
    dphi = np.array([-1 / 2, 1 / 2])
    dxdz = h/2 #h is global value in main function
    dzdx = 2/h
    return {'phi_1': int(phi_1), 'phi_2': int(phi_2), 'dphi': dphi, 
            'dxdz': int(dxdz), 'dzdx': int(dzdx)}

def element_mats(N,iee,ts):
    """ Compute K and M matrices

     Args: N and the uniform grid iee

     Returns: K and M matrices
    """
    n = np.linspace(0,1,int(ts)+1)
    Ne = N-1 #total # of 1D elements
    Kloc = np.zeros((2,2))
    Kglob = np.zeros((N,N))
    Mloc = np.zeros((2,2))
    Mglob = np.zeros((N,N))
    Fglob = np.zeros((N, 1))
    P, W = quad_points()
    b0 = basis(P[0]) # mapping to xi space with our guassian quadtrature points -.5773
    b1 = basis(P[1]) # mapping to xi space with our guassian quadtrature point .5773
    dphi = b0['dphi']
    #Phi matrix for ease of use
    phi = np.array([[b0['phi_1'], b1['phi_1']], 
                    [b0['phi_2'], b1['phi_2']]])
                    
    for i in range(Ne):
        for l in range(1):
            for m in range(1):
                
                Mloc[l][m] = (phi[l, 0] * phi[m, 0] + phi[l, 1] * phi[m, 1]) * h
                Kloc[l][m] = int(dphi[l]) * b0['dzdx'] * int(dphi[m]) * b0['dzdx'] * b0['dxdz'] * 2


        #Global assembly 
        for l in range(1):
            global_node1 = int(iee[i][l])
            for m in range(1):
                global_node2 = int(iee[i][m])
                Kglob[global_node1][global_node2] += Kloc[l][m]
                Mglob[global_node1][global_node2] += Mloc[l][m]
    
    #P[0] = -.5774 P[1] = .5774
    Fglob[i] = -(fxt((P[0]), ts) * phi[0,0] + fxt((P[1]), ts) * phi[0,1]) * (1/8)

    return Kglob, Mglob, Fglob

            
def quad_points():
    """ Function to get guassian quadraure points
     Returns: Guassian weights and points for N=2 
    """
    guassT = np.array([-.5774,.5774])#guassian quadrature points
    guassW = np.array([1,1])
    
    return guassT, guassW


def bconds(Mglob, N):
    """ Function applys dirichlet boundary conditions

     Args: Global mass Matrix and N

     Returns:  New global mass matrix
    """
    Mglob[0, :] = 0
    Mglob[:, 0] = 0
    Mglob[N - 1, :] = 0
    Mglob[:, N - 1] = 0
    Mglob[0, 0] = 1
    Mglob[N - 1, N - 1] = 1
    
    bounds = [0,0]
    DBCs = np.eye(N)
    DBCs[0,0] = bounds[0]
    DBCs[N-1,N-1] = bounds[1]
    return Mglob, DBCs
    
def solve_kuf(u_n, febe, dt, M, K, F, DBCs):
    """ Solves the Ku=f system

     Args: u_n =initial displacemnts, febe = forward backward euler, M = global mass mat, K = global K mat, 
     F = global forcing vector, dt = 1/timestep,DBCs = dirichlet bounds

     Returns: the displacement vector u_n
    """
    #inverting M and creating and inverting matrix B
    M_inv = np.linalg.inv(M)
    M_K = np.dot(M_inv, K)
    B = (1/dt) * M + K
    B_inv = np.linalg.inv(B)
    if febe == 'FE':
        for t in range(1/dt):
            u_n[:, t + 1] = u_n[:, t] - dt * M_K.dot(u_n[:, t]) + dt * M_inv.dot(F[:, t])
            u_n[:, t + 1] = DBCs.dot(u_n[:, t + 1])
    else:
        for t in range(1/dt):
            u_n[:, t + 1] = (1 / dt) * B_inv.dot(M.dot(u_n[:, t])) + B_inv.dot(F[:, t])
            u_n[:, t + 1] = DBCs.dot(u_n[:, t + 1])
    
    return u_n

def plot_soln(u, ts,febe,N):
    x = np.linspace(0, 1, N)
    xn = np.linspace(0, 1, 1000)
    sol = analytical_soln(xn,1)
    plt.plot(xn, sol, label='Analytical Solution (Exact)')
    plt.plot(x, u[:, ts], label=f'{febe} Numerical Solution with {ts} time steps')
    plt.xlabel('x')
    plt.ylabel('Solution')
    plt.title('Comparing Analytical and Numerical Solutions')
    plt.legend()
    plt.savefig(f'{febe}_{ts}_soln.png')

def main():
    """ Main entry point of the app """
    

    uin = user_in()

    global h #make h global since it doesnt change and its the easiest solution
    h = (uin['xr']-uin['xl'])/uin['N']

    print("Making Unifrom Grid\n")
    print("--------------------\n")
    grid = uni_grid(uin['N'], uin['xl'], uin['xr'])

    #appending values for specific heat equation problem inputing x0
    uin_h = user_in_heateq(uin['N'], grid['x'])


    print("Creating Mass and Stiffness Matrices and F Vector\n")
    print("------------------------------\n")
    M, K, F = element_mats(uin['N'],grid['iee'],(1/uin_h['dt']))

    print("F:\n")
    print(F)
    print("M:\n")
    print(M)
    print("K:\n")
    print(K)

    print("Applying Dirichlet Boundary Conditions\n")
    print("------------------------------\n")
    M, DBCs = bconds(M, uin['N'])

    print("New M Matrix:\n")
    print(M)

    print("Solving Ku=f System\n")
    print("------------------------------\n")
    u = solve_kuf(uin_h['u_n'], uin['emeth'],uin_h['dt'], M, K, F, DBCs)

    print("Plotting analytical vs numerical")
    print("------------------------------\n")
    plot_soln(u, (1/uin_h['dt']), uin['emeth'],uin['N'])

if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
