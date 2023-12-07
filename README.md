# Project 2
## by Trey Gower

## Overview
Solving the PDE Numerically with a 1-D FEM:
ut - uxx = f(x,t), where (x,t) ∈ (0,1) x (0,1)

### Initial conditions:
u(x,0) = sin(𝜋𝑥)
u(0,t) = u(1,t) = 0
					
Test function and analytical solution:
f(x,t) = (𝜋2 - 1)e-tsin(𝜋𝑥) & u(x,t) = e-tsin(𝜋𝑥)


## How to use the code:
For each of these enter your desired values and it will compute the solution numerically to the heat equation
Enter N the total # of global nodes:

Enter xl the left boundary:

Enter xr the right boundary:

Enter (FE) for Foward or (BE) for Backward Euler Method:

Enter initial time:

Enter finial time:

Enter number of time steps:

# Problem 2 (Forward Euler):
N = 11, Xl = 0, Xr = 1, initial time = 0, final time = 1, time steps = 551

<img src = https://github.com/TreyGower7/FEMProject/blob/main/fe_551.0_soln.png>

# Problem 3 (Backward Euler):
N = 11, Xl = 0, Xr = 1, initial time = 0, final time = 1, time steps = 551
<img src = https://github.com/TreyGower7/FEMProject/blob/main/be_551.0_soln.png>


# Conclusion:
Although I didn’t produce the correct answers fully I worked on this for a while and got far. Ultimately, I waited till the last minute to fully finish it and, unfortunately, it was not totally correct. The errors in my code arose from my K and M matrices and F vector being calculated somewhat incorrectly producing the wrong displacement vector. 

Here was what I got which produced the wrong results (My matrices and more importantly F vector was very sparse):

F:
[[0.00000000e+000]
 [0.00000000e+000]
 [0.00000000e+000]
 [0.00000000e+000]
 [0.00000000e+000]
 [0.00000000e+000]
 [0.00000000e+000]
 [0.00000000e+000]
 [0.00000000e+000]
 [3.14097022e-240]
 [0.00000000e+000]]

M:
[[11.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0. 11.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  0. 11.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0. 11.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0. 11.  0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0. 11.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0. 11.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0. 11.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0. 11.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0. 11.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]]
K:

[[0.06060867 0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.        ]
 [0.         0.06060867 0.         0.         0.         0.
  0.         0.         0.         0.         0.        ]
 [0.         0.         0.06060867 0.         0.         0.
  0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.06060867 0.         0.
  0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.06060867 0.
  0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.06060867
  0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.06060867 0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.06060867 0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.06060867 0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.06060867 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.        ]]

### Applying Dirichlet Boundary Conditions

New M Matrix:
[[ 1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0. 11.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  0. 11.  0.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0. 11.  0.  0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0. 11.  0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0. 11.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0. 11.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0. 11.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0. 11.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0. 11.  0.]
 [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.]]
