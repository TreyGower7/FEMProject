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

    return {"N":N,"xl":xl,"xr":xr,"fx":fx, "BCs": BCs}

def main():
    """ Main entry point of the app """
    uin = user_in()

    print(uin['N'])


if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()