import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def VH(x,y,z):
    """
    This function determines the potential value for a given position    
    """

    r = np.sqrt(x**2+y**2+z**2)
    V = 0.5*(3*x**2 - z**2) + 1/r 
    return V

def JH(x,y,z,v):
    """   
    This calculatesthe jacobi integral for a given point in Hill 3BP
    """
    V = VH(x,y,z)
    JH = 0.5*v**2 - V
    return JH

def allowable_region(J,xmin,xmax,ymin,ymax,zmin,zmax,n):
    """   
    This function takes in a value of J and bounds for a grid and returns a matrix 
    stating if the point is allowable. If a value is zero (ie z=0), set the min and
    max values for that variable to both zero

    Args:
        J: The jacobi intrgral value
        xmin,xmax: min and max values for the x coord
        n: number of points for each direction

    
    """

    # Creating values of the grid
    matrix = []
    if xmin and ymin:
        print('XY')
        d1 = np.linspace(xmin,xmax,n)
        d2 = np.linspace(ymin,ymax,n)
        # Nested For Loop
        for val1 in d1:  # x outer loop
            sublist = []
            for val2 in d2:  # y inner loop
                if VH(val1,val2,0) > np.abs(J):
                    sublist.append(1)
                else:
                    sublist.append(0)
            matrix.append(sublist)

    elif xmin and zmin:
        print('XZ')
        d1 = np.linspace(xmin,xmax,n)
        d2 = np.linspace(zmin,zmax,n)
        # Nested For Loop
        for val1 in d1:  # x outer loop
            sublist = []
            for val2 in d2:  # y inner loop
                if VH(val1,0,val2) > np.abs(J):
                    sublist.append(1)
                else:
                    sublist.append(0)
            matrix.append(sublist)
    elif zmin and ymin:
        print('YZ')
        d1 = np.linspace(ymin,ymax,n)
        d2 = np.linspace(zmin,zmax,n)
        # Nested For Loop
        for val1 in d1:  # x outer loop
            sublist = []
            for val2 in d2:  # y inner loop
                if VH(0,val1,val2) > np.abs(J):
                    sublist.append(1)
                else:
                    sublist.append(0)
            matrix.append(sublist)
    matrix = np.linalg.matrix_transpose(matrix)
    return matrix

def graph_contour_region(J,xmin,xmax,ymin,ymax,zmin,zmax,n,Xt,Yt,Z_tot):

    matrix = allowable_region(J,xmin,xmax,ymin,ymax,zmin,zmax,n)
    matrix = np.array(matrix,dtype=float)
    
    plt.figure(figsize=(7,6))
    from matplotlib.colors import ListedColormap
    binary_cmap = ListedColormap(["red", "green"])

    plt.imshow(matrix, cmap=binary_cmap,origin='lower',
            extent=[xmin, xmax, zmin, zmax], vmin=0, vmax=1, alpha=0.65, aspect='auto')
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title("Allowable (green) vs Unallowable (red) Regions")

    plt.show()

