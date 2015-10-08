#! /usr/bin/env python

import numpy as np
import warnings

def interpolate_cell(cell1,cell2,Lambda):
    """
    Generate interpolated set of lattice vectors R = R' + Lambda(R"-R')

    For mathematical basis, see "interpolation_scheme.pdf", included in the dvxc source directory.

    arguments:
        cell1, cell2: 3x3 numpy arrays of lattice vectors
        Lambda: scalar parameter of interpolation (usually between 0 and 1)

    returns:
        3x3 numpy array of new interpolated lattice vectors

    """

    cell = cell1 + Lambda * (cell2-cell1)
    return cell

def cell_volume(cell):
    """
    Calculate cell volume from lattice vectors

    arguments:
        cell: 3x3 numpy array of lattice vectors
    returns:
        volume: scalar volume in same units as lattice vectors

    """

    cell = np.array(cell)
    if len(cell) == 3 and cell.size == 9:
        v = np.dot(cell[0,:],np.cross(cell[1,:],cell[2,:]))
        return v
    else:
        raise Exception("Cell is not a 3x3 array")

def compute_x(cell1, cell2):
    """
    Compute the vector x of dot/cross product permutations (see supplementary information)

    arguments:
        cell1, cell2: 3x3 numpy arrays of lattice vectors

    returns:
        x: vector of dot/cross product permutations 
           a[cell1,cell2] . b[cell1,cell2] x c[cell1,cell2]

    """
    def f(R1, R2, R3):
        return np.dot(R1[0,:],np.cross(R2[1,:],R3[2,:]))

    perm = [(cell1,cell1,cell1),
            (cell1,cell1,cell2),
            (cell1,cell2,cell1),
            (cell1,cell2,cell2),
            (cell2,cell1,cell1),
            (cell2,cell1,cell2),
            (cell2,cell2,cell1),
            (cell2,cell2,cell2)]
    x = []
    for (R1, R2, R3) in perm:
        x.append(f(R1, R2, R3))
    return x

    
def solve_v(cell1, cell2, v):
    """
    Find scale factor needed to obtain target volume by solving cubic equation

    arguments:
        cell1, cell2: 3x3 numpy arrays of lattice vectors
        v: desired volume (for interpolation process this should be between the volumes of cell1 and cell2)

    returns:
        Lambda: scalar parameter of interpolation between cell1 and cell2 to obtain volume v

    raises:
        warnings: if target volume is not between the input cells, we are extrapolating!

    """

    # Check volume is in range, warn if not
    cell1_v = cell_volume(cell1)
    cell2_v = cell_volume(cell2)
    if ((v > cell1_v and v > cell2_v)
        or (v < cell1_v and v < cell2_v)):
        warnings.warn("volume {0} falls outside range of cells {1}-{2}. We are extrapolating!".format(
            v, cell1_v, cell2_v))
        
    x = compute_x(cell1,cell2)
    A = [[-1,1,1,-1, 1,-1,-1,1],
         [3,-2,-2,1,-2, 1,1,0],
         [-3,1,1, 0, 1, 0,0,0],
         [1, 0,0, 0, 0, 0,0,0]]
    cubic_coefficients = np.dot(A,x)
    # Subtract v from constant term to form cubic equation
    cubic_coefficients += [0,0,0,-v]
    roots = np.roots(cubic_coefficients)

    # Roots of a cubic are complex. We require a real, preferably small number
    real_solutions = []
    for solution in roots:
        if solution == solution.real:
            real_solutions.append(solution)
    if len(real_solutions) > 0:
        return min(real_solutions).real
    else:
        warnings.warn("No real solutions.")
        return False


def check_algebra(cell1,cell2,Lambda):
    x = compute_x(cell1,cell2)
    A = [[-1,1,1,-1, 1,-1,-1,1],
         [3,-2,-2,1,-2, 1,1,0],
         [-3,1,1, 0, 1, 0,0,0],
         [1, 0,0, 0, 0, 0,0,0]]
    y = [Lambda**3,Lambda**2,Lambda,1]
    return np.dot(np.dot(A,x),y)
    
def main():
    cell1 = np.array([[1, 0, 0],
                  [0,1,0],
                  [0,0,1]
              ])
    cell2 = np.array([[1.9,0,0],
                  [0,1.9,0],
                  [0,0,1.9]
              ])

    print "Simple interpolated volume with Lambda = 0.397"
    print cell_volume(interpolate_cell(cell1,cell2,0.397))
    print "Volume from algebraic manipulation with Lambda = 0.397"
    print check_algebra(cell1,cell2,0.397)

    print solve_v(cell1,cell2,2.5)
    
if __name__ == "__main__":
    main()
        
