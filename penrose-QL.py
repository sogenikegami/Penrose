import numpy as np
import matplotlib.pyplot as plt
import math
import Penrose_tile02

vertexdata = Penrose_tile02.pen()

def Hsite(vertexdata):#on site energy
    e_s = 1.0
    e_px = 2.0
    e_py = 2.0
    
    N = len(vertexdata)
    zeroN = np.zeros((N,N))*0j
    idenN = np.indentity(N)
    H = np.zeros((6*N,6*N))*0j

    a1 = np.hstack([idenN*e_s,zeroN,zeroN,zeroN,zeroN,zeroN])
    a2 = np.hstack([zeroN,idenN*e_px,zeroN,zeroN,zeroN,zeroN])
    a3 = np.hstack([zeroN,zeroN,idenN*e_py,zeroN,zeroN,zeroN])
    a4 = np.hstack([zeroN,zeroN,zeroN,idenN*e_s,zeroN,zeroN])
    a5 = np.hstack([zeroN,zeroN,zeroN,zeroN,idenN*e_px,zeroN])
    a6 = np.hstack([zeroN,zeroN,zeroN,zeroN,zeroN,idenN*e_py])
    hamiltonian_component = np.vstack([a1,a2,a3,a4,a5,a6])
    H += hamiltonian_component
    return H

def Hhop(vertexdata): #hopping energy
    

def main():
    """
    r1 = np.array([[1],[3]])
    r2 = np.array([[2],[4]])
    print( r1.T * r2)
    return 0
    """

if __name__ == "__main__":
    main()