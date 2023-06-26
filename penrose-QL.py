import numpy as np
import matplotlib.pyplot as plt
import math
import Penrose_tile02

vertexdata = Penrose_tile02.pen()

def Hsite(vertexdata):#on site energy
    e_s = 1.8
    e_px = -6.5
    e_py = -6.5
    
    N = len(vertexdata)
    zeroN = np.zeros((N,N))*0j
    idenN = np.identity(N)
    H = np.zeros((6*N,6*N))*0j

    a1 = np.hstack([idenN*e_s,zeroN,zeroN,zeroN,zeroN,zeroN])
    a2 = np.hstack([zeroN,idenN*e_px,zeroN,zeroN,zeroN,zeroN])
    a3 = np.hstack([zeroN,zeroN,idenN*e_py,zeroN,zeroN,zeroN])
    a4 = np.hstack([zeroN,zeroN,zeroN,idenN*e_s,zeroN,zeroN])
    a5 = np.hstack([zeroN,zeroN,zeroN,zeroN,idenN*e_px,zeroN])
    a6 = np.hstack([zeroN,zeroN,zeroN,zeroN,zeroN,idenN*e_py])
    hamiltonian_component = np.vstack([a1,a2,a3,a4,a5,a6])
    H += hamiltonian_component #6N * 6N
    return H

def hopparameter(j,bond1,bond2): #bond1,2: 1->s  2->px  3->py
    theta = 2/5*math.pi
    V_sssig = -0.4
    V_spsig = 0.9
    V_ppsig = 1.8
    V_pppi = 0.05
    l = math.cos(j*theta)
    m = math.sin(j*theta)
    t = 0
    if bond1 == 1 and bond2 == 1:   t = V_sssig
    elif bond1 == 1 and bond2 == 2: t = l * V_spsig
    elif bond1 == 1 and bond2 == 3: t = m * V_spsig
    elif bond1 == 2 and bond2 == 2: t = l*l*V_ppsig + (1-l*l)*V_pppi
    elif bond1 == 3 and bond2 == 3: t = m*m*V_ppsig + (1-m*m)*V_pppi
    elif bond1 == 2 and bond2 == 3: t = l*m*(V_ppsig - V_pppi)
    return t 


def Hhop(vertexdata): #hopping energy
    N = len(vertexdata)
    H = np.zeros((6*N,6*N))*0j
    for alpha in range(3):
        for beta in range(3):
            for i in range(N):
                neighbor_vec = vertexdata[i]["neighbor"]
                num_neighbor = len(neighbor_vec)
                for j in range(num_neighbor):
                    t = hopparameter(neighbor_vec[j][1],alpha+1,beta+1)
                    anni_pos = np.zeros((3*N,1))*0j
                    create_pos = np.zeros((3*N,1))*0j
                    create_pos[N*alpha+i] += 1
                    anni_pos[N*beta+neighbor_vec[j][0]]+= 1
                    hamiltonian_component = create_pos * anni_pos.conj().T
                    a1 = np.hstack([hamiltonian_component , np.zeros((3*N,3*N))*0j])
                    a2 = np.hstack([np.zeros((3*N,3*N))*0j , hamiltonian_component])
                    H += np.vstack([a1,a2])
    return H

def Hsoc(vertexdata):
    lam = 0.8
    N = len(vertexdata)
    H = np.zeros((6*N,6*N))*0j
    for i in range(N):
        create_vec = np.zeros((3*N,1))*0j
        anni_vec = np.zeros((3*N,1))*0j
        create_vec[2*N+i] += 1
        anni_vec[N+i] += 1
        h = 1j * lam * create_vec * anni_vec.conj().T
        hamiltonian_component = h + h.conj().T
        a1 = np.hstack([hamiltonian_component , np.zeros((3*N,3*N))*0j])
        a2 = np.hstack([np.zeros((3*N,3*N))*0j , hamiltonian_component*(-1)])
        H += np.vstack([a1,a2])
    return H
        
def Hamiltonian(vertexdata):
    H = Hsite(vertexdata) + Hhop(vertexdata) + Hsoc(vertexdata)
    return H

def eigval(vertexdata,eigval,eigvec): 
    N = len(vertexdata)
    eigenergy = []
    for i in range(len(eigval)):
        eigenergy.append(eigval[i].real)
    eigenergy.sort()
    index = list(range(1,6*N+1))
    fig = plt.figure()
    plt.plot(index,eigenergy)
    plt.savefig('/Users/sogenikegami/Documents/UT4S/non-crystal/penrose/penrose-QL_eigval')
    fig.show()


def wavefunc_plot(vertexdata,eigval,eigvec):
    N = len(vertexdata)
    R=[]
    theta = []
    wavefunc = []
    for i in range(6*N):
        if abs(eigval[i].real) < 0.005:
            for j in range(6*N):
                wavefunc.append(eigvec[j][i])
                R.append(abs(eigvec[j][i]))
                theta.append(math.atan2(wavefunc[j].imag,wavefunc[j].real))
            break
    
    fig = plt.figure()
    for i in range(N):

        plt.scatter(vertexdata[i]["pos"][0],vertexdata[i]["pos"][1],s=10)
    for i in range(6):
        for j in range(N):
            plt.scatter(vertexdata[j]["pos"][0],vertexdata[j]["pos"][1],marker='.',s=R[i*N+j]*300,c=theta[i*N+j],cmap='jet')
    plt.savefig('/Users/sogenikegami/Documents/UT4S/non-crystal/penrose/penrose-QL_wavefunc03')
    fig.show()
    
    #print(len(R))
    #print(len(theta))
    return 0
            

    

def main():
    H = Hamiltonian(vertexdata)
    eigval , eigvec = np.linalg.eig(H)
    wavefunc_plot(vertexdata,eigval,eigvec)
    #eigval(vertexdata)


if __name__ == "__main__":
    main()