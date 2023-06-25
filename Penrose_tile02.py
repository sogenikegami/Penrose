import math
import matplotlib.pyplot as plt
import numpy as np
import os

theta = 2*math.pi/5
g1 = 0.1
g2 = 0.7
g3 = -0.98
g4 = 0.43
g5 = -(g1+g2+g3+g4) 
gamma = [g1,g2,g3,g4,g5]



def d(j): #return (2,1)
    d = np.array([math.cos(j*theta),
                  math.sin(j*theta)])
    return d

def K(j,r): #return Kj input:j and position r(2,1)
    return round(np.dot(d(j),r) + gamma[j])

def project(K): #return(2,1)
    R2d = np.zeros(2)
    for j in range(5):
        R2d[0] += K[j]*d(j)[0]
        R2d[1] += K[j]*d(j)[1]
    return R2d

def findpen(r,rvalue,Kvalue,vertexvalue,Kset,vertexdata):
    Kxy = np.zeros(5)
    for k in range(5):
        Kxy[k] = K(k,r)

    flag = True
    for k in Kset:
        local_flag = True
        for l in range(5):
            if Kxy[l] != k[l]:
                local_flag = False
                break
        if local_flag:
            flag = False
            break

    if flag:
        R2d=project(Kxy)
        rvalue.append(r)
        Kvalue.append(Kxy)
        vertexvalue.append(R2d)
        Kset.append(Kxy)
        vertexdata.append({'pos':R2d,'neighbor':[],'K':Kxy,'r':r})
    return 0

def pen():
    xnum = 100
    ynum = 100
    xmax = 3 #
    ymax = 3
    radius = 3

    fig = plt.figure()
    rvalue = [] # position on 2D plane which were already searched
    vertexvalue = [] # penrose tile vertex position on 2D
    Kvalue = [] # lattice point in 5D space which were already searched
    Kset = [] #same as above
    vertexdata = []

    """
    pairs = [[0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    dr = 1/10000
    kmax = 100
    rs = np.array([[0,0],[dr,0],[0,dr],[-dr,0],[0,-dr]])

    for pair in pairs:
        i1 = pair[0]
        i2 = pair[1]
        d1 = d(i1)
        d2 = d(i2)
        M = [[math.cos(i1*theta),math.sin(i1*theta)],
             [math.cos(i2*theta),math.sin(i2*theta)]]
        Minv = np.linalg.inv(M)
        for k1 in range(-kmax,kmax):
            for k2 in range(-kmax,kmax):
                r = Minv @ np.array([k1-gamma[i1],k2-gamma[i2]])
                for ir in range(5):
                    r2 = r + rs[ir]
                    if np.linalg.norm(r2) < rmax:
                        zero = findpen(r,rvalue,Kvalue,vertexvalue,Kset)
    """
    for i in range(-xnum,xnum):
        for j in range(-ynum,ynum):
            x = i*xmax/xnum
            y = j*ymax/ynum
            if x*x + y*y < radius *radius:
                r = np.array([x,y])
                findpen(r,rvalue,Kvalue,vertexvalue,Kset,vertexdata)                
                if abs(y) < ymax/3:
                    for alpha in range(1,10):
                        for beta in range(1,10):
                            x = (i+alpha/10)*xmax/xnum
                            y = (j+beta/10)*ymax/ynum
                            r = np.array([x,y])
                            findpen(r,rvalue,Kvalue,vertexvalue,Kset,vertexdata)                     


    c='blue'
    jcount = 0
    for i in range(len(vertexdata)):
        r = vertexdata[i]['r']
        R2d = vertexdata[i]['pos']
        Kxy = vertexdata[i]['K']
        for j in range(5):
            R2dj = R2d + d(j)
            for k in range(len(vertexdata)):
                trialvec = vertexdata[k]['pos']
                if np.linalg.norm(trialvec-R2dj) < 10**(-4):
                    vertexdata[i]['neighbor'].append([k,j])
                    jcount += 1
                    plt.plot([R2d[0],R2dj[0]],[R2d[1],R2dj[1]],color=c)
                    break
            #jhop = next(filter(lambda x: np.linalg.norm(x-R2dj,ord=2) < 10**(-3),vertexvalue),"None")
           
            #if jhop != "None":


        for j in range(5):
            R2dj = R2d-d(j)
            for k in range(len(vertexdata)):
                trialvec = vertexdata[k]['pos']
                if np.linalg.norm(trialvec-R2dj) < 10**(-4):
                    vertexdata[i]['neighbor'].append([k,j])
                    jcount += 1
                    plt.plot([R2d[0],R2dj[0]],[R2d[1],R2dj[1]],color=c)
                    break
            #jhop = next(filter(lambda x: np.linalg.norm(x-R2dj,ord=2)< 10**(-3),vertexvalue),"None")
            #if jhop != "None":


    fig.show()
    plt.savefig('/Users/sogenikegami/Documents/UT4S/non-crystal/penrose/penrose06')
    print(len(vertexvalue))
    return vertexdata




def main():
    pen()


if __name__ == "__main__":
    main()