import numpy as np
import pylab as plt

L = 64
W = 2

K = np.diag(np.ones(L-1),1)
K += K.T

hit = []
hit1 = []

for arg in range(1):

    print(arg)

    H = np.array([
    [(1.75607 + 0j), (1 + 0j), (0 + 0j), (0 + 0j)],
    [(1 + 0j), (-0.173636 + 0j), (1 + 0j), (0 + 0j)],
    [(0 + 0j), (1 + 0j), (3.36858 + 0j), (1 + 0j)],
    [(0 + 0j), (0 + 0j), (1 + 0j), (-0.150378 + 0j)]])
    e,v = np.linalg.eigh(H)

    #print(H)


    psi_0 = np.zeros(L,complex)
    psi_0[L//2] += 1.


    psi = np.copy(psi_0)

    data = []
    data1 = []
    i=0
    for t in [2**x for x in np.arange(-5,20, 0.5)]:

        psi = np.dot(np.conjugate(v).T,psi_0)
        psi = np.exp(-1j*e*t)*psi
        psi = np.dot(v,psi)

        data.append(np.abs(psi[L//2])**2)
        data1.append(np.dot(np.arange(-L//2,L//2)**2,np.abs(psi)**2) - np.dot(np.arange(-L//2,L//2),np.abs(psi)**2)**2 )

        print(str(i)+":")
        print(t)
        print(np.exp(-1j*e*t))
        print(np.abs(psi[L//2])**2)
        print(np.dot(np.arange(-L//2,L//2)**2,np.abs(psi)**2) - np.dot(np.arange(-L//2,L//2),np.abs(psi)**2)**2 )
        print("________________")

        i+=1
    hit.append(data)
    hit1.append(data1)
     