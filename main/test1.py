import numpy as np
import pylab as plt

L = 4
W = 3

K = np.diag(np.ones(L-1),1)
K += K.T

hit = []
hit1 = []
hit2 = []
hit3 = []
for arg in range(1):

    print(arg)

    H = K + W*np.diag(np.random.uniform(-1,1,L))
    e,v = np.linalg.eigh(H)


    psi_0 = np.zeros(L,complex)
    psi_0[L//2] += 1.


    psi = np.copy(psi_0)

    data = []
    data1 = []

    for t in [2**x for x in np.arange(-5,20,0.5)]:

        psi = np.dot(np.conjugate(v).T,psi_0)
        psi = np.exp(-1j*e*t)*psi
        psi = np.dot(v,psi)

        data.append(np.abs(psi[L//2])**2)
        data1.append(np.dot(np.arange(-L//2,L//2)**2,np.abs(psi)**2) - np.dot(np.arange(-L//2,L//2),np.abs(psi)**2)**2 )

    hit.append(data)
    hit1.append(data1)

    ###First Method

    #############INBALANCE###############
    C = np.diag([x%2 for x in range(L)]) #Initial density matrix
    data2 = []
    for t in [2**x for x in np.arange(-5,-4,0.5)]:
        
          C1 = np.dot(np.conjugate(v).T, C)
          C1  = np.dot(np.diag(np.exp(-1j*e*t)), C1)
          C1 = np.dot(v, C1)

          C1 = np.dot(C1, v)
          C1 = np.dot(C1, np.diag(np.exp(1j*e*t)))
          C1 = np.dot(C1, np.conjugate(v).T)
          DD = np.diag(C1)
          data2.append(np.sum(DD[1::2])-np.sum(DD[::2]))

    hit2.append(data2)

