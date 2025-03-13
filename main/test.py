import numpy as np
import pylab as plt

L = 32
W = 3

K = np.diag(np.ones(L-1),1)
K += K.T

hit = []
hit1 = []

for arg in range(1):

    print(arg)

    H = K + W*np.diag(np.random.uniform(-1,1,L))
    e,v = np.linalg.eigh(H)


    psi_0 = np.zeros(L,complex)
    psi_0[L//2] += 1.


    psi = np.copy(psi_0)

    data = []
    data1 = []

    lst = [2**x for x in np.arange(-5,20,0.5)]
    vect = -1j*e
    print(np.exp(vect))

    print(np.exp(vect[0]))

    for t in [2**x for x in np.arange(-5,20,0.5)]:

        #print(t)
        psi = np.dot(np.conjugate(v).T,psi_0)
        psi = np.exp(-1j*e*t)*psi
        #print(np.exp(-1j*e*t))