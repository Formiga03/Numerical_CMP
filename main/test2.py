import numpy as np
import pylab as plt

L = 8
W = 2

K = np.diag(np.ones(L-1),1)
K += K.T

hit = []
hit1 = []

for arg in range(1):

    print(arg)

    vct = np.array([2.3, 4.3, 1.2, 1.12, 1.23, 3.4, 3.16, 2.111], dtype=complex)

    # Define the subdiagonal and superdiagonal
    aux1 = np.full(len(vct) - 1, complex(1, 0))

    # Create the tridiagonal matrix
    ham = np.diag(vct) + np.diag(aux1, k=-1) + np.diag(aux1, k=1)
    print(ham)
    e,v = np.linalg.eigh(ham)

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
        print( "1st SD:"+str(np.dot(np.arange(-L//2,L//2)**2,np.abs(psi)**2)))
        print( "2nd SD:"+str(np.dot(np.arange(-L//2,L//2),np.abs(psi)**2)**2))
        print("________________")

        i+=1
    hit.append(data)
    hit1.append(data1)
     