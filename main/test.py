import numpy as np
import pylab as plt

L = 8
W = 3

K = np.diag(np.ones(L-1),1)
K += K.T

hit = []
hit1 = []
hit2 = []
hit3 = []
hit4  =[]

for arg in range(1):
    print(arg)

    H = K + W*np.diag(np.random.uniform(-1,1,L))
    e,v = np.linalg.eigh(H)


    psi_0 = np.zeros(L,complex)
    psi_0[L//2] += 1.


    psi = np.copy(psi_0)
    #############Entaglament Entropy###############
    ######We are studying a system of free fermions and everything is encoded in the one-body density matrix C
    ### We can campute the half-partition entaglement entropy in the few line codes.

    #######Starting state, (1,0,1,0----)

    Identity = np.eye(L) #np.zeros([L,L//2], complex)

    C = Identity[:,1::2]

    print(C)


    data4 = []

    for t in [2**x for x in np.arange(-5,-4,0.5)]:

          C1 = np.dot(np.conjugate(v).T, C)
          C1  = np.dot(np.diag(np.exp(-1j*e*t)), C1)
          C1 = np.dot(v, C1)

          C_Final = np.dot(C1,np.conjugate(C1).T)

          print(len(C_Final), len(C_Final[0]))
          print(C_Final)
          print("________________________")

          ######### C_Final is the evolved density matrix at time t
          ##### now we slice it to have any the one for half-partition
          
          C_half = C_Final[:L//2,:L//2]
          print(len(C_half), len(C_half[0]))
          print("________________________")
          print("________________________")





          eigenvals = np.linalg.eigvalsh(C_half)

          eigenvals1 = eigenvals[eigenvals>10**(-20)] 
          eigenvals2 = np.ones(L//2) - eigenvals
          
          eigenvals2 = eigenvals2[eigenvals2>10**(-20)] 

          ####### Entaglement entropy = -\sum_{a} s_a\log{s_a} -\sum_a(1-s_a)\log{1-s_a) where s_a are the eigenvalues of C_half

          S = -np.dot(eigenvals1, np.log(eigenvals1)) -np.dot(eigenvals2,np.log(eigenvals2)) 
          print(S)


          data4.append(S)

    hit4.append(data4)
