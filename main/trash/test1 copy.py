import scipy.sparse as sparse
import numpy as np
#import scipy.sparse.linalg.eigen.arpack as arp
import pylab as pl
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import expm
import time
import pylab as plt
from tqdm import tqdm




def gen_spin_operators(L):
        """" Returns the spin operators sigma_x and sigma_z for L sites """

        sx = sparse.csr_matrix(np.array([[0.,1.],[1.,0.]]))
        sz = sparse.csr_matrix(np.array([[1.,0.],[0.,-1.]]))
        sy = sparse.csr_matrix(np.array([[0.,1j],[-1j,0.]]))


        d = 2

        sx_list = []
        sz_list = []
        sy_list = []


        for i_site in range(L):
                if i_site==0:
                        X=sx
                        Z=sz
                        Y=sy

                else:
                        X= sparse.csr_matrix(np.eye(d))
                        Z= sparse.csr_matrix(np.eye(d))
                        Y= sparse.csr_matrix(np.eye(d))

                for j_site in range(1,L):
                        if j_site==i_site:
                                X=sparse.kron(X,sx, 'csr')
                                Z=sparse.kron(Z,sz, 'csr')
                                Y=sparse.kron(Y,sy, 'csr')

                        else:
                                X=sparse.kron(X,np.eye(d),'csr')
                                Z=sparse.kron(Z,np.eye(d),'csr')
                                Y=sparse.kron(Y,np.eye(d),'csr')

                sx_list.append(X)
                sz_list.append(Z)
                sy_list.append(Y)
        return np.array(sx_list),np.array(sy_list), np.array(sz_list)

def entanglement_entropy(psi,L):
        " Calculate the entanglement entropy "
        block_dim =  2**(L//2)
        psi_block = np.reshape(psi,(block_dim,block_dim))
        s = np.linalg.svd(psi_block,compute_uv=0)
        s=s[s>10**(-20)]**2

        return -np.inner(np.log(s),s)

#@njit
def int_to_spin_list(state_int, length):
    # Convert integer to binary string, padded to the desired length
    binary_str = bin(state_int)[2:].zfill(length)
    # Map binary digits to spin list (-1 for 0, +1 for 1)
    spin_list = [-1 if bit == '0' else 1 for bit in binary_str]
    return spin_list

#@njit
def spin_list_to_int(spin_list):
    # Map -1 to 0, +1 to 1
    bit_list = [0 if spin == -1 else 1 for spin in spin_list]
    # Convert bit list to integer
    state_int = 0
    for bit in bit_list:
        state_int = (state_int << 1) | bit
    return state_int

def Connected_state_with_M(length):
    # Convert integer to binary string, padded to the desired length

    AA = []
    for kk in range(2**L):
        spin_list = int_to_spin_list(kk, L)
        if np.sum(spin_list)== 0.:
            AA.append(kk)

    return AA

def gen_hamiltonian(sp_list,sm_list, sz_list,L):
    H = sparse.csr_matrix((2**L,2**L))
    H_int = sparse.csr_matrix((2**L,2**L))

    s_p = sp_list
    s_m = sm_list

    for i in range(L):
        
        H += sp_list[i]*sm_list[np.mod(i+1,L)] 
        H_int += sz_list[i]*sz_list[np.mod(i+1,L)]
    return H+H.T, H_int

def gen_hamiltonian1(h, sz_list,L):
    H = sparse.csr_matrix((2**L,2**L))

    for i in range(L):    
        H += h[i]*sz_list[i]

    return H

L = 12

s_x, s_y, s_z = gen_spin_operators(L)

sp_list = s_x+1j*s_y
sm_list = s_x-1j*s_y

H, H_int = gen_hamiltonian(sp_list,sm_list, s_z,L) ## hopping and interaction
W = 8 #disorder strenght 
hit = []
Delta = 0.001 #interaction 

M_0 = Connected_state_with_M(L) ## only states with magnetization = 0--> L/2 number of particle

for arg in tqdm(range(1000)):


    h = W*np.random.uniform(-1,1,L)
    H_d = gen_hamiltonian1(h, s_z,L)

    H_tot = H + H_d + Delta*H_int
 
    M_0 = Connected_state_with_M(L)
    M_0 = np.array(M_0)

    H_tot = H_tot[M_0, :]
    H_tot = H_tot[:,M_0]

    ### random product state

    psi = np.zeros(len(M_0), complex)
    xx = spin_list_to_int([1,-1]*(L//2))
    yy = np.where(M_0==xx)[0]
    #print(yy,xx)
    psi[yy] += 1

   
    e,v = np.linalg.eigh(H_tot.toarray())

    data = []

    for t in [2**x for x in np.arange(-5,30, 0.25)]:
        psi1 = np.dot(np.conjugate(v).T,psi)
        psi1 = np.exp(-1j*e*t)*psi1
        psi1 = np.dot(v, psi1)
        Psi = np.zeros(2**L, complex)
        Psi[M_0] += psi1
        EE = entanglement_entropy(Psi,L)
        data.append(EE)

    hit.append(data)

plt.plot([2**x for x in np.arange(-5,30, 0.25)],np.mean(hit,0), "o-")
plt.semilogx()

plt.show()

