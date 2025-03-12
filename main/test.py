import numpy as np
import pylab as plt

L = 8
W = 3

K = np.diag(np.ones(L-1),1)
K += K.T

hit = []
hit1 = []

for arg in range(1):

    print(arg)

    H = K + W*np.diag(np.random.uniform(-1,1,L))
    e,v = np.linalg.eigh(H)

    print(H)