#STATIC CASE
#implementation of a Markov chain model demonstrating Zipf's law

import math
import numpy as np
from scipy.linalg import expm
import random
import matplotlib.pyplot as plt

#number of quantities
n=1000
max_quant=1000
grid_size=100
p_ii=0

#inital quantities as integers
q = [np.array([random.randint(0,max_quant) for i in range(n)])]

#generate random 2d locations for cities and compute distances
x = [(grid_size*random.random(),grid_size*random.random()) for i in range(n)]
d = [[math.dist(x[i],x[j]) for j in range(n)] for i in range(n)]

#transition probabilities are given by p_ij = 1/d_ij where d_ij
#is the distance from node i to node j.
P = [[(1/math.dist(x[i],x[j]) if i != j else p_ii) for j in range(n)] for i in range(n)]
row_sum = [sum(row) for row in P]
#normalize so that row sums (outgoing) are 1
P = np.array([[P[i][j]/row_sum[i] for j in range(n)] for i in range(n)])
#make sure row sums are 1 within numerical error
for row in P:
    assert(sum(row) - 1 < 10**-9)

#compute eigenvalues, eigenvectors, and stationary distribution
eig_val, eig_vec = np.linalg.eig(P)
pi = eig_vec[0]/sum(eig_vec[0])

#for cont. case, subtract 1 from diagonal to get zero row sum
#A is stochastic matrix
A = P - np.identity(n)

#to solve, simply exponentiate the matrix
def Q(t):
    return np.matmul(q[0],expm(t*A))

#sort list in ascending order to get rank-freq. dist
quantities=sorted(list(Q(100))).reverse()
plt.plot(quantities)

#plot the initial distribution for reference
init=sorted(q[0]).reverse()
plt.plot(init)

plt.show()
plt.clf()
