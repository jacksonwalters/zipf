#implementation of a Markov chain model demonstrating Zipf's law

import math
import numpy as np
from scipy.linalg import expm
import random
import matplotlib.pyplot as plt

#number of quantities
n=1000

#inital quantities as integers
Q_0 = [random.randint(0,1000) for i in range(n)]

#generate random 2d locations for cities and compute distances
x = [(100*random.random(),100*random.random()) for i in range(n)]
d = [[math.dist(x[i],x[j]) for j in range(n)] for i in range(n)]

#STATIC CASE
#transition probabilities are given by p_ij = 1/d_ij where d_ij
#is the distance from node i to node j.
p = [[(1/math.dist(x[i],x[j]) if i != j else 0) for j in range(n)] for i in range(n)]
row_sum = [sum(row) for row in p]
#normalize so that row sums (outgoing) are 1
p = np.array([[p[i][j]/row_sum[i] for j in range(n)] for i in range(n)])
#make sure row sums are 1 within numerical error
for row in p:
    assert(sum(row) - 1 < 10**-9)

#compute eigenvalues, eigenvectors, and stationary distribution
eig_val, eig_vec = np.linalg.eig(p)
pi = eig_vec[0]/sum(eig_vec[0])

#for cont. case, subtract 1 from diagonal to get zero row sum
#to get stochastic matrix
A = p - np.identity(n)

#to solve, simply exponentiate the matrix
def Q(t):
    return np.matmul(Q_0,expm(t*A))

quants=list(Q(100))
quants.sort()
quants.reverse()
plt.plot(quants)

Q_0.sort()
Q_0.reverse()
plt.plot(Q_0)

plt.show()
plt.clf()
