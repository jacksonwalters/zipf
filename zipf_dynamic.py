#DYNAMIC CASE
#p_ij = boole(Q_j - Q_i > 0)/d_ij
#meaning transitions from i->j only when population at Q_j exceeds Q_i
#starting with slow solution - updating transition matrix each time

import math
import numpy as np
from scipy.linalg import expm
import random
import matplotlib.pyplot as plt

#number of quantities
n=1000
max_quant=100
grid_size=100
p_ii=10

#inital quantities as integers
q = [np.array([random.randint(0,max_quant) for i in range(n)])]

#generate random 2d locations for cities and compute distances
x = [(grid_size*random.random(),grid_size*random.random()) for i in range(n)]
d = [[math.dist(x[i],x[j]) for j in range(n)] for i in range(n)]

#transition probabilities are given by p_ij = 1/d_ij where d_ij
#is the distance from node i to node j.
def p(i,j,t):
    #diagonal element p_ii determines proportion of quantity "held"
    return abs(q[t][i]-q[t][j])/math.dist(x[i],x[j]) if i != j else p_ii
def P(t):
    #probabilities at time t
    P_t = [[p(i,j,t) for j in range(n)] for i in range(n)]
    #normalize so that row sums (outgoing) are 1
    row_sum = [sum(row) for row in P_t]
    P_t = np.array([[P_t[i][j]/row_sum[i] for j in range(n)] for i in range(n)])
    return P_t

#to solve (for now) multiply matrices
#faster to find transposition points and run cont. until those
#times are hit
def Q(t):
    #if t=0, return inital quantities
    if t == 0:
        return q[0]
    #if t is exactly length of q, just reference list
    if t == len(q):
        q_t = np.matmul(q[t-1],P(t-1))
        q.append(q_t)
        return q[t]
    #otherwise compute recursively
    if t > len(q):
        q_t = np.matmul(Q(t-1),P(t-1))
        q.append(q_t)
        return q[t]
    #if already computed, just return it
    if t < len(q):
        return q[t]

#sort list in ascending order to get rank-freq. dist
quantities=sorted(list(Q(10))).reverse()
plt.plot(quantities)

#plot the initial distribution for reference
init=sorted(list(q[0])).reverse()
plt.plot(init)

plt.show()
plt.clf()
