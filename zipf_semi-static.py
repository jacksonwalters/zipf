#SEMI-STATIC CASE
#p_ij = boole(Q_j - Q_i > 0)/d_ij
#meaning transitions from i->j only when population at Q_j exceeds Q_i
#starting with slow solution - updating transition matrix each time

import math
import numpy as np
from scipy.linalg import expm
import random
import matplotlib.pyplot as plt

#number of quantities
n=5

#inital quantities as integers
q = [np.array([random.randint(0,1000) for i in range(n)])]

#generate random 2d locations for cities and compute distances
x = [(100*random.random(),100*random.random()) for i in range(n)]
d = [[math.dist(x[i],x[j]) for j in range(n)] for i in range(n)]

#transition probabilities are given by p_ij = 1/d_ij where d_ij
#is the distance from node i to node j.
def p(i,j,t):
    return bool(q[t][i]-q[t][j] < 0)/math.dist(x[i],x[j]) if i != j else 1
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
