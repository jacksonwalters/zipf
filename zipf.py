#implementation of a Markov chain model demonstrating Zipf's law

import math
import numpy as np
import scipy
import random
import matplotlib.pyplot as plt

#number of quantities
n=10

#inital quantities as integers
Q_0 = [random.randint(0,1000) for i in range(n)]

#generate random 2d locations for cities and compute distances
x = [(100*random.random(),100*random.random()) for i in range(n)]
d = [[math.dist(x[i],x[j]) for j in range(n)] for i in range(n)]

#show plot of distances between nodes
plt.scatter(*zip(*x))
plt.show()

#STATIC, CONTINUOUS CASE
#transition probabilities are given by p_ij = 1/d_ij where d_ij
#is the distance from node i to node j. normalized so that row sums
#(outgoing) are 1
p = [[(1/math.dist(x[i],x[j]) if i != j else 0) for j in range(n)] for i in range(n)]
row_sum = [sum(row) for row in p]
p = [[p[i][j]/row_sum[i] for j in range(n)] for i in range(n)]

print([sum(row) for row in p])
