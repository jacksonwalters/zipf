#rigorous Markov chain model demonstrating Zipf's law

Zipf's "law" is a rank-frequency distribution which
which occurs in various settings, e.g. for a large
corpus of natural language text, count the number of times each
word occurs (its frequency, f), and order by this
value to obtain an index (the rank, k), and plot f vs. k.
this data is well approximated by a power law f(k)=1/k^s
normalizing to get a probability distribution, we get
f(k)=1/(H_N,S * k^s) where H_N,s is the N-th generalized
harmonic number sum_{n=1}^{N} 1/n^s where N is the number
of words in the corpus. in the limit as N -> infty, this
is the Riemann zeta function.

Zipf's law also shows up in naturally occurring measures like populations of cities. To come up
with a model that explains the resulting distribution, let's assume that people are more likely
to move to a city if either

1) it is nearby
2) it is more populated

Viewing the cities as n nodes in a network (city_1,city_2,...,city_n), the probability to move from city_1 to city_2 as p_12,...,p_ij,...,p_nn (the transition matrix), we can write down a Markov-chain model where each step in the process is viewed as unit of time t_1,...,t_max. Assuming the cities lie in a 2d plane (more generally an n-dim'l vector space), we can compute the symmetric distance (or similarity) matrix d_ij. Let the population of a city (or quantity) be f_1,...,f_n.

CASE I (STATIC): p_ij does not depend on time. A person is more likely to move to a city that is closer, so

p_ij = 1/d_ij.

CASE II (SEMI-STATIC): p_ij depends on time, but is constant on intervals where the order of the
populations is constant. That is, we only care whether the population (Q_i) of one city is greater than another,
not by how much (which is changing with time, Q_i = Q_i(t)).

Let boole(statement) = {0 if statement is False, 1 if statement is True}. Then

p_ij = boole(Q_i-Q_j > 0)/d_ij

CASE III (DYNAMIC): pi_j depends on both the population difference and the distance between cities,

p_ij = (Q_i - Q_j)/d_ij

SETTING UP THE MARKOV CHAIN:
----------------------------

Since p_ij is the probability of a transition/move from city_i to city_j, the amount of people
we expect to move out of city_i and into city_j is Q_i*p_ij.

That is, we need to sum over all incoming (+) and outgoing (-) cities/nodes to get
the total population change between times t_r and t_r+1:

Q_i(t_r+1) = Q_i(t_r) + sum_k Q_k(t_r)*p_ki - sum_k Q_i(t_r)*p_ik

We can solve for this process by diagonalizing the matrix. The eigenvalues will tell you how fast things are changing, and the long term behavior will show that places with large quantities tend to attract resources from places with less resources.
