# options

This code:

(1) simulates the discounted payoffs of calls and puts in the presence of lower reflecting barrier, and shows that these converge to the analytical formulas in the paper 

(2) demonstrates that replication of a written put in the presence of the barrier always works, using EITHER the 
Black-Scholes delta, OR the barrier formula delta (the latter is preferred because it is cheaper).

This code is vectorised (i.e. all nSim simulations are performed in a single matrix, with no loop), and so runs over 50x faster than earlier unvectorised code.

