1. Option_prices_and_demo_of_replication.R

This code does the following:

(a) simulates the discounted payoffs of calls and puts in the presence of lower reflecting barrier, and shows that     these converge to the analytical formulas in the paper 

(b) demonstrates that replication of a written put in the presence of the barrier always works, using EITHER the 
Black-Scholes delta, OR the barrier formula delta (the latter is preferred because it is cheaper).

This code is vectorised (i.e. all nSim simulations are performed in a single matrix, with no loop), and so runs over 50x faster than earlier unvectorised code.


2. Direct_and_synthetic_replication_forwards_puts_calls.R

This code simulates replication paths as in Figure 6 in the paper. 

You can then use the graph-drawing file at (3.) below to produce the graphs for a single asset specimen path.


3. Draw_replication_portfolios_paths.R

Set RunNo = the run number from (2.) above that you want to graph.

To reproduce Figure 6 in the paper, graph Run No 3 after running (2.) with these parameters:

S = K = 1, b = 0.5, tau = 25, r = 0.015, q = 0.1,  sigma = 0.13, c = -0.0001, h =1, 
seed = 1930, N = 50,400, nSim = 10.

