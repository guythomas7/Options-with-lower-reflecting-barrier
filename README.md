1. Option_prices_and_demo_of_replication.R

This is the main code, and most likely to be useful for other people's initial experimentation. The code does the following:

(a) simulates the payoffs of calls and puts (discounted at the risk-free rate) in the presence of a lower reflecting barrier, and shows that these converge to the analytical formulas in the paper 

(b) demonstrates that replication of a written put (e.g. a NNEG) in the presence of the barrier always works, using EITHER barrier formula delta (preferred because it is cheaper), OR the Black-Scholes delta (as per Appendix B of the paper).

This code is vectorised (i.e. all nSim simulations are performed in a single matrix, with no loop), and so runs over 50x faster than earlier unvectorised code. The matrix contains nSim rows (simulation paths) x N columns (time steps).


2. Direct_and_synthetic_replication_forwards_puts_calls.R

This code simulates replication paths as in Figure 6 in the paper. 

You can then use the graph-drawing file at (3.) below to produce the call, put and forward paths for a single underlying asset path.


3. Draw_replication_portfolios_paths.R

Set RunNo = the run number from (2.) above that you want to graph.

To reproduce Figure 6 in the paper, I set RunNo = 3, after running (2.) with these parameters:

S = K = 1, b = 0.5, tau = 25, r = 0.015, q = 0.1,  sigma = 0.13, c = -0.0001, h =1, 
set.seed = 1930, N = 50,400, nSim = 10.

I'm using R 3.02 and 64-bit Windows 10. I'm not sure if the same seed will give the same result with different versions or hardware. 
