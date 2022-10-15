#=====================#
# DISCLAIMER
# This  code is provided "as is" and "with all faults" for educational and research use only.  
#   The author makes no representations or warranties of any kind concerning its correctness or 
#   suitability for any particular purpose.  Some numbers and graphs in published papers were
#   produced in Excel, not this code. If you find errors or ways to improve the code, please 
#   let me know at R.G.ThomasNOSPAM AT kent.ac.uk
#=====================#

#This code:
#(1) evaluates call and put prices analytically, and by Monte Carlo simulation 
#(2) demonstrates #that replication of a written put in the presence of the barrier always works, using EITHER the 
#Black-Scholes delta, OR the barrier formula delta (the latter is preferred because it is cheaper).

#This code is vectorised (i.e. all nSim asset paths are simulated in a single matrix, with no loop), 
#and so runs over 50x faster than earlier unvectorised code. 

# Setting parameters 
S <- 1.0 #stock price at time t
K <- 1.0 #strike price 
b <- 0.5 #barrier - for b = K case, set b just slightly less than K, to avoid divisions by zero problems.
tau <- 25 #time to maturity T - t (in years) 
r <- 0.015 #risk-free annual interest rate (convenient to set to zero)
q <- 0.01 #deferment rate (yield) (needs slightly different from r, to avoid division-by-zero problems in theta)
sigma <- 0.13 #annual volatility of the notional GBM price (standard deviation)

set.seed(1930) #set the seed (if not set, it's taken from the computer clock)
N <- 6300 #N is number of time steps, e.g. 252 x 25 = 6300 steps for every working day over 25 years. 
nSim <- 100 #number of simulations (paths) 


#Check validity of inputs
stopifnot(b <= min(S,K), r!=q, K!=b)

#analytic prices
# z's as in the papers
z1 <- (log(S/K) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z2 <- (log(b^2/(K*S)) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z3 <- (log(S/b) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z4 <- (log(b/S) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
theta <- 2*(r-q)/sigma^2


BS_call <- S*exp(-q*tau)*pnorm(z1) - K*exp(-r*(tau))*pnorm(z1-sigma*sqrt(tau)) #BS value call
BS_put <- -S*exp(-q*tau)*pnorm(-z1) + K*exp(-r*tau)*pnorm(-z1+sigma*sqrt(tau)) #BS value put

Analytic_barrier_call <- S*exp(-q*tau)*pnorm(z1) - K*exp(-r*(tau))*pnorm(z1-sigma*sqrt(tau))+
1/theta*(S*exp(-q*tau)*(b/S)^(1+theta)*pnorm(z2) - K*exp(-r*(tau))*(K/b)^(theta-1)*pnorm(z2-theta*sigma*sqrt(tau)))

Analytic_barrier_put <- K*exp(-r*tau)*pnorm(-z1+sigma*sqrt(tau)) - S*exp(-q*tau)*pnorm(-z1)-
b*exp(-r*tau)*pnorm(-z3+sigma*sqrt(tau)) + S*exp(-q*tau)*pnorm(-z3)+
1/theta*(b*exp(-r*tau)*pnorm(-z3+sigma*sqrt(tau)) - S*exp(-q*tau)*(b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)) - 
           K*exp(-r*tau)*(K/b)^(theta-1)*pnorm(z2-theta*sigma*sqrt(tau)))


# Compute the Monte Carlo prices 

dt <- tau/N #length of each time sub interval
Z <-  matrix(rnorm(nSim*N, mean=0, sd=1),nrow = nSim, ncol = N) #standard normal sample of N elements
dW <- Z*sqrt(dt) #Brownian motion increments (nSim simulations) x (N increments)
X <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_as_matrix <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_T <- vector(length=nSim) 

X[,1] <- S
Y[,1] <- S
b_as_matrix[,] <- b
B_running_max[,1] <- b_as_matrix[,1]/X[,1]


Option_payoff_written_put <- numeric(nSim)

Final_forward_hedge_position <- numeric(nSim)
Final_forward_profit <- numeric(nSim)
Forward_hedging_error <- numeric(nSim)

BS_final_cash <-numeric(nSim)
BS_final_stock <-numeric(nSim)
BS_hedging_error <- numeric(nSim)

Thomas_final_cash <-numeric(nSim)
Thomas_final_stock <-numeric(nSim)
Thomas_hedging_error <- numeric(nSim)


#Forward replicatione 
Forward_hedge_delta <- 1 
# Forwardhedge is static! So this is trivial, but it's included here to counter the common (wrong) 
# intuition that the barrier changes the replication cost of a forward.
 
# Forward_hedge_position <- numeric(N)
Forward_hedge_position <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
#Black-Scholes replication of barrier put  
  # NB In this program we are using Black-Scholes to hedge the put with reflecting barrier - NOT the ordinary 
  #  put without a reflecting barrier. This always works, for the reasons in Appendix B of the paper.
  # Another heuristic way of thinking about it: Black-Scholes replicates the put payoff for ANY random GBM path, 
  #  including one which happens to always turn upwards at say 0.5 (but without a barrier in place 
  #  at 0.5). The RGBM path (ie with an actual barrier at 0.5) is observationally indistinguishable 
  #  from this. 
BS_delta_written_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_before_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
BS_delta_written_put[,1] <- -exp(-q*tau)*(pnorm(z1)-1) #positive number
BS_stock_position_after_trade[,1] <- BS_delta_written_put[1] 
BS_cash_account[,1] <-  -S * BS_delta_written_put[1] - BS_put 
BS_trade_size[,1] <- 0
  
#Thomas put replication

Thomas_delta_written_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
Thomas_z1 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z2 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z3 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z4 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_delta_written_put[,1] <- -exp(-q*tau)*(pnorm(z1) - pnorm(z3) + (b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)))
#using z1,...,z4 from top of program...OK at time 1.
Thomas_stock_position_after_trade[,1] <- Thomas_delta_written_put[,1] # same as above
Thomas_cash_account[,1] <- -S * Thomas_delta_written_put[,1] - Analytic_barrier_put 
#negative number, borrowing to buy the asset. 
Thomas_trade_size[,1] <- 0
  
# Then do the loop of N time steps, with hedging at every step. This is slow, but given the recursive nature of the 
# definition of Y_t and the need to hedge at every step, I can't see how to speed it up.
 
for(j in 2:N){

    X[,j] <- X[,j-1]*exp((r - q -0.5*sigma^2)*dt + sigma*dW[,j]) 
    B[,j] <- b_as_matrix[,j]/X[,j]
    B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
    #Done as a running maxinum to avoid calculating MAX(B[,1;j]) in every interation of the loop
    Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
    
   #First forward replication

    Forward_hedge_position[,j] <- Y[,j] * Forward_hedge_delta 
    # trivial hedge, which gets ALL the benefit of the interventions. 
    # The hedge neutralises the drift, but NOT the interventions - correct for a forward!
    
   #Then for BS replication
    
    BS_delta_written_put[,j] <- -exp(-q*tau*(1-j/N))*(pnorm((log(Y[,j]/K) +(r- q +
                                 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N))))-1)
    
    
    BS_trade_size[,j] <- BS_delta_written_put[,j] - BS_delta_written_put[,j-1] 
    #if delta has gone down, need to sell some  
    BS_stock_position_after_trade[,j] <- BS_stock_position_after_trade[,j-1] + BS_trade_size[,j] 
    # trade restores position to fractional size delta, based on new asset price    
    #E.g. Selling some if price has gone up, which implies delta of a written put has fallen
    BS_cash_account[,j] <- BS_cash_account[,j-1] * exp(r*dt) + BS_stock_position_after_trade[,j-1] *Y[,j] * q *dt - 
                           BS_trade_size[,j] * Y[,j] 
    
    #Then for Thomas replication
    
    Thomas_z1[,j] <-(log(Y[,j]/K) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z2[,j] <-(log(b^2/(K*Y[,j])) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z3[,j] <-(log(Y[,j]/b) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z4[,j] <-(log(b/Y[,j]) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    
    
    if (j!=N) {
      Thomas_delta_written_put[,j] <- -exp(-q*tau*(1-j/N))*(pnorm(Thomas_z1[,j]) - pnorm(Thomas_z3[,j]) +
                                       (b/Y[,j])^(1+theta) * (pnorm(Thomas_z4[,j]) - pnorm(Thomas_z2[,j])))
    } else {
      Thomas_delta_written_put[,j] <- ifelse(K > Y[,j], 1, 0)
    }
    # Setting the final delta manually avoids the deep bug of z2, z3, z4 evaluating as 0/0 if Y_T[i,j]] = b 
    #   and time remaining = 0. 
    
    Thomas_trade_size[,j] <- Thomas_delta_written_put[,j] - Thomas_delta_written_put[,j-1]
    Thomas_stock_position_after_trade[,j] <- Thomas_stock_position_after_trade[,j-1] + Thomas_trade_size[,j] 
    Thomas_cash_account[,j] <- Thomas_cash_account[,j-1] * exp(r*dt) + 
      Thomas_stock_position_after_trade[,j-1] *Y[,j] * q *dt - Thomas_trade_size[,j] * Y[,j]
    
   } # end of N steps
  
Y_T <- Y[,N]
  
Option_payoff_written_put <- - pmax(K - Y_T,0) # minus sign as _written_ put. pmax is "parallel maximum" of 2 vectors.
  
#First for forward replication 

Final_forward_hedge_position <- Forward_hedge_position[,j] 
Final_forward_profit <- Y_T - Final_forward_hedge_position # alwqays zero!
Forward_hedging_error <- Final_forward_profit #always zero!
  
#Then for BS replication
  
BS_final_cash <- BS_cash_account[,j] 
BS_final_stock <- BS_stock_position_after_trade[,j] 
BS_hedging_error <- ((BS_final_cash + BS_final_stock*Y_T) - Option_payoff_written_put)  

#Then for Thomas replication
  
Thomas_final_cash <- Thomas_cash_account[,j] 
Thomas_final_stock <- Thomas_stock_position_after_trade[,j] 
Thomas_hedging_error <- ((Thomas_final_cash + Thomas_final_stock * Y_T) - Option_payoff_written_put) 

payoff_expiry_call <-pmax(Y_T-K,0) 
expected_payoff_call <- mean(payoff_expiry_call)
Monte_Carlo_call_price <- exp(-r*(tau))*expected_payoff_call

payoff_expiry_put <-pmax(K-Y_T,0) 
expected_payoff_put <- mean(payoff_expiry_put)
Monte_Carlo_put_price <- exp(-r*(tau))*expected_payoff_put

#First means for forward

Mean_forward_hedging_error <- mean(Forward_hedging_error) 
# Forward hedging error is always zero - I include this only to show that the usual forward hedging recipe works
#   perfectly. (People often think the barrier will change the forward price. But think: if you introduce a barrier
#   where there was none before, that will raise both spot and forward **in £ numeraire**; but the spot-to-forward
#   arbitrage is unaffected.)

#Then means for BS 
Mean_BS_final_cash <-mean(BS_final_cash)
Max_BS_final_cash <- max(BS_final_cash)
Min_BS_final_cash <- min(BS_final_cash)
Mean_BS_hedging_error <- mean(abs(BS_hedging_error))

#Then means for Thomas 

Thomas_mean_final_cash <-mean(Thomas_final_cash)
Thomas_max_final_cash <- max(Thomas_final_cash)
Thomas_min_final_cash <- min(Thomas_final_cash)
Thomas_mean_hedging_error <- mean(abs(Thomas_hedging_error))

Ratio_Call <- Monte_Carlo_call_price/Analytic_barrier_call
Ratio_Put <- Monte_Carlo_put_price/Analytic_barrier_put

hist(Forward_hedging_error)# always zero!
hist(BS_hedging_error) 
hist(Thomas_hedging_error) 

cat("Analytic B-Scholes Call :",round(BS_call, digits=5), "   Analytic B-Scholes Put :",round(BS_put, digits=5))

cat("\nAnalytic Barrier Call   :",round(Analytic_barrier_call, digits=5),
    "   Analytic Barrier Put   :",round(Analytic_barrier_put, digits=5))

cat("\nMonte-Carlo Barrier Call:",round(Monte_Carlo_call_price, digits=5), 
    " Monte-Carlo Barrier Put   :",round(Monte_Carlo_put_price, digits=5))

cat("\n    MC/Analytic for Call:",round(Ratio_Call, digits=5), 
    "        MC/Analytic for Put:",round(Ratio_Put, digits=5))

cat("\nMean B-S abs replication error (written put with barrier):",round(Mean_BS_hedging_error,digits=5),
    "Thomas ditto: ",round(Thomas_mean_hedging_error,digits=5))

#=====================#

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#For the GNU General Public License, see <https://www.gnu.org/licenses/>.

