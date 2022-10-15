#=====================#
# DISCLAIMER
# This  code is provided "as is" and "with all faults" for educational and research use only.  
#   The author makes no representations or warranties of any kind concerning its correctness or 
#   suitability for any particular purposes.  Some numbers and graphs in published papers were
#   produced in Excel, not this code. If you find errors or ways to improve the code, please 
#   let me know at  R.G.ThomasNOSPAM AT kent.ac.uk
#=====================#

# Setting parameters
# To get Figure 6 in paper, use the following parameters here, then set RunNo = 3
#in "Graphing_of_replicating_portfolios.R":

# S = K = 1, b = 0.5, tau = 25, r = 0.015, q = 0.1,  sigma = 0.13, c = -0.0001, h =1, 
# seed = 1930, N = 50,400, nSim = 10.

S <- 1.0 #stock price at time t
K <- 1.0 #strike price  # DOESN'T MATTER for Net-delta (calls) strategy - K (or X) doesn't appear in formula
b <- 0.5 #barrier - for b = K case, set b just slightly less than K, to avoid divisions by zero problems.
tau <- 25 #time to maturity T - t (in years) 
r <- 0.015 #risk-free annual interest rate (convenient to set to zero)
q <- 0.01 #deferment rate (yield) (needs slightly different from r, to avoid division-by-zero problems in theta)
sigma <- 0.13 #annual volatility of the notional GBM price (standard deviation)

c <-  -0.0001 # to investigate limit on trading close to barrier: cannot trade when Y - b < c. 
#It turns out c has only marginal effect
#If don't want limit, set c to tiny neg number is safe (tiny discrepancies otherwise, prob from "Y = b exactly" cases).  
h <- 1 # hedging frequency - hedge every h time steps. Normally I just use 1.

set.seed(1930) #set the seed (if not set, it's taken from the computer clock)
N <- 50400 #N is number of time steps, e.g. 252 x 25 = 6300 steps for every working day over 25 years. 
nSim <- 10 #number of simulations (paths) 

#Check validity of inputs
stopifnot(b <= min(S,K), r!=q, K!=b)

#analytic prices
# z's as in the paper
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
1/theta*(b*exp(-r*tau)*pnorm(-z3+sigma*sqrt(tau)) - S*exp(-q*tau)*(b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)) - K*exp(-r*tau)*(K/b)^(theta-1)*pnorm(z2-theta*sigma*sqrt(tau)))

Analytic_forward <- S * exp(-q * tau) - K * exp(-r *tau)

# z's for S = b, for lower limit for C-dashed
z1_dash <- (log(b/K) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z2_dash <- (log(b^2/(K*b)) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z3_dash <- (log(b/b) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))
z4_dash  <- (log(b/b) + (r - q + 0.5*sigma^2)*tau)/(sigma*sqrt(tau))


Lower_limit_for_C_dashed <- b * exp(-q * tau) - K * exp(-r *tau) + 
    K*exp(-r*tau)*pnorm(-z1_dash+sigma*sqrt(tau)) - b*exp(-q*tau)*pnorm(-z1_dash)-
    b*exp(-r*tau)*pnorm(-z3_dash+sigma*sqrt(tau)) + b*exp(-q*tau)*pnorm(-z3_dash)+
    1/theta*(b*exp(-r*tau)*pnorm(-z3_dash+sigma*sqrt(tau)) - 
    b*exp(-q*tau)*(b/b)^(1+theta)*(pnorm(z4_dash)-pnorm(z2_dash)) - 
    K*exp(-r*tau)*(K/b)^(theta-1)*pnorm(z2_dash-theta*sigma*sqrt(tau)))

  
# Compute the Monte Carlo prices 

dt <- tau/N #length of each time sub interval
Z <-  matrix(rnorm(nSim*N, mean=0, sd=1),nrow = nSim, ncol = N) #standard normal sample of N elements
dW <- Z*sqrt(dt) #Brownian motion increments (nSim simulations) x (N increments)
X <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
b_as_matrix <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
B_running_max <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Y_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1) 
Y_T <- vector(length=nSim) 
Y_T_log_change <- vector(length=nSim)
Intrinsic_value_of_call <-matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Intrinsic_value_of_put <-matrix(numeric(nSim*N), nrow = nSim, ncol = N)
                         
X[,1] <- S
Y[,1] <- S
b_as_matrix[,] <- b
B_running_max[,1] <- b_as_matrix[,1]/X[,1]

Option_payoff_bought_put <- numeric(nSim)

Option_exercise_payment_bought_put <- numeric(nSim)


BS_final_cash <-numeric(nSim)
BS_final_stock <-numeric(nSim)
BS_hedging_error <- numeric(nSim)

Thomas_final_cash <-numeric(nSim)
Thomas_final_stock <-numeric(nSim)
Thomas_hedging_error <- numeric(nSim)

Thomas_call_final_cash <-numeric(nSim)
Thomas_call_final_stock <-numeric(nSim)
Thomas_call_hedging_error <- numeric(nSim)

Lowest_difference_replicating_portfolios <- numeric(nSim)
Lowest_difference_forward_replicating_portfolios <- numeric(nSim)


ND_final_stock <- numeric(nSim)
ND_final_cash <- numeric(nSim)
ND_lowest_cash <- numeric(nSim)
ND_lowest_portfolio_value <- numeric(nSim)
ND_lowest_LTV <- numeric(nSim)

negative_cases <- matrix(nrow = 0, ncol=5)
barrier_touches <- matrix(nrow = 0, ncol=3)
Y_saved <-matrix(nrow=nSim, ncol=N)
ND_cash_account_saved <-matrix(nrow=nSim, ncol=N)
ND_portfolio_value_saved <-matrix(nrow=nSim, ncol=N)
ND_delta_most_recent_permitted_trade_saved <- matrix(nrow=nSim, ncol=N)
ND_LTV_saved <-matrix(nrow=nSim, ncol=N)


#Black-Scholes replication of barrier put    

BS_delta_bought_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_delta_bought_put_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_before_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
BS_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

BS_delta_bought_put[,1] <- exp(-q*tau)*(pnorm(z1)-1) #negative number
BS_stock_position_after_trade[,1] <- BS_delta_bought_put[1] 
BS_cash_account[,1] <-  -S * BS_delta_bought_put[1] + BS_put 
BS_delta_bought_put_most_recent_permitted_trade[,1] <-BS_delta_bought_put[,1] 
#opening the position always permitted
BS_replicating_portfolio[,1] <- BS_stock_position_after_trade[,1] * Y[,1] + BS_cash_account[,1]
BS_trade_size[,1] <- 0
  
#Thomas put replication 
Thomas_delta_bought_put <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_put_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Thomas_synthetic_call_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Difference_put_replicating_portfolios <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_z1 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z2 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z3 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_z4 <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
  
Thomas_delta_bought_put[,1] <- exp(-q*tau)*(pnorm(z1) - pnorm(z3) + (b/S)^(1+theta)*(pnorm(z4)-pnorm(z2)))
Thomas_stock_position_after_trade[,1] <- Thomas_delta_bought_put[,1] # same as above
Thomas_cash_account[,1] <- - S * Thomas_delta_bought_put[,1] + Analytic_barrier_put 
#first term is positive, becasue for a bought put we go short the asset. 
Thomas_replicating_portfolio[,1] <- Thomas_stock_position_after_trade[,1] * Y[,1] + 
                                   Thomas_cash_account[,1]
Thomas_delta_bought_put_most_recent_permitted_trade[,1] <-Thomas_delta_bought_put[,1] 
#opening the position always permitted
Thomas_trade_size[,1] <- 0

Thomas_synthetic_call_replicating_portfolio[,1] <-  Thomas_stock_position_after_trade[,1] *Y[,1]+
  Thomas_cash_account[,1] +
  Y[,1]*exp(-q*tau) - K*exp(-r*tau)


#Thomas call replication

Thomas_delta_bought_call <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_call_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_put_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_put_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Difference_call_replicating_portfolios <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

#Forward F' and synthetic forward F, replication

Forward_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Forward_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Synthetic_forward_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Synthetic_forward_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Difference_forward_replicating_portfolios <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_delta_bought_call <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_delta_bought_call_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_synthetic_put_replicating_portfolio <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_put_replicating_portfolio_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
Interim_rolled_up_difference_in_analytic_call_prices <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices <- 
  matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call <- 
  matrix(numeric(nSim*N), nrow = nSim, ncol = N)
Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

Thomas_delta_bought_call[,1] <- exp(-q*tau)*(pnorm(z1) - (b/S)^(1+theta)*pnorm(z2))
#using z1, z2 from top of program...OK at time 1.
Thomas_call_stock_position_after_trade[,1] <- Thomas_delta_bought_call[,1] # same as above
Thomas_call_cash_account[,1] <- -S * Thomas_delta_bought_call[,1] + Analytic_barrier_call 
#Negative number, as for BS, for low barrier. Can be positive for high barrier (when
#the delta of the Thomas call gets low).
Thomas_call_replicating_portfolio[,1] <- Thomas_call_stock_position_after_trade[,1] * Y[,1] + 
  Thomas_call_cash_account[,1]
Thomas_delta_bought_call_most_recent_permitted_trade[,1] <- Thomas_delta_bought_call[,1] 
#opening the position always permitted
Thomas_call_trade_size[,1] <- 0

Thomas_synthetic_put_replicating_portfolio[,1] <-  Thomas_call_replicating_portfolio[,1] -
                                                     (Y[,1]*exp(-q*tau) - K*exp(-r*tau))

Difference_put_replicating_portfolios[,1] <- Thomas_synthetic_put_replicating_portfolio[,1] - 
                                                Thomas_replicating_portfolio[,1]  

Interim_rolled_up_difference_in_analytic_call_prices[,1] <- Analytic_barrier_call - 
                                                           (Analytic_forward + Analytic_barrier_put)

Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices[,1] <-
  Thomas_synthetic_call_replicating_portfolio[,1] /  Interim_rolled_up_difference_in_analytic_call_prices[,1]


#Forward F and synthetic forward F*, replication
Forward_replicating_portfolio[,1] <- Y[,1]* exp(-q * tau) - K * exp(-r*tau)
Synthetic_forward_replicating_portfolio[,1] <- Thomas_call_replicating_portfolio[,1] -
                 Thomas_replicating_portfolio[,1]
Difference_forward_replicating_portfolios[,1] <- Synthetic_forward_replicating_portfolio[,1] -
                                                     Forward_replicating_portfolio[,1]
                                   

#Net delta (ND)
ND_delta <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_delta_most_recent_permitted_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_stock_position_before_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_stock_position_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_trade_size <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_cash_account <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_portfolio_value <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_stock_value_after_trade <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)
ND_stock_value_log_change <- matrix(numeric(nSim*(N-1)), nrow = nSim, ncol = N-1)
ND_LTV <- matrix(numeric(nSim*N), nrow = nSim, ncol = N)

ND_delta[,1] <- exp(-q*tau)*(1 - pnorm(z3) + (b/S)^(1+theta)*(pnorm(z4)))
ND_stock_position_after_trade[,1] <- ND_delta[,1]
ND_cash_account[,1] <-  -S *ND_delta[,1] # DON'T put initial put valuation here - just
#doing the delta hedging WITHOUT replicating the two offsetting strikes.
ND_portfolio_value[,1] <-  ND_cash_account[,1] + ND_stock_position_after_trade[,1] *Y[,1]
ND_stock_value_after_trade[,1] <-  ND_stock_position_after_trade[,1] *Y[,1]
ND_delta_most_recent_permitted_trade[,1] <- ND_delta[,1]
ND_trade_size[,1] <- 0

  
#Then do the loop - with hedging every h time steps
  
for(j in 2:N){
    
    X[,j] <- X[,j-1]*exp((r - q -0.5 * sigma^2)*dt + sigma*dW[,j]) 
    B[,j] <- b_as_matrix[,j]/X[,j]
    B_running_max[,j] <- ifelse(B[,j] > B_running_max[,j-1], B[,j], B_running_max[,j-1])
    Y[,j] <- X[,j]*pmax(1,B_running_max[,j])
    Y_log_change[,j-1] <- log(Y[,j]) - log(Y[,j-1]) 
    
    #First for BS put hedging 
   
      # If(j %% h > 0) branch (j modulo h) says: if we're not at a hedging day, 
      #  just maintain the existing cash and stock positions.
      
    
      if(j %% h > 0) {
        BS_cash_account[,j] <- BS_cash_account[,j-1] * exp(r*dt) + BS_stock_position_after_trade[,j-1] *Y[,j] * q *dt 
        BS_stock_position_after_trade[,j] <- BS_stock_position_after_trade[,j-1]
        BS_delta_bought_put_most_recent_permitted_trade[,j] <- BS_delta_bought_put_most_recent_permitted_trade[,j-1]
        
      }else{
        
        BS_delta_bought_put[,j] <- exp(-q*tau*(1-j/N))*(pnorm((log(Y[,j]/K) + 
                                     (r- q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N))))-1)
        
        BS_trade_size[,j] <- ifelse(Y[,j] - b < c, 0, BS_delta_bought_put[,j] - 
                                      BS_delta_bought_put_most_recent_permitted_trade[,j-1])
        BS_stock_position_after_trade[,j] <- BS_stock_position_after_trade[,j-1] + BS_trade_size[,j] 
        BS_cash_account[,j] <- BS_cash_account[,j-1] * exp(r*dt) + 
                                 BS_stock_position_after_trade[,j-1] *Y[,j] * q *dt - BS_trade_size[,j] * Y[,j] 
        BS_delta_bought_put_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b < c, 
                                                                BS_delta_bought_put_most_recent_permitted_trade[,j-1],
                                                                 BS_delta_bought_put[,j]) 
        BS_replicating_portfolio[,j] <- BS_stock_position_after_trade[,j] *Y[,j] + BS_cash_account[,j]
        
      }    


    #Then for Thomas put hedging
    
    Thomas_z1[,j] <-(log(Y[,j]/K) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z2[,j] <-(log(b^2/(K*Y[,j])) +(r - q +0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z3[,j] <-(log(Y[,j]/b) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
    Thomas_z4[,j] <-(log(b/Y[,j]) +(r - q + 0.5*sigma^2)*(tau*(1-j/N)))/(sigma*sqrt(tau*(1-j/N)))
 
    if(j %% h > 0) {
      Thomas_cash_account[,j] <- Thomas_cash_account[,j-1] * exp(r*dt) + 
                                   Thomas_stock_position_after_trade[,j-1] *Y[,j] * q *dt
      Thomas_stock_position_after_trade[,j] <- Thomas_stock_position_after_trade[,j-1]
      Thomas_delta_bought_put_most_recent_permitted_trade[,j] <- Thomas_delta_bought_put_most_recent_permitted_trade[,j-1]
    
       }else{
    
           if (j!=N) {
           Thomas_delta_bought_put[,j] <- exp(-q*tau*(1-j/N))*(pnorm(Thomas_z1[,j]) - pnorm(Thomas_z3[,j]) +
             (b/Y[,j])^(1+theta) * (pnorm(Thomas_z4[,j]) - pnorm(Thomas_z2[,j])))
           } else {
           Thomas_delta_bought_put[,j] <- ifelse(K < Y[,j], -1, 0)
          
            }
    
       Thomas_trade_size[,j] <- ifelse(Y[,j] - b < c, 0, Thomas_delta_bought_put[,j] - 
                                         Thomas_delta_bought_put_most_recent_permitted_trade[,j-1])
       Thomas_stock_position_after_trade[,j] <- Thomas_stock_position_after_trade[,j-1] + Thomas_trade_size[,j]
       Thomas_cash_account[,j] <- Thomas_cash_account[,j-1] * exp(r*dt) + 
                                     Thomas_stock_position_after_trade[,j-1] *Y[,j] * q *dt -
                                     Thomas_trade_size[,j] * Y[,j]
    
       Thomas_delta_bought_put_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b < c, 
             Thomas_delta_bought_put_most_recent_permitted_trade[,j-1], Thomas_delta_bought_put[,j])
    
       #Now set the delta manually if at last step
       if (j!=N) {
         Thomas_replicating_portfolio[,j] <- Thomas_stock_position_after_trade[,j] *Y[,j] +
           Thomas_cash_account[,j]
       }else{
         Thomas_replicating_portfolio[,j] <-  Thomas_delta_bought_put[,j]*Y[,j] +
           Thomas_cash_account[,j]
       }
       
       
    if (j!=N) {
      Thomas_synthetic_call_replicating_portfolio[,j] <- Thomas_stock_position_after_trade[,j]*Y[,j] +
        Thomas_cash_account[,j] + Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N))
    }else{
      Thomas_synthetic_call_replicating_portfolio[,j] <-  Thomas_delta_bought_put[,j]*Y[,j] +
        Thomas_cash_account[,j] + Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N))
      
    }
    # set the delta manually if at last step
    
    Interim_rolled_up_difference_in_analytic_call_prices[,j] <- 
         (Analytic_barrier_call - (Analytic_forward + Analytic_barrier_put))*exp(r * tau * j/N)
    
    Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices[,j] <-
        Thomas_synthetic_call_replicating_portfolio[,j] /  Interim_rolled_up_difference_in_analytic_call_prices[,j]
    
    Intrinsic_value_of_call[,j] <- pmax(Y[,j]- K,0)
    Intrinsic_value_of_put[,j] <- pmax(K - Y[,j],0)
    
    Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call[,j] <-
      Thomas_synthetic_call_replicating_portfolio[,j] - Intrinsic_value_of_call[,j]
      
    }
    
     #Then for Thomas call hedging 
    
    # Already calculated the required column of z_1[,j] and column of z_2[,j] above

    if(j %% h > 0) {
      Thomas_call_cash_account[,j] <- Thomas_call_cash_account[,j-1] * exp(r*dt) + Thomas_call_stock_position_after_trade[,j-1] *Y[,j] * q *dt
      Thomas_call_stock_position_after_trade[,j] <- Thomas_call_stock_position_after_trade[,j-1]
      Thomas_delta_bought_call_most_recent_permitted_trade[,j] <- Thomas_delta_bought_call_most_recent_permitted_trade[,j-1]    
     
       }else{
      
      
          if (j!=N) {
          Thomas_delta_bought_call[,j] <- exp(-q*tau*(1-j/N))*(pnorm(Thomas_z1[,j]) -
                                                             (b/Y[,j])^(1+theta) * pnorm(Thomas_z2[,j]))
           } else {
          Thomas_delta_bought_call[,j] <- ifelse(K < Y[,j], 1, 0)
          }
    # Setting the final delta manually avoids the deep bug of z2, z3, z4 evaluating as 0/0 if Y_T[i,j]] = b 
    #   and time remaining = 0. 
    
       
         Thomas_call_trade_size[,j] <- ifelse(Y[,j] - b < c, 0, Thomas_delta_bought_call[,j] - Thomas_delta_bought_call_most_recent_permitted_trade[,j-1])
         Thomas_call_stock_position_after_trade[,j] <- Thomas_call_stock_position_after_trade[,j-1] + Thomas_call_trade_size[,j]
         Thomas_call_cash_account[,j] <- Thomas_call_cash_account[,j-1] * exp(r*dt) + Thomas_call_stock_position_after_trade[,j-1] *Y[,j] * q *dt -
                                          Thomas_call_trade_size[,j] * Y[,j]
         
         Thomas_delta_bought_call_most_recent_permitted_trade[,j] <- ifelse(Y[,j] - b < c, 
                 Thomas_delta_bought_call_most_recent_permitted_trade[,j-1], Thomas_delta_bought_call[,j])        
         
         #Now set the delta manually if at last step
         if (j!=N) {
         Thomas_call_replicating_portfolio[,j] <- Thomas_call_stock_position_after_trade[,j] *Y[,j] +
                                                      Thomas_call_cash_account[,j]
         }else{
          Thomas_call_replicating_portfolio[,j] <-  Thomas_delta_bought_call[,j]*Y[,j] +
                                                       Thomas_call_cash_account[,j]
          }
         
         if (j!=N) {
           Thomas_synthetic_put_replicating_portfolio[,j] <- Thomas_call_replicating_portfolio[,j] -
                                                               (Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N)))
         }else{
           Thomas_synthetic_put_replicating_portfolio[,j] <-  Thomas_delta_bought_call[,j]*Y[,j] +
             Thomas_call_cash_account[,j] - (Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N)))
         } 
         #Synthetic put has to go here, where call calc has just been done
         
         Difference_put_replicating_portfolios[,j] <- Thomas_synthetic_put_replicating_portfolio[,j] - 
                                                         Thomas_replicating_portfolio[,j] 
         
         Difference_call_replicating_portfolios[,j] <- Thomas_call_replicating_portfolio[,j] - 
                                                         Thomas_synthetic_call_replicating_portfolio[,j]  
         
         
         Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call[,j] <-
           Thomas_call_replicating_portfolio[,j] - Intrinsic_value_of_call[,j]
         
         
         
       }     
    
    #Forward F and synthetic forward F*, replication
    
    Forward_replicating_portfolio[,j] <- Y[,j]*exp(-q*tau*(1-j/N)) - K * exp(-r*tau*(1-j/N))
    Synthetic_forward_replicating_portfolio[,j] <-  Thomas_call_replicating_portfolio[,j] -
                                                         Thomas_replicating_portfolio[,j]
    
    Forward_replicating_portfolio_log_change[,j-1] <- log(Y[,j]) - log(Y[,j-1])
    # only the change in the stock position counts as volatility, not the flow of yield or interest
    Synthetic_forward_replicating_portfolio_log_change[,j-1] <- 
      log(max(abs(Thomas_call_stock_position_after_trade[,j] - Thomas_stock_position_after_trade[,j]) *Y[,j], 0.0000001)) -  
      (log(max(abs(Thomas_call_stock_position_after_trade[,j-1] - Thomas_stock_position_after_trade[,j-1]) *Y[,j-1], 0.0000001)))
    #ditto  
        
    Difference_forward_replicating_portfolios[,j] <- Synthetic_forward_replicating_portfolio[,j] -
                                                         Forward_replicating_portfolio[,j]
    
    #Then for the Net Delta 

   
        if(j %% h > 0) {
          ND_cash_account[,j] <- ND_cash_account[,j-1] * exp(r*dt) + ND_stock_position_after_trade[,j-1] *Y[,j] * q *dt
          ND_stock_position_after_trade[,j] <- ND_stock_position_after_trade[,j-1]
          ND_delta_most_recent_permitted_trade[,j] <- ND_delta_most_recent_permitted_trade[,j-1]
          
        }else{
          if (j!=N) {
            ND_delta[,j] <- exp(-q*tau*(1-j/N))*(1 - pnorm(Thomas_z3[,j]) +
                                                   (b/Y[,j])^(1+theta) * (pnorm(Thomas_z4[,j])))
          } else {
            ND_delta[,j] <- 0 # The net delta at maturity must be 0
          }
        
        }   
        
        
        ND_trade_size[,j] <- ifelse(Y[,j] - b < c, 0, ND_delta[,j] - ND_delta_most_recent_permitted_trade[,j-1])
        ND_stock_position_after_trade[,j] <- ND_stock_position_after_trade[,j-1] + ND_trade_size[,j]
        ND_cash_account[,j] <- ND_cash_account[,j-1] * exp(r*dt) + ND_stock_position_after_trade[,j-1] *Y[,j] * q *dt -
                      ND_trade_size[,j] * Y[,j]
        ND_portfolio_value[,j] <- ND_cash_account[,j] + ND_stock_position_after_trade[,j]*Y[,j]
        ND_stock_value_after_trade[,j] <- ND_stock_position_after_trade[,j]*Y[,j]
        ND_stock_value_log_change[,j-1] <- 
                  log(max(ND_stock_value_after_trade[,j]*Y[,j], 0.0000001)) -
                                   log(max(ND_stock_value_after_trade[,j-1]*Y[,j-1], 0.0000001))         
        
        
        
        ND_LTV[,j] <- pmin(ND_cash_account[,j],0)/(ND_stock_position_after_trade[,j] *Y[,j])
        
        
        ND_delta_most_recent_permitted_trade[,j] <- ifelse (Y[,j] - b < c, 
                ND_delta_most_recent_permitted_trade[,j-1], ND_delta[,j])

      
} # end of N steps


#=================#
   

ND_LTV[,N] <- 0 
# Set every LTV to 0 at time N (else stock position going to zero at time N gives spuriously large LTV)
  
Y_T <- Y[,N]

Y_T_log_change <- log(Y[,N]) - log(Y[,1])
  
 
Option_payoff_bought_put <-  pmax(K - Y_T,0)# pmax is "parallel maximum" of 2 vectors.
Option_payoff_bought_call <- pmax(Y_T - K,0)
  
#First for BS replication

BS_final_cash <- BS_cash_account[,j]
BS_final_stock <-  BS_stock_position_after_trade[,j]


BS_hedging_error <- ((BS_final_cash + BS_final_stock * Y_T) - Option_payoff_bought_put)
# ("proceeds of replication scheme" - "option payoff")

#Then for Thomas bought put replication
  

Thomas_final_cash <-  Thomas_cash_account[,j]

Thomas_final_stock <-  Thomas_stock_position_after_trade[,j]

Thomas_hedging_error <- (Thomas_final_cash + Thomas_final_stock * Y_T) - Option_payoff_bought_put
  

#Then for Thomas bought call replication

Thomas_call_final_cash <- Thomas_call_cash_account[,j]
Thomas_call_final_stock <-  Thomas_call_stock_position_after_trade[,j]


Thomas_call_hedging_error <- ((Thomas_call_final_cash + Thomas_call_final_stock * Y_T) - Option_payoff_bought_call)

Lowest_forward_replicating_portfolio <- apply(Forward_replicating_portfolio,1,min)
Time_of_lowest_forward_replicating_portfolio <- apply(Forward_replicating_portfolio,1,which.min)

Lowest_Thomas_synthetic_call_replicating_portfolio <- apply(Thomas_synthetic_call_replicating_portfolio,1,min)
Time_of_lowest_Thomas_synthetic_call_replicating_portfolio <- 
  apply(Thomas_synthetic_call_replicating_portfolio,1,which.min)


Lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices <-
  apply(Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices,1,min)
Time_of_lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices <-
  apply(Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices,1,which.min)

Lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call <-
      apply(Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call,1,min)
Time_of_lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call <-
      apply(Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call,1,which.min)

Lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call <-
  apply(Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call,1,min)
Time_of_lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call <-
  apply(Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call,1,which.min)


Lowest_difference_put_replicating_portfolios <- apply(Difference_put_replicating_portfolios,1,min)
Lowest_difference_call_replicating_portfolios <- apply(Difference_call_replicating_portfolios,1,min)
Lowest_difference_forward_replicating_portfolios <- apply(Difference_forward_replicating_portfolios,1,min)

#Then for the Net Delta

ND_final_cash <- ND_cash_account[,j] + ND_stock_position_after_trade[,j] *Y[,j]


ND_lowest_cash <- apply(ND_cash_account,1,min) # returns vector: the lowest values from each row of the matrix
ND_lowest_portfolio_value <- apply(ND_portfolio_value,1,min)
Time_of_lowest_ND_portfolio_value <- apply(ND_portfolio_value,1,which.min)
ND_lowest_LTV <- apply(ND_LTV,1,min)
  

# STORING CASES WHERE BARRIER IS TOUCHED, IF WANT TO INSPECT (but it turns out there is nothing very interesting)
  
#if(sum(Y - b < 1e-6) > 0){
#
# barrier_touches <- rbind(barrier_touches, c(k, sum(Y[,]-b < 1e-6), round(ND_final_cash[k], digits=5)))
#
#  }
  
#break statement if want to stop and inspect a case - set number to the case number
#if (k+1 > 1288) break
  
#sTORING DETAILS OF NEGATIVE CASES
  
#  if (ND_final_cash[k] < 0) {
#    
#    negative_cases <- rbind(negative_cases, c(k, sum(Y-b < 1e-6), round(ND_final_cash[k], digits=5),
#                         round(ND_lowest_portfolio_value[k], digits=5), round(ND_lowest_LTV[k], digits=5)))
#    Y_saved[k,] <- Y #only way can get this to work - silly
#    ND_cash_account_saved[k,] <-ND_cash_account
#    ND_portfolio_value_saved[k,] <-ND_portfolio_value
#    ND_delta_most_recent_permitted_trade_saved[k,] <- ND_delta_most_recent_permitted_trade
#    
#    ND_LTV_saved[k,] <- ND_LTV
    
#  }
  

payoff_expiry_call <-pmax(Y_T-K,0) 
expected_payoff_call <- mean(payoff_expiry_call)
Monte_Carlo_call_price <- exp(-r*(tau))*expected_payoff_call

payoff_expiry_put <-pmax(K-Y_T,0) 
expected_payoff_put <- mean(payoff_expiry_put)
Monte_Carlo_put_price <- exp(-r*(tau))*expected_payoff_put

#First meansn for BS 
Mean_BS_final_cash <-mean(BS_final_cash)
Max_BS_final_cash <- max(BS_final_cash)
Min_BS_final_cash <- min(BS_final_cash)
Mean_BS_hedging_error <- mean(abs(BS_hedging_error))

#Then means for Thomas put
Thomas_mean_final_cash <-mean(Thomas_final_cash)
Thomas_max_final_cash <- max(Thomas_final_cash)
Thomas_min_final_cash <- min(Thomas_final_cash)
Thomas_mean_hedging_error <- mean(abs(Thomas_hedging_error))

#Then means for Thomas call
Thomas_call_mean_final_cash <-mean(Thomas_call_final_cash)
Thomas_call_max_final_cash <- max(Thomas_call_final_cash)
Thomas_call_min_final_cash <- min(Thomas_call_final_cash)
Thomas_call_mean_hedging_error <- mean(abs(Thomas_call_hedging_error))

#Then for ND
ND_mean_lowest_cash <-mean(ND_lowest_cash)
ND_min_lowest_cash <- min(ND_lowest_cash)
ND_mean_lowest_portfolio_value <- mean(ND_lowest_portfolio_value)
ND_min_lowest_portfolio_value <- min(ND_lowest_portfolio_value)
ND_mean_final_cash <- mean(ND_final_cash)
ND_mean_lowest_LTV <- mean(ND_lowest_LTV)
ND_min_lowest_LTV <- min(ND_lowest_LTV)
ND_min_final_cash <- min(ND_final_cash)
ND_count_below_zero <- sum(ND_final_cash<0)
Row_no_of_the_ND_min_lowest_portfolio_value <-  which.min(ND_lowest_portfolio_value)
Time_of_lowest_in_min_row_for_ND_portfolio_value  <-
  Time_of_lowest_ND_portfolio_value[Row_no_of_the_ND_min_lowest_portfolio_value]

Ratio_call <- Monte_Carlo_call_price/Analytic_barrier_call
Ratio_put <- Monte_Carlo_put_price/Analytic_barrier_put

Y_SD_log_change <- sqrt(1/tau * N) * mean(apply(Y_log_change,1,sd))
Thomas_replicating_portfolio_SD_log_change <- sqrt(1/tau * N)*
                                        mean(apply(Thomas_replicating_portfolio_log_change,1,sd))
Thomas_synthetic_put_replicating_portfolio_SD_log_change <- sqrt(1/tau * N)* 
                                        mean(apply(Thomas_synthetic_put_replicating_portfolio_log_change,1,sd))
Thomas_synthetic_put_multiple <- 
  Thomas_synthetic_put_replicating_portfolio_SD_log_change / Thomas_replicating_portfolio_SD_log_change


Y_T_SD_log_change <- sqrt((1/tau * N)/N)* sd(Y_T_log_change)

Mean_lowest_forward_replicating_portfolio <- mean(Lowest_forward_replicating_portfolio)
Min_lowest_forward_replicating_portfolio <- min(Lowest_forward_replicating_portfolio)
Row_no_of_the_min_forward_replicating_portfolio <-
  which.min(Lowest_forward_replicating_portfolio)
Time_of_lowest_in_min_row_for_forward_replicating_portfolio  <-
  Time_of_lowest_forward_replicating_portfolio[Row_no_of_the_min_forward_replicating_portfolio]


Mean_lowest_Thomas_synthetic_call_replicating_portfolio <- mean(Lowest_Thomas_synthetic_call_replicating_portfolio)
Min_lowest_Thomas_synthetic_call_replicating_portfolio <- min(Lowest_Thomas_synthetic_call_replicating_portfolio)
Row_no_of_the_min_Thomas_synthetic_call_replicating_portfolio <-
  which.min(Lowest_Thomas_synthetic_call_replicating_portfolio)
Time_of_lowest_Thomas_synthetic_call_replicating_portfolio  <-
  Time_of_lowest_Thomas_synthetic_call_replicating_portfolio[Row_no_of_the_min_Thomas_synthetic_call_replicating_portfolio]

Difference_with_initial_C_dashed_plus_ND_lowest_portfolio_value_for_same_row <-
  Min_lowest_Thomas_synthetic_call_replicating_portfolio -
  (ND_lowest_portfolio_value[Row_no_of_the_min_Thomas_synthetic_call_replicating_portfolio] + 
  (Analytic_forward + Analytic_barrier_put))

Mean_lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices <-
    mean(Lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices)
Min_lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices <-
    min(Lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices)
Row_no_of_the_min_ratio_for_rolled_up_difference <-
    which.min(Lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices)
Time_of_lowest_in_min_row_for_rolled_up_difference  <-
    Time_of_lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices[Row_no_of_the_min_ratio_for_rolled_up_difference]

Mean_lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call <-
  mean(Lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call)
Min_lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call <-
    min(Lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call)
Row_no_of_the_min_for_intrinsic_value <-
  which.min(Lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call)
Time_of_lowest_in_min_row_for_intrinsic_value  <-
  Time_of_lowest_Thomas_synthetic_call_replicating_portfolio_minus_intrinsic_value_of_call[Row_no_of_the_min_for_intrinsic_value]

Mean_lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call <-
  mean(Lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call)
Min_lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call <-
  min(Lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call)
Row_no_of_the_min_for_Thomas_call_intrinsic_value <-
  which.min(Lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call)
Time_of_lowest_in_min_row_for_Thomas_call_intrinsic_value  <-
  Time_of_lowest_Thomas_call_replicating_portfolio_minus_intrinsic_value_of_call[Row_no_of_the_min_for_Thomas_call_intrinsic_value]


Min_difference_put_replicating_portfolios <- min(Lowest_difference_put_replicating_portfolios)
Min_difference_call_replicating_portfolios <- min(Lowest_difference_call_replicating_portfolios)
Min_difference_forward_replicating_portfolios <- min(Lowest_difference_forward_replicating_portfolios)


#apply(matrix, margin, function) applies function to the 2margins" of the matrix, whee 1=rows, 2 = columns
#returns a vector, which I then take the mean of. Each rwo is 1 25-year simulation.

mybins <-seq(0.0, 0.1,0.001) #This can be used to make both histograms look same, with narrow bins.
mybins2 <-seq(-0.9, 0.1, 0.01) 
mybins3 <-seq(-0.6, 0, 0.01)
hist(BS_hedging_error)
hist(Thomas_hedging_error)
hist(ND_final_cash) #, breaks=mybins)
hist(ND_lowest_cash) # , breaks=mybins2)

#hist(ND_final_cash, breaks=seq(0,0.07,0.001), main="Terminal gain", xlab="", cex.axis=1.2, cex.lab=1.4, cex.main =1.8)
#hist(ND_lowest_cash, breaks=seq(-0.6,0.,0.01), main="Lowest cash", xlab="", cex.axis=1.2, cex.lab=1.4, cex.main =1.8)


#Calculate the theoretical gains from the trading - exactly offsetting the rolled-forward difference 
# in initial valuations
Rolled_up_difference_in_call_valuations <- 
  (Analytic_barrier_call - (S * exp(-q*tau) - K*exp(-r*tau) + Analytic_barrier_put)) * exp(r * tau)


Int_loss_ratio_sorted <- 
  sort(Lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices)

Row_No_of_closest_ILR_to_minus_1 <- which.min(abs(Int_loss_ratio_sorted + 1))

if(min(Int_loss_ratio_sorted) < -1){
  
  Row_No_of_closest_ILR_below_minus_1 <- Row_No_of_closest_ILR_to_minus_1
  Percentile_of_closest_ILR_below_minus_1 <- round(Row_No_of_closest_ILR_below_minus_1/nSim * 100, digits = 2)
  
}else{
  
  Row_No_of_closest_ILR_below_minus_1 <- 1
  Percentile_of_closest_ILR_below_minus_1 <- 0
}

CTE_p <- mean(Int_loss_ratio_sorted[1:Row_No_of_closest_ILR_below_minus_1])

CTE_25 <- mean(Int_loss_ratio_sorted[1:round(nSim/4, digits=0)])

cat("Analytic B-Scholes Call :",BS_call, "   Analytic B-Scholes Put :",BS_put)

cat("\nAnalytic Barrier Call   :",Analytic_barrier_call, "  Analytic Barrier Put   :",Analytic_barrier_put)

cat("\nMonte-Carlo Barrier Call:",Monte_Carlo_call_price, " Monte-Carlo Barrier Put:",Monte_Carlo_put_price)

cat("\nMean B-Scholes put abs hedging error:",round(Mean_BS_hedging_error,digits=7),  
      "Thomas ditto: ",round(Thomas_mean_hedging_error,digits=7))

cat("\nMean Thomas call abs hedging error:",round(Thomas_call_mean_hedging_error,digits=7))

cat("\nRolled-up difference in call valuations:",   Rolled_up_difference_in_call_valuations)

cat("\nLowest cash over term in net-delta strategy: MIN:",round(ND_min_lowest_cash,digits=5),  
    "MEAN:",round(ND_mean_lowest_cash,digits=5))

cat("\nLowest portfolio over term in net-delta strategy: MEAN:",round(ND_mean_lowest_portfolio_value,digits=5),
    "MIN:" ,round(ND_min_lowest_portfolio_value,digits=5),  
    "Row no of the MIN:" ,Row_no_of_the_ND_min_lowest_portfolio_value,    
    "Time of lowest in MIN row:" , Time_of_lowest_in_min_row_for_ND_portfolio_value)  

cat("\nLowest LTV over term from net-delta strategy: MIN:",round(ND_min_lowest_LTV,digits=5),
      "MEAN:",round(ND_mean_lowest_LTV,digits=5))

#Lowest LTV over term is not very meaningful - often very large near the end of the term, but this is 
# only because stock position becomes tiny. Lowest cash over term (i.e. max borrowing) is a more 
# meaningful constraint on the ND strategy.

cat("\nFinal cash from net-delta strategy: MIN:",round(ND_min_final_cash,digits=5), 
      "MEAN:",round(ND_mean_final_cash,digits=5))
    
cat("\nInitial net-delta:" ,round(ND_delta[1,1],digits=5))

cat("\nLowest forward replicating portfolio: MEAN:" ,round(Mean_lowest_forward_replicating_portfolio, digits=5),
    "MIN:" ,round(Min_lowest_forward_replicating_portfolio,digits=5),
    "Row no of the MIN:" ,Row_no_of_the_min_forward_replicating_portfolio,    
    "Time of lowest in MIN row:" , Time_of_lowest_in_min_row_for_forward_replicating_portfolio)  


cat("\nLowest synthetic call replicating portfolio: MEAN:" ,
      round(Mean_lowest_Thomas_synthetic_call_replicating_portfolio, digits=5),
      "MIN:" ,round(Min_lowest_Thomas_synthetic_call_replicating_portfolio,digits=5),
      "Row no of the MIN:" ,Row_no_of_the_min_Thomas_synthetic_call_replicating_portfolio,    
      "Time of lowest in MIN row:" , Time_of_lowest_Thomas_synthetic_call_replicating_portfolio) 

cat("\nLower limit: (F + P) at t = 0, S = b:" ,round(Lower_limit_for_C_dashed, digits=5))

#cat("\nDifference of above less (initial C' + ND_lowest_portfolio_value), for same row:" 
#      ,round(Difference_with_initial_C_dashed_plus_ND_lowest_portfolio_value_for_same_row, digits=5))
  

cat("\nLowest synthetic call replicating portfolio as multiple of interim rolled-up (C - C'):")

cat("\nMEAN:" ,round(Mean_lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices, digits=5),
    "MIN:" ,round(Min_lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices, digits=5),
    "Row no of the MIN:" ,Row_no_of_the_min_ratio_for_rolled_up_difference,    
    "Time of lowest in MIN row:" , Time_of_lowest_in_min_row_for_rolled_up_difference,
    "25th percentile:" , round(quantile(Lowest_Thomas_synthetic_call_replicating_portfolio_as_multiple_of_rolled_up_difference_in_analytic_call_prices,0.25), digits=5))  


    
cat("\nMin diff of replication portfolios (over all runs) (~0 shows one always larger):") 
    
cat("\n(P' - P):" ,round(Min_difference_put_replicating_portfolios, digits=5),
     "(C - C'):" ,round(Min_difference_call_replicating_portfolios, digits=5), 
     "(F - F'):" ,round(Min_difference_forward_replicating_portfolios, digits=5))

cat("\nPercentile of closest ILR just below -1:" , round(Percentile_of_closest_ILR_below_minus_1, digits=5), "%", 
    "  CTE of ILR  at that percentile:" , CTE_p, "  25%-CTE of ILR:" ,CTE_25)


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
