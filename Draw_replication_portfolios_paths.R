
RunNo <- 3
# Enter Run number you want to graph

linetypes <- c(1,2,3)


par(mfrow=c(2,2), oma = c(1,1,1,1) + 0.2,
    mar = c(3,3,1,1) + 1, mgp=c(1,1,0))
    
    
#SPOT ONLY:

matplot(cbind(Y[RunNo,1:N],b,0), # tempstrike),
        col=c("black", "chocolate4", "black", "black"), type="l", lty=c(1,2,3,5), lwd=0.1, xaxt = "n",
        main = "Spot", xlab="Time steps", ylab="Price\n", cex.main = 1.4, cex.axis = 1.2, cex.lab = 1.2, xlim=c((0/52*N),(52/52*N))) #, ylim=c(0.0, 0.05))  

text(2/52*N ,b+0.04,"barrier", cex=1.2, col="chocolate4")
axis(1, at=c(1,N), cex.axis = 1.2, labels = c("0", expression(italic("T"))))



#FORWARD PAYOFF:


matplot(cbind(Synthetic_forward_replicating_portfolio[RunNo,1:N],
              Forward_replicating_portfolio[RunNo,1:N],
              0), #b,
        col=c("darkolivegreen1", "springgreen4", "black", "black"), type="l", lty=linetypes, lwd=0.1, cex.axis=1.2, cex.main = 1.4, xaxt="n",
        main = "Forward", xlab="Time steps", ylab="Price\n", xlim=c((0/52*N),(52/52*N))) #, ylim=c(0.0, 0.05))  #-0.8,-0.1 0.80, 0.82

axis(1, at=c(1,N), cex.axis = 1.2, labels = c("0", expression(italic("T"))))


#CALL PAYOFF:


matplot(cbind(Thomas_call_replicating_portfolio[RunNo,1:N],
              Thomas_synthetic_call_replicating_portfolio[RunNo,1:N],
              #Intrinsic_value_of_call[RunNo,1:N],
               0), # b,
        col=c("pink", "red", "black", "black"), type="l", lty=linetypes, lwd=0.1, cex.axis = 1.2, cex.main = 1.4, xaxt="n",
        main = "Call", xlab="Time steps", ylab="Price\n") #, xlim=c((20/52*N),40/52*N)) # ylim=c(-0.02, 0.02))  #-0.8,-0.1 0.80, 0.82

axis(1, at=c(1,N), cex.axis = 1.2, labels = c("0", expression(italic("T"))))


#PUT PAYOFF:


matplot(cbind(Thomas_synthetic_put_replicating_portfolio[RunNo,1:N],
              Thomas_replicating_portfolio[RunNo,1:N],
              #BS_replicating_portfolio[RunNo,1:6300],
              #Intrinsic_value_of_put[RunNo,1:N],
              0), #b,
        col=c("lightcyan3", "blue", "black", "black"), type="l", lty=linetypes, lwd=0.1, cex.axis=1.2, cex.main = 1.4, xaxt="n",
        main = "Put", xlab="Time steps", ylab="Price\n") #,xlim=c((0/52*N),(52/52*N))) #, ylim=c(0.0, 0.05))  #-0.8,-0.1 0.80, 0.82

axis(1, at=c(1,N), cex.axis = 1.2, labels = c("0", expression(italic("T"))))

