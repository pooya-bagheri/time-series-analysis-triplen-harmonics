require(astsa)

#Building fundamental current time series
#The Data is first placed in clipboard 
Data=read.table("clipboard")
FC=ts(Data,start=1/60,frequency=60)
colnames(FC) <- "Amperes"
#Building the triplen harmonic current time series
#The Data is placed in clipboard 
Data2=read.table("clipboard")
THC=ts(Data2,,start=1/60,frequency=60)
colnames(THC) <- "Amperes"

#Saving timeseries Data
save(FC,THC,file="InputData.saved")

#Loading timeseries Data
load("InputData.saved")

#Plotting original Time Series
win.graph()
par(mfrow=c(2,1))
plot(FC,main="FC (Fundamental Current)",lwd=1.5)
plot(THC,main="THC (Triplen Harmonic Current)",lwd=1.5)

#Looking at ACF of Triplen Harmonic current:
win.graph()
acf(THC,120,main="THC")

#Looking at scatter plot and ccf for relationship between FC and THC:
win.graph()
par(mfrow=c(2,1))
plot(FC,THC,main="THC(Triplen Harmonic Current) vs. FC(Fundamental Current)")
ccf(drop(FC),drop(THC),main="CCF: THC vs. FC")

#Detrending triplen harmonic current by Regression on Fundamental
library(dynlm)  
fit = dynlm(THC ~ FC)
summary(fit)

#Investigating the residuals of time series after regression
THC_regressed=0.06199*FC+0.39704
THC_resid=THC-THC_regressed
colnames(THC_resid) <- "Amperes"
win.graph()
par(mfrow=c(2,1))
plot(THC_resid,main="Residuals of Regressing THC on FC",lwd=1.5)
acf(THC_resid,main="Residuals of Regressing THC on FC",120)

win.graph()
par(mfrow=c(2,1))
plot(FC,THC_resid)
ccf(drop(FC),drop(THC_resid),120,main="CCF: Residuals of THC vs. FC")

win.graph()
acf2(THC_resid,120)

D_THC_resid=diff(THC_resid)
win.graph()
par(mfrow=c(2,1))
plot(D_THC_resid,main="First difference of THC residuals")
acf(D_THC_resid,120,main="First difference of THC residuals")
win.graph()
acf2(D_THC_resid,120)

D8_D_THC_resid=diff(D_THC_resid,8)
win.graph()
acf2(D8_D_THC_resid,120)



#Let's try some ARIMA models (not seasonal) with different parameters: 
out = matrix(ncol = 10)
d =1
P=0;D=0;Q=0;S=0;
for (p in 0:6) { for (q in 0:4) { 
fits = sarima(THC_resid, p, d, q)
#title(sub =  paste("p,q,P,D,Q =", p,q,P,D,Q), outer = T, line = -1)
out = rbind(out, c(p, d, q, P, D, Q, S, fits$AIC, fits$AICc, fits$BIC))
}}
out1=out
out1

#And also some seasonal ARIMA models with different parameters at S=8: 
d = 1
D=1
S = 8
for (p in 0:2) { for (q in 0:2) { for (P in 0:2)   {for (Q in 0:1){
fits = sarima(THC_resid, p, d, q, P, D, Q, S)
title(sub =  paste("p,q,P,D,Q =", p,q,P,D,Q), outer = T, line = -1)
out = rbind(out, c(p, d, q, P, D, Q, S, fits$AIC, fits$AICc, fits$BIC))
}}}}
out = out[-1,]
colnames(out) = c("p", "d", "q", "P", "D", "Q", "S", "AIC", "AICc","BIC")
out

#The model Chosen by BIC:
p=5;d=1;q=2;P=0;D=0;Q=0;
win.graph()
fits1 = sarima(THC_resid, p, d, q)
title(sub =  paste("Chosen by BIC: p,d,q =", p,d,q), outer = T, line = -1)

#The model Chosen by AIC and AICc:
win.graph()
p=5;d=1;q=3;P=0;D=0;Q=0;
fits2 = sarima(THC_resid, p, d, q)
title(sub =  paste("Chosen by AIC and AICc: p,d,q=", p,d,q), outer = T, line = -1)

#Due to similar results, We choose the first one which is simpler as the final model: ARIMA(5,1,2)

##############################
#How well is the model in forcasting?
h = 60
n1 = length(THC_resid) - h
X1 = ts(THC_resid[1:n1])
X2= ts(THC_resid[(n1+1):length(THC_resid)], start = n1+1)
win.graph()
X2.pr = sarima.for(X1, n.ahead = h, p=5,d=1,q=2)
lines(X2)
abline(v=n1+.5)

THC_resid.pr=ts(X2.pr$pred,frequency=60,start=(23+1/60))

#And the prediction based on regression on fundamental current
FC2=ts(FC[n1+1:length(FC)],start=n1+1,end=length(FC))
THC.RegresPred=0.06199*FC2+0.39704

#Therefore the total prediction will be:
THC.TotalPred=THC.RegresPred+X2.pr$pred
THC.predicted=ts(THC.TotalPred,frequency=60,start=(23+1/60))
#Computing upper and lower estimation band (95% confidence interval)
THC.Upper=THC.RegresPred+X2.pr$pred+2*X2.pr$se
THC.UpperLimit=ts(THC.Upper,frequency=60,start=(23+1/60))
THC.Lower=THC.RegresPred+X2.pr$pred-2*X2.pr$se
THC.LowerLimit=ts(THC.Lower,frequency=60,start=(23+1/60))

#Plotting the final result:
win.graph()
par(mfrow=c(2,1))
plot(THC,lwd=1.5,main="Triplen Harmonic Current (Predicted vs. Actual)",xlab="Time (Hour)")
lines(THC.predicted,col='red')
points(THC.predicted,col='red',pch='o',cex=.4)
lines(THC.UpperLimit,lty=2,col='blue',lwd=1.5)
lines(THC.LowerLimit,lty=2,col='blue',lwd=1.5)
abline(v=23,lty=6)
#Plotting a zoomed version:
plot(THC,xlim=c(20,24.5),ylim=c(13,20),lwd=2,main="Triplen Harmonic Current (Predicted vs. Actual)",xlab="Time (Hour)")
lines(THC.predicted,col='red')
points(THC.predicted,col='red',pch='o',cex=.7)
lines(THC.UpperLimit,lty=2,col='blue',lwd=1.5)
lines(THC.LowerLimit,lty=2,col='blue',lwd=1.5)
abline(v=23,lty=6)
