######################## Final Project #######################

##################### Reading the Data #######################

lead <- read.table("C:/Users/Dell/Dropbox/NCSU/ST 537/Project/lead.full.txt", header = F)
colnames(lead) = c("id", "ind.age", "sex", "week", "blood", "trt")
head(lead)

Y= lead$blood
week = lead$week
ind.age = lead$ind.age
sex= lead$sex


### Indicator variable for groups
C1= as.numeric(lead$trt==1)
C2= as.numeric(lead$trt==2)
C3= as.numeric(lead$trt==3)
lead=cbind(lead,C1,C2,C3)


### Treatment Groups
lead1 = lead[C1 == 1,1:7]
lead2 = lead[C2 == 1,1:7]
lead3 = lead[C3 == 1,1:7]

## Number of rows for each treatment
table(lead$trt)

### Models 
meanform= blood  ~ -1 + C1 + C1:week + C1:sex + C1:ind.age + C1:week:ind.age+ C1:sex:week+ 
   C1:sex:ind.age+ C1:sex:week:ind.age+
  (C2+C2:week + C2:sex + C2:ind.age + C2:week:ind.age+ C2:sex:week+
   C2:sex:ind.age+ C2:sex:week:ind.age)+
  (C3+C3:week + C3:sex + C3:ind.age + C3:week:ind.age+ C3:sex:week+
   C3:sex:ind.age+ C3:sex:week:ind.age)


##Independent, where error variance does not change over weeks

library(nlme)
fit.in.ev = lme(fixed= meanform, 
           data=lead, random = ~ week|id, method="ML",
           control = lmeControl(opt='optim'))
summary(fit.in.ev)


### Independent, where error variance changes over weeks

fit.in.uv = lme(fixed= meanform, 
           data=lead, random = ~ week|id, method="ML",
           weights = varIdent(form = ~ 1 | week),
           control = lmeControl(opt='optim'))
summary(fit.in.uv )


### AR(1) correlation structure, where error variance does not change over weeks

fit.ar1.ev = lme(fixed= meanform, 
           data=lead, random = ~ week|id, method="ML",
           correlation = corAR1(form = ~ week | id),
           control = lmeControl(opt='optim'))
summary(fit.ar1.ev)


### AR(1) correlation structure, where error variance changes over weeks

fit.ar1.uv = lme(fixed= meanform, 
           data=lead, random = ~ week|id, method="ML",
           weights = varIdent(form = ~ 1 | week),
           correlation = corAR1(form = ~ week | id),
           control = lmeControl(opt='optim'))
summary(fit.ar1.uv)


### Unstructured, where error variance does not change over weeks

lead$timefact= as.numeric( factor(lead$week, labels = 1:5) )

fit.un.ev = lme(fixed= meanform, 
           data=lead, random = ~ week|id, method="ML",
           correlation = corSymm(form= ~ timefact | id ),
           control = lmeControl(opt='optim'))
summary(fit.un.ev)


### Unstructured, where error variance changes over weeks

fit.un.uv =lme(fixed= meanform, 
           data=lead, random = ~ week|id, method="ML",
           weights = varIdent(form = ~ 1 | timefact),
           correlation = corSymm(form= ~ timefact | id ),
           control = lmeControl(opt='optim'))
summary(fit.un.uv)

############## aic and bic ###############

aic = AIC(fit.in.ev,fit.in.uv,fit.ar1.ev,fit.ar1.uv,fit.un.ev, fit.un.uv )
bic = BIC(fit.in.ev,fit.in.uv,fit.ar1.ev,fit.ar1.uv,fit.un.ev, fit.un.uv )
# Display
cbind(aic, bic$BIC)

############## loglikelihood values ###############

loglik.1 = logLik(fit.in.ev)
loglik.2 = logLik(fit.in.uv)
loglik.3 = logLik(fit.ar1.ev)
loglik.4 = logLik(fit.ar1.uv)
loglik.5 = logLik(fit.un.ev)
loglik.6 = logLik(fit.un.uv)

cbind(loglik.1,loglik.2,loglik.3,loglik.4,loglik.5,loglik.6)

################# Robust covariance matrix and se of betahat ########################

library(clubSandwich)
betahat = fixed.effects(fit.in.ev)

# Robust covariance matrix and se of betahat
V.robust = vcovCR(fit.in.ev, type = "CR0")
se.robust = sqrt(diag(V.robust))

# Standard error
##SE = sqrt( diag(fit.a$varFix) )

# CI limits
df = nrow(lead) - length(betahat)
t.alpha = qt(0.05/2, df = df, lower.tail = F)
lower = betahat - t.alpha*se.robust
upper = betahat + t.alpha*se.robust

# Display the estimates
tab = cbind(betahat, se.robust, lower, upper)
tab

################# Effect of age on blood lead level ###################### 

# L matrix 
La =  matrix(0, 12, 24)
La[1,6]=1
La[2,9]=1
La[3,12]=1
La[4,13]=1
La[5,15]=1
La[6,16]=1
La[7,18]=1
La[8,19]=1
La[9,21]=1
La[10,22]=1
La[11,23]=1
La[12,24]=1
# Estimate and SE
estimate = La %*% betahat
SE = La %*% V.robust %*% t(La)
# confidence limits
df = nrow(lead) - length(betahat)
t.alpha = qt(0.05/2, df = df, lower.tail = F)
lower = estimate - t.alpha*SE
upper = estimate + t.alpha*SE
# results
tab = data.frame(estimate, SE, lower, upper)
round(tab, 4)


################### Hypothesis testing ####################### 

########## wald test #########
cc = nrow(La)
df = nrow(lead) - length(betahat)
# estimate and covariance matrix of L\beta
est = La %*% betahat
varmat = La %*% V.robust %*% t(La)
# Wald test
Wald = c( t(est) %*% solve(varmat) %*% (est) )
p.value = pchisq(q = Wald, df = cc, lower.tail=FALSE)
data.frame(Wald, p.value)


########## t stat ###############

# df
df = nrow(lead) - length(betahat)
# t stats
t.robust = betahat/se.robust
# p-values
p.value = round( 2*pt(q = abs(t.robust), df = df, lower.tail = FALSE), 4 )
# results
tab = data.frame(betahat, se.robust, t.robust, p.value)
tab


################# Effect of Gender on blood lead level ###################### 


# L matrix
Lg=  matrix(0, 12, 24)
Lg[1,5]=1
Lg[2,8]=1
Lg[3,11]=1
Lg[4,14]=1
Lg[5,15]=1
Lg[6,17]=1
Lg[7,18]=1
Lg[8,20]=1
Lg[9,21]=1
Lg[10,22]=1
Lg[11,23]=1
Lg[12,24]=1
Lg
# Estimate and SE
estimate.g = Lg %*% betahat
SE.g = Lg %*% V.robust %*% t(Lg)
# confidence limits
df = nrow(lead) - length(betahat)
t.alpha = qt(0.05/2, df = df, lower.tail = F)
lower.g = estimate.g - t.alpha*SE.g
upper.g = estimate.g + t.alpha*SE.g
# results
tab = data.frame(estimate.g, SE.g, lower.g, upper.g)
round(tab, 4)


################### Hypothesis testing ####################### 

########## wald test #########
cc = nrow(Lg)
df = nrow(lead) - length(betahat)
# estimate and covariance matrix of L\beta
est.g = Lg %*% betahat
varmat.g = Lg %*% V.robust %*% t(Lg)
# Wald test
Wald.g = c( t(est.g) %*% solve(varmat.g) %*% (est.g) )
p.value.g = pchisq(q = Wald.g, df = cc, lower.tail=FALSE)
data.frame(Wald.g, p.value.g)


########## t stat ###############

# df
df = nrow(lead) - length(betahat)
# t stats
t.robust = betahat/se.robust
# p-values
p.value = round( 2*pt(q = abs(t.robust), df = df, lower.tail = FALSE), 4 )
# results
tab = data.frame(betahat, se.robust, t.robust, p.value)
tab



################# Using L matrix for sex and age ###################### 

# L matrix
Li=  matrix(0, 6, 24)
Li[1,15]=1
Li[2,18]=1
Li[3,21]=1
Li[4,22]=1
Li[5,23]=1
Li[6,24]=1
Li
# Estimate and SE
estimate.s = Li %*% betahat
SE.s = Li %*% V.robust %*% t(Li)
# confidence limits
df = nrow(lead) - length(betahat)
t.alpha = qt(0.05/2, df = df, lower.tail = F)
lower.s = estimate.s - t.alpha*SE.s
upper.s = estimate.s + t.alpha*SE.s
# results
tab = data.frame(estimate.s, SE.s, lower.s, upper.s)
round(tab, 4)


################### Hypothesis testing ####################### 

########## wald test #########
cc = nrow(Li)
df = nrow(lead) - length(betahat)
# estimate and covariance matrix of L\beta
est.s = Li %*% betahat
varmat.s = Li %*% V.robust %*% t(Li)
# Wald test
Wald.s = c( t(est.s) %*% solve(varmat.s) %*% (est.s) )
p.value.s = pchisq(q = Wald.s, df = cc, lower.tail=FALSE)
data.frame(Wald.s, p.value.s)


################# Using L matrix for week and age ###################### 

# L matrix
Li=  matrix(0, 3, 24)
Li[1,13]=1
Li[2,16]=1
Li[3,19]=1
Li
# Estimate and SE
estimate.s = Li %*% betahat
SE.s = Li %*% V.robust %*% t(Li)
# confidence limits
df = nrow(lead) - length(betahat)
t.alpha = qt(0.05/2, df = df, lower.tail = F)
lower.s = estimate.s - t.alpha*SE.s
upper.s = estimate.s + t.alpha*SE.s
# results
tab = data.frame(estimate.s, SE.s, lower.s, upper.s)
round(tab, 4)


################### Hypothesis testing ####################### 

########## wald test #########
cc = nrow(Li)
df = nrow(lead) - length(betahat)
# estimate and covariance matrix of L\beta
est.s = Li %*% betahat
varmat.s = Li %*% V.robust %*% t(Li)
# Wald test
Wald.s = c( t(est.s) %*% solve(varmat.s) %*% (est.s) )
p.value.s = pchisq(q = Wald.s, df = cc, lower.tail=FALSE)
data.frame(Wald.s, p.value.s)


############################## Smaller model ###############################

### Models 
meanform.s= blood ~ -1 + C1 + C1:week  +  C1:ind.age+
  (C2+C2:week  +  C2:ind.age)+
  (C3+C3:week +   C3:ind.age)


##Independent, where error variance does not change over weeks

library(nlme)
fit.in.ev.s = lme(fixed= meanform.s, 
                data=lead, random = ~ week|id, method="ML",
                control = lmeControl(opt='optim'))
summary(fit.in.ev.s)


###################### Significance of Smaller model ########################
# Full model
fit.in.ev = lme(fixed= meanform, 
                data=lead, random = ~ week|id, method="ML",
                control = lmeControl(opt='optim'))
p1 = 28
# Reduced model
fit.in.ev.s = lme(fixed= meanform.s, 
                  data=lead, random = ~ week|id, method="ML",
                  control = lmeControl(opt='optim'))
p2 = 13
# log-likelihoods
loglik.full = logLik(fit.in.ev)
loglik.red = logLik(fit.in.ev.s)
# LRT and p-value
df = p1 - p2
LRT = 2*(loglik.full - loglik.red)
p.value = pchisq(LRT, df = df, lower.tail = FALSE)
# results
data.frame(L.full = loglik.full, L.reduced = loglik.red,
           LRT = LRT, df = df, p.value = p.value)


############################ Using anova.lme()########################

anova.lme(fit.in.ev, fit.in.ev.s)


################# Robust covariance matrix and se of betahat ########################

library(clubSandwich)
betahat.s = fixed.effects(fit.in.ev.s)

# Robust covariance matrix and se of betahat
V.robust.s = vcovCR(fit.in.ev.s, type = "CR0")
se.robust.s = sqrt(diag(V.robust.s))

# Standard error
##SE = sqrt( diag(fit.a$varFix) )

# CI limits
df = nrow(lead) - length(betahat.s)
t.alpha = qt(0.05/2, df = df, lower.tail = F)
lower = betahat.s - t.alpha*se.robust.s
upper = betahat.s + t.alpha*se.robust.s

# Display the estimates
tab = cbind(betahat.s, se.robust.s, lower, upper)
tab

###################### Comparison of mean trends #######################

### Placebo and Succimer- Low
anova.lme(fit.in.ev.s, L = c(1,-1,0,0,0,0,0,0,0), adjustSigma = TRUE)
anova.lme(fit.in.ev.s, L = c(0,0,0,1,0,-1,0,0,0), adjustSigma = TRUE)
anova.lme(fit.in.ev.s, L = c(0, 0,0,0,1,0,-1,0,0), adjustSigma = TRUE)

### Placebo and Succimer- High
anova.lme(fit.in.ev.s, L = c(1, 0,-1,0,0,0,0,0,0), adjustSigma = TRUE)
anova.lme(fit.in.ev.s, L = c(0,0,0,1,0,0,0,-1,0), adjustSigma = TRUE)
anova.lme(fit.in.ev.s, L = c(0, 0,0,0,1,0,0,0,-1), adjustSigma = TRUE)

### Succimer- Low and Succimer -High
anova.lme(fit.in.ev.s, L = c(0, 1,-1,0,0,0,0,0,0), adjustSigma = TRUE)
anova.lme(fit.in.ev.s, L = c(0,0,0,0,0,1,0,-1,0), adjustSigma = TRUE)
anova.lme(fit.in.ev.s, L = c(0, 0,0,0,0,0,1,0,-1), adjustSigma = TRUE)


############################# Mean Trends ##############################

cf= fit.in.ev.s$coefficients$fixed
response <- function(c1,c2,c3,t,ag){
  
    u= c1*cf[1]+c1*cf[4]*t+c1*cf[5]*ag+
      c2*cf[2]+c2*cf[6]*t+c2*cf[7]*ag+
      c3*cf[3]+c3*cf[8]*t+c3*cf[9]*ag
 u
}   ### function to get estimated value


response <- Vectorize(response)
wk <- c(0,2,4,6,8)

P.Male.A0 = response(1,0,0,wk,0)
P.Male.A1 = response(1,0,0,wk,1)
SL.Male.A0 = response(0,1,0,wk,0)
SL.Male.A1 = response(0,1,0,wk,1)
SH.Male.A0 = response(0,0,1,wk,0)
SH.Male.A1 = response(0,0,1,wk,1)


matplot(c(0,2,4,6,8), P.Male.A0, type = "b", lty = 2, col = "green", pch="0",
        xlab = "week", ylab = "blood lead level", main = 'Estimated mean trends-Males',
        ylim = c(10,42))
lines(wk, P.Male.A1 , type = "b", lty = 2, col = "green", pch="1",)
lines(wk, SL.Male.A0, type = "b", lty = 2, col = "red", pch="0",)
lines(wk, SL.Male.A1, type = "b", lty = 2, col = "red", pch="1",)
lines(wk, SH.Male.A0, type = "b", lty = 2, col = "blue", pch="0",)
lines(wk, SH.Male.A1, type = "b", lty = 2, col = "blue", pch="1",)

legend(5,42,legend=c('Placebo','Succimer-Low','Succimer-High'),
       col=c("green","red", "blue"), lty=c(1,1))



############################# Model Diagnostics ##################################


################### Checking equal variance assumption ################

############ Subject level predictions - Fitted values ################

sub.pred = predict(fit.in.ev.s, level = 1)
P.pred = sub.pred[1:174]
SL.pred = sub.pred[175:349]
SH.pred = sub.pred[350:526]

############# Pearson Residuals ##############

res = resid(fit.in.ev.s, level = 1, type = "pearson")
P.res = res[1:174]
SL.res = res[175:349]
SH.res = res[350:526]

par(mfrow=c(1,3))

matplot(P.pred, P.res, col = "green", pch="P",
        xlab = "fitted values", ylab = "Standardised residuals", main = 'Placebo')
matplot(SL.pred, SL.res, col = "red", pch="L",
        xlab = "fitted values", ylab = "Standardised residuals", main = 'Succimer-Low')
matplot(SH.pred, SH.res, col = "blue", pch="H",
        xlab = "fitted values", ylab = "Standardised residuals", main = 'Succimer-High')


################### Checking assumption of normal errors ################

par(mfrow=c(1,3))

qqnorm(P.res, main = "Normal Q-Q Plot-Placebo",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(P.res)

qqnorm(SL.res, main = "Normal Q-Q Plot-Succimer Low",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(SL.res)

qqnorm(SH.res, main = "Normal Q-Q Plot-Succimer High",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(SH.res)


################ Quality of fit ###################


matplot(sub.pred, lead$blood, col = "orange", pch=16,
        xlab = "Fitted values", ylab = "Observed Response", main = 'Quality of fit')
abline(lm(lead$blood~ sub.pred))


################ Prediction of random effects ###################

b.hat = random.effects(fit.in.ev.s,)
b.hat[1:3, ]

#################### Q-Q plots for random effects ################

#################### Q-Q plots for intercepts ################

par(mfrow=c(1,3))
qqnorm(b.hat[1:40,1], main = "Placebo- Intercept",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(b.hat[1:40,1])

qqnorm(b.hat[41:80,1], main = "Succimer Low- Intercept",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(b.hat[41:80,1])

qqnorm(b.hat[81:120,1], main = "Succimer High - Intercept",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(b.hat[81:120,1])

#################### Q-Q plots for weeks ########################

par(mfrow=c(1,3))

qqnorm(b.hat[1:40,2], main = "Placebo- Week",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(b.hat[1:40,2])

qqnorm(b.hat[41:80,2], main = "Succimer Low- Week",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(b.hat[41:80,2])

qqnorm(b.hat[81:120,2], main = "Succimer High - Week",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(b.hat[81:120,2])

############### Scatter plots for random effects ######################

par(mfrow=c(1,3))

plot(b.hat[1:40,1],b.hat[1:40,2],main = "Placebo",
     xlab = "intercept", ylab = "week",col= "green",pch = 19,
     plot.it = TRUE, datax = FALSE)

plot(b.hat[41:80,1],b.hat[41:80,2],main = "Succimer-Low",
     xlab = "intercept", ylab = "week",col= "red",pch=19,
     plot.it = TRUE, datax = FALSE)

plot(b.hat[81:120,1],b.hat[81:120,2],main = "Succimer-High",
     xlab = "intercept", ylab = "week",col= "blue",pch=19,
     plot.it = TRUE, datax = FALSE)




#################### New Smaller model ###############################


### Models 
meanform.s2= blood ~   C1:week  + ind.age+ C2:week  +C3:week 


##Independent, where error variance does not change over weeks

library(nlme)
fit.in.ev.s2 = lme(fixed= meanform.s2, 
                  data=lead, random = ~ week|id, method="ML",
                  control = lmeControl(opt='optim'))
summary(fit.in.ev.s2)


###################### Significance of Smaller model ####################
# Full model
fit.in.ev = lme(fixed= meanform, 
                data=lead, random = ~ week|id, method="ML",
                control = lmeControl(opt='optim'))

p1 <- 28
# Reduced model
fit.in.ev.s2 = lme(fixed= meanform.s2, 
                  data=lead, random = ~ week|id, method="ML",
                  control = lmeControl(opt='optim'))
p2 <- 9
# log-likelihoods
loglik.full <- logLik(fit.in.ev)
loglik.red <- logLik(fit.in.ev.s2)
# LRT and p-value
df <- p1 - p2
LRT <- 2*(loglik.full - loglik.red)
p.value <- pchisq(LRT, df = df, lower.tail = FALSE)
# results
data.frame(L.full = loglik.full, L.reduced = loglik.red,
           LRT = LRT, df = df, p.value = p.value)


######### Using anova.lme()
anova.lme(fit.in.ev, fit.in.ev.s2)


############################# Mean Trends ##############################

cf= fit.in.ev.s2$coefficients$fixed
response <- function(c1,c2,c3,t,ag){
  
  u= cf[1]+ag*cf[2]+c1*cf[3]*t+c2*cf[4]*t+c3*cf[5]*t
  u
}   ### function to get estimated value


response <- Vectorize(response)
wk <- c(0,2,4,6,8)

P.Male.A0 = response(1,0,0,wk,0)
P.Male.A1 = response(1,0,0,wk,1)
SL.Male.A0 = response(0,1,0,wk,0)
SL.Male.A1 = response(0,1,0,wk,1)
SH.Male.A0 = response(0,0,1,wk,0)
SH.Male.A1 = response(0,0,1,wk,1)


matplot(c(0,2,4,6,8), P.Male.A0, type = "b", lty = 2, col = "green", pch="0",
        xlab = "week", ylab = "blood lead level", main = 'Estimated mean trends',
        ylim = c(10,42))
lines(wk, P.Male.A1 , type = "b", lty = 2, col = "green", pch="1",)
lines(wk, SL.Male.A0, type = "b", lty = 2, col = "red", pch="0",)
lines(wk, SL.Male.A1, type = "b", lty = 2, col = "red", pch="1",)
lines(wk, SH.Male.A0, type = "b", lty = 2, col = "blue", pch="0",)
lines(wk, SH.Male.A1, type = "b", lty = 2, col = "blue", pch="1",)

legend(5,42,legend=c('Placebo','Succimer-Low','Succimer-High'),
       col=c("green","red", "blue"), lty=c(1,1))


