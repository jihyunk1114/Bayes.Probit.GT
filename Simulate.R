N = 5000                                  # Number of in-sample individuals
d = 5                                     # Pool size
Se.true = c(0.95,0.98)                    # Sensitivity
Sp.true = c(0.98,0.99)                    # Specificity
X<-cbind(rnorm(N,0,0.5),rbinom(N,1,0.5))  # Covariate
c<-runif(N,0.1,2)                         # test time

## Model ##
g<-function(t){log(t/12)+t^2/2}           # function
theta = c(-0.5,0.5)                       # regression coefficient

#########################################
# Generate true statuses
set.seed(12345)
p = pnorm(g(c)+X%*%theta)
Y.true = rbinom(N,1,p)
mean(Y.true)

library(Rcpp)
Rcpp::sourceCpp("SampLatent.cpp")
source("Testing Functions.txt")
source("Bayes.Probit.GT.R")

#########################################
# Simulate master pooling group testing
Pool.Data = Pool.test(Y.true,0.95,0.98,d)
Z1 = Pool.Data$Z
Y1 = Pool.Data$Y

mp = Bayes.Probit.GT(Z1,X,Y1,c,Se=0.95,Sp=0.98,na=1)
mp$summary.theta

########################################
# Simulate Dorfman group testing
Dorf.Data = Dorfman.decode.diff.error(Y.true, Se.true, Sp.true, d)
Z2 = Dorf.Data$Z
Y2 = Dorf.Data$Y

df = Bayes.Probit.GT(Z2,X,Y2,c,Se=Se.true,Sp=Sp.true,na=2)
df$summary.theta

# unknown assay accuracies
df.unknown = Bayes.Probit.GT(Z2,X,Y2,c,na=2,err.est=TRUE)
df.unknown$summary.theta

########################################
# Simulate Array testing
Array.Data = Array.decode.diff.error(Y.true, Se.true, Sp.true, d)
Z3 = Array.Data$Z
Y3 = Array.Data$Y

at = Bayes.Probit.GT(Z3,X,Y3,c,Se=Se.true,Sp=Sp.true,na=2)
at$summary.theta

# unknown assay accuracies
at.unknown = Bayes.Probit.GT(Z3,X,Y3,c,na=2,err.est=TRUE)
at.unknown$summary.theta

