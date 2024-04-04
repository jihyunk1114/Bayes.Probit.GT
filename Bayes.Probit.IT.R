#Y: individual status by group test result (length: sample size)
#X: covariates in matrix form
#c: test time
#grid: can define grid and find corresponding baseline survival function
#theta0, Sigma0, m0, v0, a0, b0: For priors
#Se: Sensitivity
#Sp: Specificity
#order: order for I splines (usually 3 or 4)
#knots: If NULL default, can set up own knots
#maxiter: maximum iteration number
#burn.in: burn-in number
#make a note that m will be redefined when knots are specified.
#warning if both Y and Yi not provided
#m: number of interior knots
#' @export
Bayes.Probit.IT = function(X, Y, c, grid = NULL, n.grid = 100, init.theta = numeric(ncol(X)),
          eta = 1, gam0 = -3, gam = NULL, theta0 = NULL, Sigma0 = NULL, m0 = 0, v0 = 0.1,
          a0 = 1, b0 = 1, ae = 1, be = 1, ap = 1, bp = 1, Se = NULL, Sp = NULL, order = 3,
          knots = NULL, m = 10, quantile = TRUE, maxiter = 5000, burn.in = 2000, err.est = FALSE) {

  #make a note that m will be re-defined when knots are specified.
  if(is.null(knots)==F) {print("m will be redefined when knots are specified")}
  if(err.est == FALSE & (is.null(Se) | is.null(Sp))) {print("Se and Sp should be specified if err.est = FALSE")}

  #Ispline function
  Ispline<-function(x,order,knots){
  # get I spline matrix with order
  # x is a row vector
  # k is the order of I spline
  # knots are a sequence of increasing points
  # the number of free parameters in M spline is the length of knots plus 1.


  m=length(knots)
  ### get Ispline bases ###
  k=order+1
  n=m-2+k # number of parameters
  t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

  yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
  for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
  }

  yytem1=yy1
  for (ii in 1:order){
    yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
    for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
    }
    yytem1=yytem2
  }


  index=rep(0,length(x))
  for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
  }

  ibases=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

  if (order==1){
    for (i in 2:n){
      ibases[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
    }
  }else{
    for (j in 1:length(x)){
      for (i in 2:n){
        if (i<(index[j]-order+1)){
          ibases[i-1,j]=1
        }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
          ibases[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
        }else{
          ibases[i-1,j]=0
        }
      }
    }
  }
  return(ibases)
}

  # when user do not specify knots or grid:
  if (is.null(knots)) {
    if (quantile == TRUE) {knots = c(0, quantile(c, seq(0.1, 0.9, length.out = m)), max(c) + 0.01)
    } else{knots = seq(0, max(c)+ 0.01, length.out = m + 2)}
  } else{m = length(knots) - 2}
  if (is.null(grid)) {grid = seq(0, max(c)+0.01, length.out = n.grid)}

  n.knots = length(knots); n.grid = length(grid)

  if(err.est == TRUE) {
    #print("Warning: Unknown assay accuracies are not applicable for master pooling.")
    Se.mat = rep(0,nrow=maxiter);Sp.mat = rep(0,nrow=maxiter)
    Se = 0.9; Sp = 0.9
  }
  
  N = nrow(X); p = ncol(X); lambda = Se + Sp - 1
  I = Ispline(c, order, knots = knots); I1 = Ispline(grid, order, knots = knots)
  Y.true = numeric(N); ui = Y.true

  #prior
  if(is.null(theta0)) {theta0 = numeric(ncol(X))}
  if(is.null(Sigma0)) {Sigma0 = nrow(X)*solve(t(X)%*%X)}

  k = m + order
  #init
  theta = init.theta; gam = rep(0.1,k)

  # To save chains
  theta.mat = matrix(0, maxiter, k + p + 1)
  g.mat = matrix(0, maxiter, n.grid)

  for (iter in 1:maxiter) {
    u = t(gam0 + gam %*% I) + X %*% theta
    
    #Calculate probability
    m1 = 1 - pnorm(u)
    prY1 = Se - lambda * m1; prD1 = 1 - m1
    prY0 = 1 - prY1; prD0 = 1 - prD1
    prob1 = Se * prD1 / prY1; prob2 = 1 - Sp * prD0 / prY0
    prob1[prob1 > 1] = 1; prob2[prob2 > 1] = 1
    
    #Derive true group status
    Y.true[Y == 1] = rbinom(sum(Y==1), 1, prob1[Y == 1])
    Y.true[Y == 0] = rbinom(sum(Y==0), 1, prob2[Y == 0])

    ui[Y.true == 0] = truncnorm::rtruncnorm(sum(Y.true==0),-Inf,0,u[Y.true==0],1)
    ui[Y.true == 1] = truncnorm::rtruncnorm(sum(Y.true==1),0, Inf,u[Y.true==1],1)
    
    #gam0 posterior
    W0 = v0 + N; E0 = (v0 * m0 + sum(unlist(ui) - t(I) %*% gam - X %*% theta)) / W0
    gam0 = rnorm(1, E0, sqrt(1 / W0))

    #gam posterior
    W = rowSums(I ^ 2)
    for (l in 1:length(W)) {
      if (W[l] == 0) {
        gam[l] = rexp(1, eta)
      } else{
        El = (sum(I[l,] * (ui - gam0 - t(I[-l,]) %*% gam[-l] - X %*% theta)) - eta) / W[l]
        gam[l] = truncnorm::rtruncnorm(1, 0, Inf, El, sqrt(1 / W[l]))
      }
    }

    #theta posterior
    Sigmatilde = solve(solve(Sigma0) + t(X) %*% X)
    thetatilde = c(Sigmatilde %*% (solve(Sigma0) %*% theta0 + t(X) %*% (ui - gam0 - t(gam %*% I))))
    theta = c(mvtnorm::rmvnorm(1, thetatilde, Sigmatilde))

    #eta posterior
    eta = rgamma(1, a0 + k, b0 + sum(gam))
    
    #Se and Sp predicting
    if(err.est == TRUE){
      Se = rbeta(1,ae+sum(Y*Y.true),be+sum((1-Y)*Y.true))
      Sp = rbeta(1,ap+sum((1-Y)*(1-Y.true)),bp+sum(Y*(1-Y.true)))
      Se.mat[iter] = Se; Sp.mat[iter] = Sp
      lambda = Se + Sp - 1
    }
    
    #Save Gibbs result
    theta.mat[iter, 1:p] = theta
    theta.mat[iter, p + 1] = gam0
    theta.mat[iter,-(1:(p + 1))] = gam

    g.mat[iter,] = t(gam0 + gam %*% I1)

    if (iter %% 100 == 0) {cat(paste("iteration", iter, "complete\n"))}

  }

  #Baseline survival function Plot
  St = colMeans(1 - pnorm(g.mat[-(1:burn.in),]))
  plot(grid, St, xlab = "t", ylab = expression(S[0](t)), lwd = 2, col = 2, type = "l")

  #Summary table
  theta.mat = theta.mat[-(1:burn.in),]
  mean = colMeans(theta.mat[, 1:p])
  q.025 = apply(theta.mat[, 1:p], 2, function(x) {quantile(x, 0.025)})
  q.975 = apply(theta.mat[, 1:p], 2, function(x) {quantile(x, 0.975)})
  rname = NULL
  for (ii in 1:p) {rname = c(rname, paste0("theta", ii))}
  summary.theta = cbind(mean, q.025, q.975)
  
  if(err.est == TRUE){
    Se.mat = Se.mat[-(1:burn.in)]; Sp.mat = Sp.mat[-(1:burn.in)]
    new = rbind(c(mean(Se.mat), quantile(Se.mat, c(0.025,0.975))),
                c(mean(Sp.mat), quantile(Sp.mat, c(0.025,0.975))))
    rname = c(rname,"Se1","Sp1")
    theta.mat = cbind(theta.mat,Se.mat,Sp.mat)

    summary.theta = rbind(summary.theta,new)
  }
  
  summary.theta = as.data.frame(summary.theta, row.names = rname)

  #Results to show
  return(list(theta.mat = theta.mat, g.mat = g.mat, grid = grid, St = St, summary.theta = summary.theta))

  #summary (list: mean/ quantiles / baseline cdf /grid and S(t)) # plot grid(t) and S(t)

}

