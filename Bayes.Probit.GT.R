# Input:
#Z: A matrix of testing responses whose jth row is of the form (col1=Z_j, col2=dj, col3=assay used, col4:col(4+dj-1)=indices of the individuals assigned to the jth pool)
#X: Covariate matrix in which the ith row contains the covariate information pertaining to the ith individual
#Y: matrix whose ith row is of the form (col1=0, col2=number of pools ith individual involved in, col3:col(3+col2-1)=pools ith individual was assigned to)
#c: censoring/testing time for each individual
#grid: Define grid and find corresponding baseline survival function. Default is NULL.
#n.grid: length of grid.
#theta0, Sigma0, m0, v0, a0, b0, ae, be, ap, bp: For priors
#Se: vector of sensitivity values, if known
#Sp: vector of specificity values, if known
#order: order for I splines (usually 3 or 4)
#knots: NULL is default, but you can set up own knots
#m: number of interior knots
#quantile: If knots = NULL and quantile = TRUE, knots are created based on quantiles.
#          If knots = NULL and quantile = FALSE, equally spaced knots are created.
#maxiter: maximum iteration number
#burn.in: burn-in number
#na: number of arrays
#err.est: Default is FALSE. If assay accuracies are unknown change this to TRUE.
#' @import Rcpp
#' @import mvtnorm
#' @import truncnorm
#'
#' @export
Bayes.Probit.GT = function(Z, X, Y, c, grid = NULL, n.grid = 100,
                 init.theta = numeric(ncol(X)), eta = 1, gam0 = -3, gam = NULL,
                 theta0 = NULL, Sigma0 = NULL, m0 = 0, v0 = 0.1, a0 = 1, b0 = 1,
                 ae = 1, be = 1, ap = 1, bp = 1, Se = NULL, Sp = NULL, order = 3,
                 knots = NULL, m = 10, quantile = TRUE,
                 maxiter = 5000, burn.in = 2000, err.est = FALSE) {

  #make a note that m will be re-defined when knots are specified.
  if(is.null(knots)==FALSE) {print("m will be redefined when knots are specified")}
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

  #package
  Rcpp::sourceCpp("SampLatent.cpp")

  # when user do not specify knots or grid:
  if (is.null(knots)) {
    if (quantile == TRUE) {knots = c(0, quantile(c, seq(0.1, 0.9, length.out = m)), max(c) + 0.01)
    } else{knots = seq(0, max(c)+ 0.01, length.out = m + 2)}
  } else{m = length(knots) - 2}
  if (is.null(grid)) {grid = seq(0, max(c)+0.01, length.out = n.grid)}

  n.knots = length(knots); n.grid = length(grid)

  N = nrow(X); nx = ncol(X); na = max(Z[,3])

  # monotone splines
  I = Ispline(c, order, knots = knots); I1 = Ispline(grid, order, knots = knots)

  # latent varirable ui
  ui = numeric(N)

  #
  if(err.est == TRUE) {
    #print("Warning: Unknown assay accuracies are not applicable for master pooling.")
    Se.mat = matrix(-99,nrow=maxiter,ncol=na);Sp.mat = matrix(-99,nrow=maxiter,ncol=na)
    Se = rep(0.9,na); Sp = rep(0.9,na)
  }

  #prior
  if(is.null(theta0)) {theta0 = numeric(nx)}
  if(is.null(Sigma0)) {Sigma0 = nrow(X)*solve(t(X)%*%X)}

  k = m + order
  #init
  theta = init.theta; gam = rep(0.1,k)

  # To save chains
  theta.mat = matrix(0, maxiter, k + nx + 1)
  g.mat = matrix(0, maxiter, n.grid)

  for (iter in 1:maxiter) {
    m = t(gam0 + gam %*% I) + X %*% theta
    p = pnorm(t(gam0 + gam %*% I) + X %*% theta)

    u = runif(N)
    Y[,1] = SampLatent(N, p, Y=Y, Z=Z, U=u, se=Se,sp=Sp,na)  #na=number of assays

    ui[Y[,1] == 0] = truncnorm::rtruncnorm(sum(Y[,1]==0),-Inf,0,m[Y[,1]==0],1)
    if (any(Y[, 1] == 1)) {
        ui[Y[, 1] == 1] = truncnorm::rtruncnorm(sum(Y[, 1] == 1), 0, Inf, m[Y[, 1] == 1], 1)
    }

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
      if(na == 1) {
        flag = sapply(1:nrow(Z),function(x) Y[setdiff(Z[x,-(1:3)],-99),1])
        Z.true = unlist(lapply(flag,function(x){ifelse(sum(x)>0,1,0)}))
        Se = rbeta(1,ae+sum(Z[,1]*Z.true),be+sum((1-Z[,1])*Z.true))
        Sp = rbeta(1,ap+sum((1-Z[,1])*(1-Z.true)),bp+sum(Z[,1]*(1-Z.true)))
      } else{
        flag1 = sapply(1:nrow(Z),function(x) Y[setdiff(Z[x,-(1:3)],-99),1])
        Z.true1 = unlist(lapply(flag1,function(x){ifelse(sum(x)>0,1,0)}))
        Z.list = sapply(1:na, function(i) Z[Z[,3]==i,1])
        Z.true.list = sapply(1:na, function(i) Z.true1[Z[,3]==i])
        Se = sapply(1:na,function(i){rbeta(1,ae+sum(Z.list[[i]]*Z.true.list[[i]]),
                                   be+sum((1-Z.list[[i]])*Z.true.list[[i]]))})
        Sp = sapply(1:na,function(i){rbeta(1,ap+sum((1-Z.list[[i]])*(1-Z.true.list[[i]])),
                                   bp+sum(Z.list[[i]]*(1-Z.true.list[[i]])))})
      }
      Se.mat[iter,] = Se; Sp.mat[iter,] = Sp
    }

    #Save Gibbs result
    theta.mat[iter, 1:nx] = theta
    theta.mat[iter, nx + 1] = gam0
    theta.mat[iter,-(1:(nx + 1))] = gam
    g.mat[iter,] = t(gam0 + gam %*% I1)

    if (iter %% 100 == 0) {cat(paste("iteration", iter, "complete\n"))}

  }

  #Baseline survival function Plot
  St = colMeans(1 - pnorm(g.mat[-(1:burn.in),]))
  plot(grid, St, xlab = "t", ylab = expression(S[0](t)), lwd = 2, col = 2, type = "l")

  #Summary table
  theta.mat = theta.mat[-(1:burn.in),]
  mean = colMeans(theta.mat[, 1:nx])
  q.025 = apply(theta.mat[, 1:nx], 2, function(x) {quantile(x, 0.025)})
  q.975 = apply(theta.mat[, 1:nx], 2, function(x) {quantile(x, 0.975)})
  rname = NULL
  for (ii in 1:nx) {rname = c(rname, paste0("theta", ii))}
  summary.theta = cbind(mean, q.025, q.975)

  if(err.est == TRUE){
    if(na == 1) {
      Se.mat = Se.mat[-(1:burn.in)]; Sp.mat = Sp.mat[-(1:burn.in)]
      new = rbind(c(mean(Se.mat), quantile(Se.mat, c(0.025,0.975))),
      c(mean(Sp.mat), quantile(Sp.mat, c(0.025,0.975))))
      rname = c(rname,"Se1","Sp1")
      } else{
        Se.mat = Se.mat[-(1:burn.in),]; Sp.mat = Sp.mat[-(1:burn.in),]
        new = rbind(cbind(colMeans(Se.mat), apply(Se.mat, 2, function(x) {quantile(x, 0.025)}),
        apply(Se.mat, 2, function(x) {quantile(x, 0.975)})),
        cbind(colMeans(Sp.mat), apply(Sp.mat, 2, function(x) {quantile(x, 0.025)}),
        apply(Sp.mat, 2, function(x) {quantile(x, 0.975)})))
        for (ii in (nx+1):(nx+na)) {rname = c(rname, paste0("Se", ii-nx))}
        for (ii in (nx+na+1):(nx+2*na)) {rname = c(rname, paste0("Sp", ii-nx-na))}
      }
    
    theta.mat = cbind(theta.mat,Se.mat,Sp.mat)
    summary.theta = rbind(summary.theta,new)
  }

  summary.theta = as.data.frame(summary.theta, row.names = rname)

  #Results to show
  return(list(theta.mat = theta.mat, g.mat = g.mat, grid = grid, St = St, summary.theta = summary.theta))

  #summary (list: mean/ quantiles / baseline cdf /grid and S(t)) # plot grid(t) and S(t)

}

