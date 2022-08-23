Euler3D <- function(N =1000,M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt,driftx,diffx,drifty,diffy,
                    driftz,diffz,type=c("ito","str"),corr=NULL,...)
{
  if (type=="ito"){
    Ax    <- function(t,x,y,z) eval(driftx)
    Ay    <- function(t,x,y,z) eval(drifty)
    Az    <- function(t,x,y,z) eval(driftz)}else{
      if (is.null(corr)){
        driftstrx <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]]))))
        driftstry <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]]))))
        driftstrz <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2), list(e1 = driftz[[1]], e2 = diffz[[1]], de2 = Deriv(diffz,"z",cache.exp=FALSE)[[1]]))))
      }else{
        C <- stats::cov2cor(corr)
        driftstrx <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho1 *e3 * de3 + 0.5 * rho2 *e4 * de4), list(e1 = driftx[[1]], e2 = diffx[[1]], de2 = Deriv(diffx,"x",cache.exp=FALSE)[[1]],
                                                                                                                                  e3 = diffy[[1]], de3 = Deriv(diffx,"y",cache.exp=FALSE)[[1]],e4 = diffz[[1]], de4 = Deriv(diffx,"z",cache.exp=FALSE)[[1]],roh1=C[1,2],rho2=C[1,3]))))
        driftstry <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho1 *e3 * de3 + 0.5 * rho3 *e4 * de4), list(e1 = drifty[[1]], e2 = diffy[[1]], de2 = Deriv(diffy,"y",cache.exp=FALSE)[[1]],
                                                                                                                                  e3 = diffx[[1]], de3 = Deriv(diffy,"x",cache.exp=FALSE)[[1]],e4 = diffz[[1]], de4 = Deriv(diffy,"z",cache.exp=FALSE)[[1]],roh1=C[1,2],rho3=C[2,3]))))			  
        driftstrz <- eval(Simplify(substitute(expression(e1 + 0.5 * e2 * de2 + 0.5 * rho2 *e3 * de3 + 0.5 * rho3 *e4 * de4), list(e1 = driftz[[1]], e2 = diffz[[1]], de2 = Deriv(diffz,"z",cache.exp=FALSE)[[1]],
                                                                                                                                  e3 = diffx[[1]], de3 = Deriv(diffz,"x",cache.exp=FALSE)[[1]],e4 = diffy[[1]], de4 = Deriv(diffz,"y",cache.exp=FALSE)[[1]],roh2=C[1,3],rho3=C[2,3]))))			  
      }
      Ax <- function(t,x,y,z) eval(driftstrx) 
      Ay <- function(t,x,y,z) eval(driftstry)
      Az <- function(t,x,y,z) eval(driftstrz)
    }
  Sx <- function(t,x,y,z) eval(diffx)
  Sy <- function(t,x,y,z) eval(diffy) 
  Sz <- function(t,x,y,z) eval(diffz)
  if (missing(Dt)) {
    t <- seq(t0, T, by=Dt)
  } else {
    t <- c(t0, t0 + cumsum(rep(Dt, N)))
    T <- t[N + 1]
  } 
  Dt <- (T - t0)/N
  Z1 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
  Z2 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
  Z3 <- matrix(rnorm(N * M, 0, sqrt(Dt)), N, M)
  if (!is.null(corr)){
    C <- stats::cov2cor(corr)
    W1 <- Z1
    W2 <- C[1,2] * Z1 + sqrt(1-C[1,2]^2) * Z2	
    W3 <- C[1,3] * Z1 + ((C[2,3]-C[1,2]*C[1,3])/(sqrt(1-C[1,2]^2))) * Z2 + ( sqrt(1-C[1,3]^2 - ((C[2,3]-C[1,2]*C[1,3])/(sqrt(1-C[1,2]^2)))^2 )  ) * Z3
  }else{
    W1 <- Z1
    W2 <- Z2
    W3 <- Z3
  }	
  X <- matrix(x0, N+1, M)
  Y <- matrix(y0, N+1, M)	
  Z <- matrix(z0, N+1, M)	
  for (i in 1L:N) {
    X[i + 1L,] <- X[i,] + Ax(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sx(t[i],X[i,],Y[i,],Z[i,]) * W1[i,]  
    Y[i + 1L,] <- Y[i,] + Ay(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sy(t[i],X[i,],Y[i,],Z[i,]) * W2[i,]
    Z[i + 1L,] <- Z[i,] + Az(t[i],X[i,],Y[i,],Z[i,]) * Dt + Sz(t[i],X[i,],Y[i,],Z[i,]) * W3[i,]}
  name <- c("X","Y","Z")
  name <- if(M > 1) c(paste(name[1],1:M,sep=""),paste(name[2],1:M,sep=""),paste(name[3],1:M,sep=""))
  X <- ts(X, start = t0, deltat = Dt, names=name[1:M])
  Y <- ts(Y, start = t0, deltat = Dt, names=name[(M+1):(2*M)])
  Z <- ts(Z, start = t0, deltat = Dt, names=name[(2*M+1):(3*M)])
  return(list(X=X,Y=Y,Z=Z))
} 

fx <- expression(2*(y>0)-2*(z<=0) , -2*y, -2*z)
gx <- rep(expression(0.2),3)

Euler3D(N =1000, M=1,x0=0,y0=0,z0=0,t0=0,T=1,Dt,driftx = 0.2, diffx = expression(2*(y>0)-2*(z<=0)), 
        drifty = 0.2,diffy = expression(-2*y),
        driftz = 0.2, diffz = expression(-2*z), type=c("ito"),corr=NULL)
