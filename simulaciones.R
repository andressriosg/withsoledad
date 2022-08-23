library(bbmle)
library(deSolve)

intervalo <- list(1:150, 150:177, 177:205, 205:270, 270:312, 312:330, 330:363, 363:385)
x0 <- list(c(7749210, 0, 1, 0), NULL, NULL, NULL, NULL, NULL, NULL, NULL)
ajuste <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
parameter <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
maximaver <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
theta <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
likelihood <- list()
beta <- seq(0, 0.999, by = 0.01)
solutions <- list()
for (i in 1:length(intervalo)) {
  solutions[[i]] <- matrix(NA, ncol = 5, nrow = length(intervalo[[i]]))
}

maximaver[[1]] <- function(lbeta, lnu, lgamma) {
  weeks <- intervalo[[1]]
  parms <- c(beta=plogis(lbeta), nu=plogis(lnu), gamma = plogis(lgamma))
  x0 <- x0[[1]]
  out <- ode(y=x0, weeks, seir, parms)
  SD1 <- sqrt(sum((estimacion$X1[weeks]-out[,4])^2)/length(weeks))
  - sum(dnorm(estimacion$X1[weeks], mean=out[,4], sd=SD1, log=TRUE))
}

ajuste[[1]] <- mle2(maximaver[[1]],
                    start=list(lbeta=qlogis(0.00001),
                               lnu = qlogis(.02),
                               lgamma=qlogis(.02)),  
                    method="Nelder-Mead",
                    control=list(maxit=1E5,trace=0),
                    trace=FALSE)

theta[[1]] <- as.numeric(c(plogis(coef(ajuste[[1]])[1:3])))
parameter[[1]] <- c(beta= theta[[1]][1], nu = theta[[1]][2], gamma = theta[[1]][3])

solutions[[1]] <- ode(y=x0[[1]], intervalo[[1]], seir, parameter[[1]])
colnames(solutions[[1]]) <- c("time","S", "E", "I","R")

for (j in 2:length(intervalo)) {
  
  x0[[j]] <- c(solutions[[j-1]][length(intervalo[[j-1]]), 2],
               solutions[[j-1]][length(intervalo[[j-1]]), 3],
               estimacion$X1[intervalo[[j-1]][length(intervalo[[j-1]])]],
               estimacion1$X1[intervalo[[j-1]][length(intervalo[[j-1]])]])
  
  maximaver[[j]] <- function(lbeta, lnu, lgamma) {
    weeks <- intervalo[[j]]
    parms <- c(beta=plogis(lbeta), nu=plogis(lnu), gamma = plogis(lgamma))
    x0 <- x0[[j]]
    out <- ode(y=x0, weeks, seir, parms)
    SD1 <- sqrt(sum((estimacion$X1[weeks]-out[,4])^2)/length(weeks))
    - sum(dnorm(estimacion$X1[weeks], mean=out[,4], sd=SD1, log=TRUE))
    + sqrt(sum((estimacion1$X1[weeks]-out[,5])^2)/length(weeks))
    - sum(dnorm(estimacion1$X1[weeks], mean=out[,5], sd=SD1, log=TRUE))
  }
  
  ajuste[[j]] <- mle2(maximaver[[j]],
                      start=list(lbeta=qlogis(0.01),
                                 lnu=qlogis(.002),
                                 lgamma=qlogis(.02)),  
                      method="Nelder-Mead",
                      control=list(maxit=1E5,trace=0),
                      trace=FALSE)
  
  theta[[j]] <- as.numeric(c(plogis(coef(ajuste[[j]])[1:3])))
  parameter[[j]] <- c(beta= theta[[j]][1], nu = theta[[j]][2], gamma = theta[[j]][3])
  
  solutions[[j]] <- ode(y=x0[[j]], intervalo[[j]], seir, parameter[[j]])
  colnames(solutions[[1]]) <- c("time","S", "E", "I","R")
}

plot(tiempo, estimacion$X1, col = "blue", lwd = 2, type = "l", ylab = "Number of cases",
     xlab = "Number of days")
for (k in 1:length(intervalo)) {
  lines(intervalo[[k]], solutions[[k]][, 4], col = k, type = "l", lwd = 2)
}

plot(1:315, rep(0, 315), col = 1, type = "l", lwd = 2,
     ylab = "Number of susceptibles", xlab = "Number of days",
     ylim = c(0, max(solutions[[1]])))
for (k in 1:length(intervalo)) {
  lines(intervalo[[k]], solutions[[k]][, 2], col = k, type = "l", lwd = 2)
}

plot(1:315, rep(0, 315), col = "white", type = "l", lwd = 2,
     ylab = "Number of exposed", xlab = "Number of days",
     ylim = c(0, max(solutions[[1]]))) # Dividir
for (k in 1:length(intervalo)) {
  lines(intervalo[[k]], solutions[[k]][, 3], col = k, type = "l", lwd = 2)
}