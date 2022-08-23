install.packages("bbmle")
library(bbmle)
library(deSolve)

mu = 4/1000 
eta = (9.42/1000)*mean(proyeccion_bogota_final)
upsilon = 1/((5+6)/2)

intervalo <- list(1:150, 150:177, 177:205, 205:270, 270:312, 312:330, 330:363, 363:385)

x0 <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)

susceptibles_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
expuestos_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
infectados_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
recuperados_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
susceptibles1_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
expuestos1_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
infectados1_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
recuperados1_b <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)

for (j in 1:length(recuperados_b)) {
  susceptibles_b[[j]] <- susceptibles_bogota[intervalo[[j]]] 
  expuestos_b[[j]] <- expuestos_bogota[intervalo[[j]]] 
  infectados_b[[j]] <- suavizado_infectados$X1[intervalo[[j]]] 
  recuperados_b[[j]] <- suavizado_recuperados$X1[intervalo[[j]]] 
}

for (j in 1:length(recuperados_b)) {
  susceptibles1_b[[j]] <- susceptibles_bogota[intervalo[[j]]+1] 
  expuestos1_b[[j]] <- expuestos_bogota[intervalo[[j]]+1] 
  infectados1_b[[j]] <- suavizado_infectados$X1[intervalo[[j]]+1] 
  recuperados1_b[[j]] <- suavizado_recuperados$X1[intervalo[[j]]+1] 
}

parametros_a <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
for (j in 1:length(recuperados_b)) {
  parametros_a[[j]] <- c(eta, 
    ((1-mu)*sum(((susceptibles_b[[j]])^2*infectados_b[[j]])) - 
                sum(susceptibles1_b[[j]]*susceptibles_b[[j]]*infectados_b[[j]]) + 
                eta*sum(susceptibles_b[[j]]*infectados_b[[j]]) - 
                (1 - upsilon - mu)*sum(susceptibles_b[[j]]*expuestos_b[[j]]*infectados_b[[j]]) + 
                sum(expuestos1_b[[j]]*susceptibles_b[[j]]*
                      infectados_b[[j]]))/(2*sum((susceptibles_b[[j]])^2*(infectados_b[[j]])^2)), 
    upsilon, 
    ((1-mu)*(sum((infectados_b[[j]])^2) - sum(recuperados_b[[j]]*infectados_b[[j]])) - 
       sum(infectados1_b[[j]]*infectados_b[[j]]) + 
       sum(infectados1_b[[j]]*recuperados_b[[j]]) + 
       upsilon*sum(expuestos_b[[j]]*infectados_b[[j]]))/(2*sum((infectados_b[[j]])^2)), 
    mu)
}

parametros_a

solutions <- list()

seir <- function(t,x,parameter){
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  with(as.list(parameter),
       {
         dS <- parameter[1] - parameter[2]*S*I - parameter[5]*S # https://datosmacro.expansion.com/
         dE <- parameter[2]*S*I - parameter[3]*E - parameter[5]*E #
         dI <- parameter[3]*E - parameter[4]*I - parameter[5]*I
         dR <- parameter[4]*I - parameter[5]*R
         res <- c(dS, dE, dI, dR)
         list(res)
       })
}

for (i in 1:length(intervalo)) {
  library(deSolve)
  x0[[i]] <- c(susceptibles_b[[i]][1], expuestos_b[[i]][1], infectados_b[[i]][1], recuperados_b[[i]][1])
  solutions[[i]] <- ode(y = x0[[i]], intervalo[[i]], seir, parametros[[i]])
}

soluciones <- rbind(solutions[[1]], solutions[[2]], solutions[[3]], solutions[[4]], solutions[[5]], 
                    solutions[[6]], solutions[[7]], solutions[[8]])

plot(soluciones[,4], type = "p", col = "red", lwd = 2, pch = 16)
lines(expuestos_bogota, col = "blue", pch = 15, lwd = 3)

parametros # Los valores de gamma son demasiado altos 

parametermv <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
theta <- list()
solucionesmv <- list()
seir1 <- list()

for (j in 1:length(intervalo)) { 
  for (j in 1:length(intervalo)) {
    seir1[[j]] <- function(t,x,parameter){
      S <- x[1]
      E <- x[2]
      I <- x[3]
      R <- x[4]
      with(as.list(parameter),
           {
             dS <- parametros[[j]][1] - parametros[[j]][2]*S*I - parametros[[j]][5]*S
             dE <- parametros[[j]][2]*S*I - parametros[[j]][3]*E - parametros[[j]][5]*E 
             dI <- parametros[[j]][3]*E - parameter*I - parametros[[j]][5]*I
             dR <- parameter*I - parametros[[j]][5]*R
             res <- c(dS, dE, dI, dR)
             list(res)
           })
    }
  }
  
  maximaver[[j]] <- function(lgamma) {
    parameter <- c(gamma = plogis(lgamma))
    out <- ode(y= x0[[j]], intervalo[[j]], seir1[[j]], parameter)
    SD <- sqrt(sum((infectados_b[[j]][intervalo[[j]]] -out[,4])^2)/length(intervalo[[j]])) 
    - sum(dnorm(infectados_b[[j]][intervalo[[j]]], mean=out[,4], sd=SD, log=TRUE)) 
    + sqrt(sum((recuperados_b[[j]][intervalo[[j]]] -out[,5])^2)/length(intervalo[[j]]))
    -  sum(dnorm(recuperados_b[[j]][intervalo[[j]]], mean=out[,5], sd=SD, log=TRUE))
  }
  
  ajuste[[j]] <- mle2(maximaver[[j]],
                      start=list(lgamma= 1e-5),  
                      method="Brent",
                      control=list(maxit=1E5,trace=0),
                      trace=FALSE, lower = 0, upper = 1)
  
  theta[[j]] <- as.numeric(c(coef(ajuste[[j]])[1]))
  parametermv[[j]] <- c(gamma = theta[[j]][1])
  
  solucionesmv[[j]] <- ode(y = x0[[j]], intervalo[[j]], seir1[[j]], as.numeric(parametermv[[j]]))
}

soluciones_mv <- NULL

for (k in 1:8) {
  soluciones_mv <- rbind(solucionesmv[[k]])
}
plot(soluciones_mv[,3], type = "p", col = "red", lwd = 2, pch = 16)
lines(suavizado_infectados$X1[intervalo[[1]]], col = "blue", pch = 15, lwd = 3)


# Calculo del parámetro gamma por mínimos cuadrados ordinarios 

seir1 <- list()
for (j in 1:length(intervalo)) {
  seir1[[j]] <- function(t,x,parameter){
    S <- x[1]
    E <- x[2]
    I <- x[3]
    R <- x[4]
    with(as.list(parameter),
         {
           dS <- parametros[[j]][1] - parametros[[j]][2]*S*I - parametros[[j]][5]*S
           dE <- parametros[[j]][2]*S*I - parametros[[j]][3]*E - parametros[[j]][5]*E 
           dI <- parametros[[j]][3]*E - parameter*I - parametros[[j]][5]*I
           dR <- parameter*I - parametros[[j]][5]*R
           res <- c(dS, dE, dI, dR)
           list(res)
         })
  }
}

library(deSolve)
minimun <- list()

minimun[[1]] <- function(data, parameter){
    return(sum((data[, 1] - ode(y=x0[[1]], intervalo[[1]], seir1[[1]], parameter)[, 2])^2) + 
             sum((data[, 2] - ode(y=x0[[1]], intervalo[[1]], seir1[[1]], parameter)[, 3])^2) + 
             sum((data[, 3] - ode(y=x0[[1]], intervalo[[1]], seir1[[1]], parameter)[, 4])^2) + 
             sum((data[, 4] - ode(y=x0[[1]], intervalo[[1]], seir1[[1]], parameter)[, 5])^2))
}
minimun[[2]] <- function(data, parameter){
  return(sum((data[, 1] - ode(y=x0[[2]], intervalo[[2]], seir1[[2]], parameter)[, 2])^2) + 
           sum((data[, 2] - ode(y=x0[[2]], intervalo[[2]], seir1[[2]], parameter)[, 3])^2) + 
           sum((data[, 3] - ode(y=x0[[2]], intervalo[[2]], seir1[[2]], parameter)[, 4])^2) + 
           sum((data[, 4] - ode(y=x0[[2]], intervalo[[2]], seir1[[2]], parameter)[, 5])^2))
}
minimun[[3]] <- function(data, parameter){
  return(sum((data[, 1] - ode(y=x0[[3]], intervalo[[3]], seir1[[3]], parameter)[, 2])^2) + 
           sum((data[, 2] - ode(y=x0[[3]], intervalo[[3]], seir1[[3]], parameter)[, 3])^2) + 
           sum((data[, 3] - ode(y=x0[[3]], intervalo[[3]], seir1[[3]], parameter)[, 4])^2) + 
           sum((data[, 4] - ode(y=x0[[3]], intervalo[[3]], seir1[[3]], parameter)[, 5])^2))
}
minimun[[4]] <- function(data, parameter){
  return(sum((data[, 1] - ode(y=x0[[4]], intervalo[[4]], seir1[[4]], parameter)[, 2])^2) + 
           sum((data[, 2] - ode(y=x0[[4]], intervalo[[4]], seir1[[4]], parameter)[, 3])^2) + 
           sum((data[, 3] - ode(y=x0[[4]], intervalo[[4]], seir1[[4]], parameter)[, 4])^2) + 
           sum((data[, 4] - ode(y=x0[[4]], intervalo[[4]], seir1[[4]], parameter)[, 5])^2))
}
minimun[[5]] <- function(data, parameter){
  return(sum((data[, 1] - ode(y=x0[[5]], intervalo[[5]], seir1[[5]], parameter)[, 2])^2) + 
           sum((data[, 2] - ode(y=x0[[5]], intervalo[[5]], seir1[[5]], parameter)[, 3])^2) + 
           sum((data[, 3] - ode(y=x0[[5]], intervalo[[5]], seir1[[5]], parameter)[, 4])^2) + 
           sum((data[, 4] - ode(y=x0[[5]], intervalo[[5]], seir1[[5]], parameter)[, 5])^2))
}
minimun[[6]] <- function(data, parameter){
  return(sum((data[, 1] - ode(y=x0[[6]], intervalo[[6]], seir1[[6]], parameter)[, 2])^2) + 
           sum((data[, 2] - ode(y=x0[[6]], intervalo[[6]], seir1[[6]], parameter)[, 3])^2) + 
           sum((data[, 3] - ode(y=x0[[6]], intervalo[[6]], seir1[[6]], parameter)[, 4])^2) + 
           sum((data[, 4] - ode(y=x0[[6]], intervalo[[6]], seir1[[6]], parameter)[, 5])^2))
}
minimun[[7]] <- function(data, parameter){
  return(sum((data[, 1] - ode(y=x0[[7]], intervalo[[7]], seir1[[7]], parameter)[, 2])^2) + 
           sum((data[, 2] - ode(y=x0[[7]], intervalo[[7]], seir1[[7]], parameter)[, 3])^2) + 
           sum((data[, 3] - ode(y=x0[[7]], intervalo[[7]], seir1[[7]], parameter)[, 4])^2) + 
           sum((data[, 4] - ode(y=x0[[7]], intervalo[[7]], seir1[[7]], parameter)[, 5])^2))
}
minimun[[8]] <- function(data, parameter){
  return(sum((data[, 1] - ode(y=x0[[8]], intervalo[[8]], seir1[[8]], parameter)[, 2])^2) + 
           sum((data[, 2] - ode(y=x0[[8]], intervalo[[8]], seir1[[8]], parameter)[, 3])^2) + 
           sum((data[, 3] - ode(y=x0[[8]], intervalo[[8]], seir1[[8]], parameter)[, 4])^2) + 
           sum((data[, 4] - ode(y=x0[[8]], intervalo[[8]], seir1[[8]], parameter)[, 5])^2))
}

min <- list()

for (j in 1:length(intervalo)) {
  min[[j]] <- minimun[[j]](data = data.frame(cbind(susceptibles_bogota[intervalo[[j]]], 
                                       expuestos_bogota[intervalo[[j]]], 
                                       suavizado_infectados$X1[intervalo[[j]]], 
                                       suavizado_recuperados$X1[intervalo[[j]]])), parameter = 0.000001)
}

result <- list()
for (k in 1:length(intervalo)) {
  result[[k]] <- optim(par = c(0), minimun[[k]], data = data.frame(cbind(susceptibles_bogota[intervalo[[k]]], 
                                                                         expuestos_bogota[intervalo[[k]]], 
                                                                         suavizado_infectados$X1[intervalo[[k]]], 
                                                                         suavizado_recuperados$X1[intervalo[[k]]])), 
                       method = "Brent", lower = 0, upper = 1) 
}

soluciones2 <- list()

for (w in 1:length(intervalo)) {
  soluciones2[[w]] <- ode(y=x0[[w]], intervalo[[w]], seir1[[w]], result[[w]]$par)
}

soluciones_2 <- rbind(soluciones2[[1]], soluciones2[[2]], soluciones2[[3]], soluciones2[[4]], soluciones2[[5]], 
                       soluciones2[[6]], soluciones2[[7]], soluciones2[[8]])
plot(soluciones_2[, 3], lwd = 2, type = "l")
points(suavizado_recuperados$X1, pch = 16, col = "red")
1.642049e+01
3.570945e+01