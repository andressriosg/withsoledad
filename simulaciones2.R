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

x0 <- c(7743953, 0, 1, 0)
times = 1:419
parametros = c((149090/7743955)*7743955, 4.246413e-06, 1/5.1, 15.0715, 84643/7743955)  

# install.packages("deSolve")
library(deSolve)
solution_1 <- ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, 0.2, 0.1, 0.3, 0.2))

mu = 0.2
eta = 4
upsilon = 0.1
beta = ((1-mu)*sum(((solution_1[1:99, 2])^2*solution_1[1:99, 4])) - sum(solution_1[2:100, 2]*solution_1[1:99, 2]*solution_1[1:99, 4])  + eta*sum(solution_1[1:99, 2]*solution_1[1:99, 4]) - (1 - upsilon - mu)*sum(solution_1[1:99, 2]*solution_1[1:99, 3]*solution_1[1:99, 4]) + sum(solution_1[2:100, 3]*solution_1[1:99, 2]*solution_1[1:99, 4]))/(2*sum((solution_1[1:99, 2])^2*(solution_1[1:99, 4])^2))

gamma = ((1-mu)*(sum((solution_1[1:99, 4])^2) - sum(solution_1[1:99, 5]*solution_1[1:99, 4])) - sum(solution_1[1:99, 4]*solution_1[2:100, 4]) + sum(solution_1[1:99, 4]*solution_1[2:100, 5]) + upsilon*sum(solution_1[1:99, 3]*solution_1[1:99, 4]))/(2*sum((solution_1[1:99, 4])^2))

plot(solution_1[, 4], col = "red", type = "p", lwd = 2, pch = 16)
lines(ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, beta, 0.1, gamma, 0.2))[,4], col = "blue", type = "l", lwd = 2)

plot(solution_1[, 5], col = "green", type = "p", lwd = 2, pch = 16)
lines(ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, beta, 0.1, gamma, 0.2))[,5], col = "blue", type = "l", lwd = 2)

# Conclusión: El método de estimación funciona perfecto. El problema está en los datos 

# Veamos con un modelo con perturbaciones aleatorias 

library(deSolve)
solution_1 <- ode(y= c(4, 0, 0.1, 0), 0:101, seir, c(4, 0.2, 0.1, 0.3, 0.2))

solution_2 <- matrix(0, nrow = 102, ncol = 5)
solution_2[, 1:3] = solution_1[, 1:3]
for (j in 4:5) {
  for (i in 1:51) {
    #    solution_2[2*i, 2] = solution_1[2*i, 2] + 1
    #    solution_2[2*i, 3] = solution_1[2*i, 3] + 0.5
    solution_2[2*i, j] = solution_1[2*i, j] + 0.2
    #    solution_2[(2*i - 1), 2] = solution_1[(2*i - 1), 2] - 1
    #    solution_2[(2*i - 1), 3] = solution_1[(2*i - 1), 3] - 0.5
    solution_2[(2*i - 1), j] = solution_1[(2*i - 1), j] - 0.2
  } 
}
solution_2[1, 2:5] <- c(4, 0, 0.1, 0)
solution_2[1, ] = solution_1[1, ]

plot(solution_2[, 4], type = "l", col = "red", lwd = 2)
plot(solution_2[, 5], type = "l", col = "darkgreen", lwd = 2)

eta = 4
mu = 0.2
upsilon = 0.1

beta_1 = ((1-mu)*sum(((solution_2[1:99, 2])^2*solution_2[1:99, 4])) - sum(solution_2[2:100, 2]*solution_2[1:99, 2]*solution_2[1:99, 4])  + eta*sum(solution_2[1:99, 2]*solution_2[1:99, 4]) - (1 - upsilon - mu)*sum(solution_2[1:99, 2]*solution_2[1:99, 3]*solution_2[1:99, 4]) + sum(solution_2[2:100, 3]*solution_2[1:99, 2]*solution_2[1:99, 4]))/(2*sum((solution_2[1:99, 2])^2*(solution_2[1:99, 4])^2))

gamma_1 = ((1-mu)*(sum((solution_2[1:99, 4])^2) - sum(solution_2[1:99, 5]*solution_2[1:99, 4])) - sum(solution_2[1:99, 4]*solution_2[2:100, 4]) + sum(solution_2[1:99, 4]*solution_2[2:100, 5]) + upsilon*sum(solution_2[1:99, 3]*solution_2[1:99, 4]))/(2*sum((solution_2[1:99, 4])^2))

plot(solution_2[, 4], col = "red", type = "p", lwd = 2, pch = 16)
lines(ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, beta_1, 0.1, gamma_1, 0.2))[,4], col = "blue", type = "l", lwd = 2)

plot(solution_2[, 5], col = "green", type = "p", lwd = 2, pch = 16)
lines(ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, beta_1, 0.1, gamma_1, 0.2))[,5], col = "blue", type = "l", lwd = 2)

plot(solution_2[, 2], col = "green", type = "p", lwd = 2, pch = 16)
lines(ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, beta_1, 0.1, gamma_1, 0.2))[,2], col = "blue", type = "l", lwd = 2)

solution_d = ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, beta_1, 0.1, gamma_1, 0.2))

# ¿Estima bien el parámetro sigma?

# PARA ESTIMAR TODOS LOS PARÁMETROS DEL MODELO

# Primero escribimos los datos dados por Sj e Ij: 

S = NULL 
E = NULL 
I = NULL 
S[1] = solution_2[1, 2]
E[1] = solution_2[1, 3] 
I[1] = solution_2[1, 4]
for(j in 1:99){
  S[j+1] = solution_2[j, 2] + eta - beta_1*solution_2[j, 2]*solution_2[j, 4] - mu*solution_2[j, 2]
  E[j+1] = solution_2[j, 3] + beta_1*solution_2[j, 2]*solution_2[j, 4] - upsilon*solution_2[j, 3] - mu*solution_2[j, 3]
}

plot(solution_2[, 2], col = "blue", type = "p", lwd = 2, pch = 16)
lines(S, col = "cyan", lwd = 2, pch = 16)

plot(solution_2[, 3], col = "orange", type = "p", lwd = 2, pch = 18)
points(E, col = "blue", lwd = 2, pch = 18)

# PRIMERO: CONSIDERANDO QUE SE CONOCE LA SOLUCION REAL
sigma2 = (1/(2*99))*(sum(((solution_2[2:100, 2]-solution_1[2:100, 2])^2)/((solution_1[1:99, 2])^2*(solution_1[1:99, 4])^2)) + sum(((solution_2[2:100, 3]-solution_1[2:100, 3])^2)/((solution_1[1:99, 2]*solution_1[1:99, 4])^2)))
sigma = sqrt(sigma2)

SN = list() 
EN = list()
IN = list()
RN = list()
z = list()

for (k in 1:5) {
  SN[[k]] = rep(0, 100)
  EN[[k]] = rep(0, 100)
  IN[[k]] = rep(0, 100)
  RN[[k]] = rep(0, 100)
  z[[k]] = as.vector(rnorm(100, 0, 1))
}

for (k in 1:5) {
  SN[[k]][1] = solution_2[1, 2]
  EN[[k]][1] = solution_2[1, 3]
  IN[[k]][1] = solution_2[1, 4]
  RN[[k]][1] = solution_2[1, 5]
}

for (k in 1:5) {
  for(j in 1:99){
    SN[[k]][j+1] = SN[[k]][j] + (eta - beta_1*SN[[k]][j]*IN[[k]][j] - mu*SN[[k]][j]) - sigma*SN[[k]][j]*IN[[k]][j]*z[[k]][j]
    EN[[k]][j+1] = EN[[k]][j] + (beta_1*SN[[k]][j]*IN[[k]][j] - upsilon*EN[[k]][j] - mu*EN[[k]][j]) + sigma*SN[[k]]*IN[[k]]*z[[k]][j]
    IN[[k]][j+1] = IN[[k]][j] + (upsilon*EN[[k]][j] - gamma_1*IN[[k]][j] - mu*IN[[k]][j])
    RN[[k]][j+1] = RN[[k]][j] + (gamma_1*IN[[k]][j] - mu*RN[[k]][j])
}}

warnings()
plot(solution_2[, 5], col = "green", type = "p", lwd = 2, pch = 18)
lines(RN, col = "blue", lwd = 3, pch = 18)

plot(solution_2[, 2], col = "blue", type = "p", lwd = 2, pch = 18)
lines(SN, col = "cyan", lwd = 3, pch = 18)

plot(solution_2[, 4], col = "red", type = "p", lwd = 2, pch = 18)
lines(IN, col = "purple", lwd = 3, pch = 18)

plot(solution_2[, 3], col = "black", type = "p", lwd = 3, pch = 18, ylim = c(0, 9))
lines(EN, col = "orange", lwd = 3, pch = 18)

install.packages("ggplot2")
library(ggplot2)

graf_suscp <- data.frame(solution_2[, 1:2], SN[[1]], SN[[2]], SN[[3]], SN[[4]], SN[[5]]) 
graf_expue <- data.frame(solution_2[, 1], solution_2[, 3], EN[[1]], EN[[2]], EN[[3]], EN[[4]], EN[[5]])
graf_infec <- data.frame(solution_2[, 1], solution_2[, 4], IN[[1]], IN[[2]], IN[[3]], IN[[4]], IN[[5]])
graf_recup <- data.frame(solution_2[, 1], solution_2[, 5], RN[[1]], RN[[2]], RN[[3]], RN[[4]], RN[[5]])

p1 <- ggplot(data = graf_suscp) + geom_point(aes(x = graf_suscp[, 1], y = graf_suscp[, 2]), col = 'darkblue', 
      size = 2) + labs(x = 'Time (days)', y = 'Number of susceptible population') + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 3]), col = 'dodgerblue', size = 1) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 4]), col = 'dodgerblue1', size = 1) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 5]), col = 'deepskyblue', size = 1) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 6]), col = 'deepskyblue1', size = 1) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 7]), col = 'cyan3', size = 1)

p2 <- ggplot(data = graf_expue) + geom_point(aes(x = graf_expue[, 1], y = graf_expue[, 2]), col = 'darkorange4', 
      size = 2) + labs(x = 'Time (days)', y = 'Number of exposed population') + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 3]), col = 'chocolate1', size = 1) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 4]), col = 'chocolate2', size = 1) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 5]), col = 'goldenrod2', size = 1) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 6]), col = 'goldenrod', size = 1) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 7]), col = 'darkorange', size = 1)

p3 <- ggplot(data = graf_infec) + geom_point(aes(x = graf_infec[, 1], y = graf_infec[, 2]), col = 'red4', 
      size = 2) + labs(x = 'Time (days)', y = 'Number of infected population') + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 3]), col = 'firebrick1', size = 1) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 4]), col = 'firebrick2', size = 1) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 5]), col = 'orangered1', size = 1) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 6]), col = 'orangered2', size = 1) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 7]), col = 'brown2', size = 1)

p4 <- ggplot(data = graf_recup) + geom_point(aes(x = graf_recup[, 1], y = graf_recup[, 2]), col = 'seagreen4', 
      size = 2) + labs(x = 'Time (days)', y = 'Number of recovered population') + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 3]), col = 'seagreen2', size = 1) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 4]), col = 'lawngreen', size = 1) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 5]), col = 'palegreen2', size = 1) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 6]), col = 'springgreen2', size = 1) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 7]), col = 'lightgreen', size = 1)

# install.packages("gridExtra")
library(gridExtra)
grid.arrange(p1, p2, p3, p4)

# Estimando los expuestos y los susceptibles 

pob <- NULL 
for (j in 1:length(solution_1[, 4])) {
  pob[j] = sum(solution_1[j, 2:5])
}

delta_pob <- NULL
for (j in 1:length(solution_2[, 4])) {
  delta_pob[j] = pob[j+1] - pob[j]
}

delta_inf <- NULL
for (j in 1:length(solution_2[, 4])) {
  delta_inf[j] = solution_1[, 4][j+1] - solution_1[, 4][j]
}

delta_rec <- NULL
for (j in 1:length(solution_2[, 4])) {
  delta_rec[j] = solution_1[, 5][j+1] - solution_1[, 5][j]
}

mu = 0.2
eta = 4
upsilon = 0.1
expuestos_e <- ((delta_inf[1:99] + delta_rec[1:99] - delta_pob[1:99])/upsilon) + (eta/upsilon) - 
  ((mu/upsilon)*(pob[1:99] - solution_2[1:99, 4] - solution_2[1:99, 5]))
plot(solution_2[, 3], type = "l", lwd = 2, col = "orange") ## NO SIRVE 
points(expuestos_e, type = "p", lwd = 2, col = "darkorange",  pch = 16)

susceptibles_s <- pob[1:99] - solution_2[1:99, 4] - solution_2[1:99, 5] - expuestos_e 
plot(solution_2[, 2], type = "l", lwd = 2, col = "blue") ## NO SIRVE 
points(susceptibles_s, type = "p", lwd = 2, col = "darkblue",  pch = 16)

library(ggplot2)

graf_expuestos <- data.frame(1:99, solution_2[1:99, 3], expuestos_e)
graf_susceptibles <- data.frame(1:99, solution_2[1:99, 2], susceptibles_s)


p5 <- ggplot(data = graf_susceptibles) + geom_point(aes(y = graf_susceptibles[, 2], x = graf_susceptibles[, 1]), col = 'cyan', 
                                                 size = 2) + labs(x = 'Time (days)', y = 'Number of susceptible population') + 
  geom_line(aes(y = graf_susceptibles[, 3], x = graf_susceptibles[, 1]), col = 'blue', size = 1)
p6 <- ggplot(data = graf_expuestos) + geom_point(aes(y = graf_expuestos[, 2], x = graf_expuestos[, 1]), col = 'orange', 
                                                 size = 2) + labs(x = 'Time (days)', y = 'Number of exposed population') + 
  geom_line(aes(y = graf_expuestos[, 3], x = graf_expuestos[, 1]), col = 'darkorange', size = 1)

library(gridExtra)
grid.arrange(p5, p6)
