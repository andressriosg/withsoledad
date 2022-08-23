library(bbmle)
library(deSolve)

intervalo <- list(1:150, 150:177, 177:205, 205:270, 270:312, 312:330, 330:363, 363:385)
x0 <- list(c(susceptibles_bogota[1], expuestos_bogota[1], 1, 0), 
           NULL, NULL, NULL, NULL, NULL, NULL, NULL)
ajuste <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
parametermv <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
theta <- list()
solucionesmv <- list()
seir1 <- list()

for (j in 2:length(intervalo)) {
  x0[[j]] <- c(susceptibles_bogota[length(intervalo[[j-1]])],
               expuestos_bogota[length(intervalo[[j-1]])],
               suavizado_infectados$X1[intervalo[[j-1]][length(intervalo[[2-1]])]],
               suavizado_reuperados$X1[intervalo[[j-1]][length(intervalo[[2-1]])]])
  }

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

  maximaver[[j]] <- function(lgamma) {
    parameter <- c(gamma = plogis(lgamma))
    weeks = intervalo[[j]]
    out <- ode(y= x0[[j]], intervalo[[j]], seir1[[j]], 0.1)
    SD1 <- sqrt(sum((suavizado_infectados$X1[weeks]-out[,4])^2)/length(weeks)) 
  }
  
  ajuste[[j]] <- mle2(maximaver[[j]],
                      start=list(lgamma= 1),  
                      method="Brent",
                      control=list(maxit=1E5,trace=0),
                      trace=FALSE, lower = 0, upper = 3)
  
  theta[[j]] <- as.numeric(c(coef(ajuste[[j]])[1]))
  parametermv[[j]] <- c(gamma = theta[[j]][1])
  
  solucionesmv[[j]] <- ode(y = x0[[j]], intervalo[[j]], seir1[[j]], as.numeric(parametermv[[j]]))
}

soluciones_mv <- rbind(solucionesmv[[1]], solucionesmv[[2]], solucionesmv[[3]], solucionesmv[[4]], solucionesmv[[5]], 
                    solucionesmv[[6]], solucionesmv[[7]], solucionesmv[[8]])

plot(soluciones_mv[, 2], type = "p", col = "red", lwd = 2, pch = 16)
plot(suavizado_recuperados$X1, col = "blue", pch = 15, lwd= 2)

### Máxima de máxima verosimilitud por intervalos 

# Estimacion de los par?metros por intervalos 

library(bbmle)
library(deSolve)

intervalo <- list(1:150, 150:205, 205:270, 270:312, 312:363, 363:385)
x0 <- list(c(susceptibles_bogota[1], expuestos_bogota[1], 1, 0),
           NULL, NULL, NULL, NULL, NULL)
ajuste <- list(NULL, NULL, NULL, NULL, NULL, NULL)
parameter <- list(NULL, NULL, NULL, NULL, NULL, NULL)
maximaver <- list(NULL, NULL, NULL, NULL, NULL, NULL)
theta <- list(NULL, NULL, NULL, NULL, NULL, NULL)
likelihood <- list()
beta <- seq(0, 0.999, by = 0.01) 
solutions <- list()
for (i in 1:length(intervalo)) {
  solutions[[i]] <- matrix(NA, ncol = 5, nrow = length(intervalo[[i]])) 
}

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

require(bbmle)
maximaver[[1]] <- function(lbeta, lnu, lgamma) {
  parms <- c(eta, beta=plogis(lbeta), nu=plogis(lnu), gamma = plogis(lgamma), mu)
  x0 <- x0[[1]]
  out <- ode(y=x0, intervalo[[1]], seir, parms)
  SD <- sqrt(sum((suavizado_infectados$X1[intervalo[[1]]]-out[,4])^2)/length(intervalo[[1]]))
  - sum(dnorm(suavizado_infectados$X1[intervalo[[1]]], mean=out[,4], sd=SD, log=TRUE))
  + sqrt(sum((suavizado_recuperados$X1[intervalo[[1]]]-out[,5])^2)/length(intervalo[[1]]))
  - sum(dnorm(suavizado_recuperados$X1[intervalo[[1]]], mean=out[,5], sd=SD, log=TRUE))
}

ajuste[[1]] <- mle2(maximaver[[1]], 
            start=list(lbeta=qlogis(0.05), 
                       lnu=qlogis(.2),
                       lgamma=qlogis(.2)),  
            method="Nelder-Mead",
            control=list(maxit=1E5,trace=0),
            trace=FALSE)

summary(ajuste[[1]])

for (j in 2:length(intervalo)) {
  
  x0[[j]] <- c(susceptibles_bogota[intervalo[[j-1]][length(intervalo[[j-1]])]], 
               expuestos_bogota[intervalo[[j-1]][length(intervalo[[j-1]])]], 
               suavizado_infectados$X1[intervalo[[j-1]][length(intervalo[[j-1]])]], 
               suavizado_recuperados$X1[intervalo[[j-1]][length(intervalo[[j-1]])]])
  
  maximaver[[j]] <- function(lbeta, lnu, lgamma) {
    parms <- c(eta, beta=plogis(lbeta), nu=plogis(lnu), gamma = plogis(lgamma), mu)
    x0 <- x0[[j]]
    out <- ode(y=x0, intervalo[[j]], seir, parms)
    SD <- sqrt(sum((suavizado_infectados$X1[intervalo[[j]]]-out[,4])^2)/length(intervalo[[j]]))
    - sum(dnorm(suavizado_infectados$X1[intervalo[[j]]], mean=out[,4], sd=SD, log=TRUE))
    + sqrt(sum((suavizado_recuperados$X1[intervalo[[j]]]-out[,5])^2)/length(intervalo[[j]]))
    - sum(dnorm(suavizado_recuperados$X1[intervalo[[j]]], mean=out[,5], sd=SD, log=TRUE))
  }
  
  ajuste[[j]] <- mle2(maximaver[[j]], 
                      start=list(lbeta=qlogis(0.01), 
                                 lnu=qlogis(.5), 
                                 lgamma=qlogis(.1)),   
                      method="Nelder-Mead",
                      control=list(maxit=1E5,trace=0),
                      trace=FALSE)
  
  theta[[j]] <- as.numeric(c(plogis(coef(ajuste[[j]])[1:3])))
  parameter[[j]] <- c(eta, beta= theta[[j]][1], nu = theta[[j]][2], gamma = theta[[j]][3], mu)
  
  solutions[[j]] <- ode(y=x0[[j]], intervalo[[j]], seir, parameter[[j]])
  colnames(solutions[[j]]) <- c("time","S", "E", "I","R")
}

theta[[1]] <- as.numeric(c(plogis(coef(ajuste[[1]])[1:3])))
parameter[[1]] <- c(eta, beta= theta[[1]][1], nu = theta[[1]][2], gamma = theta[[1]][3], mu)
solutions[[1]] <- ode(y=x0[[1]], intervalo[[1]], seir, parameter[[1]])

ajuste[[5]] <- mle2(maximaver[[5]], 
                    start=list(lbeta=qlogis(0.1), 
                               lnu=qlogis(.5), 
                               lgamma=qlogis(.8)),   
                    method="Nelder-Mead",
                    control=list(maxit=1E5,trace=0),
                    trace=FALSE)

theta[[5]] <- as.numeric(c(plogis(coef(ajuste[[5]])[1:3])))
parameter[[5]] <- c(eta, beta= theta[[5]][1], nu = theta[[5]][2], gamma = theta[[5]][3], mu)

solutions[[5]] <- ode(y=x0[[5]], intervalo[[5]], seir, parameter[[5]])

plot(tiempo[1:385], suavizado_infectados$X1[1:385], col = "blue", lwd = 2, type = "l", ylab = "Number of cases", 
     xlab = "Number of days")
points(tiempo[1:385], casos$Infectados[1:385], pch = 20)
for (k in 1:length(intervalo)) {
  lines(intervalo[[k]], solutions[[k]][, 4], col = k, type = "l", lwd = 2) 
}

solutions_1 <- rbind(solutions[[1]], solutions[[2]], solutions[[3]], solutions[[4]], solutions[[5]], solutions[[6]])
time = c(1:150, 150:205, 205:270, 270:312, 312:363, 363:385)
solo_infectados <- data.frame(time = c(tiempo[1:150], tiempo[150:205], tiempo[205:270], tiempo[270:312], tiempo[312:363], 
                                       tiempo[363:385]), 
                              solutions_1[, 2:5], suavizado_infectados$X1[time], casos$Infectados[time], 
                              suavizado_recuperados$X1[time], casos$recuperados[time])

r0_solo <- list()
for (j in 1:length(intervalo)) {
  r0_solo[[j]] = (eta*parameter[[j]][2]*parameter[[j]][3])/(parameter[[j]][5]*(parameter[[j]][3]-parameter[[j]][5])*(parameter[[j]][4]-parameter[[j]][5]))
}


### Sólo sobre los infectados 
colors1 <- c("Smoothed infected" = "firebrick3", "Estimated model" = "black", "Infected data" = "indianred2")
colors2 <- c('Smoothed recovered' = "darkgreen",  "Estimated model" = "black", "Recovered data" = "chartreuse3")

p_s_i <- ggplot(data = solo_infectados) + 
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 7], color = 'Infected data'), size = 1.5) + 
  geom_line(aes(x = solo_infectados[, 1], y = solo_infectados[, 6], color = 'Smoothed infected'), size = 1.5) +
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 4], color = 'Estimated model'), size = 0.8)  +
  labs(x = 'Time (days)', y = 'No. of infected population') + scale_color_manual(values = colors1, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
    linetype = c("solid", "solid", "blank"),
    shape = c(NA, NA, 16)))) + 
  geom_vline(xintercept = tiempo[150]) +
  geom_vline(xintercept = tiempo[205]) +
  geom_vline(xintercept = tiempo[270]) +
  geom_vline(xintercept = tiempo[312]) +
  geom_vline(xintercept = tiempo[363]) +
  theme(legend.position="right") + scale_x_date(breaks = "2 months", date_labels = "%d-%m-%Y")

p_s_r <- ggplot(data = solo_infectados) + 
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 9], color = 'Recovered data'), size = 1.5) + 
  geom_line(aes(x = solo_infectados[, 1], y = solo_infectados[, 8], color = 'Smoothed recovered'), size = 1.5) +
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 5], color = 'Estimated model'), size = 0.8)  +
  labs(x = 'Time (days)', y = 'No. of recovered population') + scale_color_manual(values = colors2, name = "", 
                                guide = guide_legend(override.aes = list(linetype = c("solid", "solid", "blank"),
                                                                    shape = c(NA, NA, 16)))) + 
  geom_vline(xintercept = tiempo[150]) +
  geom_vline(xintercept = tiempo[205]) +
  geom_vline(xintercept = tiempo[270]) +
  geom_vline(xintercept = tiempo[312]) +
  geom_vline(xintercept = tiempo[363]) +
  theme(legend.position="right") + scale_x_date(breaks = "2 months", date_labels = "%d-%m-%Y")

grid.arrange(p_s_i, p_s_r, ncol = 1)  

susceptibles_solo_i <- data.frame(tiempo[1:150], solutions_1[1:150, 2:3])

colors3 <- c("Susceptible population" = "dodgerblue", 'Exposed population' = "chocolate1")
p_s_s <- ggplot(data = susceptibles_solo_i) + 
  geom_line(aes(x = susceptibles_solo_i[, 1],  y = susceptibles_solo_i[, 2], color = 'Susceptible population'), size = 1.1) +
  geom_line(aes(x = susceptibles_solo_i[, 1],  y = susceptibles_solo_i[, 3], color = 'Exposed population'), size = 1.1) + 
  labs(x = 'Time (days)', y = 'Number of population') + scale_color_manual(values = colors3, name = "") + 
  theme(legend.position="right") + scale_x_date(breaks = "1 month", date_labels = "%d-%m-%Y")
grid.arrange(p_s_s, ncol = 1)  

max(susceptibles_solo_i)
mean(parameter[[1]][4], parameter[[2]][4], parameter[[3]][4], parameter[[4]][4], parameter[[5]][4], parameter[[6]][4])
mean(parameter[[1]][3], parameter[[2]][3], parameter[[3]][3], parameter[[4]][3], parameter[[5]][3], parameter[[6]][3])
beta_solo_i <- seq(0, 1, by = 0.01)
max_beta <- NULL
for (i in 1:length(beta_solo_i)) {
  max_beta[i] <- maximaver[[4]](lbeta = qlogis(beta_solo_i[i]), lnu = qlogis(parameter[[2]][3]), lgamma = qlogis(parameter[[2]][4]))
}

plot(beta_solo_i[2:100], max_beta[2:100], type = "l")

beta_grafica <- data.frame(beta_solo_i, max_beta)
p_beta_s <- ggplot(data = beta_grafica) + 
  geom_line(aes(x = beta_grafica[, 1],  y = beta_grafica[, 2]), col = "darkviolet", size = 1.2)  + 
  labs(x = expression(beta), y = 'Negative log-likelihood') 
grid.arrange(p_beta_s, ncol = 1)  

intervalo1 <- list(1:150, 150:177, 177:205, 205:270, 270:312, 312:330, 330:363, 363:385)
ajuste1 <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
parameter1 <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
maximaver1 <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
theta1 <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
x01 <- list(c(susceptibles_bogota[1], expuestos_bogota[1], 1, 0),
           NULL, NULL, NULL, NULL, NULL, NULL, NULL)

require(bbmle)
maximaver1[[1]] <- function(lbeta, lnu, lgamma) {
  parms <- c(eta, beta=plogis(lbeta), nu=plogis(lnu), gamma = plogis(lgamma), mu)
  x0 <- x01[[1]]
  out <- ode(y=x0, intervalo1[[1]], seir, parms)
  SD <- sqrt(sum((suavizado_infectados$X1[intervalo1[[1]]]-out[,4])^2)/length(intervalo1[[1]]))
  - sum(dnorm(suavizado_infectados$X1[intervalo1[[1]]], mean=out[,4], sd=SD, log=TRUE))
  + sqrt(sum((suavizado_recuperados$X1[intervalo1[[1]]]-out[,5])^2)/length(intervalo1[[1]]))
  - sum(dnorm(suavizado_recuperados$X1[intervalo1[[1]]], mean=out[,5], sd=SD, log=TRUE))
}
maximaver1[[1]](lbeta = 0.1, lnu = 0.1, lgamma = 0.1)

library(deSolve)
ajuste1[[1]] <- mle2(maximaver1[[1]], 
                    start=list(lbeta=qlogis(0.1), 
                               lnu=qlogis(.2),
                               lgamma=qlogis(.2)),  
                    method="Nelder-Mead",
                    control=list(maxit=1E5,trace=0),
                    trace=FALSE)

summary(ajuste1[[1]])

for (j in 2:length(intervalo1)) {
  x01[[j]] <- c(susceptibles_bogota[intervalo1[[j-1]][length(intervalo1[[j-1]])]], 
               expuestos_bogota[intervalo1[[j-1]][length(intervalo1[[j-1]])]], 
               suavizado_infectados$X1[intervalo1[[j-1]][length(intervalo1[[j-1]])]], 
               suavizado_recuperados$X1[intervalo1[[j-1]][length(intervalo1[[j-1]])]])
  
  maximaver1[[j]] <- function(lbeta, lnu, lgamma) {
    parms <- c(eta, beta=plogis(lbeta), nu=plogis(lnu), gamma = plogis(lgamma), mu)
    x0 <- x01[[j]]
    out <- ode(y=x0, intervalo1[[j]], seir, parms)
    SD <- sqrt(sum((suavizado_infectados$X1[intervalo1[[j]]]-out[,4])^2)/length(intervalo1[[j]]))
    - sum(dnorm(suavizado_infectados$X1[intervalo1[[j]]], mean=out[,4], sd=SD, log=TRUE))
    + sqrt(sum((suavizado_recuperados$X1[intervalo1[[j]]]-out[,5])^2)/length(intervalo1[[j]]))
    - sum(dnorm(suavizado_recuperados$X1[intervalo1[[j]]], mean=out[,5], sd=SD, log=TRUE))
  }
  ajuste1[[j]] <- mle2(maximaver1[[j]], 
                      start=list(lbeta=qlogis(0.01), 
                                 lnu=qlogis(.5), 
                                 lgamma=qlogis(.1)),   
                      method="Nelder-Mead",
                      control=list(maxit=1E5,trace=0),
                      trace=FALSE)
}

ajuste1[[7]] <- mle2(maximaver1[[7]], 
                     start=list(lbeta=qlogis(0.99999999), 
                                lnu=qlogis(0.7), 
                                lgamma=qlogis(0.05)),   
                     method="Nelder-Mead",
                     control=list(maxit=1E5,trace=0),
                     trace=FALSE)

solutions1 <- list()
for (i in 1:length(intervalo)) {
  solutions1[[i]] <- matrix(NA, ncol = 5, nrow = length(intervalo[[i]])) 
}
for (j in 1:length(intervalo1)) {
  theta1[[j]] <- as.numeric(c(plogis(coef(ajuste1[[j]])[1:3])))
  parameter1[[j]] <- c(eta, beta= theta1[[j]][1], nu = theta1[[j]][2], gamma = theta1[[j]][3], mu)
  solutions1[[j]] <- ode(y=x01[[j]], intervalo1[[j]], seir, parameter1[[j]])
  colnames(solutions1[[j]]) <- c("time","S", "E", "I","R")
}

solutions_2 <- rbind(solutions1[[1]], solutions1[[2]], solutions1[[3]], solutions1[[4]], solutions1[[5]], solutions1[[6]], 
                     solutions1[[7]], solutions1[[8]])
time1 <- c(1:150, 150:177, 177:205, 205:270, 270:312, 312:330, 330:363, 363:385)
solo_infectados_1 <- data.frame(time1 = c(tiempo[1:150], tiempo[150:177], tiempo[177:205], tiempo[205:270],
                                          tiempo[270:312], tiempo[312:330], tiempo[330:363], tiempo[363:385]), 
                              solutions_2[, 2:5], suavizado_infectados$X1[time1], casos$Infectados[time1], 
                              suavizado_recuperados$X1[time1], casos$recuperados[time1])

colors1 <- c("Smoothed infected" = "firebrick3", "Estimated model" = "black", "Infected data" = "indianred2")
colors2 <- c('Smoothed recovered' = "darkgreen",  "Estimated model" = "black", "Recovered data" = "chartreuse3")

p_s_i_1 <- ggplot(data = solo_infectados_1) + 
  geom_point(aes(x = solo_infectados_1[, 1],  y = solo_infectados_1[, 7], color = 'Infected data'), size = 1.5) + 
  geom_line(aes(x = solo_infectados_1[, 1], y = solo_infectados_1[, 6], color = 'Smoothed infected'), size = 1.5) +
  geom_point(aes(x = solo_infectados_1[, 1],  y = solo_infectados_1[, 4], color = 'Estimated model'), size = 0.8)  +
  labs(x = 'Time (days)', y = 'No. of infected population') + scale_color_manual(values = colors1, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
                                                                                      linetype = c("solid", "solid", "blank"),
                                                                                      shape = c(NA, NA, 16)))) + 
  geom_vline(xintercept = tiempo[150]) +
  geom_vline(xintercept = tiempo[177]) +
  geom_vline(xintercept = tiempo[205]) +
  geom_vline(xintercept = tiempo[270]) +
  geom_vline(xintercept = tiempo[312]) +
  geom_vline(xintercept = tiempo[330]) +
  geom_vline(xintercept = tiempo[363]) + scale_x_date(breaks = "2 months", date_labels = "%d-%m-%Y")

p_s_r_1 <- ggplot(data = solo_infectados_1) + 
  geom_point(aes(x = solo_infectados_1[, 1],  y = solo_infectados_1[, 9], color = 'Recovered data'), size = 1.5) + 
  geom_line(aes(x = solo_infectados_1[, 1], y = solo_infectados_1[, 8], color = 'Smoothed recovered'), size = 1.5) +
  geom_point(aes(x = solo_infectados_1[, 1],  y = solo_infectados_1[, 5], color = 'Estimated model'), size = 0.8)  +
  labs(x = 'Time (days)', y = 'No. of infected population') + scale_color_manual(values = colors2, name = "", 
                                                                                      guide = guide_legend(override.aes = list(
                                                                                        linetype = c("solid", "solid", "blank"),
                                                                                        shape = c(NA, NA, 16)))) +  
  geom_vline(xintercept = tiempo[150]) +
  geom_vline(xintercept = tiempo[177]) +
  geom_vline(xintercept = tiempo[205]) +
  geom_vline(xintercept = tiempo[270]) +
  geom_vline(xintercept = tiempo[312]) +
  geom_vline(xintercept = tiempo[330]) +
  geom_vline(xintercept = tiempo[363]) + scale_x_date(breaks = "2 months", date_labels = "%d-%m-%Y")

grid.arrange(p_s_i_1, p_s_r_1, ncol = 1)  

max1_beta <- NULL
for (i in 1:length(beta_solo_i)) {
  max1_beta[i] <- maximaver1[[7]](lbeta = qlogis(parameter1[[7]][2]), lnu = qlogis(parameter1[[7]][3]), lgamma = qlogis(beta_solo_i[i]))
}

beta1_grafica <- data.frame(beta_solo_i[2:100], max1_beta[2:100])
p1_beta_s <- ggplot(data = beta1_grafica) + 
  geom_line(aes(x = beta1_grafica[, 1],  y = beta1_grafica[, 2]), col = "darkviolet", size = 1.2)  + 
  labs(x = expression(beta), y = 'Negative log-likelihood') 
grid.arrange(p1_beta_s, ncol = 1)  


### Sólo sobre los infectados 
colors1 <- c("Smoothed infected" = "firebrick3", "Estimated model" = "black", "Infected data" = "deeppink3")
colors2 <- c("Recovered data" = "green3", 'Smoothed recovered' = "forestgreen", "Estimated model" = "black")

p_s_i <- ggplot(data = solo_infectados) + 
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 7], color = 'Infected data'), size = 1.5) + 
  geom_line(aes(x = solo_infectados[, 1], y = solo_infectados[, 6], color = 'Smoothed infected'), size = 1.5) +
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 4], color = 'Estimated model'), size = 0.5)  +
  labs(x = 'Time (days)', y = 'Number of infected population') + scale_color_manual(values = colors1) + 
  scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p_s_r <- ggplot(data = solo_infectados) + 
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 9], color = 'Recovered data'), size = 1.5) + 
  geom_line(aes(x = solo_infectados[, 1], y = solo_infectados[, 8], color = 'Smoothed recovered'), size = 1.5) +
  geom_point(aes(x = solo_infectados[, 1],  y = solo_infectados[, 5], color = 'Estimated model'), size = 0.5)  +
  labs(x = 'Time (days)', y = 'Number of infected population') + scale_color_manual(values = colors2) + 
  scale_x_date(breaks = "1 months", date_labels = "%d-%m-%Y")

grid.arrange(p_s_i,nrow = 1)  

susceptibles_solo_i_1 <- data.frame(tiempo[10:150], solutions_2[10:150, 2:3])

colors4 <- c("Susceptible population" = "dodgerblue")
colors5 <- c('Exposed population' = "chocolate1")
p_s_s_2 <- ggplot(data = susceptibles_solo_i_1) + 
  geom_line(aes(x = susceptibles_solo_i_1[, 1],  y = susceptibles_solo_i_1[, 2], color = 'Susceptible population'), size = 1.1) + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colors4, name = " ") +
  theme(legend.position="bottom")  + scale_x_date(breaks = "45 days", date_labels = "%d-%m-%Y")

p_s_e_2 <- ggplot(data = susceptibles_solo_i_1) + 
  geom_line(aes(x = susceptibles_solo_i_1[, 1],  y = susceptibles_solo_i_1[, 3], color = 'Exposed population'), size = 1.1) + 
  labs(x = 'Time (days)', y = 'No. of exposed population') + scale_color_manual(values = colors5, name = " ") + 
  theme(legend.position="bottom") +  scale_x_date(breaks = "45 days", date_labels = "%d-%m-%Y")

grid.arrange(p_s_s_2, p_s_e_2, nrow = 1)  

max(susceptibles_solo_i_1[, 3])
min(susceptibles_solo_i_1[, 2])
  
ajuste1[[1]] <- mle2(maximaver1[[1]], 
                      start=list(lbeta=qlogis(1.738857e-05), 
                                 lnu=qlogis(1.383661e-05), 
                                 lgamma=qlogis(1.633080e-02)),   
                      method="Nelder-Mead",
                      control=list(maxit=1E5,trace=0),
                      trace=FALSE)
ajuste1[[2]] <- mle2(maximaver1[[2]], 
                     start=list(lbeta=qlogis(9.882397e-01), 
                                lnu=qlogis(1.201782e-11), 
                                lgamma=qlogis(5.492948e-03)),   
                     method="Nelder-Mead",
                     control=list(maxit=1E5,trace=0),
                     trace=FALSE)
ajuste1[[3]] <- mle2(maximaver1[[3]], 
                     start=list(lbeta=qlogis(0.9), 
                                lnu=qlogis(0.00005), 
                                lgamma=qlogis(0.001)),   
                     method="Nelder-Mead",
                     control=list(maxit=1E5,trace=0),
                     trace=FALSE)
ajuste1[[4]] <- mle2(maximaver1[[4]], 
                     start=list(lbeta=qlogis(0.01), 
                                lnu=qlogis(0.2), 
                                lgamma=qlogis(0.2)),   
                     method="Nelder-Mead",
                     control=list(maxit=1E5,trace=0),
                     trace=FALSE)
ajuste1[[5]] <- mle2(maximaver1[[5]], 
                     start=list(lbeta=qlogis(0.1), 
                                lnu=qlogis(0.5), 
                                lgamma=qlogis(0.8)),   
                     method="Nelder-Mead",
                     control=list(maxit=1E5,trace=0),
                     trace=FALSE)
ajuste1[[6]] <- mle2(maximaver1[[6]], 
                     start=list(lbeta=qlogis(7.040308e-01), 
                                lnu=qlogis(4.579072e-06), 
                                lgamma=qlogis(3.334993e-03)),   
                     method="Nelder-Mead",
                     control=list(maxit=1E5,trace=0),
                     trace=FALSE)
  
for (j in 1:length(intervalo1)) {
  theta1[[j]] <- as.numeric(c(plogis(coef(ajuste1[[j]])[1:3])))
  parameter1[[j]] <- c(eta, beta= theta1[[j]][1], nu = theta1[[j]][2], gamma = theta1[[j]][3], mu)
  solutions1[[j]] <- ode(y=x0[[j]], intervalo1[[j]], seir, parameter[[j]])
  colnames(solutions1[[j]]) <- c("time","S", "E", "I","R")
  }

plot(tiempo[1:385], suavizado_infectados$X1[1:385], col = "blue", lwd = 2, type = "l", ylab = "Number of cases", 
     xlab = "Number of days")
points(tiempo[1:385], casos$Infectados[1:385], pch = 20)
for (k in 1:length(intervalo1)) {
  lines(intervalo[[k]], solutions1[[k]][, 4], col = k, type = "l", lwd = 2) 
}

plot(tiempo[4:100], solutions1[[1]][, 2][4:100], col = "blue", lwd = 2, type = "l", ylab = "Number of cases", 
     xlab = "Number of days")
plot(tiempo[1:100], solutions1[[1]][, 3][1:100], col = "red")
points(tiempo[1:385], casos$Infectados[1:385], pch = 20)
for (k in 2:length(intervalo)) {
  lines(intervalo[[k]], solutions1[[k]][, 2], col = k, type = "l", lwd = 2) 
}
