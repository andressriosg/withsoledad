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
plot(solution_2[, 4])

# Estimación de los susceptibles y expuestos 

pob <- NULL # Población total: N(t) = S(t) + E(t) + I(t) + R(t)
for (j in 1:length(solution_1[, 4])) {
  pob[j] = sum(solution_1[j, 2:5])
}

delta_pob <- NULL
for (j in 1:length(solution_2[, 4])) {
  delta_pob[j] = pob[j+1] - pob[j]
}

delta_inf <- NULL
for (j in 1:length(solution_2[, 4])) {
  delta_inf[j] = solution_2[, 4][j+1] - solution_2[, 4][j]
}

delta_rec <- NULL
for (j in 1:length(solution_2[, 4])) {
  delta_rec[j] = solution_2[, 5][j+1] - solution_2[, 5][j]
}

mu_1 = 0.2
eta_1 = 4
upsilon_1 = 0.1
expuestos_e <- ((delta_inf[1:101] + delta_rec[1:101] - delta_pob[1:101])/upsilon_1) + (eta_1/upsilon_1) - 
  ((mu_1/upsilon_1)*(pob[1:101] - solution_2[1:101, 4] - solution_2[1:101, 5]))

susceptibles_s <- pob[1:101] - solution_2[1:101, 4] - solution_2[1:101, 5] - expuestos_e
plot(susceptibles_s)
plot(solution_2[, 2], type = "l", lwd = 2, col = "blue")
points(susceptibles_s, type = "p", lwd = 2, col = "darkblue",  pch = 16)

library(ggplot2)

graf_expuestos <- data.frame(0:100, solution_2[1:101, 3], expuestos_e)
graf_susceptibles <- data.frame(0:100, solution_2[1:101, 2], susceptibles_s)

colors4 <- c("Estimated susceptible" = "royalblue1", "Modeled susceptible" = "blue4")
colors5 <- c("Estimated exposed" = "gold2", 'Modeled exposed' = "darkorange3")
p5 <- ggplot(data = graf_susceptibles) + geom_point(aes(x = graf_susceptibles[, 1], y = graf_susceptibles[, 3], 
                                                       color = "Estimated susceptible"), size = 1.5) + 
      labs(x = 'Time (days)', y = 'No. of susceptible population') +
      geom_line(aes(x = graf_susceptibles[, 1], y = graf_susceptibles[, 2], color = "Modeled susceptible"), size = 1) +  
      scale_color_manual(values = colors4, name = "", guide = guide_legend(override.aes = list(linetype = c("blank", "solid"),
      shape = c(16, NA))))  + theme(legend.position="bottom")  
p6 <- ggplot(data = graf_expuestos) + geom_point(aes(y = graf_expuestos[, 3], x = graf_expuestos[, 1], color = "Estimated exposed"), 
      size = 1.5) + labs(x = 'Time (days)', y = 'No. of exposed population') + 
      geom_line(aes(y = graf_expuestos[, 2], x = graf_expuestos[, 1], color = 'Modeled exposed'), size = 1) +  
      scale_color_manual(values = colors5, name = "", guide = guide_legend(override.aes = list(linetype = c("blank", "solid"),
      shape = c(16, NA))))  + theme(legend.position="bottom")

library(gridExtra)
grid.arrange(p5, p6, ncol = 2)

# Suavizamiento de las funciones usando regresión no paramétrica 

tiempo = 0:(length(solution_2[, 4])-1)
library(npregfast)
suavizado_i_e <- frfast(solution_2[, 4] ~ tiempo, model = "np", smooth = "kernel", kbin = 102, p =3) # Suavizamos infectados 
suavizado_i_e <- data.frame(suavizado_i_e$p)
suavizado_i_e = suavizado_i_e$X1

suavizado_r_e <- frfast(solution_2[, 5] ~ tiempo, model = "np", smooth = "kernel", kbin = 102, p =3) # Suavizamos recuperados 
suavizado_r_e <- data.frame(suavizado_r_e$p)
suavizado_r_e = suavizado_r_e$X1

colors_i <- c("Smoothed infected" = "firebrick3", "Simulated infected data" = "indianred2")
colors_r <- c('Smoothed recovered' = "darkgreen", "Simulated recovered data" = "chartreuse3")

graf_infectados_recuperados <- data.frame(0:100, solution_2[1:101, 4:5], suavizado_i_e[1:101], suavizado_r_e[1:101])

library(ggplot2)
library(gridExtra)

p10 <- ggplot(data = graf_infectados_recuperados ) + 
       geom_point(aes(x = graf_infectados_recuperados[, 1], y = graf_infectados_recuperados[, 2], color = "Simulated infected data"), size = 1.5) + 
       geom_line(aes(x = graf_infectados_recuperados[, 1], y = graf_infectados_recuperados[, 4], color = "Smoothed infected"), size = 1.2) +
       labs(x = 'Time (days)', y = 'No. of infected population') + scale_color_manual(values = colors_i, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("solid", "blank"),
                                                                                   shape = c(NA, 16)))) + 
       theme(legend.position="bottom")
p11 <- ggplot(data = graf_infectados_recuperados ) +
       geom_point(aes(x = graf_infectados_recuperados[, 1], y = graf_infectados_recuperados[, 3], color = "Simulated recovered data"), size = 1.5) + 
       geom_line(aes(x = graf_infectados_recuperados[, 1], y = graf_infectados_recuperados[, 5], color = 'Smoothed recovered'), size = 1.2)  +
       labs(x = 'Time (days)', y = 'No. of recovered population') + scale_color_manual(values = colors_r, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("solid", "blank"),
                                                                                   shape = c(NA, 16)))) + 
       theme(legend.position="bottom")

grid.arrange(p10, p11, ncol = 2)

delta_inf_s <- NULL
for (j in 1:length(suavizado_i_e)) {
  delta_inf_s[j] = suavizado_i_e[j+1] - suavizado_i_e[j]
}

delta_rec_s <- NULL
for (j in 1:length(suavizado_r_e)) {
  delta_rec_s[j] = suavizado_r_e[j+1] - suavizado_r_e[j]
}

expuestos_e_s <- ((delta_inf_s[1:101] + delta_rec_s[1:101] - delta_pob[1:101])/upsilon_1) + (eta_1/upsilon_1) - 
  ((mu_1/upsilon_1)*(pob[1:101] - suavizado_i_e[1:101] - suavizado_r_e[1:101]))

susceptibles_s_s <- pob[1:101] - suavizado_i_e[1:101] - suavizado_r_e[1:101] - expuestos_e_s 
plot(solution_2[, 3], type = "l", lwd = 2, col = "blue")
points(expuestos_e_s, type = "p", lwd = 2, col = "darkblue",  pch = 16)

library(ggplot2)

graf_expuestos_s <- data.frame(0:100, solution_2[1:101, 3], expuestos_e_s)
graf_susceptibles_s <- data.frame(0:100, solution_2[1:101, 2], susceptibles_s_s)
colors4 <- c("Estimated susceptible" = "royalblue1", "Modeled susceptible" = "blue4")
colors5 <- c("Estimated exposed" = "gold2", 'Modeled exposed' = "darkorange3")

p7 <- ggplot(data = graf_susceptibles_s) + geom_point(aes(x = graf_susceptibles_s[, 1], y = graf_susceptibles_s[, 3], 
                                                          color = "Estimated susceptible"), size = 1.5) + 
      geom_line(aes(x = graf_susceptibles_s[, 1], y = graf_susceptibles_s[, 2], color = "Modeled susceptible"), size = 1) + 
      labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colors4, name = "", 
                                                                                        guide = guide_legend(override.aes = list(
                                                                                          linetype = c( "blank", "solid"),
                                                                                          shape = c(16, NA))))  +
      theme(legend.position="bottom")  
p8 <- ggplot(data = graf_expuestos_s) +
      geom_point(aes(y = graf_expuestos_s[, 3], x = graf_expuestos_s[, 1], color = "Estimated exposed"), size = 1.5) +
      geom_line(aes(y = graf_expuestos_s[, 2], x = graf_expuestos_s[, 1], color = 'Modeled exposed'), size = 1) +
      labs(x = 'Time (days)', y = 'No. of exposed population') + scale_color_manual(values = colors5, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
                                                                                      linetype = c("blank", "solid"),
                                                                                      shape = c(16, NA))))   + 
      theme(legend.position="bottom")

library(gridExtra)
grid.arrange(p7, p8, ncol = 2)

# Estimación sin suavización 

mu_1 = 0.2
eta_1 = 4
upsilon_1 = 0.1
beta_e = ((1-mu_1)*sum(((susceptibles_s[1:100])^2*solution_2[1:100, 4])) - sum(susceptibles_s[2:101]*susceptibles_s[1:100]*solution_2[1:100, 4])  + eta_1*sum(susceptibles_s[1:100]*solution_2[1:100, 4]) - (1 - upsilon_1 - mu_1)*sum(susceptibles_s[1:100]*expuestos_e[1:100]*solution_2[1:100, 4]) + sum(expuestos_e[2:101]*susceptibles_s[1:100]*solution_2[1:100, 4]))/(2*sum((susceptibles_s[1:100])^2*(solution_2[1:100, 4])^2))
beta_e
gamma_e = ((1-mu_1)*(sum((solution_2[1:100, 4])^2) - sum(solution_2[1:100, 5]*solution_2[1:100, 4])) - sum(solution_2[1:100, 4]*solution_2[2:101, 4]) + sum(solution_2[1:100, 4]*solution_2[2:101, 5]) + upsilon*sum(expuestos_e[1:100]*solution_2[1:100, 4]))/(2*sum((solution_1[1:100, 4])^2))
gamma_e


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
library(deSolve)
solution_3 <- ode(y= solution_2[1, 2:5], 1:100, seir, c(4, beta_e, 0.1, gamma_e, 0.2))

estimacion_final_s <- data.frame(susceptibles_s[1:100], expuestos_e[1:100], solution_2[1:100, 4], 
                                 solution_2[1:100, 5], solution_3[1:100, ], solution_1[1:100, 2:5])

colors4 <- c("Estimated susceptible" = "deepskyblue", "Solution 1" = "blue4", "Solution 2" = "royalblue1")
colors5 <- c("Estimated exposed" = "gold2", "Solution 1" = "darkorange3", "Solution 2" = "orange2")
colors1 <- c("Infected data" = "indianred2",  "Solution 1" = "firebrick3", "Solution 2" = "magenta3")
colors2 <- c("Recovered data" = "chartreuse3", "Solution 1" = "darkgreen", "Solution 2" = "springgreen3")

p12 <- ggplot(data = estimacion_final_s) + geom_point(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[, 1]
                                                          , color = "Estimated susceptible"), size = 2.3) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[,6], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[,10], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colors4, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
                                                                                      linetype = c("blank", "solid", "solid"),
                                                                                      shape = c(16, NA, NA))))  +
  theme(legend.position="bottom") 
p13 <- ggplot(data = estimacion_final_s) + geom_point(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[, 2], 
                                                          color = "Estimated exposed"), size = 2.3) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[,7], color =  "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[,11], color =  "Solution 2"), size = 1.2, linetype = "dashed") +
  labs(x = 'Time (days)', y = 'No. of exposed population')  + scale_color_manual(values = colors5, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("blank", "solid", "solid"),
                                                                                   shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")
p14 <- ggplot(data = estimacion_final_s) + geom_point(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[, 3], 
                                                          color = "Infected data"), size = 2.3) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[, 8], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[,12], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of infected population')  + scale_color_manual(values = colors1, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")
p15 <- ggplot(data = estimacion_final_s) + geom_point(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[, 4], 
                                                          color = 'Recovered data'), size = 2.3) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[, 9], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_s[, 5], y = estimacion_final_s[,13], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + scale_color_manual(values = colors2, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom") 
library(gridExtra)
grid.arrange(p12, p13, p14, p15)

# Estimación con suavización 

mu_1 = 0.2
eta_1 = 4
upsilon_1 = 0.1
beta_s = ((1-mu_1)*sum(((susceptibles_s_s[1:100])^2*suavizado_i_e[1:100])) - sum(susceptibles_s_s[2:101]*susceptibles_s_s[1:100]*suavizado_i_e[1:100])  + eta_1*sum(susceptibles_s_s[1:100]*suavizado_i_e[1:100]) - (1 - upsilon_1 - mu_1)*sum(susceptibles_s_s[1:100]*expuestos_e_s[1:100]*suavizado_i_e[1:100]) + sum(expuestos_e_s[2:101]*susceptibles_s_s[1:100]*suavizado_i_e[1:100]))/(2*sum((susceptibles_s_s[1:100])^2*(suavizado_i_e[1:100])^2))
beta_s
gamma_s = ((1-mu_1)*(sum((suavizado_i_e[1:100])^2) - sum(suavizado_r_e[1:100]*suavizado_i_e[1:100])) - sum(suavizado_i_e[1:100]*suavizado_i_e[2:101]) + sum(suavizado_i_e[1:100]*suavizado_r_e[1:100]) + upsilon*sum(expuestos_e_s[1:100]*suavizado_i_e[1:100]))/(2*sum((suavizado_i_e[1:100])^2))
gamma_s

library(deSolve)
solution_4 <- ode(y= solution_2[1,2:5], 1:101, seir, c(4, beta_s, 0.1, gamma_s, 0.2))

estimacion_final_c <- data.frame(susceptibles_s_s[1:100], expuestos_e_s[1:100], suavizado_i_e[1:100], 
                                 suavizado_r_e[1:100], solution_4[1:100, ], solution_1[1:100, 2:5], solution_2[1:100, 4], 
                                 solution_2[1:100, 5])
head(estimacion_final_c)

colors4 <- c("Estimated susceptible" = "royalblue1", "Solution 1" = "blue4", "Solution 2" = "deepskyblue")
colors5 <- c("Estimated exposed" = "orange2", "Solution 1" = "darkorange3", "Solution 2" = "gold2")
colors1 <- c("Infected data" = "indianred2",  "Solution 1" = "firebrick3", "Solution 2" = "magenta3")
colors2 <- c("Recovered data" = "chartreuse3", "Solution 1" = "darkgreen", "Solution 2" = "springgreen3")

p16 <- ggplot(data = estimacion_final_c) + geom_point(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[, 1]
                                                          , color = "Estimated susceptible"), size = 2.3) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[,6], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[,10], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colors4, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid", "solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom") 

p17 <- ggplot(data = estimacion_final_c) + geom_point(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[, 2], 
                                                          color = "Estimated exposed"), size = 2.3) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[,7], color =  "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[,11], color =  "Solution 2"), size = 1.2, linetype = "dashed") +
  labs(x = 'Time (days)', y = 'No. of exposed population')  + scale_color_manual(values = colors5, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("blank", "solid", "solid"),
                                                                                   shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")
p18 <- ggplot(data = estimacion_final_c) + geom_point(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[, 14], 
                                                         color = "Infected data"), size = 2.3) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[, 8], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[,12], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of infected population')  + scale_color_manual(values = colors1, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p19 <- ggplot(data = estimacion_final_c) + geom_point(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[, 15], 
                                                      color = 'Recovered data'), size = 2.3) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[, 9], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c[, 5], y = estimacion_final_c[,13], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + scale_color_manual(values = colors2, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")  

grid.arrange(p16, p17, p18, p19)

#### Cuando únicamente se tienen datos para infectados 

gamma_1 = 0.3

recuperados_sim <- NULL
recuperados_sim[1] <- 0 # Sin suavización de datos 
for (k in 2:101) {
  recuperados_sim[k] = recuperados_sim[k-1] + (gamma_1*solution_2[, 4][k-1] - mu_1*recuperados_sim[k-1]) 
}

expuestos_e_sin_inf <- NULL 
expuestos_e_sin_inf <- ((delta_inf[1:101] - delta_pob[1:101])/upsilon_1) + (eta_1/upsilon_1) - 
  ((mu_1/upsilon_1)*(pob[1:101] - solution_2[1:101, 4] + gamma_1*solution_2[1:101, 4]))
plot(expuestos_e_sin_inf)
lines(solution_2[,3])

susceptibles_e_sin_inf <-pob[1:101] - expuestos_e_sin_inf - solution_2[1:101, 4] - recuperados_sim

beta_e_sin = ((1-mu_1)*sum(((susceptibles_e_sin_inf[1:100])^2*solution_2[1:100, 4])) - 
            sum(susceptibles_e_sin_inf[2:101]*susceptibles_e_sin_inf[1:100]*solution_2[1:100, 4])  + 
            eta_1*sum(susceptibles_e_sin_inf[1:100]*solution_2[1:100, 4]) - (1 - upsilon_1 - mu_1)*
            sum(susceptibles_e_sin_inf[1:100]*expuestos_e_sin_inf[1:100]*solution_2[1:100, 4]) + 
            sum(expuestos_e_sin_inf[2:101]*susceptibles_e_sin_inf[1:100]*solution_2[1:100, 4]))/
            (2*sum((susceptibles_e_sin_inf[1:100])^2*(solution_2[1:100, 4])^2))

solution_33 <- ode(y= solution_2[1, 2:5], 1:100, seir, c(4, beta_e_sin, 0.1, 0.3, 0.2))

estimacion_final_sin_rec <- data.frame(susceptibles_e_sin_inf[1:100], expuestos_e_sin_inf[1:100], 
                                       solution_2[1:100, 4], recuperados_sim[1:100], 
                                       solution_33[1:100, ], solution_1[1:100, 2:5])

colors4 <- c("Estimated susceptible" = "deepskyblue", "Solution 1" = "blue4", "Solution 2" = "royalblue1")
colors5 <- c("Estimated exposed" = "gold2", "Solution 1" = "darkorange3", "Solution 2" = "orange2")
colors1 <- c("Infected data" = "indianred2",  "Solution 1" = "firebrick3", "Solution 2" = "magenta3")
colors2 <- c('Estimated recovered' = "springgreen3", "Solution 1" = "chartreuse3", "Solution 2" = "darkgreen")

p12_12 <- ggplot(data = estimacion_final_sin_rec) + geom_point(aes(x = estimacion_final_sin_rec[, 5], 
                                                                y = estimacion_final_sin_rec[, 1]
                                                          , color = "Estimated susceptible"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[,6], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[,10], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colors4, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
                                                                                      linetype = c("blank", "solid", "solid"),
                                                                                      shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p13_13 <- ggplot(data = estimacion_final_sin_rec) + geom_point(aes(x = estimacion_final_sin_rec[, 5], 
                                                                y = estimacion_final_sin_rec[, 2], 
                                                          color = "Estimated exposed"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[,7], color =  "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[,11], color =  "Solution 2"), size = 1.2, linetype = "dashed") +
  labs(x = 'Time (days)', y = 'No. of exposed population')  + scale_color_manual(values = colors5, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("blank", "solid", "solid"),
                                                                                   shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p14_14 <- ggplot(data = estimacion_final_sin_rec) + geom_point(aes(x = estimacion_final_sin_rec[, 5], 
                                                                   y = estimacion_final_sin_rec[, 3], 
                                                          color = "Infected data"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[, 8], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[,12], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of infected population')  + scale_color_manual(values = colors1, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p15_15 <- ggplot(data = estimacion_final_sin_rec) + geom_point(aes(x = estimacion_final_sin_rec[, 5], 
                                                                   y = estimacion_final_sin_rec[, 4], 
                                                          color = 'Estimated recovered'), size = 1.5) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[, 9], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_sin_rec[, 5], y = estimacion_final_sin_rec[,13], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + scale_color_manual(values = colors2, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom") 
library(gridExtra)
grid.arrange(p12_12, p13_13, p14_14, p15_15)

recuperados_sim_1 <- NULL
recuperados_sim_1[1] <- 0 # Con suavización de datos 
for (k in 2:101) {
  recuperados_sim_1[k] = recuperados_sim_1[k-1] + (gamma_1*suavizado_i_e[k-1] - mu_1*recuperados_sim_1[k-1]) 
}

expuestos_e_sin_inf_suave <- NULL 
expuestos_e_sin_inf_suave <- ((delta_inf_s[1:101] - delta_pob[1:101])/upsilon_1) + (eta_1/upsilon_1) - 
  ((mu_1/upsilon_1)*(pob[1:101] - suavizado_i_e[1:101] + gamma_1*suavizado_i_e[1:101]))
plot(expuestos_e_sin_inf_suave)
lines(solution_2[,3])

susceptibles_e_sin_inf_suave <- pob[1:101] - expuestos_e_sin_inf_suave - suavizado_i_e[1:101] - recuperados_sim_1

beta_e_sin_suave = ((1-mu_1)*sum(((susceptibles_e_sin_inf_suave[1:100])^2*suavizado_i_e[1:100])) - 
                sum(susceptibles_e_sin_inf_suave[2:101]*susceptibles_e_sin_inf_suave[1:100]*suavizado_i_e[1:100])  + 
                eta_1*sum(susceptibles_e_sin_inf_suave[1:100]*suavizado_i_e[1:100]) - (1 - upsilon_1 - mu_1)*
                sum(susceptibles_e_sin_inf_suave[1:100]*expuestos_e_sin_inf_suave[1:100]*suavizado_i_e[1:100]) + 
                sum(expuestos_e_sin_inf_suave[2:101]*susceptibles_e_sin_inf_suave[1:100]*suavizado_i_e[1:100]))/
                (2*sum((susceptibles_e_sin_inf_suave[1:100])^2*(suavizado_i_e[1:100])^2))

solution_44 <- ode(y= solution_2[1, 2:5], 1:100, seir, c(4, beta_e_sin_suave, 0.1, 0.3, 0.2))

estimacion_final_c_sin_recu <- data.frame(susceptibles_e_sin_inf_suave[1:100], expuestos_e_sin_inf_suave[1:100], 
                                          suavizado_i_e[1:100], recuperados_sim_1[1:100], 
                                          solution_44[1:100, ], solution_1[1:100, 2:5], solution_2[1:100, 4])

colors4 <- c("Estimated susceptible" = "royalblue1", "Solution 1" = "blue4", "Solution 2" = "deepskyblue")
colors5 <- c("Estimated exposed" = "orange2", "Solution 1" = "darkorange3", "Solution 2" = "gold2")
colors1 <- c("Infected data" = "indianred2",  "Solution 1" = "firebrick3", "Solution 2" = "magenta3")
colors2 <- c('Estimated recovered' = "chartreuse3", "Solution 1" = "springgreen3", "Solution 2" = "darkgreen")

p16_16 <- ggplot(data = estimacion_final_c_sin_recu) + geom_point(aes(x = estimacion_final_c_sin_recu[, 5], 
                                                                   y = estimacion_final_c_sin_recu[, 1]
                                                          , color = "Estimated susceptible"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[,6], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[,10], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colors4, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
                                                                                      linetype = c("blank", "solid", "solid"),
                                                                                      shape = c(16, NA, NA))))  +
  theme(legend.position="bottom") 

p17_17 <- ggplot(data = estimacion_final_c_sin_recu) + geom_point(aes(x = estimacion_final_c_sin_recu[, 5], 
                                                                   y = estimacion_final_c_sin_recu[, 2], 
                                                          color = "Estimated exposed"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[,7], color =  "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[,11], color =  "Solution 2"), size = 1.2, linetype = "dashed") +
  labs(x = 'Time (days)', y = 'No. of exposed population')  + scale_color_manual(values = colors5, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("blank", "solid", "solid"),
                                                                                   shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p18_18 <- ggplot(data = estimacion_final_c_sin_recu) + geom_point(aes(x = estimacion_final_c_sin_recu[, 5], 
                                                                      y = estimacion_final_c_sin_recu[, 14], 
                                                          color = "Infected data"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[, 8], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[,12], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of infected population')  + scale_color_manual(values = colors1, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p19_19 <- ggplot(data = estimacion_final_c_sin_recu) + geom_point(aes(x = estimacion_final_c_sin_recu[, 5], 
                                                                   y = estimacion_final_c_sin_recu[, 4], 
                                                          color = 'Estimated recovered'), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[, 9], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu[, 5], y = estimacion_final_c_sin_recu[,13], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + scale_color_manual(values = colors2, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")  

grid.arrange(p16_16, p17_17, p18_18, p19_19)

#########################
expuestos_e_s_final <- ((delta_inf_s[1:101] + recuperados_sim_1[2:101] - recuperados_sim_1[1:100] - delta_pob[1:101])/upsilon_1) + (eta_1/upsilon_1) - 
  ((mu_1/upsilon_1)*(pob[1:101] - suavizado_i_e[1:101] - recuperados_sim_1[1:101]))
plot(expuestos_e_s_final)
lines(solution_1[, 3])

susceptibles_e_sin_inf_suave_mas <- pob[1:101] - expuestos_e_s_final - suavizado_i_e[1:101] - recuperados_sim_1

beta_e_sin_suave_fin = ((1-mu_1)*sum(((susceptibles_e_sin_inf_suave_mas[1:100])^2*suavizado_i_e[1:100])) - 
                      sum(susceptibles_e_sin_inf_suave_mas[2:101]*susceptibles_e_sin_inf_suave_mas[1:100]*suavizado_i_e[1:100])  + 
                      eta_1*sum(susceptibles_e_sin_inf_suave_mas[1:100]*suavizado_i_e[1:100]) - (1 - upsilon_1 - mu_1)*
                      sum(susceptibles_e_sin_inf_suave_mas[1:100]*expuestos_e_s_final[1:100]*suavizado_i_e[1:100]) + 
                      sum(expuestos_e_s_final[2:101]*susceptibles_e_sin_inf_suave_mas[1:100]*suavizado_i_e[1:100]))/
  (2*sum((susceptibles_e_sin_inf_suave_mas[1:100])^2*(suavizado_i_e[1:100])^2))

solution_44_44 <- ode(y= solution_2[1, 2:5], 1:100, seir, c(4, beta_e_sin_suave_fin, 0.1, 0.3, 0.2))

estimacion_final_c_sin_recu_suave <- data.frame(susceptibles_e_sin_inf_suave_mas[1:100], expuestos_e_s_final[1:100], 
                                          suavizado_i_e[1:100], recuperados_sim_1[1:100], 
                                          solution_44_44[1:100, ], solution_1[1:100, 2:5], solution_2[1:100, 4])

colors4 <- c("Estimated susceptible" = "royalblue1", "Solution 1" = "blue4", "Solution 2" = "deepskyblue")
colors5 <- c("Estimated exposed" = "orange2", "Solution 1" = "darkorange3", "Solution 2" = "gold2")
colors1 <- c("Infected data" = "indianred2",  "Solution 1" = "firebrick3", "Solution 2" = "magenta3")
colors2 <- c('Estimated recovered' = "chartreuse3", "Solution 1" = "springgreen3", "Solution 2" = "darkgreen")

p16_16_16 <- ggplot(data = estimacion_final_c_sin_recu_suave) + geom_point(aes(x = estimacion_final_c_sin_recu_suave[, 5], 
                                                                      y = estimacion_final_c_sin_recu_suave[, 1]
                                                                      , color = "Estimated susceptible"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[,6], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[,10], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + scale_color_manual(values = colors4, name = "", 
                                                                                    guide = guide_legend(override.aes = list(
                                                                                      linetype = c("blank", "solid", "solid"),
                                                                                      shape = c(16, NA, NA))))  +
  theme(legend.position="bottom") 

p17_17_17 <- ggplot(data = estimacion_final_c_sin_recu_suave) + geom_point(aes(x = estimacion_final_c_sin_recu_suave[, 5], 
                                                                      y = estimacion_final_c_sin_recu_suave[, 2], 
                                                                      color = "Estimated exposed"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[,7], color =  "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[,11], color =  "Solution 2"), size = 1.2, linetype = "dashed") +
  labs(x = 'Time (days)', y = 'No. of exposed population')  + scale_color_manual(values = colors5, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("blank", "solid", "solid"),
                                                                                   shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p18_18_18 <- ggplot(data = estimacion_final_c_sin_recu_suave) + geom_point(aes(x = estimacion_final_c_sin_recu_suave[, 5], 
                                                                      y = estimacion_final_c_sin_recu_suave[, 14], 
                                                                      color = "Infected data"), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[, 8], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[,12], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of infected population')  + scale_color_manual(values = colors1, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")

p19_19_19 <- ggplot(data = estimacion_final_c_sin_recu_suave) + geom_point(aes(x = estimacion_final_c_sin_recu_suave[, 5], 
                                                                      y = estimacion_final_c_sin_recu_suave[, 4], 
                                                                      color = 'Estimated recovered'), size = 1.5) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[, 9], color = "Solution 1"), size = 1.2) + 
  geom_line(aes(x = estimacion_final_c_sin_recu_suave[, 5], y = estimacion_final_c_sin_recu_suave[,13], color = "Solution 2"), size = 1.2, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + scale_color_manual(values = colors2, name = "", 
                                                                                  guide = guide_legend(override.aes = list(
                                                                                    linetype = c("blank", "solid","solid"),
                                                                                    shape = c(16, NA, NA))))  +
  theme(legend.position="bottom")  

grid.arrange(p16_16_16, p17_17_17, p18_18_18, p19_19_19)


# Estimación de la volatilidad por MCO: 

z <- NULL 
for (i in 1:100) {
  z[i] <- rnorm(1)
}
sigma_s = (sum((solution_4[2:101, 2] - susceptibles_s[2:101] + expuestos_e[2:101] 
               - solution_4[2:101, 3])*susceptibles_s[1:100]*solution_2[1:100, 4]*z))/(2*
               sum((susceptibles_s[1:100]*solution_2[1:100, 4])^2*z^2))
sigma_s

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
  SN[[k]][1] = solution_1[1, 2]
  EN[[k]][1] = solution_1[1, 3]
  IN[[k]][1] = solution_1[1, 4]
  RN[[k]][1] = solution_1[1, 5]
}

for (k in 1:5) {
  for(j in 1:99){
    SN[[k]][j+1] = SN[[k]][j] + (eta - beta_s*SN[[k]][j]*IN[[k]][j] - mu*SN[[k]][j]) - sigma_s*SN[[k]][j]*IN[[k]][j]*z[[k]][j]
    EN[[k]][j+1] = EN[[k]][j] + (beta_s*SN[[k]][j]*IN[[k]][j] - upsilon*EN[[k]][j] - mu*EN[[k]][j]) + sigma_s*SN[[k]]*IN[[k]]*z[[k]][j]
    IN[[k]][j+1] = IN[[k]][j] + (upsilon*EN[[k]][j] - gamma_s*IN[[k]][j] - mu*IN[[k]][j])
    RN[[k]][j+1] = RN[[k]][j] + (gamma_s*IN[[k]][j] - mu*RN[[k]][j])
  }}

graf_suscp <- data.frame(solution_2[1:100, 1:2], SN[[1]], SN[[2]], SN[[3]], SN[[4]], SN[[5]]) 
graf_expue <- data.frame(solution_2[1:100, 1], solution_2[1:100, 3], EN[[1]], EN[[2]], EN[[3]], EN[[4]], EN[[5]])
graf_infec <- data.frame(solution_2[1:100, 1], solution_2[1:100, 4], IN[[1]], IN[[2]], IN[[3]], IN[[4]], IN[[5]])
graf_recup <- data.frame(solution_2[1:100, 1], solution_2[1:100, 5], RN[[1]], RN[[2]], RN[[3]], RN[[4]], RN[[5]])

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

# Estimación no paramétrica de gamma 

gama <- (delta_recuperados[1:385] + 0.004*suavizado_recuperados$X1[1:385])/(suavizado_infectados$X1[1:385])
mean(gama)
plot(gama, type = "l")

no_recuperados <- NULL
no_recuperados_1 <- NULL 
no_recuperados[1] <- suavizado_recuperados$X1[1]
no_recuperados_1[1] <- suavizado_recuperados$X1[1]
no_recuperados[2:386] <- suavizado_recuperados$X1[1:385] + (gama*suavizado_infectados$X1[1:385] - 0.004*suavizado_recuperados$X1[1:385])
no_recuperados_1[2:386] <- suavizado_recuperados$X1[1:385] + (mean(gama)*suavizado_infectados$X1[1:385] - 0.004*suavizado_recuperados$X1[1:385])

plot(no_recuperados, type = "l", col = "blue")
no_recuperados - no_recuperados_1
points(no_recuperados_1)
points(suavizado_recuperados$X1)

upsilon_no = (delta_infectados[1:385] + (0.004 + gama)*suavizado_infectados$X1[1:385])/(expuestos_bogota[1:385])
plot(upsilon_no)
hist(upsilon_no, breaks = 10)
mean(upsilon_no)

no_infectados <- NULL 
no_infectados[1] <- suavizado_infectados$X1[1]
no_infectados[2:386] <- suavizado_infectados$X1[1:385] + upsilon_no*expuestos_bogota[1:385] - (0.004 + gama)*suavizado_infectados$X1[1:385]
plot(no_infectados, type = "l", col = "red")
points(suavizado_infectados$X1)
sum(no_infectados - suavizado_infectados$X1[1:386])

beta_no <- (expuestos_bogota[2:386] - expuestos_bogota[1:385] + (upsilon_no + 0.004)*expuestos_bogota[1:385])/
           (susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385])
plot(beta_no, ylim = c(0, 0.00001))
hist(beta_no, breaks = 500, xlim = c(0, 4e-06))
mean(beta_no)

beta1_no <- (- susceptibles_bogota[2:386] + susceptibles_bogota[1:385] + 73660.35 - 0.004*susceptibles_bogota[1:385])/ 
            (susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385])
hist(beta1_no, breaks = 100)
mean(beta1_no)

no1_infectados <- NULL 
no1_infectados[1] <- suavizado_infectados$X1[1]
no1_infectados[2:386] <- suavizado_infectados$X1[1:385] + (1/5.2)*expuestos_bogota[1:385] - (0.004 + gama)*suavizado_infectados$X1[1:385]
plot(no1_infectados[2:386], type = "l", col = "red", ylim = c(0, max(no1_infectados[2:386])))
lines(suavizado_infectados$X1)
no1_infectados - suavizado_infectados$X1[1:386]

library(ggplot2)
no_sirve_estimacion <- data.frame(tiempo[2:385], no1_infectados[2:385], suavizado_infectados$X1[2:385], casos$Infectados[2:385])

colors <- c("Infected data" = "indianred2", "Smoothed infected" = "firebrick3", "Data update estimation" = "darkorange4")
pnosirve <- ggplot(data = no_sirve_estimacion) + geom_point(aes(x = no_sirve_estimacion[, 1], y = no_sirve_estimacion[, 4]
                                                               , color = "Infected data"), size = 1.5) + 
  geom_line(aes(x = no_sirve_estimacion[, 1], y = no_sirve_estimacion[, 3], color = "Smoothed infected"), size = 1.2) + 
  geom_line(aes(x = no_sirve_estimacion[, 1], y = no_sirve_estimacion[,2], color = "Data update estimation"), size = 1.2) + 
  labs(x = 'Time (days)', y = 'No. of infected population') +  scale_color_manual(values = colors, name = "", 
                                                                            guide = guide_legend(override.aes = list(
                                                                              linetype = c("blank", "solid", "solid"),
                                                                              shape = c(16, NA, NA)))) + 
  theme(legend.position="right") + scale_x_date(breaks = "2 months", date_labels = "%d-%m-%Y")

grid.arrange(pnosirve)
mean(upsilon_no)

expuestos_no <- NULL
expuestos_no[1] <- expuestos_bogota[1]
expuestos_no[2:386] <- expuestos_bogota[1:385] + beta_no*susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385]- (upsilon_no + 0.004)*expuestos_bogota[1:385]
plot(expuestos_no, type = "l", col = "red")
points(expuestos_bogota)
expuestos_no - expuestos_bogota

expuestos1_no <- NULL
expuestos1_no[1] <- expuestos_bogota[1]
expuestos1_no[2:386] <- expuestos_bogota[1:385] + beta1_no*susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385]- (upsilon + mu)*expuestos_bogota[1:385]
plot(expuestos1_no[1:385], type = "l", col = "red", ylim = c(min(expuestos_bogota), max(expuestos1_no[2:385])))
lines(expuestos_bogota)

susceptibles_no <- susceptibles_bogota[1:385] + 73660.35 - 
                   beta_no*susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385] - 0.004*susceptibles_bogota[1:385]
susceptibles_no1 <- susceptibles_bogota[1:385] + 73660.35 - 
  beta1_no*susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385] - 0.004*susceptibles_bogota[1:385]
plot(susceptibles_bogota, type = "p", ylim = c(7400000, 7800000))
lines(susceptibles_no, lwd = 2, col = "blue")
lines(susceptibles_no1, lwd = 2, col = "red")
susceptibles_no <- proyeccion_bogota_final[1:385] - expuestos_no[1:385] - suavizado_infectados$X1[1:385] - suavizado_recuperados$X1[1:385]

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

x02 <- c(susceptibles_bogota[1], expuestos_bogota[1], suavizado_infectados$X1[1], suavizado_recuperados$X1[1])
parametros2 = c(eta, mean(beta_no), mean(upsilon_no), mean(gama), mu) # Mala estimación 

library(deSolve)
plot(ode(y=x02, 1:385, seir, parametros2)[, 2], type = "l")
points(1:385, susceptibles_bogota[1:385], pch = 16)
plot(ode(y=x02, 1:385, seir, parametros2)[, 3], type = "l")
points(1:385, expuestos_bogota[1:385], pch = 16)
plot(ode(y=x02, 1:385, seir, parametros2)[, 4], type = "l")
points(1:385, suavizado_infectados$X1[1:385], pch = 16)
plot(ode(y=x02, 1:385, seir, parametros2)[, 5], type = "l")
points(1:385, suavizado_recuperados$X1[1:385], pch = 16)

susceptibles_bogota[1] + (9.42/1000)*mean(proyeccion_bogota_final) - 
  beta_no[1]*susceptibles_bogota[1]*suavizado_infectados$X1[1] - 
  0.004*susceptibles_bogota[1]
susceptibles_bogota[2] + (9.42/1000)*mean(proyeccion_bogota_final) - 
  beta_no[2]*susceptibles_bogota[2]*suavizado_infectados$X1[2] - 
  0.004*susceptibles_bogota[2]
plot(no_no_susceptibles[[1]])

plot(beta_no, col = "blue", lwd = 2, type = "l", ylim = c(min(beta_no), max(beta_no)))
lines(beta1_no, col = "red", lwd = 2)

# Si lo hacemos por intervalos 

gamma_no_no <- list()
upsilon_no_no <- list()
beta_no_no <- list()
no_parametrico <- matrix()
no_no_recuperados <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
no_no_infectados <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
no_no_expuestos <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
no_no_susceptibles <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
for (j in 1:8) {
  gamma_no_no[[j]] <- (delta_recuperados[intervalo[[j]]] + 0.004*suavizado_recuperados$X1[intervalo[[j]]])/(suavizado_infectados$X1[intervalo[[j]]])
  upsilon_no_no[[j]] <- (delta_infectados[intervalo[[j]]] + (0.004 + gamma_no_no[[j]])*suavizado_infectados$X1[intervalo[[j]]])/(expuestos_bogota[intervalo[[j]]])
  beta_no_no[[j]] <- (expuestos_bogota[intervalo[[j]]+1] - expuestos_bogota[intervalo[[j]]] + 
                     (upsilon_no_no[[j]] + 0.004)*expuestos_bogota[intervalo[[j]]])/
                     (susceptibles_bogota[intervalo[[j]]]*suavizado_infectados$X1[intervalo[[j]]])
  no_no_recuperados[[1]][1] <- suavizado_recuperados$X1[intervalo[[1]]][1]
  no_no_infectados[[1]][1] <- suavizado_infectados$X1[intervalo[[1]]][1]
  no_no_expuestos[[1]][1] <- expuestos_bogota[intervalo[[1]]][1]
  no_no_susceptibles[[1]][1] <- susceptibles_bogota[intervalo[[1]]][1]
  no_no_recuperados[[j]][2:(length(intervalo[[j]])+1)] <- suavizado_recuperados$X1[intervalo[[j]]] + 
                                                (gamma_no_no[[j]]*suavizado_infectados$X1[intervalo[[j]]] - 
                                                   0.004*suavizado_recuperados$X1[intervalo[[j]]])
  no_no_infectados[[j]][2:(length(intervalo[[j]])+1)] <- suavizado_infectados$X1[intervalo[[j]]] + 
                                               upsilon_no_no[[j]]*expuestos_bogota[intervalo[[j]]] - 
                                               (0.004 + gamma_no_no[[j]])*suavizado_infectados$X1[intervalo[[j]]]
  no_no_expuestos[[j]][2:(length(intervalo[[j]])+1)] <- expuestos_bogota[intervalo[[j]]] + 
                                              beta_no_no[[j]]*susceptibles_bogota[intervalo[[j]]]*
                                              suavizado_infectados$X1[intervalo[[j]]]- 
                                              (upsilon_no_no[[j]] + 0.004)*expuestos_bogota[intervalo[[j]]]
  no_no_susceptibles[[j]][2:(length(intervalo[[j]])+1)] <- susceptibles_bogota[intervalo[[j]]] + (9.42/1000)*mean(proyeccion_bogota_final) - 
                                                 beta_no_no[[j]]*susceptibles_bogota[intervalo[[j]]]*
                                                 suavizado_infectados$X1[intervalo[[j]]] - 
                                                 0.004*susceptibles_bogota[intervalo[[j]]]
}
susceptibles_bogota
final_s <- NULL 
final_s <- na.omit(c(no_no_susceptibles[[1]][1:length(intervalo[[1]])], 
             no_no_susceptibles[[2]][1:length(intervalo[[2]])], 
             no_no_susceptibles[[3]][1:length(intervalo[[3]])], 
             no_no_susceptibles[[4]][1:length(intervalo[[4]])], 
             no_no_susceptibles[[5]][1:length(intervalo[[5]])], 
             no_no_susceptibles[[6]][1:length(intervalo[[6]])], 
             no_no_susceptibles[[7]][1:length(intervalo[[7]])], 
             no_no_susceptibles[[8]][1:length(intervalo[[8]])]))
par(mfrow=c(1,1))
plot(2:385, final_s[2:385], ylim = c(6000000, max(final_s)), type = "l")
lines(1:385, susceptibles_bogota[1:385])

final_e <- NULL 
final_e <- na.omit(c(no_no_expuestos[[1]][1:length(intervalo[[1]])], 
                     no_no_expuestos[[2]][1:length(intervalo[[2]])], 
                     no_no_expuestos[[3]][1:length(intervalo[[3]])], 
                     no_no_expuestos[[4]][1:length(intervalo[[4]])], 
                     no_no_expuestos[[5]][1:length(intervalo[[5]])], 
                     no_no_expuestos[[6]][1:length(intervalo[[6]])], 
                     no_no_expuestos[[7]][1:length(intervalo[[7]])], 
                     no_no_expuestos[[8]][1:length(intervalo[[8]])]))
plot(2:385, final_e[2:385], type = "l", col = "red", ylim = c(min(expuestos_bogota), max(final_e)))
points(1:385, expuestos_bogota[1:385])

final_i <- NULL 
final_i <- na.omit(c(no_no_infectados[[1]][1:length(intervalo[[1]])], 
                     no_no_infectados[[2]][1:length(intervalo[[2]])], 
                     no_no_infectados[[3]][1:length(intervalo[[3]])], 
                     no_no_infectados[[4]][1:length(intervalo[[4]])], 
                     no_no_infectados[[5]][1:length(intervalo[[5]])], 
                     no_no_infectados[[6]][1:length(intervalo[[6]])], 
                     no_no_infectados[[7]][1:length(intervalo[[7]])], 
                     no_no_infectados[[8]][1:length(intervalo[[8]])]))
plot(2:385, final_i[2:385], type = "l", ylim = c(min(suavizado_infectados$X1), max(final_i)), col = "red")
points(1:385, suavizado_infectados$X1[1:385])

final_r <- NULL 
final_r <- na.omit(c(no_no_recuperados[[1]][1:length(intervalo[[1]])], 
                     no_no_recuperados[[2]][1:length(intervalo[[2]])], 
                     no_no_recuperados[[3]][1:length(intervalo[[3]])], 
                     no_no_recuperados[[4]][1:length(intervalo[[4]])], 
                     no_no_recuperados[[5]][1:length(intervalo[[5]])], 
                     no_no_recuperados[[6]][1:length(intervalo[[6]])], 
                     no_no_recuperados[[7]][1:length(intervalo[[7]])], 
                     no_no_recuperados[[8]][1:length(intervalo[[8]])]))
plot(2:385, final_r[2:385], type = "l", ylim = c(min(suavizado_recuperados$X1), max(final_r)), col = "red")
points(1:385, suavizado_recuperados$X1[1:385])

colors4 <- c("Estimated\n susceptible" = "royalblue1", "Data update\n susceptible estimation" = "blue4")
colors5 <- c("Estimated\n exposed" = 'tan1', "Data update\n exposed estimation" = "darkorange3")
colors1 <- c("Smoothed\n infected" = "firebrick3", "Data update\n infected estimation" 
             = "magenta3", "Infected\n data" = "indianred2")
colors2 <- c("Smoothed\n recovered" = "springgreen3", "Data update\n recovered estimation" = "darkgreen", 
             "Recovered\n data" = "chartreuse3")
tiempo <- seq(from = as.Date("2020-03-06"), to = as.Date("2021-03-29"), by=1)
grafica_no_parametrica <- data.frame(tiempo[1:385], final_s[1:385], final_e[1:385], final_i[1:385], final_r[1:385], susceptibles_bogota[1:385], 
                                     expuestos_bogota[1:385], suavizado_infectados$X1[1:385], suavizado_recuperados$X1[1:385], 
                                     casos$Infectados[1:385], casos$recuperados[1:385])

p_ns <- ggplot(data = grafica_no_parametrica) + geom_line(aes(x = grafica_no_parametrica[, 1], y = grafica_no_parametrica[, 6],
                                                               color = "Estimated\n susceptible"), size = 1.2) + 
  geom_point(aes(x = grafica_no_parametrica[, 1], y = grafica_no_parametrica[,2], color = "Data update\n susceptible estimation"),
            size = 0.8) + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + # ylim(0, max(grafica_no_parametrica[, 2])) + 
  scale_color_manual(values = colors4, name = "", guide = guide_legend(override.aes = list(linetype = 
                      c("solid", "solid"), shape = c(NA, NA)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p_ne <- ggplot(data = grafica_no_parametrica) + geom_point(aes(x = grafica_no_parametrica[, 1], y = grafica_no_parametrica[, 7], 
                                                      color = "Estimated\n exposed"), size = 1.5) + 
  geom_line(aes(x =grafica_no_parametrica[, 1], y = grafica_no_parametrica[,3], color = "Data update\n exposed estimation", 
                ), size = 1.1, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of exposed population') + 
  scale_color_manual(values = colors5, name = "", guide = guide_legend(override.aes = list(linetype = 
                     c("blank", "solid"), shape = c(16, NA)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p_ni <- ggplot(data = grafica_no_parametrica) +   geom_point(aes(x = grafica_no_parametrica[, 1], 
        y = grafica_no_parametrica[,10], color =  "Infected\n data"), size = 1.5) +
  geom_line(aes(x = grafica_no_parametrica[, 1], y = grafica_no_parametrica[, 8], color = "Smoothed\n infected"), 
            size = 1.5) + 
  geom_line(aes(x = grafica_no_parametrica[, 1], y = grafica_no_parametrica[,4], 
                color =  "Data update\n infected estimation"), size = 1.5, linetype = "4C88C488") + 
  labs(x = 'Time (days)', y = 'No. of infected population') +
  scale_color_manual(values = colors1, name = "", guide = guide_legend(override.aes = list(linetype = 
                     c("solid", "solid", "blank"), shape = c(NA, NA, 16)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p_nr <- ggplot(data = grafica_no_parametrica) + geom_point(aes(x = grafica_no_parametrica[, 1], 
  y = grafica_no_parametrica[, 11], color = "Recovered\n data"), size = 1.5) + 
  geom_point(aes(x = grafica_no_parametrica[, 1], y = grafica_no_parametrica[, 9], 
             color = "Smoothed\n recovered"), size = 1.5) + 
  geom_line(aes(x = grafica_no_parametrica[, 1], y = grafica_no_parametrica[,5], 
                color = "Data update\n recovered estimation"), size = 1.5, linetype = "4C88C488") + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + 
  scale_color_manual(values = colors2, name = "", guide = guide_legend(override.aes = list(linetype = 
                     c("solid", "solid", "blank"), shape = c(NA, NA, 16)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")
library(ggplot2)
library(gridExtra)
grid.arrange(p_ns, p_ne, p_ni, p_nr) 

grafica_susceptibles_mejorado <- data.frame(tiempo[1:385], proyeccion_bogota_final[1:385] - final_e[1:385] - final_i[1:385] - final_r[1:385], 
                                            susceptibles_bogota[1:385])

p_ns_m <- ggplot(data = grafica_susceptibles_mejorado) + 
  geom_point(aes(x = grafica_susceptibles_mejorado[, 1], y = grafica_susceptibles_mejorado[, 3],
                 color = "Estimated\n susceptible"), size = 1.5) + 
  geom_line(aes(x = grafica_susceptibles_mejorado[, 1], y = grafica_susceptibles_mejorado[,2], 
                 color = "Data update\n susceptible estimation"), size = 1.5, linetype = "dashed") + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + 
  scale_color_manual(values = colors4, name = "", guide = guide_legend(override.aes = list(linetype = 
                     c("dashed", "solid"), shape = c(NA, NA)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")
plot(p_ns_m)

mean(susceptibles_bogota[1:385] - final_s)

library(gridExtra)
grid.arrange(p_i, p_r, ncol = 2)

matriz_1 <- matrix(0, nrow = length(intervalo), ncol = 6)

for (i in 1:length(intervalo)) {
  matriz_1[i, 1] <- mean(beta_no[intervalo[[i]]])
  matriz_1[i, 5]<- mean(gama[intervalo[[i]]])
  matriz_1[i, 3] <- mean(upsilon_no[intervalo[[i]]])
}

for (i in 1:length(intervalo)) {
  matriz_1[i, 2] <- sd(beta_no[intervalo[[i]]])
  matriz_1[i, 6] <- sd(gama[intervalo[[i]]])
  matriz_1[i, 4] <- sd(upsilon_no[intervalo[[i]]])
}

matriz_1

library(deSolve)
x0_np <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)
solutions_np <- list()
for (i in 1:length(intervalo)) {
  library(deSolve)
  x0_np[[i]] <- c(susceptibles_b[[i]][1], expuestos_b[[i]][1], infectados_b[[i]][1], recuperados_b[[i]][1])
  solutions_np[[i]] <- ode(y = x0_np[[i]], intervalo[[i]], seir, c(eta, matriz_1[i, 1], matriz_1[i, 3], 
                                                                   matriz_1[i, 5], mu))
}

soluciones_np <- rbind(solutions_np[[1]], solutions_np[[2]], solutions_np[[3]], solutions_np[[4]], solutions_np[[5]], 
                    solutions_np[[6]], solutions_np[[7]], solutions_np[[8]])
plot(soluciones_np[, 4], col = "blue", type = "p", ylim = c(0, max(suavizado_infectados$X1)), pch = 16)
lines(suavizado_infectados$X1)

plot(soluciones_np[, 5], col = "blue", type = "l", ylim = c(0, max(suavizado_recuperados$X1)), pch = 16)
lines(suavizado_recuperados$X1)

plot(soluciones_np[, 2], col = "blue", type = "l", ylim = c(0, max(soluciones_np[, 2])), pch = 16)
lines(susceptibles_bogota)

### MÁS DEL ENFOQUE NO PARAMÉTRICO

upsilon_espuestos = (delta_infectados[1:385] + (mu + gama)*suavizado_infectados$X1[1:385])
plot(upsilon_espuestos)
lambda = upsilon_espuestos[1:385] + mu*s_e[1:385] - (s_e[2:386] - s_e[1:385])
plot(upsilon_espuestos + mu*s_e - (s_e[2:386] - s_e[1:385]))

# No tenemos eta: 

s_e = proyeccion_bogota_final[1:386] - suavizado_recuperados$X1[1:386] -  suavizado_infectados$X1[1:386]
plot(s_e)

expuestos = (eta - mu*s_e[2:386] - (s_e[2:386] - s_e[1:385]))/upsilon
plot(expuestos)
lines(expuestos_bogota)

expuestos1 = delta_infectados[1:385] + (mu + gama[1:385])*suavizado_infectados$X1[1:385]
plot(expuestos1)

tiempo = 1:385
library(npregfast)
beta_no_s <- frfast(beta_no ~ tiempo, model = "np", smooth = "kernel", kbin = 385, p =3)
beta_no_s <- data.frame(beta_no_s$p)
gamma_no_s <- frfast(gama ~ tiempo, model = "np", smooth = "kernel", kbin = 385, p =3)
gamma_no_s <- data.frame(gamma_no_s$p)
upsilon_no_s <- frfast(upsilon_no ~ tiempo, model = "np", smooth = "kernel", kbin = 385, p =3)
upsilon_no_s <- data.frame(upsilon_no_s$p)

grafica_parametros <- data.frame(tiempo[1:385], beta_no, upsilon_no, gama, beta_no_s$X1, upsilon_no_s$X1, gamma_no_s$X1)

colors_para1 <- c("Data update\n estimated beta" = 'mediumorchid1', "Smoothed data\n update estimated\n for beta" = 'mediumorchid4')  
colors_para2 <- c("Data update\n estimated upsilon" = "goldenrod1", "Smoothed data\n update estimated\n for upsilon" = 'goldenrod4')
colors_para3 <- c("Data update\n estimated gamma" = "cadetblue2", "Smoothed data\n update estimated\n for gamma" = 'cyan4')
p_beta <- ggplot(data = grafica_parametros) + geom_point(aes(x = grafica_parametros[, 1], y = grafica_parametros[, 2], 
                                                             color = "Data update\n estimated beta"), size = 2)  +
  geom_line(aes(x = grafica_parametros[, 1], y = grafica_parametros[, 5], color = "Smoothed data\n update estimated\n for beta"), 
            size = 1.2, linetype = "twodash") + 
  labs(x = 'Time (days)', y = expression(hat(beta)[j])) + 
  scale_color_manual(values = colors_para1 , name = "", guide = guide_legend(override.aes = list(linetype = 
       c("blank", "twodash"), shape = c(16, NA)))) + 
  theme(legend.position="right") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")
p_upsilon <- ggplot(data = grafica_parametros) +  geom_point(aes(x = grafica_parametros[, 1], y = grafica_parametros[, 3], 
                                                                 color = "Data update\n estimated upsilon"), size = 2)  +
  geom_line(aes(x = grafica_parametros[, 1], y = grafica_parametros[, 6], color = "Smoothed data\n update estimated\n for upsilon"), 
            size = 1.2, linetype = "twodash") +
  labs(x = 'Time (days)', y = expression(hat(upsilon)[j])) + 
  scale_color_manual(values = colors_para2 , name = "", guide = guide_legend(override.aes = list(linetype = 
   c("blank", "twodash"), shape = c(16, NA)))) + 
  theme(legend.position = "right") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")
p_gamma <- ggplot(data = grafica_parametros)  +  geom_point(aes(x = grafica_parametros[, 1], y = grafica_parametros[, 4], 
                                                                color = "Data update\n estimated gamma"), size = 2)  +
  geom_line(aes(x = grafica_parametros[, 1], y = grafica_parametros[, 7], color = "Smoothed data\n update estimated\n for gamma"), 
            size = 1.2, linetype = "twodash") +
  labs(x = 'Time (days)', y = expression(hat(gamma)[j])) + 
  scale_color_manual(values = colors_para3 , name = "", guide = guide_legend(override.aes = list(linetype = 
  c("blank", "twodash"), shape = c(16, NA)))) + 
  theme(legend.position = "right") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

grid.arrange(p_beta, p_upsilon, p_gamma, ncol = 1)

mu = 4/1000 # https://saludata.saludcapital.gov.co/osb/index.php/datos-de-salud/demografia/tm-bruta/
# 4 muertes por cada 1000 habitantes 
# mean(759/proyeccion_bogota_final) https://datosmacro.expansion.com/demografia/mortalidad/colombia 759 muertos dia
eta = (9.42/1000)*mean(proyeccion_bogota_final)  # max(as.numeric(delta_poblacion[1:420])) # https://datosmacro.expansion.com/
# 9.42 nacimientos por cada 1000 habitantes https://saludata.saludcapital.gov.co/osb/index.php/datos-de-salud/demografia/natalidad/

r0 = (eta*beta_no*upsilon_no)/(mu*(upsilon_no + mu)*(gama + mu))

tiempo = 1:385
library(npregfast)
r0_s <- frfast(r0 ~ tiempo, model = "np", smooth = "kernel", kbin = 385, p =3)
r0_s <- data.frame(r0_s$p)

grafica_r0 <- data.frame(tiempo[1:385], r0, r0_s$X1)

color_r0 <- c("Data update\n basic reproduction\n number" = "navajowhite3", "Smoothed data\n update estimated\n for the basic\n reproduction number" = 'navajowhite4')
library(ggplot2)
p_r0 <- ggplot(data = grafica_r0) + geom_point(aes(x = grafica_r0[, 1], y = grafica_r0[, 2], 
                                              color = "Data update\n basic reproduction\n number"), size = 2)  + 
  geom_line(aes(x = grafica_r0[, 1], y = grafica_r0[, 3], 
             color = "Smoothed data\n update estimated\n for the basic\n reproduction number"), size = 1.2) + 
  labs(x = 'Time (days)', y = expression(hat(R)[0[j]]))  +  ylim(quantile(r0, 0.01), quantile(r0, 0.99)) +
  scale_color_manual(values = color_r0 , name = "", guide = guide_legend(override.aes = list(linetype = 
                     c("blank", "solid"), shape = c(16, NA)))) +
  theme(legend.position = "right") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")
plot(p_r0)

matriz_2 <- matrix(0, nrow = length(intervalo), ncol = 4)

for (i in 1:length(intervalo)) {
  matriz_2[i, 1] <- mean(r0[intervalo[[i]]])
  matriz_2[i, 2] <- sd(r0[intervalo[[i]]])
  matriz_2[i, 3] <- mean(r0_s$X1[intervalo[[i]]])
  matriz_2[i, 4] <- sd(r0_s$X1[intervalo[[i]]])
}

matriz_2
matriz_1
