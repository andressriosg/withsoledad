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

tiempo = 0:(length(solution_2[, 4])-1)
library(npregfast)
suavizado_i_e <- frfast(solution_2[, 4] ~ tiempo, model = "np", smooth = "kernel", kbin = 102, p =3) # Suavizamos infectados 
suavizado_i_e <- data.frame(suavizado_i_e$p)
suavizado_i_e = suavizado_i_e$X1

suavizado_r_e <- frfast(solution_2[, 5] ~ tiempo, model = "np", smooth = "kernel", kbin = 102, p =3) # Suavizamos recuperados 
suavizado_r_e <- data.frame(suavizado_r_e$p)
suavizado_r_e = suavizado_r_e$X1

mu = 0.2
upsilon = 0.1
eta = 4

delta_inf_s <- NULL
for (j in 1:length(suavizado_i_e)) {
  delta_inf_s[j] = suavizado_i_e[j+1] - suavizado_i_e[j]
}

delta_rec_s <- NULL
for (j in 1:length(suavizado_r_e)) {
  delta_rec_s[j] = suavizado_r_e[j+1] - suavizado_r_e[j]
}

expuestos_e_s <- ((delta_inf_s[1:101] + delta_rec_s[1:101] - delta_pob[1:101])/upsilon) + (eta/upsilon) - 
  ((mu/upsilon)*(pob[1:101] - suavizado_i_e[1:101] - suavizado_r_e[1:101]))

susceptibles_s_s <- pob[1:101] - suavizado_i_e[1:101] - suavizado_r_e[1:101] - expuestos_e_s 

plot(susceptibles_s_s)
plot(suavizado_i_e)
lines(solution_2[, 4])

mu = 0.2
eta = 4
upsilon = 0.1
beta_s = ((1-mu)*sum(((susceptibles_s_s[1:100])^2*suavizado_i_e[1:100])) - sum(susceptibles_s_s[2:101]*susceptibles_s_s[1:100]*suavizado_i_e[1:100])  + eta*sum(susceptibles_s_s[1:100]*suavizado_i_e[1:100]) - (1 - upsilon - mu)*sum(susceptibles_s_s[1:100]*expuestos_e_s[1:100]*suavizado_i_e[1:100]) + sum(expuestos_e_s[2:101]*susceptibles_s_s[1:100]*suavizado_i_e[1:100]))/(2*sum((susceptibles_s_s[1:100])^2*(suavizado_i_e[1:100])^2))
beta_s
gamma_s = ((1-mu)*(sum((suavizado_i_e[1:100])^2) - sum(suavizado_r_e[1:100]*suavizado_i_e[1:100])) - sum(suavizado_i_e[1:100]*suavizado_i_e[2:101]) + sum(suavizado_i_e[1:100]*suavizado_r_e[1:100]) + upsilon*sum(expuestos_e_s[1:100]*suavizado_i_e[1:100]))/(2*sum((suavizado_i_e[1:100])^2))
gamma_s

library(deSolve)
solution_fin <- ode(y= c(4, 0, 0.1, 0), 1:100, seir, c(4, beta_s, 0.1, gamma_s, 0.2))
plot(solution_fin[, 3])
lines(expuestos_e_s)

z = NULL 
for (l in 1:99) {
  z[l] = mean(rnorm(10000, 0, 1))
}

sigma_3 = (1/2)*(sum(((solution_fin[7:100, 2]-susceptibles_s_s[7:100] + expuestos_e_s[7:100]-solution_fin[7:100, 3])*susceptibles_s_s[6:99]*suavizado_i_e[6:99]*z[6:99])))/
  (sum((susceptibles_s_s[6:99]*suavizado_i_e[6:99]*z[6:99])^2))
sigma_3 # 0.07527228

sigma_3 = (1/2)*(sum(((solution_fin[7:100, 2]-susceptibles_s_s[7:100] + expuestos_e_s[7:100]-solution_fin[7:100, 3])*solution_fin[6:99, 2]*solution_fin[6:99, 4]*z[6:99])))/
  (sum((solution_fin[6:99, 2]*solution_fin[6:99, 4]*z[6:99])^2))
sigma_3

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
  z[[k]] = as.vector(rnorm(1000, 0, 1))
}

for (k in 1:5) {
  SN[[k]][1] = solution_fin[1, 2]
  EN[[k]][1] = solution_fin[1, 3]
  IN[[k]][1] = solution_fin[1, 4]
  RN[[k]][1] = solution_fin[1, 5]
}

for (k in 1:5) {
  for(j in 1:99){
    SN[[k]][j+1] = SN[[k]][j] + (eta - beta_s*SN[[k]][j]*IN[[k]][j] - mu*SN[[k]][j]) - sigma_3*SN[[k]][j]*IN[[k]][j]*z[[k]][j]
    EN[[k]][j+1] = EN[[k]][j] + (beta_s*SN[[k]][j]*IN[[k]][j] - upsilon*EN[[k]][j] - mu*EN[[k]][j]) + sigma_3*SN[[k]]*IN[[k]]*z[[k]][j]
    IN[[k]][j+1] = IN[[k]][j] + (upsilon*EN[[k]][j] - gamma_s*IN[[k]][j] - mu*IN[[k]][j])
    RN[[k]][j+1] = RN[[k]][j] + (gamma_s*IN[[k]][j] - mu*RN[[k]][j])
  }}

graf_suscp <- data.frame(solution_2[1:100, 1], solution_fin[1:100, 2], susceptibles_s_s[1:100], 
                         SN[[1]], SN[[2]], SN[[3]], SN[[4]], SN[[5]]) 
graf_expue <- data.frame(solution_2[1:100, 1], solution_fin[1:100, 3], expuestos_e_s[1:100], 
                         EN[[1]], EN[[2]], EN[[3]], EN[[4]], EN[[5]])
graf_infec <- data.frame(solution_2[1:100, 1], solution_fin[1:100, 4], solution_2[1:100, 4], IN[[1]], 
                         IN[[2]], IN[[3]], IN[[4]], IN[[5]])
graf_recup <- data.frame(solution_2[1:100, 1], solution_fin[1:100, 5], solution_2[1:100, 5], RN[[1]], 
                         RN[[2]], RN[[3]], RN[[4]], RN[[5]])

library(ggplot2)

colors7 <- c("Estimated\n susceptible" = 'darkblue', "Deterministic\n solution" = "deepskyblue3", 
             "Stochastic\n paths" = 'cadetblue2')
colors8 <- c("Estimated\n exposed" = 'darkorange4', "Deterministic\n solution" = 'chocolate1', 
             "Stochastic\n paths" = 'lightsalmon')
colors9 <- c("Infected\n data" = 'red4', "Deterministic\n solution" = 'firebrick1', 
             "Stochastic\n paths" = 'palevioletred1')
colors10 <- c("Recovered\n data" = 'springgreen4', "Deterministic\n solution" = 'green3', 
              "Stochastic\n paths" = 'lightgreen')
p1 <- ggplot(data = graf_suscp)  + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 8], color = "Stochastic\n paths"), size = 0.8) +  
  geom_point(aes(x = graf_suscp[, 1], y = graf_suscp[, 3], color = "Estimated\n susceptible"), size = 1.7) + 
  geom_line(aes(x = graf_suscp[, 1], y = graf_suscp[, 2], color = "Deterministic\n solution"), size = 1.2) + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + 
  scale_color_manual(values = colors7, name = "", guide = guide_legend(override.aes = list(linetype = 
  c("blank", "solid", "solid"), shape = c(16, NA, NA)))) + 
  theme(legend.position="bottom")

p2 <- ggplot(data = graf_expue) + labs(x = 'Time (days)', y = 'No. of exposed population') + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_point(aes(x = graf_expue[, 1], y = graf_expue[, 3], color = "Estimated\n exposed"), size = 1.7) +
  geom_line(aes(x = graf_expue[, 1], y = graf_expue[, 2], color = "Deterministic\n solution"), size = 1.2) + 
  labs(x = 'Time (days)', y = 'No. of exposed population') + 
  scale_color_manual(values = colors8, name = "", guide = guide_legend(override.aes = list(linetype = 
  c("blank", "solid", "solid"), shape = c(16, NA, NA)))) + 
  theme(legend.position="bottom")

p3 <- ggplot(data = graf_infec) + labs(x = 'Time (days)', y = 'No. of infected population') + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_point(aes(x = graf_infec[, 1], y = graf_infec[, 3], color = "Infected\n data"), size = 1.7) + 
  geom_line(aes(x = graf_infec[, 1], y = graf_infec[, 2], color = "Deterministic\n solution"), size = 1.2)+ 
  labs(x = 'Time (days)', y = 'No. of infected population') + 
  scale_color_manual(values = colors9, name = "", guide = guide_legend(override.aes = list(linetype = 
  c("blank", "solid", "solid"), shape = c(16, NA, NA)))) + 
  theme(legend.position="bottom")

p4 <- ggplot(data = graf_recup) + labs(x = 'Time (days)', y = 'No. of recovered population') + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_point(aes(x = graf_recup[, 1], y = graf_recup[, 3], color = "Recovered\n data"), size = 1.7) +
  geom_line(aes(x = graf_recup[, 1], y = graf_recup[, 2], color = "Deterministic\n solution"), size = 1.2) + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + 
  scale_color_manual(values = colors10, name = "", guide = guide_legend(override.aes = list(linetype = 
  c("blank", "solid", "solid"), shape = c(16, NA, NA)))) + 
  theme(legend.position="bottom")

# install.packages("gridExtra")
library(gridExtra)
grid.arrange(p1, p2, p3, p4)

w = NULL 
for (l in 1:385) {
  w[l] = mean(rnorm(10000, 0, 1))
}

sigma_4 = (1/2)*(sum(((final_s[2:385] - susceptibles_bogota[2:385] + expuestos_bogota[2:385]- final_e[2:385])*susceptibles_bogota[1:384]*suavizado_infectados$X1[1:384]*w[1:384])))/
  (sum((susceptibles_bogota[1:384]*suavizado_infectados$X1[1:384]*w[1:384])^2))
sigma_4

sigma_5 = (1/2)*(sum(((proyeccion_bogota_final[2:385] - final_e[2:385] - final_i[2:385] - final_r[2:385] 
                       - susceptibles_bogota[2:385] + expuestos_bogota[2:385]- final_e[2:385])*susceptibles_bogota[1:384]*suavizado_infectados$X1[1:384]*w[1:384])))/
  (sum((susceptibles_bogota[1:384]*suavizado_infectados$X1[1:384]*w[1:384])^2))
sigma_5

a = NULL 
for (l in 1:385) {
  a[l] = mean(rnorm(10000, 0, 1))
}

SN1 = list() 
EN1 = list()
IN1 = list()
RN1 = list()
z1 = list()

for (k in 1:5) {
  SN1[[k]] = rep(0, 386)
  EN1[[k]] = rep(0, 386)
  IN1[[k]] = rep(0, 386)
  RN1[[k]] = rep(0, 386)
  z1[[k]] = as.vector(rnorm(386, 0, 1))
}

for (k in 1:5) {
  SN1[[k]][1] = susceptibles_bogota[1]
  EN1[[k]][1] = expuestos_bogota[1]
  IN1[[k]][1] = casos$Infectados[1]
  RN1[[k]][1] = casos$recuperados[1]
}


for (k in 1:5) {
#  SN1[[k]][2:386] = final_s[1:385] - sigma_4*susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385]*z1[[k]][1:385]
  SN1[[k]][2:386] = grafica_susceptibles_mejorado[1:385, 2] - sigma_4*susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385]*z1[[k]][1:385]
  EN1[[k]][2:386] = final_e[1:385] + sigma_4*susceptibles_bogota[1:385]*suavizado_infectados$X1[1:385]*z1[[k]][1:385]
  IN1[[k]][2:386] = final_i[1:385]  
  RN1[[k]][2:386] = final_r[1:385] 
  }


graf_suscp1 <- data.frame(tiempo[1:385], grafica_susceptibles_mejorado[1:385, 2], susceptibles_bogota[1:385], 
                         SN1[[1]][1:385], SN1[[2]][1:385], SN1[[3]][1:385], SN1[[4]][1:385], SN1[[5]][1:385]) 
graf_expue1 <- data.frame(tiempo[1:385], final_e[1:385], expuestos_bogota[1:385], 
                         EN1[[1]][1:385], EN1[[2]][1:385], EN1[[3]][1:385], EN1[[4]][1:385], EN1[[5]][1:385])
graf_infec1 <- data.frame(tiempo[1:385], final_i[1:385], suavizado_infectados$X1[1:385], IN1[[1]][1:385], 
                         IN1[[2]][1:385], IN1[[3]][1:385], IN1[[4]][1:385], IN1[[5]][1:385])
graf_recup1 <- data.frame(tiempo[1:385], final_r[1:385], suavizado_recuperados$X1[1:385], RN1[[1]][1:385],
                         RN1[[2]][1:385], RN1[[3]][1:385], RN1[[4]][1:385], RN1[[5]][1:385])

colors11 <- c("Estimated\n susceptible" = 'darkblue', "Data update\n solution" = "deepskyblue3", 
             "Stochastic\n paths" = 'cadetblue2')
colors21 <- c("Estimated\n exposed" = 'darkorange4', "Data update\n solution" = 'chocolate1', 
             "Stochastic\n paths" = 'lightsalmon')
colors31 <- c("Infected\n data" = 'red4', "Deterministic\n solution" = 'firebrick1', 
             "Stochastic\n paths" = 'palevioletred1')
colors41 <- c("Recovered\n data" = 'springgreen4', "Deterministic\n solution" = 'green3', 
              "Stochastic\n paths" = 'lightgreen')
p111 <- ggplot(data = graf_suscp1) + 
  geom_line(aes(x = graf_suscp1[, 1], y = graf_suscp1[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp1[, 1], y = graf_suscp1[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp1[, 1], y = graf_suscp1[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp1[, 1], y = graf_suscp1[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_suscp1[, 1], y = graf_suscp1[, 8], color = "Stochastic\n paths"), size = 0.8) +  
  geom_point(aes(x = graf_suscp1[, 1], y = graf_suscp1[, 3], color = "Estimated\n susceptible"), size = 2)  + 
  geom_line(aes(x = graf_suscp1[, 1], y = graf_suscp1[, 2], color = "Data update\n solution"), size = 1) + 
  labs(x = 'Time (days)', y = 'No. of susceptible population') + 
  scale_color_manual(values = colors11, name = "", guide = guide_legend(override.aes = list(linetype = 
  c("blank", "solid", "solid"), shape = c(16, NA, NA)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p21 <- ggplot(data = graf_expue1) + 
  geom_line(aes(x = graf_expue1[, 1], y = graf_expue1[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue1[, 1], y = graf_expue1[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue1[, 1], y = graf_expue1[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue1[, 1], y = graf_expue1[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_expue1[, 1], y = graf_expue1[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_point(aes(x = graf_expue1[, 1], y = graf_expue1[, 3], color = "Estimated\n exposed"), size = 2) +
  geom_line(aes(x = graf_expue1[, 1], y = graf_expue1[, 2], color = "Data update\n solution"), size = 1) + 
  labs(x = 'Time (days)', y = 'No. of exposed population') + 
  scale_color_manual(values = colors21, name = "", guide = guide_legend(override.aes = list(linetype = 
  c("blank", "solid", "solid"), shape = c(16, NA, NA)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p31 <- ggplot(data = graf_infec1) +  
  geom_line(aes(x = graf_infec1[, 1], y = graf_infec1[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec1[, 1], y = graf_infec1[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec1[, 1], y = graf_infec1[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec1[, 1], y = graf_infec1[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_infec1[, 1], y = graf_infec1[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_point(aes(x = graf_infec1[, 1], y = graf_infec1[, 3], color = "Infected\n data"), size = 1.7) + 
  geom_line(aes(x = graf_infec1[, 1], y = graf_infec1[, 2], color = "Deterministic\n solution"), size = 1.2) + 
  scale_color_manual(values = colors31, name = "",
                     guide = guide_legend(override.aes = list(linetype = c("blank", "solid", "solid"),
                                                              shape = c(16, NA, NA)))) + 
  theme(legend.position = "bottom") + guides(col = guide_legend(ncol = 3)) + 
  scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p41 <- ggplot(data = graf_recup1) + 
  geom_line(aes(x = graf_recup1[, 1], y = graf_recup1[, 4], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup1[, 1], y = graf_recup1[, 5], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup1[, 1], y = graf_recup1[, 6], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup1[, 1], y = graf_recup1[, 7], color = "Stochastic\n paths"), size = 0.8) + 
  geom_line(aes(x = graf_recup1[, 1], y = graf_recup1[, 8], color = "Stochastic\n paths"), size = 0.8) + 
  geom_point(aes(x = graf_recup1[, 1], y = graf_recup1[, 3], color = "Recovered\n data"), size = 1.7) +
  geom_line(aes(x = graf_recup1[, 1], y = graf_recup1[, 2], color = "Deterministic\n solution"), size = 1.2) + 
  labs(x = 'Time (days)', y = 'No. of recovered population') + 
  scale_color_manual(values = colors41, name = "", guide = guide_legend(override.aes = list(linetype = 
  c("solid", "solid", "blank"), shape = c(NA, NA, 16)))) + 
  theme(legend.position="bottom") + guides(col = guide_legend(ncol = 3)) + 
  scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

# install.packages("gridExtra")
library(gridExtra)
grid.arrange(p111, p21, ncol = 2)

max(c(SN1[[1]][1:385] - grafica_susceptibles_mejorado[1:385, 2], 
      SN1[[2]][1:385] - grafica_susceptibles_mejorado[1:385, 2], 
      SN1[[3]][1:385] - grafica_susceptibles_mejorado[1:385, 2], 
      SN1[[4]][1:385] - grafica_susceptibles_mejorado[1:385, 2], 
      SN1[[5]][1:385] - grafica_susceptibles_mejorado[1:385, 2]))

max(c(EN1[[1]][1:385] - final_e[1:385], 
      EN1[[2]][1:385] - final_e[1:385], 
      EN1[[3]][1:385] - final_e[1:385], 
      EN1[[4]][1:385] - final_e[1:385], 
      EN1[[5]][1:385] - final_e[1:385]))

IN2 = list(NULL, NULL, NULL, NULL, NULL)
RN2 = list(NULL, NULL, NULL, NULL, NULL)
IN2[[1]][1:385] = proyeccion_bogota_final[1:385] - SN1[[1]][1:385] - EN1[[1]][1:385] - final_r[1:385] 
RN2[[1]][1:385] = proyeccion_bogota_final[1:385] - SN1[[1]][1:385] - EN1[[1]][1:385] - final_i[1:385]
IN2[[2]][1:385] = proyeccion_bogota_final[1:385] - SN1[[2]][1:385] - EN1[[2]][1:385] - final_r[1:385] 
RN2[[2]][1:385] = proyeccion_bogota_final[1:385] - SN1[[2]][1:385] - EN1[[2]][1:385] - final_i[1:385]
IN2[[3]][1:385] = proyeccion_bogota_final[1:385] - SN1[[3]][1:385] - EN1[[3]][1:385] - final_r[1:385] 
RN2[[3]][1:385] = proyeccion_bogota_final[1:385] - SN1[[3]][1:385] - EN1[[3]][1:385] - final_i[1:385]
IN2[[4]][1:385] = proyeccion_bogota_final[1:385] - SN1[[4]][1:385] - EN1[[4]][1:385] - final_r[1:385] 
RN2[[4]][1:385] = proyeccion_bogota_final[1:385] - SN1[[4]][1:385] - EN1[[4]][1:385] - final_i[1:385]
IN2[[5]][1:385] = proyeccion_bogota_final[1:385] - SN1[[5]][1:385] - EN1[[5]][1:385] - final_r[1:385] 
RN2[[5]][1:385] = proyeccion_bogota_final[1:385] - SN1[[5]][1:385] - EN1[[5]][1:385] - final_i[1:385]

graf_infec2 <- data.frame(tiempo[1:385], final_i[1:385], suavizado_infectados$X1[1:385], IN2[[1]][1:385], 
                          IN2[[2]][1:385], IN2[[3]][1:385], IN2[[4]][1:385], IN2[[5]][1:385], casos$Infectados[1:385])
graf_recup2 <- data.frame(tiempo[1:385], final_r[1:385], suavizado_recuperados$X1[1:385], RN2[[1]][1:385],
                          RN2[[2]][1:385], RN2[[3]][1:385], RN2[[4]][1:385], RN2[[5]][1:385], casos$recuperados[1:385])

colors32 <- c("Infected\n data" = "indianred2", "Smoothed\n infected" = 'red4', "Data update\n solution" = 'firebrick1', 
              "Stochastic\n paths" = 'palevioletred1')
colors42 <- c("Recovered data" = "chartreuse3", "Smoothed\n recovered" = 'springgreen4', "Data update\n solution" = 'green3', 
              "Stochastic\n paths" = 'lightgreen')
library(ggplot2)
p32 <- ggplot(data = graf_infec2) + 
  geom_point(aes(x = graf_infec2[, 1], y = graf_infec2[, 9], color = "Infected\n data"), size = 1.7) +
  geom_point(aes(x = graf_infec2[, 1], y = graf_infec2[, 3], color = "Smoothed\n infected"), size = 1.7) + 
  geom_line(aes(x = graf_infec2[, 1], y = graf_infec2[, 2], color = "Data update\n solution"), size = 1, linetype = "dashed") + 
  geom_line(aes(x = graf_infec2[, 1], y = graf_infec2[, 4], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_infec2[, 1], y = graf_infec2[, 5], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_infec2[, 1], y = graf_infec2[, 6], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_infec2[, 1], y = graf_infec2[, 7], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_infec2[, 1], y = graf_infec2[, 8], color = "Stochastic\n paths"), size = 1) + 
  labs(x = 'Time (days)', y = 'No. of infected population') + 
  scale_color_manual(values = colors32, name = "", guide = guide_legend(ncol = 2, override.aes = list(linetype = 
  c("blank", "solid", "dashed", "solid"), shape = c(16, NA, NA, NA)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

p42 <- ggplot(data = graf_recup2) + labs(x = 'Time (days)', y = 'No. of recovered population') + 
  geom_point(aes(x = graf_recup2[, 1], y = graf_recup2[, 9], color = "Recovered data"), size = 1) +
  geom_point(aes(x = graf_recup2[, 1], y = graf_recup2[, 3], color = "Smoothed\n recovered"), size = 1.7) +
  geom_line(aes(x = graf_recup2[, 1], y = graf_recup2[, 2], color = "Data update\n solution"), size = 1, linetype = "dashed") + 
  geom_line(aes(x = graf_recup2[, 1], y = graf_recup2[, 4], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_recup2[, 1], y = graf_recup2[, 5], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_recup2[, 1], y = graf_recup2[, 6], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_recup2[, 1], y = graf_recup2[, 7], color = "Stochastic\n paths"), size = 1) + 
  geom_line(aes(x = graf_recup2[, 1], y = graf_recup2[, 8], color = "Stochastic\n paths"), size = 1) +
  labs(x = 'Time (days)', y = 'No. of recovered population') + 
  scale_color_manual(values = colors42, name = "", guide = guide_legend(ncol = 2, override.aes = list(linetype = 
  c("blank", "solid", "dashed", "solid"), shape = c(16, NA, NA, NA)))) + 
  theme(legend.position="bottom") + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")

library(gridExtra)
grid.arrange(p32, p42, ncol = 2)

max(c(IN2[[1]][1:385] - final_i[1:385], 
      IN2[[2]][1:385] - final_i[1:385], 
      IN2[[3]][1:385] - final_i[1:385], 
      IN2[[4]][1:385] - final_i[1:385], 
      IN2[[5]][1:385] - final_i[1:385]))

max(c(RN2[[1]][1:385] - final_r[1:385], 
      RN2[[2]][1:385] - final_r[1:385], 
      RN2[[3]][1:385] - final_r[1:385], 
      RN2[[4]][1:385] - final_r[1:385], 
      RN2[[5]][1:385] - final_r[1:385]))