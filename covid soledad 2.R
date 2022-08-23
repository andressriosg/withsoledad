### Estimaci?n de los susceptibles y de los infectados

## Estimaci?n del par?metros sigma por m?mimos cuadrados ordinarios

# Determinar los puntos de inflexion. Una vez los determinamos procedemos a estimar
# Los valores para cada uno de los puntos

install.packages("npregfast")
library(npregfast)
suavizadoi <- frfast(casos$Infectados ~ tiempo, model = "np", smooth = "kernel", kbin = length(casos$Infectados), 
                     p =3)
plot(suavizadoi)
suavizado_infectados <- data.frame(suavizadoi$p)
sum(suavizado_infectados$X1)

suavizador <- frfast(casos$recuperados ~ tiempo, model = "np", smooth = "kernel", kbin = length(casos$recuperados), 
                     p =3)
plot(suavizador)
suavizado_recuperados <- data.frame(suavizador$p)
sum(suavizado_recuperados$X1)
 
# Proyecci?n de las poblaciones: 

library(readxl)
View(anexo_proyecciones_poblacion_bogota_desagreacion_loc_2018_2035_UPZ_2018_2024)
anexo_proyecciones_poblacion_bogota_desagreacion_loc_2018_2035_UPZ_2018_2024 <- read_excel("C:/Users/User/Downloads/anexo-proyecciones-poblacion-bogota-desagreacion-loc-2018-2035-UPZ-2018-2024.xlsx")
proyeccion_bog <- data.frame(cbind(anexo_proyecciones_poblacion_bogota_desagreacion_loc_2018_2035_UPZ_2018_2024$...3, 
                             anexo_proyecciones_poblacion_bogota_desagreacion_loc_2018_2035_UPZ_2018_2024$...4,
                             anexo_proyecciones_poblacion_bogota_desagreacion_loc_2018_2035_UPZ_2018_2024$...310))
proyeccion_bog <- na.omit(proyeccion_bog)

proyeccion_bog1 = NULL
for (i in 1:2035-2017) {
  proyeccion_bog1[i] <- sum(as.numeric(proyeccion_bog[proyeccion_bog[, 2] == i + 2017 & 
                                                      proyeccion_bog[, 1] != "Total", 3]))
}

par(mfrow=c(1,1))
proyeccion_bogota <- data.frame(2018:2035, proyeccion_bog1)
plot(proyeccion_bogota[, 1], proyeccion_bogota[, 2], type = "b", col = 3, lwd = 2)

install.packages("np")
library(np)
bogota_np <- npreg(proyeccion_bogota[, 2]~ proyeccion_bogota[, 1], eydat = seq(2018, 2035, by = 1/365), 
                   exdat = seq(2018,2035, by = 1/365), bws = 0.5)

plot(bogota_np$bws, lwd = 2, col = "blue")
points(proyeccion_bogota[, 1], proyeccion_bogota[, 2], col = 5, pch = 16)

library(npregfast)
suavizadopb <- frfast(proyeccion_bogota[, 2]~ proyeccion_bogota[, 1], model = "np", smooth = "kernel", 
                     kbin = 6206, p =3)
suavizadopb <- data.frame(suavizadopb$p)

par(mfrow=c(1,1))
plot(seq(2018, 2035, by = 1/365), suavizadopb$X1, type = "l",lwd = 2, col = "blue", xlim = c(2018, 2021)
     , ylim = c(7400000, 7800000))
lines(proyeccion_bogota[, 1], proyeccion_bogota[, 2], col = 4, lwd = 2)

proyeccion_bogota_final <- suavizadopb$X1[(365*2+ 58): (365*2+ 58 + length(casos$Infectados))] # Total de la poblacion por unidad de tiempo 
plot(proyeccion_bogota_final, type = "l", lwd = 2, col = "darkblue")

# Estimaci?n de los expuestos y los susceptibles 

delta_poblacion <- NULL 
for (j in 1:length(proyeccion_bogota_final)) {
  delta_poblacion[j] = proyeccion_bogota_final[j+1] - proyeccion_bogota_final[j]
}
plot(delta_poblacion, type = "l")

delta_infectados <- NULL 
for (j in 1:421) {
  delta_infectados[j] = suavizado_infectados$X1[j+1] - suavizado_infectados$X1[j]
}
plot(delta_infectados)

delta_recuperados <- NULL 
for (j in 1:421) {
  delta_recuperados[j] = suavizado_recuperados$X1[j+1] - suavizado_recuperados$X1[j]
}
plot(delta_recuperados)

# https://www.dane.gov.co/index.php/estadisticas-por-tema/demografia-y-poblacion/nacimientos-y-defunciones
mu = 4/1000 # https://saludata.saludcapital.gov.co/osb/index.php/datos-de-salud/demografia/tm-bruta/
# 4 muertes por cada 1000 habitantes 
# mean(759/proyeccion_bogota_final) https://datosmacro.expansion.com/demografia/mortalidad/colombia 759 muertos dia
eta = (9.42/1000)*mean(proyeccion_bogota_final)  # max(as.numeric(delta_poblacion[1:420])) # https://datosmacro.expansion.com/
# 9.42 nacimientos por cada 1000 habitantes https://saludata.saludcapital.gov.co/osb/index.php/datos-de-salud/demografia/natalidad/
upsilon = 1/5.2 # https://elmedicointeractivo.com/el-periodo-de-incubacion-de-covid-19-se-situa-en-51-dias/

expuestos_bogota <- ((delta_infectados[1:420] + delta_recuperados[1:420] - delta_poblacion[1:420])/upsilon) + (eta/upsilon) - 
                    ((mu/upsilon)*(proyeccion_bogota_final[1:420] - suavizado_infectados$X1[1:420]  - suavizado_recuperados$X1[1:420]))
susceptibles_bogota = (proyeccion_bogota_final[1:420] - expuestos_bogota - suavizado_recuperados$X1[1:420] - 
                         suavizado_infectados$X1[1:420])

plot(proyeccion_bogota_final, type = "l", lwd = 2, col = "brown", ylim = c(min(susceptibles_bogota), max(proyeccion_bogota_final)))
lines(susceptibles_bogota, type = "l", lwd = 2, col = "blue")

plot(expuestos_bogota, type = "l", lwd = 2, col = "orange")

# susceptibles_bogota <- 0.2*(proyeccion_bogota_final[1:420] - expuestos_bogota - suavizado_recuperados$X1[1:420] - suavizado_infectados$X1[1:420])
# plot(susceptibles_bogota, type = "l", lwd = 2, col = "green")

#### REVISAR LA ESTIMACIÓN DE LOS SUSCEPTIBLES. ANOTAR UN CRECIMIENTO LOGÍSTICO 

susceptibles_expuestos <- (proyeccion_bogota_final[1:420] - suavizado_recuperados$X1[1:420] - suavizado_infectados$X1[1:420])
plot(susceptibles_expuestos)

(- delta_poblacion + 2*delta_recuperados + 2*delta_infectados)/(2*upsilon) + eta/(2*upsilon) - mu/(2*upsilon)*(proyeccion_bogota_final[1:420] - 2*suavizado_recuperados$X1[1:420] - 2*suavizado_infectados$X1[1:420])


# Estimaci?n por m?nimos cuadrados ordinarios

beta = ((1-mu)*sum(((susceptibles_bogota[1:419])^2*suavizado_infectados$X1[1:419])) - sum(susceptibles_bogota[2:420]*susceptibles_bogota[1:419]*suavizado_infectados$X1[1:419])  + eta*sum(susceptibles_bogota[1:419]*suavizado_infectados$X1[1:419]) - (1 - upsilon - mu)*sum(susceptibles_bogota[1:419]*expuestos_bogota[1:419]*suavizado_infectados$X1[1:419]) + sum(expuestos_bogota[2:420]*susceptibles_bogota[1:419]*suavizado_infectados$X1[1:419]))/(2*sum((susceptibles_bogota[1:419])^2*(suavizado_infectados$X1[1:419])^2))
gamma = ((1-mu)*(sum((suavizado_infectados$X1[1:419])^2) - sum(suavizado_recuperados$X1[1:419]*suavizado_infectados$X1[1:419])) - sum(suavizado_infectados$X1[1:419]*suavizado_infectados$X1[2:420]) + sum(suavizado_infectados$X1[1:419]*suavizado_recuperados$X1[2:420]) + upsilon*sum(expuestos_bogota[1:419]*suavizado_infectados$X1[1:419]))/(2*sum((suavizado_infectados$X1[1:419])^2))

# Estimar los susceptibles vía ecuación logísitica 
# Usando las estimaciones dadas: 

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

library("xlsx")
write.csv(proyeccion_bogota_final, file = "proyeccion_bogota_final.csv")

x0 <- c(7743953, 0, 1, 0)
times = 1:419
parametros = c(eta, beta, 1/5.1, gamma, mu)  

install.packages("deSolve")
library(deSolve)
solution <- ode(y=x0, times, seir, parametros)

plot(solution[, 2], col = 2, lwd = 2, type ="l")
lines(1:420, susceptibles_bogota, lwd = 2, col = "blue")

plot(solution[, 3], col = 2, lwd = 2, type ="l")
lines(1:420, expuestos_bogota, lwd = 2, col = "blue")

plot(solution[, 4], col = 2, lwd = 2, type ="l")
lines(1:420, suavizado_infectados$X1[1:420], lwd = 2, col = "blue")

plot(solution[, 5], col = 2, lwd = 2, type ="l")
lines(1:420, suavizado_recuperados$X1[1:420], lwd = 2, col = "blue")

plot(1:420, susceptibles_bogota, lwd = 2, col = "blue")
plot(1:420, expuestos_bogota, lwd = 2, col = "blue")

par(mfrow=c(2,2))
plot(susceptibles_bogota)
plot(expuestos_bogota)
plot(suavizado_infectados$X1)
plot(suavizado_recuperados$X1)

sum(suavizado_infectados$X1[1:385])
sum(suavizado_recuperados$X1[1:400]) # La base de datos está teniendo problemas. Hay más recuperados que infectados 
# Se están bajo-estimando los infectados 

grafica_poblacion <- data.frame(tiempo[1:385], proyeccion_bogota_final[1:385])
grafica_susceptibles <- data.frame(tiempo[1:385], susceptibles_bogota[1:385])
grafica_expuestos <- data.frame(tiempo[1:385], expuestos_bogota[1:385])

p_p <- ggplot(data = grafica_poblacion) + geom_line(aes(x = grafica_poblacion[, 1], y = grafica_poblacion[, 2]), 
                                                      col = 'brown', size = 1) +
  labs(x = 'Time (days)', y = 'Number of total population')  + scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")
p_s <- ggplot(data = grafica_susceptibles) + geom_line(aes(x = grafica_susceptibles[, 1], y = grafica_susceptibles[, 2]), 
                                                     col = 'darkblue', size = 1) +
       labs(x = 'Time (days)', y = 'Number of susceptible population') + 
  scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y")
p_e <- ggplot(data = grafica_expuestos) + geom_line(aes(x = grafica_expuestos[, 1], y = grafica_expuestos[, 2]), 
                                                       col = 'darkorange', size = 1) + 
  labs(x = 'Time (days)', y = 'Number of exposed population') +  
  scale_x_date(breaks = "3 months", date_labels = "%d-%m-%Y") 

library(gridExtra)
grid.arrange(p_p, p_s, p_e, ncol = 2)
