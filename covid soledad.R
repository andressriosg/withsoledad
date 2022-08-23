library(readr)
# Casos_positivos_de_COVID_19_en_Colombia <- read_csv("C:/Users/User/Downloads/Casos_positivos_de_COVID-19_en_Colombia.csv")
Salida_Datos_Abiertos <- read_csv("C:/Users/User/Downloads/Salida_Datos_Abiertos.csv")
head(Salida_Datos_Abiertos)
filas_bogota = which(Salida_Datos_Abiertos$Ciudad_municipio_nom == "BOGOTA")
datos_bogota <- Salida_Datos_Abiertos[filas_bogota, ]
head(datos_bogota)
count(datos_bogota)
casos_infectados <- datos_bogota$Fecha_diagnostico
head(casos_infectados)
infectados <- data.frame(transform(table(casos_infectados)))

install.packages("tidyverse")
install.packages("lubridate")
install.packages("datos")
library(tidyverse)
library(lubridate)
library(datos)

infectados$casos_infectados <- as.Date(infectados$casos_infectados, "%d/%m/%Y")
infectados$casos_infectados[order(infectados$casos_infectados)]

infectados <- infectados[order(infectados$casos_infectados), ]
infectados

plot(1:length(infectados$Freq), infectados$Freq, type = "h", col = "red")

casos_recuperados <- datos_bogota$Fecha_recuperado
head(casos_recuperados)
recuperados <- data.frame(transform(table(casos_recuperados)))

recuperados$casos_recuperados <- as.Date(recuperados$casos_recuperados, "%d/%m/%Y")
recuperados$casos_recuperados[order(recuperados$casos_recuperados)]

recuperados <- recuperados[order(recuperados$casos_recuperados), ]
recuperados 

plot(1:length(recuperados$Freq), recuperados$Freq, type = "h", col = "green")

infectados$casos_infectados = as.character(infectados$casos_infectados)
recuperados$casos_recuperados = as.character(recuperados$casos_recuperados)

comunes = intersect(infectados$casos_infectados,
                    recuperados$casos_recuperados)

diferentes <- setdiff(infectados$casos_infectados, recuperados$casos_recuperados)

casos <- data.frame(Fecha = infectados$casos_infectados,
                    Infectados = infectados$Freq,
                    recuperados = NA)

for (i in 1:length(diferentes)) {
  casos[which(casos$Fecha == diferentes[i]), 3] = 0
}
casos
for (i in 1:length(recuperados$casos_recuperados)) {
  casos[which(casos$Fecha == recuperados$casos_recuperados[i]), 3] =
    recuperados$Freq[i]
}
casos
sum(casos$Infectados) # Se tuvo que tomar la fecha de diagnóstico como día en que fueron infectados 
sum(casos$recuperados)

tiempo <- 1:length(casos$Infectados)

casos$Fecha[385]

plot(tiempo, casos$Infectados, col = "red", pch = 16, ylim = c(0, max(casos$recuperados)))
points(tiempo, casos$recuperados, col = "green", pch = 16)
lines(tiempo, suavizado_infectados$X1, col = "darkred", lwd = 2)
lines(tiempo, suavizado_recuperados$X1, col = "darkgreen", lwd = 2)
sum(suavizado_infectados)
sum(suavizado_recuperados)

library(ggplot2)

grafica_infectados <- data.frame(1:385, casos$Infectados[1:385], suavizado_infectados$X1[1:385])
grafica_recuperados <- data.frame(1:385, casos$recuperados[1:385], suavizado_recuperados$X1[1:385])

p_i <- ggplot(data = grafica_infectados) + geom_point(aes(x = grafica_infectados[, 1], y = grafica_infectados[, 2]), 
                                                      col = 'firebrick1', size = 2) + 
  geom_line(aes(x = grafica_infectados[, 1], y = grafica_infectados[, 3]), col = 'black', size = 1) + 
  labs(x = 'Time (days)', y = 'Number of infected population') 
p_r <- ggplot(data = grafica_recuperados) + geom_point(aes(x = grafica_recuperados[, 1], y = grafica_recuperados[, 2]), 
                                                      col = 'darkgreen', size = 2) + 
  geom_line(aes(x = grafica_recuperados[, 1], y = grafica_recuperados[, 3]), col = 'black', size = 1) + 
  labs(x = 'Time (days)', y = 'Number of recovered population')  

library(gridExtra)
library(ggplot2)
grid.arrange(p_i, p_r, ncol = 2)

tiempo <- seq(from = as.Date("2020-03-06"), to = as.Date("2021-03-29"), by=1)
grafica_intervalos <- data.frame(tiempo[1:385], suavizado_infectados$X1[1:385], suavizado_recuperados$X1[1:385], casos$Infectados[1:385], 
                                 casos$recuperados[1:385])

colors <- c("Smoothed infected" = "firebrick3", 
            # "Smoothed recovered" = "darkgreen", 
            "Infected data" = "indianred2" 
            #, "Recovered data" = "chartreuse3"
            )
library(ggplot2)
p_inter <- ggplot(data = grafica_intervalos) + 
  # geom_point(aes(x = grafica_intervalos[, 1],  y = grafica_intervalos[, 5], color = 'Recovered data'), size = 1.5) +
  geom_point(aes(x = grafica_intervalos[, 1],  y = grafica_intervalos[, 4], color = 'Infected data'), size = 1.5)  + 
  # geom_line(aes(x = grafica_intervalos[, 1],  y = grafica_intervalos[, 3], color = 'Smoothed recovered'), size = 1.1)  +
  geom_line(aes(x = grafica_intervalos[, 1], y = grafica_intervalos[, 2], color = 'Smoothed infected'), size = 1.1) +
  #geom_vline(xintercept = tiempo[150]) +
  #geom_vline(xintercept = tiempo[205]) +
  #geom_vline(xintercept = tiempo[270]) +
  #geom_vline(xintercept = tiempo[312]) +
  #geom_vline(xintercept = tiempo[363]) +
  labs(x = 'Time (days)', y = 'Number of population') +  scale_color_manual(values = colors, name = "", 
                                                                                 guide = guide_legend(override.aes = list(
                                                                                   linetype = c("solid", 
                                                                                                # "solid", "blank", 
                                                                                                "blank"),
                                                                                   shape = c(NA, 
                                                                                             # NA, 16, 
                                                                                             16)))) + 
  theme(legend.position="right") + scale_x_date(breaks = "2 months", date_labels = "%d-%m-%Y")

library(gridExtra)
grid.arrange(p_inter, ncol = 1)  

plot(1:420, suavizado_recuperados$X1[1:420], lwd = 3, type = "l", col = "green")
lines(1:420, suavizado_infectados$X1[1:420], lwd = 3, col = "red")
abline(v = 150, col = "purple", lwd = 2)
abline(v = 177, col = "purple", lwd = 2)
abline(v = 205, col = "purple", lwd = 2)
abline(v = 270, col = "purple", lwd = 2)
abline(v = 312, col = "purple", lwd = 2)
abline(v = 330, col = "purple", lwd = 2)
abline(v = 363, col = "purple", lwd = 2)
abline(v = 385, col = "black", lwd = 4)
