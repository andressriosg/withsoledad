# Vamos a tomar los residuales 

residuales_susceptibles_1 <- graf_suscp1[, 4] - susceptibles_bogota[1:385]
residuales_susceptibles_2 <- graf_suscp1[, 5] - susceptibles_bogota[1:385]
residuales_susceptibles_3 <- graf_suscp1[, 6] - susceptibles_bogota[1:385]
residuales_susceptibles_4 <- graf_suscp1[, 7] - susceptibles_bogota[1:385]
residuales_susceptibles_5 <- graf_suscp1[, 8] - susceptibles_bogota[1:385]
residuales_expuestos_1 <- graf_expue1[, 4] - expuestos_bogota[1:385]
residuales_expuestos_2 <- graf_expue1[, 5] - expuestos_bogota[1:385]
residuales_expuestos_3 <- graf_expue1[, 6] - expuestos_bogota[1:385]
residuales_expuestos_4 <- graf_expue1[, 7] - expuestos_bogota[1:385]
residuales_expuestos_5 <- graf_expue1[, 8] - expuestos_bogota[1:385]
residuales_infectados_1 <- graf_infec1[, 4] - casos$Infectados[1:385]
residuales_infectados_2 <- graf_infec1[, 5] - casos$Infectados[1:385]
residuales_infectados_3 <- graf_infec1[, 6] - casos$Infectados[1:385]
residuales_infectados_4 <- graf_infec1[, 7] - casos$Infectados[1:385]
residuales_infectados_5 <- graf_infec1[, 8] - casos$Infectados[1:385]
residuales_recuperados_1 <- graf_recup1[, 4] - casos$recuperados[1:385]
residuales_recuperados_2 <- graf_recup1[, 5] - casos$recuperados[1:385]
residuales_recuperados_3 <- graf_recup1[, 6] - casos$recuperados[1:385]
residuales_recuperados_4 <- graf_recup1[, 7] - casos$recuperados[1:385]
residuales_recuperados_5 <- graf_recup1[, 8] - casos$recuperados[1:385]

plot(susceptibles_bogota[1:385], residuales_susceptibles) # Para ver homocedasticidad 
residuales_susceptibles <- data.frame(susceptibles_bogota[1:385], residuales_susceptibles_1, residuales_susceptibles_2, 
                                      residuales_susceptibles_3, residuales_susceptibles_4, 
                                      residuales_susceptibles_5)
residuales_expuestos <- data.frame(expuestos_bogota[1:385], residuales_expuestos_1, residuales_expuestos_2, 
                                      residuales_expuestos_3, residuales_expuestos_4, 
                                      residuales_expuestos_5)
residuales_infectados <- data.frame(casos$Infectados[1:385], residuales_infectados_1, residuales_infectados_2, 
                                   residuales_infectados_3, residuales_infectados_4, 
                                   residuales_infectados_5)
residuales_recuperados <- data.frame(casos$recuperados[1:385], residuales_recuperados_1, residuales_recuperados_2, 
                                    residuales_recuperados_3, residuales_recuperados_4, 
                                    residuales_recuperados_5)

colors_residuales_sup <- c("Residuals\n 1st path" = 'royalblue1', "Residuals\n 2nd path" = "royalblue2", 
              "Residuals\n 3rd path" = 'royalblue3', "Residuals\n 4th path" = "steelblue3", 
              "Residuals\n 5th path" = "steelblue4")
plot_resi_sus <- ggplot(data = residuales_susceptibles) + 
  geom_point(aes(x = residuales_susceptibles[, 1], y = residuales_susceptibles[, 2], color = "Residuals\n 1st path"), size = 1.2) + 
  geom_point(aes(x = residuales_susceptibles[, 1], y = residuales_susceptibles[, 3], color = "Residuals\n 2nd path"), size = 1.2) + 
  geom_point(aes(x = residuales_susceptibles[, 1], y = residuales_susceptibles[, 4], color = "Residuals\n 3rd path"), size = 1.2) + 
  geom_point(aes(x = residuales_susceptibles[, 1], y = residuales_susceptibles[, 5], color = "Residuals\n 4th path"), size = 1.2) + 
  geom_point(aes(x = residuales_susceptibles[, 1], y = residuales_susceptibles[, 6], color = "Residuals\n 5th path"), size = 1.2) + 
  geom_hline(yintercept = 0, color = "blue", size = 0.8) + 
  labs(x = 'Estimated susceptible data', y = 'Residuals') +
  theme(legend.position="bottom") + 
  scale_color_manual(values = colors_residuales_sup, name = "", guide = guide_legend(ncol = 3, 
                     override.aes = list(linetype = c(rep("solid", 5)), shape = c(rep(16, 5))))) 
colors_residuales_exp <- c("Residuals\n 1st path" = 'darkorange', "Residuals\n 2nd path" = "darkorange1", 
                           "Residuals\n 3rd path" = 'tan2', "Residuals\n 4th path" = "tan1", 
                           "Residuals\n 5th path" = "darkorange2")
plot_resi_exp <- ggplot(data = residuales_expuestos) + 
  geom_point(aes(x = residuales_expuestos[, 1], y = residuales_expuestos[, 2], color = "Residuals\n 1st path"), size = 1.2) + 
  geom_point(aes(x = residuales_expuestos[, 1], y = residuales_expuestos[, 3], color = "Residuals\n 2nd path"), size = 1.2) + 
  geom_point(aes(x = residuales_expuestos[, 1], y = residuales_expuestos[, 4], color = "Residuals\n 3rd path"), size = 1.2) + 
  geom_point(aes(x = residuales_expuestos[, 1], y = residuales_expuestos[, 5], color = "Residuals\n 4th path"), size = 1.2) + 
  geom_point(aes(x = residuales_expuestos[, 1], y = residuales_expuestos[, 6], color = "Residuals\n 5th path"), size = 1.2) + 
  geom_hline(yintercept = 0, color = "chocolate", size = 0.8) +
  labs(x = 'Estimated exposed data', y = 'Residuals') + 
  theme(legend.position="bottom") + 
  scale_color_manual(values = colors_residuales_exp, name = "", guide = guide_legend(ncol = 3, 
                     override.aes = list(linetype = c(rep("solid", 5)), shape = c(rep(16, 5))))) 
colors_residuales_inf <- c("Residuals\n 1st path" = 'brown1', "Residuals\n 2nd path" = "brown2", 
                           "Residuals\n 3rd path" = 'brown3', "Residuals\n 4th path" = "tomato1", 
                           "Residuals\n 5th path" = "tomato2")
plot_resi_inf <- ggplot(data = residuales_infectados) + 
  geom_point(aes(x = residuales_infectados[, 1], y = residuales_infectados[, 2], color = "Residuals\n 1st path"), size = 1.2) + 
  geom_point(aes(x = residuales_infectados[, 1], y = residuales_infectados[, 3], color = "Residuals\n 2nd path"), size = 1.2) + 
  geom_point(aes(x = residuales_infectados[, 1], y = residuales_infectados[, 4], color = "Residuals\n 3rd path"), size = 1.2) + 
  geom_point(aes(x = residuales_infectados[, 1], y = residuales_infectados[, 5], color = "Residuals\n 4th path"), size = 1.2) + 
  geom_point(aes(x = residuales_infectados[, 1], y = residuales_infectados[, 6], color = "Residuals\n 5th path"), size = 1.2) +  
  geom_hline(yintercept = 0, color = "red", size = 0.8) +
  labs(x = 'Infected data', y = 'Residuals') + 
  theme(legend.position="bottom") + 
  scale_color_manual(values = colors_residuales_inf, name = "", guide = guide_legend(ncol = 3, 
                     override.aes = list(linetype = c(rep("solid", 5)), shape = c(rep(16, 5))))) 
colors_residuales_rec <- c("Residuals\n 1st path" = 'green3', "Residuals\n 2nd path" = "limegreen", 
                           "Residuals\n 3rd path" = 'springgreen3', "Residuals\n 4th path" = "seagreen3", 
                           "Residuals\n 5th path" = "chartreuse3")
plot_resi_rec <- ggplot(data = residuales_recuperados) + 
  geom_point(aes(x = residuales_recuperados[, 1], y = residuales_recuperados[, 2], color = "Residuals\n 1st path"), size = 1.2) + 
  geom_point(aes(x = residuales_recuperados[, 1], y = residuales_recuperados[, 3], color = "Residuals\n 2nd path"), size = 1.2) + 
  geom_point(aes(x = residuales_recuperados[, 1], y = residuales_recuperados[, 4], color = "Residuals\n 3rd path"), size = 1.2) + 
  geom_point(aes(x = residuales_recuperados[, 1], y = residuales_recuperados[, 5], color = "Residuals\n 4th path"), size = 1.2) + 
  geom_point(aes(x = residuales_recuperados[, 1], y = residuales_recuperados[, 6], color = "Residuals\n 5th path"), size = 1.2) +  
  geom_hline(yintercept = 0, color = "green", size = 0.8) +
  labs(x = 'Recovered data', y = 'Residuals') + 
  theme(legend.position="bottom") + 
  scale_color_manual(values = colors_residuales_rec, name = "", guide = guide_legend(ncol = 3, 
                     override.aes = list(linetype = c(rep("solid", 5)), shape = c(rep(16, 5))))) 

grid.arrange(plot_resi_sus, plot_resi_exp, plot_resi_inf, plot_resi_rec, ncol = 2) # Heterocedasticidad 

c(shapiro.test(residuales_susceptibles_1)$p.value,  # Normalidad 
shapiro.test(residuales_susceptibles_2)$p.value, 
shapiro.test(residuales_susceptibles_3)$p.value, 
shapiro.test(residuales_susceptibles_4)$p.value, 
shapiro.test(residuales_susceptibles_5)$p.value)

c(shapiro.test(residuales_expuestos_1)$p.value,  # Normalidad 
  shapiro.test(residuales_expuestos_2)$p.value, 
  shapiro.test(residuales_expuestos_3)$p.value, 
  shapiro.test(residuales_expuestos_4)$p.value, 
  shapiro.test(residuales_expuestos_5)$p.value)

c(shapiro.test(residuales_infectados_1)$p.value,  # Normalidad 
  shapiro.test(residuales_infectados_2)$p.value, 
  shapiro.test(residuales_infectados_3)$p.value, 
  shapiro.test(residuales_infectados_4)$p.value, 
  shapiro.test(residuales_infectados_5)$p.value)

c(shapiro.test(residuales_recuperados_1)$p.value,  # Normalidad 
  shapiro.test(residuales_recuperados_2)$p.value, 
  shapiro.test(residuales_recuperados_3)$p.value, 
  shapiro.test(residuales_recuperados_4)$p.value, 
  shapiro.test(residuales_recuperados_5)$p.value)
