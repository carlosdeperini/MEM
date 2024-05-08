#Preparación del Entorno--------------------------------------------------------
install.packages("xtable") #Para convertir matrices a tablas de latex
library(xtable)
library("robustbase")
library("RobStatTM")
library("MASS")
library("fastmatrix") #Esto requiere L1pack según su documentación
library("L1pack")

rm(list = ls())

#Ejercicio 10-------------------------------------------------------------------
###Item a)----

#Parámetro que se puede ajustar
epsilon <- 0.20 #Este parámetro se cambia entre 0.01, 0.05, 0.10 y 0.20

#Definicion de parametros que no varían
n <- 10000 #Tamaño de la muestra
nRep <- 500 #Cantidad de repeticiones
tau <- 3:20 #Desvío de la normal contaminante
mu <- 0 #Media de ambas normales

#Definición de datos
matriz_resultado <- matrix(nrow = length(tau), ncol = 2)
colnames(matriz_resultado) <- c("n*Var_media", "n*Var_mediana")
set.seed(123) #Semilla para repetibilidad

#Construcción de la tabla
for (tau_contador in tau) {

  #Se guardan los resultados
  matriz_resultado[tau_contador-2, 1] <- (1-epsilon) + epsilon*tau_contador^2 #Varianza del promedio teórica
  matriz_resultado[tau_contador-2, 2] <- pi/(2*(1-epsilon + epsilon/tau_contador)^2) #Varianza de la mediana teórica
}

print(xtable(matriz_resultado), include.rownames = FALSE)

###Item b)----

rm(list = ls())

#Funciones
generador_contaminado <- function(n, mu, tau, epsilon){
  var_binomial <- rbinom(1, 1, epsilon) #Vale 1 si está contaminada
  normal_centrada <- rnorm(n, mu, 1) #Normal conocida
  normal_contaminante <- rnorm(n, mu, tau) #Normal que contamina
  
  muestra <- c()
  if (var_binomial == 1) { #Si da 1, quiere decir que está contaminada
    muestra <- normal_contaminante
  } else{
    muestra <- normal_centrada
  }
  return(muestra)
}

#Definicion de parametros que no varían
n <- 10000 #Tamaño de la muestra
nRep <- 500 #Cantidad de repeticiones
tau <- 10 #Desvío de la normal contaminante
mu <- 0 #Media de ambas normales
epsilon <- c(0.01, 0.05, 0.10, 0.20)

#Definición de datos
matriz_auxiliar <- matrix(nrow = nRep, ncol = 1)
matriz_resultado <- matrix(nrow = length(epsilon), ncol = 3)
set.seed(123) #Semilla para repetibilidad
contador <- 1

#Construcción de la tabla
for (eps in epsilon) {
  
  #Se genera la matriz de muestras contaminadas
  for (i in 1:nRep) {
    muestra_contaminada <- generador_contaminado(n, mu, tau, eps)
    matriz_auxiliar[i, 1] <- mean(muestra_contaminada, trim = 0.05)
  }
  
  #Se guardan los resultados
  matriz_resultado[contador, 1] <- (1-eps) + eps*tau^2 #Varianza del promedio teórica
  matriz_resultado[contador, 2] <- pi/(2*(1-eps + eps/tau)^2) #Varianza de la mediana teórica
  matriz_resultado[contador, 3] <- n*var(matriz_auxiliar[,1]) #Varianza de la mediana
  
  contador <- contador + 1
}

print(xtable(matriz_resultado), include.rownames = FALSE)

#Ejercicio 21-------------------------------------------------------------------

rm(list = ls())

###Item a)----
data("mineral") #Carga de datos

#Gráfico
plot(mineral$copper, mineral$zinc, col = "red", pch = 20,
     main = "Rectas de Regresión y Puntos del DF Mineral",
     xlab = "Cobre", ylab = "Zinc")
grid()

###Item b)----

#Se ajustan los modelos
ajuste_LS <- lm(mineral$zinc~mineral$copper)
ajuste_M <- rlm(mineral$zinc~mineral$copper, method="M")
ajuste_MM <- rlm(mineral$zinc~mineral$copper, method="MM")

cobre_ext <- cbind(rep(1, length(mineral$copper)), mineral$copper) #Hay que agregar columnas de 1.
ajuste_L1 <- lmrob.lar(cobre_ext,mineral$zinc)

#Se realizan los gráficos de las rectas obtenidas

#Ajuste LS
y_hat_ajuste_LS <- ajuste_LS$fitted.values
lines(mineral$copper, y_hat_ajuste_LS, col="blue", lwd = 2)

#Ajuste L1
zinc_hat_L1 <- mineral$zinc - ajuste_L1$residuals
lines(mineral$copper, zinc_hat_L1, col="#BF3EFF", lwd = 2)

#Ajuste M
y_hat_ajuste_M <- ajuste_M$fitted.values
lines(mineral$copper, y_hat_ajuste_M, col="black", lwd = 2)

#Ajuste MM
y_hat_ajuste_MM <- ajuste_MM$fitted.values
lines(mineral$copper, y_hat_ajuste_MM, col="green", lwd = 2)

legend("topleft",
       legend=c("LS", "L1", "M-estimador Tukey", "MM-estimador Tukey"),
       col=c("blue", "#BF3EFF", "black", "green"),
       lwd=2,
       cex=0.8) 


###Item c)----

#QQPlot
qqnorm(residuals(ajuste_LS), main = "QQ-plot de residuos (LS)")
qqline(residuals(ajuste_LS), col = "steelblue")
grid()

#Boxplot
boxplot(ajuste_LS$residuals, col = "lightblue", 
        main = "Boxplot de los residuos del modelo LS",
        ylab = "Residuos",
        outcol = "red",pch = 19)
grid()

###Item d)----

#QQPlot
qqnorm(residuals(ajuste_MM), main = "QQ-plot de residuos (MM)")
qqline(residuals(ajuste_MM), col = "steelblue")
grid()


#Boxplot
boxplot(ajuste_MM$residuals, col = "lightblue", 
        main = "Boxplot de los residuos del modelo MM",
        ylab = "Residuos",
        outcol = "red",pch = 19)
grid()

###Item e)----

#Eliminación del outlier y creación de nuevo data frame
maximo_cobre <- which.max(mineral$copper) #Posición en el df del outlier
mineral_sin_outlier <- mineral[-maximo_cobre,] #Data sin Outlier extremo
rownames(mineral_sin_outlier) <- NULL #Se reinicia la numeración del df

#Ajuste sin el outlier
ajuste_LS_sin_outlier <- lm(zinc~copper, data = mineral_sin_outlier)

#Gráficos
x_grilla <- seq(min(mineral$copper),max(mineral$copper), length.out = 100)
y_hat_ajuste_LS_sin_outliers <- predict(ajuste_LS_sin_outlier, newdata = data.frame(copper = x_grilla))
lines(x_grilla, y_hat_ajuste_LS_sin_outliers, col="brown", lwd = 2)

legend("topleft",
       legend=c("LS", "L1", "M-estimador Tukey", "MM-estimador Tukey", "LS sin outlier"),
       col=c("blue", "#BF3EFF", "black", "green", "brown"),
       lwd=2,
       cex=0.8) 





