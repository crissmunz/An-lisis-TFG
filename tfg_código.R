install.packages("dplyr")
install.packages("readr")
install.packages("tseries")
install.packages("Kendall")
install.packages("e1071")
install.packages("psych")
install.packages("ggplot2")
install.packages("FSA")
install.packages("writexl")

library(dplyr)
library(readr)
library(tseries)
library(Kendall)
library(e1071)
library(psych)
library(ggplot2)
library(FSA)
library(writexl)
#citaciones de las librerias utilizadas
citation("dplyr")
citation("readr")
citation("tseries")
citation("Kendall")
citation("e1071")
citation("psych")
citation("ggplot2")
citation("FSA")
citation("writexl")

#Cargar la base de datos
data<-read.csv("/Users/cristinamunoz/Downloads/heart_failure_clinical_records_dataset.csv")

#Cálculo de los estadísticos descriptivos
str(data)
head(data)
descriptive_stats<-describe(data)
descriptive_stats<-as.data.frame(descriptive_stats)
descriptive_stats$Variable<-rownames(descriptive_stats)
print(names(descriptive_stats))

descriptive_stats$Q1<-apply(data,2,quantile,probs = 0.25, na.rm=TRUE)
descriptive_stats$Q3<-apply(data,2,quantile,probs = 0.75, na.rm=TRUE)
descriptive_stats$Rango_Intercuartílico<-descriptive_stats$Q3-descriptive_stats$Q1

descriptive_stats<-descriptive_stats %>%
  rename(
         Media = mean, 
         Mediana = median,
         Desciación_típica = sd,
         Mínimo = min,
         Máximo = max,
         Asimetría = skew)
descriptive_stats<-descriptive_stats%>%
  select(Variable,Media,Mediana,Desciación_típica, Mínimo, Máximo, descriptive_stats$Rango_Intercuartílico, Asimetría)
print(descriptive_stats)
#descarga de la tabla del análisis descriptivo 
write.csv(descriptive_stats,"/Users/cristinamunoz/Downloads/descriptive_statistics.csv", row.names = FALSE)
# Calcular el número de pacientes con diabetes y/o anemia
num_pacientes_diabetes_anemia <- sum(data$diabetes == 1 | data$anaemia == 1)
total_pacientes <- nrow(data)
porcentaje_diabetes_anemia <- (num_pacientes_diabetes_anemia / total_pacientes) * 100
porcentaje_diabetes_anemia


#PRUEBA DE RACHAS
# No normalidad de los datos
shapiro.test(datos$serum_sodium)
# Seleccionar la variable 'serum_sodium' para el test de rachas
mediana <- median(datos$serum_sodium, na.rm = TRUE)
variable_binario <- ifelse(datos$serum_sodium > mediana, 1, 0)
# Función para contar el número total de rachas
cuenta_rachas <- function(series) {
  n <- length(series)
  if (n == 0) return(0)
  rachas <- 1
  for (i in 2:n) {
    if (series[i] != series[i - 1]) {
      rachas <- rachas + 1
    }
  }
  return(rachas)
}

total_rachas <- cuenta_rachas(variable_binario)
n <- sum(variable_binario == 1)
m <- sum(variable_binario == 0)
N <- n + m

esperanza_rachas <- 1+((2 * n * m) / N )
varianza_rachas <- (2 * n * m * (2 * n * m - N)) / (N^2 * (N - 1))

z_valor <- (total_rachas - esperanza_rachas) / sqrt(varianza_rachas)
z_valor

p_valor <- 2 * pnorm(abs(z_valor), lower.tail = FALSE)

list(
  R = total_rachas,
  Esperanza_Rachas = esperanza_rachas,
  Varianza_Rachas = varianza_rachas,
  Z_Valor = z_valor,
  p_Valor= p_valor
)
#gráfico
datos <- datos %>%
  mutate(suero_sod_bin = ifelse(serum_sodium > mediana, "Supera la mediana", "No supera la mediana"))

ggplot(datos, aes(x = 1:nrow(datos), y = serum_sodium, color = suero_sod_bin)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("Supera la mediana" = "blue", "No supera la mediana" = "red")) +
  labs(title = "Secuencia de Valores de Sodio Sérico",
       x = "Índice de Observación",
       y = "Sodio Sérico (mEq/L)") +
  theme_minimal()


#MANN-WHITNEY
data_clean1 <- data %>% filter(!is.na(creatinine_phosphokinase) & !is.na(anaemia))
#definición de las variables a comprarar
tratamiento1 <- data_clean$creatinine_phosphokinase[data_clean1$anaemia == 1]
tratamiento2 <- data_clean$creatinine_phosphokinase[data_clean1$anaemia == 0]
#No normalidad datos
shapiro.test(tratamiento1)
shapiro.test(tratamiento2)
#prueba
resultado_mann_whitney <- wilcox.test(tratamiento1, tratamiento2)
resultado_mann_whitney
W<-resultado_mann_whitney$statistic
#Valor Z aproximado
n1_mann_whitney <- length(tratamiento1)
n2_mann_whitney <- length(tratamiento2)
mu_mann_whitney <- n1_mann_whitney * n2_mann_whitney / 2
sigma_mann_whitney <- sqrt(n1_mann_whitney * n2_mann_whitney * (n1_mann_whitney + n2_mann_whitney + 1) / 12)
Z_mann_whitney <- (W - mu_mann_whitney) / sigma_mann_whitney
Z_mann_whitney
#tabla/gráfico
data$anaemia_factor <- as.factor(data$anaemia)
ggplot(data, aes(x = anaemia_factor, y = creatinine_phosphokinase)) +
  geom_boxplot() +
  labs(title = "Niveles de CPK según Presencia de Anemia",
       x = "Anemia (0 = No, 1 = Sí)",
       y = "Niveles de CPK (U/L)")


#KENDALL
data_clean2 <- data %>% filter(!is.na(age) & !is.na(ejection_fraction))
#definición variables
variable1 <- data_clean2$age
variable2 <- data_clean2$ejection_fraction
#No normalidad datos
shapiro.test(variable1)
shapiro.test(variable2)


#prueba
resultado_kendall <- cor.test(variable1, variable2, method = "kendall")
resultado_kendall
#gráfico
ggplot(data_clean2, aes(x = age, y = ejection_fraction)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Relación entre Edad y Fracción de Eyección", x = "Edad", y = "Fracción de Eyección (%)") +
  theme_minimal()



#KRUSKAL-WALLIS
data$grupos_edad <- cut(data$age, breaks = c(0,40, 50, 60, 70, 80, 100), 
                      labels = c("0-40","41-50", "51-60", "61-70", "71-80", "81-100"))
#No normalidad datos
shapiro.test(data$serum_creatinine)
shapiro.test(data$age)
# Realizar el test de Kruskal-Wallis
resultados_KW <- kruskal.test(serum_creatinine ~ grupos_edad, data = data)
H<-resultados_KW$statistic
print(resultados_KW)
print(H)
#gráfico
ggplot(data, aes(x = grupos_edad, y = serum_creatinine)) +
  geom_boxplot() +
  labs(title = "Niveles de Creatinina Sérica según Rangos de Edad",
       x = "Grupo de Edad",
       y = "Creatinina Sérica (mg/dL)")
#Prueba post hoc
dunn_results <- dunnTest(serum_creatinine ~ grupos_edad, data = data, method = "bonferroni")
print(resultados_dunn)
tabla_resultados_dunn <- resultados_dunn$res
write_xlsx(tabla_resultados_dunn, "tabla_resultados_dunn.xlsx")
