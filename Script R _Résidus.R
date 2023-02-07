####script de Mathieu Mortz 2022#####

library(readxl)
library(ggplot2)
#vous avez pas besoin de toutes ces librairies, je sais juste pas lesquels sont les bonnes ;)

#Allometry <- read_excel(file.choose(),1)

Allometry$logO2=log10(Allometry$O2) #changer en log
Allometry$logPoids=log10(Allometry$Poids)

model <- lm(logPoids~logO2, data=Allometry) #fonction linéaire
summary(model)
res <- resid(model)
res #regarde tes résidues si ça fit
plot(fitted(model), res)
abline(0,0)

Allometry$predicted <- predict(model) #valeurs prédites
Allometry$residuals <- residuals(model) #rentrer tes résidues dans ton tableau initial
Allometry

par(mfrow = c(2, 2)) #split 
plot(model) #plus ou moins utile pour nous

#deuxième méthode pour un visuel avec l'équation de la droite de régression
lm_eqn <- function(df){
  m <- lm(logO2 ~ logPoids, df);
  eq <- substitute(italic(logO2) == a + b %.% italic(logPoids)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

rhodo <- ggplot(Allometry, aes(x = logPoids, y = logO2)) + 
  geom_point() + geom_smooth(method="lm")  + geom_text(x = 0, y = 2.2, label = lm_eqn(Allometry), parse = TRUE)
rhodo

