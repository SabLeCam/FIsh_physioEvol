# Analyses contrastes phylogénétiques sous R
À ce stade, vous devriez avoir un arbre phylogénétique pour 25 espèces de poissons (arbre_16S_poissons), ainsi qu’un tableau de données avec, pour chaque espèce, des valeurs de longévité, de température et de nombre d’oeufs par portée et/ou taille des oeufs. La suite de l’analyse va consister à estimer le degré de corrélation entre ces traits en considérant l’impact dû aux relations phylogénétiques entre les espèces. 
Nous utilisons pour cela des modèles généralisés (GLS pour Generalized Least Squares), pour lesquels on va spécifier la structure de covariance due à la phylogénie dans les résidus du modèle. Pour estimer cette structure de covariance, on utilise un modèle appelé « Lambda de Pagel », qui va tenter de quantifier le degré de signal phylogénétique dans les données (ici les résidus). Une valeur de 1 correspond à un signal fort (les espèces proches sont significativement plus semblables que des espèces éloignées), et une valeur de 0 indique qu’il n’y a aucune structure phylogénétique. Et les valeurs intermédiaires illustreront différents niveaux de signal brownien. 
Je vous présente ici les différentes étapes de l’analyse sur R avec les commandes à réaliser.
Tout d’abord charger les 2 packages nécessaires à l’analyse :
```r
library(ape)
library(nlme)
```
Charger l’arbre phylogénétique :
```r
tree<-read.tree("arbre_16S_poissons.tre")

Visualiser l’arbre :
plot(tree);axisPhylo()
Charger les traits biologiques (taille, longévité, température max, nombre d’œufs par portée et/ou taille des oeufs) :
traits_log<-read.table("traits_log_fishes.txt",h=T)
Ylog <- traits_log$Length
Xlog <- traits_log$Longevity
Zlog <- traits_log$Temp_max
Nous avons réalisé au préalable une transformation log des données brutes, en faisant l’hypothèse qu’il existerait une relation allométrique (et donc log linéaire) entre les traits biologiques étudiés. 
On estime ensuite la structure de covariance due à la phylogénie par un modèle de Pagel :
corStr <- corPagel(1, tree)
Et on réalise ensuite un modèle linéaire généralisé en spécifiant la structure de covariance estimé dans l’étape précédente (corStr) :
fit1 <- gls(Ylog ~ Xlog, correlation=corStr)
Lors de l’exécution de cette commande, vous devriez avoir ce message d’erreur :
In Initialize.corPhyl(X[[i]], ...) :
  Rownames in data frame do not match tree tip names; data taken to be in the same order as in tree
C’est parce que le nom des espèces n’est pas indiqué dans le fichier des traits. Par défaut, chaque donnée est attribuée à une espèce selon l’ordre indiqué par l’arbre phylogénétique. Il faudra donc classer les données dans l’ordre des espèces indiqué par l’arbre que vous avez réalisé, chose à faire également si vous êtes amené à tester d’autres jeux de données ou d’autres arbres phylogénétiques. 
À noter aussi qu’on regarde ici la relation entre le log de la taille (Ylog) et le log de la longévité (Xlog). Vous pouvez tout à fait remplacer les deux paramètres étudiés dans la ligne pour tester d’autre trait, comme la température max (Zlog).
Pour afficher un résumé des résultats :
summary(fit1)
Vous devriez voir apparaître ce résultat :
 

On observe tout d’abord un signal phylogénétique plutôt fort (lambda d’environ 0.755), avec une relation allométrique très forte entre la taille et la longévité (p-value=0). On peut noter les valeurs de la pente de la courbe (0.679) et de l’intercept avec l’axe des ordonnées (1.043), autrement dit l’équation de la courbe linéaire (ici Ylog = 0.679 Xlog + 1.043).

Pour réaliser une figure avec la ligne de régression :
plot(Xlog,Ylog)
abline(fit1, col="red")

 
On voit bien ici la corrélation forte entre taille et longévité, légèrement gommé lorsqu’on considère la phylogénie des espèces étudiées. Pour illustrer cela, on peut également réaliser un modèle sans considéré la structure de covariance due à la phylogénie :
fit1_ols <- lm(Ylog ~ Xlog)
Et ajouter la ligne de régression de ce modèle au graphe réalisé juste avant :
abline(fit1_ols, col="green")
 
La pente de la courbe verte est effectivement plus importante que celle de la courbe rouge, ce qui nous indique que la phylogénie intervient dans la relation allométrique positive que l’on observe entre taille et longévité chez les poissons. Mais la phylogénie n’explique pas à elle seule cette relation puisque l’on observe tout de même une corrélation lorsqu’on la considère dans notre modèle. D’autres mécanismes devraient donc intervenir ici pour expliquer cette relation (lien taille-métabolisme vs métabolisme-longévité ?)
Pour finir, voici d’autres commandes que vous pourriez être amené à utiliser selon les analyses que vous effectuerez. 
Si vous souhaitez recréer un arbre qu’avec certaines espèces :
species<-c("Anguilla_rostrata","Cyprinus_carpio","Carassius_auratus "…) 
pour spécifier le nom des espèces entre guillemets, séparés par une virgule. Puis recréer l’arbre :
tree$tip.label[tree$tip.label%in%species]
tree2 <- drop.tip(tree, tip=tree$tip.label[!tree$tip.label%in%species])
Pour afficher le nouvel arbre :
plot(tree2);axisPhylo()

Vous pouvez également réaliser ce genre d’analyse en étudiant les relations existantes entre un trait et les résidus du modèle réalisé pour les deux autres. Exemple ici pour étudier la relation entre température maximale et les résidus du modèle taille-longévité :
fit1_ols <-lm(Xlog ~ Ylog)
res_ fit1_ols <-resid(fit1_ols)
fit2<-gls(Zlog~ res_ fit1_ols,correlation=corStr)
summary(fit3)
plot(res_longevity,Zlog)
abline(fit3,col="red")
fit4<-lm(Zlog~ res_ fit1_ols)
summary(fit4)
abline(fit4,col="green")
