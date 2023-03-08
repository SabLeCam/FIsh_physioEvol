# PHYLOGÉNÉTIQUE COMPARATIVE
[source](https://static1.squarespace.com/static/5459da8ae4b042d9849b7a7b/t/57ea64eae58c62718aa34769/1474979059782/Nesin_Winternitz_Practical_1and2.pdf)(Introduction to the phylogenetic (comparative) method, J. Wintternitz 2016)

Ici vous trouverez comment réaliser l'enesemble des analyses de phylogénétique comparative sous R

## package R nécessaires pour ces analyses

```r
adegenet
ape
geiger
phytools
picante
stringr
msa
caper
```
installer les packages necessaires et les appeler

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("msa")

install.packages("caper") ##continue with following package names
library(caper) ##continue with following packages

adegenet
ape
geiger
phytools
picante
stringr
msa
treetools




```

## Importer des données

1. les données phénotypiques

sauver les données des traits phénotypiques en txt ou csv au préalable

```r
#read data as comma-delimitated (.csv) file with headers on columns
traits<-read.csv("traits.csv",header=TRUE)
#or read data as tab-delimitated (.txt) files
traits<- read.delim("traits.txt", header=TRUE)
#check if your data are properly loaded
traits
str(traits)
```

2. les données moléculaires

```r
seq <-read.dna("fish_16S.fasta", format="fasta")
##check your sequences
image(seq)

path<-'PATH_TO_YOUR_FASTA_FILE'

seq <-Biostrings::readDNAStringSet(path)# load file
seq

```

## Aligner les sequences

```r
Aln <- msa(seq, method = "ClustalOmega", order="aligned")
#export alignment as DNAbin object
seq_align <- msaConvert(Aln, type="seqinr::alignment")
```

## Construire des phylogénies à partir des données moléculaires

Calculer une matrice de distance basée sur le nombre de différences nucléotidiques
```r
?dist.dna()###regarder les différents modèles d'évolution moléculaire disponibles
distxj<-dist.dna(seq_align)
```


# Construire un arbre avec la méthode du neighbor-joining
```r
tree <- nj(distxj)
#plot a basic tree
plot.phylo(tree, type="phylogram")
```
![image](https://user-images.githubusercontent.com/20643860/219545168-19ba9230-fc18-4eab-aa21-e6464b8dba6f.png)

Exemples de représentations de l'arbre
```r
par(mfrow=c(2,2)) #plot 2 rows by 2 columns
plot.phylo(tree, type="phylogram")
plot.phylo(x=tree, type="cladogram", edge.width=2)
plot.phylo(x=tree, type="fan", edge.width=2, edge.lty=2, cex=.8)
plot.phylo(x=tree, type="radial", edge.color="red", edge.width=2, edge.lty=3, cex=.8)
```
![image](https://user-images.githubusercontent.com/20643860/219545266-2ea7772e-f342-40ce-8913-8e4dd23d133e.png)

Bootstrap pour evaluer la robustesse des noeuds

```r
TipLabels(tree) # connaitre la position de notre outgroup
outgroup <- 2  
#fonction pour enraciner l'arbre:
  foo <- function(xx) root(nj(dist.dna(xx)), outgroup)
tr <- foo(seq_align) 
bp <- boot.phylo(tr, seq_align, foo, B=1000) #1000 bootstraps

#représentation graphique de l'arbre avec les valeur de boostrap pour les noeuds
  plot(tr)
nodelabels(round(bp/10))
```
![image](https://user-images.githubusercontent.com/20643860/219698076-aa703587-2dc4-4748-b337-de87759cd0ba.png)



Représenter les familles plutôt que les noms scientifiques et utiliser des symboles

```r
plot.phylo(x=tr, type="cladogram", show.tip=FALSE, lwd=3, main="Neighbour-Joining tree")
#add axis with distances
axisPhylo()
#use symbols (pch) for tip labels
tiplabels(frame="none", pch=rep(x=0:7, times=c(2, 1,1,1,1,1,1,1)), lwd=2, cex=2)
 legend(x="topleft", legend=c("Acipenseridae","Rajidae","Cyprinidae","Carcharhinidae","Gadidae", 
 "Percidae", "Gasterosteidae", "Pleuronectidae"), border="black", pch=0:7, pt.lwd=2, pt.cex=1.5, bty="o", bg="lightgrey", box.lwd=1, cex=1.2, title="Famille")
 ```
![image](https://user-images.githubusercontent.com/20643860/219698249-4b62de44-d42b-459c-bc00-e0ee93f1602b.png)
 
 ##Représenter les traits sur l'arbre phylogénétique
 
pour un trait continu
 
 ```r
rownames(traits)<-traits$BINOMIAL.NAME
body.size<-as.matrix(traits)[,6]
mode(body.size)<-'numeric' 

body.size<-body.size[tree$tip.label]#on assignt le trait aux tiplabels

#on peut construire l'arbre
obj<-contMap(tree,body.size,plot=FALSE)
plot(obj,type="phylogram",leg.txt="Max length",lwd=6, mar=c(4,2,4,2))
title(main="fish phylogenetic tree")
axis(1)
title(xlab="Time from the root")
```
![image](https://user-images.githubusercontent.com/20643860/219698443-e9b3116c-776f-432a-8a2a-0ce0aedc3640.png)


pour un trait discret

```r
#create a discrete trait vector

traits$habitat<-c("Freshwater","Both","Marine","Freshwater","Marine","Both","Marine","Freshwater","Both","Marine")

habitat<-as.data.frame(traits$habitat)
rownames(habitat)<-tree$tip.label
x<-as.matrix(habitat)[,1] 



tree<-tr #rooted tree
plot(tree,type="phylogram", cex=0.9, lwd=2, label.offset=.002)
title("Habitat")


cols<-setNames(c("White","Green","Blue"), sort(unique(x)))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.5)
#prompt=TRUE means you click to draw the legend
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8,
                  vertical=TRUE, shape="square")
                  
```
![image](https://user-images.githubusercontent.com/20643860/219698627-3c4bd929-0c9d-4bb1-8c34-452399f18912.png)

On peut aussi reconstruire l'état ancestral du trait , méthode de maximum de vraisemblance

```r
#vérifier que l'arbre est bien dichotomique
tree<-multi2di(tree)
is.binary(tree) #on veut un TRUE

#l'arbre doit être enraciné
dtree<-tree
dtree$edge.length[dtree$edge.length==0]<-max(nodeHeights(tree))*1e-6

ans<-ace(x,dtree, model="ER", type="discrete")
ans

```
<img width="575" alt="image" src="https://user-images.githubusercontent.com/20643860/219701999-8c85f086-1213-4d50-83af-80104024f236.png">
Dans ce cas l'estimation n'est pas très bonne car on a une vraisemblance équivalente du caractère à la racine.
Si on regarde ces valeurs par noeud, le constat est le même.

```r
ans$lik.anc
```
La cause probable est que le caractère n'est pas bien ségrégé par la phylogénie, on retrouve les différents états du caractère dans chaque clade.

autre méthode: inférence bayésienne
```r
#SImulation de 300 carte stochastique du caractère
  trees<-make.simmap(tree, x ,model="ER", nsim=300)
#obtenir l'état ancestral par nooeud
aa<-describe.simmap(trees,plot=FALSE)
aa
```

visualisation
```r
plotSimmap(trees[[1]],cols,pts=FALSE, lwd=3, setEnv=TRUE, offset=.5)
nodelabels(pie = aa$ace, piecol = cols, cex = 0.6)
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
#prompt=TRUE means you click to draw the legend
add.simmap.legend(colors=cols,prompt=TRUE,fsize=0.8,vertical=TRUE, shape="square")
```

![image](https://user-images.githubusercontent.com/20643860/219701552-8f3619fe-daa8-4a16-b1a5-2e6b550a734c.png)


# Tester des modèles évolutifs en prenant en compte la phylogénie

Vérifier que l'arbre est dichotomique (pleinement résolu):
```r
is.binary.tree(tree)
```
Vérifier que l'arbre est enraciné et que la racine a une longueur non nulle:
```r
is.rooted(tree) #check if tree is rooted
#IF FALSE
tree<-root(tree, outgroup, node = NULL, resolve.root = TRUE)
is.rooted(tree) #check if tree is rooted
```
Définir les branche de longueur null à 1/10 000e de la taille de l'arbre
```r
tree$edge.length[tree$edge.length==0]<-max(nodeHeights(tree))*1e-4
```
Vérifiez que les noms des feuilles correspondent EXACTEMENT à ceux du dataset
```r
name.check(tree, dataset)
#if necessary
tree$tip.label<-gsub(" ", "_", tree$tip.label) #replace any whitespace with underscore
rownames(dataset)<-dataset$Species #name the data rows
```
Regardons les données:
```r
par(mfrow=c(1,2))
hist(traits$AM)
hist(traits$MAX.AGE)
 dev.off()
 ```
 
 On peut  normaliser les données avec une tranformation log10, afin de se rapprocher de l'hypothèse de normalité
 
 ```r
 par(mfrow=c(1,1))
hist(log10(traits$AM))
traits$log.AM<-log10(traits$AM)
```

On peut maintenant réaliser une régression linéaire pour voir la corrélation entre ces traits
```r
model1<-lm(MAX.AGE~log.AM, data=traits)
summary(model1)
```
<img width="461" alt="image" src="https://user-images.githubusercontent.com/20643860/223836405-de93b103-0a5f-4e94-a293-5cba95847652.png">

```r
plot(MAX.AGE~log.AM, data=traits)
abline(model1, col="Red")
```
![image](https://user-images.githubusercontent.com/20643860/223836812-17c46e6b-79cf-4dd2-ac5d-397401471d13.png)

Il y a une corrélation significative entre ces traits mais à ce stade on ne sait pas si cette relation est indépendante de la phylogénie.

**1) Calculer différentes mesures du signal phylogénétique**

<img width="411" alt="image" src="https://user-images.githubusercontent.com/20643860/223840055-d2043619-214a-4f2d-bca1-425cd05939e4.png">
Figure issue de *"Meireles, José Eduardo, Brian O’Meara, and Jeannine Cavender-Bares. "Linking leaf spectra to the plant tree of life." Remote sensing of plant biodiversity (2020): 155-172.*


- <ins>Signal phylogénétique avec le λ (lambda) de Pagel</ins> (Pagel, 1999)

Lambda est une facteur de transformation des branches de l'arbre afin que les données obtenues 'fit' une modèle de mouvement brownien (i.e. dérive génétique). Si le lambda estimé = 0, les traits inférés ont évolué indépendamment de la phylogénie. Lambda = 1 correspond à un modèle brownien (l'évolution du trait est entièrement dû un signal génétique). Lambda varie généralement entre 0 et 1 mais Lambda peut être supérieur à 1 lorsque les covariances phylogénétiques sont plus fortes qu'attendes sous le modèle brownien.

```r
tree0 <- rescale(tree, model = "lambda", 0)
tree1 <- rescale(tree, model = "lambda", 1)
tree.3 <- rescale(tree, model = "lambda", 0.3)
par(mfcol = c(1, 4))
plot(tree)
title("original tree")
plot(tree1) #complete brownian motion
title("lambda=1")
plot(tree.3) #some weak signal
title("lambda=0.3")
plot(tree0) #star phylogeny (no phylo signal)
title("lambda=0")
```
![image](https://user-images.githubusercontent.com/20643860/223847559-c8990344-1df7-4b60-879f-c42c9a85fe79.png)

La fonction *phylosig* (package *phytools*) permet de tester si lambda est significativement différent de 0 pour les deux traits
```r
logAM <- traits$log.AM
names(logAM) <- rownames(traits)
phylosig(tree, logAM, method = "lambda", test = TRUE)
```
<img width="285" alt="image" src="https://user-images.githubusercontent.com/20643860/223848628-8f41a6d5-e168-4f44-9ba9-1853567a66a0.png">

Il semble y avoir un signal phylogénétique derrière l'evolution de l'âge à la maturité. c'est corroborré par les deux méthode d'estimation de lambda.

On peut aussi estimer lambda avec les fonctions *fitDiscrete* ou *fitContinuous* du package *geiger*, selon si on a des traits continus ou discrets. On génère l'estimation de lambda avec le maximum de vraisemblance.

```r
lambda.AM<-fitContinuous(phy = tree, dat = logAM, model = "lambda")
```
<img width="387" alt="image" src="https://user-images.githubusercontent.com/20643860/223850011-0532d2a7-c53f-4af4-9bfd-ab16866efd5f.png">

On procède de même avec la longévité:
```r
MAXAGE <- traits$MAX.AGE
names(MAXAGE) <- rownames(traits)
phylosig(tree, MAXAGE, method = "lambda", test = TRUE)

lambda.MAXAGE<-fitContinuous(phy = tree, dat = MAXAGE, model = "lambda")
```
Il n'y a pas de signal phylogénétique pour la longévité


- <ins>Signal phylogénétique avec le K de Blomberg</ins> (Blomberg et al. 2003)

K peut être considéré comme la proportion de covariance dûe à la phylogénie. Varie généralement entre O et 1 mais peut être supérieur à 1 si les covariances phylogénétiques sont plus fortes qu'attendues sous le modèle brownien.
```r
phylosig(tree, logAM, method = "K", test = T)
```
<img width="344" alt="image" src="https://user-images.githubusercontent.com/20643860/223858933-797fc79e-5e09-4eea-9a7f-85b04e5a02a7.png">
```r
phylosig(tree, MAXAGE, method = "K", test = T)
```
<img width="344" alt="image" src="https://user-images.githubusercontent.com/20643860/223859443-b105555c-7c54-4421-9f39-a35d800d373b.png">



**2) Analyses des contrastes indépendants**

les contrastes sont les différences de valeurs de traits entre espèces (et noeuds). On utilise la méthode des contrastes indépendant afin d'estimer les corrélations évolutives entre les caractères.

```r
#regrouper les données phénotypiques et phylogénétiques en un objet

comp<- comparative.data(tree, traits, names.col=BINOMIAL.NAME, vcv=TRUE,
                        na.omit=FALSE, warn.dropped=TRUE)
#avec la fonction crunch(), on fitte les modèles de régression des contrats indépendants
comp.ic<- crunch(MAX.AGE ~ log.AM, data=comp, equal.branch.length=TRUE)
summary(comp.ic)
```
<img width="462" alt="image" src="https://user-images.githubusercontent.com/20643860/223866208-32ce305e-4204-433f-a014-d04674331dc1.png">


On visualise les résultats et les compare avec les corrélations ne tenant pas compte du signal phylogénétique
```r
comp_contr<-caic.table(comp.ic)

par(mfrow=c(1,2))
plot(MAX.AGE ~ log.AM , data=comp_contr)
abline(comp.ic, col="red")
title("IC regression")
plot(MAX.AGE ~ log.AM, data=traits)
abline(model1, col="red")
title("original regression")
```
![image](https://user-images.githubusercontent.com/20643860/223866358-34b71fe1-2fc8-495c-bc57-3daa075d9b16.png)


**3) Réaliser des analyses PGLS (phylogenetic generalized least squares)**

Cette méthode est plus flexible que celle des contrastes indépendants. Le modèle d'évolution peut s'écarter du modèle brownien strcit (Lambda ou K=1). Plusieurs paramètres d'ajustement peuvent être incoporé dans l'analyse (voir help du package *caper* pour des options) afin d'améliorer l'estimation des corrélations entre les traits. De plus, l'intercept de la régression ne doit pas être fixée à 0.

```r
comp.pgls<-pgls(MAX.AGE ~ log.AM, data=comp, lambda="ML")
summary(comp.pgls)
```
<img width="459" alt="image" src="https://user-images.githubusercontent.com/20643860/223879490-b8b382bf-38cd-4121-ba60-63c880b72021.png">
La corrélation reste significative mais les paramètres de l'équation sont différents (comparaison avec modèle 1)

```r
plot(MAX.AGE ~ log.AM , data=traits)
abline(comp.pgls, col="red",lty=1)
abline(model1, col="red",lty=2)
legend()
legend('topleft', legend=c("PGLS", "Original"), cex=0.8, col="red", lty=c(1,2))
```

![image](https://user-images.githubusercontent.com/20643860/223879323-7df79950-3597-403f-bb2d-60f152c857a8.png)


























