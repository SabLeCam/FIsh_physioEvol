# PHYLOGÉNÉTIQUE COMPARATIVE
source(Introduction to the phylogenetic (comparative) method, J. Wintternitz)

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














