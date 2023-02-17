# PHYLOGÉNÉTIQUE COMPARATIVE

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
babette



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
seq_align <- msaConvert(Aln, type="ape::DNAbin")
```

## Construire des phylogénies à partir des données moléculaires

Calculer une matrice de distance basée sur le nombre de différences nucléotidiques
```r
?dist.dna()
distxj<-dist.dna(seq_align)

#calculates the tree by neighbor joining
tree <- nj(distxj)
#plot a basic tree
plot.phylo(tree, type="phylogram")
```
Calculer une matrice de distance basée sur le nombre de différences nucléotidiques





