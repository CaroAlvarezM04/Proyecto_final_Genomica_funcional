#Librerias
#Esophagus es una base de datos precargada en phyloseq
library("phyloseq")
#Es una libreria para hacer mas amigables las graficas y plots generados
library("ggplot2")

#Conjunto de temas del paquete ggplot2
theme_set(theme_bw())

#Para cargar datos de ejemplo en el entorno de trabajo
data(esophagus)

#Los datos engloban ciertas tablas, los siguientes comandos son para familiarizarnos y explorar dichos datos
View(esophagus)
#Cada renglón es una taxa y cada columna es una muestra
View(otu_table(esophagus))
#Guardar la tabla en un nuevo objeto para poder manipularlo
otu <- otu_table(esophagus)

#Analsis de riqueza 
riqueza <- plot_richness(otu)
riqueza
riqueza2 <- riqueza + geom_boxplot(data=riqueza$data, aes(color="samples"), alpha=0.1)
riqueza2

#Filogenia que agrupa los OTUS presentes en las tres muestras de acuerdo a la abundancia
arbol <- plot_tree(esophagus, shape="OTU", size="Abundance")
arbol
#Cargar la matriz
matriz <- read.csv("matriz.csv")
matriz

#Matriz de presencia/ ausencia
arbol_matriz <- read.csv("arbol_matriz.csv")
arbol_matriz

#calculo de la matriz de distancia
dist_mat <- dist(arbol_matriz, method = "euclidean")
#clusterizacion
model <- hclust(dist_mat,"average")
model
#dendograma a partir del metodo de clusterizacion
dendograma <- as.dendrogram(model)
plot(dendograma)

matriz_abun <- read.csv("matriz_abun.csv")
head(matriz_abun)
# generate number of colours equal to number of phyla
install.packages("viridis")
library("viridis")
colours <- viridis::viridis(length(unique(matriz_abun$Genero)))
color.indices <- match(matriz_abun$Genero, unique(matriz_abun$Genero))
# generate a vector with color for each taxon
colvec <- colours[color.indices]

generos <- barplot(as.matrix(matriz_abun[,-1]), col=colvec)+
legend('bottomright', fill=colours,legend = matriz_abun$Genero)

generos2 <- heatmap(as.matrix(matriz_abun[,-1]), col=colvec)+
  legend('bottomright', fill=colours,legend = matriz_abun$Genero)


#_________________________________________

# Cómo volver a crear el conjunto de datos del esófago usando la función import_mothur
mothlist <- system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
mothgroup <- system.file("extdata", "esophagus.good.groups.gz", package="phyloseq")
mothtree <- system.file("extdata", "esophagus.tree.gz", package="phyloseq")
mothtree
show_mothur_cutoffs(mothlist)
cutoff <- "0.10"
esophman <- import_mothur(mothlist, mothgroup, mothtree, cutoff)

esophman

#________

taxotutable <- phyloseq::psmelt(esophagus)
taxotutable
write.csv(taxotutable, file = "OTU.csv", sep = ",")


OTU <- read.csv("OTU.csv")
#hacer una matriz de adyacencia
library("igraph")
ad_mat <- graph_from_adjacency_matrix(as.matrix(OTU[-1,-1]))
plot(ad_mat)

BiocManager::install("RCy3")
library(RCy3)
createNetworkFromIgraph(ad_mat)