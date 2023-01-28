######## Smart mobility in Milan, Italy. A cluster analysis at district level.
# October 2022
# Italian Journal of Applied Statistics

### Caricamento dati e pacchetti utili all'analisi
library(readr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(sf)

'%notin%' <- Negate('%in%')

setwd("C:/Users/paulm/OneDrive/Documenti/GitHub/NC_MS_PM_PMC_MilanoMobility")
Ind_grezzi <- read_excel("MobMi3.xlsx")
Ind_grezzi <- Ind_grezzi %>%
  rename(IdNIL = IdNil)

Index <- read_excel("Index_new.xlsx")

# colnames(Ind_grezzi)[1] = 'IdNIL'
Index <- inner_join(Index, Ind_grezzi, by="IdNIL")
Index_original <- Index
# Escludiamo i quartieri poco abitati
Index <- Index %>%
  filter(IdNIL %notin% c(3,8,39,47,75,86,87)) %>%
  select(-c(Nil))
write_csv(x = Index_original, file = "NILS_DataCluster.csv", col_names = T)

##############################################################
########## Cluster analisys non gerarchica: K-means ##########
##############################################################

library(cluster)
library(stats)
library(fpc)
library(factoextra) 

IndexClust <- Index %>% 
  column_to_rownames(.,'NIL') %>%     # Assegno il nome alle righe usando la variabile codice 
  scale(.,center=T,scale=T) %>%       # Standardizzo Z = (X - mu)/sigma
  data.frame(.)

##### Matrice delle distanze tra le unità usando JJI
IndexClust_JJI <- IndexClust %>%
  dplyr::select(JJI)
distance_JJI <- get_dist(IndexClust_JJI)
fviz_dist(distance_JJI, gradient = list(low="white", mid="grey", high="black"))

##### Matrice delle distanze tra le unità usando AMPI
IndexClust_AMPI <- IndexClust %>%
  dplyr::select(AMPI)
distance_AMPI <- get_dist(IndexClust_AMPI)
fviz_dist(distance_AMPI, gradient = list(low="white", mid="grey", high="black"))

##### Matrice delle distanze tra le unità usando JJI + AMPI
IndexClust_JJIAMPI <- IndexClust %>%
  dplyr::select(AMPI,JJI)
distance_JJIAMPI <- get_dist(IndexClust_JJIAMPI)
fviz_dist(distance_JJIAMPI, gradient = list(low="white", mid="grey", high="black"))



###################################################################################################
########## Identificazione (statistica) del numero ottimale di gruppi: package 'NbClust' ##########
###################################################################################################

library(factoextra)
library(NbClust)
k <- 4

########################
##### Jevons Index #####
########################

## Calcolo degli indicatori di performance per k-means e Ward
nb_kmeans_JJI <- NbClust(IndexClust_JJI, distance = "euclidean", min.nc = 2, max.nc = 6,
                         method = "kmeans", 
                         index ="all")
## "Regola della maggioranza"
p1 <- fviz_nbclust(IndexClust_JJI, kmeans, nstart = 25,  method = "gap_stat", nboot = 100) +
  labs(subtitle = "Gap statistic method", title = "JJI index")
p2 <- fviz_nbclust(nb_kmeans_JJI,method = "silhouette") + 
  labs(subtitle = "Silhouette", title = "JJI index") + 
  scale_x_discrete(name ="Dose (mg)", limits=factor(c(3,4)))
p3 <- fviz_nbclust(nb_kmeans_JJI,method = "wss") + 
  labs(subtitle = "Majority rule-of-thumb", title = "JJI index") + 
  scale_x_discrete(name ="Dose (mg)", limits=factor(c(3,4)))

## Numero di clusters ottimali per ogni metodo
nb_kmeans_JJI$Best.nc
## "Regola della maggioranza"
fviz_nbclust(nb_kmeans_JJI)
## suddivisione ottimale delle unitÃ 
nb_kmeans_JJI$Best.partition

### Optimal number of groups: 5
clust_k5_JJI <- kmeans(IndexClust_JJI, centers=k, nstart = 50) # 5 gruppi

## Optimal clusters
output_JJI <- cbind(k_kmeans_opt_JJI = clust_k5_JJI$cluster)
output <- cbind(IdNIL = Index$IdNIL, output_JJI)


######################
##### AMPI Index #####
######################

## Calcolo degli indicatori di performance per k-means e Ward
nb_kmeans_AMPI <- NbClust(IndexClust_AMPI, distance = "euclidean", min.nc = 2, max.nc = 6,
                          method = "kmeans", index ="all")

## Numero di clusters ottimali per ogni metodo
nb_kmeans_AMPI$Best.nc
## "Regola della maggioranza"
p4 <- fviz_nbclust(IndexClust_AMPI, kmeans, nstart = 25,  method = "gap_stat", nboot = 100) +
  labs(subtitle = "Gap statistic method", title = "AMPI index")
p5 <- fviz_nbclust(nb_kmeans_AMPI,method = "silhouette") + 
  labs(subtitle = "Silhouette", title = "AMPI index") + 
  scale_x_discrete(name ="Dose (mg)", limits=factor(c(3,4)))
p6 <- fviz_nbclust(nb_kmeans_AMPI,method = "wss") + 
  labs(subtitle = "Majority rule-of-thumb", title = "AMPI index") + 
  scale_x_discrete(name ="Dose (mg)", limits=factor(c(3,4)))

## suddivisione ottimale delle unità
nb_kmeans_AMPI$Best.partition

### Optimal number of groups: 3 --> use 5 for completeness
clust_k5_AMPI <- kmeans(IndexClust_AMPI, centers=k, nstart = 50) # 5 gruppi

## Optimal clusters
output <- cbind(output, k_kmeans_opt_AMPI = clust_k5_AMPI$cluster)

############################
##### JJI + AMPI Index #####
############################

## Calcolo degli indicatori di performance per k-means e Ward
nb_kmeans_JJIAMPI <- NbClust(IndexClust_JJIAMPI, distance = "euclidean", min.nc = 2, max.nc = 6,
                             method = "kmeans", index ="all")
## "Regola della maggioranza"
p7<- fviz_nbclust(IndexClust_JJIAMPI, kmeans, nstart = 25,  method = "gap_stat", nboot = 100) +
  labs(subtitle = "Gap statistic method", title = "JJI+AMPI index")
p8 <- fviz_nbclust(nb_kmeans_JJIAMPI,method = "silhouette") + 
  labs(subtitle = "Silhouette", title = "JJI+AMPI index") + 
  scale_x_discrete(name ="Dose (mg)", limits=factor(c(3,4)))
p9 <- fviz_nbclust(nb_kmeans_JJIAMPI,method = "wss") + 
  labs(subtitle = "Majority rule-of-thumb", title = "JJI+AMPI index") + 
  scale_x_discrete(name ="Dose (mg)", limits=factor(c(3,4)))

## Numero di clusters ottimali per ogni metodo
nb_kmeans_AMPI$Best.nc
## "Regola della maggioranza"
fviz_nbclust(nb_kmeans_JJIAMPI)
## suddivisione ottimale delle unità
nb_kmeans_AMPI$Best.partition

### Optimal number of groups: 3 --> use 5 for completeness
clust_k5_JJIAMPI <- kmeans(IndexClust_JJIAMPI, centers=k, nstart = 50) # 5 gruppi

## Optimal clusters
output <- cbind(output, k_kmeans_opt_JJIAMPI = clust_k5_JJIAMPI$cluster)

png(filename = "Cluster_metrics.png", width = 800, height = 1000)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol = 3,nrow = 3)
dev.off()


############################################
########## Mapping of the results ##########
############################################

## Load shape Milano
library(sf)
Quartieri <- read_sf("Quartieri_Milano/Quartieri_Milano.shp", 
                     quiet = FALSE) %>% st_transform(crs = 4326) %>%
  rename(IdNIL = ID_NIL)

output <- data.frame(output)
Quartieri_cluster <- full_join(Quartieri, output, by="IdNIL")
Quartieri_cluster <- Quartieri_cluster %>%
  mutate(Group_JJI = case_when(k_kmeans_opt_JJI == 3 ~ "Medium-low",
                               k_kmeans_opt_JJI == 4 ~ "High",
                               k_kmeans_opt_JJI == 1 ~ "Low",
                               k_kmeans_opt_JJI == 2 ~ "Medium-high",
                               is.na(k_kmeans_opt_JJI) ~ "Outlier"),
         Group_AMPI = case_when(k_kmeans_opt_AMPI == 3 ~ "Low",
                                k_kmeans_opt_AMPI == 4 ~ "High",
                                k_kmeans_opt_AMPI == 2 ~ "Medium-high",
                                k_kmeans_opt_AMPI == 1 ~ "Medium-low",
                                is.na(k_kmeans_opt_AMPI) ~ "Outlier"),
         Group_JJIAMPI = case_when(k_kmeans_opt_JJIAMPI == 3 ~ "Medium-low",
                                   k_kmeans_opt_JJIAMPI == 1 ~ "Low",
                                   k_kmeans_opt_JJIAMPI == 4 ~ "High",
                                   k_kmeans_opt_JJIAMPI == 2 ~ "Medium-high",
                                   is.na(k_kmeans_opt_JJIAMPI) ~ "Outlier")) %>%
  select(IdNIL, NIL, Group_JJI, Group_AMPI, Group_JJIAMPI)

Quartieri_cluster_df <- Quartieri_cluster
st_geometry(Quartieri_cluster_df) <- NULL

library(xtable)
write_csv(x = Quartieri_cluster_df, file = "NILS_clusters.csv", col_names = T)
Quartieri_cluster_df <- Quartieri_cluster_df %>%
  arrange(IdNIL) %>%
  select(-c(IdNIL))
Cluster_grp <- xtable(Quartieri_cluster_df,
                      caption = paste0("Cluster results using Jevons only, AMPI only and Jevons+AMPI, respectively"), 
                      align=c("l","l","c","c","c"))
Cluster_grp

# Graphical parameters
fills <- c("Outlier" = "white",
           "High" = "#1F3457",
           "Medium-high"="#F1625F",
           "Medium-low"="#9D9A99", 
           "Low"="#CBE8DD")

#Generate the maps for the 3 clustering
#JJI cluster map
Quartieri_cluster$Group_JJI <- factor(Quartieri_cluster$Group_JJI, 
                                      levels = c("High", "Medium-high", "Medium-low",
                                                 "Low", "Outlier"))
JJI_cluster <- ggplot()+
  geom_sf(data = Quartieri_cluster, aes(fill = Group_JJI), show.legend = "point")+
  theme_void()+
  ggtitle('\n1. K-means: Jevons static index (K=4)') +
  theme(plot.title = element_text(size = 20, face = "bold")) +
  scale_fill_manual(values=fills, labels=c("Outlier","High", "Medium-high", "Medium-low", "Low"), name="Clusters:") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5))) +
  theme(legend.position="bottom") +
  theme(legend.title = element_text(face="bold", size=15)) +
  theme(legend.key.size = unit(0.7, "cm")) +
  theme(legend.text=element_text(size=13)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(aspect.ratio=3/4)
JJI_cluster  


#AMPI cluster map
Quartieri_cluster$Group_AMPI <- factor(Quartieri_cluster$Group_AMPI, 
                                       levels = c("High", "Medium-high", "Medium-low",
                                                  "Low", "Outlier"))
AMPI_cluster <- ggplot()+
  geom_sf(data = Quartieri_cluster, aes(fill = Group_AMPI), show.legend = "point")+
  theme_void()+
  ggtitle('\n2. K-means: AMPI index (K=4)') + 
  theme(plot.title = element_text(size = 20, face = "bold")) +
  scale_fill_manual(values=fills, labels=c("Outlier","High", "Medium-high", "Medium-low", "Low"), name="Clusters:") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5))) +
  theme(legend.position="bottom") +
  theme(legend.title = element_text(face="bold", size=15)) +
  theme(legend.key.size = unit(0.7, "cm")) +
  theme(legend.text=element_text(size=13)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(aspect.ratio=3/4) 
AMPI_cluster


#JJI AMPI cluster map
Quartieri_cluster$Group_JJIAMPI <- factor(Quartieri_cluster$Group_JJIAMPI, 
                                          levels = c("High", "Medium-high", "Medium-low",
                                                     "Low", "Outlier"))
JJI_AMPI_cluster <- ggplot()+
  geom_sf(data = Quartieri_cluster, aes(fill = Group_JJIAMPI), show.legend = "point")+
  theme_void()+
  ggtitle('\n3. K-means: JJI+AMPI index (K=4)') + 
  theme(plot.title = element_text(size = 20, face = "bold")) +
  scale_fill_manual(values=fills, labels=c("Outlier","High", "Medium-high", "Medium-low", "Low"), name="Clusters:") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5))) +
  theme(legend.position="bottom") +
  theme(legend.title = element_text(face="bold", size=15)) +
  theme(legend.key.size = unit(0.7, "cm")) +
  theme(legend.text=element_text(size=13)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(aspect.ratio=3/4)
JJI_AMPI_cluster 

#combine the maps together
library(cowplot)
Map_cluster <- plot_grid(
  JJI_cluster, AMPI_cluster, JJI_AMPI_cluster, ncol = 1)
#Save the maps
ggsave(Map_cluster, file= "Map_cluster.png", width = 32, height = 60, units = "cm")





#################################################################################
png("Correlation.png", width = 40, height = 20, res = 800, units = 'cm')
ggplot(data = IndexClust, aes(x=JJI, y=AMPI)) + 
  geom_point(lwd=3) + 
  geom_smooth(method = "lm", se = FALSE, lwd=1.5) + 
  geom_text(mapping = aes(x = 3, y=5), size=7,label = "Linear correlation", col="#000099") +
  geom_text(mapping = aes(x = 3, y=4.5), size=7,
            label = expression(paste(rho,"=0.4615")), col="#000099") +
  labs(title = "Adjusted Mazziota-Pareto Index VS Jevons Static Index",
       subtitle = "Estimated values and linear correlation",
       x="Standardized JJI", y="Standardized AMPI") + 
  theme(
    axis.title.y=element_text(color="black", size=14),
    axis.text.y=element_text(color="black", size=12),
    axis.title.x=element_text(color="black", size=14),
    axis.text.x=element_text(color="black"),
    title = element_text(size=16))
dev.off()