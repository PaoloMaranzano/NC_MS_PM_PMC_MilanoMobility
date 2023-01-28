library(readr)
library(tidyverse)
library(sf)
library(spatialEco)

setwd("C:/Users/paulm/OneDrive/Documenti/GitHub/NC_MS_PM_PMC_MilanoMobility")

##### Dataset used for K-means clustering
DataClusters <- read_csv("NILS_DataCluster.csv")
DataClusters <- DataClusters %>%
  select(IdNIL,JJI,Media01,IMG,AMPI)
##### K-means output
Clusters <- read_csv("NILS_clusters.csv")
Clusters <- Clusters %>%
  select(-c(NIL))
##### Services data
Services <- read_csv("Services.csv", col_types = cols(...1 = col_skip()))

##### Join
Data <- full_join(x = Services, y = Clusters, by = "IdNIL")
Data <- full_join(x = Data, y = DataClusters, by = "IdNIL")

##### Shape file
Quartieri <- read_sf("Quartieri_Milano/Quartieri_Milano.shp", 
                     quiet = FALSE) %>% 
  st_transform(crs = 4326) %>%
  select(IdNIL = ID_NIL)

Data <- inner_join(Quartieri, Data, by="IdNIL")

##### Extract coordinates and Euclidean distance
# Centr_coords <- Data %>%
#   st_geometry() %>%
#   st_centroid() %>%
#   st_coordinates()
Data_point <- st_point_on_surface(x = Data)
Centr_coords <- st_coordinates(Data_point)
mtx_distance <- st_distance(Data_point, Data_point)
mtx_distance_pearson <- matrix(data = 1/dim(Centr_coords)[1],
                               nrow = nrow(mtx_distance),ncol = ncol(mtx_distance))

##### Spatial cross-correlation
X_names <- c("Area","Population","Pop_density","N_Hospitals","N_Hospital_beds",
             "N_Stations","N_Uni_campus","N_Uni_headquarters","N_Schools","N_Schools_students",
             "N_Hotels","N_Hotels_rooms","N_Hotels_beds")
Y_names <- c("AMPI","JJI","Media01","IMG")
Cor_list <- vector("list", length = length(Y_names))
Cor_mat <- matrix(data = NA, ncol = 5, nrow = length(X_names))
rownames(Cor_mat) <- X_names
colnames(Cor_mat) <- c("Index","Variable","Rc","Rp","R0")
Cor_mat <- data.frame(Cor_mat)
for (i in 1:length(Y_names)) {
  for (j in 1:length(X_names)) {
    print(paste0("Y = ",Y_names[i]," e X = ",X_names[j]))
    r_temp <- crossCorrelation(x = Data_point[[X_names[j]]], y = Data_point[[Y_names[i]]],
                               type =  "LSCI",coords = Centr_coords,
                               k = 999, scale.partial = T, scale.matrix = T,
                               dist.function = "inv.power",  clust = T)
    x <- as.vector(scale(Data_point[[X_names[j]]],center = T,scale = T))
    y <- as.vector(scale(Data_point[[Y_names[i]]],center = T,scale = T))
    R0 <- cor(x,y)
    Rc <- r_temp$I
    Rp <- R0 - Rc
    Cor_mat[j,] <- cbind(Y_names[i],X_names[j],Rc,Rp,R0)
  }
  Cor_list[[i]] <- Cor_mat
}
SpatCorr <- bind_rows(Cor_list)
rownames(SpatCorr) <- NULL
SpatCorr <- SpatCorr %>%
  mutate(across(c(Rc,Rp,R0),as.numeric),
         across(c(Rc,Rp,R0),~ round(.x,digits = 3))) %>%
  arrange(R0)
SpatCorr %>%
  View()


library(xtable)
SpatCorr_latex <- xtable(SpatCorr,
                         caption = paste0("Spatial cross-correlation results."), 
                         align=c("l","l","l","c","c","c"))
SpatCorr_latex



x <- as.vector(scale(Data_point[["N_Hotels"]],center = T,scale = T))
y <- as.vector(scale(Data_point[["JJI"]],center = T,scale = T))
# ## 2. Filtro di Hampel
# soglia <- 1.4826*3*median(abs((x - median(x))))
# # Posizione dei valori anomali
# pos_out_Hampel <- which(abs(x - median(x)) > soglia)
spcor_hotel <- crossCorrelation(x = x, y = y,
                           type =  "LSCI",coords = Centr_coords,
                           k = 999, scale.partial = T, scale.matrix = T,
                           dist.function = "inv.power",  clust = T)
plot(spcor_hotel$SCI[,1],spcor_hotel$SCI[,2])
abline(h = 0,col="red", lwd=2)
abline(v = 0, col="red", lwd=2)
