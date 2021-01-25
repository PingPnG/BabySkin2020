#############################################################
# Load the tidyverse and readxl
#https://www.datanovia.com/en/blog/types-of-clustering-methods-overview-and-quick-start-r-code/
#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
#https://www.datacamp.com/community/tutorials/pca-analysis-r
###################################################################
rm(list = ls())
library(tidyverse)
library(readxl)

# Install and load the "here" package
#install.packages('here')
library(here)
#setwd("C:/Users/hu.p/OneDrive - Procter and Gamble/Desktop/Holomonas/8_1_2020")
# Try to read in data with defaults
#test1 <- read_excel(here("Data", "halo_Prokka_gene.xlsx"), sheet = "Sheet3")
test1 <- read_excel(here("Data", "Halomonas_KO.xlsx"), sheet = "Sheet2")
head(test1)
P<-data.frame(test1[,2:199])

rownames(P)=unlist(test1[,1])
colnames(P)=colnames(test1)[2:199]
#summary(P)
strains=colnames(P)
pca_P <- prcomp(P, scale. = TRUE)
groups<-c(rep("Halo", 185), rep("PG", 13))
summary(pca_P)
gcolor<-c(rep("Blue", 185), rep("red", 13))
gcolor[strains=="H_elongata_DSM_2581_type_strain:_DSM_2581.GCF_000196875.2"]="pink"
gcolor[strains=="H_sp_GFAJ-1.GCF_000236625.2"]="orange"
gcolor[strains=="H_subglaciescola_ACAM_12.GCF_900142895.1"]="black"

#install.packages("ggfortify")
library(ggfortify)
P2=as.data.frame(t(P))
dim(P2)
pca_res <- prcomp(P2, scale. = TRUE)
jpeg("Halo_KO_PCA.jpg", width = 350, height = 350)
autoplot(pca_res, colour = gcolor)
dev.off()

jpeg("Halo_KO_PCA_AllwName.jpg", width = 1000, height = 1000)
autoplot(pca_res, colour = gcolor, shape = FALSE, label.size = 3)
dev.off()
####outlier: H. sp MG34 and H.utahensis 
P2<-P2[-c(131,173,174),]
P2<-P2[, colSums(P2 != 0) > 0]
pca_res <- prcomp(P2)
gcolor<-c(rep("Blue", 182), rep("red", 13))
strains=rownames(P2)
gcolor[strains=="H_elongata_DSM_2581_type_strain:_DSM_2581.GCF_000196875.2"]="pink"
gcolor[strains=="H_sp_GFAJ-1.GCF_000236625.2"]="orange"
gcolor[strains=="H_subglaciescola_ACAM_12.GCF_900142895.1"]="black"

jpeg("Halo_KO_PCA_Clean3wName.jpg", width = 1000, height = 1000)
autoplot(pca_res, colour = gcolor, shape = FALSE, label.size = 3)
dev.off()
###########A good and clean cluster formed#########################
jpeg("Halo_KO_PCA_Clean3.jpg", width = 350, height = 350)
autoplot(pca_res, colour = gcolor)
dev.off()
#########

jpeg("Halo_KO_PCA_Clean3_biplot.jpg", width = 1000, height = 1000)
autoplot(pca_res, data=P2,colour = gcolor,
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)
dev.off()
#########################################
####Plotting k-mean, 8 cluster looks very nice seperate P&G strains################
set.seed(1)

jpeg("Halo_KO_PCA_Kmean8.jpg", width = 1000, height = 1000)
autoplot(kmeans(P2, 8), data = P2, shape = FALSE, label.size = 3, frame=TRUE)
dev.off()

jpeg("Halo_KO_PCA_Kmean8_frame.jpg", width = 1000, height = 1000)
autoplot(kmeans(P2, 8), data = P2, shape = FALSE, label.size = 3, frame=TRUE)
dev.off()

jpeg("Halo_KO_PCA_Kmean8_frame_noname.jpg", width = 1000, height = 1000)
autoplot(kmeans(P2, 8), data = P2,  frame=TRUE)
dev.off()

######################################
#library(cluster)

#autoplot(clara(P2, 8))
#autoplot(fanny(P2, 3), frame = TRUE) ### this method not work for 8 clusters####always show 3
#autoplot(pam(iris[-5], 3), frame = TRUE, frame.type = 'norm') ###oval frame

##################method not work so far. I think we can use machine learning algorithm to test hypothesis
###########about which feature is the most predictive about the contamination####
#####need to come back to this.
#install.packages("lfda")
#library(lfda)
#groups<-as.factor(c(rep("PubHalo", 182), rep("PGHalo", 13)))
# Local Fisher Discriminant Analysis (LFDA)
#model <- lfda(P2, groups, r = 2, metric="plain")
#autoplot(model, data = P2, frame = TRUE, frame.colour = gcolor)
# Kernel Local Fisher Discriminant Analysis (KLFDA)
#model <- klfda(kmatrixGauss(iris[-5]), iris[, 5], 4, metric="plain")
#autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')
# Semi-supervised Local Fisher Discriminant Analysis (SELF)
#model <- self(iris[-5], iris[, 5], beta = 0.1, r = 3, metric="plain")
#autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')

###########################################################################################
library("cluster")
library("factoextra")
library("magrittr")

res.dist <- get_dist(P2, stand = TRUE)
jpeg("Halo_KO_eu_dist.jpg", width = 3500, height = 3500)
fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
dev.off()

#######################

#####determine optimal cluster
jpeg("rplot1.jpg", width = 350, height = 350)
fviz_nbclust(P2, kmeans, method = "gap_stat")
dev.off()
#####did not run anything below8-5-2020
set.seed(123)
km.res <- kmeans(P2, 3, nstart = 25)
jpeg("rplot2.jpg", width = 350, height = 350)
fviz_cluster(km.res, data = P2,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())
dev.off()
##################
pam.res <- pam(P2, 3)
jpeg("rplot3.jpg", width = 350, height = 350)
fviz_cluster(pam.res)
dev.off()
#########################
# Compute hierarchical clustering
res.hc <- P2 %>%
  scale() %>%                    # Scale the data
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

# Visualize using factoextra
# Cut in 4 groups and color by groups
jpeg("rplot4.jpg", width = 3500, height = 3500)
fviz_dend(res.hc, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)
dev.off()

################################################
set.seed(123)

# Compute
library("NbClust")
res.nbclust <- P2 %>%
  scale() %>%
  NbClust(distance = "euclidean",
          min.nc = 2, max.nc = 10, 
          method = "complete", index ="all") 
# Visualize
library(factoextra)
jpeg("rplot1.jpg", width = 350, height = 350)
fviz_nbclust(res.nbclust, ggtheme = theme_minimal())
dev.off()











#####have to install from file####
#library(devtools)
#install_github("vqv/ggbiplot")
#install.packages("C:/Users/hu.p/OneDrive - Procter and Gamble/Desktop/Holomonas/8_1_2020/ggbiplot-master.zip", repos = NULL, type = "win.binary")
##########################
library(ggbiplot)
jpeg("rplot.jpg", width = 350, height = 350)
ggbiplot(pca_P, obs.scale = 1, var.scale = 1,
         #groups = groups,
         ellipse = TRUE, circle = TRUE) 
  #+scale_color_discrete(name = '') +
  #theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()


#autoplot(pca_res)


autoplot(pca_res, data = iris, colour = 'Species')
autoplot(pca_res, data = iris, colour = 'Species', loadings = TRUE)

# Ward Hierarchical Clustering
#d <- dist(P, method = "euclidean") # distance matrix
#fit <- hclust(d, method="ward")
#plot(fit) # display dendogram
#groups <- cutree(fit, k=5) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
#rect.hclust(fit, k=5, border="red")

library("cluster")
#install.packages("factoextra")
library("factoextra")
library("magrittr")

my_data <- P %>%
  na.omit() %>%          # Remove missing values (NA)
  scale()                # Scale variables
head(my_data, n = 3)
#res.dist <- get_dist(P, stand = TRUE, method = "pearson")
#fviz_dist(res.dist, 
#          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#library("factoextra")
#fviz_nbclust(my_data, kmeans, method = "gap_stat")