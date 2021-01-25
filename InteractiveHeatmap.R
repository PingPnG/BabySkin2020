rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]
##########################
#install.packages("heatmaply")
#library("heatmaply")

filename="keratin.data"
A<-read.table(filename, sep="\t", header=TRUE)
d<-dim(A)
B=A[1:d[1], 5:d[2]]
rownames(B)= paste0(A[,2], ".", A[,1])
Cname=colnames(A)[5:d[2]]

Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[.]")
Group1=rep(NA, Clen)
Group2=rep(NA,Clen)
Group3=rep(NA, Clen)
GSSID=rep(NA, Clen)
for(mm in  1:Clen ){
       Group1[mm]=splitname[[mm]][1]
       Group2[mm]=splitname[[mm]][2]
       Group3[mm]=splitname[[mm]][3]
      GSSID[mm]=splitname[[mm]][4]
}

C=B[rowMeans(B)>100, Group1=="Baby"]
#Group1 =="Baby"
D=C[rowMeans(C)>100,]
library("heatmaply")
df <- normalize(C)
heatmaply(df)

###Generate static heatmap
p1<- ggheatmap(B)
p1

p2<-gplots::heatmap.2(
  as.matrix(df),
  trace = "none",
  col = viridis(100),
  key = FALSE
)


p3<-heatmaply(
  as.matrix(df),
  seriate = "mean", 
  row_dend_left = TRUE,
  plot_method = "plotly"
)

p3

p4<-heatmaply(
  df,
  k_col = 2, 
  k_row = 2
)

p4

p5<-heatmaply(
  df,
  colors = viridis(n = 256,  option = "magma"),
  k_col = 2, 
  k_row = 2
)

p5
library(RColorBrewer)
p6<-heatmaply(
  df,
  colors = colorRampPalette(brewer.pal(3, "RdBu"))(256),
  k_col = 2, 
  k_row = 2
)
gradient_col <- ggplot2::scale_fill_gradient2(
  low = "blue", high = "red", 
  midpoint = 0.5, limits = c(0, 1)
)
p7<-heatmaply(
  df,
  scale_fill_gradient_fun = gradient_col
)

library(dendextend)
# Create dendrogram for rows
mycols <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")
row_dend <-  df %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  set("branches_lwd", 1) %>% 
  set("branches_k_color", mycols[1:3], k = 3) 

# Create dendrogram for columns
col_dend <-  df %>%
  t() %>%
  dist() %>%
  hclust() %>%
  as.dendrogram() %>%
  set("branches_lwd", 1) %>% 
  set("branches_k_color", mycols[1:2], k = 2)

# Visualize the heatmap
heatmaply(
  df,
  Rowv = row_dend,
  Colv = col_dend
)

heatmaply(
  df[, -c(8, 9)],
  col_side_colors = c(rep(0, 5), rep(1, 4)),
  row_side_colors = df[, 8:9]
)

heatmaply(df, cellnote = mtcars)

mat <- df
mat[] <- paste("This cell is", rownames(mat))
mat[] <- lapply(colnames(mat), function(colname) {
  paste0(mat[, colname], ", ", colname)
})

heatmaply(
  df,
  custom_hovertext = mat
)
dir.create("folder")
heatmaply(mtcars, file = "folder/heatmaply_plot.html")
browseURL("folder/heatmaply_plot.html")
dir.create("folder")
heatmaply(mtcars, file = "folder/heatmaply_plot.png")
browseURL("folder/heatmaply_plot.png")
tmp <- heatmaply(mtcars, file = "folder/heatmaply_plot.png")
rm(tmp)

