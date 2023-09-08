# bladder cancer dataset- chen-naturecomm2020

load("all_cell_filtered_CellView.Rds")
tsne.data<-readRDS("tsne.data.rds")

library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(Seurat)

tsne_plot <- data.frame(x = tsne.data[,1], y = tsne.data[,2], col = tsne.data[,4])

tsne.df<-rownames_to_column(tsne.data, var = "rowname")
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 1, "Myeloid"))
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 2, "Epithelial"))
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 3, "mCAF"))
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 4, "iCAF"))
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 5, "Endothelial"))
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 6, "T cell"))
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 7, "B cell"))
tsne.df<- tsne.df %>%mutate(dbCluster = replace(dbCluster, dbCluster == 8, "Mast cell"))

ggplot(tsne.df) + geom_point(aes(x=tsne.df[,2], y=tsne.df[,3], color=tsne.df[,5]), size=0.1)
cell.type<-tsne.df %>% dplyr::select(rowname, dbCluster)

# take all cells and make seurat object
allcells<- tsne.df
allcells.ptid<- allcells%>% separate(rowname, c("id", "second") , remove = FALSE)
allcells.ptid<-dplyr::select(allcells.ptid, c(1, 2, 4, 5, 7))

log2cpm.hgnc<-readRDS("log2cpm.hgnc.rds")

# creating metadata file
sample.info <- data.frame(cell.names = colnames(log2cpm.hgnc)) 
sample.info<-sample.info %>% separate(cell.names, c("sampleID", "second") , remove = FALSE)
sample.info<-dplyr::select(sample.info, c(1, 2))

allcells.ptid.meta<-allcells.ptid[, c("rowname", "dbCluster")]
colnames(allcells.ptid.meta)<- c("cell.names", "cell.type")
sample.meta<-full_join(sample.info, allcells.ptid.meta)
row.names(sample.meta) <- sample.meta$cell.names 
sample.meta$cell.names <- as.character(sample.meta$cell.names) 

# create and save seurat object
scrna <- CreateSeuratObject(counts = log2cpm, meta.data = sample.meta)

# standard seurat processing follows