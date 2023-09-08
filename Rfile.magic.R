library(Rmagic)
library(phateR)
library(readr)
library(ggplot2)
library(viridis)
library(Seurat)

scrna<-readRDS("Seuratobject.rds") 
magic.dat <- magic(scrna)
magic.dat@active.assay <- 'MAGIC_RNA'