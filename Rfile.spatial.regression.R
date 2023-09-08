set.seed(1)

library(Matrix)
library(spdep)
library(raster)
library(eSDM)
library(spatialreg)
library(sp)
library(ade4)
library(adegenet)
library(adegraphics)
library(adespatial)
library(raster)
library(spatialreg)
library(dplyr)
library(Seurat)

library(Seurat)

# run seurat for 10x spatial file

# load data
load("barcode_by_coor_matrix.RDS") # spatial coordinates file (spatial coord file)
load("barcode_by_cell_type_abundance_matrix.RDS") # cell type deconvolution file (decon file)
load("barcode_by_pathway_score_matrix.sct.RDS")  # pathway expression metadata file (pathway file)
load("barcode_by_expression_score_matrix.sct.RDS") # gene expression metadata file (expresssion file)

#spatial_coor
pathway_score_combined <-  pathway_score_combined %>% 
  tibble::rownames_to_column(var = "barcode") %>%
  tidyr::separate(col= barcode, 
                  into=c("barcodeID", "sampleID"), sep="_", remove =F)

decon_raw_combined <-  decon_raw_combined %>% 
  tibble::rownames_to_column(var = "barcode") %>%
  tidyr::separate(col= barcode, 
                  into=c("barcodeID", "sampleID"), sep="_", remove =F)

exp_score_combined <-  exp_score_combined %>% 
  tibble::rownames_to_column(var = "barcode") %>%
  tidyr::separate(col= barcode, 
                  into=c("barcodeID", "sampleID"), sep="_", remove =F)

exp_score_combined<- exp_score_combined %>%
  subset(select= c(barcode, barcodeID, sampleID, 
                   SMAD3, 
                   KLF4, 
                   PPARG)) # include relevant genes

#add barcodes
spatial_coor<-spatial_coor %>% tibble::rownames_to_column(var = "barcode")

#merge  scores and barcodes
pathway_score_combined<-dplyr::full_join(pathway_score_combined, spatial_coor, by = "barcode")
colnames(pathway_score_combined) 
decon_raw_combined<-dplyr::full_join(decon_raw_combined, spatial_coor, by = "barcode")
colnames(decon_raw_combined)
exp_score_combined<-dplyr::full_join(exp_score_combined, spatial_coor, by = "barcode")
colnames(exp_score_combined)

#list
coords.list=list()
coords.list[[1]]<-dplyr::filter(pathway_score_combined, sampleID == 1) %>%
  subset(select = c(row, col)) %>% dplyr::rename(X1=row) %>% dplyr::rename(X2=col)
coords.list[[2]]<-dplyr::filter(pathway_score_combined, sampleID == 2) %>%
  subset(select = c(row, col)) %>% dplyr::rename(X1=row) %>% dplyr::rename(X2=col)
coords.list[[3]]<-dplyr::filter(pathway_score_combined, sampleID == 3) %>%
  subset(select = c(row, col)) %>% dplyr::rename(X1=row) %>% dplyr::rename(X2=col)
coords.list[[4]]<-dplyr::filter(pathway_score_combined, sampleID == 4) %>%
  subset(select = c(row, col)) %>% dplyr::rename(X1=row) %>% dplyr::rename(X2=col)

pathway.list=list()
pathway.list[[1]]<-dplyr::filter(pathway_score_combined, sampleID == 1) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
pathway.list[[2]]<-dplyr::filter(pathway_score_combined, sampleID == 2) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
pathway.list[[3]]<-dplyr::filter(pathway_score_combined, sampleID == 3) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
pathway.list[[4]]<-dplyr::filter(pathway_score_combined, sampleID == 4) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))

decon.list=list()
decon.list[[1]]<-dplyr::filter(decon_raw_combined, sampleID == 1) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
decon.list[[2]]<-dplyr::filter(decon_raw_combined, sampleID == 2) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
decon.list[[3]]<-dplyr::filter(decon_raw_combined, sampleID == 3) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
decon.list[[4]]<-dplyr::filter(decon_raw_combined, sampleID == 4) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))

exp.list=list()
exp.list[[1]]<-dplyr::filter(exp_score_combined, sampleID == 1) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
exp.list[[2]]<-dplyr::filter(exp_score_combined, sampleID == 2) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
exp.list[[3]]<-dplyr::filter(exp_score_combined, sampleID == 3) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))
exp.list[[4]]<-dplyr::filter(exp_score_combined, sampleID == 4) %>%
  subset(select = -c(barcode,barcodeID,sampleID,row,col,imagerow,imagecol))

for(i in c(1:4)) {
  sp.tme=coords.list[[i]]
  path.tme=pathway.list[[i]]
  cell.tme=decon.list[[i]]
  
  coordinates(sp.tme) <- ~ X1 + X2
  tme.poly=pts2poly_centroids(data.frame(coordinates(sp.tme)), 0.5)
  neib=poly2nb(tme.poly)
  neib.listw <- nb2listw(neib, zero.policy = TRUE)
  
  moran.lm.stat<-list()
  df=cbind(path.tme$pathway,cell.tme)
  df=log10(df)
  moran.lm.stat[[i]]=lagsarlm(path.tme$pathway ~ ., data =df, neib.listw, zero.policy = TRUE)

  summary.p.value<-data.frame(pathway=summary(moran.lm.stat[[i]])$Coef[,4])
  
  summary.p.value<--log10(as.matrix(summary.p.value))
  summary.p.value <- summary.p.value[-c(1,2),]
  summary.p.value<- as.data.frame(summary.p.value)
  summary.p.value$tumorname<-as.factor(c(rep(i, nrow(summary.p.value))))
  assign(paste0("summary.p.value.", i), summary.p.value)
}

summary.p.value.list<-list(summary.p.value.1, summary.p.value.2, summary.p.value.3, summary.p.value.4)

# all samples combine
summary.p.value<-rbind(rownames_to_column(summary.p.value.1, var = "Cell.type"),
                       rownames_to_column(summary.p.value.2, var = "Cell.type"),
                       rownames_to_column(summary.p.value.3, var = "Cell.type"),
                       rownames_to_column(summary.p.value.4, var = "Cell.type"))
# plots
ggplot(summary.p.value, aes(x = tumorname, y = Cell.type, fill = pathway)) +
  scale_fill_gradient(low = "light yellow", high = "dark red") +
  geom_tile()

for (i in 2:(length(summary.p.value)-1)){
  ggplot(summary.p.value, aes(x = tumorname, y = Cell.type, fill = names(summary.p.value)[i])) +
  scale_fill_gradient(low = "light yellow", high = "dark red") +
  geom_tile()
}

for (i in names(summary.p.value)[2:(length(summary.p.value)-1)]){
  print(ggplot(summary.p.value, aes(x = tumorname, y = Cell.type, fill = .data[[i]])) +
    scale_fill_gradient(low = "light yellow", high = "dark red") +
    geom_tile())
}





