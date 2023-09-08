library(ggfortify)
library(dplyr)
library(tidyverse)
library(ggplot2)
#######
set.seed(1)
filename<-xx
cellpro.res<-read.csv("MyExpt_IdentifyPrimaryObjects.csv", header=T)
cellpro.dat$ImageNumber<-as.factor(cellpro.dat$ImageNumber)
cellpro.dat.des<-cellpro.dat %>%arrange(desc(ImageNumber)) 
col.man<-c("black","darkgreen","orange", "red", "pink") #when only 5 images

# plotting x, y coordinates
ggplot(dplyr::filter(cellpro.dat, ImageNumber %in% c(1, 2, 3, 4, 5)),
       aes(x = AreaShape_Center_X, y = AreaShape_Center_Y, col = ImageNumber)) + 
  geom_point(shape = ".") +
  facet_grid( . ~ ImageNumber)+
  scale_y_continuous(trans = "reverse")+
  scale_colour_manual(values = col.man)

cellpro.dat.des.xy$axis.x<-as.factor (cellpro.dat.des.xy$axis.x)
cellpro.dat.des.xy$axis.y<-as.factor (cellpro.dat.des.xy$axis.y)
cellpro.dat.des.xy$axis.xy<-as.factor (cellpro.dat.des.xy$axis.xy)

#for clustering
cellpro.df <- dplyr::select(cellpro.dat.des.xy, -c("ImageNumber","ObjectNumber","AreaShape_Center_X","AreaShape_Center_Y", "axis.x", "axis.y", "axis.xy"))
cellpro.meta <- dplyr::select(cellpro.dat.des.xy, c("ImageNumber","ObjectNumber","AreaShape_Center_X","AreaShape_Center_Y", "axis.x", "axis.y", "axis.xy"))

# pca
cellpro.pca<-scale(cellpro.df)
cellpro.pca <- prcomp(cellpro.pca, center = TRUE, scale = T)
dims.pca <- cellpro.pca$x[, 1:2]
colnames(dims.pca) <- c("dimension.x", "dimension.y")
# extract imagenumber info
image.num <- cellpro.meta$ImageNumber
# extract axis info
xaxis <- cellpro.meta$axis.x
yaxis <- cellpro.meta$axis.y
xyaxis <- cellpro.meta$axis.xy
dims.pca <- cbind(as.data.frame(dims.pca), image.num, xaxis, yaxis, xyaxis)
# generate plot using ggplot
ggplot(dims.pca, aes(x = dimension.x, y = dimension.y, color = image.num)) + 
  geom_point(size = 0.2) + 
  labs(x = "dimension x", y = "dimension y") + 
  theme_bw()+
  ggtitle("PCA visualization")+
  theme(plot.title = element_text(hjust = 0.5))
# generate plot using autoplot
p0<-autoplot(cellpro.pca, data = cellpro.meta, 
                           size = 0.2,
                           colour = "ImageNumber",
                           loadings = TRUE, loadings.colour = 'blue',
                           loadings.label = TRUE, loadings.label.size = 5, loadings.label.vjust = 1.2)
p0
# top loadings
library(dplyr)
library(tibble)
rot.top20<- cellpro.pca$rotation %>% as.data.frame %>% rownames_to_column%>% dplyr::select(rowname, PC1, PC2) %>% arrange(desc(PC1^2+PC2^2)) %>% head(20)
# umap
library(umap)
cellpro.umap<-umap::umap(scale(cellpro.df)[,c(rot.top20$rowname)])
dims.umap<-cellpro.umap$layout
colnames(dims.umap) <- c("dimension.x", "dimension.y")
dims.umap <- cbind(as.data.frame(dims.umap), image.num, xaxis, yaxis, xyaxis)
# make umap and color with image number
ggplot(dims.umap, aes(x = dimension.x, y = dimension.y, color = image.num)) + 
  geom_point(size = 0.2) + 
  labs(x = "dimension x", y = "dimension y") + 
  theme_bw()+
  ggtitle("UMAP visualization")+
  theme(plot.title = element_text(hjust = 0.5))

#### DESCRIPTION OF CLASSIFIERS
# solidity - minimum shape-1-solid object no indentations, <1-object with irregular boundary or holes.
# FormFactor-Equals 1 for a perfectly circular object, more spindle-like objects closer to 0. 
# eccentricity- between 0 and 1. 0-a circle, 1-a line segment
# compactness-filled circle will have a compactness of 1, with irregular objects or objects with holes-greater than 1.
# zernike-characterize the distribution of intensity across the object,higher order carry less information


# assigning Holoclone (mesenchymal-M) and Meroclone (epithelial-E)  phenotypic states
cellpro.assign <- cellpro.df %>% 
  mutate(radius = ifelse(AreaShape_MeanRadius >= mean(AreaShape_MeanRadius), "E" , "M"))%>%
  mutate(area = ifelse(AreaShape_Area >= mean(AreaShape_Area), "E" , "M"))%>%
  mutate(perimeter = ifelse(AreaShape_Perimeter >= mean(AreaShape_Perimeter), "E" , "M"))%>%
  mutate(solidity = ifelse(AreaShape_Solidity >= mean(AreaShape_Solidity), "M" , "E"))%>%
  mutate(formfactor = ifelse(AreaShape_FormFactor >= mean(AreaShape_FormFactor), "M" , "E"))%>%
  mutate(eccentricity = ifelse(AreaShape_Eccentricity >= mean(AreaShape_Eccentricity), "E" , "M"))%>%
  mutate(compactness = ifelse(AreaShape_Compactness >= mean(AreaShape_Compactness), "E" , "M"))

cellpro.assign$radius<-as.factor (cellpro.assign$radius)
cellpro.assign$area<-as.factor (cellpro.assign$area)
cellpro.assign$perimeter<-as.factor (cellpro.assign$perimeter)
cellpro.assign$solidity<-as.factor (cellpro.assign$solidity)
cellpro.assign$formfactor<-as.factor (cellpro.assign$formfactor)
cellpro.assign$eccentricity<-as.factor (cellpro.assign$eccentricity)
cellpro.assign$compactness<-as.factor (cellpro.assign$compactness)
head(cellpro.assign)

# plot colored by classifier- eg-solidity
ggplot(dims.umap, aes(x = dimension.x, y = dimension.y, color = cellpro.assign$solidity)) + 
  geom_point(size = 0.2) + 
  labs(x = "dimension x", y = "dimension y") + 
  theme_bw()+
  ggtitle("UMAP visualization")+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(paste0(filename,".hm.solidity.png"),plot = last_plot(),)


########final assignment
cellpro.assign <- cellpro.assign %>% 
  mutate(class = ifelse((perimeter == "M" & solidity == "M" & formfactor == "M" & compactness == "M"), "M" , "E"))
cellpro.assign <- cbind(cellpro.assign, cellpro.meta)

#or
cellpro.assign <- cellpro.assign %>% 
  mutate(class = ifelse((perimeter == "E" & solidity == "E" & formfactor == "E" & compactness == "E"), "E" , "M"))
cellpro.assign <- cbind(cellpro.assign, cellpro.meta)

#for scratch

#table
cellpro.assign.scratch.cells <- cellpro.assign %>% 
  group_by (ImageNumber, axis.xy, class) %>% 
  summarise(n = n()) %>% 
  mutate(phenotypepercent = (n / sum(n))*100)

#count
cellpro.assign.scratch<-cellpro.assign %>% count(axis.xy, class)

#count by class
cellpro.assign.classpercent <- cellpro.assign %>% 
  group_by (ImageNumber, class) %>% 
  summarise(n = n()) %>% 
  mutate(classpercent = (n / sum(n))*100)

#count by classaxis.xy
cellpro.assign.classpercentEaxis <- cellpro.assign %>% 
  filter (axis.xy == "E") %>% 
  group_by (ImageNumber, class) %>% 
  summarise(n = n()) %>% 
  mutate(classpercent = (n / sum(n))*100)

cellpro.assign.classpercentMaxis <- cellpro.assign %>% 
  filter (axis.xy == "M") %>% 
  group_by (ImageNumber, class) %>% 
  summarise(n = n()) %>% 
  mutate(classpercent = (n / sum(n))*100)

#scratch assay area plots
cellpro.area<-read.csv("MyExpt_Image.csv", header=T)
cellpro.area<-dplyr::select(cellpro.area, c(ImageNumber, 
                                     AreaOccupied_AreaOccupied_IdentifyPrimaryObjects,
                                     AreaOccupied_TotalArea_IdentifyPrimaryObjects,
                                     Count_IdentifyPrimaryObjects))

cellpro.area <- cellpro.area %>% 
  mutate(Time.Hours= 2*(row_number())) %>%
  mutate(max.scratch.area=max(AreaOccupied_AreaOccupied_IdentifyPrimaryObjects)-min(AreaOccupied_AreaOccupied_IdentifyPrimaryObjects)) %>% 
  mutate(Percent.area.covered=((AreaOccupied_AreaOccupied_IdentifyPrimaryObjects/max(AreaOccupied_AreaOccupied_IdentifyPrimaryObjects))*100)) %>% 
  mutate(Percent.scratch.area.covered=(((AreaOccupied_AreaOccupied_IdentifyPrimaryObjects-min(AreaOccupied_AreaOccupied_IdentifyPrimaryObjects))/max.scratch.area)*100)) %>% 
  mutate(Percent.scratch.area=100-Percent.scratch.area.covered)

#area.image number
ggplot(data=cellpro.area, aes(x=ImageNumber, y=AreaOccupied_AreaOccupied_IdentifyPrimaryObjects)) +
  geom_line()+
  geom_point()+
  coord_cartesian(expand = TRUE)

#Percent.area.covered
ggplot(data=cellpro.area, aes(x=Time.Hours, y=Percent.area.covered)) +
  geom_line()+
  geom_point()+
  coord_cartesian(ylim = c(NA, 100), expand = TRUE)

#Percent.scratch.area
ggplot(data=cellpro.area, aes(x=Time.Hours, y=Percent.scratch.area)) +
  geom_line()+
  geom_point()+
  coord_cartesian(ylim = c(NA, 100), expand = TRUE)

#Percent.scratch.area.covered
ggplot(data=cellpro.area, aes(x=Time.Hours, y=Percent.scratch.area.covered)) +
  geom_line()+
  geom_point()+
  coord_cartesian(ylim = c(0, 100), expand = TRUE)

# mean of three readings 
mean<-mean %>% group_by(Time.Hours) %>% 
  summarise(mean.Percent.area.covered = mean(Percent.area.covered),
            sd.Percent.area.covered = sd(Percent.area.covered),
            sem.Percent.area.covered = (sd(Percent.area.covered))/sqrt(length(Percent.area.covered)),
            mean.Percent.scratch.area.covered = mean(Percent.scratch.area.covered),
            sd.Percent.scratch.area.covered = sd(Percent.scratch.area.covered),
            sem.Percent.scratch.area.covered = (sd(Percent.scratch.area.covered))/sqrt(length(Percent.scratch.area.covered)),
            mean.Percent.scratch.area = mean(Percent.scratch.area),
            sd.Percent.scratch.area = sd(Percent.scratch.area),
            sem.Percent.scratch.area = (sd(Percent.scratch.area))/sqrt(length(Percent.scratch.area)))

scratch.mean<-rbind(mesenchymal.mean, epithelial.mean)
scratch.mean <- scratch.mean %>% 
  mutate(Time.Hours=Time.Hours-2)


library(ggpubr)
pval<-compare_means(c(mean.Percent.area.covered, mean.Percent.scratch.area, mean.Percent.scratch.area.covered) ~ phenotype, scratch.mean)
pval

#Percent.area.covered - barplot
ggplot(data=scratch.mean, aes(x=Time.Hours, y=mean.Percent.area.covered, fill =phenotype)) +
  geom_bar(stat="identity", position = position_dodge())+
  scale_color_manual(values=c('Red','Blue'))+
  geom_errorbar(aes(ymin=mean.Percent.area.covered-sem.Percent.area.covered, 
                    ymax=mean.Percent.area.covered+sem.Percent.area.covered,
                    color=phenotype), position=position_dodge(),
                size=0.5,
                width=1,
                linetype=4,
                alpha=0.9)+
  coord_cartesian(ylim = c(NA, 100), expand = TRUE)+
  geom_text(aes(label = paste("p=", round(pval$p[1],2)), x = 40, y = 10),col = "black", position=position_dodge())

#Percent.area.covered
ggplot(data=scratch.mean, aes(x=Time.Hours, y=mean.Percent.area.covered)) +
  geom_line(aes(color=phenotype))+
  geom_point()+
  scale_color_manual(values=c('Red','Blue'))+
  geom_errorbar(aes(ymin=mean.Percent.area.covered-sem.Percent.area.covered, 
                    ymax=mean.Percent.area.covered+sem.Percent.area.covered,
                    color=phenotype),
                size=0.5,
                width=1,
                linetype=4,
                alpha=0.9)+
  coord_cartesian(ylim = c(NA, 100), expand = TRUE)+
  scale_x_discrete(limits = c(0, 12, 24, 36, 48))+
  geom_text(aes(label = paste("p=", round(pval$p[1],2)), x = 40, y = 10),col = "black")


#Percent.scratch.area
ggplot(data=scratch.mean, aes(x=Time.Hours, y=mean.Percent.scratch.area)) +
  geom_line(aes(color=phenotype))+
  geom_point()+
  scale_color_manual(values=c('Red','Blue'))+
  geom_errorbar(aes(ymin=mean.Percent.scratch.area-sem.Percent.scratch.area, 
                    ymax=mean.Percent.scratch.area+sem.Percent.scratch.area,
                    color=phenotype),
                size=0.5,
                width=1,
                linetype=4,
                alpha=0.7)+
  coord_cartesian(ylim = c(NA, 100), expand = TRUE)+
  geom_text(aes(label = paste("p=", round(pval$p[2],2)), x = 40, y = 100),col = "black")

#Percent.scratch.area.covered
ggplot(data=scratch.mean, aes(x=Time.Hours, y=mean.Percent.scratch.area.covered)) +
  geom_line(aes(color=phenotype))+
  geom_point()+
  scale_color_manual(values=c('Red','Blue'))+
  geom_errorbar(aes(ymin=mean.Percent.scratch.area.covered-sem.Percent.scratch.area.covered, 
                    ymax=mean.Percent.scratch.area.covered+sem.Percent.scratch.area.covered,
                    color=phenotype),
                size=0.5,
                width=1,
                linetype=4,
                alpha=0.9)+
  coord_cartesian(ylim = c(NA, 100), expand = TRUE)+
  geom_text(aes(label = paste("p=", round(pval$p[3],2)), x = 40, y = 10),col = "black")
