library(MASS)
set.seed(1)
# cell profiler output and phenotype assignment
x<-read.csv("MyExpt_IdentifyPrimaryObjects.csv", header=T)
cell<-read.csv("cellpro.assign.csv", sep = ",", header=T)

x1=x[,c("ImageNumber","AreaShape_Center_X","AreaShape_Center_Y")]
colnames(x1)=c("image","x","y")
state=round(runif(nrow(x1),0,1))
state=cell$class

# dummy state assignment. In reality, it should represent the mesenchymal/epithelial phenotypic state classification
input=cbind(x1,state)
#for mixed population c1
input=input[which(input$image > 15 & input$image < 50),]
#for pure population c2, c3 - 50 or 70
input=input[which(input$image > 15 & input$image < 70),]

dir<-"transition/"
cell.id="C"
plate="24"
for(i in 2:6) {
  ## assigning cells to coarse grids
  #2:6- when 2-large grid size, less no. of grids; 6-smaller grid size, more no. of grids
  gridz= i   ## represent n x n grids
  binz=gridz-1 ## since binz counts up 0 - n-1
  x_adj=round((input$x-min(input$x))*binz/(max(input$x)-min(input$x)))
  y_adj=round((input$y-min(input$y))*binz/(max(input$y)-min(input$y)))
  # normalizing data
  xy_adj=paste0(x_adj,"_",y_adj,sep="")
  mat=cbind(input,xy_adj)
  mat$state=as.factor(mat$state)
  mat$image=as.factor(mat$image)
  mat.list=split(mat, mat$xy_adj)
  #mat.list
  # splitting the data by xy_adj
  stat.list=lapply(mat.list,function(x) { table(x$state,x$image) })
  stat.list
  names(stat.list)
  # applying function to mat.list, finding how many of each class in each image,
  # by the xy coords
  
  # making into a data frame, columns are values of table, rows are xy_coord grops
  stat.df <- data.frame(matrix(unlist(stat.list), nrow=length(stat.list), byrow=T))
  head(stat.df)
  rownames(stat.df)=names(mat.list)
  colnames(stat.df)[seq(2,ncol(stat.df),2)] <- (rep("_C2", length(seq(2,ncol(stat.df),2)))) # rename even columns-Mero
  colnames(stat.df)[seq(1,ncol(stat.df),2)] <- (rep("_C1", length(seq(1,ncol(stat.df),2)))) # rename odd columns-Holo
  #where C1-mesenchymal, C2-epithelial, .1-source, .2-destination
  #mesenchymal-source, epithelial-source, mesenchymal-destination, epithelial-destination
  
  #stat.df.filtered
  ## filters to exclude over-crowded or sparse grids - or other filters as applicable
  stat.df.filtered=stat.df[which((apply(stat.df, 1, FUN=sum) >10) & (apply(stat.df, 1, FUN=sum) <5000)),]
  apply(stat.df.filtered, 1, FUN=sum)
  ## solves the redundant system of equations to determine the cell state transition rates
  for (j in (seq(1,(ncol(stat.df)-3),2))){
    t1<-j
    t2<-j+1
    t3<-j+2
    t4<-j+3
    A=as.matrix(stat.df.filtered[,c(t1,t2)]) #source
    b=as.matrix(stat.df.filtered[,c(t3,t4)]) #destination
    
    r.Ct=ginv(A) %*% b
    rownames(r.Ct)=c("C1_","C2_") #source
    colnames(r.Ct)=c("_C1","_C2") #destination
    sum <- apply(stat.df.filtered, 1, FUN=sum)
    write.table(as.data.frame(round(r.Ct, 3)),file = paste0(dir,"r.Ct.",cell.id,".",plate,".",i,".csv"), append=TRUE,sep = ",",col.names = F, row.names = T)
    }
  library(reshape2)
  imagenumber<-colnames(dplyr::select(reshape2::dcast(reshape2::melt(stat.list, 
                                                                     value.name = "classsum", 
                                                                     varnames = c("class", "image")),
                                                      L1~image), -1))
  imagenumber<-imagenumber[1:(((ncol(stat.df))/2)-1)] #image number
  gridsize<-rep(i, length(imagenumber))
  mincells = rep(min(sum),length(imagenumber)) 
  maxcells = rep(max(sum),length(imagenumber)) 
  grid.cells<-data.frame(gridsize, mincells, maxcells)
  write.table(as.data.frame(grid.cells),file = paste0(dir,"grid.cells.",cell.id,".",plate,".",i,".csv"), sep = ",",col.names = T, row.names = F)
  write(imagenumber, file = paste0(dir,"imagenumber.",cell.id,".",plate,".",i,".txt"))
  }

#read back 
library(tidyverse) # for group_by, slice
library(dplyr) #for bind_cols and summarise
library(stringr) #for replacing string

#read each grid
r.Ct.well.plate.grid # read transition rate file for each grid
grid.cells # read grid cells
imagenumber # read image number file

colnames(r.Ct.well.plate.grid)<-c("class", "C1","C2")
r.Ct.well.plate.grid<-r.Ct.well.plate.grid %>%
  dplyr::mutate(C2C1=C1 , C2C2=C2) %>%
  dplyr::mutate_at(c("C2C1"), list(lead), n = 1 ) %>%
  dplyr::mutate_at(c("C2C2"), list(lead), n = 1 ) %>%
  dplyr::filter(row_number() %% 2 == 1) %>%
  rename(C1C1=C1,
         C1C2=C2)

dim(r.Ct.well.plate.grid) ##check if nrow is same
length(imagenumber) #check if nrow is same
r.Ct.well.plate.grid$imagenumber<-imagenumber
r.Ct.well.plate.grid<-r.Ct.well.plate.grid[,-1]
r.Ct.well.plate.grid$phenotype<-rep("M or E", nrow(r.Ct.well.plate.grid)) # assign phenotype
r.Ct.well.plate.grid<-bind_cols(r.Ct.well.plate.grid,grid.cells)

#assign and combine all grids
r.Ct.well.plate.grid.2<-r.Ct.well.plate.grid
r.Ct.well.plate.grid.3<-r.Ct.well.plate.grid
r.Ct.well.plate.grid.4<-r.Ct.well.plate.grid
r.Ct.well.plate.grid.5<-r.Ct.well.plate.grid
r.Ct.well.plate.grid.6<-r.Ct.well.plate.grid

r.Ct.well.plate.grid.all<-bind_rows(list(r.Ct.well.plate.grid.2,
                                      r.Ct.well.plate.grid.3,
                                      r.Ct.well.plate.grid.4,
                                      r.Ct.well.plate.grid.5,
                                      r.Ct.well.plate.grid.6))

grid.phenotype.meansd<-r.Ct.well.plate.grid.all %>% group_by(gridsize, phenotype) %>%
  summarise(across(C1C1:C2C2, list(mean=mean,sd=sd)))
grid.phenotype.meansd[,-c(1:2)] <-round(grid.phenotype.meansd[,-c(1:2)],3)

#for plot
grid.phenotype.mean<-grid.phenotype.meansd %>%
  pivot_longer(!c(gridsize,phenotype), names_to = "class", values_to = "mean") %>% 
  filter(row_number() %% 2 == 1) %>% 
  mutate_at("class", str_replace, "_mean", "")
grid.phenotype.sd<-grid.phenotype.meansd %>%
  pivot_longer(!c(gridsize,phenotype), names_to = "class", values_to = "sd") %>% 
  filter(row_number() %% 2 == 0) %>% 
  mutate_at("class", str_replace, "_sd", "")
grid.phenotype.meansd.long<-left_join(grid.phenotype.mean, grid.phenotype.sd)
well<-rep(2, nrow(grid.phenotype.meansd.long))
grid.phenotype.meansd.long<-bind_cols(list(grid.phenotype.meansd.long, well=well))

#plot rate transition
ggplot(grid.phenotype.meansd.long, aes(x = mean, y = class, color = class)) + 
  geom_point(aes(size = gridsize))+
  scale_size(trans = 'reverse')+
  geom_line() +
  labs(x = "mean", y = "class") + 
  theme_bw()+
  ggtitle(well)
