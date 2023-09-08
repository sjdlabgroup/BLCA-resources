# get TCGA data
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("GSVA")
library("ComplexHeatmap")

TCGAbiolinks:::getProjectSummary("TCGA-BLCA")

#gene names
query.exp.blca <- GDCquery(
  project = "TCGA-BLCA",
  legacy = T,
  data.category = "Gene expression", 
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq",
  experimental.strategy = "RNA-Seq",
  file.type = "results", #imp for nodups,
  sample.type = "Primary Tumor",
  access = "open")

#downlaod files
GDCdownload(query.exp.blca)
exp.blca <- GDCprepare(query = query.exp.blca, summarizedExperiment = T)

# data preprocessing
# deseq dataobject
dge = DGEList(counts=assay(exp.blca),
              samples=colData(exp.blca),
              genes=as.data.frame(rowData(exp.blca)))
  
group<-factor(make.names(paste(exp.blca@colData$vital_status,
                               exp.blca@colData$tissue_or_organ_of_origin,
                               exp.blca@colData$gender,
                               exp.blca@colData$primary_site,
                               sep="."), unique=F))
keep = filterByExpr(dge,group = group)
dge = dge[keep,keep.lib.sizes=FALSE]

# Normalization (TMM followed by voom)
dge = calcNormFactors(dge)

# log-CPM values are used for exploratory plots
logcpm<-cpm(dge, log=T) #for analysis other than differential expression, voom output and logcpm are similar
logcpmscale<- t(scale(t(logcpm))) #logzscore

# differential expression analysis
# voom transformation and calculation of variance weights
# design matrix
design<-model.matrix(~0+group)
colnames(design)<-levels(group)
v = voom(dge, design, plot=TRUE)

# Fit model in limma using weights frm voom and comparison between groups
fit = lmFit(v, design)
contr <- makeContrasts(Alive.Bladder..NOS.female.Bladder - Dead.Bladder..NOS.female.Bladder, 
                       levels = colnames(design))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)

# Show top genes
topGenes = topTable(tmp, number=100, sort.by="p")
head(topGenes)

limma_res<-list(
  logcpm=logcpm,
  logcpmscale=logcpmscale,
  voomObj=v, # normalized data
  fit=tmp, # fitted model and statistics
  topGenes=topGenes)

#pca
library(PCAtools)
limma_res<-readRDS("limma_res.RDS")
mat<-(limma_res$voomObj$E)
metadata<-(limma_res$voomObj$targets)
p <- pca(limma_res$voomObj$E, metadata = limma_res$voomObj$targets, removeVar = 0.1)
pca_res<-prcomp(t(limma_res$voomObj$E))

ggplot(as.data.frame(pca_res$x), aes(x = PC1, y = PC2, col = limma_res$voomObj$targets$vital_status)) +
  geom_point() +
  labs(title = "PCA analysis of voom-normalized expression")


# gsva for gene set enrichment
# input-expression dataset
library(tidyr)
library(tidyverse)
mat<-as.data.frame(limma_res$voomObj$E) %>% rownames_to_column("rowname") %>% 
  separate (rowname, c("rowname", NA), sep = "\\|", convert = T) %>% 
  distinct(rowname, .keep_all= TRUE) %>%
  column_to_rownames(var = "rowname") %>%
  as.matrix()

library(Biobase)
eSet<-ExpressionSet(as.matrix(mat),
                    phenoData=AnnotatedDataFrame(limma_res$voomObj$targets),
                    featureData = AnnotatedDataFrame(limma_res$voomObj$genes))

# this used as dup row names
eSet<-ExpressionSet(as.matrix(mat),
                    phenoData=AnnotatedDataFrame(limma_res$voomObj$targets))

# input-geneset collection
genesetlist<-readRDS("genesetlist.rds")
head(genesetlist)
names(genesetlist)

for (i in 1: length(genesetlist)){
  genesetname<-names(genesetlist[i])
  assign(paste0("gs", i), GSEABase::GeneSet(as.character(unlist(genesetlist[[i]])), setName=genesetname))
}

geneSets <- GSEABase::GeneSetCollection(gs1, ...)

# run gsva
# default, method gsea
es <- gsva(eSet, geneSets)
esplot<-exprs(es)
cols_group = cutree(hclust(dist(t(esplot))), k=5) 

Heatmap(esplot, show_row_dend = FALSE, row_names_side = "left", 
        row_order = geneset_order, show_column_names = F, column_dend_reorder = T,
        column_names_rot = 0, column_names_side = "top",
        column_names_gp = grid::gpar(fontsize = 20),
        row_names_gp = grid::gpar(fontsize = 20))

# corr plot
library(corrplot)
round(cor(t(esplot)), digits = 2)
corrplot(cor(t(esplot)),method = "number",type = "upper")
corrplot(cor(t(esplot)),method = "color",type = "upper")
corrplot(cor(t(esplot)), method = "color",type = "upper", addCoef.col = "black", diag=FALSE )

# for umap
eSet<-readRDS("eSet.RDS")

# pca
library(ggfortify)
pca_results <- prcomp((t(exprs(eSet))))
pca_plot_df <- data.frame(pca_results$x) %>%
  tibble::rownames_to_column("sampleID") %>%
  dplyr::inner_join(meta, by = "sampleID")
ggplot(pca_plot_df,aes(x = PC1,y = PC2)) +
  geom_point() 
ggplot(pca_plot_df,aes(x = PC1,y = PC2,color = emt))+ 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))+
  geom_point() 

library(tidyverse)
pca_loadings <- pca_results$rotation
pca_loadings <- pca_loadings %>% 
  as_tibble(rownames = gene)
 
top_genes <- pca_loadings %>% 
 dplyr::select(gene, PC1, PC2) %>%
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  group_by(PC) %>% 
  arrange(desc(abs(loading))) %>%
  slice(1:10) %>%
  pull(gene) %>% 
  unique()

top_loadings <- pca_loadings %>% 
  filter(gene %in% top_genes)

ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "black") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))

# correlation plot cell type
corrplot(cor(t(esplot.cor)),method = "color",type = "upper")

# scatterpie plot
library(ggplot2)
library(scatterpie)
umap_plot_df_x<-data.matrix(umap_plot_df[c("X1", "X2")])
rownames(umap_plot_df_x)<- umap_plot_df$sampleID
umap_plot_df_y<-data.matrix(umap_plot_df[c(coumn_names)])
rownames(umap_plot_df_y)<- umap_plot_df$sampleID
umap_plot_df_pie<-as.data.frame(cbind(umap_plot_df_x, umap_plot_df_y)) 
celltype.percent<-umap_plot_df_pie %>% dplyr::select(c(barcode, type, value, value_total))
celltype.percent<-epi.percent %>% mutate(percent=(value/value_total)*100)
ggplot(umap_plot_df_pie_long_dat) + 
  geom_arc_bar(aes(x0 = X1, y0 = X2, r0 = 0, r=0.04,
                   #r = 0.1
                   start = start_angle, end = end_angle, fill = type),
               color= NA) +
  coord_fixed() + theme(aspect.ratio=1)

# survival plots
library(survival)
# relevant metadata
clin_df <- eSet@phenoData@data[, c("patient",
                                   "vital_status",
                                   "days_to_death",
                                   "days_to_last_follow_up")]


# create boolean column deceased and overall survival
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
clin_df$gene <- t(exprs(eSet))[rownames(clin_df), gene]
median_value <- median(clin_df$gene)
clin_df$med = ifelse(clin_df$gene >= median_value, "UP", "DOWN")
clin_df$med <- factor(clin_df$gene, levels = c("UP","DOWN"))

fit = survfit(Surv(overall_survival, deceased) ~ med, data=clin_df)

# Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, pval=T, ggtheme = theme(aspect.ratio = 1))

# hazard ratio
res.cox <- coxph(Surv(overall_survival, deceased) ~ gene, data =  clin_df)
ggforest(res.cox)