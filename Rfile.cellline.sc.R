library("Seurat")
library("sctransform")
library("dplyr")
library("ggplot2")

# read cell ranger output-feature-bardcode matrices
data.10x = list()
data.10x[[i]] <- Read10X(data.dir = xx)

# make seurat objects
scrna = CreateSeuratObject(counts = data.10x)

# get cell cycle genes
cell.cycle.tirosh <- read.csv("http://genomedata.org/rnaseq-tutorial/scrna/CellCycleTiroshSymbol2ID.csv", header=TRUE); # read in the list of genes
s.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G1/S")]; # create a vector of S-phase genes
g2m.genes = cell.cycle.tirosh$Gene.Symbol[which(cell.cycle.tirosh$List == "G2/M")]; # create a vector of G2/M-phase genes

#QC plots

scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)

plot1 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.1)
plot2 <- FeatureScatter(object = scrna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.1)
plot1+plot2

# calculate cell cycle score for each cell
scrna <- CellCycleScoring(scrna, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)

# filter cells
scrna <- subset(scrna, subset = nFeature_RNA > Feature01  & nCount_RNA < Count99 & percent.mt < 10)

# normalize - normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result
scrna <- NormalizeData(scrna)
scrna <- FindVariableFeatures(scrna)

# scale and center data - mean expression across cells=0, variance across cells =1
scrna <- ScaleData(scrna, features = rownames(scrna)) 

# reduce dimensionality
scrna <- RunPCA(scrna)
DimPlot(scrna, reduction = "pca")
scrna<-RunTSNE(scrna, dims = 1:10, reduction = "pca")
DimPlot(scrna, reduction = "tsne") 
scrna<-RunUMAP(scrna, dims = 1:10, reduction = "pca")
DimPlot(scrna, reduction = "umap") 
# find clusters
scrna <- FindNeighbors(object=scrna, dims = 1:10)
scrna <- FindClusters(object=scrna, resolution = 0.5)

# composition plot
ggplot(scrna@meta.data, aes(x=seurat_clusters , fill = DataSet)) + 
  geom_bar()

# compute module score and plot
library(RColorBrewer) #for brewer.pal
scrna <- AddModuleScore(object = scrna,
                        features = list(features))
FeaturePlot(scrna, features = features) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

scrna.dat<-as.data.frame(t(GetAssayData(scrna)))

# violin plots
scrna.dat %>% rownames_to_column() %>% select(rowname, orig.ident, gene1, gene2) %>% pivot_longer(cols=c(gene1, gene2), names_to = "metadata", values_to = "count") %>%
  mutate(meta = factor(metadata, levels=c(gene1, gene2))) %>%
  ggplot(aes(x = orig.ident, y = count, fill = metadata)) +
  geom_violin() +
  stat_compare_means(aes(group=!!ensym(x)), label = "p.signif", method="t.test", comparisons = combn(1:length(unique(scrna.dat$orig.ident)), 2, FUN = list)) +
  stat_summary(fun ="median", geom = "text", aes(label=round(..y..,2), color=metadata, vjust=-5), size=5)+
  theme_bw() 

#############
# Differential expression
library("enrichR")
DEenrichRPlot(object = non, ident.1 = i,ident.2 = j, enrich.database = "MSigDB_Oncogenic_Signatures", max.genes = 50)

# Differential expression to pathways
library(presto) #presto fast Wilcoxon rank sum test for gsea
library(dplyr)
library(tibble)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)

de.genes <- wilcoxauc(non, 'seurat_clusters')
cluster.genes<- de.genes %>% dplyr::filter(group == "1") %>% arrange(desc(auc)) %>% dplyr::select(feature, auc)
ranks<- deframe(cluster.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "H")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

# for fast gsea
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaRes<- fgseaRes[-22,]
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(abs(NES)))
fgseaResTidy %>% dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% arrange(padj) %>% head()

# only plot the top pathways
ggplot(fgseaResTidy %>% head(n= 10) , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = NES>0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways") + 
theme_minimal()