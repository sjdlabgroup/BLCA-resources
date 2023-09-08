library(infercnv)

# read seurat file
library(Seurat)
scrna #seurat object

#make infercnv object
setwd("infercnv")
counts_matrix = GetAssayData(scrna, slot="counts")

# creating annotation file
anno <- data.frame(cell.names = colnames(counts_matrix)) 

# assigning normal patients as normal
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix, 
                                    annotations_file="anno_ptid.tsv",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names=c("normal"))
# run
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir="output_dir",  
                             cluster_by_groups=T, 
                             analysis_mode = "subclusters",
                             denoise=T,
                             HMM=T)