### Load packages --------------------------------------------------------------------------------------------------------------------------------------
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(fgsea)
library(devtools)
library(msigdbr)
library(org.Hs.eg.db)
library(biomaRt)
library(GEOquery)
library(ggrepel)
library(gprofiler2)
library(ggbreak)
library(GO.db)
library(reshape2)
library(clusterProfiler)
library(factoextra)
library(tidyr)

### Define functions --------------------------------------------------------------------------------------------------------------------------------------

#load pathways for fgsea
LoadGSEApathways <- function(gsea_db,
                             species = "Homo sapiens",
                             genenames = "ensembl"){
  
  # load correct geneset
  if(gsea_db == "H" | gsea_db == "C3"){
    msigdb <- msigdbr(species = species, category = gsea_db)
  }else{
    msigdb <- msigdbr(species = species, subcategory = gsea_db)
  }
  
  # split on name or symbol
  if(genenames == "ensembl"){
    pathways <- split(msigdb$ensembl_gene, msigdb$gs_name)
  }else{
    pathways <- split(msigdb$gene_symbol, msigdb$gs_name)
  }
  
  #return
  return(pathways)
}


# plot results from fgsea based on NES
plotGSEA <- function(gseaRes, 
                     ntop = 10,
                     showgridlines = T,
                     textsize = 15,
                     title = "",
                     col_sig = "cadetblue",
                     col_ns = "gray",
                     showlegend = T){
  
  gsea_db_tidy <- strsplit(gseaRes[[1,1]], split = "_")[[1]][1]
  
  topPathwaysUp <- gseaRes %>% 
    filter(ES > 0) %>% 
    arrange(desc(NES)) %>%
    head(ntop)
  
  topPathwaysDown <- gseaRes %>% 
    filter(ES < 0) %>% 
    arrange(NES) %>%
    head(ntop)
  
  gseaplot <- rbind(topPathwaysUp, topPathwaysDown) %>%
    mutate(pathway_tidy = gsub(paste0(gsea_db_tidy, "_"), "", pathway)) %>%
    transform(pathway_tidy = gsub("_", " ", pathway_tidy)) %>%
    ggplot(aes(reorder(pathway_tidy, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="", 
         y="Normalized enrichment score") + 
    theme_minimal(base_size = textsize) +
    scale_fill_manual(values = c(col_ns, col_sig))
  
  if(!showgridlines){
    gseaplot <- gseaplot +
      theme(panel.grid = element_blank())
  }
  
  if(!showlegend){
    gseaplot <- gseaplot +
      theme(legend.position = "")
  }
  
  if(title == ""){
    gseaplot <- gseaplot +
      ggtitle(gsea_db_tidy) +
      theme(plot.title = element_text(hjust = .5))
  }else{
    gseaplot <- gseaplot +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = .5))
  }
  
  print(gseaplot)
}




### Load & organize data --------------------------------------------------------------------------------------------------------------------------------------

##### metadata and counts #####

# metadata
SampleInfo <- read.csv("SampleInfo.csv", row.names = 1)

# count data
Cts <- read.csv("Count_table.csv", row.names = 1)

# ensure order
SampleInfo <- SampleInfo[colnames(Cts), ]

# Pre-filtering of low count genes 
Cts = Cts[rowSums(Cts) >= 10, ]




##### Optional: select protein coding genes only #####
mart <- useMart("ensembl","hsapiens_gene_ensembl") 

GeneBiotypes <- getBM(c("ensembl_gene_id",
                       "hgnc_symbol",
                       "gene_biotype"), 
                     mart = mart)

ProteinCoding <- GeneBiotypes %>%
  filter(gene_biotype == "protein_coding")

Cts <- Cts %>%
  rownames_to_column(var = "gene_id") %>%
  filter(gene_id %in% ProteinCoding[,"ensembl_gene_id"]) %>%
  column_to_rownames(var = "gene_id")




##### Create DESeq2 object #####

# create object
dds <- DESeqDataSetFromMatrix(countData = Cts,
                             colData = SampleInfo,
                             design= ~ Condition) 
dds <- DESeq(dds)

# save normalized counts in file
nCts <- counts(dds, normalized=T)
write.table(nCts, "Normalized_counts.csv", sep = ",")





### Data exploration  --------------------------------------------------------------------------------------------------------------------------------------

##### Gene dispersion plot #####
plotDispEsts(dds)
ggsave("Gene_dispersion_plot.png")


##### PCA and top loadings #####

# transform data
rld <- rlogTransformation(dds)

# select top 500 variable features and perform PCA
rv <- matrixStats::rowVars(assay(rld), useNames = TRUE)
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]

pca <- prcomp(t(assay(rld)[select, ]))

# plot PCA
autoplot(pca,
         x = 1,
         y = 2,
         data = SampleInfo,
         colour = "Condition",
         main = "",
         legend = T,
         size = 2) +
  theme_classic(base_size = 20) +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_blank(),
        legend.title = element_blank()) +
  scale_colour_manual(values = c("#0072b2", "#d55e00"))
ggsave("PCA_condition.png")


# extract top laodings of PC1 --> change if your group separates over another PC!!!
pca_contrib <- pca$rotation %>% as.data.frame() %>% dplyr::arrange(PC1)
top_contrib <- c(head(rownames(pca_contrib), 10), tail(rownames(pca_contrib), 10))

# plot top loadings
pca_contrib[top_contrib, ] %>% 
  rownames_to_column(var = "Feature") %>%
  ggplot(aes(reorder(Feature, PC1), PC1)) +
    geom_col(aes(fill=PC1<=0)) +
    coord_flip() +
    theme_minimal(base_size = 15) +
    scale_fill_manual(values = c("#0072b2", "#d55e00")) +
    ggtitle("Top loadings component 1") +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          legend.title = element_blank(),
          legend)
ggsave("Top_loadings_PC1.png")  




### Differential expression  --------------------------------------------------------------------------------------------------------------------------------------

# results for padj <= .05
res <- results(dds, alpha = 0.05)
write.table(res, "DEA_results.csv", sep = ",")


# volcano plot
res %>%
  as.data.frame() %>%
  filter(!is.na(padj)) %>%
  mutate(colour = ifelse(padj > .05, "Not significant", 
                         ifelse(log2FoldChange > 0, "Upregulated",
                                "Downregulated"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = colour)) + 
    geom_point(size = 1, show.legend = T) +
    theme_classic(base_size = 20) +
    #geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
    scale_color_manual(values = c("blue", "grey", "red")) +
    xlab("log2(Fold Change)") +
    ylab("-log10(padj)") +
    theme(legend.title = element_blank())
ggsave("Volcano_plot.png")

# MA - plot
plotMA(res, ylim=c(-2,2))
ggsave("MA_plot.png")




### GSEA  --------------------------------------------------------------------------------------------------------------------------------------

##### order genelist #####
GeneList <- setNames(res$log2FoldChange, rownames(res)) 
GeneList  <-  GeneList[order(GeneList, decreasing = T)]

##### load pathways #####
databases <- c("GO:BP", 
               "GO:MF", 
               "GO:CC", 
               "H")
names(databases) <- c("GOBP", 
                      "GOMF", 
                      "GOCC", 
                      "Hallmark")
  
  
pathways <- list()

for (i in 1:length(databases)){
  pathways[[i]] <- LoadGSEApathways(databases[i])
  names(pathways)[i] <- names(databases)[i]
}



##### run #####
gseaRes <- list()

for (i in 1:length(pathways)){
  set.seed(42)
  
  gseaRes[[i]] <- fgsea(pathways[[i]],
                        GeneList, 
                        minSize = 15, 
                        maxSize = 300,
                        nPermSimple = 10000)
  names(gseaRes)[i] <- names(pathways)[i]
}


##### plot #####

for(i in 1:length(gseaRes)){
  
  if(length(gseaRes[[i]]) > 0){
    
    plotGSEA(gseaRes[[i]],
             ntop = 10,
             showgridlines = T,
             textsize = 15)
    ggsave(paste0("GSEA_", names(gseaRes)[i], ".png"))
    
  }
}




