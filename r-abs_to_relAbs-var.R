#!/usr/bin/env Rscript

library(dendextend)
library(ape)
library(dplyr)
library(phylogram)


groups<-c("agalactiae", "salivarius_vestibularis", "mitis_infantis_oralis_pseudopneumoniae_pneumoniae",
             "gallolyticus_lutetiensis_infantarius_equinus", "intermedius_constellatus_anginosus_gordonii",
             "rubneri_australis_ilei", "sinensis_cristatus_gordonii_sanguinis", "sobrinus_criceti", "urinalis",
             "equi", "canis", "pyogenes", "dysgalactiae", "iniae_uberis_porcinus_halichoeri", "mutans", "parasanguinis",
             "massiliensis")


for (group in groups) {
  
  # Step 1: load trees
  PANtree<-ape::read.tree(paste0("/workspace/jmarkwelchlab/SPECIES_LEVEL_PANGENOMES/Streptococcus_species_pan/16_LAYERS/",group,"_GC_freq_tree"))
  MLtree<-ape::read.tree(paste0("/workspace/jmarkwelchlab/SPECIES_LEVEL_PANGENOMES/Streptococcus_species_pan/09_PHYLOGENOMICS/",group,"_Bac71_fasta.clean.fa.contree"))
  
  # Step 2: convert to dendrogram objects
  dend_MLtree <- ape::chronos(MLtree)
  dend_MLtree_2 <- phylogram::as.dendrogram(dend_MLtree)
  dend_PANtree <- ape::chronos(PANtree)
  dend_PANtree_2 <- phylogram::as.dendrogram(dend_PANtree)
  
  # Step 3:  make dendrogram list
  SCG_PAN <- dendlist(dend_MLtree_2, dend_PANtree_2)
  
  # Step 4: Call the pdf command to start the plot
  rows <- as.data.frame(MLtree$tip.label)
  height=nrow(rows) *0.0625
  pdf(file = paste0("/workspace/jmarkwelchlab/SPECIES_LEVEL_PANGENOMES/Streptococcus_species_pan/17_PHY_vs_PAN/", group,"_tanglegram.pdf"),
      width = 10,
      height = height) 
  
  # Step 5: plot
  SCG_PAN %>% dendextend::untangle(method="step2side") %>%
    dendextend::tanglegram(common_subtrees_color_lines=FALSE,
                           highlight_distinct_edges=FALSE,
                           common_subtrees_color_branches=FALSE,
                           highlight_branches_lwd=FALSE,
                           lwd=1,color_lines="black", 
                           lab.cex = 0.5, edge.lwd = 1, 
                           margin_inner = 20, 
                           columns_width = c(1, 0.5, 1), 
                           axes=FALSE)
  # Step 6: Run dev.off() to create the file!
  dev.off()
}