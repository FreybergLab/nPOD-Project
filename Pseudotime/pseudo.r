
library(tidyverse)
packageVersion("tidyverse") ## ‘2.0.0’
library(monocle3)
packageVersion("monocle3") # ‘1.4.26’
library(Seurat)
packageVersion("Seurat") ## ‘5.3.1’
library(SeuratWrappers)
packageVersion("SeuratWrappers") ## ‘0.4.0’
library(igraph)
packageVersion("igraph") ## ‘2.1.4’
library(ggpubr)
source("Pseudotime/R/irlba_small_matrix_fix.R")   # lets preprocess_cds() run on the 5-feature matrix; see that file
packageVersion("ggpubr") ## ‘0.6.1’


list.files(getwd())

##### input
mtx <- read.delim("./Data/quad_data.csv", header = T, sep = ',')
dim(mtx)

mtx$rownumber <- 1:nrow(mtx)

features <- c("Count.ins", "Count.glu", "Count.soma", "Count.PP", "Count.others")
group_levels <- c("ND", "0-5\nyears\nT1D", "≥6\nyears\nT1D")
comparisons <- list(c(group_levels[1], group_levels[3]), c(group_levels[2], group_levels[3]), c(group_levels[1], group_levels[2]))


com1 <-  mtx%>%dplyr::select(rownumber,Total.Cells,Count.ins,Count.glu,Count.soma,Count.PP)%>%
  mutate("Count.others" = Total.Cells - (Count.ins+Count.glu+Count.soma+Count.PP))%>%
  mutate(cells = paste("cells", rownumber, sep = '_'))%>%
  dplyr::select(-c(Total.Cells, rownumber))%>%
  select(cells, Count.ins, Count.glu, Count.soma, Count.PP, Count.others)%>%
  column_to_rownames(var = "cells")%>%t()%>%data.frame()


com_meta <- mtx%>%dplyr::select(rownumber, Sex, AgeCat, Diagnosis, Disease.Duration, ImageID, Donor, Region)%>%
  mutate(cells = paste("cells", rownumber, sep = '_'))%>%column_to_rownames(var = "cells")

#### create a obj
obj1 <- CreateSeuratObject(com1)
obj1 <- AddMetaData(obj1, com_meta)
obj1
# An object of class Seurat 
# 5 features across 202585 samples within 1 assay 
# Active assay: RNA (5 features, 0 variable features)
# 1 layer present: counts


######## Monocle3 STAR!!!!
cds <- as.cell_data_set(obj1)
rowData(cds)$gene_short_name <- c("Count.ins", "Count.glu", "Count.soma", "Count.PP", "Count.others")

### add group 
cds$durationcat <- " "
cds$durationcat[cds$Disease.Duration < 6] <- "0-5 years"
cds$durationcat[cds$Disease.Duration >= 6] <- "≥6 years"

cds$durationcattag <- paste0('(',cds$Disease.Duration, ' yr)')
cds$durationcattag[cds$Disease.Duration < 1] <- paste0('(<1 yr)')
cds$durationcattag[cds$Diagnosis == "ND"] <- ''

cds$Donortime <- paste0(obj1$Donor, '\n', cds$Diagnosis, cds$durationcattag)

cds <- preprocess_cds(cds, num_dim = 5)
cds <- reduce_dimension(cds, umap.min_dist = 0.35) ## default: PCA

### pseudo analysis func
analysisfun <- function(cds_sub, root_clusters, min_component_size = 40) {
  ### sub cluster
  cds_sub <- cluster_cells(cds_sub, k = 15)
  
  # Learn trajectory graph
  cds_sub <- learn_graph(cds_sub, learn_graph_control = list(minimal_branch_len = 30), close_loop = FALSE)
  
  # Select root cells
  clusters <- cds_sub@clusters$UMAP$clusters
  
  root_cells <- names(clusters)[as.character(clusters) %in% as.character(root_clusters)]
  
  
  # Order cells
  cds_sub <- order_cells(cds_sub, root_cells = root_cells)
  
  # Prune small graph components
  g <- principal_graph(cds_sub)[["UMAP"]]
  comp <- igraph::components(g)
  
  small_comps <- which(comp$csize < min_component_size)
  
  if (length(small_comps) > 0) {
    
    nodes_to_remove <- igraph::V(g)[comp$membership %in% small_comps]
    g_pruned <- igraph::delete_vertices(g, nodes_to_remove)
    cds_sub@principal_graph[["UMAP"]] <- g_pruned
  }
  
  # Extract pseudotime
  pse_df <- data.frame(
    cells = names(pseudotime(cds_sub)),
    pseudotime = as.numeric(pseudotime(cds_sub)))
  
  pse_df <- com_meta %>%
    tibble::rownames_to_column("cells") %>%
    dplyr::inner_join(pse_df, by = "cells")
  
  # Disease duration groups
  pse_df$durationcat <- NA_character_
  
  pse_df$durationcat[
    pse_df$Disease.Duration < 6
  ] <- "0-5\nyears\nT1D"
  
  pse_df$durationcat[
    pse_df$Disease.Duration >= 6
  ] <- "≥6\nyears\nT1D"
  
  pse_df$durationcat[
    pse_df$Diagnosis == "ND"
  ] <- "ND"
  
  pse_df$durationcat <- factor(
    pse_df$durationcat,
    levels = c(
      "ND",
      "0-5\nyears\nT1D",
      "≥6\nyears\nT1D"
    )
  )
  
  return(list(
    cds = cds_sub,
    pse_df = pse_df
  ))
}

### pseudo analysis plot func
plotfun <- function(res) {
  
  cds_sub <- res$cds
  pse_df <- res$pse_df
  
  # Trajectory plot
  trajectory_fig <- plot_cells(
    cds = cds_sub,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    label_roots = FALSE,
    cell_size = 0.4)+
    scale_color_gradientn(colors = c("#4575B4", "#91BFDB", "#E0F3F8","#FFFFBF", "#FEE090", "#FC8D59","#D73027"))+
    scale_x_reverse()
  
  # Comparisons
  comparisons <- list(
    c("ND", "≥6\nyears\nT1D"),
    c("0-5\nyears\nT1D", "≥6\nyears\nT1D"),
    c("ND", "0-5\nyears\nT1D"))
  
  # Violin plot
  violin_fig <- ggplot(pse_df,aes(x = durationcat, y = pseudotime))+
    geom_violin(aes(fill = Diagnosis))+
    geom_boxplot(width = 0.04,fill = "white",alpha = 0.6,outlier.shape = NA)+
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test",
      label = "p.signif",
      format.args = list(digits = 2),
      tip.length = 0,
      symnum.args = list(cutpoints = c(0, 0.01, 1),symbols = c("***", "ns")))+
    scale_fill_manual(values = c("#5C94C8","#F4553D","#F4553D"))+
    theme_classic() +
    labs(x = "", y = "Pseudotime") +
    theme(axis.text = element_text(size = 13, colour = "black"),
          axis.title = element_text(size = 16, colour = "black"))
  
  return(list(
    trajectory = trajectory_fig,
    violin = violin_fig))
}


############### sub cluster
### Region
cds_age_head <- cds[, colData(cds)$Region == "Head"] ## 6 88072
cds_age_body <- cds[, colData(cds)$Region %in% c("Body", "Tail")] ## 6 116172

#### age 
cds_age_y <- cds[, colData(cds)$AgeCat == "<20_years"] ## 6 106296
cds_age_m <- cds[, colData(cds)$AgeCat == "20-40_years"] ## 6 97948



### -----------------------
### figure
### -----------------------

## fig 9A
figure9A <- plot_cells(cds, genes=c("Count.ins", "Count.glu", "Count.soma", "Count.PP"),
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_grid(~feature_label,scales = "free",
             labeller = as_labeller(c(
               "Count.ins" = "\u03B2-cells",
               "Count.glu" = "\u03B1-cells",
               "Count.soma" = "\u03B4-cells",
               "Count.PP" = "PP-cells"
             )))+
  theme(legend.title = element_text(hjust = 0.3))+
  scale_x_reverse()
print(figure9A)


## fig 9B
#### set order
cds$Region <- factor(cds$Region, levels = c("Head", "Body", "Tail"))

figure9B.ins <- plot_cells(cds, genes=c("Count.ins"),
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+Region, scales = "free", ncol = 3)+
  scale_x_reverse()
print(figure9B.ins)


## fig S21A
#### set order
cds$Region <- factor(cds$Region, levels = c("Head", "Body", "Tail"))

figureS21A.glu <- plot_cells(cds, genes=c("Count.glu"), 
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+Region, scales = "free", ncol = 3)+
  scale_x_reverse()
print(figureS21A.glu)

## fig S21C
#### set order
cds$Region <- factor(cds$Region, levels = c("Head", "Body", "Tail"))

figureS21C.soma <- plot_cells(cds, genes=c("Count.soma"), 
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+Region, scales = "free", ncol = 3)+
  scale_x_reverse()
print(figureS21C.soma)

## fig S21E
#### set order
cds$Region <- factor(cds$Region, levels = c("Head", "Body", "Tail"))

figureS21E.PP <- plot_cells(cds, genes=c("Count.PP"), 
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+Region, scales = "free", ncol = 3)+
  scale_x_reverse()
print(figureS21E.PP)


## fig 9C&D
############ head
head_res <- analysisfun(
  cds_age_head,
  root_clusters = c(30, 13, 15, 14, 24),
  min_component_size = 40)

head_fig <- plotfun(head_res)

### fig 9C&D head
print(head_fig$trajectory)
print(head_fig$violin)


############ body_tail
body_res <- analysisfun(
  cds_age_body,
  root_clusters = c(45, 46, 20, 11, 9, 30,3, 32, 37, 15, 35))

body_fig <- plotfun(body_res)

### fig9 C&D body
print(body_fig$trajectory)
print(body_fig$violin)


## fig 9E age
figure9E.ins.age <- plot_cells(cds, genes=c("Count.ins"), 
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+AgeCat, scales = "free", ncol = 2)+
  scale_x_reverse()
print(figure9E.ins.age)


## fig S21B
figureS21B.glu.age <- plot_cells(cds, genes=c("Count.glu"), 
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+AgeCat, scales = "free", ncol = 2)+
  scale_x_reverse()
print(figureS21B.glu.age)


## fig S21D
figureS21D.soma.age <- plot_cells(cds, genes=c("Count.soma"), 
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+AgeCat, scales = "free", ncol = 2)+
  scale_x_reverse()
print(figureS21D.soma.age)

## fig S21F
figureS21F.PP.age <- plot_cells(cds, genes=c("Count.PP"), 
           show_trajectory_graph = FALSE,
           label_cell_groups = F,
           label_groups_by_cluster = F,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE)+
  facet_wrap(~Diagnosis+AgeCat, scales = "free", ncol = 2)+
  scale_x_reverse()
print(figureS21F.PP.age)

## fig 9F&G
############ AGE < 20
age_y_res <- analysisfun(
  cds_age_y,
  root_clusters = c(40, 23, 17, 13, 28, 38))

age_y_fig <- plotfun(age_y_res)

## fig 9F&G < 20
print(age_y_fig$trajectory)
print(age_y_fig$violin)


############ AGE 20–40
age_m_res <- analysisfun(
  cds_age_m,
  root_clusters = c(17, 7, 24, 42, 18, 35))

age_m_fig <- plotfun(age_m_res)

## fig 9F&G 20-40
print(age_m_fig$trajectory)
print(age_m_fig$violin)


## fig S20
### set order
cds$Donortime <- factor(cds$Donortime, levels = c("6516\nND", "6476\nND", "6530\nND", "6416\nND", "6501\nND", "6439\nND", "6434\nND", "6453\nND", "6492\nND", "6525\nND", "6481\nND", "6468\nND", "6531\nND", "6487\nND", "6417\nND", "6401\nND", "6482\nND", "6470\nND", "6448\nND", "6493\nND", "6440\nND", "6535\nND", "6431\nND", "6430\nND", "6445\nND", "6415\nND", "6471\nND", "6462\nND", "6467\nND", "6486\nND", "6511\nND", "6461\nND", "6474\nND", "6465\nND", "6442\nND", "6552\nND", "6583\nND", "6548\nND", "6560\nND", "6454\nND", "6426\nND", "6446\nND", "6425\nND", "6479\nND", "6537\nND", "6420\nND", "6557\nND", "6412\nND", "6413\nND", "6460\nND", "6546\nND", "6544\nND", "6543\nND", "6547\nND", "6234\nND", "6427\nND", "6178\nND", "6179\nND", "6162\nND",
                                                  "6533\nT1D(<1 yr)", "6528\nT1D(<1 yr)", "6534\nT1D(<1 yr)", "6520\nT1D(<1 yr)", "6551\nT1D(<1 yr)", "6405\nT1D(<1 yr)", "6456\nT1D(<1 yr)", "6563\nT1D(<1 yr)", "6362\nT1D(<1 yr)", "HPAP078\nT1D(<1 yr)", "6526\nT1D(1 yr)", "6593\nT1D(1 yr)", "5000\nT1D(2 yr)", "6436\nT1D(3 yr)", "6523\nT1D(3 yr)", "6536\nT1D(4 yr)", "6441\nT1D(5 yr)", "HPAP030\nT1D(5 yr)", "6466\nT1D(6 yr)", "6473\nT1D(6 yr)", "6529\nT1D(7 yr)", "6506\nT1D(8 yr)", "6443\nT1D(9 yr)", "6519\nT1D(10 yr)", "6459\nT1D(10 yr)", "6438\nT1D(10 yr)", "HPAP048\nT1D(10 yr)", "6418\nT1D(11 yr)", "6180\nT1D(11 yr)", "HPAP041\nT1D(12 yr)", "HPAP068\nT1D(13 yr)", "6435\nT1D(14.75 yr)", "6510\nT1D(18 yr)", "6422\nT1D(21 yr)"))


figureS20.ins <- plot_cells(cds, genes=c("Count.ins"),
                            show_trajectory_graph = FALSE,
                            label_cell_groups = F,
                            label_groups_by_cluster = F,
                            label_leaves=FALSE,
                            label_branch_points=FALSE,
                            label_roots = FALSE)+
  facet_wrap(~Donortime,scales = "free", ncol = 10)+
  scale_x_reverse()
print(figureS20.ins)

## fig S22
figureS22.glu <- plot_cells(cds, genes=c("Count.glu"),
                            show_trajectory_graph = FALSE,
                            label_cell_groups = F,
                            label_groups_by_cluster = F,
                            label_leaves=FALSE,
                            label_branch_points=FALSE,
                            label_roots = FALSE)+
  facet_wrap(~Donortime,scales = "free", ncol = 10)+
  scale_x_reverse()
print(figureS22.glu)


