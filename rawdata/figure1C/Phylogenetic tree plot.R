if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

library(tidyverse)
library(ggtree)
library(treeio)
library(cowplot)

nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

setwd("C:/Users/cf466/OneDrive - University of Exeter/Documents/Documents - CF laptop/PhD/Phage Cocktails/Host Genomes (Main)")

nwk <- read.newick("core-genes-full.nwk", node.label='support')
df_tip_data <- read.csv("panel_strains.csv") %>% 
  mutate(Origin = ifelse(panel == "Y", "Panel Strain", "This Study")) %>% 
  select(-panel)

p1 <- ggtree(nwk, branch.length = "none") %<+% df_tip_data
p2 <- rotate(p1, 64) +
  geom_tiplab(size=6, aes(color=Origin)) +
  scale_colour_manual(values = c("#fdb462", "#a35eb5")) +
  #geom_nodepoint(color="black", size=0.1) +
  geom_nodelab(geom='label', aes(label=support, subset=support < 0.9), hjust = 1.8, size = 4) +
  geom_hilight(node=69, fill="steelblue", alpha=.4) +
  labs(caption = "Confidence values only shown if less than 0.9") +
  geom_hilight(node=c(109,85, 82, 83), fill="darkgreen", to.bottom = TRUE, alpha=.4) +
  hexpand(.04, direction = 1) +
  theme(legend.position = "bottom")
p2

scale_plot <- ggtree(nwk, layout = "equal_angle") %<+% df_tip_data +
  geom_treescale(offset = 0.0005, width = 0.01, linesize = 1) +
  geom_tippoint(aes(color=Origin), shape=16, size=3) +
  scale_colour_manual(values = c("#fdb462", "#a35eb5")) +
  geom_hilight(node=69, fill="steelblue", alpha=.4) +
  geom_hilight(node=c(109,85,83,82), fill="darkgreen", to.bottom = TRUE, alpha=.4) +
  theme(plot.margin = grid::unit(c(-100, 15, -100, 15), "mm"))
scale_plot

plot_grid(p2, scale_plot, labels = "AUTO", nrow = 2, ncol = 1, relwidths = c(1,1))
