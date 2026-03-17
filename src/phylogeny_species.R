# ==============================================================
#                   PHYLOGENY ANALYSIS - SPECIES
#                CLUSTER AND NUMBER OF INTAGRATION
# ==============================================================


# Load necessary libraries for data manipulation, plotting, and phylogenetic analysis
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)
library(ggh4x)
library(aplot)
library(ggdendro)
library(glue)
library(data.table)

#Path directories
dir_In_processed <- "data/processed/"
dir_Output <- "output/"
dir_Input <- "data/raw/"
dir_fig <- "figures/"

# Load phylogenetic tree from file
myTree <- ape::read.tree(paste0(dir_Input,"Insect_phylogeny.treefile"))
# Load cluster data (e.g. sequences assigned to clusters)
#clust <- fread(paste0(dir_In_processed,"clustered_hit_IV27_J2_0.95_0.7_cluster.tsv"), header=FALSE)
# Load the average number of integrations per host based on circles
nb_hit_genome_sp <- fread(paste0(dir_Output,"nb_hit_genome_sp.txt"), header=TRUE)
# Load the final table of hits representing integrations in hosts genome
dt_hits_cons_assi <- fread(paste0(dir_Output,"dt_integration_final.txt"), header=TRUE)
# Load PDV names with wasp species IV and BV
PDV <- fread(paste0(dir_Input,"PDV_wasp_sp.txt")) %>% distinct(segment, .keep_all = TRUE)

# Root the tree using a specific outgroup (Allacma_fusca)
rooted.tree <- root(myTree, which(myTree$tip.label %in% c("Allacma_fusca")))

# Remove specific species from the tree (not needed for this analysis)
myTree_2 <- drop.tip(rooted.tree, c("Campodea_augens", "Thermobia_domestica", "Ischnura_elegans"))

# Plot the tree with node labels for visual inspection
plot_test <- ggtree(myTree_2) + geom_text(aes(label=node), hjust=-.3) + theme_tree2()

# Create a metadata dataframe from the tree tips
metadata <- data.frame(genspe=rooted.tree$tip.label) %>%
  # Replace species name for consistency
  mutate(genspe=if_else(genspe=="Lochmaea_capreae", "Lochmaea_crataegi", genspe)) %>%
  # Split genus and species
  separate(genspe, into=c('genus','sp'), sep="_", remove=FALSE) %>%
  # Create formatted label for italic genus/species names
  mutate(lab = glue("italic({genus})~italic({sp})")) %>%
  # Join with external metadata (e.g., hits per genome)
  left_join(nb_hit_genome_sp, by=c("genspe")) %>%
  # Flag whether the species has a defined family
  mutate(TH=if_else(is.na(family), "n", "y"))

# Create a dataframe indicating which nodes represent main insect orders for highlighting
d <- data.frame(node=c(92, 71, 127, 121), type=c("Coleoptera", "Hymenoptera", "Phasmoptera", "Orthoptera"))

# Join metadata with tree tips (used later for coloring tips)
d1 <- data.frame(genspe = rep(myTree$tip.label)) %>%
  left_join(nb_hit_genome_sp, by=c("genspe")) %>%
  mutate(TH=if_else(is.na(family), "n", "y"))

# Same as above, but without the TH column
d2 <- data.frame(genspe = rep(myTree$tip.label)) %>%
  left_join(nb_hit_genome_sp, by=c("genspe"))

# Build a phylogenetic tree plot
ptree <- ggtree(myTree_2) %<+% metadata +
  geom_tiplab(aes(label=lab, color = TH), parse=T, align=TRUE, size=2.7) +
  #geom_tippoint(aes(color=TH), size=2, alpha=.6) +
  scale_color_manual(values=c("black","brown3"))+
  geom_hilight(data=d, aes(node=node, fill=type), align="left", type = "rect") + 
  scale_fill_manual(values=c("plum3", "lightpink","bisque", "lightsalmon")) +
  geom_cladelabel(node=92, label=c("Coleoptera"), fontsize = 3.6 ,barsize=NA, hjust=2.3, vjust=7, color="plum4")+
  geom_cladelabel(node=71, label=c("Hymenoptera"),fontsize = 3.6, barsize=NA, hjust=1.6, vjust=5.8, color="lightpink4")+
  geom_cladelabel(node=127, label=c("Phasmatodea"), fontsize = 3.6 ,barsize=NA, hjust=1.4, vjust=-1.7, color="lightsalmon3")+
  geom_cladelabel(node=121, label=c("Orthoptera"),fontsize = 3.6, barsize=NA, hjust=1.8, vjust=-2.5, color="bisque4")+
  xlim_tree(2) + theme(legend.position="none", plot.title =element_text(size=9)) + 
  ggtitle("A")

# Plot bar chart: mean number of integrations per genome
p2 <- ggplot(d2, aes(y=genspe, x=mean_nb_hits, fill=circles)) + 
  geom_col() +
    scale_fill_manual(values=c( "#023858","#045a8d","#0570b0","#3690c0","#74a9cf","#a6bddb","#d0d1e6","#ece7f2", 
                                         "#fff7f3","#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1","#dd3497","#7a0177","#49006a"
                                )) +
#  scale_fill_manual(values=c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1",
#                                      "#c6dbef", "#deebf7","#f7fbff","#f7fcb9","#d9f0a3","#addd8e", 
#                                      "#78c679", "#41ab5d", "#238443","#006837","#004529")) +
                                        xlim(0, 150) + 
  theme_classic() +
  theme(plot.title =element_text(size=9),
        axis.title.y=element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y =element_blank(),
        axis.ticks.x =element_line(),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=9), 
        axis.line.x.bottom = element_line(colour = "black"), 
        axis.line.y  = element_blank()) +
  ggtitle("C") +
  xlab("Mean number of PDV \n integrations per genome")

# ==============================================================
#                         MAIN FUNCTION
# ==============================================================

create_phylogenetic_analysis <- function(clust_filename, output_suffix) {
  
  cat("\n=== Processing:", clust_filename, "===\n")
  
  # Read cluster file
  clust <- fread(paste0(dir_In_processed, clust_filename), header=FALSE)
  clust <- setNames(clust, c('cluster_name', 'seq_name'))
  
# Prepare and format integration data
dt_hit_CLUS <- filter(dt_hits_cons_assi, order!="Lepidoptera") %>%
  mutate(FinalStart=round(FinalStart, digits = -3), FinalEnd=round(FinalEnd, digits = -3)) %>% 
  separate(query_id , into = c("query_id", "lepi"), sep="-") %>%
  distinct(subject_id, genspe, query_id , FinalStart, FinalEnd)

# Extract start and end positions from cluster sequences
clust_1 <- separate(clust, seq_name, into=c("ref_seq", "plus"), sep=":", remove = FALSE) %>% 
  separate(plus, into=c("start", "end"), sep="-", remove = TRUE) %>%
  mutate(FinalStart=round(as.numeric(start), digits = -3),
         FinalEnd=round(as.numeric(end), digits = -3)) %>%
  select(-start, -end)

stat_clust <- group_by(clust_1, cluster_name) %>% tally() %>% ungroup() %>% 
  summarise(mean_clus=mean(n), max_clus=max(n), min_clus=min(n), median_clus=median(n))


# Join cluster and hit data
clust_2 <- left_join(dt_hit_CLUS, clust_1, by=c('subject_id'='ref_seq', "FinalEnd"))

# Count number of hits per cluster and species
family_clust <- filter(clust_2, !is.na(cluster_name)) %>%
  group_by(cluster_name, genspe, query_id ) %>%
  tally() %>% ungroup()

# Reconstruct cluster table across tree species
d6 <- data.frame(genspe= rep(myTree$tip.label)) %>%
  left_join(family_clust, by=c("genspe"))

# Format for heatmap: log transform, pivot, convert to matrix
df_1 <- mutate(d6, n=log(n+1)) %>%
  filter(!is.na(d6$query_id )) %>%
  pivot_wider(id_cols = cluster_name, names_from = genspe, values_from = n , values_fn = sum, values_fill = as.numeric(0)) %>%
  mutate_if(is.integer,as.numeric) %>%
  column_to_rownames("cluster_name")

# Create dendrogram for clustering clusters
otter_dendro <- as.dendrogram(hclust(d = dist(x = data.matrix(df_1)), method="ward.D"))
dendro_plot <- ggdendrogram(data = otter_dendro) + 
  theme(axis.text.y =element_blank(),
        axis.text.x =element_blank(),
        plot.margin = unit(c(0.2,0,0.2,0.2), "cm"))

# Save dendrogram as SVG
ggsave(paste0(dir_fig, "dendro_plot_", output_suffix,".svg"), 
       dendro_plot,
       height = 1.7, 
       width = 10)

# Extract order of clusters from dendrogram
otter_order <- order.dendrogram(otter_dendro)

# Reorder cluster levels accordingly
d6$cluster_name <- factor(x = d6$cluster_name,
                          levels = tibble::rownames_to_column(df_1, "cluster_name")$cluster_name[otter_order], 
                          ordered = TRUE)

# Create heatmap of integration counts per cluster and species
p4 <- ggplot(d6, aes(cluster_name, genspe, fill= n)) + 
  geom_tile() +
  scale_fill_gradient(low="gray87", high="black") +
  theme_classic() +
  theme(plot.title =element_text(size=9),
        legend.title = element_text(size = 10),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y =element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        axis.line.y  = element_blank()) +
  ggtitle("B") + labs(fill = "Number of \n integrations \n per cluster") +
  xlab("Cluster") 

# Load wasp segment data and match to clusters


d7 <- d6 %>% filter(!is.na(query_id)) %>% 
  left_join(PDV, by=c("query_id"="segment")) %>%
  mutate(y=1)

# Create tile plot to show which wasp species contributed to each cluster
p6 <- ggplot(d7) + 
  geom_tile(aes(cluster_name, y, fill= wasp_species)) +
  theme_void() +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 10),
        axis.ticks.x =element_blank(),
        axis.title.y=element_blank(),
        axis.text.y =element_blank(),
        axis.ticks.y =element_blank(),
        legend.key.size = unit(0.3, 'cm')) +
  labs(x = "Cluster") + 
  labs(fill = "Wasp species") 

# Combine main plots: tree, heatmap, bar plot
p <- p4 %>% insert_left(ptree, width=1.7) %>% insert_right(p2, width=0.8)

# Save figure as SVG
output_file <- paste0(dir_fig, "phylogenie_species_", output_suffix, ".svg")
ggsave(output_file, 
       p,
       height = 7, 
       width = 10)

cat("✓ Saved:", output_file, "\n")

cat("Statistic:", output_file, "\n")
return(stat_clust)
}

# ==============================================================
#                    EXECUTE ANALYSES
# ==============================================================

# Process each cluster file
create_phylogenetic_analysis("clustered_hit_J1_0.8_0.7_cluster.tsv", "J1_0.8_0.7")
create_phylogenetic_analysis("clustered_hit_J2_0.8_0.7_cluster.tsv", "J2_0.8_0.7")
create_phylogenetic_analysis("clustered_hit_IV27_J2_0.95_0.7_cluster.tsv", "J2_IV27_0.95_0.7")
create_phylogenetic_analysis("clustered_hit_IV27_J1_0.95_0.7_cluster.tsv", "J1_IV27_0.95_0.7")
create_phylogenetic_analysis("clustered_hit_J1_0.95_0.7_cluster.tsv", "J1_0.95_0.7")
create_phylogenetic_analysis("clustered_hit_J2_0.95_0.7_cluster.tsv", "J2_0.95_0.7")

cat("\n=== All analyses complete! ===\n")
