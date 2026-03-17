# ==========================================================================================================
#                   PHYLOGENY ANALYSIS - Genomic Integrations at the Phylum, Order, and Family
# ==========================================================================================================

# ----------------------
# Load Required Libraries
# ----------------------
library(ggplot2)        # For creating plots
library(ggtree)         # For visualizing phylogenetic trees
library(ape)            # For handling phylogenetic trees
library(dplyr)          # For data manipulation
library(tidyr)          # For data tidying
library(data.table)     # For fast data loading and manipulation
require(treeio)         # For importing tree data
library(ggtreeExtra)    # For enhanced ggtree plots (e.g. adding extra data layers)
library(aplot)          # For insert an associated plot of a main plot

# ----------------------
# Define Input/Output Paths
# ----------------------
dir_In_processed <- "data/processed/"
dir_fig <- "figures/"
dir_Output <- "output/"

# ----------------------
# Load Integration Data
# ----------------------
genome_taxo_all <- fread(paste0(dir_Output, "genome_taxo_all.txt"), header = TRUE)        # Genomic and taxonomic data
dt_hits_cons_assi <- fread(paste0(dir_Output, "dt_integration_final.txt"), header = TRUE) # Integration results data

# ==============================================================

### PHYLUM-LEVEL PHYLOGENY ###

# Load the phylum phylogenetic tree created using TimeTree (https://timetree.org/)
myTree_p <- ape::read.tree(paste0(dir_In_processed, "phylo_phylum.nwk"))

# Drop unwanted taxa from the tree
drop_taxa <- c("Porifera", "Ctenophora", "Hemichordata", "Echinodermata", "Xenacoelomorpha", "Chordata", 
               "Nematoda-2", "Chaetognatha", "Brachiopoda-2", "Brachiopoda-3", "Sipuncula", 
               "Entoprocta", "Cnidaria", "Placozoa")
myTree_p <- drop.tip(myTree_p, drop_taxa)

# Create the base tree plot
tree <- ggtree(myTree_p) + geom_tiplab(size=7) + xlim_tree(1000) + theme_tree2() +
  scale_x_continuous(breaks = seq(-700, 0, by = 350))

ptree <- revts(tree) +
  xlab("Mya") + theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=22))
# Count number of unique species per phylum in both datasets
number_sp <- distinct(genome_taxo_all, genspe, .keep_all = TRUE) %>% group_by(phylum) %>% tally()
number_sp_int <- distinct(dt_hits_cons_assi, genspe, .keep_all = TRUE) %>% group_by(phylum) %>% summarise(n_int = n())

# Merge counts with the tree tip labels
d1 <- data.frame(phylum = myTree_p$tip.label) %>% 
  left_join(number_sp, by = "phylum") %>% 
  left_join(number_sp_int, by = "phylum")

# Barplot showing number of species with and without integration per phylum
p1 <- ggplot(d1, aes(y = phylum, x = n, fill = 'Without integration')) + 
  geom_col() +
  geom_col(aes(y = phylum, x = n_int, fill = 'With integration')) +
  scale_fill_manual(values = c('firebrick', 'lightblue')) +
  theme_classic() +
  theme(
    plot.title = element_text(size=9),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size=22),
    axis.text.x = element_text(size=15), 
    axis.line.x.bottom = element_line(colour = "black"), 
    axis.line.y = element_blank(), 
    legend.text = element_text(size=14)
  ) +
  xlab("Number of species")

# Combine the barplot with the tree
p <- p1 %>% insert_left(ptree, width = 0.7)

# Save the phylum-level plot
ggsave(paste0(dir_fig, "phylogeny_phylum.svg"), p, height = 7, width = 10)

# ==============================================================

### INSECT ORDER-LEVEL PHYLOGENY ###

# Load the insect order tree (also from TimeTree)
myTree_i <- ape::read.tree(paste0(dir_In_processed, "phylo_insecct_orders.nwk"))

# Drop redundant or placeholder taxa
drop_taxa_i <- c("Zoraptera", "Psocoptera-2", "Phthiraptera-2", "Psocoptera-3", "Psocoptera-4", "Psocoptera-5", 
                 "Psocoptera-6", "Psocoptera-7", "Psocoptera-8", "Psocoptera-9", "Phthiraptera-3", "Phthiraptera-4", 
                 "Phthiraptera-5", "Phthiraptera-6", "Phthiraptera-7", "Phthiraptera-8", "Phthiraptera-9", 
                 "Phthiraptera-10", "Phthiraptera-11", "Phthiraptera-12", "Phthiraptera-13")
myTree_i <- drop.tip(myTree_i, drop_taxa_i)

# Plot the tree
tree <- ggtree(myTree_i) + geom_tiplab(size=6) + xlim_tree(600) + theme_tree2() +
  scale_x_continuous(breaks = seq(-400, 0, by = 200)) 
           
  
ptree <- revts(tree) +
  xlab("Mya") + theme(axis.text.x = element_text(size=15), axis.title.x = element_text(size=22))

# Count unique species per order (Insecta only)
number_sp <- filter(genome_taxo_all, class == "Insecta") %>% 
  distinct(genspe, .keep_all = TRUE) %>% group_by(order) %>% tally()
number_sp_int <- distinct(dt_hits_cons_assi, genspe, .keep_all = TRUE) %>% group_by(order) %>% summarise(n_int = n())

# Merge data with tree tips
d2 <- data.frame(order = myTree_i$tip.label) %>% 
  left_join(number_sp, by = "order") %>% 
  left_join(number_sp_int, by = "order")

# Barplot with species per order
p2 <- ggplot(d2, aes(y = order, x = n, fill = 'Without integration')) + 
  geom_col() +
  geom_col(aes(y = order, x = n_int, fill = 'With integration')) +
  scale_fill_manual(values = c('firebrick', 'lightblue')) +
  theme_classic() +
  theme(
    plot.title = element_text(size=9),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.title.x = element_text(size=22),
    axis.text.x = element_text(size=15),
    axis.line.x.bottom = element_line(colour = "black"),
    axis.line.y = element_blank(),
    legend.text = element_text(size=12)
  ) +
  xlab("Number of species")

# Combine with tree
p_i  <- p2 %>% insert_left(ptree, width = 0.8)

# Save insect order-level plot
ggsave(paste0(dir_fig, "phylogeny_insect_order.svg"), p_i, height = 7, width = 10)

# ==============================================================

### FAMILY-LEVEL PHYLOGENY ###

# Read mapping table for species name corrections
equi_sp <- read.delim(paste0(dir_In_processed, "equi_family_for_R.txt"), header = FALSE)
equi_sp_1 <- unite(equi_sp, "old_sp", c("V1", "V2"), sep = "_") %>% 
  mutate(V4 = if_else(V4 == "sp.", "", V4)) %>% 
  unite("new_sp", c("V3", "V4"), sep = " ", remove = FALSE) %>% 
  mutate(new_sp = if_else(V4 == "", V3, new_sp))

# filter to keep only orders with integration 
sp <- filter(genome_taxo_all, order=="Hymenoptera" | order=="Lepidoptera" | order=="Orthoptera" | 
               order=="Coleoptera" | order=="Phasmatodea") %>% 
  distinct(family, .keep_all = TRUE) %>% select(order, family, genspe)

# Replace old species names with updated ones
sp_plus <- left_join(sp, equi_sp_1, by = c("genspe" = "old_sp")) %>% 
  mutate(genspe = if_else(!is.na(new_sp), new_sp, genspe)) %>% 
  select(-new_sp, -V3, -V4)
sp_plus$genspe <- gsub(" ", "_", sp_plus$genspe)

# Load the final family-level tree
myTree <- ape::read.tree(paste0(dir_In_processed, "sp_family_final.nwk"))

# Merge species with metadata
d <- data.frame(genspe = myTree$tip.label) %>% left_join(sp_plus, by = "genspe")

# Define nodes to highlight (corresponding to insect orders)
d1 <- data.frame(node = c(253, 206, 162, 154, 156), type = c("Lepidoptera", "Coleoptera", "Hymenoptera", "Phasmoptera", "Orthoptera"))

# Count species per family (for selected orders)
number_sp <- filter(genome_taxo_all, order %in% c("Hymenoptera", "Lepidoptera", "Orthoptera", "Coleoptera", "Phasmatodea")) %>% 
  distinct(genspe, .keep_all = TRUE) %>% group_by(family) %>% tally()
number_sp_int <- distinct(dt_hits_cons_assi, genspe, .keep_all = TRUE) %>% group_by(family) %>% summarise(n_int = n())

# Merge counts and prepare data for plotting
d2 <- data.frame(genspe = myTree$tip.label) %>% 
  left_join(sp_plus, by = "genspe") %>% 
  left_join(number_sp, by = "family") %>% 
  left_join(number_sp_int, by = "family") %>%
  mutate(int = if_else(!is.na(n_int), "yes", "no")) %>%
  filter(family != "Ichneumonidae" | family != "Braconidea") %>%
  mutate(n_100 = if_else(n > 100, n - 100, NA)) %>%
  mutate(n = if_else(n > 100, 100, n),
         n = if_else(!is.na(n_int) & n >= n_int, n - n_int, n))

# Convert to long format for plotting
d3 <- d2 %>% pivot_longer(cols = c('n', 'n_int'), names_to = 'p', values_to = 'nb')

# Plot the phylogenetic tree in fan layout with colored points and highlighted clades
tree <- ggtree(myTree, layout = "fan", open.angle = 10) %<+% d2 + 
  geom_hilight(data = d1, aes(node = node, fill = type), type = "rect", alpha=0.75) +
  theme(legend.position = "none")

# Add barplots to the tree (number of species and integrations per family)
tree_f <- tree + 
  geom_fruit(data = d3, geom = geom_col, mapping = aes(y = genspe, x = nb, fill = p),
             pwidth = 1, orientation = "y", size = 0.02, offset = 0.02,
             axis.params = list(axis = "x", text.size = 1.8, hjust = 1, vjust = 1, nbreak = 4),
             grid.params = list()) +
  scale_fill_manual(values = c("#a9511a", "#467575", "#5b005bbf", 'lightblue', 'firebrick', "#cfa100", "#db5a76"))

tree_f
# Save the final tree
ggsave(paste0(dir_fig, "phylogenie_family.svg"), tree_f, height = 7, width = 10)
