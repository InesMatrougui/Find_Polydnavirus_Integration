# ==============================================================
#                   PHYLOGENY ANALYSIS - J1 & J2
# ==============================================================

# ----------------------
# Load Required Libraries
# ----------------------
library(ggplot2)
library(ggtree)
library(tidyr)
library(tidytree)
library(dplyr)
library(ape)
library(phytools)
library(data.table)
library(ggrepel)
library(ggtext)

# ----------------------
# Define Input/Output Paths
# ----------------------
dir_In_processed <- "data/processed/"
dir_Output <- "output/"
dir_Input <- "data/raw/"
dir_fig <- "figures/"

# ----------------------
# Load Integration Data
# ----------------------
dt_hits_cons_assi <- fread(paste0(dir_Output,"dt_integration_final.txt"), header = TRUE)

# ----------------------
# Function to Plot Phylogenetic Tree
# ----------------------
plot_phylogeny <- function(tree_filename, add_data, suffix, tip_colors, breaks, output_name) {
  
  # Read the tree
  myTree <- ape::read.tree(paste0(dir_In_processed, tree_filename))
  
  # Keep unique combinations of subject_id/order/family/genspe from integration data
  Family_order <- distinct(dt_hits_cons_assi, subject_id, .keep_all = TRUE) %>% 
    select(subject_id, order, family, genspe)
  
  # Bind additional data about parasitoid wasp cercle (e.g., Hd27, DsIV)
  Family_order <- rbind(Family_order, add_data)
  
  
  # Generate metadata from tip labels
  metadata <- data.frame(subject_id = myTree$tip.label) %>%
    separate(col = "subject_id", into = c("subject_id", "plus"), sep = "_", fill = "right") %>%
    mutate(plus = if_else(plus == "DsIV", "DsIV_15", plus)) %>%
    left_join(Family_order, by = "subject_id") %>%
    unite("subject_id_NA", c("subject_id", "plus"), sep = "_", remove = FALSE) %>%
    mutate(subject_id = if_else(is.na(plus), subject_id, subject_id_NA)) %>%
    select(subject_id, order, family, genspe)
  
  # Compute bootstrap values (only keep >70)
  bootstrap_vector <- c(rep(NA, ((length(myTree$edge) / 2) - length(myTree$node.label)) + 1),
                        ifelse(as.numeric(myTree$node.label) < 70, NA, as.numeric(myTree$node.label)))
  
  node_data <- data.frame(node = 1:length(bootstrap_vector),
                          has_bootstrap = !is.na(bootstrap_vector))
    
    p <- ggtree(myTree, branch.length = 'none', layout = "equal_angle", aes(color = genspe)) %<+% metadata %<+% node_data +
      #geom_tiplab(aes(color = genspe), size = 2, hjust = -0.1, ) +
      geom_nodepoint(aes(subset = has_bootstrap, fill = ">70"), 
                     size = 1.2,
                     shape = 21,
                     color = "grey50",
                     stroke = 0.3,
                     show.legend = TRUE) +
      coord_cartesian(clip = 'off') +
      geom_treescale() +
      theme_tree2(plot.margin = margin(80, 80, 80, 80)) +
     # geom_tippoint(aes(shape = order, color=genspe), size = 3, alpha=0.3) +
      theme(legend.position = "right") +
      scale_color_manual(name= "Species", breaks = breaks, values = tip_colors, 
                         labels = function(x)
                           ifelse(x == "parasitoid_wasp",
                                  gsub("_", " ", x),
                                  parse(text = paste0("italic('", gsub("_", " ", x), "')")))) +
      scale_fill_manual(name = "Bootstrap", 
                        values = c(">70" = "grey50"),
                        guide = guide_legend(override.aes = list(shape = 20, color = "grey50", stroke = 0.3))) +
      theme(legend.text = element_text(),
        legend.title = element_text()) 
    
    print(p) 
    
  # Save the plot
  ggsave(paste0(dir_fig, output_name, ".svg"), 
         plot = p, height = 8, width = 10)
}

# ----------------------
# Run for J1
# ----------------------
add_data_J1 <- data.frame(
  subject_id = c("Hd27", "KF156228.1", "MZ129252.1", "NC.008952.1"),
  order = rep("wasp", 4),
  family = rep("wasp", 4),
  genspe = rep("parasitoid_wasp", 4)
)

color_breaks_J1 <- c("parasitoid_wasp", "Bacillus_rossius", "Clitarchus_hookeri", "Dryococelus_australis", "Crioceris_asparagi", "Lochmaea_crataegi",
                     "Neodiprion_pinetum", "Neodiprion_lecontei", "Neodiprion_virginianus", "Neodiprion_fabricii", "Diprion_similis", "Abia_candens", "Athalia_rosae",
                     "Tenthredo_mesomela", "Tenthredo_livida", "Tenthredo_notha", "Rhogogaster_chlorosoma", "Macrophya_annulata", "Ametastegia_equiseti")

color_values_J1 <- c("black", "#c51b8a", "#fa9fb5","#fde0dd", "#d95f0e", "#fec44f",  
                      "#004a49", "#005e5d","#126b69", "#257675", "#3e8685", "#579695",
                      "#70a6a5", "#89b6b5", "#9ec4c3", "#b4d2d1", "#c4dddc", "#d4e7e6", "#f0f7f6", "#f0f7f9")
                            
plot_phylogeny(
  tree_filename = "J1_IV27_alignment_trimmed_unique.fasta.treefile",
  add_data = add_data_J1,
  suffix = "J1",
  tip_colors = color_values_J1,
  breaks = color_breaks_J1,
  output_name = "phylogeny_J1_Hd27_Ds15_b"
)

# ----------------------
# Run for J2
# ----------------------
add_data_J2 <- data.frame(
  subject_id = c("Hd27", "KF156228.1", "MZ129252.1", "NC.008952.1"),
  order = rep("wasp", 4),
  family = rep("wasp", 4),
  genspe = rep("parasitoid_wasp", 4)
)

color_breaks_J2 <- c("parasitoid_wasp", "Bacillus_rossius", "Clitarchus_hookeri", "Dryococelus_australis",  "Lochmaea_crataegi",
                     "Neodiprion_pinetum", "Neodiprion_lecontei", "Neodiprion_virginianus", "Neodiprion_fabricii", "Diprion_similis", "Abia_candens", "Athalia_rosae",
                     "Tenthredo_mesomela", "Tenthredo_notha", "Rhogogaster_chlorosoma", "Macrophya_annulata", "Macrophya_alboannulata" )

color_values_J2 <- c("black", "#c51b8a", "#fa9fb5","#fde0dd", "#d95f0e",  
                            "#004a49", "#005e5d","#126b69", "#257675", "#3e8685", "#579695",
                            "#70a6a5", "#89b6b5", "#9ec4c3", "#b4d2d1", "#c4dddc", "#d4e7e6", "#f0f7f6")
                            
plot_phylogeny(
  tree_filename = "J2_IV27_alignment_trimmed_unique.fasta.treefile",
  add_data = add_data_J2,
  suffix = "J2",
  tip_colors = color_values_J2,
  breaks = color_breaks_J2,
  output_name = "phylogeny_J2_Hd27_Ds15_b"
)

# ==============================================================

