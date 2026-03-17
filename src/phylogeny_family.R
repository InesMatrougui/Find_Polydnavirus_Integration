setwd("C:/Users/inesm/Documents/cle_usb/Lecteur USB/PDV_arthopodes/14-02-25")

library("dplyr")
library("tidyr")
library("data.table")
require(treeio)
require(ggtree)
library(ggtreeExtra)


sp <- filter(genome_taxo_all, order=="Hymenoptera" | order=="Lepidoptera" | order=="Orthoptera" | 
               order=="Coleoptera" | order=="Phasmatodea") %>% 
  distinct(family, .keep_all = TRUE) %>% select(order, family, genspe)

supp_sp <- read.delim("equi_family_sp_supp.txt", header=FALSE)
supp_sp_1 <- unite(supp_sp, "sp", c("V1", "V2"))

old_sp <- read.table("test_family_red.txt", header=FALSE)

old_sp_1 <- filter(old_sp, !(V1 %in% supp_sp_1$sp))

new_sp <- read.delim("sp_family.txt", header=FALSE)

sp_new_old <- rbind(new_sp, old_sp_1)


#write.table(rbind(new_sp, old_sp_1), "sp_family_final.txt", quote = FALSE, row.names = F, col.names = F)

equi_sp <- read.delim("equi_family_for_R.txt", header=FALSE)
equi_sp_1 <- unite(equi_sp, "old_sp", c("V1","V2"), sep="_") %>% 
              mutate(V4=if_else(V4=="sp.", "", V4)) %>% 
              unite("new_sp", c("V3", "V4"), sep=" ", remove = FALSE) %>% 
                mutate(new_sp=if_else(V4=="", V3, new_sp))


sp_plus <- left_join(sp, equi_sp_1, by=c("genspe"="old_sp")) %>% 
            mutate(genspe=if_else(!is.na(new_sp), new_sp, genspe)) %>% 
              select(-new_sp, -V3, -V4)

sp_plus$genspe=gsub(pattern = " ", replacement = "_", x = sp_plus$genspe)


myTree <- ape::read.tree("sp_family_final.nwk")

d <- data.frame(genspe=myTree$tip.label) %>% 
  left_join(sp_plus, by="genspe")
d1 <- data.frame(node=c(253, 206, 162, 154, 156), type=c("Lepidoptera", "Coleoptera", "Hymenoptera", "Phasmoptera", "Orthoptera"))


number_sp <- filter(genome_taxo_all, order=="Hymenoptera" | order=="Lepidoptera" | order=="Orthoptera" | 
                      order=="Coleoptera" | order=="Phasmatodea") %>% distinct(genspe, .keep_all = T) %>% group_by(family) %>% tally() 
number_sp_int <- distinct(dt_hits_cons_assi, genspe, .keep_all = T) %>% group_by(family) %>% summarise(n_int= n())

d2 <- data.frame(genspe=myTree$tip.label) %>% left_join(sp_plus, by="genspe") %>% 
  left_join(number_sp, by="family") %>% left_join(number_sp_int, by="family") %>% 
  mutate(int=if_else(!is.na(n_int),"yes", "no")) %>% 
  filter(family!="Ichneumonidae" | family!="Braconidea") %>%
  mutate(n_100=if_else(n>100, n-100, NA)) %>% 
  mutate(n=if_else(n>100,100, n),
         n=if_else(!is.na(n_int) & n>=n_int ,n-n_int, n))
d3 <- d2 %>% pivot_longer(cols=c('n', 'n_int'), names_to='p', values_to='nb')


tree <- ggtree(myTree, layout="fan", open.angle=10) %<+% d2 + geom_tippoint(mapping=aes(color=int), 
                             size=1) + scale_color_manual(values=c("lightgrey","red"))+ 
  geom_hilight(data=d1, aes(node=node, fill=type), type = "rect") + theme(legend.position="none") 

tree_f <-  tree +  geom_fruit(data=d3, geom=geom_col, mapping=aes(y=genspe, x=nb, fill=p), pwidth=1, orientation="y", size=0.02, offset = 0.04,
                                                         axis.params=list(axis= "x",text.size= 1.8,hjust= 1,vjust= 1,nbreak= 4),
                                                         grid.params=list()) +
#   geom_fruit(data=d2, geom=geom_col, mapping=aes(y=genspe, x=n_100), fill="lightblue", pwidth=1, orientation="y", size=0.02, offset = 0.04,
#                axis.params=list(axis= "x",text.size= 1.8,hjust= 1,vjust= 1,nbreak= 4),
#                grid.params=list()) +
  scale_fill_manual(values=c("plum3", "lightpink", "cadetblue3",'lightblue','firebrick',"bisque", "lightsalmon"))

tree_f

ggsave("C:/Users/inesm/Documents/cle_usb/Lecteur USB/PDV_arthopodes/PDV_arthro/Figures/phylogenie_family_test_2.svg", 
       tree_f,
       height = 7, 
       width = 10)
