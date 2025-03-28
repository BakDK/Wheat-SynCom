#Plotting figure 4

# calculate the dissimilarity (Bray Curtis) between samples
library(ampvis2)
library(tidyverse)
library(vegan)
library(ggplot2)


combi_amp<-readRDS("total_ampvis.rds")

# Samples are rows
synnat<-combi_amp %>% amp_subset_samples(Treatment %in% c("SynCom","Control")) %>% amp_subset_samples(Compartment %in% "Rhizoplane") 

synnat_df<-t(synnat$abund)

# Calculate disssimilarity after rarefying
SNrare_dist <-avgdist(synnat_df,dmethod = "bray", sample= min(colSums(synnat$abund))) # min is 20,040, iterations = 100 default

SNrare_dist_tibble<-SNrare_dist %>%
  as.matrix() %>%
  as_tibble(rownames = "Sample") %>%
  pivot_longer(-Sample) %>%
  filter(name< Sample)

# Use dissimilarities out and  make a plot comparing only within day and treatment
# get the metadata first
synnat_met<-synnat$metadata %>% as_tibble() %>% select(-Medium,-Compartment,-Batch,-Pot_number)

synnat_rare_dist<-inner_join(SNrare_dist_tibble, synnat_met, by = c("Sample" ="Sample_ID")) %>%
  inner_join(.,synnat_met, by = c("name" = "Sample_ID")) %>%
  mutate(TT_sam =paste(Treatment.x,Timepoint.x,sep = "_"), TT_name =paste(Treatment.y,Timepoint.y,sep = "_")) %>%
  group_by(TT_sam) %>%
  filter(TT_sam == TT_name) %>%
  ungroup()

#Order factors
synnat_rare_dist$Treatment.x<-factor(synnat_rare_dist$Treatment.x, levels = c("SynCom","Control"))

# Plot the NatCom and SynCom dissimilarities

SynNat_syn_bray<-synnat_rare_dist %>% filter(Treatment.x %in% c("Control","SynCom")) %>%
  ggplot(aes(x = as.factor(Days.x), y = value, color = Treatment.x))+ geom_point(position = position_dodge(width = 0.75)) +
  geom_boxplot(aes(fill = Treatment.x),alpha = 0.1)+
  scale_color_manual(values=c("#1ab1cc","#725d4c"),labels = c("SynCom", "NatCom"))+
  scale_fill_manual(values=c("#1ab1cc","#725d4c"),labels = c("SynCom","NatCom"))+
  ylab("Bray Curtis dissimilarity \n among replicates")+xlab("Days after sowing")+theme_bw()+
  labs(color ="Treatment", fill = "Treatment")+theme(legend.position = "bottom")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#PCoA only NatCom and SynCom

pca_SN_rare<-ecodist::pco(SNrare_dist)
# Determine how much the first two Principal components explain

PC_SN<-sum(pca_SN_rare$values)
PC1_SN<-pca_SN_rare$values[1]/PC_SN*100 # 25.1%
PC2_SN<-pca_SN_rare$values[2]/PC_SN*100 # 14.8%
SN_bray_curtis_pcoa_df <- data.frame(pcoa1 = pca_SN_rare$vectors[,1], 
                                  pcoa2 = pca_SN_rare$vectors[,2])

SN_bray_curtis_pcoa_df$samples<-rownames(pca_SN_rare$vectors)

#Order factors
synnat_met$Treatment<-factor(synnat_met$Treatment, levels = c("SynCom","Control"))

# Create a plot
PCoA_pl_SN<-inner_join(SN_bray_curtis_pcoa_df,synnat_met, by = c("samples" ="Sample_ID"))%>%
  ggplot(aes(x=pcoa1, y=pcoa2)) + geom_point(aes(shape = Timepoint, color = Treatment))+
  labs(x = "PCoA 1 25%" ,
       y = "PCoA 2 14%" ) +
  theme(title = element_text(size = 10)) +
  scale_color_manual(values=c("#1ab1cc","#725d4c"),labels = c("SynCom","NatCom"), name = "Inoculum")+
  scale_shape_discrete(labels=c("7","14","21","28"))+labs(shape = "Days")

#Improving aesthetics
PCoA_SN_2<-PCoA_pl_SN + theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(-1,-1,-1,-1))

SynNat_syn_bray+theme(legend.position = "none")
# Plot together
library(ggpubr)

legende_til_plot<-get_legend(PCoA_SN_2)
falles_plot<-ggarrange(SynNat_syn_bray+theme(legend.position = "none"),PCoA_SN_2+theme(legend.position = "none"),
                       labels = c("A","B"))
ggarrange(falles_plot,as_ggplot(legende_til_plot), ncol =1, heights = c(9,2))

# Plot figure 4
x11(width =8, height = 5)
ggarrange(falles_plot,as_ggplot(legende_til_plot), ncol =1, heights = c(9,2))
ggsave("Plots/BC_dis_PCoA_comb.png",width = 8, height = 5, dpi = 300)
dev.off()


## Use only field samples to determine the effect size of each factor for driving community composition

Field_only <- combi_amp %>% amp_subset_samples(Treatment %in% "Field") %>% amp_subset_samples(Compartment %in% "Rhizoplane")
library(permute)
h1 <- with(data.frame(Field_only$metadata), how(nperm = 999, blocks = Plot)) # Accounting for block effect
set.seed(5)
adonis2(rare_dist_Field~Variety*Days, data = data.frame(Field_only$metadata), by = "terms", permutations = h1)


