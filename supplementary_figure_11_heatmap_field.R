# Heatmap for the supplementary files
library(ampvis2)
library(tidyverse)
#Load data
comb_amp<-readRDS("total_ampvis.rds")
comb_amp$metadata
# Heatmap of the field
heatmap_over_times<-comb_amp %>% 
  amp_subset_samples(Treatment %in% "Field") %>%
  amp_subset_samples(!Days %in% "0") %>%
  amp_heatmap(tax_show = 15,
              tax_aggregate = "Genus",
              facet_by = "Variety",
              group_by ="Days")

ggsave("Plots/heatmap_supp_field.png",heatmap_over_times, width = 6, height = 5)
