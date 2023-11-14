
# indicator species ####

library(indicspecies)
# community matrix as df with sites in rows
t(asv_table_r) %>% as.data.frame() -> sp.df
# make a vector with sampling sites
sp.df %>% rownames_to_column("sample_name") %>% 
  dplyr::select("sample_name") %>% 
  left_join(.,env_lobau_final_r, by = c("sample_name" = "sample_id")) %>% 
  dplyr::select("well_id.y") %>% 
  table()

c(rep(1, 17), rep(2, 23) ) -> group_vect


indval <- multipatt(sp.df, 
                    group_vect,
                    func="IndVal.g",  
                    duleg=TRUE,  
                    control = how(nperm=999)) 
indval$sign
# the NA p values are the species present in both sites





# adjust with fdr:

# see how many are significant before the fdr
indval[["sign"]][["p.value"]] %>% as.data.frame() %>% 
  filter(.<0.05) # 908

# after fdr
p.adjust(indval[["sign"]][["p.value"]], method="BH") %>% 
  as.data.frame() %>% 
  filter(.<0.05) # 597


# add species names 
p.adjust(indval[["sign"]][["p.value"]], method="BH") %>% 
  as.data.frame() %>% 
  dplyr::rename("adj.p" = ".") %>% 
  # after correcting add the names from the original output  
  cbind(indval$sign) %>% 
  # filter now only the truly sign ones
  dplyr::filter(adj.p < 0.05) %>% 
  dplyr::select("adj.p", "s.1", "s.2", "stat") %>% 
  
  rownames_to_column("ASVs") %>% 
  left_join(., tax_table %>% rownames_to_column("ASVs"), 
            by = "ASVs") -> species.df







# plot df ####

## bubble plot ####

# I remove the adj p because I dont use it for plotting.
```{r}
species.df %>% 
  mutate (Sample_site = ifelse(s.1 == 1, "D15", "ESW") ) %>% 
  dplyr::select("Division", "Family", "Sample_site", "stat") %>%
  dplyr::filter(Family != "") %>% 
  arrange(desc(stat) ) %>% 
  slice_max(stat, n= 70) %>% 
  ggplot() +
  geom_point(aes(x = as.factor(Sample_site), y= Family, 
                 color = Division, cex = stat )
             #shape = 21, color = "black"
  ) + 
  theme_bw() + 
  scale_color_manual(values = 
                       #distinct_palette(n = NA, pal = "brewerPlus") 
                       c("#A6CEE3",  "#1F78B4",  "#B2DF8A", "#33A02C",  "#E31A1C",   "#FDBF6F",
                         "#CAB2D6",   "#6A3D9A",  "#1ff8ff",  "#83d5af" ,  "#E7298A",
                         "#66A61E",  "#A6761D",  "#4b6a53",   "#b249d5",    "#7edc45",
                         "#5c47b8",   "#cfd251",   "#ff69b4", "#69c86c")) -> plot_ind

plot_ind

#"#cd3e50"    "#da6130"   "#5e79b2"   "#c29545"   "#532a5a"   "#5f7b35"  
#"#c497cf"   "#773a27"   "#7cb9cb"   "#594e50"   "#d3c4a8"   "#c17e7f" ))




## lollipop plot ####

species.df %>% 
  mutate (Sample_site = ifelse(s.1 == 1, "D15", "ESW") ) %>% 
  dplyr::select("Division", "Family", "Sample_site", "stat") %>%
  dplyr::filter(Family != "") %>% 
  arrange(desc(stat) ) %>% 
  slice_max(stat, n= 70) %>% 
  ggplot() +
  geom_point(aes(x = stat, y= Family, 
                 color = Sample_site),
             #shape = 21, color = "black"
             size = 3 ) + 
  geom_segment(aes(x=0.8, xend= stat, y=Family, yend =Family)) +
  theme_bw()


## try to cut sp down more  ####

# This function allows reducing drastically the number of species combinations to be retained for a given target site group.

sc <- indicators(X=sp.df, cluster=group_vect, 
                 group=2, 
                 max.order = 1, 
                 verbose=TRUE, 
                 At=0.5, Bt=0.2)


pruneindicators(sc, At = 0.2,
                Bt = 0.3,
                sqrtIVt = 0,
                alpha = 1,
                max.indicators = 4,
                verbose = FALSE)
