
asv_table_r %>% 
  mutate_if(is.integer, as.numeric) %>% 
  mutate(class = tax_table$Class) %>%
  pivot_longer(cols = -class, 
               values_to = "abundance", 
               names_to = "sample") %>%
  left_join(.,env_lobau_final_r, by = c("sample" = "sample_id")) %>%  
  dplyr::select ("abundance", "class", "well_id.y") %>%
  mutate(well_id.y = as.factor(well_id.y),
         class = as.factor(class)) %>%  
  
  # to find relative abundance within each category  
  dplyr::group_by(well_id.y) %>%
  mutate (rel_abund = abundance / sum(abundance)) %>%
  ungroup() %>%
  
  # to sum up all the counts of one class from all the samples  
  dplyr::group_by(well_id.y, class) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  ungroup() %>%
  
  # those very rare categorize as Diverse others 
  dplyr::group_by(well_id.y, class) %>%
  mutate(class = case_when(rel_abund < 0.02 ~ "Diverse others", 
                           TRUE ~ class)) %>%
  ungroup() %>%
  
  # final regrouping after categorizing for Diverse others  
  dplyr::group_by(well_id.y, class) %>%
  summarise(rel_abund = sum (rel_abund)) %>%
  ungroup() %>%
  
  left_join(., tax_table[,c("Division", "Class" )] %>% 
              distinct(Class, .keep_all=TRUE) , 
            by = c("class" = "Class") ) %>% 
  
  dplyr::mutate(active_comm = paste("DNA")) %>% 
  
  as.data.frame() -> data_DNA






asv_table_rna_r %>% 
  mutate_if(is.integer, as.numeric) %>% 
  mutate(class = tax_table$Class) %>%
  pivot_longer(cols = -class, 
               values_to = "abundance", 
               names_to = "sample") %>%
  left_join(.,env_lobau_final_rna_r, by = c("sample" = "sample_id")) %>%  
  dplyr::select ("abundance", "class", "well_id.y") %>%
  mutate(well_id.y = as.factor(well_id.y),
         class = as.factor(class)) %>%  
  
  # to find relative abundance within each category  
  dplyr::group_by(well_id.y) %>%
  mutate (rel_abund = abundance / sum(abundance)) %>%
  ungroup() %>%
  
  # to sum up all the counts of one class from all the samples  
  dplyr::group_by(well_id.y, class) %>%
  summarise(rel_abund = sum(rel_abund)) %>%
  ungroup() %>%
  
  # those very rare categorize as Diverse others 
  dplyr::group_by(well_id.y, class) %>%
  mutate(class = case_when(rel_abund < 0.02 ~ "Diverse others", 
                           TRUE ~ class)) %>%
  ungroup() %>%
  
  # final regrouping after categorizing for Diverse others  
  dplyr::group_by(well_id.y, class) %>%
  summarise(rel_abund = sum (rel_abund)) %>%
  ungroup() %>%
  
  left_join(., tax_table[,c("Division", "Class" )] %>% 
              distinct(Class, .keep_all=TRUE) , 
            by = c("class" = "Class") ) %>% 
  
  dplyr::mutate(active_comm = paste("RNA")) %>% 
  
  as.data.frame() -> data_RNA





rbind(data_DNA, data_RNA) %>% 
  mutate(class = case_when(class == "" ~ "Unclassified Ciliophora", 
                           TRUE ~ class) ) %>% 
  mutate(Division =  replace_na(Division, "Others") ) %>% 
  ggplot(aes( x=well_id.y, y= rel_abund, fill = class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~active_comm) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = distinct_palette(n = 26, 
                                              pal = "brewerPlus",
                                              add = "lightgrey"))


# Define the order of fill colors
rbind(data_DNA, data_RNA) %>% 
  mutate(class = case_when(class == "" ~ "Unclassified Ciliophora", 
                           TRUE ~ class) ) %>% 
  mutate(Division =  replace_na(Division, "Others") ) %>% 
  dplyr::select("Division") %>%  unique()


# rbind(data_DNA, data_RNA) %>% mutate(class = case_when(class == "" ~ "Unclassified Ciliophora", TRUE ~ class) ) %>%  mutate(Division =  replace_na(Division, "Others") ) %>% dplyr::filter(Division== "Ciliophora") %>%  dplyr::select("class") %>% unique()

fill_order <- c(#1. Ciliophora
                "Unclassified Ciliophora",
                # "Litostomatea",
                "Oligohymenophorea", 
                "Spirotrichea",
                # "Colpodea",
                # "Nassophorea", "Phyllopharyngea", 
  
                #2. Conosa,
                 "Archamoebea",
                
                
                #3. Ochrophyta,
                "Bacillariophyta" , "Chrysophyceae","Eustigmatophyceae",
                 "Dictyochophyceae", "Synurophyceae",
                
                
                #4. Centroheliozoa,
                # "Centroheliozoa_X",
                
                
                #5. Cryptophyta,
                "Cryptophyceae",
                
                
                #6. Dinoflagellata,
                "Dinophyceae",
                
                
                #7. Discoba,
                "Euglenozoa",
                
                
                #8. Cercozoa,
                 "Filosa-Sarcomonadea", "Filosa-Thecofilosea",
                
                
                #9. Katablepharidophyta,
                 "Katablepharidaceae",
                
                
                #10. Pseudofungi,
                 "Oomycota",
                
                
                #11. Opalozoa,
                 "Opalinata",
                
                
                #12. Chlorophyta,
                 "Chlorophyceae", 
                # "Trebouxiophyceae",
                
                
                #13. Haptophyta,
                # "Prymnesiophyceae",
                
                
                #14. Streptophyta
                # "Embryophyceae",
                
                
                "Diverse others"
  )






# Convert Category to a factor with a specific order
rbind(data_DNA, data_RNA) %>% 
  mutate(class = case_when(class == "" ~ "Unclassified Ciliophora", 
                           TRUE ~ class) ) %>% 
  mutate(Division =  replace_na(Division, "Others") ) -> dataset_bar
dataset_bar$class <- factor(dataset_bar$class, levels = fill_order)



rbind(data_DNA, data_RNA) %>% 
  mutate(class = case_when(class == "" ~ "Unclassified Ciliophora", 
                           TRUE ~ class) ) %>% 
  mutate(Division =  replace_na(Division, "Others") ) -> dataset_bar
dataset_bar$class <- factor(dataset_bar$class, levels = fill_order)

# Create a bar plot with ordered fill colors
ggplot(dataset_bar, aes( x=well_id.y, y= rel_abund, fill = class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~active_comm) +
  scale_fill_manual(values = c( 
"Unclassified Ciliophora" = "#990000"  , "Litostomatea" = "#CC3300",  
"Oligohymenophorea" = "#FF0000" , "Spirotrichea" = "#FF3333",  "Colpodea" =  "#FF6666", 
"Nassophorea" = "#FF9999"  , "Phyllopharyngea" =  "#FFCCCC",


"Archamoebea" =  "#FFCC00" ,


"Bacillariophyta" = "#009999", "Chrysophyceae" = "#00CCCC", "Eustigmatophyceae" = "#00FFFF", 
"Dictyochophyceae" = "#66FFFF", "Synurophyceae" = "#99FFFF" ,

#"#3366CC" 


"Centroheliozoa_X" = "#990099"  ,

"Cryptophyceae" ="#33CC33" ,

"Dinophyceae" = "#660009",

"Euglenozoa" = "#CC9933",


"Filosa-Sarcomonadea" = "#0066CC" , "Filosa-Thecofilosea" = "#66B2FF",  #"#669966",
#"#336633"


"Katablepharidaceae" = "#FF33FF",

"Oomycota" = "#c0d65c",

"Opalinata" = "#FF6600" ,


"Chlorophyceae" = "navyblue" , "Trebouxiophyceae" = "blue",


"Prymnesiophyceae" = "#666699", 

"Embryophyceae" = "#CC66CC" ,

"Diverse others" = "gray")) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 13, face = "bold"),
        strip.text = element_text(size =14)) +
  guides(fill = guide_legend(ncol = 1)) + 
  labs(y="Relative abundance %", x = NULL, fill = NULL) +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) -> comp_plot




# format legend


dataset_bar %>% dplyr::select("rel_abund", "well_id.y", "Division", "class") %>% 
  dplyr::rename("Abundance" = "rel_abund",
                "Sample"= "well_id.y",
                "Class" = "class") %>% 
  mutate(Top_Division = Division, 
          Top_Class = Class) %>% 
  mutate(group = paste(Top_Division, Top_Class, sep="-")) -> mdf_GP


hex <-c (

#"Unclassified Ciliophora"
 "#990000"  , 
#"Archamoebea" 
"#FFCC00" ,
#"Bacillariophyta"  
"#009999", 
#"Chrysophyceae"  
"#00CCCC", 
#"Cryptophyceae"    
"#33CC33" ,
#"Dinophyceae" 
"#660009",
#"Diverse others" 
"gray",
#"Euglenozoa"  
"#CC9933",
#"Eustigmatophyceae"  
"#00FFFF",
#"Filosa-Sarcomonadea"
"#0066CC" , 
#"Katablepharidaceae"  
"#FF33FF",
#"Oligohymenophorea"  
"#FF0000" , 
# "Opalinata"  
"#FF6600" ,
# "Spirotrichea"
 "#FF3333",  
# "Unclassified Ciliophora"
"#990000",  
# "Bacillariophyta" 
"#009999", 
# "Chrysophyceae"  
"#00CCCC", 
# "Cryptophyceae"  
"#33CC33" ,
# "Dictyochophyceae"
"#66FFFF", 
# "Dinophyceae"
"#660009",
# "Diverse others"   
"gray",
# "Filosa-Thecofilosea" 
 "#66B2FF",  
# "Oomycota"   
"#c0d65c",
# "Spirotrichea"  
 "#FF3333",  
# "Synurophyceae"  
"#99FFFF" ,
# "Unclassified Ciliophora"
"#990000",  
# "Bacillariophyta"  
"#009999", 
# "Chrysophyceae" 
"#00CCCC", 
# "Cryptophyceae" 
"#33CC33" ,
# "Dinophyceae" 
"#660009",
# "Diverse others" 
"gray",
# "Euglenozoa"
"#CC9933",
# "Eustigmatophyceae" 
"#00FFFF",
# "Filosa-Sarcomonadea"    
"#0066CC" , 
# "Katablepharidaceae"    
"#FF33FF",
# "Oligohymenophorea"    
"#FF0000" , 
# "Spirotrichea"  
 "#FF3333",  
# "Unclassified Ciliophora"
"#990000",  
# "Bacillariophyta" 
"#009999", 
# "Chlorophyceae"     
 "navyblue" , 
# "Chrysophyceae"  
"#00CCCC", 
# "Cryptophyceae" 
"#33CC33" ,
# "Dinophyceae"   
"#660009",
# "Diverse others" 
"gray",
#"Katablepharidaceae"  
"#FF33FF",
#"Spirotrichea"  
 "#FF3333"

)




dataset_bar %>% 
  dplyr::select("Division", "class") %>% 
  dplyr::rename("Top_Division" = "Division", 
         "Top_Class" = "class") %>% 
  mutate(group = paste(Top_Division, Top_Class, sep="-")) %>% 
  mutate(hex=hex) -> cdf_GP


custom_legend(mdf_GP, 
              cdf_GP, 
              group_level = "Division", 
              subgroup_level = "Class",
              legend_key_size = 1,
              legend_text_size = 12)

# To make two columns
# the package does not do that 
# so make two legends and combine them

mdf_GP %>%  dplyr::filter(Division %in% c("Ciliophora", "Conosa", "Ochrophyta", 
                                     "Cryptophyta", "Dinoflagellata", "Others" )) -> mdf_sub1

cdf_GP %>%  dplyr::filter(Top_Division %in% c("Ciliophora", "Conosa", "Ochrophyta", 
                                         "Cryptophyta", "Dinoflagellata", "Others" )) -> cdf_sub1




mdf_GP %>%  dplyr::filter(!Division %in% c("Ciliophora", "Conosa", "Ochrophyta", 
                                          "Cryptophyta", "Dinoflagellata", "Others" )) -> mdf_sub2

cdf_GP %>%  dplyr::filter(!Top_Division %in% c("Ciliophora", "Conosa", "Ochrophyta", 
                                              "Cryptophyta", "Dinoflagellata", "Others" )) -> cdf_sub2


custom_legend(mdf_sub1, 
              cdf_sub1, 
              group_level = "Division", 
              subgroup_level = "Class",
              legend_key_size = 1,
              legend_text_size = 14) -> leg_plot1

custom_legend(mdf_sub2, 
              cdf_sub2, 
              group_level = "Division", 
              subgroup_level = "Class",
              legend_key_size = 1,
              legend_text_size = 14) -> leg_plot2

comp_plot + leg_plot1 + leg_plot2 + plot_layout(widths =  c(1,1,1))
