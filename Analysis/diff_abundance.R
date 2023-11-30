




# differential abundance

# DeSeq2 ####
library(DESeq2)

# disputable since I manually added pseudocount 1 to each value in asv table in order to be able to compute geometric means...

# for this daa we need unrarified data

# remove those with few reads
asv_table[,colSums(asv_table) > 4000] -> asv_table_unrar
# remove them from env data
env_lobau_final[env_lobau_final$sample_id %in%  colnames(asv_table),] -> env_lobau_final_unrar



as.data.frame(env_lobau_final_unrar) %>% column_to_rownames("sample_id") -> env.data_unrar

ps_dna_unrar <- phyloseq(otu_table(
  ( as.matrix(asv_table_unrar) + 1) , 
  taxa_are_rows=TRUE), 
  tax_table(as.matrix(tax_table) ) , 
  sample_data(env.data_unrar ) )


# make a deseq object
dds.data = phyloseq_to_deseq2( ps_dna_unrar, ~ well_id.y)


dds = DESeq(dds.data)
dds$well_id.y
res <- results(dds)  

res <- res[order(res$padj, na.last=NA), ]
sigtab <- res %>%  as.data.frame()  %>% 
  dplyr::filter(padj < 0.05)  %>% 
  rownames_to_column("ASV") %>% 
  left_join(. , tax_table %>% rownames_to_column("ASV"), 
            by = "ASV")

# Bayesian estimation of dispersion
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)


ggplot(sigtab, aes(x=log2FoldChange, y=Class, color=Division)) + 
  geom_jitter(size=2, width = 0.2) +  
  scale_color_manual(values = distinct_palette(n = NA, pal = "brewerPlus")  ) 



# MA plots
# "Log fold-change vs. mean shows how well the homoskedasticity assumption holds,
# and identifies unusual OTUs where fold-change is at a limit."
plotMA(res)




### heatmap ####

select <- sigtab$ASV
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select, ] 

df <-  colData(dds) %>%  as.data.frame() %>%  dplyr::select(well_id.y)
df2 <- log2.norm.counts %>% as.data.frame() %>% 
  rownames_to_column("ASV") %>% 
  left_join(. , tax_table %>% rownames_to_column("ASV"), 
            by = "ASV") %>% 
  column_to_rownames("ASV") %>% 
  dplyr::select("Division") %>% 
  mutate(Division = as.factor(Division))

library(pheatmap)
pheatmap::pheatmap(log2.norm.counts %>% t(), 
                   annotation_col=df2,
                   annotation_row = df,
                   main="" ,  
                   show_colnames = F,
                   show_rownames = T) +
  scale_color_manual(values = microViz::distinct_palette(n = NA, pal = "brewerPlus")) + 
  scale_fill_manual(values = microViz::distinct_palette(n = NA, pal = "brewerPlus"))


# ANCOM ####
install.packages("ANCOM")
library(ANCOM)




# SIMPER ####
library(vegan)

# Cruaud says: 
# 25 OTUs that explained the most of the dissimilarity observed
# between the eukaryotic community structures of cold and warm
# season samples were selected for further analyses.


## step 1 ####
# just for groundwater and rna
asv_table_rna_r %>% 
  mutate(across(everything(), ~ ./sum(asv_table_rna_r))) %>% 
.[, 1:19] -> asv_simper

# disp
disp_s <- betadisper( 
  vegdist(t(asv_simper),
          method = "bray"), 
  env_lobau_final_rna_r[env_lobau_final_rna_r$well_id.y == "D15",]$season ) 

anova(disp_s) #no sign. differences in dispersion


# per,anova
permanova_s <- adonis2(t(asv_simper) ~ 
            env_lobau_final_rna_r[env_lobau_final_rna_r$well_id.y == "D15",]$season , 
                        method = "bray" , 
                        permutations = 999)
permanova_s # significant difference





## step 2 ####
# find the asvs that contribute the most to the seasonal differences

# asv_table in format taxa as columns, sample as rows 
data_matrix <- asv_simper %>% t() %>% as.data.frame()
group_my_data <- env_lobau_final_rna_r[env_lobau_final_rna_r$well_id.y == "D15",]$season %>% as.factor()

# Perform SIMPER analysis
simper_result <- vegan::simper(data_matrix, group = group_my_data, permutations = 999)
summary(simper_result)

# Get the 30 most pronounced contributors
rbind ( simper_result$summer_fall %>% 
          as.data.frame() %>% 
          dplyr::filter(p<0.05) %>% 
          dplyr::mutate(diff = 1), 
        
        # summer - winter
        simper_result$summer_winter %>% 
          as.data.frame() %>% 
          dplyr::filter(p<0.05) %>% 
          dplyr::mutate(diff = 2),
        
        # summer - spring
        simper_result$summer_spring %>% 
        as.data.frame() %>% 
          dplyr::filter(p<0.05) %>% 
          dplyr::mutate(diff = 3),
        
        # fall - winter
        simper_result$fall_winter %>% 
          as.data.frame() %>% 
          dplyr::filter(p<0.05) %>% 
          dplyr::mutate(diff = 4),
        
        # fall - spring
        simper_result$fall_spring %>% 
          as.data.frame() %>% 
          dplyr::filter(p<0.05) %>% 
          dplyr::mutate(diff = 5),
        
        # winter - spring
        simper_result$winter_spring %>% 
          as.data.frame() %>% 
          dplyr::filter(p<0.05) %>% 
          dplyr::mutate(diff = 6)
        ) -> top_contributors

top_contributors[!duplicated(top_contributors$species) | 
                   duplicated(top_contributors$species, fromLast = TRUE), ] %>% slice_max(order_by = average, n = 70) -> top_contributors

# 30 species
top_contributors %>% dim()



# Plot the abundances of the 30 most differential species 
asv_simper[ rownames(asv_simper) %in% top_contributors$species, ] %>% 
  rownames_to_column("asv") %>% 
  pivot_longer(cols = -asv) %>% 
  left_join(., env_lobau_final_rna_r[env_lobau_final_rna_r$well_id.y == "D15",],
            by = c("name" = "sample_id") ) %>% 
  dplyr::select("asv", "value", "Date") %>% 
  left_join(., tax_table %>% rownames_to_column("asv"), 
            by = "asv") %>% 
  ggplot(aes(x=Date, y=asv)) +
  geom_point(aes( alpha = value), size = 4) +
  theme_bw() + 
  geom_vline(xintercept = c(5.5, 13.5), linetype = "dashed", color = "gray") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45)) 
  scale_alpha_continuous(range = c(0, 0.004))
  #scale_colour_gradient(low = "green", high = "red") 

# not great..