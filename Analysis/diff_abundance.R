




# differential abundance ####

 # Highly disputable since I needed to add 1 to each value in asv table in order to be able to compute geometric means...

# for daa we need unrarified data

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




# heatmap ####

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


pheatmap::pheatmap(log2.norm.counts %>% t(), 
                   annotation_col=df2,
                   annotation_row = df,
                   main="" ,  
                   show_colnames = F,
                   show_rownames = T) +
  scale_color_manual(values = distinct_palette(n = NA, pal = "brewerPlus")) + 
  scale_fill_manual(values = distinct_palette(n = NA, pal = "brewerPlus"))


