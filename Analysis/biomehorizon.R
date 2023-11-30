library(biomehorizon) 

# asv table  
# The first column contains microbial taxon IDs (or OTUs for 16S data), and
## all other columns are samples. Values represent sample reads per microbe within a given sample.
## Though in this case values are integer sample reads, they can also be represented as
## proportions or percentages of the total sample. Columns do not need defined names.

asv_table_rna_r %>% rownames_to_column("taxon_id")


# Metadata format
# Must include sample IDs that match the column names of otusample,
## subject IDs, and collection dates in either date or numeric format. 
## The columns with sample IDs, collection dates and subject names must be named 
## "sample", "collection_date" and "subject" respectively. 

env_lobau_final_rna_r %>% 
  dplyr::select("sample_id", "Date", "well_id.y") %>% 
  dplyr::rename("subject" = "well_id.y",
                "sample" = "sample_id") %>% as.data.frame() -> met_hor


# Taxonomydata format
# You can supply a vector of strings each with the entire taxonomy of a microbe,
## with levels separated by semicolons, or a table with columns for each taxonomic level where the first
## column is the OTU ID. Columns do not need defined names.
## Supports classification up to Subspecies (8 levels)

head(taxonomysample_diet)
tax_table %>% rownames_to_column("taxon_id")  %>% 
  mutate(div_class = paste(Division, Class, sep = "_")) -> tax_tab2






paramList <- prepanel(otudata = asv_table_rna_r %>% 
                        rownames_to_column("taxon_id"), 
                      metadata = met_hor,
                      taxonomydata = tax_tab2$div_class,
                      subj = "D15", 
                      facetLabelsByTaxonomy = TRUE,
                      thresh_prevalence = 40, thresh_abundance = 0.4,
                      maxGap = 35,
                      # order asvs
                      otulist = c(
                        # ciliophora
                        "ASV_233_pes", "ASV_3j4_6ar", "ASV_3o1_gt4", "ASV_4ky_nm7", 
                        "ASV_4m4_yrm",  "ASV_cs8_xfz", "ASV_92d_2wr", 
                        "ASV_e78_p1f", "ASV_edd_ug4", "ASV_g9j_ozs", "ASV_h7w_on6", 
                        "ASV_lky_7t9", "ASV_nmq_mou", "ASV_pcq_5ah", " ASV_q5h_df7", 
                        "ASV_qac_2d8", 
                        # dinoflagellata
                        "ASV_51o_7wt",
                        # chlorophyta
                        "ASV_536_3gb", "ASV_o3m_ak8", 
                        # ocrophyta,
                        "ASV_77l_gxf", "ASV_f0o_3b2", "ASV_fvs_jzl", "ASV_jh1_uc5",
                        "ASV_occ_rqd", 
                        #cryptophyta
                        "ASV_7i3_75n", "ASV_q3c_ijc",  "ASV_k26_zi8", 
                        "ASV_qfc_n4v", "ASV_bm6_pkn",
                        # lobosa
                        "ASV_t8z_hiq"
                        )
                      )

#c( 18472 18507 18541 18569 18605 18639 18660 18690 18751 18814 18907 18935 18970 19004 19032 19058 19110 19206 19235)

horizonplot(paramList) +
  ggplot2::scale_x_continuous(expand = c(0,0),
                              breaks =  seq(1,19,1),
                              labels = as.Date(env_lobau_final_rna_r[env_lobau_final_rna_r$well_id.y == "D15",]$Date) ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(5.5, 13.5), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(1.5, 4.5,
                            7.5, 9.5,
                            10.5, 12.5,
                            15.5, 17.5), linetype = "dashed", color = "gray") 
