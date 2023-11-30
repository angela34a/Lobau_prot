# continuation of the complete_analysis script and the nmds analysis
 
disp_season <- betadisper( 
 sp.dist, 
  env_lobau_final_r$season ) 

anova(disp_season)

disp_season$distances %>% as.data.frame() %>% 
  dplyr::rename("distance" = ".") %>% 
  cbind(.,  disp_season$group  %>%  as.data.frame()  ) %>% 
  dplyr::rename("Season" = ".") %>% 
  ggplot( aes(x= Season, y = distance)) +
  geom_boxplot() +
  annotate(geom="text", x=1.1, y=0.35, label="p = 0.0526",
           color="black") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_x_discrete(limits =c("fall", "winter", "spring", "summer") ) +
  labs(y = "Distance to centroid")
  
