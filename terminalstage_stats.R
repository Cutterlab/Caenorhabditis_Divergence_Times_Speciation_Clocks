

data <- read.csv("./nls_and_terminalstage_data.csv")


#Non-phylogenetically independent data----------------------------------------------------------

#using a Kruskal-Wallis test to look for differences in div time based on latest terminal stage
kruskal.test(div_time_mgens_adjusted ~ terminal_stage, data = data) #p = 1.194e-15


#pairwise post-hoc tests
pairwise.wilcox.test(data$div_time_mgens_adjusted, data$terminal_stage, p.adj = "bonf", paired = FALSE, exact = FALSE)


#computing spearman's rank correlation
cor.test(x = data$F1_stage, y = data$div_time_mgens_adjusted, method = "spearman") 
#rho=0.427 #p<1e-11



#Phylogenetically independent data-------------------------------------------------------------


PIC_data <- data %>%
  #remove unnecessary columns
  select(!c(14:16),) %>%
  #arrange adjusted divergence times from smallest to largest (mgens)
  arrange(div_time_mgens_adjusted) %>%
  #remove NA values in divergence time data
  filter(!is.na(div_time_mgens_adjusted)) %>%
  #create new columns that identify each species in a pair individually
  separate(species_combo, into = c("species1", "species2"), sep = "-", remove = FALSE) %>%
  #filter out conspecific pairs from new columns 
  filter(species1 != species2)

#removing observations without terminal life stage values
vs_data_PIC <- PIC_data %>% 
  filter(!is.na(terminal_stage))


#algorithm for selecting the maximum number of pairs of species that are phylogenetically independent
visited <- list()
good_entries <- list()
j <- 1
for (i in 1:nrow(vs_data_PIC)) {
  if (!vs_data_PIC$species1[i] %in% visited & !vs_data_PIC$species2[i] %in% visited) {
    visited <- append(visited, c(vs_data_PIC$species1[i], vs_data_PIC$species2[i]))
    good_entries[[j]] <- list(
      vs_data_PIC$species1[i], 
      vs_data_PIC$species2[i],
      vs_data_PIC$div_time_mgens_adjusted[i]
    )
    j <- j + 1
  }
}

#add the returned combinations of species to dataframe
vs_data_PIC <- as.data.frame(do.call("rbind", good_entries)) %>% 
  rename(species1 = V1, 
         species2 = V2, 
         divergence_time_mgens_adjusted = V3) %>% 
  
  #change class of div time data from list to numeric
  mutate(divergence_time_mgens_adjusted = as.numeric(divergence_time_mgens_adjusted)) %>%
  
  #create columns for each reciprocal cross
  unite(species_combo, species1, species2, sep = "-", remove = FALSE) %>%
  unite(species_combo_recip, species2, species1, sep = "-", remove = TRUE) %>%
  
  #add in vs RI data from original dataset
  left_join(., PIC_data, by = "species_combo") %>%
  rename(vs_cross1 = terminal_stage) %>%
  left_join(., data, join_by("species_combo_recip" == "species_combo")) %>%
  rename(vs_cross2 = terminal_stage) %>%
  
  #move improtant columns to the front and remove unecessary columns from left joins
  relocate(c(vs_cross1, vs_cross2), .after = divergence_time_mgens_adjusted) %>%
  select(!c(6:29,)) %>%
  
  #fill NA values for F1 RI with data from reciprocal cross
  mutate(vs_cross2 = coalesce(vs_cross2, vs_cross1))  


#Fill in missing values and gather all the latest stages out of both reciprocal crosses into one column (using the vs_cross2 column) and rename it to "latest_stage".
vs_data_PIC[8,10] <- 3
vs_data_PIC$latest_stage <- vs_data_PIC$vs_cross2
vs_data_PIC[5,9] <- "sterile_adult"
vs_data_PIC[8,9] <- "larvae"



#reorder terminal stage factor levels
vs_data_PIC$latest_stage <- factor(vs_data_PIC$latest_stage, levels = c("adult", 
                                                                        "sterile_adult", 
                                                                        "larvae", 
                                                                        "emb", 
                                                                        "ambig_emb", 
                                                                        "no fert"))


#using a Kruskal-Wallis test to look for differences in div time based on latest terminal stage
kruskal.test(divergence_time_mgens_adjusted ~ latest_stage, data = vs_data_PIC) #p=0.06884

#computing spearman's rank correlation
cor.test(x = vs_data_PIC$F1_stage, y = vs_data_PIC$divergence_time_mgens_adjusted, method = "spearman") 
#p=0.008267, rho=0.653




