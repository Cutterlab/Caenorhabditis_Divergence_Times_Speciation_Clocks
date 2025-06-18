library(tidyverse)
library(nlstools)

#------------------------------------------------------------------------------
#Generating NLS models using all available species-pair data
#------------------------------------------------------------------------------

data <- read.csv("./nls_and_terminalstage_data.csv")


#premating reproductive isolation-----------------------------------------------
premating_data <- data %>% #removing observations without premating RI values
  filter(!is.na(premating_RI))

#including adjusted conspecific coalescence times
prematingRI.nls <- nls(premating_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_adjusted))), 
                       data = premating_data, 
                       start = list(b = 0.005, a = 960))
summary(prematingRI.nls)

#conspecific coalescence times set to 0
prematingRI_conzero.nls <- nls(premating_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens))), 
                       data = premating_data, 
                       start = list(b = 0.005, a = 960))
summary(prematingRI_conzero.nls)


#excluding conspecific crosses
prematingRI_nocons.nls <- nls(premating_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_nocons))), 
                               data = premating_data, 
                               start = list(b = 0.005, a = 960))
summary(prematingRI_nocons.nls)


#prezygotic reproductive isolation---------------------------------------------

#including adjusted conspecific coaelscence times
prezygotic_data <- data %>% #removing observations without prezygotic RI values
  filter(!is.na(prezyg_RI))

prezygoticRI.nls <- nls(prezyg_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_adjusted))), 
                        data = prezygotic_data, 
                        start = list(b = 0.0028, a = 7))
summary(prezygoticRI.nls)


#conspecific coalescence times set to 0
prezygoticRI_conzero.nls <- nls(prezyg_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens))), 
                        data = prezygotic_data, 
                        start = list(b = 0.0028, a = 7))
summary(prezygoticRI_conzero.nls)


#excluding conspecific crosses
prezygoticRI_nocons.nls <- nls(prezyg_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_nocons))), 
                                data = prezygotic_data, 
                                start = list(b = 0.0028, a = 7))
summary(prezygoticRI_nocons.nls)




#postzygotic reproductive isolation---------------------------------------------
postzygotic_data <- data %>% #removing observations without postzygotic RI values
  filter(!is.na(postzyg_RI))


#including adjusted conspecifi coalescence times
postzygoticRI.nls <- nls(postzyg_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_adjusted))), 
                         data = postzygotic_data, 
                         start = list(b = 0.2, a = 3e04))
summary(postzygoticRI.nls)


#conspecific coalescence times set to 0
postzygoticRI_conzero.nls <- nls(postzyg_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens))), 
                         data = postzygotic_data, 
                         start = list(b = 0.2, a = 3e04))
summary(postzygoticRI_conzero.nls)


#excluding conspecific crosses
postzygoticRI_nocons.nls <- nls(postzyg_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_nocons))), 
                                 data = postzygotic_data, 
                                 start = list(b = 0.2, a = 3e04))
summary(postzygoticRI_nocons.nls)




#F1 reproductive isolation------------------------------------------------------
F1_data <- data %>% #removing observations without F1 RI values
  filter(!is.na(F1_RI))

#using adjusted conspecific coalescence times
F1RI.nls <- nls(F1_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_adjusted))), 
                data = F1_data, 
                start = list(b = 0.1, a = 595))

summary(F1RI.nls)


#conspecific coalescence times set to 0
F1RI_conzero.nls <- nls(F1_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens))), 
                data = F1_data, 
                start = list(b = 0.1, a = 595))

summary(F1RI_conzero.nls)


#excluding conspecific crosses
F1RI_nocons.nls <- nls(F1_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_nocons))), 
                        data = F1_data, 
                        start = list(b = 0.1, a = 595))

summary(F1RI_nocons.nls)



#F2 reproductive isolation------------------------------------------------------
F2_data <- data %>% #removing observations without F2 RI values
  filter(!is.na(F2_RI))

#using adjusted conspecific coalescence times
F2RI.nls <- nls(F2_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_adjusted))), 
                data = F2_data, 
                start = list(b = 0.5310, a = 3.825))
summary(F2RI.nls)


#conspecific coalescence times set to 0
F2RI_conzero.nls <- nls(F2_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens))), 
                data = F2_data, 
                start = list(b = 0.611, a = 4.15e11),
                nls.control(maxiter = 1000, minFactor = 0.00001)) #does not converge


#excluding conspecific coalescence times
F2RI_nocons.nls <- nls(F2_RI ~ 1 / (1 + a*exp(-b * (div_time_mgens_nocons))), 
                        data = F2_data, 
                        start = list(b = 0.611, a = 4.15e11),
                        nls.control(maxiter = 1000, minFactor = 0.00001)) #does not converge




#-------------------------------------------------------------------------------
#Selecting species pairs to use in NLS model of phylogenetic independent contrasts
#-------------------------------------------------------------------------------


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


#premating reproductive isolation--------------------------------------------------------------------

premating_data_PIC <- PIC_data %>% #removing observations without premating RI values
  filter(!is.na(premating_RI))

visited <- list()
good_entries <- list()
j <- 1
for (i in 1:nrow(premating_data_PIC)) {
  if (!premating_data_PIC$species1[i] %in% visited & !premating_data_PIC$species2[i] %in% visited) {
    visited <- append(visited, c(premating_data_PIC$species1[i], premating_data_PIC$species2[i]))
    good_entries[[j]] <- list(
      premating_data_PIC$species1[i], 
      premating_data_PIC$species2[i],
      premating_data_PIC$div_time_mgens_adjusted[i]
    )
    j <- j + 1
  }
}


#add the returned combinations of species to dataframe
premating_data_PIC <- as.data.frame(do.call("rbind", good_entries)) %>% 
  rename(species1 = V1, 
         species2 = V2, 
         divergence_time_mgens_adjusted = V3) %>% 
  
  #change class of div time data from list to numeric
  mutate(divergence_time_mgens_adjusted = as.numeric(divergence_time_mgens_adjusted)) %>%
  
  #create columns for each reciprocal cross
  unite(species_combo, species1, species2, sep = "-", remove = FALSE) %>%
  unite(species_combo_recip, species2, species1, sep = "-", remove = TRUE) %>%
  
  #add in prezygotic RI data from original dataset
  left_join(., PIC_data, by = "species_combo") %>%
  rename(premating_RI_cross1 = premating_RI) %>%
  left_join(., data, join_by("species_combo_recip" == "species_combo")) %>%
  rename(premating_RI_cross2 = premating_RI) %>%
  
  #move improtant columns to the front and remove unecessary columns from left joins
  relocate(c(premating_RI_cross1, premating_RI_cross2), .after = divergence_time_mgens_adjusted) %>%
  select(!c(6:29,)) %>%
  
  #fill NA values for prezygotic RI with data from reciprocal cross
  mutate(premating_RI_cross2 = coalesce(premating_RI_cross2, premating_RI_cross1))

#add column that averages the RI between reciprocal crosses
premating_data_PIC <- premating_data_PIC %>%
  mutate(mean_premating_RI = rowMeans(select(premating_data_PIC, c(premating_RI_cross1, premating_RI_cross2))))
  

prematingRI_PIC.nls <- nls(mean_premating_RI ~ 1 / (1 + a*exp(-b * (divergence_time_mgens_adjusted))), 
                       data = premating_data_PIC, 
                       start = list(b = 0.005, a = 960),
                       nls.control(maxiter = 1000))
#nls fails to converge 


#prezygotic reproductive isolation------------------------------------------------------------------

#removing observations without prezygotic RI values
prezygotic_data_PIC <- PIC_data %>% 
  filter(!is.na(prezyg_RI))


#algorithm for selecting the maximum number of pairs of species that are phylogenetically independent
visited <- list()
good_entries <- list()
j <- 1
for (i in 1:nrow(prezygotic_data_PIC)) {
  if (!prezygotic_data_PIC$species1[i] %in% visited & !prezygotic_data_PIC$species2[i] %in% visited) {
    visited <- append(visited, c(prezygotic_data_PIC$species1[i], prezygotic_data_PIC$species2[i]))
    good_entries[[j]] <- list(
      prezygotic_data_PIC$species1[i], 
      prezygotic_data_PIC$species2[i],
      prezygotic_data_PIC$div_time_mgens_adjusted[i]
    )
    j <- j + 1
  }
}

#add the returned combinations of species to dataframe
prezygotic_data_PIC <- as.data.frame(do.call("rbind", good_entries)) %>% 
  rename(species1 = V1, 
         species2 = V2, 
         divergence_time_mgens_adjusted = V3) %>% 
  
  #change class of div time data from list to numeric
  mutate(divergence_time_mgens_adjusted = as.numeric(divergence_time_mgens_adjusted)) %>%
  
  #create columns for each reciprocal cross
  unite(species_combo, species1, species2, sep = "-", remove = FALSE) %>%
  unite(species_combo_recip, species2, species1, sep = "-", remove = TRUE) %>%
  
  #add in prezygotic RI data from original dataset
  left_join(., PIC_data, by = "species_combo") %>%
  rename(prezygotic_RI_cross1 = prezyg_RI) %>%
  left_join(., data, join_by("species_combo_recip" == "species_combo")) %>%
  rename(prezygotic_RI_cross2 = prezyg_RI) %>%
  
  #move important columns to the front and remove unecessary columns from left joins
  relocate(c(prezygotic_RI_cross1, prezygotic_RI_cross2), .after = divergence_time_mgens_adjusted) %>%
  select(!c(6:29,)) %>%
  
  #fill NA values for prezygotic RI with data from reciprocal cross
  mutate(prezygotic_RI_cross2 = coalesce(prezygotic_RI_cross2, prezygotic_RI_cross1))
  
#add column that averages the RI between reciprocal crosses
prezygotic_data_PIC <- prezygotic_data_PIC %>%
  mutate(mean_prezygotic_RI = rowMeans(select(prezygotic_data_PIC, c(prezygotic_RI_cross1, prezygotic_RI_cross2))))

#nls
prezygoticRI_PIC.nls <- nls(mean_prezygotic_RI ~ 1 / (1 + a*exp(-b * (divergence_time_mgens_adjusted))), 
                        data = prezygotic_data_PIC, 
                        start = list(b = 0.0028, a = 7))
summary(prezygoticRI_PIC.nls)



#postzygotic reproductive isolation----------------------------------------------------------------------

#removing observations without postzygotic RI values
postzygotic_data_PIC <- PIC_data %>% 
  filter(!is.na(postzyg_RI))


#algorithm for selecting the maximum number of pairs of species that are phylogenetically independent
visited <- list()
good_entries <- list()
j <- 1
for (i in 1:nrow(postzygotic_data_PIC)) {
  if (!postzygotic_data_PIC$species1[i] %in% visited & !postzygotic_data_PIC$species2[i] %in% visited) {
    visited <- append(visited, c(postzygotic_data_PIC$species1[i], postzygotic_data_PIC$species2[i]))
    good_entries[[j]] <- list(
      postzygotic_data_PIC$species1[i], 
      postzygotic_data_PIC$species2[i],
      postzygotic_data_PIC$div_time_mgens_adjusted[i]
    )
    j <- j + 1
  }
}

#add the returned combinations of species to dataframe
postzygotic_data_PIC <- as.data.frame(do.call("rbind", good_entries)) %>% 
  rename(species1 = V1, 
         species2 = V2, 
         divergence_time_mgens_adjusted = V3) %>% 
  
  #change class of div time data from list to numeric
  mutate(divergence_time_mgens_adjusted = as.numeric(divergence_time_mgens_adjusted)) %>%
  
  #create columns for each reciprocal cross
  unite(species_combo, species1, species2, sep = "-", remove = FALSE) %>%
  unite(species_combo_recip, species2, species1, sep = "-", remove = TRUE) %>%
  
  #add in postzygotic RI data from original dataset
  left_join(., PIC_data, by = "species_combo") %>%
  rename(postzygotic_RI_cross1 = postzyg_RI) %>%
  left_join(., data, join_by("species_combo_recip" == "species_combo")) %>%
  rename(postzygotic_RI_cross2 = postzyg_RI) %>%
  
  #move important columns to the front and remove unecessary columns from left joins
  relocate(c(postzygotic_RI_cross1, postzygotic_RI_cross2), .after = divergence_time_mgens_adjusted) %>%
  select(!c(6:29,)) %>%
  
  #fill NA values for postzygotic RI with data from reciprocal cross
  mutate(postzygotic_RI_cross2 = coalesce(postzygotic_RI_cross2, postzygotic_RI_cross1))

#add column that averages the RI between reciprocal crosses
postzygotic_data_PIC <- postzygotic_data_PIC %>%
  mutate(mean_postzygotic_RI = rowMeans(select(postzygotic_data_PIC, c(postzygotic_RI_cross1, postzygotic_RI_cross2))))


#nls
postzygoticRI_PIC.nls <- nls(mean_postzygotic_RI ~ 1 / (1 + a*exp(-b * (divergence_time_mgens_adjusted))), 
                            data = postzygotic_data_PIC, 
                            start = list(b = 0.2, a = 3e04))
summary(postzygoticRI_PIC.nls)




#F1 reproductive isolation----------------------------------------------------------------------------

#removing observations without F1 RI values
F1_data_PIC <- PIC_data %>% 
  filter(!is.na(F1_RI))


#algorithm for selecting the maximum number of pairs of species that are phylogenetically independent
visited <- list()
good_entries <- list()
j <- 1
for (i in 1:nrow(F1_data_PIC)) {
  if (!F1_data_PIC$species1[i] %in% visited & !F1_data_PIC$species2[i] %in% visited) {
    visited <- append(visited, c(F1_data_PIC$species1[i], F1_data_PIC$species2[i]))
    good_entries[[j]] <- list(
      F1_data_PIC$species1[i], 
      F1_data_PIC$species2[i],
      F1_data_PIC$div_time_mgens_adjusted[i]
    )
    j <- j + 1
  }
}

#add the returned combinations of species to dataframe
F1_data_PIC <- as.data.frame(do.call("rbind", good_entries)) %>% 
  rename(species1 = V1, 
         species2 = V2, 
         divergence_time_mgens_adjusted = V3) %>% 
  
  #change class of div time data from list to numeric
  mutate(divergence_time_mgens_adjusted = as.numeric(divergence_time_mgens_adjusted)) %>%
  
  #create columns for each reciprocal cross
  unite(species_combo, species1, species2, sep = "-", remove = FALSE) %>%
  unite(species_combo_recip, species2, species1, sep = "-", remove = TRUE) %>%
  
  #add in F1 RI data from original dataset
  left_join(., PIC_data, by = "species_combo") %>%
  rename(F1_RI_cross1 = F1_RI) %>%
  left_join(., data, join_by("species_combo_recip" == "species_combo")) %>%
  rename(F1_RI_cross2 = F1_RI) %>%
  
  #move improtant columns to the front and remove unecessary columns from left joins
  relocate(c(F1_RI_cross1, F1_RI_cross2), .after = divergence_time_mgens_adjusted) %>%
  select(!c(6:29,)) %>%
  
  #fill NA values for F1 RI with data from reciprocal cross
  mutate(F1_RI_cross2 = coalesce(F1_RI_cross2, F1_RI_cross1))

#add column that averages the RI between reciprocal crosses
F1_data_PIC <- F1_data_PIC %>%
  mutate(mean_F1_RI = rowMeans(select(F1_data_PIC, c(F1_RI_cross1, F1_RI_cross2))))

#nls
F1_PIC.nls <- nls(mean_F1_RI ~ 1 / (1 + a*exp(-b * (divergence_time_mgens_adjusted))), 
                             data = F1_data_PIC, 
                             start = list(b = 0.1, a = 595))
summary(F1_PIC.nls)



#F2 reproductive isolation------------------------------------------------------------------------------

#removing observations without F2 RI values
F2_data_PIC <- PIC_data %>% 
  filter(!is.na(F2_RI))


#algorithm for selecting the maximum number of pairs of species that are phylogenetically independent
visited <- list()
good_entries <- list()
j <- 1
for (i in 1:nrow(F2_data_PIC)) {
  if (!F2_data_PIC$species1[i] %in% visited & !F2_data_PIC$species2[i] %in% visited) {
    visited <- append(visited, c(F2_data_PIC$species1[i], F2_data_PIC$species2[i]))
    good_entries[[j]] <- list(
      F2_data_PIC$species1[i], 
      F2_data_PIC$species2[i],
      F2_data_PIC$div_time_mgens_adjusted[i]
    )
    j <- j + 1
  }
}

#add the returned combinations of species to dataframe
F2_data_PIC <- as.data.frame(do.call("rbind", good_entries)) %>% 
  rename(species1 = V1, 
         species2 = V2, 
         divergence_time_mgens_adjusted = V3) %>% 
  
  #change class of div time data from list to numeric
  mutate(divergence_time_mgens_adjusted = as.numeric(divergence_time_mgens_adjusted)) %>%
  
  #create columns for each reciprocal cross
  unite(species_combo, species1, species2, sep = "-", remove = FALSE) %>%
  unite(species_combo_recip, species2, species1, sep = "-", remove = TRUE) %>%
  
  #add in F2 RI data from original dataset
  left_join(., PIC_data, by = "species_combo") %>%
  rename(F2_RI_cross1 = F2_RI) %>%
  left_join(., data, join_by("species_combo_recip" == "species_combo")) %>%
  rename(F2_RI_cross2 = F2_RI) %>%
  
  #move improtant columns to the front and remove unecessary columns from left joins
  relocate(c(F2_RI_cross1, F2_RI_cross2), .after = divergence_time_mgens_adjusted) %>%
  select(!c(6:29,)) %>%
  
  #fill NA values for F1 RI with data from reciprocal cross
  mutate(F2_RI_cross2 = coalesce(F2_RI_cross2, F2_RI_cross1))

#add column that averages the RI between reciprocal crosses
F2_data_PIC <- F2_data_PIC %>%
  mutate(mean_F2_RI = rowMeans(select(F2_data_PIC, c(F2_RI_cross1, F2_RI_cross2))))

#nls
F2_PIC.nls <- nls(mean_F2_RI ~ 1 / (1 + a*exp(-b * (divergence_time_mgens_adjusted))), 
                  data = F2_data_PIC, 
                  start = list(b = 0.5310, a = 3.825))
summary(F2_PIC.nls)

#nls gradient is singular













