# First run the R scripts: "peds_data_cleaning.R" and "rvs_data_cleaning.R"

# Load the required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(here)

# Load the workspace from the .RData file
load(here("data", "baja_data_cleaning.RData"))
load(here("data", "reef_data_cleaning.RData"))

# Add lat, lon from REEF to eDNA file
finalcleaned_eDNA_df <- left_join(finalcleaned_eDNA_df, 
                                  merged_reef_survey %>% 
                                    select(Site_ID, lat, lon, lat_group) %>% 
                                    distinct(), 
                                  by = "Site_ID")

# Create metadata df without ASVs or taxa
edna_metadata_df <- finalcleaned_eDNA_df %>%
  select(-c(Label:Species))

# Condense into one row per sample
condensed_metadata <- edna_metadata_df %>%
  group_by(sample) %>%
  summarize_all(list(~ first(.)))

condensed_metadata <- condensed_metadata %>%
  select(Sample_name, sample, everything())

# Write as csv
folder <- here("data")
path1 <- here::here(folder, "baja_edna_metadata.csv")
write.csv(condensed_metadata, file=path1, row.names=FALSE)