# Load the required libraries
library(dplyr)
library(ggplot2)
library(here)

# Read REEF csv files
reefsurvey_species <- read.csv(here("data", "reef_data", "TEPspecies-5.csv"))
reefsurvey_geo <- read.csv(here("data", "reef_data", "Bajageog.csv"))
reefsurvey <- read.csv(here("data", "reef_data", "ExploreBaja22.csv"))

# Combine data frames by species
merged1_reef_survey <- left_join(reefsurvey, reefsurvey_species, by = "species")

# Add z's to the one location in REEF Surveys not in eDNA data so that it will get assigned the last site number in analysis
# Rename sites to match REEF site names
merged2_reef_survey <- merged1_reef_survey %>%
  mutate(Location = case_when(
    Location == "Chayo's Cave (San Pedro Martir)" ~ "zz Chayo's Cave (San Pedro Martir)",
    TRUE ~ Location  # Keep the existing value if none of the conditions match
  ))

# Add "Site 1" etc to each site by alphabetical order
merged3_reef_survey <- merged2_reef_survey %>%
  mutate(Site_ID = paste0("Site ", dense_rank(Location)))

# Rename species
merged_reef_survey <- merged3_reef_survey %>%
  mutate(scientificname = case_when(
    scientificname == "Paranthias colonus" ~ "Cephalopholis colonus",
    scientificname == "Hermosilla azurea" ~ "Kyphosus azureus",
    TRUE ~ scientificname # Keep the existing name if none of the conditions match
  ))

# Add habitat description to each site
reef_habitat_types <- read.csv(here("data", "reef_data", "REEF habitat types.csv"))
merged_reef_survey$Description <- reef_habitat_types$Description[match(merged3_reef_survey$Habitat, reef_habitat_types$Habitat)]

# Write to .csv
folder <- here("data")
path3 <- here::here(folder, "merged_reef.csv")
write.csv(merged_reef_survey, file=path3, row.names=FALSE)

# Save the workspace to a .RData file
path4 <- here::here(folder, "reef_data_cleaning.RData")
save.image(file=path4)

