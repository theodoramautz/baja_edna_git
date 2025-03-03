# Some pre-analysis data cleaning steps

# Load the required libraries
library(dplyr)
library(ggplot2)
library(here)

# Read REEF csv files
reefsurvey_species <- read.csv(here("data", "reef_data", "TEPspecies-5.csv"))
reefsurvey_geo <- read.csv(here("data", "reef_data", "Bajageog.csv"))
reefsurvey <- read.csv(here("data", "reef_data", "ExploreBaja22.csv"))
sites_lat_lon <- read.csv(here("data","reef_data", "Bajageog Google maps.csv"))

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

# Add lat and lon to sites
merged_reef_survey <- merged_reef_survey %>%
  left_join(sites_lat_lon, by = "Site_ID")

# Convert lat and lon to decimal-degrees
separate_degrees_minutes <- function(coord) {
  split_coord <- strsplit(coord, " ")
  degrees <- as.numeric(sapply(split_coord, `[`, 1))
  minutes_decimal <- as.numeric(sapply(split_coord, `[`, 2))
  return(data.frame(deg = degrees, min = minutes_decimal))
}

merged_reef_survey <- cbind(
  merged_reef_survey,
  separate_degrees_minutes(merged_reef_survey$lat),
  separate_degrees_minutes(merged_reef_survey$lon)
)

# Rename the new columns
names(merged_reef_survey)[(ncol(merged_reef_survey)-3):
                            ncol(merged_reef_survey)] <- c("lat_deg", 
                                                           "lat_min", 
                                                           "lon_deg", 
                                                           "lon_min")

# Convert to decimal-degrees
merged_reef_survey <- merged_reef_survey %>%
  dplyr::mutate(lat.decimal = lat_deg + lat_min/60,
                lon.decimal = lon_deg - lon_min/60)

# Remove original lat / lon columns
merged_reef_survey <- merged_reef_survey[, -c(29:34)]

# Rename those columns to "lat" and "lon"
merged_reef_survey <- merged_reef_survey %>%
  rename(lat = lat.decimal) %>%
  rename(lon = lon.decimal)

# Assign latitude groups
merged_reef_survey <- merged_reef_survey %>%
  mutate(lat_group = case_when(
    lat >= 28 & lat <= 30 ~ "28-30",
    lat >= 25.7 & lat <= 27.9 ~ "25.7-27.9",
    lat >= 24 & lat <= 25.6 ~ "24-25.6",
    TRUE ~ NA_character_  # For values outside of specified ranges
  ))

# Add Family name
family <- read.table(here("data", "reef_data", "TEPfamily0823.txt"), 
                     header = TRUE, sep = "\t")
colnames(family)[1] <- "Family"
colnames(family)[2] <- "family_scientific"
family <- family[, 1:2]
merged_reef_survey <- merge(merged_reef_survey, family, by = "Family", all.x = TRUE)
merged_reef_survey <- merged_reef_survey[, c(2:length(merged_reef_survey), 1)]

# Write to .csv
folder <- here("data")
path3 <- here::here(folder, "merged_reef.csv")
write.csv(merged_reef_survey, file=path3, row.names=FALSE)

# Save the workspace to a .RData file
path4 <- here::here(folder, "reef_data_cleaning.RData")
save.image(file=path4)