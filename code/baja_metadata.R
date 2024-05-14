# eDNA Metadata file from raw data (no samples filtered out, unlike in
# "baja_analysis.R" file)

library(dplyr)
library(here)

# Load the workspace from the .RData file
load(here("data", "baja_data_cleaning.RData"))
load(here("data", "reef_data_cleaning.RData"))

# Add habitat type, lat, lon from REEF to eDNA file
# Find the most common Habitat for each Site_ID in merged_reef_survey
most_common_hab <- merged_reef_survey %>%
  group_by(Site_ID, Habitat) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  slice(1) %>%
  ungroup() %>%
  select(Site_ID, Most_Common_Habitat = Habitat)

# Remove blanks and empty rows from metadata_columns
metadata_columns <- metadata_columns[metadata_columns$Site != "" & metadata_columns$Diver != "blank", ]

# Rename sites to match REEF site names
metadata_columns <- metadata_columns %>%
  mutate(Site = case_when(
    Site == "Arecife del Pardito (near Isla Coyote)" ~ "Parrotfish City (Reef East of Isla Coyote)",
    Site == "Cabeza de Gorila (Isla Coronado)" ~ "Cabeza de Gorila / Gorilla Head (Isla Coronado)",
    Site == "Cayo Island" ~ "Cayo Island (near Isla Lobos) - Isla San Jose",
    Site == "Ensenada Grande / Mobula Dive Night Dive" ~ "Ensenada Grande Mobula Dive",
    Site == "Honeymoon (Isla Danzante) Night Dive" ~ "Honeymoon (Isla Danzante)",
    Site == "Lolo's Cove / Los Nidos" ~ "Lolo's Cove / Los Nidos Sud",
    Site == "Sea Lion Rock/La Cueva/The Cave (Isla Las Animas)" ~ "Sea Lion Rock / La Cueva / The Cave (Isla Las Animas)",
    Site == "Silhoutte (East side Isla Danzante)" ~ "Silhouette (East side Isla Danzante)",
    TRUE ~ Site  # Keep the existing value if none of the conditions match
  ))

# Add Site ID to Baja metadata
metadata_columns <- metadata_columns %>%
  mutate(Site_ID = paste0("Site ", dense_rank(Site)))

# Join most_common_hab with metadata_columns to include the most common Description
metadata_df <- left_join(metadata_columns, 
                                  most_common_hab, 
                                  by = "Site_ID")

# Join merged_reef_survey to include lat and lon
metadata_df <- left_join(metadata_df, merged_reef_survey %>% 
                                    select(Site_ID, lat, lon) %>% 
                                    distinct(), 
                                  by = "Site_ID")

# Convert lat and lon to decimal-degrees
# Function to separate degrees and minute-decimals
separate_degrees_minutes <- function(coord) {
  # Split the coordinate into degrees and minute-decimals
  split_coord <- strsplit(coord, " ")
  # Extract degrees and minute-decimals
  degrees <- as.numeric(sapply(split_coord, `[`, 1))
  minutes_decimal <- as.numeric(sapply(split_coord, `[`, 2))
  # Create a data frame with separated degrees and minute-decimals
  return(data.frame(deg = degrees, min = minutes_decimal))
}

# Apply the function to separate degrees and minute-decimals for latitude and longitude
metadata_df <- cbind(
  metadata_df,
  separate_degrees_minutes(metadata_df$lat),
  separate_degrees_minutes(metadata_df$lon)
)

# Rename the new columns
names(metadata_df)[(ncol(metadata_df)-3):ncol(metadata_df)] <- c("lat_deg", "lat_min", "lon_deg", "lon_min")

# Convert to decimal-degrees
metadata_df <- metadata_df %>%
  dplyr::mutate(lat.decimal = lat_deg + lat_min/60,
                lon.decimal = lon_deg - lon_min/60)

# Remove original lat / lon columns
metadata_df <- metadata_df[, -c(10:15)]

# Rename those columns to "lat" and "lon"
metadata_df <- metadata_df %>%
  rename(lat = lat.decimal) %>%
  rename(lon = lon.decimal)

# Assign latitude groups
metadata_df <- metadata_df %>%
  mutate(lat_group = case_when(
    lat >= 28 & lat <= 30 ~ "28-30",
    lat >= 25.7 & lat <= 27.9 ~ "25.7-27.9",
    lat >= 24 & lat <= 25.6 ~ "24-25.6",
    TRUE ~ NA_character_  # For values outside of specified ranges
  ))

# Add 'Seq_ID' column to metadata_df
metadata_df <- merge(metadata_df, Baja_ID_column, by = "sample")

# Change "sample" values from characters to numbers
metadata_df$sample <- as.numeric(metadata_df$sample)

# Order by Baja sample ID
metadata_df <- metadata_df[order(metadata_df$sample),]

# Move "Seq_ID" to first column and rename to "Sample_name"
metadata_df <- metadata_df[, c(13, 1:12)]
metadata_df <- metadata_df %>%
  rename('sample-id' = Seq_ID)

# Write metadata_df to a .tsv file with tab-separated values
write.table(metadata_df, file = "Baja_QIIME2_metadata.tsv", sep = "\t", row.names = FALSE)

