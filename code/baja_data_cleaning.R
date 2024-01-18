# Load the required libraries
library(dplyr)
library(ggplot2)
library(here)

# Read the csv files into a data frame
baja_asv_df <- read.csv(here("edna_data", "BAJALIB_MFU_ASV_table.csv"))
baja_taxonomy_df <- read.csv(here("edna_data", "BAJALIB_MFU_taxonomy_output.csv"))

# Combine the data frames based on the common column
combined_baja_df <- left_join(baja_asv_df, baja_taxonomy_df, by = "Hash")

# Sort combined_baja_df alphabetically by Sample_name
combined_baja_sorted_df <- combined_baja_df[order(combined_baja_df$Sample_name),]

# Find every unique Sample name
sample_name <- combined_baja_sorted_df$Sample_name
unique_samples <- unique(sample_name)
print(unique_samples)

# Read UD index names csv file
ud_index_df <- read.csv(here("edna_data", "Bajalib_nextera_ud_indexes.xlsx - Baja library index.csv"))

# Remove extra parts of name from Sample_name column in combined_baja_sorted_df
Sample_name <- combined_baja_sorted_df$Sample_name
combined_baja_sorted_df$Sample_name <- gsub("MFU-|-d1-1$", "", combined_baja_sorted_df$Sample_name)

# Append proper Baja sample ID column to combined_baja_sorted_df
Baja_ID_column <- ud_index_df[c("sample", "Seq_ID")]
appended_combined_baja_sorted_df <- merge(combined_baja_sorted_df, Baja_ID_column, by.x = "Sample_name", by.y = "Seq_ID", all.x = TRUE)

# Change "sample" values from characters to numbers
appended_combined_baja_sorted_df$sample <- as.numeric(appended_combined_baja_sorted_df$sample)

# Order by Baja sample ID
Baja_ID_sorted_df <- appended_combined_baja_sorted_df[order(appended_combined_baja_sorted_df$sample),]

# Add locations to Baja sample ID spreadsheet
# First manually rename first column of 'eDNA Explore Baja REEF 2022 - All Data.csv' to "sample"
# Also manually find and replace "El Lavadero (Las Animas) " with a space at the end to get rid of space at the end
baja_metadata_df <- read.csv(here("edna_data", "eDNA Explore Baja REEF 2022 - All Data.csv"))
location_column <- baja_metadata_df[c("sample", "Site")]
# Remove BRPF- and leading 0's from those columns
location_column$sample <- sub("BRPF-", "", location_column$sample)
location_column$sample <- sub("^0+", "", location_column$sample)
appended_locations_df <- merge(Baja_ID_sorted_df, location_column, all.x = TRUE)

# Remove terrestrial species
terrestrial_species <- c("Homo sapiens", "Sus scrofa", "Bos taurus", "Felis catus")
firstcleaned_eDNA_df <- appended_locations_df %>% filter(!Species %in% terrestrial_species)

# Rename sites to match REEF site names
secondcleaned_eDNA_df <- firstcleaned_eDNA_df %>%
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

# Add "Site 1" etc to each site by alphabetical order
finalcleaned_eDNA_df <- secondcleaned_eDNA_df %>%
  mutate(Site_ID = paste0("Site ", dense_rank(Site)))

# Write as csv file
write.csv(finalcleaned_eDNA_df, "merged_baja.csv", row.names=FALSE)

# Save the workspace to a .RData file
save.image("baja_data_cleaning.RData")
