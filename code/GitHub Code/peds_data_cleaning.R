# Some pre-analysis data cleaning steps

# Load the required libraries
library(dplyr)
library(ggplot2)
library(here)

# Read the csv files into a data frame
baja_asv_df <- read.csv(here("data", "sequencing_data", "rerun_bioinformatics_newseqs", "Baja_MiFish_ASVs_v4.csv"))
baja_taxonomy_df <- read.csv(here("data", "sequencing_data", "rerun_bioinformatics_newseqs", "Baja_MiFish_taxa_table_newseqs_v4.csv"))
baja_hash_df <- read.csv(here("data", "sequencing_data", "rerun_bioinformatics_newseqs", "Baja_MiFish_hash_key_v4.csv"))
baja_reads_df <- read.csv(here("data", "sequencing_data","rerun_bioinformatics_newseqs", "Baja_MiFish_reads_track_v4.csv"))

# Add BAJA_ to baja_asv_df
baja_asv_df$Sample_name <- paste0("BAJA_", baja_asv_df$Sample_name)
folder <- here("data", "sequencing_data")

# Add BAJA_ and remove leading 0s to baja_reads_df
baja_reads_df <- baja_reads_df %>%
  rename(Sample_name = X)
baja_reads_df$Sample_name <- sub("^0+", "", baja_reads_df$Sample_name)
baja_reads_df$Sample_name <- paste0("BAJA_", baja_reads_df$Sample_name)

# Add hashes to taxonomy df
baja_taxonomy_df <- left_join(baja_taxonomy_df, baja_hash_df, by = "Sequence")

# Combine the data frames based on the common column
combined_baja_df <- left_join(baja_asv_df, baja_taxonomy_df, by = "Hash")

# Find every unique Sample name
sample_name <- combined_baja_df$Sample_name
unique_samples <- unique(sample_name)
print(unique_samples)

# Read UD index names csv file
ud_index_df <- read.csv(here("data", "edna_data", "Bajalib_nextera_ud_indexes.xlsx - Baja library index.csv"))

# Append proper Baja sample ID column to combined_baja_df
Baja_ID_column <- ud_index_df[c("sample", "Seq_ID")]
Baja_ID_column$Seq_ID <- gsub("BAJA_", "", Baja_ID_column$Seq_ID) # remove "BAJA_"
Baja_ID_column$Seq_ID <- sub("^0+", "", Baja_ID_column$Seq_ID) # remove leading 0s
Baja_ID_column$Seq_ID <- paste0("BAJA_", "", Baja_ID_column$Seq_ID) # re-add "BAJA_"

# Delete rows after BAJA_135 (which are no longer Baja MiFish samples)
index <- which(Baja_ID_column$Seq_ID == "BAJA_136")
Baja_ID_column <- Baja_ID_column[1:(index-1), ]
appended_combined_baja_df <- merge(combined_baja_df, Baja_ID_column, by.x = "Sample_name", by.y = "Seq_ID", all.x = TRUE)

# Add locations, dive time, transport time, time of collection
baja_metadata_df <- read.csv(here("data", "edna_data", "eDNA Explore Baja REEF 2022 - All Data.csv"))
metadata_columns <- baja_metadata_df[c("sample", "Site", "Diver", "TimeofCollection_localtime", "BottomExposureTime_minutes", "MaxDepth", "TransportTime")]

# Remove BRPF- and leading 0's from those columns
metadata_columns$sample <- sub("^BRPF-", "", metadata_columns$sample)
metadata_columns$sample <- sub("^0+", "", metadata_columns$sample)

# Merge by sample
# Merge the two data frames based on the "sample" column
appended_locations_df <- merge(appended_combined_baja_df, metadata_columns, by = "sample", all.x = TRUE)

# Reorder based on sample
appended_locations_df$sample <- as.numeric(appended_locations_df$sample)
appended_locations_df <- appended_locations_df[order(appended_locations_df$sample), ]

# Rename sites to match REEF site names
secondcleaned_eDNA_df <- appended_locations_df %>%
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

# Specify folder to write files to
folder <- here("data")
path1 <- here::here(folder, "merged_baja.csv")

# Write as csv file
write.csv(finalcleaned_eDNA_df, file=path1, row.names=FALSE)

# Save the workspace to a .RData file
path2 <- here::here(folder, "baja_data_cleaning.RData")
save.image(file=path2)
