library(dplyr)

# Load the data
taxa_table <- read.csv("/Users/theodoramautz/Desktop/SIO/Research/baja_edna/baja_edna_git/data/sequencing_data/Baja_MiFish_taxa_table_newseqs.csv")
asvs <- read.csv("/Users/theodoramautz/Desktop/SIO/Research/baja_edna/baja_edna_git/data/sequencing_data/Baja_MiFish_ASVs.csv")
hash_key <- read.csv("/Users/theodoramautz/Desktop/SIO/Research/baja_edna/baja_edna_git/data/sequencing_data/Baja_MiFish_hash_key.csv")

# Merge hash_key with taxa_table to get species information
merged_data <- merge(hash_key, taxa_table, by = "Sequence")

# Merge the resulting data with asvs to get nReads and Sample information
final_data <- merge(merged_data, asvs, by = "Hash")

# Calculate the number of samples and total reads for each species
species_rarity <- final_data %>%
  group_by(Species) %>%
  summarize(total_samples = n_distinct(Sample_name),
            total_reads = sum(nReads)) %>%
  arrange(total_samples, total_reads)

# Extract the three rarest species
rarest_species <- head(species_rarity, 10)

# Find the ASVs with the highest nReads for the three rarest species
result <- final_data %>%
  filter(Species %in% rarest_species$Species) %>%
  group_by(Species, Hash) %>%
  summarize(max_reads = max(nReads)) %>%
  arrange(Species, desc(max_reads))

print(result)
