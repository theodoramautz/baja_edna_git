############################################################################################
# PRELIMINARY
############################################################################################
# Load the required libraries
library(dplyr)
library(ggplot2)
library(ggvenn)
library(grid)
library(gridExtra)
library(lattice)
library(vegan)
library(tidyr)
library(tidyverse)
library(here)

# Load the workspace from the .RData file
load(here("data", "baja_data_cleaning.RData"))
load(here("data", "reef_data_cleaning.RData"))

############################################################################################
# EDNA SPECIES
############################################################################################
# List of unique eDNA species
eDNA_species <- finalcleaned_eDNA_df$Species
eDNA_species_df <- as.data.frame(eDNA_species)
eDNA_species_number <- na.omit(eDNA_species_df)
eDNA_species_list <- as.data.frame(unique(eDNA_species_number))

# Make a bar chart with ASV counts per sample
sample_asv_counts <- table(finalcleaned_eDNA_df$sample)
sample_asv_counts_df <- as.data.frame(sample_asv_counts)
o <- ggplot(sample_asv_counts_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  labs(y = "Number of ASVs", x = "Sample") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Number of ASVs Per eDNA Sample")
print(o)

# Count how many unique samples had ASVs at all
total_unique_samples <- finalcleaned_eDNA_df %>%
  distinct(sample) %>%
  count()
print(total_unique_samples)

# Print bar graph of sighted species in eDNA
t <- ggplot(eDNA_species_number, aes(x=reorder(eDNA_species, -table(eDNA_species)[eDNA_species]))) +
  geom_bar() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  labs(y = "Number of Detections", x = "Species") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("eDNA Detection Counts of Species")+
  coord_flip()
print(t)

############################################################################################
# REEF SPECIES
############################################################################################
# Make a frequency table for abundance of species sightings in REEF surveys 
frequency_table <- table(merged_reef_survey$scientificname)
species_frequency_df <- as.data.frame(frequency_table)
species_frequency_df <- species_frequency_df[order(species_frequency_df$Freq), ]

# List of unique REEF species
REEF_species <- species_frequency_df$Var1
REEF_species_df <- as.data.frame(REEF_species)

# Bar chart of REEF survey observations per all species
r <- ggplot(data=merged_reef_survey, aes(x=reorder(scientificname, -table(scientificname)[scientificname]))) +
    geom_bar() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + 
    labs(y = "Number of observations", x = "Species") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Taxa Counts Per Site, REEF Survey") +
    coord_flip()
print(r)
  
# Filter species to those observed >200 times
species_frequency_filtered <- subset(species_frequency_df, Freq > 200)

# Plot only species observed >200 times in REEF surveys
species_frequency_filtered$Var1 <- factor(species_frequency_filtered$Var1,
  levels = species_frequency_filtered$Var1[order(species_frequency_filtered$Freq, decreasing = TRUE)])

s <- ggplot(species_frequency_filtered, aes(x=Freq, y=reorder(Var1, -Freq))) +
  geom_bar(stat = "identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
    labs(y = "Species", x = "Frequency of Observation") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("REEF Surveys Observation Counts for Frequently Observed Species")
print(s)

############################################################################################
# SPECIES OVERLAP
############################################################################################
# See how many species sighted in REEF surveys were detected in eDNA samples
common_species <- intersect(eDNA_species, REEF_species)
common_species_df <- as.data.frame(common_species)
length(common_species)

# See which species are in eDNA that were not in REEF surveys
in_eDNA_notin_REEF_df <- anti_join(eDNA_species_list, common_species_df, by = c("eDNA_species"="common_species"))

# See how many of those common species are in the most abundant REEF species list
abundant_REEF_species <- species_frequency_filtered$Var1
abundant_REEF_species_df <- as.data.frame(abundant_REEF_species)
common_abundant_species <- intersect(eDNA_species, abundant_REEF_species)
common_abundant_species_df <- as.data.frame(common_abundant_species)
length(common_abundant_species)

# Make Venn Diagram of species of overlap
eDNA_list <- unique(eDNA_species_number$eDNA_species)
REEF_list <- as.character(unique(REEF_species_df$REEF_species))
both_species <- list(eDNA = eDNA_list, REEF = REEF_list)
venn_plot <- ggvenn(both_species, columns = c("eDNA", "REEF"), 
                    stroke_size = 1, 
                    set_name_color = "black",
                    set_name_size = 8,
                    text_color = "black",
                    text_size = 5, fill_color = c("seagreen3", "plum"))
print(venn_plot)
############################################################################################
# SITES
############################################################################################
# Make a bar chart with site taxa counts
# eDNA samples
p <- ggplot(data=finalcleaned_eDNA_df, aes(x = reorder(Site, -table(Site)[Site]))) +
  geom_bar() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + 
  labs(y = "Number of Taxa", x = "Site") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxa Counts Per Site, eDNA") +
  coord_flip()
print(p)

# REEF surveys
q <- ggplot(data=merged_reef_survey, aes(x=reorder(Location, -table(Location)[Location]))) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  labs(y = "Number of Taxa", x = "Site") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Taxa Counts Per Site, REEF Survey") +
  coord_flip()
print(q)

# Double check eDNA and REEF have same site names
unique_eDNA_sites <- sort(unique(finalcleaned_eDNA_df$Site))
unique_REEF_sites <- sort(unique(merged_reef_survey$Location))

# Bar plot of site data for both eDNA and REEF side by side
grid.arrange(p, q, ncol = 2)

# Make data frames with frequency of sightings / detections by site
eDNA_sitefreq_df <- finalcleaned_eDNA_df %>%
  group_by(Site_ID) %>%
  summarise(frequency = n())

REEF_sitefreq_df <- merged_reef_survey %>%
  group_by(Site_ID) %>%
  summarise(frequency = n())

# Bar plot combining these two things
site_barplot <- ggplot() +
  geom_bar(data = eDNA_sitefreq_df, aes(x = Site_ID, y = frequency, fill = "eDNA"),
           position = position_dodge(width = 0.8), stat = "identity") +
  geom_bar(data = REEF_sitefreq_df, aes(x = Site_ID, y = -frequency, fill = "REEF"),
           position = position_dodge(width = 0.8), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  labs(y = "Frequency", x = "Site_ID") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ASV and Observed Species Frequencies per Site") +
  scale_fill_manual(values = c("eDNA" = "blue", "REEF" = "red")) +
  coord_flip()

print(site_barplot)

# Plotting %s of frequencies to standardize y-axis scales
# Calculate total frequencies
total_frequency_eDNA <- sum(eDNA_sitefreq_df$frequency)
total_frequency_REEF <- sum(REEF_sitefreq_df$frequency)

# Calculate percentages
eDNA_sitefreq_df$percentage <- (eDNA_sitefreq_df$frequency / total_frequency_eDNA) * 100
REEF_sitefreq_df$percentage <- (REEF_sitefreq_df$frequency / total_frequency_REEF) * 100

# Plotting
site_barplot_percentages <- ggplot() +
  geom_col(data = eDNA_sitefreq_df, aes(x = Site_ID, y = -percentage, fill = "eDNA"),
           position = position_dodge(width = 0.8)) +
  geom_col(data = REEF_sitefreq_df, aes(x = Site_ID, y = percentage, fill = "REEF"),
           position = position_dodge(width = 0.8)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  labs(y = "Percentage", x = "Site_ID") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ASV and Observed Species Frequencies per Site, by % of Frequencies") +
  scale_fill_manual(values = c("eDNA" = "blue", "REEF" = "red")) +
  coord_flip()

print(site_barplot_percentages)

# Reorder from most to least %s in REEF data
# Order Site_ID based on REEF frequencies
order_sites <- order(REEF_sitefreq_df$frequency)

# Reorder Site_ID factor levels
eDNA_sitefreq_df$Site_ID <- factor(eDNA_sitefreq_df$Site_ID, levels = eDNA_sitefreq_df$Site_ID[order_sites])

# Plotting
site_barplot_percentages_ordered <- ggplot() +
  geom_col(data = REEF_sitefreq_df, aes(x = reorder(Site_ID, percentage), y = percentage, fill = "REEF"),
           position = position_dodge(width = 0.8)) +
  geom_col(data = eDNA_sitefreq_df, aes(x = reorder(Site_ID, -percentage), y = -percentage, fill = "eDNA"),
           position = position_dodge(width = 0.8)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
  labs(y = "Percentage", x = "Site_ID") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("ASV and Observed Species Frequencies per Site, by % of Frequencies") +
  scale_fill_manual(values = c("eDNA" = "blue", "REEF" = "red")) +
  coord_flip()

print(site_barplot_percentages_ordered)

############################################################################################
# NMDS
############################################################################################
### eDNA ###
# First make a dataframe where each row is a sample and each column is an ASV, and numbers in
# each cell are nReads
nmds_eDNA_subset <- finalcleaned_eDNA_df[, c("sample", "Hash", "nReads")]
nmds_matrix_eDNA <- pivot_wider(data = nmds_eDNA_subset,
                                names_from = Hash,
                                values_from = nReads,
                                values_fill = 0)  # replace NAs with 0s

# Transform data using Wisconsin double-standardization
nmds_matrix_eDNA_standardized <- decostand(nmds_matrix_eDNA, "max", na.rm=FALSE)

# Transform data manually
max_matrix <- max(nmds_matrix_eDNA)
normalized_nmds_matrix_eDNA <- nmds_matrix_eDNA / max_matrix

# Calculate the dissimilarity matrix
diss_matrix_eDNA <- vegdist(normalized_nmds_matrix_eDNA, method = "bray")

# Perform & plot NMDS
nmds_result_eDNA <- metaMDS(diss_matrix_eDNA, k = 2)
plot_nmds_eDNA <- plot(nmds_result_eDNA)

# Add metadata to plot
# Make new dataframe with sample & site info
nmds_site_eDNA <- finalcleaned_eDNA_df %>%
  select(sample, Site_ID) %>%
  distinct() %>%
  arrange(sample) %>%
  mutate(sample = as.character(sample))
# nmds_coords <- scores(nmds_result, display = "sites")

scores(nmds_result_eDNA) %>%
  as_tibble(rownames="sample") %>%
  inner_join(., nmds_site_eDNA, by = "sample") %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=Site_ID)) +
    geom_point()+
  theme_classic()

### REEF ###
# Take the data I need for NMDS
# nmds_REEF_subset <- merged_reef_survey[, c("Form", "scientificname", "Abundance")]

# Filter data to retain rows with the highest Abundance within each group (bc there are duplicates
# for scientificname and Form due to separate entries for juveniles / adults of certain species)
# nmds_REEF_subset <- merged_reef_survey %>%
#   group_by(Form, scientificname) %>%
#   arrange(desc(Abundance)) %>%
#   slice(1) %>%
#   ungroup() %>%
#   select(Form, scientificname, Abundance)
# 
# # Convert Form and Abundance to appropriate data types
# nmds_REEF_subset$Form <- as.numeric(nmds_REEF_subset$Form)
# nmds_REEF_subset$Abundance <- as.integer(nmds_REEF_subset$Abundance)
# 
# # Now make a dataframe where each row is a form, each column is a species, and numbers in
# # each cell are Abundance
# nmds_matrix_REEF <- pivot_wider(data = nmds_REEF_subset,
#                                 names_from = scientificname,
#                                 values_from = Abundance,
#                                 values_fill = 0)  # replace NAs with 0s
# 
# # Calculate the dissimilarity matrix
# diss_matrix_REEF <- vegdist(nmds_matrix_REEF, method = "bray")
# 
# # Perform & plot NMDS
# nmds_result_REEF <- metaMDS(diss_matrix_REEF, k = 2)
# plot_nmds_REEF <- plot(nmds_result_REEF)
# 
# # Add metadata to plot
# # Make new dataframe with sample & site info
# nmds_site_REEF <- merged_reef_survey %>%
#   select(Form, Site_ID) %>%
#   distinct() %>%
#   arrange(Form) %>%
#   mutate(Form = as.character(Form))
# # nmds_coords <- scores(nmds_result, display = "sites")
# 
# scores(nmds_result_REEF) %>%
#   as_tibble(rownames="Form") %>%
#   inner_join(., nmds_site_REEF, by = "Form") %>%
#   ggplot(aes(x=NMDS1, y=NMDS2, color=Site_ID)) +
#   geom_point()+
#   theme_classic()

############################################################################################
# MAP
############################################################################################

# # Make a map of sites
# data("coastlineWorldFine")
# coast = as.data.frame(coastlineWorldFine@data)
# bathy.raw <- getNOAA.bathy(Lon1 = -123.3, Lon2 = 117,
#                              lat1 = 31.4, lat2 = 34.9, resolution = 1)
# bathy <- as_tibble(fortify.bathy(bathy.raw)) #add bathymetry (Lily's code)
# 
# library(ggplot2)
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthdata)
# 
# # Load coastline data
# world <- ne_countries(scale = "medium", returnclass = "sf")
# baja_coast <- world[world$name == "Mexico" & world$admin == "Baja California", ]
# 
# # Specify bounding box for bathymetry data
# lon1 <- -123.3
# lon2 <- 117
# lat1 <- 31.4
# lat2 <- 34.9
# 
# # Get bathymetry data
# bathy.raw <- getNOAA.bathy(lon1, lon2, lat1, lat2, resolution = 1)
# bathy.df <- fortify(bathy.raw)
# 
# # Plot coastline and bathymetry
# ggplot() +
#   geom_sf(data = baja_coast, fill = NA, color = "black") +
#   geom_contour(data = bathy.df, aes(x = x, y = y, z = z)) +
#   scale_fill_viridis_c(name = "Depth") +  # You can choose a different color scale
#   labs(title = "Baja California Coastline with Bathymetry") +
#   theme_minimal()
