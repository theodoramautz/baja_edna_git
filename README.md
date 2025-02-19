# baja_edna_git
Baja eDNA project git repository for code. 

In the code folder:
"baja_data_cleaning.R" "reef_data_cleaning.R" "baja_metadata.R" and "baja_analysis.R" were early analysis attempts and mostly go through some preliminary data cleaning steps. Run those first to create the .RData files referenced in "theodora_phyloseq_newseqs.Rmd"
"dada2_tutorial-v3-Baja.Rmd" is the code to run all the bioinformatics.
"theodora_phyloseq_newseqs.Rmd" is the main coding file with all analyses used in the paper.

In the data folder:
"baja_edna_metadata.csv" is the metadata file.
"eDNA_data" folder has:
"Bajalib_nextera_ud_indexes.xlsx - Baja library index.csv" which connects "sample" to "Seq_ID" which is how the sequencer identified each sample.
"diver_number_key.xlsx" helps associate a diver number with a diver name for divers who carried PEDS on their dives.
"eDNA Explore Baja REEF 2022 - All Data.csv" has metadata for the PEDS.

"reef_data" folder has:
"ExploreBaja22.csv" which is all raw RVS data.
Feel free to read the "readme.txt" file in that folder for further info on other data files.

"reference_databases" folder has:
"12S_efc_derep_and_clean-Baja-v3.fa" - the final reference database used to assign taxonomy.

"rerun_bioinformatics_newseqs" folder has:
All files for phyloseq analysis from final bioinformatics analyses (with updated reference database after adding fish tissue sequences).
