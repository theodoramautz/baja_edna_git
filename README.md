# baja_edna_git
Baja eDNA project git repository for code. 

In the code folder:
There are 2 folders. "mifish" has the code for analyses using MiFish primers. "dloop" has the code for analyses using D-loop primers.
In "mifish":
"peds_data_cleaning.R" and "rvs_data_cleaning.R" should be run first. These cleaned data prior to analyses.
"metadata_creation.R" should be run next to create the metadata file.
"mifish_bioinformatics.Rmd" is the code to run all the bioinformatics.
"mifish_analysis.Rmd" contains all analyses used in the paper for MiFish primers. This should be run after all other code files are run because it uses cleaned / processed data created from the others.

In the data folder:
"baja_edna_metadata.csv" is the cleaned metadata file.

"edna_data" folder has:
"Bajalib_nextera_ud_indexes.xlsx - Baja library index.csv" which connects "sample" to "Seq_ID" which is how the sequencer identified each sample.
"eDNA Explore Baja REEF 2022 - All Data.csv" is the original metadata file for the PEDS.
"MiFish" folder has:
All processed data from the bioinformatics workflow found in "mifish_bioinformatics.Rmd" from the code folder.
"dloop" folder has:
All processed data from the bioinformatics workflow found in "Baja_dloop_dada2.R" from the code folder.

"reef_data" folder has:
"ExploreBaja22.csv" which is all raw RVS data.
Feel free to read the "readme.txt" file in that folder for further info on other data files.

"reference_databases" folder has:
"12S_efc_derep_and_clean-Baja-v3.fa" - the final reference database used to assign taxonomy.
