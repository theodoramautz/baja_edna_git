##############################################################################
##                         Welcome to The Kelly Lab's                       ##
##              Metabarcoding Taxonomic Assignment Pipeline                 ##
## This pipeline will take you from de-multiplexed fastq files to a matrix  ##
## with taxon names and counts. It should work for linux or mac. PC users,  ##
## look into WSL2.                                                          ##
##############################################################################
##    Before you start, ensure you have access to the following files.      ##
##            Many can be found in the OneDrive under:                      ##
##             "2. KellyLab/Projects/AnnotationDatabases/code/"             ##
##                                                                          ##
## Required Files:                                                          ##
## 1. Your Fastqs, of course. Their names must start with the primer name   ##
## 2. CEG_BLAST_function.R - Function to run BLAST on the CEG cluster       ##
## 3. LCA_function.R - Function to run LCA                                  ##
## 4. Pipeline_code.R - This very code                                      ##
## 5. MURI_taxids_to_exclude.txt #optional                                  ##
## 6. MURIblast_*_template.sh - Code for taxonomic assignment               ##
##    depending on delimitation                                             ##
## 7. make_primer_shell_script.R - A script that writes another script      ##
## 8. primer_data.csv - Sheet with primer information like sequence, etc    ##
##                                                                          ##
## Required Programs:                                                       ##
## 1. Taxonkit: https://bioinf.shenwei.me/taxonkit/                         ##
## 2. Cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html ##
##############################################################################

#Restart your Rstudio session before running each time.
.rs.restartR()
#You may also clean your environment if you want to. 
#This solves many bugs you may have by rerunning the pipeline back to back.
##Be careful when running this
rm(list = ls())

##put the path to this very code here.
here::i_am("pipeline_code_3.0-NVP-2.0-Baja-DL.R")

#load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(dada2))
suppressMessages(library(digest))
suppressMessages(library(seqinr))
suppressMessages(library(sys))
suppressMessages(library(ShortRead)) 
suppressMessages(library(here))

### SET DIRECTORY PATHS

# Your R working directory. This will change again downstream 
setwd(here("/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Baja/DL/DL-MURI/"))

# Original source code from which scripts can be copied and modified
SOURCE_CODE_LOCATION <- "/Users/nastassiapatin/OneDrive - SCCWRP/Bioinformatics/MURI metabarcoding/Pipeline/"

#the folder in which a subfolder, called `Fastq` lives, and which contains only fastq or fastq.gz files
PARENT_LOCATION <- "/Volumes/CalCOFI_OMICS/Sequence_data/GGBC_NextSeq_111324"  
  if (strsplit(PARENT_LOCATION,"")[[1]][nchar(PARENT_LOCATION)] != "/"){PARENT_LOCATION = paste0(PARENT_LOCATION, "/")} #ensures a trailing "/" as the last character on the path here

#then, ultimately, where do you want the small handful of relevant output files to live?
PROCESSED_LOCATION <- "/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Baja/DL/" # laptop

### SET VARIABLES

RUN_NAME <- "BAJA-DL" #basename(PARENT_LOCATION)

PRIMERNAME <- "DL" #used to look up primer information, so needs to be one of a set number of known primers

BLASTSCOPE <- "eukaryote" #  Default is 97% identity for vertebrate, 90% for eukaryote. Euk can take the masking, vertebrate cannot

#read in primer data
primer.data <- read.csv(paste0(SOURCE_CODE_LOCATION, "primer.data.csv"))

# primer trimming variables
PRIMERSEQ_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_f)
PRIMERSEQ_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(seq_r)
PRIMERLENGTH_F <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_f)
PRIMERLENGTH_R <- primer.data %>% filter(name == PRIMERNAME) %>% pull(primer_length_r)
MAX_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(max_amplicon_length)
MIN_AMPLICON_LENGTH <- primer.data %>% filter(name == PRIMERNAME) %>% pull(min_amplicon_length)
OVERLAP <- primer.data %>% filter(name == PRIMERNAME) %>% pull(overlap)

### IDENTIFY FILES

##File with the list of taxa to be excluded from Blast search (optional)
# Read the taxids from the file, assuming one taxid per line
NEGATIVE_TAXIDS_FILE <- paste0(SOURCE_CODE_LOCATION, "MURI_taxids_to_exclude.txt")

#where the hash database is. This is to check if sequences have been seen before and avoid double effort in taxon assignment
DATABASE_LOCATION <- "/Volumes/CalCOFI_OMICS/Databases/MURI hash DBs/"   #local
  if (strsplit(DATABASE_LOCATION,"")[[1]][nchar(DATABASE_LOCATION)] != "/"){DATABASE_LOCATION = paste0(DATABASE_LOCATION, "/")} #ensures a trailing "/" as the last character on the path here

#### SET UP FILE SYSTEM

#system2("mkdir", shQuote(paste0(PROCESSED_LOCATION))) #create folder for code, etc
FASTQ_LOCATION <- paste0(PARENT_LOCATION, "RAW") #folder within Parent Location
CODE_LOCATION <- paste0(PROCESSED_LOCATION, "code_etc") #folder within Parent Location
  system2("mkdir", shQuote(CODE_LOCATION)) #create folder for code, etc
TRIMMED_LOCATION <- paste0(PROCESSED_LOCATION, "cutadapt") #folder within Parent Location
FILTERED_LOCATION <- paste0(PROCESSED_LOCATION, "dada2_filtered") #folder within Parent Location
OUTPUT_LOCATION <- paste0(PROCESSED_LOCATION, "outputs") #folder within Parent Location
  system2("mkdir", args=shQuote(FILTERED_LOCATION)) #make folder for processed reads
  system2("mkdir", args=shQuote(OUTPUT_LOCATION)) #make folder for pipeline outputs

  
### EXTRA CODE
  
#check dependencies and get path for taxonkit
TAXONKIT_PATH = "/opt/miniconda3/bin/taxonkit"  # laptop

#check if taxonkit is working
system2(TAXONKIT_PATH, args=shQuote("version"))

# identify LCA source code
LCA_CODE <- paste0(SOURCE_CODE_LOCATION, "LCA_function_3.0-NVP.R")

### COPY FILES TO SCRIPTS FOLDER

write.csv(primer.data, paste0(CODE_LOCATION, "/primer.data.csv"))
system2("cp", args = c(shQuote(here(NEGATIVE_TAXIDS_FILE)), 
                       shQuote(CODE_LOCATION))) #taxa masked from blast search
system2("cp", args = c(shQuote(here(LCA_CODE)), shQuote(CODE_LOCATION))) # main pipeline
  source(paste0(CODE_LOCATION, "/", "LCA_function_3.0-NVP.R"))

if (BLASTSCOPE == "vertebrate"){
  system2("cp", args = c(shQuote(here("MURIblast_vertebrate_template.sh")), shQuote(SOURCE_CODE_LOCATION))) # Copy vertebrate template  
}
if (BLASTSCOPE == "eukaryote"){
  system2("cp", args = c(shQuote(here("MURIblast_eukaryote_template.sh")), shQuote(SOURCE_CODE_LOCATION))) # Copy eukaryote template  
}

#-----------------------------------------------------------------------------------------------
#Downstream of here, no further user input is required. 
#Therefore, you may select all lines up to "ANNOTATION" and run all at once if you wish to. 
#Good luck!-------------------------------------------------------------------------------------

###RUN DADA2---------------------------------------------------------------------

filelist <- system2("ls", args = shQuote(TRIMMED_LOCATION), stdout = TRUE)

fnFs <- filelist[str_detect(filelist, "_R1")] #filelist[grep(pattern="_R1_001.fastq", filelist)]
fnRs <- filelist[str_detect(filelist, "_R2")] #filelist[grep(pattern="_R2_001.fastq", filelist)]

sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- paste(sample.names1, sample.names2, sep = "_")

### Name filtered files in filtered/subdirectory ----------------------------------

filtFs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### Filter out Empty Samples -----------------------------------------------------
# if we don't filter out empty samples, an error happens during finding qual trimming length

setwd(TRIMMED_LOCATION)
file.empty <- function(filenames) file.info(filenames)$size == 20
empty_files <- file.empty(fnFs) | file.empty(fnRs)
fnFs <- fnFs[!empty_files]
fnRs <- fnRs[!empty_files]
filtFs <- filtFs[!empty_files]
filtRs <- filtRs[!empty_files]
sample.names <- sample.names[!empty_files]

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimRight = c(15, 15),
                     # truncLen = MAX_AMPLICON_LENGTH,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE) #try single thread if multithread gives errors

### Filter ---------------------------------------------------------------

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
sample.names <- sample.names[exists]

##Learn error rates---------------------------------------------------------------

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

### Dereplicate and Learn Error Rates ---------------------------------------------------------------

dadaFs <- dada(filtFs, err=errF, selfConsist = TRUE, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, selfConsist = TRUE, multithread=TRUE)

### Merge Paired Reads ---------------------------------------------------------------------

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, 
                      minOverlap = OVERLAP, verbose=TRUE, trimOverhang = TRUE)

### Construct sequence table ---------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)

### Remove chimeras ------------------------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Filter by Size -------------------------------------------------------------------------

indexes.to.keep <- which((nchar(colnames(seqtab.nochim)) <= MAX_AMPLICON_LENGTH) & 
                           ((nchar(colnames(seqtab.nochim))) >= MIN_AMPLICON_LENGTH))
cleaned.seqtab.nochim <- seqtab.nochim[,indexes.to.keep]
filteredout.seqtab.nochim <- seqtab.nochim[,!indexes.to.keep]
write.csv(filteredout.seqtab.nochim,paste0(OUTPUT_LOCATION,
                                           "/filtered_out_asv.csv"))

### Track reads through pipeline -----------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out[exists,], sapply(dadaFs, getN), 
               sapply(dadaRs, getN), sapply(mergers, getN), 
               rowSums(seqtab.nochim),
               rowSums(cleaned.seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","length_filter")
rownames(track) <- sample.names
head(track)
write.csv(track,paste0(OUTPUT_LOCATION,"/tracking_reads.csv"))

### Create Hashing  ------------------------------------------------------------------------

# define output files
conv_file <- paste0(OUTPUT_LOCATION, "/", 
                    paste0(RUN_NAME,"_",PRIMERNAME,"_hash_key.csv"))
conv_file.fasta <- file.path(OUTPUT_LOCATION, "/",
                             paste0(RUN_NAME,"_",PRIMERNAME,"_hash_key.fasta"))
ASV_file <-  file.path(OUTPUT_LOCATION, "/",
                       paste0(RUN_NAME,"_",PRIMERNAME,"_ASV_table.csv"))
taxonomy_file <- file.path(OUTPUT_LOCATION, "/",
                           paste0(RUN_NAME,"_",PRIMERNAME,"_taxonomy_output.csv"))
bootstrap_file <- file.path(OUTPUT_LOCATION, "/",
                            paste0(RUN_NAME,"_",PRIMERNAME,"_tax_bootstrap.csv"))

# create ASV table and hash key 
print(paste0("creating ASV table and hash key...", Sys.time()))
seqtab.nochim.df <- as.data.frame(cleaned.seqtab.nochim)
conv_table <- tibble( Hash = "", Sequence ="")
Hashes <- map_chr (colnames(seqtab.nochim.df), ~ digest(.x, algo = "sha1", serialize = F, skip = "auto"))
conv_table <- tibble (Hash = Hashes,
                      Sequence = colnames(seqtab.nochim.df))

write_csv(conv_table, conv_file) # write hash key into a csv
write.fasta(sequences = as.list(conv_table$Sequence), # write hash key into a fasta
            names     = as.list(conv_table$Hash),
            file.out = conv_file.fasta)
sample.df <- tibble::rownames_to_column(seqtab.nochim.df,"Sample_name")
sample.df <- data.frame(append(sample.df,c(Label=PRIMERNAME), after = 1))
current_asv <- bind_cols(sample.df %>%
                           dplyr::select(Sample_name, Label),
                         seqtab.nochim.df)
current_asv <- current_asv %>%
  pivot_longer(cols = c(- Sample_name, - Label),
               names_to = "Sequence",
               values_to = "nReads") %>%
  filter(nReads > 0)
current_asv <- merge(current_asv,conv_table, by="Sequence") %>%
  select(-Sequence) %>%
  relocate(Hash, .after=Label)

write_csv(current_asv, ASV_file) # write asv table into a csv

###ANNOTATION------------------------------------------------------------------------------------------

##Looking at the existing hash database to check if any observed sequences have been seen before.
if (file.exists(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))){
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
  seen <- which(conv_table$Hash %in% db$Hash) #hashes in this run that already occur in the database 
  notseen <- which(!conv_table$Hash %in% db$Hash) #hashes in this run that have not previously been annotated
  
  #write new (as-yet-unannotated sequences) to fasta
  write.fasta(sequences = as.list(conv_table$Sequence[notseen]), # write hash key into a fasta
              names     = as.list(conv_table$Hash[notseen]),
              file.out = paste0(OUTPUT_LOCATION, "/seqs_to_annotate.fasta"))
#if you want to blast everything instead, only run these lines below
  } else {
  write.fasta(sequences = as.list(conv_table$Sequence), # write hash key into a fasta
              names     = as.list(conv_table$Hash),
              file.out = paste0(OUTPUT_LOCATION, "/seqs_to_annotate.fasta"))
}

### Local BLAST against nt_euk database ###
system2("sh", args = shQuote(paste0(CODE_LOCATION, 
                                    "/MURIblast_vertebrate_template-NVP-localnteuk.sh")))

# Import BLAST outputs and sequence fasta
BLASTOUTPUT = paste0(OUTPUT_LOCATION, "/", "BLAST_output.txt")
FASTA = paste0(OUTPUT_LOCATION, "/", "seqs_to_annotate.fasta")

##Taxonomic assignment
##There may be warnings when taxa is not in the taxonomic database. You may ignore.
if (file.size(BLASTOUTPUT) == 0L){
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
} else if (file.exists(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))){
  #LCA table for output, and add to database
  LCA(BLASTOUTPUT = BLASTOUTPUT,
      FASTA = FASTA,
      DB_PATH_IN = paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"),
      DB_PATH_OUT = paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))
  
  # re-load newly updated db
  db <- read.csv(paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"), row.names = 1)
#if you have blasted everything again instead, only run these lines below  
  } else {
    db <- LCA(BLASTOUTPUT,
              FASTA)
    write.csv(db %>% distinct(), paste0(DATABASE_LOCATION, PRIMERNAME, "_database.csv"))
}

##Writing the output taxon tables
tax_table <- current_asv %>%
  left_join(db %>% dplyr::select(Hash, BestTaxon, Class)) %>%
  group_by(Sample_name, BestTaxon, Class) %>%
  summarise(nReads = sum(nReads))

# Tax table with full taxonomic lineage
tax_table_full <- current_asv %>%
  left_join(db %>% dplyr::select(Hash, BestTaxon, Kingdom, Phylum, 
                                 Class, Order, Family, Genus, Species)) %>%
  group_by(Sample_name, BestTaxon, Kingdom, Phylum, 
           Class, Order, Family, Genus, Species) %>%
  summarise(nReads = sum(nReads))

write_csv(tax_table, file = paste0(OUTPUT_LOCATION, "/", "taxon_table.csv"))
write_csv(tax_table_full, 
          file = paste0(OUTPUT_LOCATION, "/", "taxon_table_lineage.csv"))

tax_table %>%
  pivot_wider(id_cols = c(BestTaxon, Class), names_from = Sample_name, values_from = nReads, values_fill = 0) %>% 
  write_csv(file = paste0(OUTPUT_LOCATION, "/", "taxon_table_wide.csv"))

tax_table_full %>%
  pivot_wider(id_cols = c(BestTaxon, Kingdom, Phylum, 
                          Class, Order, Family, Genus, Species), 
              names_from = Sample_name, values_from = nReads, 
              values_fill = 0) %>% 
  write_csv(file = paste0(OUTPUT_LOCATION, "/", "taxon_table_lineage_wide.csv"))

#A taxon table but that lists variants/haplotypes separately
merged_data <- current_asv %>%
  left_join(db, by = "Hash")
merged_data <- merged_data %>%
  mutate(BestTaxon_Haplotype = paste(BestTaxon, HaplotypeNumber, sep = "_"))
summarized_data <- merged_data %>%
  group_by(BestTaxon_Haplotype, Sample_name, Class) %>%
  summarise(nReads = sum(nReads, na.rm = TRUE), .groups = 'drop')
reshaped_data <- summarized_data %>%
  pivot_wider(id_cols = c(BestTaxon_Haplotype, Class), names_from = Sample_name, values_from = nReads, values_fill = 0)
  
write.csv(reshaped_data, paste0(OUTPUT_LOCATION, "/", "haplotype_table.csv"), row.names = FALSE)

message("Pipeline complete. Outputs are now available in the Processed Location.")
