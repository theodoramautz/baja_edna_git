# This workflow was used to run Baja dloop amplicons thorugh dada2 
# This code assumes you already have primers trimmed off in a folder named "cutadapt"
# Nastassia V. Patin
# February 20 2025

#load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(dada2))
suppressMessages(library(digest))
suppressMessages(library(seqinr))
suppressMessages(library(sys))
suppressMessages(library(ShortRead)) 
suppressMessages(library(here))

### SET DIRECTORY PATHS

# Your R working directory
setwd(here("/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Baja/DL/DL-MURI/"))

### SET VARIABLES

RUN_NAME <- "BAJA-DL"
PRIMER_NAME <- "DL"

#### SET UP FILE SYSTEM

# directory for processed files

PROCESSED_LOCATION <- "/Users/nastassiapatin/OneDrive - SCCWRP/Experiments/Baja/DL/" # laptop
CODE_LOCATION <- paste0(PROCESSED_LOCATION, "code_etc") #folder within Parent Location
  system2("mkdir", shQuote(CODE_LOCATION)) #create folder for code, etc
TRIMMED_LOCATION <- paste0(PROCESSED_LOCATION, "cutadapt") #folder within Parent Location
FILTERED_LOCATION <- paste0(PROCESSED_LOCATION, "dada2_filtered") #folder within Parent Location
OUTPUT_LOCATION <- paste0(PROCESSED_LOCATION, "outputs") #folder within Parent Location
  system2("mkdir", args=shQuote(FILTERED_LOCATION)) #make folder for processed reads
  system2("mkdir", args=shQuote(OUTPUT_LOCATION)) #make folder for pipeline outputs

### RUN DADA2---------------------------------------------------------------------

filelist <- system2("ls", args = shQuote(TRIMMED_LOCATION), stdout = TRUE)

fnFs <- filelist[str_detect(filelist, "_R1")] 
fnRs <- filelist[str_detect(filelist, "_R2")] 

sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- paste(sample.names1, sample.names2, sep = "_")

### Name filtered files in filtered/subdirectory ----------------------------------

filtFs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(FILTERED_LOCATION, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### Filter out empty samples -----------------------------------------------------

setwd(TRIMMED_LOCATION)
file.empty <- function(filenames) file.info(filenames)$size == 20
empty_files <- file.empty(fnFs) | file.empty(fnRs)
fnFs <- fnFs[!empty_files]
fnRs <- fnRs[!empty_files]
filtFs <- filtFs[!empty_files]
filtRs <- filtRs[!empty_files]
sample.names <- sample.names[!empty_files]

### Filter ---------------------------------------------------------------

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimRight = c(15, 15),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE) 

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
sample.names <- sample.names[exists]

### Learn error rates---------------------------------------------------------------

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

### Generate ASVs ---------------------------------------------------------------

dadaFs <- dada(filtFs, err=errF, selfConsist==TRUE, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, selfConsist=TRUE, multithread=TRUE)

### Merge Paired Reads ---------------------------------------------------------------------

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, 
                      minOverlap=12, verbose=TRUE, trimOverhang=TRUE)

### Construct sequence table ---------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)

### Remove chimeras ------------------------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Filter by Size -------------------------------------------------------------------------

# these values filter the amplicons by maximum and minimum lengths 
indexes.to.keep <- which((nchar(colnames(seqtab.nochim)) <= 600) & 
                           ((nchar(colnames(seqtab.nochim))) >= 240))
cleaned.seqtab.nochim <- seqtab.nochim[, indexes.to.keep]
filteredout.seqtab.nochim <- seqtab.nochim[, !indexes.to.keep]
write.csv(filteredout.seqtab.nochim, paste0(OUTPUT_LOCATION,
                                           "/filtered_out_asv.csv"))

### Track reads through pipeline -----------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out[exists,], sapply(dadaFs, getN), 
               sapply(dadaRs, getN), sapply(mergers, getN), 
               rowSums(seqtab.nochim),
               rowSums(cleaned.seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", 
                     "denoisedR", "merged", "nonchim","length_filter")
rownames(track) <- sample.names
head(track)
write.csv(track, paste0(OUTPUT_LOCATION, "/tracking_reads.csv"))

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
seqtab.nochim.df <- as.data.frame(cleaned.seqtab.nochim)
conv_table <- tibble( Hash = "", Sequence ="")
Hashes <- map_chr (colnames(seqtab.nochim.df), ~ digest(.x, algo = "sha1", serialize = F, skip = "auto"))
conv_table <- tibble (Hash = Hashes,
                      Sequence = colnames(seqtab.nochim.df))

write_csv(conv_table, conv_file) # write hash key into a csv
write.fasta(sequences = as.list(conv_table$Sequence), # write hash key into a fasta
            names     = as.list(conv_table$Hash),
            file.out = conv_file.fasta)
sample.df <- tibble::rownames_to_column(seqtab.nochim.df, "Sample_name")
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
