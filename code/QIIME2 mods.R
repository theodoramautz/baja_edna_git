# First removed leading 0s from "BAJA_001" etc separately

# Read the .tsv file
data <- read.table("Baja_QIIME2_manifest.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Rename the column from "sample.id" to "sample-id"
colnames(data)[which(colnames(data) == "sample.id")] <- "sample-id"
colnames(data)[which(colnames(data) == "forward.absolute.filepath")] <- "forward-absolute-filepath"
colnames(data)[which(colnames(data) == "reverse.absolute.filepath")] <- "reverse-absolute-filepath"

# Replace the path string in each column
data$'forward-absolute-filepath' <- gsub("/your/absolute/file/path/", "/Users/theodoramautz/Desktop/SIO/Research/baja_edna/baja_edna_git/data/baja_raw_sequences/", data$'forward-absolute-filepath')
data$'reverse-absolute-filepath' <- gsub("/your/absolute/file/path/", "/Users/theodoramautz/Desktop/SIO/Research/baja_edna/baja_edna_git/data/baja_raw_sequences/", data$'reverse-absolute-filepath')

# Get rid of BAJA_165, BAJA_166, BAJA_167, BAJA_168
sample_ids_to_remove <- c("BAJA_165", "BAJA_166", "BAJA_167", "BAJA_168")
data <- data[!data$'sample-id' %in% sample_ids_to_remove, ]


# Write the modified data to a new .tsv file
write.table(data, file="Baja_QIIME2_manifest.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
