## eDNA Analysis in QIIME2
#### March 2024
Before you start, please look at these resources on the QIIME2 website:

https://docs.qiime2.org/2024.2/getting-started/  
https://docs.qiime2.org/2024.2/tutorials/overview/#useful-points-for-beginners


### 1. Installation
Back to our old friend conda! Open a Terminal window and update your miniconda:

```
conda update conda
```
Install wget (might already be installed)

```
conda install wget
```
Download the QIIME2 yaml file and install it in a new conda environment

```
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-osx-conda.yml
conda env create -n qiime2-amplicon-2024.2 --file qiime2-amplicon-2024.2-py38-osx-conda.yml
```

Activate the conda environment and test your installation

```
conda activate qiime2-amplicon-2024.2
qiime --help
```

### 2. Prepare your data!

This is the most important part of running QIIME2 smoothly. If you prepare your data and metadata, everything else will be a breeze!

#### Sequence files
Your sequence data files should be in fastq.gz format and live in their own folder.
#### Manifest
This document will tell QIIME2 where to find your raw fastq files. It should contain three columns: 'sample-id', 'forward-absolute-filepath', and 'reverse-absolute-filepath'. The 'filepath' columns should point to the full path of each corresponding forward or reverse read. The manifest file needs to be a tab-separated file (.tsv extension). You can make a .csv file in Excel, then open it in a text editor like Atom or BBEdit (NOT MS Word!), replace the commas with tabs, and save it with the .tsv extension. Make sure there are no weird characters at the beginning or end of the lines.  
For more information, read the section called '"Fastq manifest" formats' on the QIIME2 website: https://docs.qiime2.org/2024.2/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq
#### Metadata
The metadata file needs to be tab-separated and include the column 'sample-id' whose contents should match the 'sample-id' column in the manifest.

### 3. Import your sequence data

QIIME2 works with zipped files it calls "artifacts." These files have .qza extensions, or .qzv if they are meant to be visualized.

In this command, you will use the manifest file to provide the fastq file locations. The output will be a single file ("artifact") with a .qza extension. You could unzip this file and see its individual contents if you wanted to. 

This step can take a while, so be patient.

```
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Baja_QIIME2_manifest.tsv \
  --input-format PairedEndFastqManifestPhred33V2 \
  --output-path demux-paired-end-12S.qza
```

### 4. Trim primers from the sequences
We will use our old friend cutadapt to trim the reads, but this time we are using it within QIIME2. Just like with cutadapt alone, we will set the primer and reverse-complement primer sequences as variables.

```
primerF="GCCGGTAAAACTCGTGCCAGC" 
primerR="CATAGTGGGGTATCTAATCCCAGTTTG"
```

```
revcomp_primerF=`echo $primerF | \
tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni \
TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
revcomp_primerR=`echo $primerR | \
tr ACGTUWSMKRYBDHVNIacgtuwsmkrybdhvni \
TGCAAWSKMYRVHDBNNtgcaawskmyrvhdbnn | rev`
```
Now we can run cutadapt!

```
qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-12S.qza \
  --p-front-f "${primerF};required...${revcomp_primerR};optional" \
  --p-front-r "${primerR};required...${revcomp_primerF};optional" \
  --verbose --p-error-rate 0.1 --p-minimum-length 50 --o-trimmed-sequences demux-trimmed-12S.qza
```
** Bonus step: you can export the trimmed reads, run them through fastqc and multiqc, and see how they compare to the plots you made when you ran the reads through cutadapt as a standalone program.

To export files from QIIME2 artifacts:

```
qiime tools export --input-path demux-trimmed-12S.qza --output-path 12S-trimmed-reads
```

### 5. Denoise with DADA2

```
qiime dada2 denoise-paired --i-demultiplexed-seqs demux-trimmed-12S.qza \
  --p-trunc-len-f 110 --p-trunc-len-r 110 \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-max-ee-f 2 --p-max-ee-r 2 \
  --p-trunc-q 2 \
  --p-n-threads 4 \
  --o-table table.qza \
  --o-representative-sequences repseqs.qza \
  --o-denoising-stats dada2-stats.qza \
  --verbose
```
Convert the qiime artifact (.qza extension) to a visualization (.qzv extension) 

```
qiime metadata tabulate --m-input-file dada2-stats.qza \
  --o-visualization dada2-stats.qzv
```

Take a look at the stats

```
qiime tools view dada2-stats.qzv
```
This table is interactive, so you can sort by, e.g., percentage of non-chimeric reads. You can also download the table in .tsv format and use that table in R to make plots!

When you're done looking at the table, you can close the window and hit 'q' or 'Ctrl+C' to quit and return to the qiime2 commands. 

### 6. Explore the read DADA2 results

First look at the ASV table (or "feature table")

```
qiime feature-table summarize --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file Baja_QIIME2_metadata.tsv
```

```
qiime tools view table.qzv
```

As we saw when we ran DADA2 in R, the denoising and merging process cut out A LOT of reads! Let's get rid of samples with fewer than 100 reads in them.

```
qiime feature-table filter-samples \
  --i-table table.qza \
  --p-min-frequency 100 \
  --o-filtered-table table-filt.qza 
```
You can summarize and view this table like we did with the original one. How many samples did you lose with this filter?

You can also look at the individual ASV sequences.

```
qiime feature-table tabulate-seqs \
  --i-data repseqs.qza \
  --o-visualization repseqs.qzv
```

```
qiime tools view repseqs.qzv
```
Using both the ASV table and the ASV sequences (repseqs) can you figure out what taxon the most common ASV across all samples belongs to?

### 7. Assign taxonomy
Here you will use a QIIME2-formatted version of the 12S database we used in DADA2 in R. You can download the '12S_efc_rCRUX_20240307-classifier.qza' database file from the Google Drive.

```
qiime feature-classifier classify-sklearn \
  --i-classifier 12S_efc_rCRUX_20240307-classifier.qza \
  --i-reads repseqs.qza \
  --o-classification taxonomy.qza
```

```
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```
In the browser, you can click a button to download the taxonomy table as a tsv file if you wish.

### 8. Phylogeny
QIIME2 will run a suite of diversity analyses for you; some of these rely on a phylogeny of all ASVs, so we first need to align the sequences and make a tree.

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences repseqs.qza \
  --o-alignment aligned-repseqs.qza \
  --o-masked-alignment masked-aligned-repseqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

### 9. Alpha and beta diversity

```
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table-filt.qza \
  --p-sampling-depth 100 \
  --m-metadata-file Baja_QIIME2_metadata.tsv \
  --output-dir core-metrics-results
```
Many of these outputs aren't immediately available as visualizations, but some of the beta diversity plots are! 

```
qiime tools view core-metrics-results/bray_curtis_emperor.qzv
```