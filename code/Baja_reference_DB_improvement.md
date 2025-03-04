# Baja reference sequence database improvement
### Nastassia V. Patin
### January 2024

#####1. Manually searched GenBank for 12S rRNA gene sequences belonging to fish species missing from the rCRUX database and retrieval of associated NCBI Accession numbers. Keywords used: [species name], '12S rRNA'

#####2. Downloaded sequences associated with NCBI accession numbers (42 total)

#####3. Quality control of downloaded sequences with QIIME2 [1] and RESCRIPt [2]

* Generated matching taxonomy file for new sequences to add using entrez_qiime.py
	
	```
     # approx 37MB
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -zxvf taxdump.tar.gz
    # approx 900MB
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    gunzip nucl_gb.accession2taxid.gz
	
	```
	
	```
	 python /Users/Home/GitHub/entrez_qiime/entrez_qiime.py \
	 -i new_baja_seqs.fa -o new_baja_seqs_taxonomy.tsv \
	 -n /Users/Home/databases/NCBI/TAXO \
	 -a /Users/Home/databases/NCBI/nucl_gb.accession2taxid \
	 -r superkingdom,phylum,class,order,family,genus,species
	```
	
* Manually added 'Feature ID' and 'Taxon' as column headers to new\_baja\_seqs\_taxonomy.tsv and import the tsv file into QIIME2.	
	
	```
	qiime tools import --type 'FeatureData[Taxonomy]' \
	--input-path new_baja_seqs_taxonomy.tsv --output-path new_baja_seqs_taxonomy.qza
	```
* Dereplicated the new Baja seqs.
	
	```
	qiime rescript dereplicate --i-sequences new_baja_seqs.qza \
	--i-taxa new_baja_seqs_taxonomy.qza --p-mode 'uniq' --p-threads 4 \
	--o-dereplicated-sequences new_baja_seqs_derep.qza \
	--o-dereplicated-taxa new_baja_seqs_taxonomy_derep.qza
	```
* Quality control of dereplicated sequences.
	
	```
	qiime rescript cull-seqs --i-sequences new_baja_seqs_derep.qza \
	--p-n-jobs 4 --p-num-degenerates 1 --p-homopolymer-length 8 \
	--o-clean-sequences new_baja_seqs_derep_cull.qza
	```
	This brought the sequence number from 42 down to 41.
	
* Extract sequence regions flanked by MiFish primers  
	
	```
	qiime feature-classifier extract-reads \
	--i-sequences new_baja_seqs_derep_cull.qza \
	--p-f-primer GCCGGTAAAACTCGTGCCAGC --p-r-primer CATAGTGGGGTATCTAATCCCAGTTTG \
	--p-min-length 140 --p-max-length 400 \
	--o-reads new_baja_seqs_derep_cull-trimmed.qza
	```
	 This extracted only 5 of the 41 sequences!  
	 
* Extracted sequence regions with 70% ID match to the rCRUX database ('12S\_efc\_derep\_and\_clean.fasta')
	
	```
	qiime tools import --type 'FeatureData[Sequence]' \
	--input-path 12S_derep_and_clean.fasta --output-path 12S_derep_and_clean.qza
	```
	
	```
	qiime rescript extract-seq-segments \
	--i-input-sequences new_baja_seqs_derep_cull.qza \
	--i-reference-segment-sequences 12S_derep_and_clean.qza --p-perc-identity 0.7 \
	--p-min-seq-len 140 --p-threads 4 \
	--o-extracted-sequence-segments new_baja_seqs-rescript.qza \
	--o-unmatched-sequences new_baja_seqs-unmatched.qza --verbose
   ```	
   This extracted 10 of the 41 sequences. However, there were sequences that were shorter than the (theoretically) minimum 140 bp from the above command. I used seqkit to filter out short sequences.
   
   ```
   seqkit seq -m 140  new_baja_seqs-rescript.fa  > new_baja_seqs-rescript-filt.fa
   ```
   This left us with 7 sequences!
   
#####4. Changed headers in new fasta file from NCBI Accession numbers to full taxonomy

* Use the tab-delimited taxonomy file with the Python script replace\_headers.py to make the 	replacement (make sure the 'Feature ID' column entires start with '>').
	
	```
	python ../replace_fasta_headers.py -i new_baja_seqs-rescript-filt.fa \
	-l new_baja_seqs_taxonomy_v2.tsv -o new_baja_seqs-rescript-filt-MURI.fa
   ```
   
#####5. Add the new sequences to the rCRUX database
* Combine the new filtered sequence fasta and the existing reference DB fasta files
	
	```
	cat new_baja_seqs-rescript-filt-MURI.fa 12S_efc_derep_and_clean-MURI.fasta > \ 
	12S_efc_derep_and_clean-MURI-v2.fa
	```
	
#####6. References

[1]. Bolyen et al. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology. 37: 832-857. https://doi.org/10.1038/s41587-019-0209-9
[2]. Robeson et al. 2021. RESCRIPt: Reproducible sequence taxonomy reference database management. PLoS Computational Biology 17 (11): e1009581. https://doi.org/10.1371/journal.pcbi.1009581