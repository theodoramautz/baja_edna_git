import csv
from Bio import Entrez

# Set your email (required by NCBI)
Entrez.email = "thmautz@ucsd.edu"

# Define the specific taxonomic levels you want
desired_ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

def get_taxonomic_info(species_name):
    try:
        # Search NCBI taxonomy database for the species name
        handle = Entrez.esearch(db="taxonomy", term=species_name)
        record = Entrez.read(handle)
        
        if not record["IdList"]:
            return f"TaxID not found for {species_name}"
        
        # Extract the TaxID from the search result
        taxid = record["IdList"][0]
        
        # Fetch the full taxonomic information using the TaxID
        handle = Entrez.efetch(db="taxonomy", id=taxid)
        record = Entrez.read(handle)
        
        # Filter the lineage by the desired taxonomic levels
        lineage_dict = {item["Rank"]: item["ScientificName"] for item in record[0]["LineageEx"] if item["Rank"] in desired_ranks}
        lineage = [lineage_dict.get(rank, "") for rank in desired_ranks]
        lineage.append(record[0]["ScientificName"])  # Add species name
        
        # Join the filtered lineage into a single string separated by semicolons
        return ";".join(lineage)
    except Exception as e:
        return f"Error retrieving taxonomy for {species_name}: {str(e)}"

def process_species_from_csv(input_csv, output_csv):
    # Open the input CSV file and read species names
    with open(input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        species_list = [row[0] for row in reader]
    
    # Open the output CSV file to write the results
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Species Name", "Taxonomic Information"])
        
        # Loop through the species and fetch their taxonomic information
        for species in species_list:
            taxonomy = get_taxonomic_info(species)
            writer.writerow([species, taxonomy])
            print(f"Processed: {species}")

# File paths
input_csv = "species.csv"  # Replace with your input CSV file path
output_csv = "taxonomy_output2.csv"  # Replace with your desired output CSV file path

# Process species from the CSV file
process_species_from_csv(input_csv, output_csv)
