import csv
import sys
sys.path.append('/Users/theodoramautz/Library/Python/3.9/lib/python/site-packages')
from Bio import SeqIO

def search_fasta(fasta_file, species_list):
    found_species = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        for species in species_list:
            if species.lower() in record.description.lower():
                found_species.add(species)
    return found_species

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python search_fasta_with_csv.py <fasta_file> <species_csv>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    species_csv = sys.argv[2]

    # Read species names from the CSV file
    species_list = []
    with open(species_csv, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            species_list.extend(row)

    found_species = search_fasta(fasta_file, species_list)

    if found_species:
        print("Found the following species in the .fasta file:")
        for species in found_species:
            print(species)
        print()

        # Print species from CSV file not found in .fasta file
        not_found_species = set(species_list) - found_species
        if not_found_species:
            print("The following species from the CSV file are not found in the .fasta file:")
            for species in not_found_species:
                print(species)
        else:
            print("All species from the CSV file are found in the .fasta file.")
    else:
        print("No species found in the .fasta file.")
