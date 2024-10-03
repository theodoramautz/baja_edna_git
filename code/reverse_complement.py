from Bio import SeqIO
from Bio.Seq import Seq

def create_reverse_complement_fasta(input_fasta, output_fasta):
    with open(output_fasta, 'w') as output_handle:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            # Generate reverse complement
            rev_comp_seq = record.seq.reverse_complement()
            # Create a new record with the reverse complement
            rev_comp_record = record
            rev_comp_record.seq = rev_comp_seq
            # Write the reverse complement to the output file
            SeqIO.write(rev_comp_record, output_handle, 'fasta')

# Example usage
create_reverse_complement_fasta('Sequenced Fish v3.fasta', 'Sequenced Fish Rev Complement v2.fasta')