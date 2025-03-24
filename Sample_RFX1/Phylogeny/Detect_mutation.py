from Bio import AlignIO
from collections import Counter

# Load MSA file
msa_file = "DBD_MSA_fixed.fasta"
alignment = AlignIO.read(msa_file, "fasta")

# Get sequence length
seq_length = alignment.get_alignment_length()

# Compute consensus sequence (most common residues at each position)
consensus_seq = ""
for i in range(seq_length):
    column = alignment[:, i]  # Extract residues at position i
    freq = Counter(column)     # Count occurrences of each residue
    most_common_residue, _ = freq.most_common(1)[0]  # Get most common residue
    consensus_seq += most_common_residue

# Store conservation data
conserved = []
variable_sites = {}
classified_mutations = []

# Define amino acid properties
bulky_residues = {"W": "Tryptophan (Bulky, Aromatic, Hydrophobic)",
                  "F": "Phenylalanine (Bulky, Aromatic, Hydrophobic)",
                  "Y": "Tyrosine (Bulky, Aromatic, Polar)",
                  "L": "Leucine (Bulky, Aliphatic, Hydrophobic)",
                  "I": "Isoleucine (Bulky, Aliphatic, Hydrophobic)"}

positive_residues = {"K": "Lysine (Positively charged, Basic, Hydrophilic)",
                     "R": "Arginine (Positively charged, Basic, Hydrophilic)"}

negative_residues = {"D": "Aspartic Acid (Negatively charged, Acidic, Hydrophilic)",
                     "E": "Glutamic Acid (Negatively charged, Acidic, Hydrophilic)"}

polar_residues = {"S": "Serine (Polar, Uncharged, Hydrophilic)",
                  "T": "Threonine (Polar, Uncharged, Hydrophilic)",
                  "N": "Asparagine (Polar, Uncharged, Hydrophilic)",
                  "Q": "Glutamine (Polar, Uncharged, Hydrophilic)",
                  "H": "Histidine (Polar, Weakly Basic, Can be Positively Charged)"}

# Analyze each column (position in MSA)
for i in range(seq_length):
    column = alignment[:, i]  # Extract residues at position i
    freq = Counter(column)     # Count occurrences of each residue
    consensus_residue = consensus_seq[i]
    
    if len(freq) == 1:  # All sequences have the same residue
        conserved.append(consensus_residue)
    else:  # Variable site (mutation)
        variable_sites[i + 1] = dict(freq)  # Store site & mutations
        
        for residue in freq.keys():
            if residue != consensus_residue:
                mutation_type = "General"
                
                if consensus_residue in bulky_residues and residue not in bulky_residues:
                    mutation_type = "Steric Clash"
                elif consensus_residue in positive_residues and residue in negative_residues:
                    mutation_type = "Electrostatic Repulsion"
                elif consensus_residue in negative_residues and residue in positive_residues:
                    mutation_type = "Electrostatic Repulsion"
                elif consensus_residue in polar_residues and residue not in polar_residues:
                    mutation_type = "Unpaired Polar Atom"
                
                classified_mutations.append(f"Position {i+1}: {consensus_residue} -> {residue} ({mutation_type})")

# Output conserved residues
print("\n✅ Conserved Residues:")
print("".join(conserved))

# Output variable mutation sites
print("\n🔬 Variable Residues (Mutations Across Species):")
for pos, residues in variable_sites.items():
    print(f"Position {pos}: {residues}")

# Output classified mutations
print("\n🔍 Classified Mutations:")
for mutation in classified_mutations:
    print(mutation)

# Save mutation sites to a file
with open("mutation_list.txt", "w") as f:
    f.write("Position\tMutation\tType\n")
    for mutation in classified_mutations:
        f.write(mutation + "\n")

print("\n📄 Mutation analysis saved to mutation_list.txt")
