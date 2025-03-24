# PyMOL script to visualize mutations on a TF structure
import sys
sys.path.append("/home/zero/pymol/lib/python3.8/site-packages")  # Ensure PyMOL's Python path is included
from pymol import cmd

def visualize_mutations(pdb_file, mutation_list_file):
    # Load the TF structure
    cmd.load(pdb_file, "TF_structure")
    
    # Print available residue numbers in PyMOL for debugging
    print("🔍 Available Residue Numbers in PDB:")
    cmd.iterate("all", "print(resi)")
    
    # Read mutation sites from mutation_sites.txt
    mutation_sites = []
    with open(mutation_list_file, "r") as f:
        for line in f.readlines()[1:]:  # Skip header
            line = line.strip()
            if line.startswith("Position"):
                parts = line.split(":")  # Extract the part after "Position"
                if len(parts) > 1:
                    pos = parts[0].replace("Position", "").strip()  # Extract number
                    if pos.isdigit():  # Ensure it's a valid number
                        mutation_sites.append(pos)
    
    # Check if valid mutation sites were extracted
    if not mutation_sites:
        print("⚠️ No valid mutation sites found! Please check mutation_sites.txt format.")
        return
    
    # Select and color mutation sites
    for site in mutation_sites:
        print(f"🔹 Attempting to select: resi {site}")
        cmd.select(f"mutation_site_{site}", f"resi {site}")
        cmd.show("sticks", f"mutation_site_{site}")
        cmd.color("red", f"mutation_site_{site}")
    
    # Zoom into the mutations
    cmd.zoom("mutation_site_*")
    
    # Save image automatically
    cmd.png("mutation_visualization.png", width=1200, height=900, dpi=300)
    print("📸 Image saved as mutation_visualization.png")
    
    # Save modified PDB with mutations highlighted
    cmd.save("TF_structure_mutations.pdb")
    print("📄 Modified PDB saved as TF_structure_mutations.pdb")
    
    print("✅ Mutations visualized in PyMOL!")
    print("📌 Save the image manually or export the PDB from PyMOL if needed.")
    
# Example usage
if __name__ == "__main__":
    sys.path.append("/home/zero/pymol/bin")  # Add PyMOL binary path
    pdb_file = "TF_structure.pdb"  # Make sure this file exists
    mutation_list_file = "mutation_list.txt"  # Ensure this file is generated from Step 3
    visualize_mutations(pdb_file, mutation_list_file)
