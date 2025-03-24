from ete3 import Tree, TextFace

# Load phylogenetic tree
tree_file = "DBD_phylogenetic_tree.nwk"
tree = Tree(tree_file, format=1)

# Load mutation sites from mutation analysis
mutation_file = "mutation_list.txt"
mutations = {}

with open(mutation_file, "r") as f:
    for line in f:
        if line.startswith("##") or line.strip() == "" or "Mutation Type" in line:
            continue  # Skip headers and empty lines
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            position = parts[1]  # Mutation position
            mutation_desc = parts[2]  # Mutation details
            mutations[position] = mutation_desc

# Annotate tree nodes with mutations
for leaf in tree.iter_leaves():
    species = leaf.name  # Extract species ID
    if species in mutations:
        mutation_text = mutations[species]
        leaf.add_face(TextFace(f" {mutation_text}", fsize=10, fgcolor="red"), column=0)

# Save annotated tree as an image
tree.render("mutation_tree.png", w=800, units="px")
print("📄 Mutation tree saved as mutation_tree.png")
