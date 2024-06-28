from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
import matplotlib.pyplot as plt
from Bio import SeqIO

## READ FROM FASTA file and do multiple sequence alignment
# Specify the path to your FASTA file
fasta_file = 'alignment.fasta'

# Read sequences from the file
seq_records = list(SeqIO.parse(fasta_file, 'fasta'))

# Create a Multiple Sequence Alignment (if needed)
alignment = MultipleSeqAlignment(seq_records)

# Calculate genetic distances
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)
print(dm)

# Construct the tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Write the tree distances and structure to a text file
def write_tree_details(tree, distance_matrix, filename):
    with open(filename, 'w') as f:
        # f.write("Distance Matrix:\n")
        # for row in distance_matrix.matrix:
        #     f.write("\t".join(f"{dist:.4f}" for dist in row) + "\n")
        f.write("\nTree Structure:\n")
        Phylo.write(tree, f, 'newick')

write_tree_details(tree, dm, 'tree_details.txt')

# # Draw and save a figure of the tree
# def plot_tree(tree, filename):
#     fig = plt.figure(figsize=(10, 10))
#     ax = fig.add_subplot(1, 1, 1)
#     Phylo.draw(tree, do_show=False, axes=ax)
#     plt.savefig(filename)

# plot_tree(tree, 'tree.png')

##Same function as plot_tree but it prints on the tree the distances of the species and nodes to 3 demcimal places
def plot_tree(tree, filename):
    def format_branch_length(branch):
        if branch.branch_length is not None:
            return f"{branch.branch_length:.3f}"
        return None

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax, branch_labels=format_branch_length)
    plt.savefig(filename)
plot_tree(tree, 'tree.png')