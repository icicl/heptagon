# HEPTAGON - Heuristic Phylogenetic Tree Algorithm for GPU-Optimized Neighbor-Joining

### HEPTAGON is a software designed to efficiently construct phylogenetic trees for pairwise distance matrices, by performing neighbor joining on NVIDIA GPUs. It takes as input a file containing a distance matrix in Phylip format, and outputs a Newick formatted phylogenetic tree.

## Installation
To compile, simply run `make` to use the included makefile. This will produce the executable binary file `./heptagon`.

## Usage
`./heptagon --in <input_file> [--unflatten] [--skeleton-size n]`
The `--in` argument is required, and must be the path to a distance matrix in phylip format. It may have any extension.
The optional `--unflatten` argument will cause the output to print each taxon to a new line. If omitted, the output will be a single line.
The optional `--skeleton-size` argument specifies the number of taxa *k* that should be used to build the skeleton. If the given file has more than *k* taxa, a skeleton will constructed from *k* of them using neighbor joining, and then the remaining taxa will be added to this skeleton tree using a nearest-neighbor heuristic. Default value is 20000.  Larger values will be more accurate, however you will need a GPU with roughly **(4⋅k²) Bytes** of VRAM.

Debug info is printed to stderr, so you can pipe the output tree into a file to save it.

## Example
`./heptagon --in ./data/sim100k.dist.phylip --skeleton-size 30000 > sim100k.tree`