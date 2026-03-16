# HEPTAGON - Heuristic Phylogenetic Tree Algorithm for GPU-Optimized Neighbor-Joining

### HEPTAGON is a software designed to efficiently construct phylogenetic trees for pairwise distance matrices, by performing neighbor joining on NVIDIA GPUs. It takes as input a file containing a distance matrix in Phylip format, and outputs a Newick formatted phylogenetic tree.

## Installation
To compile, simply download and run `make` to use the included makefile. This will produce the executable binary file `./heptagon`.

```
git clone https://github.com/icicl/heptagon.git
cd heptagon/
make
```

## Usage
`./heptagon --in <input_file> [--unflatten] [--skeleton-size n]`
The `--in` argument is required, and must be the path to a distance matrix in phylip format. It may have any extension.
The optional `--unflatten` argument will cause the output to print each taxon to a new line. If omitted, the output will be a single line.
The optional `--skeleton-size` argument specifies the number of taxa *k* that should be used to build the skeleton. If the given file has more than *k* taxa, a skeleton will constructed from *k* of them using neighbor joining, and then the remaining taxa will be added to this skeleton tree using a nearest-neighbor heuristic. Default value is 20000.  Larger values will be more accurate, however you will need a GPU with roughly **(4⋅k²) Bytes** of VRAM.

Debug info is printed to stderr, so you can pipe the output tree into a file to save it.

## Example
`./heptagon --in ./data/sim100k.dist.phylip --skeleton-size 30000 > sim100k.tree`


## Testing

I have included some test files in `test/`. There you may download test files for a larger number of taxa *N*:
```
pip install gdown
gdown --folder https://drive.google.com/drive/folders/1uHfUkZ9JqKQ4S4YB4C-CDasb5oqsm110
```

Extract the test files: 
```
gunzip test/*.gz
```

I compared my results against [Quicktree](https://github.com/khowe/quicktree) and [DIPPER](https://github.com/TurakhiaLab/DIPPER), both of which also implement neighbor joining. See DIPPER's GitHub for how to install it.
You can install Quicktree with the following:
```
git clone https://github.com/khowe/quicktree.git
cd quicktree/
# HEPTAGON - Heuristic Phylogenetic Tree Algorithm for GPU-Optimized Neighbor-Joining

### HEPTAGON is a software designed to efficiently construct phylogenetic trees for pairwise distance matrices, by performing neighbor joining on NVIDIA GPUs. It takes as input a file containing a distance matrix in Phylip format, and outputs a Newick formatted phylogenetic tree.

## Installation
To compile, simply download and run `make` to use the included makefile. This will produce the executable binary file `./heptagon`.

```
git clone https://github.com/icicl/heptagon.git
cd heptagon/
make
```

## Usage
`./heptagon --in <input_file> [--unflatten] [--skeleton-size n]`
The `--in` argument is required, and must be the path to a distance matrix in phylip format. It may have any extension.
The optional `--unflatten` argument will cause the output to print each taxon to a new line. If omitted, the output will be a single line.
The optional `--skeleton-size` argument specifies the number of taxa *k* that should be used to build the skeleton. If the given file has more than *k* taxa, a skeleton will constructed from *k* of them using neighbor joining, and then the remaining taxa will be added to this skeleton tree using a nearest-neighbor heuristic. Default value is 20000.  Larger values will be more accurate, however you will need a GPU with roughly **(4⋅k²) Bytes** of VRAM.

Debug info is printed to stderr, so you can pipe the output tree into a file to save it.

## Example
`./heptagon --in ./data/sim100k.dist.phylip --skeleton-size 30000 > sim100k.tree`


## Testing

I have included some test files in `test/`. There you may download test files for a larger number of taxa *N*:
```
pip install gdown
gdown --folder https://drive.google.com/drive/folders/1uHfUkZ9JqKQ4S4YB4C-CDasb5oqsm110
```

Extract the test files: 
```
gunzip test/*.gz
```

I compared my results against [Quicktree](https://github.com/khowe/quicktree) and [DIPPER](https://github.com/TurakhiaLab/DIPPER), both of which also implement neighbor joining. See DIPPER's GitHub for how to install it.
You can install Quicktree with the following:
```
git clone https://github.com/khowe/quicktree.git
cd quicktree/
make
cd ..
mv quicktree/quicktree qtree
rm -rf quicktree/
mv qtree quicktree
```


To compare HEPTAGON's computed tree against the reference baseline tree, I use normalized Robinson-Foulds distance (nRF), and normalized quartet distance (nQD). I used [MAPLE](https://github.com/NicolaDM/MAPLE) and [tqDist](https://www.birc.au.dk/~cstorm/software/tqdist/), respectively, to compute these metrics.

The linked websites have install instructions, or you can install them as follows:

```
# Download and build tqDist
wget https://www.birc.au.dk/~cstorm/software/tqdist/files/tqDist-1.0.2.tar.gz
tar -xzf tqDist-1.0.2.tar.gz && rm tqDist-1.0.2.tar.gz
cd tqDist-1.0.2/
cmake .
make
make test
cd ..
ln tqDist-1.0.2/bin/quartet_dist 
```
```
# Download MAPLE
wget https://github.com/NicolaDM/MAPLE/blob/main/previousReleases/MAPLEv0.7.5.py -O MAPLE.py
```
Once you have downloaded these packages, run a test for a given number of taxa *N* by running `bash test/run_test.sh N`. For example:
```
bash test/run_test.sh 2000
```
This will run HEPTAGON, and if they are present in the current directory, Quicktree and DIPPER.

Note that the GitHub only has `N=1000,2000` tests, and the Google Drive only has tests for `N=1000,2000,5000,10000,15000,20000,25000`.

