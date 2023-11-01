# SCOPE

## Background

We study local subgraph counting queries, Q = (p, o), to count how many times a given k-node pattern graph p appears around every node v in a data graph G when the given center node o in p maps to v. 

## Compile

1. Compile and link to the [nauty](https://pallini.di.uniroma1.it) library.  The nauty library is used to compute automorphisms and symmetry-breaking rules. We include a copy of the nauty library in /utility/automorphisms and show the steps.

```shell
cd utility/automorphism/
./configure
make -j
mv nauty.a libnauty.a
```

2. Compile and link to the [GLPK](https://www.gnu.org/software/glpk/) library. The GLPK library is used to compute fractional edge covers. Edit the paths in CMakeLists.txt accordingly.

3. Build the project.

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

During the compilation, if it reports the error "relocation R_X86_64_32 against `.rodata.str1.8' can not be used when making a shared object; recompile with -fPIC", you should edit the makefile of nauty and recompile. 

```shell
cd utility/automorphism/
sed -i '6s/$/ -fPIC/' makefile
make clean
make -j
mv nauty.a libnauty.a
```

## Input format

The data graph should start with 'n, m' where n is the number of nodes and m is the number of undirected edges, followed by the edge list. The node id should be consecutive and should start from 0.

Example:

```
3 2
0 1
1 2
```

The query graph file has an additional line, '1 id', where 'id' is this pattern's orbit(center node).

Example:

```
3 2
0 1
1 2
1 0
```

The queries are in the ./exp/pattern_graph directory, and the data graphs can be downloaded from [SNAP](https://snap.stanford.edu/data/index.html) or the [Network Repository](https://networkrepository.com).

## Preprocessing (optional)

Usage:

```shell
./preprocess.out <data graph input path> <reordered graph output path> <intersection cache output path>
```

## Execution and output

Command line options:

| Option | Description                                                  |
| ------ | ------------------------------------------------------------ |
| -q     | the query graph path (single query) or directory (batch query) |
| -d     | the data graph path                                          |
| -t     | the intersection cache path, optional                        |
| -r     | the result path (single query) or directory (batch query), optional |
| -b     | with -b: batch query, without -b: single query               |
| -share | with -share: enable sharing (only for 5-node queries), without -share: disable sharing |

Example:

```
./build/executable/scope.out -q ./exp/pattern_graph/5node -d ./exp/data_graph/web-spam.txt -t ./exp/data_graph/web-spamt.bin -r ./result/5node/web-spam -b -share
./build/executable/scope.out -q ./exp/pattern_graph/5node/34.txt -d ./exp/data_graph/web-spam.txt -r ./result/5node/web-spam/34.txt
```

In the output, the $i$-th line shows the local subgraph count of the data node $i-1$.

