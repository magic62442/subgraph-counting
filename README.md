# SCOPE

## Background

We study local subgraph counting queries, Q = (p, o), to count how many times a given k-node pattern graph p appears around every node v in a data graph G when the given node orbit o in p maps to v. 

## Compile

1. Compile and link to the [nauty](https://pallini.di.uniroma1.it) library.  The nauty library is used to compute automorphisms and symmetry-breaking rules. We include a copy of the nauty library in /utility/automorphisms and show the steps.

```shell
cd utility/automorphism/
./configure
make
mv nauty.a libnauty.a
```

If it complains, "relocation R_X86_64_32 against `.rodata.str1.1' can not be used when making a shared object; recompile with the "-fPIC" option. 

```shell
cd utility/automorphism/
vim makefile
# add -fPIC to the end of line 6.
make
mv nauty.a libnauty.a
```

2. Compile and link to the [GLPK](https://www.gnu.org/software/glpk/) library. The GLPK library is used to compute fractional edge covers. Edit the paths in CMakeLists.txt accordingly.

3. Build the project.

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Input format

The data graph should start with 'n, m' where n is the number of nodes and m is the number of undirected edges, followed by the edge list. The node id should be consecutive and should start from 0.

Example:

```
3 2
0 1
1 2
```

The query graph file has an additional line, '1 id', where 'id' is this pattern's orbit(representative).

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
| -share | with -share: enable sharing, without -share: disable sharing |

Example:

```
./build/executable/scope.out -q ./exp/pattern_graph/q1 -d ./exp/data_graph/web-spam.txt -t ./exp/data_graph/web-spamt.bin -r ./result/q1/web-spam -b -share
./build/executable/scope.out -q ./exp/pattern_graph/q1/34.txt -d ./exp/data_graph/web-spam.txt -r ./result/q1/web-spam/34.txt -share
```

In the output, the $i$-th line shows the local subgraph count of the data node $i-1$.

