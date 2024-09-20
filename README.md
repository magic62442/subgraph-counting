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

The data graph should start with 'n, m' where n is the number of nodes and m is the number of undirected edges, followed by the edge list. The node id should be consecutive and should start from 0. Each edge only appears once, in the format of 'smaller_vertex_id larger_vertex_id'. The separator is a ' '.

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

## Preprocessing

Usage:

```shell
./preprocess.out <data graph input path> <reordered graph output path> <intersection cache output path>
```

We preprocessed the data graphs in our experiments by indexing triangles as a simple intersection cache. It is static, and there is no swap. This can make SCOPE faster. Note that EVOKE also index triangles., and DISC has a more advanced intersection cache with swapping.

## Execution and output

### scope.out:

| Option | Description                                                  |
| ------ | ------------------------------------------------------------ |
| -q     | the query graph path (single query) or directory (batch query) |
| -d     | the data graph path                                          |
| -t     | the intersection cache path, optional                        |
| -r     | the result path (single query) or directory (batch query), optional |
| -b     | with -b: batch query, without -b: single query               |
| -share | with -share: enable sharing the computation of Cov pattern counts if there are multiple queries. without -share: disable sharing |

Example:

```
./build/executable/scope.out -q ./exp/pattern_graph/5voc -d ./exp/data_graph/web-spam.txt -t ./exp/data_graph/web-spamt.bin -r ./result/5voc/web-spam -b -share
./build/executable/scope.out -q ./exp/pattern_graph/5voc/62.txt -d ./exp/data_graph/web-spam.txt -r ./result/5voc/web-spam/62.txt
```

In the output, the $i$-th line shows the local subgraph count of the data node $i-1$.

### batch.out

This version uses the precomputed plan for queries. We use it for GNN datasets, where there are thousands of relatively small data graphs to count. You need to specify the query graph path and intersection cache path. We provided generated plans in ./exp/plan for all 5-node and 6-node queries. To precompute plans for other query sets, you need to modify the function 'generatePlan' in 'batch.cpp'.

### 5voc.out

We further optimized the code for orbit counting for 5-vertex queries by hand, following the generated plans. This is not included in the paper. It is about 2-3 times faster than EVOKE. You need to specify the data graph path and intersection cache path.
