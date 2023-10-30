//
// Created by Qiyan LI on 2022/8/30.
//

#include "command.h"
#include "execution.h"

int main(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string queryGraphPath = cmd.getQueryGraphPath();
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    std::string trianglePath = cmd.getTrianglePath();
    bool batchQuery = cmd.getBatchQuery();
    bool singleAggregation = cmd.getSingleAggregation();
    bool shareNode = cmd.getShareNode();
    bool useTriangle = !trianglePath.empty();
    std::cout << "query graph path: " << queryGraphPath << std::endl;
    std::cout << "data graph path: " << dataGraphPath << std::endl;
    std::cout << "result path: " << resultPath << std::endl;
    std::cout << "using batch query: " << batchQuery << std::endl;
    std::cout << "using multi aggregation: " << !singleAggregation << std::endl;
    std::cout << "sharing nodes computation: " << shareNode << std::endl;
    std::cout << "using triangle: " << useTriangle << std::endl;
    std::cout << "set intersection type: " << SI << std::endl;
    DataGraph dun = DataGraph();
    dun.loadDataGraph(dataGraphPath);
    const DataGraph din = constructDirectedDataGraph(dun, false);
    const DataGraph dout = constructDirectedDataGraph(dun, true);
    specialsparse *sg = (specialsparse *)malloc(sizeof(specialsparse));
    dout.initSpecialSparse(sg);
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    Triangle triangle;
    if (useTriangle) {
        triangle.load(trianglePath, m / 2);
        std::cout << "triangle memory cost: " << triangle.memoryCost() << std::endl;
    }
    auto start = std::chrono::steady_clock::now();
    EdgeID *outID = buildInID2OutID(din, dout);
    EdgeID *unID = buildUnID2OutID(dun, dout);
    EdgeID *reverseID = buildReverseUnID(dun);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::cout << "finished building edge ID mappings. time: " << elapsedSeconds.count() << "s" << std::endl;
    std::vector<std::string> files;
    std::vector<PatternGraph> patternGraphs = loadPatternGraph(queryGraphPath, batchQuery, files);
    int patternSize = patternGraphs[0].getNumVertices();
    bool forestShare = patternSize <= 5;
    HashTable factorSum = new Count[m + 1];
    HashTable ht[MAX_NUM_NODE];
    for (auto & h: ht)
        h = new Count[m + 1];
    ui *candPos = new ui[MAX_PATTERN_SIZE];
    memset(candPos, 0, sizeof(ui) * (MAX_PATTERN_SIZE));
    VertexID **patternV = new VertexID *[MAX_NUM_NODE];
    VertexID **dataV = new VertexID *[MAX_NUM_NODE];
    EdgeID **startOffset = new EdgeID *[MAX_NUM_NODE];
    bool **visited = new bool *[MAX_NUM_NODE];
    for (int i = 0; i < MAX_NUM_NODE; ++i) {
        patternV[i] = new VertexID[MAX_PATTERN_SIZE];
        dataV[i] = new VertexID[MAX_PATTERN_SIZE];
        startOffset[i] = new VertexID[MAX_PATTERN_SIZE];
        visited[i] = new bool[n];
        memset(visited[i], false, sizeof(bool) * n);
    }
    VertexID *tmp = new VertexID[n];
    VertexID *allV = new VertexID[n];
    for (VertexID i = 0; i < n; ++i) {
        allV[i] = i;
    }
    VertexID *cliqueVertices = new VertexID[MAX_PATTERN_SIZE + 1];
    std::set<CanonType> visitedNode;

    double totalPlanTime = 0.0, totalExeTime = 0.0;
    ui totalNumPatterns = 0, totalNodes = 0;
    double averageNodeSize = 0.0;
    int numVertexTable = 0, numEdgeTable = 0;
    if (!shareNode) {
        std::vector<HashTable> mathCalH(patternGraphs.size());
        std::vector<int> orbitTypes(patternGraphs.size());
#ifdef DEBUG
        HashTable hPattern = new Count[n];
        memset(hPattern, 0, sizeof(Count) * dun.getNumVertices());
#endif
        for (int i = 0; i < patternGraphs.size(); ++i) {
            double exeTime = 0.0;
            std::map<int, std::vector<Pattern>> patterns;
            std::map<int, std::vector<std::vector<Tree>>> trees;
            start = std::chrono::steady_clock::now();
            ConNode cn;
            if (!patternGraphs[i].isClique()) {
                genEquation(patternGraphs[i], patterns, trees, cn, useTriangle, true, true);
            }
            end = std::chrono::steady_clock::now();
            elapsedSeconds = end - start;
            totalPlanTime += elapsedSeconds.count();
            start = std::chrono::steady_clock::now();
            HashTable H;
            int orbitType = patternGraphs[i].getOrbitType();
            orbitTypes[i] = orbitType;
            if (orbitType == 0) {
                H = new Count[1];
                H[0] = 0;
            }
            else if (orbitType == 1) {
                H = new Count[dun.getNumVertices()];
                memset(H, 0, sizeof(Count) * n);
            }
            else {
                H = new Count[m];
                memset(H, 0, sizeof(Count) * m);
            }
            if (cn.num != 0) {
                Pattern p(patternGraphs[i]);
                // prepare hash tables for cn
                if (orbitType == 1) {
                    for (int j = 0; j < cn.aggrePos.size(); ++j)
                        cn.hashTables[j] = H;
                }
                else {
                    for (int j = 0; j < cn.aggrePos.size() / 2; ++j)
                        cn.hashTables[j] = H;
                }
                executeConNode(cn, din, dout, dun, useTriangle, triangle, p, outID, unID, startOffset[0], patternV[0],
                               dataV[0], visited[0], candPos, tmp, allV);
                if (!resultPath.empty()) {
                    int divideFactor = cn.divideFactor;
                    if (orbitType == 0) H[0] /= divideFactor;
                    else if (orbitType == 1)
                        for (VertexID l = 0; l < n; ++l)
                            H[l] /= divideFactor;
                    else {
                        for (EdgeID l = 0; l < m ; ++l)
                            H[l] /= divideFactor;
                    }
                }
            }
            if (patternGraphs[i].isClique()) {
                start = std::chrono::steady_clock::now();
                int k = patternGraphs[i].getNumVertices();
                mkspecial(sg, k);
                kclique(k, k, sg, cliqueVertices, H, orbitType);
                freesub(sg, k);
                end = std::chrono::steady_clock::now();
                elapsedSeconds = end - start;
                totalExeTime += elapsedSeconds.count();
            }
            for (auto it = patterns.begin(); it != patterns.end(); ++it) {
                int divideFactor = it->first;
                memset(factorSum, 0, sizeof(Count) * (m + 1));
                for (int j = 0; j < it->second.size(); ++j) {
                    for (int j2 = 0; j2 < trees[divideFactor][j].size(); ++j2) {
                        for (int l = 0; l < trees[divideFactor][j][0].getNumNodes(); ++l) {
                            memset(ht[l], 0, sizeof(Count) * m);
                        }
                        int k = patterns[divideFactor][j].u.getNumVertices();
                        if (patterns[divideFactor][j].u.isClique() && k >= 4) {
                            HashTable h = ht[trees[divideFactor][j][j2].getRootID()];
                            int aggreWeight = trees[divideFactor][j][j2].getAggreWeight()[0];
                            mkspecial(sg, k);
                            kclique(k, k, sg, cliqueVertices, h, orbitType);
                            freesub(sg, k);
                            if (aggreWeight != 1) {
                                if (orbitType == 0) h[0] *= aggreWeight;
                                else if (orbitType == 1) {
                                    for (VertexID v = 0; v < n; ++v)
                                        h[v] *= aggreWeight;
                                }
                                else {
                                    for (EdgeID e = 0; e < m; ++e)
                                        h[e] *= aggreWeight;
                                }
                            }
                        }
                        else {
                            const Tree &t = trees[divideFactor][j][j2];
                            if (t.getExecuteMode()) {
                                executeTree(t, din, dout, dun, useTriangle, triangle, patterns[divideFactor][j],
                                            ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
                            }
                            else {
                                multiJoinTree(t, din, dout, dun, useTriangle, triangle, patterns[divideFactor][j],
                                              ht, outID, unID, reverseID, startOffset, patternV, dataV, visited, tmp, allV);
                            }
                        }
                        HashTable h = ht[trees[divideFactor][j][j2].getRootID()];
                        int multiFactor = trees[divideFactor][j][j2].getMultiFactor();
                        if (!resultPath.empty()) {
                            if (orbitType == 0) factorSum[0] += h[0];
                            else if (orbitType == 1)
                                for (VertexID l = 0; l < n; ++l) {
                                    factorSum[l] += h[l] * multiFactor;
                                }
                            else
                                for (EdgeID l = 0; l < m + 1; ++l) {
                                    factorSum[l] += h[l] * multiFactor;
                                }
                        }
#ifdef DEBUG
                        if (divideFactor == 2 && patterns[divideFactor][j].u.getCanonValue() == 281875) {
                        for (VertexID l = 0; l < n; ++l) {
                            hPattern[l] += h[l];
                        }
                    }
#endif
                        for (VertexID nID = 0; nID < trees[divideFactor][j][j2].getNumNodes(); ++nID) {
                            ++totalNodes;
                            averageNodeSize += trees[divideFactor][j][j2].getNode(nID).numVertices;
                            visitedNode.insert(trees[divideFactor][j][j2].getNode(nID).canonValue);
                            if (nID != trees[divideFactor][j][j2].getPostOrder().back()) {
                                if (trees[divideFactor][j][j2].getNode(nID).keySize == 1)
                                    ++numVertexTable;
                                else
                                    ++numEdgeTable;
                            }
                        }
                    }
                }
                if (!resultPath.empty()) {
                    if (orbitType == 0) H[0] += factorSum[0] / divideFactor;
                    else if (orbitType == 1)
                        for (VertexID l = 0; l < n; ++l) {
                            H[l] += factorSum[l] / divideFactor;
                        }
                    else
                        for (EdgeID l = 0; l < m + 1; ++l) {
                            H[l] += factorSum[l] / divideFactor;
                        }
                }
            }
            end = std::chrono::steady_clock::now();
            elapsedSeconds = end - start;
            exeTime = elapsedSeconds.count();
            if (batchQuery)
                std::cout << "file: " << files[i] << ", ";
            std::cout << "execution time: " << exeTime << "s";
            ui numPatterns = 0;
            for (auto it = patterns.begin(); it != patterns.end(); ++it)
                numPatterns += it->second.size();
            std::cout << ", number of patterns: " << numPatterns << std::endl;
            totalExeTime += exeTime;
            totalNumPatterns += numPatterns;
            mathCalH[i] = H;
        }
        if (!resultPath.empty()) saveCount(resultPath, mathCalH, dun, batchQuery, files, orbitTypes);
        for (int i = 0; i < patternGraphs.size(); ++i)
            delete[] mathCalH[i];
        std::cout << "compute shrinakge time: " << gShrinkTime << ", " << gShrinkTime / totalPlanTime << std::endl;
        std::cout << "planning time: " << totalPlanTime << ", total execution time: "
                  << totalExeTime << ", total number of patterns: " << totalNumPatterns << "\ntotal number of nodes: "
                  << totalNodes << ", average size of nodes: "<< averageNodeSize / totalNodes
                  << ", number of distinct nodes: " << visitedNode.size() << std::endl
                  << "number of vertex hash tables (except for root node): " << numVertexTable
                  << ", number of edge hash tables (except for root node): " << numEdgeTable << std::endl;
        std::cout << "number of match: " << gNumMatch << ", number of intersect: " << gNumIntersect
                  << ", number of edge ID: " << gNumEdgeID << ", number of update: " << gNumUpdate << std::endl;
    }

    if (shareNode) {
        std::vector<ConNode> conNodes;
        std::vector<Pattern> conPatterns;
        std::vector<int> conFactors;
        std::vector<int> conIDs;
        Forest fu, fuShrink;
        std::vector<Forest> directedF;
        std::vector<CanonType> shrinkCanon;
        std::vector<int> shrinkOrbit;
        std::vector<int> shrinkRootOrbit;
        std::vector<HashTable> result(patternGraphs.size());
        std::vector<HashTable> shrinkageH;
        std::vector<HashTable> uH;
        std::vector<std::vector<HashTable>> dH;
        // for the ith Pattern, jth shrinkage, its position in the shrinkageH and the factor
        std::vector<std::vector<int>> shrinkPos(patternGraphs.size());
        // for the ith pattern, position j, the multiply factor of that shrinkage
        std::vector<std::vector<int>> shrinkMultiFactors(patternGraphs.size());
        std::vector<int> divideFactors(patternGraphs.size());
        std::vector<bool> isUndirected(patternGraphs.size(), false);
        std::vector<HashTable> allocatedHashTable;
        allocatedHashTable.reserve(10000);
        int nTable = 0, mTable = 0;
        for (int i = 0; i < patternGraphs.size(); ++i) {
            const PatternGraph &pg = patternGraphs[i];
            int orbitType = pg.getOrbitType();
            if (orbitType == 0) {
                result[i] = new Count[1];
                result[i][0] = 0;
            }
            else if (orbitType == 1) {
                result[i] = new Count[n];
                memset(result[i], 0, sizeof(Count) * n);
                ++nTable;
            }
            else {
                if (pg.isEOrbitDir()) {
                    result[i] = new Count[m / 2];
                    memset(result[i], 0, sizeof(Count) * m / 2);
                    ++mTable;
                }
                else {
                    result[i] = new Count[m];
                    memset(result[i], 0, sizeof(Count) * m);
                    mTable += 2;
                }
            }
        }
        for (int i = 0; i < patternGraphs.size(); ++i) {
            const PatternGraph &pg = patternGraphs[i];
            int orbitType = pg.getOrbitType();
            start = std::chrono::steady_clock::now();
            if (patternGraphs[i].isClique()) {
                start = std::chrono::steady_clock::now();
                int k = patternGraphs[i].getNumVertices();
                mkspecial(sg, k);
                kclique(k, k, sg, cliqueVertices, result[i], orbitType);
                freesub(sg, k);
                end = std::chrono::steady_clock::now();
                elapsedSeconds = end - start;
                totalExeTime += elapsedSeconds.count();
            }
            else {
                std::map<int, std::vector<Pattern>> patterns;
                std::map<int, std::vector<std::vector<Tree>>> trees;
                ConNode cn;
                bool directed = genEquation(pg, patterns, trees, cn, useTriangle, true, true);
                if (cn.num != 0) {
                    if (orbitType == 0) cn.hashTables[0] = result[i];
                    else if (orbitType == 1) {
                        for (int j = 0; j < cn.aggrePos.size(); ++j)
                            cn.hashTables[j] = result[i];
                    }
                    else {
                        for (int j = 0; j < cn.aggrePos.size() / 2; ++j)
                            cn.hashTables[j] = result[i];
                    }
                    bool exists = false;
                    for (int j = 0; j < conNodes.size(); ++j) {
                        if (conNodes[j].num == cn.num && conNodes[j].canonValue == cn.canonValue && conNodes[j].edgeKey == cn.edgeKey) {
                            conNodes[j].merge(cn, conPatterns[j].u, patternGraphs[i]);
                            exists = true;
                            break;
                        }
                    }
                    if (!exists) {
                        conNodes.push_back(cn);
                        conPatterns.emplace_back(pg);
                    }
                    conIDs.push_back(i);
                    conFactors.push_back(cn.divideFactor);
                }
                else if (directed) {
                    std::cout << files[i] << " uses DAG" << std::endl;
                    if (pg.getNumVertices() >= 6 || directedF.empty()) {
                        Forest fd;
                        std::vector<HashTable> h;
                        h.push_back(result[i]);
                        dH.push_back(h);
                        fd.loadQuery(patterns, trees, true);
                        directedF.push_back(fd);
                    }
                    else {
                        directedF[0].loadQuery(patterns, trees, true);
                        dH[0].push_back(result[i]);
                    }
                }
                else {
                    isUndirected[i] = true;
                    // put the first pattern to fu, the remaining patterns to fShrinkage
                    int divideFactor = patterns.begin() -> first;
                    divideFactors[i] = divideFactor;
                    const std::vector<Pattern> &allPattern = patterns.begin()->second;
                    std::vector<std::vector<Tree>> &allTree = trees.begin()->second;
                    std::map<int, std::vector<Pattern>> rootPattern;
                    std::map<int, std::vector<std::vector<Tree>>> rootTree;
                    rootPattern[divideFactor].push_back(allPattern[0]);
                    rootTree[divideFactor].push_back(allTree[0]);
                    shrinkPos[i] = std::vector<int>(allPattern.size() - 1);
                    shrinkMultiFactors[i] = std::vector<int>(allPattern.size() - 1 + shrinkCanon.size());
                    for (int j = 1; j < allPattern.size(); ++j) {
                        bool exists = false;
                        int rootOrbit = allTree[j][0].getNode(allTree[j][0].getRootID()).v2o[0];
                        for (int k = 0; k < shrinkCanon.size(); ++k) {
                            if (allPattern[j].u.getCanonValue() == shrinkCanon[k] && allPattern[j].u.getOrbit(0) == shrinkOrbit[k]
                               && rootOrbit == shrinkRootOrbit[k]) {
                                shrinkPos[i][j - 1] = k;
                                shrinkMultiFactors[i][k] = allTree[j][0].getMultiFactor();
                                exists = true;
                                break;
                            }
                        }
                        if (!exists) {
                            HashTable h;
                            if (orbitType == 0) {
                                h = new Count[1];
                                h[0] = 0;
                            }
                            else if (orbitType == 1) {
                                h = new Count[n];
                                memset(h, 0, sizeof(Count) * n);
                                ++nTable;
                            }
                            else {
                                if (pg.isEOrbitDir()) {
                                    h = new Count[m / 2];
                                    memset(h, 0, sizeof(Count) * m / 2);
                                    ++mTable;
                                }
                                else {
                                    h = new Count[m];
                                    memset(h, 0, sizeof(Count) * m);
                                    mTable += 2;
                                }
                            }
                            shrinkPos[i][j - 1] = shrinkCanon.size();
                            shrinkMultiFactors[i][shrinkCanon.size()] = allTree[j][0].getMultiFactor();
                            allTree[j][0].setMultiFactor(1);
                            shrinkCanon.push_back(allPattern[j].u.getCanonValue());
                            shrinkOrbit.push_back(allPattern[j].u.getOrbit(0));
                            shrinkRootOrbit.push_back(rootOrbit);
                            std::map<int, std::vector<Pattern>> shrinkPattern;
                            std::map<int, std::vector<std::vector<Tree>>> shrinkTree;
                            shrinkPattern[1].push_back(allPattern[j]);
                            shrinkTree[1].push_back(allTree[j]);
                            fuShrink.loadQuery(shrinkPattern, shrinkTree, forestShare);
                            shrinkageH.push_back(h);
                        }
                    }
                    uH.push_back(result[i]);
                    fu.loadQuery(rootPattern, rootTree, forestShare);
                }
            }
            end = std::chrono::steady_clock::now();
            elapsedSeconds = end - start;
            totalPlanTime += elapsedSeconds.count();
        }
        double tableSize = (double)(nTable * n + mTable * m / 2) * sizeof(Count) / 1e9;
//        std::cout << "allocated memory for hash tables in the main function: size = " <<  tableSize << "GB" << std::endl;
//        std::cout << "number of vertex tables: " << nTable << ", number of edge tables: " << mTable << std::endl;
#ifndef ONLY_PLAN
        // 1. execute ConNodes
        start = std::chrono::steady_clock::now();
        for (int i = 0; i < conNodes.size(); ++i) {
            executeConNode(conNodes[i], din, dout, dun, useTriangle, triangle, conPatterns[i], outID, unID,
                           startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
        }
        for (int i = 0; i < conIDs.size(); ++i) {
            int divideFactor = conFactors[i];
            int orbitType = patternGraphs[i].getOrbitType();
            bool eOrbitDir = patternGraphs[i].isEOrbitDir();
            int id = conIDs[i];
            if (orbitType == 0) result[id][0] /= divideFactor;
            else if (orbitType == 1) {
                for (int k = 0; k < n; ++k)
                    result[id][k] /= divideFactor;
            }
            else {
                if (eOrbitDir) {
                    for (int k = 0; k < m / 2; ++k)
                        result[id][k] /= divideFactor;
                }
                else {
                    for (int k = 0; k < m; ++k)
                        result[id][k] /= divideFactor;
                }
            }
        }
        end = std::chrono::steady_clock::now();
        elapsedSeconds = end - start;
        totalExeTime += elapsedSeconds.count();
        std::cout << "finished compressed node. time: " << totalExeTime << "s" << std::endl;
        // 2. execute undirected shrinkCanon, undirected root patterns and directed patterns
        if (patternSize <= 5) {
            start = std::chrono::steady_clock::now();
            executeForest(fuShrink, shrinkageH, din, dout, dun, useTriangle, triangle, allocatedHashTable,
                          outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV, sg, cliqueVertices);
            for (HashTable h: allocatedHashTable) delete[] h;
            allocatedHashTable.clear();
            end = std::chrono::steady_clock::now();
            elapsedSeconds = end - start;
            totalExeTime += elapsedSeconds.count();
            std::cout << "finished undirected shrinkages. time: " << totalExeTime << std::endl;
            start = std::chrono::steady_clock::now();
            executeForest(fu, uH,  din, dout, dun, useTriangle, triangle, allocatedHashTable,
                          outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV, sg, cliqueVertices);
            for (HashTable h: allocatedHashTable) delete[] h;
            allocatedHashTable.clear();
            end = std::chrono::steady_clock::now();
            elapsedSeconds = end - start;
            totalExeTime += elapsedSeconds.count();
            std::cout << "finished undirected patterns. time: " << totalExeTime << std::endl;
        }
        else {
            start = std::chrono::steady_clock::now();
            for (int i = 0; i < fuShrink.numQuery; ++i) {
                const Tree &t = fuShrink.allTree[i];
                const Pattern &p = fuShrink.allPattern[i];
                for (int l = 0; l < t.getNumNodes(); ++l) memset(ht[l], 0, sizeof(Count) * m);
                if (t.getExecuteMode()) {
                    executeTree(t, din, dout, dun, useTriangle, triangle, p, ht, outID, unID, reverseID,
                                startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
                }
                else {
                    multiJoinTree(t, din, dout, dun, useTriangle, triangle, p, ht, outID, unID, reverseID, startOffset,
                                  patternV, dataV, visited, tmp, allV);
                }
                HashTable h = ht[t.getRootID()];
                int orbitType = patternGraphs[i].getOrbitType();
                if (orbitType == 0) {
                    shrinkageH[i][0] = h[0];
                }
                else if (orbitType == 1) {
                    for (int j = 0; j < n; ++j) {
                        shrinkageH[i][j] = h[j];
                    }
                }
            }
            end = std::chrono::steady_clock::now();
            elapsedSeconds = end - start;
            totalExeTime += elapsedSeconds.count();
            std::cout << "finished undirected shrinkages. time: " << totalExeTime << std::endl;
            start = std::chrono::steady_clock::now();
            for (int i = 0; i < fu.numQuery; ++i) {
                const Tree &t = fu.allTree[i];
                const Pattern &p = fu.allPattern[i];
                for (int l = 0; l < t.getNumNodes(); ++l) memset(ht[l], 0, sizeof(Count) * m);
                if (t.getExecuteMode()) {
                    executeTree(t, din, dout, dun, useTriangle, triangle, p, ht, outID, unID, reverseID,
                                startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
                }
                else {
                    multiJoinTree(t, din, dout, dun, useTriangle, triangle, p, ht, outID, unID, reverseID, startOffset,
                                  patternV, dataV, visited, tmp, allV);
                }
                HashTable h = ht[t.getRootID()];
                int multiFactor = t.getMultiFactor();
                int orbitType = patternGraphs[i].getOrbitType();
                if (orbitType == 0) {
                    uH[i][0] = h[0] * multiFactor;
                }
                else if (orbitType == 1) {
                    for (int j = 0; j < n; ++j) {
                        uH[i][j] = h[j] * multiFactor;
                    }
                }
            }
            end = std::chrono::steady_clock::now();
            elapsedSeconds = end - start;
            totalExeTime += elapsedSeconds.count();
            std::cout << "finished undirected patterns. time: " << totalExeTime << std::endl;
        }
        start = std::chrono::steady_clock::now();
        for (int i = 0; i < directedF.size(); ++i) {
            executeForest(directedF[i], dH[i], din, dout, dun, useTriangle, triangle, allocatedHashTable,
                          outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV, sg, cliqueVertices);
            for (HashTable h: allocatedHashTable) delete[] h;
            allocatedHashTable.clear();
        }
        end = std::chrono::steady_clock::now();
        elapsedSeconds = end - start;
        totalExeTime += elapsedSeconds.count();
        std::cout << "finished directed patterns. time : " << totalExeTime << std::endl;
        // 3. subtract shrinkage count for undirected patterns
        for (int i = 0; i < patternGraphs.size(); ++i) {
            if (!isUndirected[i]) continue;
            const PatternGraph &pg = patternGraphs[i];
            int orbitType = pg.getOrbitType();
            int divideFactor = divideFactors[i];
            for (int j = 0; j < shrinkPos[i].size(); ++j) {
                HashTable h = shrinkageH[shrinkPos[i][j]];
                if (orbitType == 0)
                    result[i][0] += h[0] * shrinkMultiFactors[i][shrinkPos[i][j]];
                else if (orbitType == 1) {
                    for (int k = 0; k < n; ++k)
                        result[i][k] += h[k] * shrinkMultiFactors[i][shrinkPos[i][j]];
                }
                else if (orbitType == 2 && pg.isEOrbitDir()) {
                    for (int k = 0; k < m / 2; ++k)
                        result[i][k] += h[k] * shrinkMultiFactors[i][shrinkPos[i][j]];
                }
                else {
                    for (int k = 0; k < m; ++k)
                        result[i][k] += h[k] * shrinkMultiFactors[i][shrinkPos[i][j]];
                }
            }
            if (orbitType == 0)
                result[i][0] /= divideFactor;
            else if (orbitType == 1) {
                for (int k = 0; k < n; ++k)
                    result[i][k] /= divideFactor;
            }
            else if (orbitType == 2 && pg.isEOrbitDir()) {
                for (int k = 0; k < m / 2; ++k)
                    result[i][k] /= divideFactor;
            }
            else {
                for (int k = 0; k < m; ++k)
                    result[i][k] /= divideFactor;
            }
        }
#endif
        std::cout << "compute shrinakge time: " << gShrinkTime << ", " << gShrinkTime / totalPlanTime << std::endl;
        std::cout << "total planning time: " << totalPlanTime << ", total execution time: " << totalExeTime
                  << ", total time: " << totalPlanTime + totalExeTime << std::endl;
        std::cout << "number of match: " << gNumMatch << ", number of intersect: " << gNumIntersect
                  << ", number of edge ID: " << gNumEdgeID << ", number of update: " << gNumUpdate << std::endl;
        int orbitType = patternGraphs[0].getOrbitType();
        std::vector<int> orbitTypes(files.size(), orbitType);
        if (!resultPath.empty()) saveCount(resultPath, result, dun, batchQuery, files, orbitTypes);
    }

    return 0;
}