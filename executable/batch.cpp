#include "command.h"
#include "execution.h"

void generatePlan() {
    std::string patternDir = "../../exp/pattern_graph";
    std::string planDir = "../../exp/plan";
    std::vector<std::string> queries = {"", "", "2voc", "3voc", "4voc", "5voc", "6voc"};
    std::vector<std::vector<PatternGraph>> patternGraphs;
    std::vector<std::vector<std::string>> files;
    patternGraphs.resize(8);
    files.resize(8);
    for (int i = 2; i < 8; ++i) {
        std::string patternSubDir = pathJoin(patternDir, queries[i]);
        std::string planSubDir = pathJoin(planDir, queries[i]);
        patternGraphs[i] = loadPatternGraph(patternSubDir, true, files[i]);
        for (std::string& str : files[i]) {
            str.erase(str.length() - 4);
        }
        for (int j = 0; j < patternGraphs[i].size(); ++j) {
            std::map<int, std::vector<Pattern>> patterns;
            std::map<int, std::vector<std::vector<Tree>>> trees;
            ConNode cn;
            std::string planPath = pathJoin(planSubDir, files[i][j] + "-tree");
            std::string equationPath = pathJoin(planSubDir, files[i][j] + "-equation");
            std::string conNodePath = pathJoin(planSubDir, files[i][j] + "-connode");
            genEquation(patternGraphs[i][j], patterns, trees, cn, true, true, true, false);
            int divideFactor = patterns.begin() -> first;
            if (cn.num != 0) {
                std::ofstream outFile(conNodePath);
                cn.writeToStream(outFile);
            }
            else {
                // write the equation
                std::ofstream outFile(equationPath);
                outFile.write(reinterpret_cast<const char*>(&divideFactor), sizeof(divideFactor));
                std::vector<ui> patternSizes;
                std::vector<std::string> patternIDs;
                std::vector<int> multiFactors;
                for (int k = 1; k < patterns[divideFactor].size(); ++k) {
                    const PatternGraph &p = patterns[divideFactor][k].u;
                    ui patternSize = p.getNumVertices();
                    patternSizes.push_back(patternSize);
                    multiFactors.push_back(trees[divideFactor][k][0].getMultiFactor());
                    for (int l = 0; l < patternGraphs[patternSize].size(); ++l) {
                        if (patternGraphs[patternSize][l].getCanonValue() == p.getCanonValue()
                        && patternGraphs[patternSize][l].getOrbit(0) == p.getOrbit(0))
                            patternIDs.push_back(files[patternSize][l]);
                    }
                }
                writeVectorToStream(outFile, patternSizes);
                writeVectorToStream(outFile, patternIDs);
                writeVectorToStream(outFile, multiFactors);
                outFile.close();
            }
            // write the td
            trees[divideFactor][0][0].writeToFile(planPath);
        }
        std::cout << "finish k = " << i << std::endl;
    }
}

void readPlan(int k, const std::vector<PatternGraph> &patternGraphs, const std::vector<std::string> &files,
              std::vector<ConNode> &conNodes, std::vector<Pattern> &conPatterns, std::vector<int> &conFactors,
              std::vector<int> &conIDs, std::vector<std::vector<int>> &coverPos, std::vector<Tree> &coverTrees,
              std::vector<Pattern> &coverPatterns, std::vector<std::vector<int>> &coverMultiFactors,
              const std::vector<HashTable> &result,
              std::vector<Tree> &trees, std::vector<Pattern> &patterns,std::vector<int> &divideFactors) {
    std::string patternDir = "../../exp/pattern_graph";
    std::string planDir = "../../exp/plan";
    std::vector<std::string> queries = {"", "", "2voc", "3voc", "4voc", "5voc", "6voc", "new7voc"};
    std::vector<ui> coverSizes;
    std::vector<std::string> coverFiles;

    std::string patternSubDir = pathJoin(patternDir, queries[k]);
    std::string planSubDir = pathJoin(planDir, queries[k]);
    // read the equation and plan for each query
    for (int i = 0; i < patternGraphs.size(); ++i) {
        const PatternGraph &pg = patternGraphs[i];
        std::string treePath = pathJoin(planSubDir, files[i].substr(0, files[i].length() - 4) + "-tree");
        std::string equationPath = pathJoin(planSubDir, files[i].substr(0, files[i].length() - 4) + "-equation");
        std::string conNodePath = pathJoin(planSubDir, files[i].substr(0, files[i].length() - 4) + "-connode");
        std::ifstream equationStream(equationPath), conNodeStream(conNodePath);
        if (!equationStream.good()) {
            ConNode cn;
            cn.readFromStream(conNodeStream);
            cn.hashTables = std::vector<HashTable>(cn.aggrePos.size(), nullptr);
            for (int j = 0; j < cn.aggrePos.size(); ++j)
                cn.hashTables[j] = result[i];
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
        else {
            // read td
            trees[i].readFromFile(treePath);
            // read equation and the td of covers
            equationStream.read(reinterpret_cast<char*>(&divideFactors[i]), sizeof(divideFactors[i]));
            std::vector<ui> patternSizes;
            std::vector<std::string> patternIDs;
            std::vector<int> multiFactors;
            readVectorFromStream(equationStream, patternSizes);
            readVectorFromStream(equationStream, patternIDs);
            readVectorFromStream(equationStream, multiFactors);
            coverPos[i] = std::vector<int>(patternSizes.size());
            coverMultiFactors[i] = multiFactors;
            for (int j = 0; j < patternSizes.size(); ++j) {
                bool exists = false;
                std::string id = patternIDs[j];
                for (int l = 0; l < coverSizes.size(); ++l) {
                    if (coverSizes[l] == patternSizes[j] && coverFiles[l] == id) {
                        coverPos[i][j] = l;
                        exists = true;
                        break;
                    }
                }
                if (!exists) {
                    // load a cover pattern and its td
                    std::vector<std::string> tmp;
                    std::string coverPatternPath = pathJoin(patternDir, queries[patternSizes[j]]);
                    coverPatternPath = pathJoin(coverPatternPath, patternIDs[j] + ".txt");
                    coverPatterns.emplace_back(loadPatternGraph(coverPatternPath, false, tmp)[0]);
                    std::string coverTreePath = pathJoin(planDir, queries[patternSizes[j]]);
                    coverTreePath = pathJoin(coverTreePath, patternIDs[j] + "-tree");
                    Tree t;
                    t.readFromFile(coverTreePath);
                    coverTrees.push_back(t);
                    coverPos[i][j] = coverFiles.size();
                    coverSizes.push_back(patternSizes[j]);
                    coverFiles.push_back(patternIDs[j]);
                }
            }
        }
    }
}

int main(int argc, char **argv) {
//    generatePlan();

    Command cmd(argc, argv);
    std::string resultPath = cmd.getResultPath();
    std::string queryGraphPath = cmd.getQueryGraphPath();
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string trianglePath = cmd.getTrianglePath();
    DataGraph dun = DataGraph();
    dun.loadDataGraph(dataGraphPath);
    const DataGraph din = constructDirectedDataGraph(dun, false);
    const DataGraph dout = constructDirectedDataGraph(dun, true);
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    Triangle triangle;
    triangle.load(trianglePath, m / 2);
    EdgeID *outID = buildInID2OutID(din, dout);
    EdgeID *unID = buildUnID2OutID(dun, dout);
    EdgeID *reverseID = buildReverseUnID(dun);
    std::vector<std::string> files;
    std::vector<PatternGraph> patternGraphs = loadPatternGraph(queryGraphPath, true, files);
    int orbitType = patternGraphs[0].getOrbitType();
    int patternSize = patternGraphs[0].getNumVertices();
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
    std::vector<ConNode> conNodes;
    std::vector<Pattern> conPatterns;
    std::vector<int> conFactors;
    std::vector<int> conIDs;
    std::vector<Tree> coverTrees;
    std::vector<Pattern> coverPatterns;
    std::vector<std::vector<int>> coverPos(patternGraphs.size());
    std::vector<std::vector<int>> coverMultiFactors(patternGraphs.size());
    std::vector<HashTable> result(patternGraphs.size());
    for (int i = 0; i < patternGraphs.size(); ++i) {
        result[i] = new Count[n];
        memset(result[i], 0, sizeof(Count) * n);
    }
    std::vector<HashTable> coverH;
    std::vector<HashTable> uH;
    std::vector<int> divideFactors(patternGraphs.size());
    std::vector<Tree> trees(patternGraphs.size());
    std::vector<Pattern> patterns(patternGraphs.size());
    // read plan from files
    readPlan(patternSize, patternGraphs, files, conNodes, conPatterns, conFactors, conIDs, coverPos, coverTrees,
             coverPatterns, coverMultiFactors, result, trees, patterns, divideFactors);
    for (int i = 0; i < conNodes.size(); ++i) {
        executeConNode(conNodes[i], din, dout, dun, true, triangle, conPatterns[i], outID, unID,
                       startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
    }
    for (int i = 0; i < conIDs.size(); ++i) {
        int divideFactor = conFactors[i];
        int id = conIDs[i];
        for (int k = 0; k < n; ++k)
            result[id][k] /= divideFactor;
    }
    // covers
    coverH.resize(coverTrees.size());
    for (int i = 0; i < coverTrees.size(); ++i) {
        coverH[i] = new Count[n];
        memset(coverH[i], 0, sizeof(Count) * n);
    }
    for (int i = 0; i < coverTrees.size(); ++i) {
        const Tree &t = coverTrees[i];
        const Pattern &p = coverPatterns[i];
        for (int l = 0; l < t.getNumNodes(); ++l) memset(ht[l], 0, sizeof(Count) * m);
        if (t.getExecuteMode()) {
            executeTree(t, din, dout, dun, true, triangle, p, ht, outID, unID, reverseID,
                        startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
        }
        else {
            multiJoinTree(t, din, dout, dun, true, triangle, p, ht, outID, unID, reverseID, startOffset,
                          patternV, dataV, visited, tmp, allV);
        }
        HashTable h = ht[t.getRootID()];
        for (int j = 0; j < n; ++j) {
            coverH[i][j] = h[j];
        }
    }
    // patterns
    for (int i = 0; i < patternGraphs.size(); ++i) {
        if (std::find(conIDs.begin(), conIDs.end(), i) != conIDs.end()) continue;
        const Tree &t = trees[i];
        const Pattern &p = patterns[i];
        for (int l = 0; l < t.getNumNodes(); ++l) memset(ht[l], 0, sizeof(Count) * m);
        if (t.getExecuteMode()) {
            executeTree(t, din, dout, dun, true, triangle, p, ht, outID, unID, reverseID,
                        startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
        }
        else {
            multiJoinTree(t, din, dout, dun, true, triangle, p, ht, outID, unID, reverseID, startOffset,
                          patternV, dataV, visited, tmp, allV);
        }
        HashTable h = ht[t.getRootID()];
        int multiFactor = t.getMultiFactor();
        for (int j = 0; j < n; ++j) {
            result[i][j] = h[j] * multiFactor;
        }
    }
    // 3. subtract cover count
    for (int i = 0; i < patternGraphs.size(); ++i) {
        if (std::find(conIDs.begin(), conIDs.end(), i) != conIDs.end()) continue;
        const PatternGraph &pg = patternGraphs[i];
        int divideFactor = divideFactors[i];
        for (int j = 0; j < coverPos[i].size(); ++j) {
            HashTable h = coverH[coverPos[i][j]];
            for (int k = 0; k < n; ++k)
                result[i][k] += h[k] * coverMultiFactors[i][j];
        }
        for (int k = 0; k < n; ++k)
            result[i][k] /= divideFactor;
    }
    std::vector<int> orbitTypes(files.size(), orbitType);
    if (!resultPath.empty()) saveCount(resultPath, result, dun, true, files, orbitTypes);

    return 0;
}