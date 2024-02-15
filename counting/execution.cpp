//
// Created by anonymous author on 2022/8/29.
//

#include "execution.h"

Count gNumIntersect = 0;
Count gNumMatch = 0;
Count gNumIntermediate = 0;
Count gNumEdgeID = 0;
Count gNumUpdate = 0;

void saveCount(const std::string &resultPath, std::vector<HashTable> H, const DataGraph &d, bool batchQuery,
               const std::vector<std::string> &files, const std::vector<int> &orbitTypes) {
    if (batchQuery) {
        for (int i = 0; i < files.size(); ++i)
            saveCount(resultPath + "/" + files[i], H[i], d, orbitTypes[i]);
    }
    else saveCount(resultPath, H[0], d, orbitTypes[0]);
}

void saveCount(const std::string &resultPath, HashTable h, const DataGraph &d, int orbitType) {
    std::ofstream out(resultPath);
    if (orbitType == 0) out << h[0] << std::endl;
    else if (orbitType == 1) {
        for (VertexID i = 0; i < d.getNumVertices(); ++i)
            out << h[i] << std::endl;
    } else {
        for (VertexID u = 0; u < d.getNumVertices(); ++u) {
            ui nbrCnt = 0;
            VertexID *outNeighbors = d.getNeighborsLargerID(u, nbrCnt);
            for (EdgeID j = 0; j < nbrCnt; ++j) {
                VertexID v = outNeighbors[j];
                out << u << "\t\t\t" << v << "\t\t\t" << h[d.getEdgeID(u, v)];
            }
        }
    }
//    std::cout << "save result to " << resultPath << std::endl;
}

// give the current mapping, the (directed) pattern graph, the (directed) data graph, the next pattern vertex(pos),
// return the data vertices that can map to the pattern vertex
// inPos: positions of vertices in dataV that are in-neighbor of dataV[pos]
// note that we should avoid swapping neighbor with tmp
void generateCandidate(
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        VertexID *dataV,
        int pos,
        const std::vector<int> &inPos,
        const std::vector<int> &outPos,
        const std::vector<int> &unPos,
        VertexID **candidate,
        ui *&candCount,
        VertexID *&tmp
) {
    int num = inPos.size() + outPos.size() + unPos.size();
    ui tmpCount;
    if (num > 2) {
        VertexID **arrays = new VertexID *[num];
        ui *counts = new ui[num];
        int p = 0;
        for (int position: inPos) {
            arrays[p] = din.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        for (int position: outPos) {
            arrays[p] = dout.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        for (int position: unPos) {
            arrays[p] = dun.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        ComputeSetIntersection::LeapfrogJoin(arrays, counts, num, tmp, tmpCount);
        candCount[pos] = tmpCount;
        std::swap(candidate[pos], tmp);
        delete[] counts;
        delete[] arrays;
        return;
    }
    VertexID *firstNeighbor;
    ui firstCnt;
    bool firstFlag = true, secondFlag = true;
    for (int position: inPos) {
        ui nbrCnt;
        VertexID *neighbor = din.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstNeighbor = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
    for (int position: outPos) {
        ui nbrCnt;
        VertexID *neighbor = dout.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstNeighbor = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
    for (int position: unPos) {
        ui nbrCnt;
        VertexID *neighbor = dun.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstNeighbor = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstNeighbor, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
}

// generate candidate with the help of triangles stored
void generateCandidateT(
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        VertexID *dataV,
        int pos,
        bool isTriangle,
        EdgeID e,
        int type,
        const std::vector<int> &inPos,
        const std::vector<int> &outPos,
        const std::vector<int> &unPos,
        VertexID **candidate,
        ui *&candCount,
        VertexID *&tmp
) {
    if (isTriangle) {
        candidate[pos] = tri.getTriend(e, type, candCount[pos]);
        return;
    }
    int num = inPos.size() + outPos.size() + unPos.size() + (type != 0);
    ui tmpCount;
    if (num > 2) {
        VertexID **arrays = new VertexID *[num];
        ui *counts = new ui[num];
        int p = 0;
        if (type != 0) {
            arrays[p] = tri.getTriend(e, type, counts[p]);
            ++p;
        }
        for (int position: inPos) {
            arrays[p] = din.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        for (int position: outPos) {
            arrays[p] = dout.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        for (int position: unPos) {
            arrays[p] = dun.getNeighbors(dataV[position], counts[p]);
            ++p;
        }
        ComputeSetIntersection::LeapfrogJoin(arrays, counts, num, tmp, tmpCount);
        candCount[pos] = tmpCount;
        std::swap(candidate[pos], tmp);
        delete[] counts;
        delete[] arrays;
        return;
    }
    VertexID *firstList;
    ui firstCnt;
    bool firstFlag = true, secondFlag = true;
    if (type != 0) {
        firstList = tri.getTriend(e, type, firstCnt);
        firstFlag = false;
    }
    for (int position: inPos) {
        ui nbrCnt;
        VertexID *neighbor = din.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstList = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstList, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
    for (int position: outPos) {
        ui nbrCnt;
        VertexID *neighbor = dout.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstList = neighbor;
        }
        else {
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstList, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
    for (int position: unPos) {
        ui nbrCnt;
        VertexID *neighbor = dun.getNeighbors(dataV[position], nbrCnt);
        if (firstFlag) {
            firstFlag = false;
            firstCnt = nbrCnt;
            firstList = neighbor;
        }
        else {
            ui tmpCount;
            if (secondFlag) {
                secondFlag = false;
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(firstList, firstCnt, neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntersect;
#endif
                ComputeSetIntersection::ComputeCandidates(candidate[pos], candCount[pos], neighbor, nbrCnt, tmp, tmpCount);
                candCount[pos] = tmpCount;
                std::swap(candidate[pos], tmp);
            }
        }
    }
}

ui computeEdgeKey(
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        int mappingSize,
        ui *pos,
        int edgeType,
        VertexID src,
        VertexID dst
) {
    ui edgeID = 0;
    switch (edgeType) {
        case 1:
            edgeID = startOffset[mappingSize] + pos[mappingSize] - 1;
            break;
        case 2:
            edgeID = outID[startOffset[mappingSize] + pos[mappingSize] - 1];
            break;
        case 3:
#ifdef COLLECT_STATISTICS
            ++gNumEdgeID;
#endif
            edgeID = dout.getEdgeID(src, dst);
            break;
        case 4:
#ifdef COLLECT_STATISTICS
            ++gNumEdgeID;
#endif
            edgeID = dun.getUndirectedEID(src, dst);
            break;
        case 5:
            edgeID = unID[startOffset[mappingSize] + pos[mappingSize] - 1];
            break;
        case 6:
            edgeID = reverseID[startOffset[mappingSize] + pos[mappingSize] - 1];
            break;
    }
    return edgeID;
}

void executeNode(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        bool isRoot,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();

    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                              nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }

    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;

            if (depth == tau.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    cnt *= H[cID][dataV[childKeyPos[j][0]]];
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else if (orbitType == 1) {
                        for (int j = 0; j < aggreV.size(); ++j) {
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            VertexID key = dataV[aggrePos[j]];
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
#ifdef COLLECT_STATISTICS
                    ++gNumUpdate;
#endif
                    h[dataV[aggrePos[0]]] += cnt;
                    if (keyPosSize < sizeBound) {
                        keyPos[keyPosSize] = dataV[aggrePos[0]];
                        ++keyPosSize;
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                                      nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
}

void executeNodeT(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        bool isRoot,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    const std::vector<std::pair<int, int>> &nodeTriPos = t.getNodeTriPos(nID);
    const std::vector<int> &triEdgeType = t.getTriEdgeType(nID);
    const std::vector<int> &triEndType = t.getTriEndType(nID);
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();

    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            EdgeID e = 0;
            if (triEdgeType[mappingSize] != 0) {
                int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                if (triEdgeType[mappingSize] == 1)
                    e = startOffset[pos2] + pos[pos2] - 1;
                else if (triEdgeType[mappingSize] == 2)
                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                else if (triEdgeType[mappingSize] == 5)
                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                else {
                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                    ++gNumEdgeID;
#endif
                    e = dout.getEdgeID(key1, key2);
                }
            }
            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                               nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            if (depth == tau.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    cnt *= H[cID][dataV[childKeyPos[j][0]]];
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else if (orbitType == 1) {
                        for (int j = 0; j < aggreV.size(); ++j) {
                            VertexID key = dataV[aggrePos[j]];
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
#ifdef COLLECT_STATISTICS
                    ++gNumUpdate;
#endif
                    h[dataV[aggrePos[0]]] += cnt;
                    if (keyPosSize < sizeBound) {
                        keyPos[keyPosSize] = dataV[aggrePos[0]];
                        ++keyPosSize;
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    EdgeID e = 0;
                    if (triEdgeType[mappingSize] != 0) {
                        int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                        if (triEdgeType[mappingSize] == 1)
                            e = startOffset[pos2] + pos[pos2] - 1;
                        else if (triEdgeType[mappingSize] == 2)
                            e = outID[startOffset[pos2] + pos[pos2] - 1];
                        else if (triEdgeType[mappingSize] == 5)
                            e = unID[startOffset[pos2] + pos[pos2] - 1];
                        else{
                            VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                            ++gNumEdgeID;
#endif
                            e = dout.getEdgeID(key1, key2);
                        }
                    }
                    generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                                       nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
}

void executeNodeEdgeKey(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        bool isRoot,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    const std::vector<std::vector<int>> &posChildEdge = t.getPosChildEdge(nID);
    const std::vector<std::vector<int>> &posAggreEdge = t.getPosAggreEdge(nID);
    const std::vector<int> &childEdgeType = t.getChildEdgeType(nID);
    const std::vector<int> &aggreEdgeType = t.getAggreEdgeType(nID);
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();

    std::vector<EdgeID> childKey(child.size());
    std::vector<EdgeID> aggreKey;
    if (isRoot && orbitType == 2) aggreKey = std::vector<EdgeID>(aggreWeight.size());
    else if (tau.keySize == 2) aggreKey.push_back(0);
    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                              nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    for (int i = 1; i < mappingSize; ++i) {
        // for child edge keys, compute the key value
        for (int j: posChildEdge[i]) {
            childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
        }
        // for aggregation edge keys, compute the key value
        for (int j: posAggreEdge[i]) {
            aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
        }
    }
    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            // for child edge keys, compute the key value
            for (int j: posChildEdge[mappingSize]) {
                childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
            }
            // for aggregation edge keys, compute the key value
            for (int j: posAggreEdge[mappingSize]) {
                aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
            }
            if (depth == tau.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    if (childKeyPos[j].size() == 1) {
                        cnt *= H[cID][dataV[childKeyPos[j][0]]];
                    }
                    else {
                        cnt *= H[cID][childKey[j]];
                    }
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else if (orbitType == 1) {
                        for (int j = 0; j < aggreV.size(); ++j) {
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            VertexID key = dataV[aggrePos[j]];
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                    else {
                        for (int j = 0; j < aggreKey.size(); ++j) {
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[aggreKey[j]] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
                    if (tau.keySize < 2) {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePos[0]]] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = dataV[aggrePos[0]];
                            ++keyPosSize;
                        }
                    }
                    else {
                        EdgeID e = aggreKey[0];
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[e] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = e;
                            ++keyPosSize;
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                                      nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
}

void executeNodeEdgeKeyT(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID **candidate,
        ui *candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        bool isRoot,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        ui *keyPos,
        ui &keyPosSize,
        ui sizeBound,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    const std::vector<std::vector<int>> &posChildEdge = t.getPosChildEdge(nID);
    const std::vector<std::vector<int>> &posAggreEdge = t.getPosAggreEdge(nID);
    const std::vector<int> &childEdgeType = t.getChildEdgeType(nID);
    const std::vector<int> &aggreEdgeType = t.getAggreEdgeType(nID);
    const std::vector<std::pair<int, int>> &nodeTriPos = t.getNodeTriPos(nID);
    const std::vector<int> &triEdgeType = t.getTriEdgeType(nID);
    const std::vector<int> &triEndType = t.getTriEndType(nID);
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    std::vector<EdgeID> childKey(child.size());
    std::vector<EdgeID> aggreKey;
    if (isRoot && orbitType == 2) aggreKey = std::vector<EdgeID>(aggreWeight.size());
    else if (tau.keySize == 2) aggreKey.push_back(0);
    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            EdgeID e = 0;
            if (triEdgeType[mappingSize] != 0) {
                int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                if (triEdgeType[mappingSize] == 1)
                    e = startOffset[pos2] + pos[pos2] - 1;
                else if (triEdgeType[mappingSize] == 2)
                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                else if (triEdgeType[mappingSize] == 5)
                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                else {
                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                    ++gNumEdgeID;
#endif
                    e = dout.getEdgeID(key1, key2);
                }
            }
            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                               nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    for (int i = 1; i < mappingSize; ++i) {
        // for child edge keys, compute the key value
        for (int j: posChildEdge[i]) {
            childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
        }
        // for aggregation edge keys, compute the key value
        for (int j: posAggreEdge[i]) {
            aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
        }
    }
    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            // for child edge keys, compute the key value
            for (int j: posChildEdge[mappingSize]) {
                childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
            }
            // for aggregation edge keys, compute the key value
            for (int j: posAggreEdge[mappingSize]) {
                aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
            }
            if (depth == tau.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    if (childKeyPos[j].size() == 1) {
                        cnt *= H[cID][dataV[childKeyPos[j][0]]];
                    }
                    else {
                        cnt *= H[cID][childKey[j]];
                    }
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else if (orbitType == 1) {
                        for (int j = 0; j < aggreV.size(); ++j) {
                            VertexID key = dataV[aggrePos[j]];
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                    else {
                        for (int j = 0; j < aggreKey.size(); ++j) {
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[aggreKey[j]] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
                    if (tau.keySize < 2) {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePos[0]]] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = dataV[aggrePos[0]];
                            ++keyPosSize;
                        }
                    }
                    else {
                        EdgeID e = aggreKey[0];
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[e] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = e;
                            ++keyPosSize;
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    EdgeID e = 0;
                    if (triEdgeType[mappingSize] != 0) {
                        int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                        if (triEdgeType[mappingSize] == 1)
                            e = startOffset[pos2] + pos[pos2] - 1;
                        else if (triEdgeType[mappingSize] == 2)
                            e = outID[startOffset[pos2] + pos[pos2] - 1];
                        else if (triEdgeType[mappingSize] == 5)
                            e = unID[startOffset[pos2] + pos[pos2] - 1];
                        else {
                            VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                            ++gNumEdgeID;
#endif
                            e = dout.getEdgeID(key1, key2);
                        }
                    }
                    generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                                       nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
}

void executePartition(
        VertexID pID,
        const Tree &t,
        VertexID ***candidate,
        ui **candCount,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool useTriangle,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    const std::vector<std::vector<VertexID>> &globalOrder = t.getGlobalOrder();
    const std::vector<std::vector<std::vector<VertexID>>> &nodesAtStep = t.getNodesAtStep();
    const std::vector<VertexID> &partitionOrder = globalOrder[pID];
    const std::vector<std::vector<VertexID>> &child = t.getChild();
    const std::vector<VertexID> &postOrder = t.getPostOrder();
    const std::vector<int> &partitionPos = t.getPartitionPos();
    const std::vector<std::vector<int>> &partitionInPos = t.getPartitionInPos(pID);
    const std::vector<std::vector<int>> &partitionOutPos = t.getPartitionOutPos(pID);
    const std::vector<std::vector<int>> &partitionUnPos = t.getPartitionUnPos(pID);
    const std::vector<bool> &partitionInterPos = t.getPartitionInterPos(pID);
    const std::vector<std::vector<int>> &greaterPos = t.getPartitionGreaterPos(pID);
    const std::vector<std::vector<int>> &lessPos = t.getPartitionLessPos(pID);
    const std::vector<bool> &partitionCandPos = t.getPartitionCandPos(pID);
    const std::vector<std::pair<int, int>> &partitionTriPos = t.getPartitionTriPos(pID);
    const std::vector<int> &triEdgeType = t.getPartitionEdgeType(pID);
    const std::vector<int> &triEndType = t.getPartitionEndType(pID);
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    EdgeID *unNbors = dun.getNbors();
    int startPos;
    if (pID == 0) startPos = 0;
    else startPos = partitionPos[pID - 1] + 1;
    int endPos = partitionPos[pID] + 1;
    bool isRoot = endPos == postOrder.size();
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    int mappingSize = 0;
    const std::vector<std::vector<VertexID>> &allChild = t.getChild();
    ui argument = 0;
    if (partitionOrder.empty()) {
        for (int i = startPos; i < endPos; ++i) {
            VertexID nID = postOrder[i];
            isRoot = endPos == postOrder.size() && i == endPos - 1;
            if (!t.nodeEdgeKey(nID) && !useTriangle) {
                executeNode(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, p, isRoot, outID, unID, reverseID,
                            startOffset, patternV, dataV, 0, visited, pos, nullptr, argument, argument, tmp, allV);
            }
            else if (!t.nodeEdgeKey(nID) && useTriangle) {
                executeNodeT(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                             startOffset, patternV, dataV, 0, visited, pos, nullptr, argument, argument, tmp, allV);

            }
            else if (t.nodeEdgeKey(nID) && !useTriangle) {
                executeNodeEdgeKey(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, p, isRoot, outID, unID, reverseID,
                                   startOffset, patternV, dataV, 0, visited, pos, nullptr, argument, argument, tmp, allV);
            }
            else {
                executeNodeEdgeKeyT(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                                    startOffset, patternV, dataV, 0, visited, pos, nullptr, argument, argument, tmp, allV);
            }
        }
        return;
    }
    // let candidate points to a memory whose length is the number of vertices in the data graph
    VertexID **partitionCandidate = new VertexID *[partitionOrder.size()];
    for (int i = 0; i < partitionOrder.size(); ++i) {
        if (partitionCandPos[i])
            partitionCandidate[i] = new VertexID[dout.getNumVertices()];
    }
    ui *partitionCandCount = new ui[partitionOrder.size()];
    // for each node, store the keys whose values are not 0
    ui **keyPos = new ui *[t.getNumNodes()];
    for (int i = startPos; i < endPos; ++i) {
        VertexID nID = postOrder[i];
        if (t.getNode(nID).keySize == 0)
            keyPos[nID] = new ui[1];
        else if (t.getNode(nID).keySize == 1)
            keyPos[nID] = new ui[n / 8 + 1];
        else if (t.getNode(nID).keySize == 2)
            keyPos[nID] = new ui[m / 8 + 1];
    }

    ui *keyPosSize = new ui[t.getNumNodes()];
    memset(keyPosSize, 0, sizeof(ui) * (t.getNumNodes()));
    partitionCandidate[0] = allV;
    partitionCandCount[0] = n;
    pos[0] = 0;
    // depth iterate through the pattern vertices in the partition order
    int depth = 0;
    while (depth >= 0) {
        while (pos[depth] < partitionCandCount[depth]) {
            VertexID v = partitionCandidate[depth][pos[depth]];
            ++pos[depth];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[depth] = partitionOrder[depth];
            dataV[depth] = v;
            for (VertexID nID: nodesAtStep[pID][depth]) {
                if (nID == postOrder[endPos - 1]) {
                    if (!t.nodeEdgeKey(nID) && !useTriangle) {
                        executeNode(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, p, isRoot, outID, unID, reverseID,
                                    startOffset, patternV, dataV, mappingSize + 1, visited, pos, nullptr, argument, argument, tmp, allV);
                    }
                    else if (!t.nodeEdgeKey(nID) && useTriangle) {
                        executeNodeT(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                                     startOffset, patternV, dataV, mappingSize + 1, visited, pos, nullptr, argument, argument, tmp, allV);
                    }
                    else if (t.nodeEdgeKey(nID) && !useTriangle) {
                        executeNodeEdgeKey(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, p, isRoot, outID, unID, reverseID,
                                           startOffset, patternV, dataV, mappingSize + 1, visited, pos, nullptr, argument, argument, tmp, allV);
                    }
                    else {
                        executeNodeEdgeKeyT(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                                            startOffset, patternV, dataV, mappingSize + 1, visited, pos, nullptr, argument, argument, tmp, allV);
                    }
                }
                else if (t.getNode(nID).keySize < 2) {
                    if (keyPosSize[nID] < n / 8 + 1) {
                        for (int j = 0; j < keyPosSize[nID]; ++j) {
                            H[nID][keyPos[nID][j]] = 0;
                        }
                    }
                    else
                        memset(H[nID], 0, sizeof(Count) * n);
                    keyPosSize[nID] = 0;
                    ui sizeBound = 1;
                    if (t.getNode(nID).keySize == 1) sizeBound = n / 8 + 1;
                    if (t.getNode(nID).keySize == 2) sizeBound = m / 8 + 1;
                    if (!t.nodeEdgeKey(nID) && !useTriangle) {
                        executeNode(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, p, false, outID, unID, reverseID,
                                    startOffset, patternV, dataV, mappingSize + 1, visited, pos, keyPos[nID], keyPosSize[nID], sizeBound, tmp, allV);
                    }
                    else if (!t.nodeEdgeKey(nID) && useTriangle) {
                        executeNodeT(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, false, outID, unID, reverseID,
                                     startOffset, patternV, dataV, mappingSize + 1, visited, pos, keyPos[nID], keyPosSize[nID], sizeBound, tmp, allV);
                    }
                    else if (t.nodeEdgeKey(nID) && !useTriangle) {
                        executeNodeEdgeKey(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, p, false, outID, unID, reverseID,
                                           startOffset, patternV, dataV, mappingSize + 1, visited, pos, keyPos[nID], keyPosSize[nID], sizeBound, tmp, allV);
                    }
                    else {
                        executeNodeEdgeKeyT(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, false, outID, unID, reverseID,
                                            startOffset, patternV, dataV, mappingSize + 1, visited, pos, keyPos[nID], keyPosSize[nID], sizeBound, tmp, allV);
                    }
                }
                else {
                    if (keyPosSize[nID] < m / 8) {
                        for (int j = 0; j < keyPosSize[nID]; ++j) {
                            H[nID][keyPos[nID][j]] = 0;
                        }
                    }
                    else
                        memset(H[nID], 0, sizeof(Count) * m);
                    keyPosSize[nID] = 0;
                    if (!useTriangle) {
                        executeNodeEdgeKey(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, p, false, outID, unID, reverseID, startOffset,
                                           patternV, dataV, mappingSize + 1, visited, pos, keyPos[nID], keyPosSize[nID], m / 8, tmp, allV);
                    }
                    else {
                        executeNodeEdgeKeyT(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, false, outID, unID, reverseID, startOffset,
                                            patternV, dataV, mappingSize + 1, visited, pos, keyPos[nID], keyPosSize[nID], m / 8, tmp, allV);
                    }
                }
            }
            if (depth == partitionOrder.size() - 1) {
                visited[dataV[mappingSize]] = false;
            }
            else {
                ++mappingSize;
                ++depth;
                pos[depth] = 0;
                if (!partitionInterPos[depth]) {
                    if (!partitionOutPos[depth].empty()) {
                        VertexID w = dataV[partitionOutPos[depth][0]];
                        startOffset[mappingSize] = outOffset[w];
                        partitionCandidate[depth] = outNbors + outOffset[w];
                        partitionCandCount[depth] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!partitionInPos[depth].empty()){
                        VertexID w = dataV[partitionInPos[depth][0]];
                        startOffset[mappingSize] = inOffset[w];
                        partitionCandidate[depth] = inNbors + inOffset[w];
                        partitionCandCount[depth] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[partitionUnPos[depth][0]];
                        startOffset[mappingSize] = unOffset[w];
                        partitionCandidate[depth] = unNbors + unOffset[w];
                        partitionCandCount[depth] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    if (useTriangle) {
                        EdgeID e = 0;
                        if (triEdgeType[mappingSize] != 0) {
                            int pos1 = partitionTriPos[mappingSize].first, pos2 = partitionTriPos[mappingSize].second;
                            if (triEdgeType[mappingSize] == 1)
                                e = startOffset[pos2] + pos[pos2] - 1;
                            else if (triEdgeType[mappingSize] == 2)
                                e = outID[startOffset[pos2] + pos[pos2] - 1];
                            else if (triEdgeType[mappingSize] == 5)
                                e = unID[startOffset[pos2] + pos[pos2] - 1];
                            else{
                                VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                                ++gNumEdgeID;
#endif
                                e = dout.getEdgeID(key1, key2);
                            }
                        }
                        generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !partitionCandPos[mappingSize], e, triEndType[mappingSize],
                                           partitionInPos[mappingSize], partitionOutPos[mappingSize], partitionUnPos[mappingSize], partitionCandidate, partitionCandCount, tmp);

                    }
                    else {
                        generateCandidate(din, dout, dun, dataV, depth, partitionInPos[depth], partitionOutPos[depth],
                                          partitionUnPos[depth], partitionCandidate, partitionCandCount, tmp);
                    }
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(partitionCandidate[mappingSize], 0, partitionCandCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    partitionCandCount[mappingSize] = firstPosGreaterThan(partitionCandidate[mappingSize], 0, partitionCandCount[mappingSize], minTarget);
                }

            }
        }
        --depth;
        --mappingSize;
        if (depth >= 0)
            visited[dataV[depth]] = false;
    }
    for (int j = 0; j < partitionOrder.size(); ++j) {
        if (partitionCandPos[j])
            delete[] partitionCandidate[j];
    }
    delete[] partitionCandidate;
    delete[] partitionCandCount;
    for (int j = startPos; j < endPos; ++j) {
        VertexID nID = postOrder[j];
        delete[] keyPos[postOrder[j]];
    }
    delete[] keyPos;
    delete[] keyPosSize;
}

void executeTree (
        const Tree &t,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool useTriangle,
        const Triangle &tri,
        const Pattern &p,
        HashTable *H,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    int numNodes = (int)t.getNumNodes();
    VertexID ***candidate = new VertexID **[numNodes];
    ui **candCount = new ui *[numNodes];
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        const Node &tau = t.getNode(nID);
        int partOrderLength = int(tau.nodeOrder.size() - tau.localOrder.size());
        candCount[nID] = new ui[tau.nodeOrder.size()];
        candidate[nID] = new VertexID *[tau.nodeOrder.size()];
        const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
        for (int i = partOrderLength; i < nodeCandPos.size(); ++i) {
            if (nodeCandPos[i])
                candidate[nID][i] = new VertexID[dout.getNumVertices()];
        }
    }

    // for each partition, calls executePartition to build a full hash table
    for (VertexID pID = 0; pID < t.getPartitionPos().size(); ++pID) {
        executePartition(pID, t, candidate, candCount, H, din, dout, dun, useTriangle, tri,
                         p, outID, unID, reverseID, startOffset, patternV, dataV, visited, pos, tmp, allV);
    }
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        const Node &tau = t.getNode(nID);
        int partOrderLength = int(tau.nodeOrder.size() - tau.localOrder.size());
        const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
        for (int i = partOrderLength; i < nodeCandPos.size(); ++i) {
            if (nodeCandPos[i])
                delete[] candidate[nID][i];
        }
        delete[] candidate[nID];
        delete[] candCount[nID];
    }
    delete[] candCount;
    delete[] candidate;
}

void multiJoin(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    const std::vector<std::vector<VertexID>> &nodesAtStep = t.getNodesAtStep()[nID];
    const std::vector<std::vector<int>> &prefixPos = t.getPrefixPos();
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    VertexID **candidate = candidates[nID];
    ui *candCount = candCounts[nID];
    ui *pos = poses[nID];
    VertexID *patternV = patternVs[nID];
    VertexID *dataV = dataVs[nID];
    EdgeID *startOffset = startOffsets[nID];
    bool *visited = visits[nID];
    ui *keyPos = keyPoses[nID];
    ui &keyPosSize = keyPosSizes[nID];
    ui sizeBound = sizeBounds[nID];
    bool isRoot = nID == t.getRootID();
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    if (tau.prefixSize != 0) {
        for (int i = 0; i < nodesAtStep[0].size(); ++i) {
            VertexID cID = nodesAtStep[0][i];
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                patternVs[cID][j] = patternV[prefixPos[cID][j]];
                dataVs[cID][j] = dataV[prefixPos[cID][j]];
                startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                poses[cID][j] = pos[prefixPos[cID][j]];
                visits[cID][dataVs[cID][j]] = true;
            }
            if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
            else if (t.getNode(cID).keySize == 1) {
                if (keyPosSizes[cID] < n / 8 + 1) {
                    for (int j = 0; j < keyPosSizes[cID]; ++j) {
                        H[cID][keyPoses[cID][j]] = 0;
                    }
                }
                else memset(H[cID], 0, sizeof(Count) * n);
            }
            else {
                if (keyPosSizes[cID] < m / 8 + 1) {
                    for (int j = 0; j < keyPosSizes[cID]; ++j) {
                        H[cID][keyPoses[cID][j]] = 0;
                    }
                }
                else memset(H[cID], 0, sizeof(Count) * m);
            }
            multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                visits[cID][dataVs[cID][j]] = false;
            }
        }
    }
    int mappingSize = tau.prefixSize;
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[0] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                              nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    pos[0] = 0;
    while (mappingSize >= (int)tau.prefixSize) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            for (int i = 0; i < nodesAtStep[mappingSize].size(); ++i) {
                VertexID cID = nodesAtStep[mappingSize][i];
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    patternVs[cID][j] = patternV[prefixPos[cID][j]];
                    dataVs[cID][j] = dataV[prefixPos[cID][j]];
                    startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                    poses[cID][j] = pos[prefixPos[cID][j]];
                    visits[cID][dataVs[cID][j]] = true;
                }
                if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
                else if (t.getNode(cID).keySize == 1) {
                    if (keyPosSizes[cID] < n / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * n);
                }
                else {
                    if (keyPosSizes[cID] < m / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * m);
                }
                multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    visits[cID][dataVs[cID][j]] = false;
                }
            }
            if (mappingSize == tau.numVertices - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    if (childKeyPos[j].empty()) cnt *= H[cID][0];
                    else cnt *= H[cID][dataV[childKeyPos[j][0]]];
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else {
                        for (int j = 0; j < aggreV.size(); ++j) {
                            VertexID key = dataV[aggrePos[j]];
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
                    if (tau.keySize == 0) h[0] += cnt;
                    else {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePos[0]]] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = dataV[aggrePos[0]];
                            ++keyPosSize;
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                                      nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --mappingSize;
        if (mappingSize >= (int)tau.prefixSize) {
            visited[dataV[mappingSize]] = false;
        }
    }
}

void multiJoinT(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    const std::vector<std::pair<int, int>> &nodeTriPos = t.getNodeTriPos(nID);
    const std::vector<int> &triEdgeType = t.getTriEdgeType(nID);
    const std::vector<int> &triEndType = t.getTriEndType(nID);
    const std::vector<std::vector<VertexID>> &nodesAtStep = t.getNodesAtStep()[nID];
    const std::vector<std::vector<int>> &prefixPos = t.getPrefixPos();
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    VertexID **candidate = candidates[nID];
    ui *candCount = candCounts[nID];
    ui *pos = poses[nID];
    VertexID *patternV = patternVs[nID];
    VertexID *dataV = dataVs[nID];
    EdgeID *startOffset = startOffsets[nID];
    bool *visited = visits[nID];
    ui *keyPos = keyPoses[nID];
    ui &keyPosSize = keyPosSizes[nID];
    ui sizeBound = sizeBounds[nID];
    bool isRoot = nID == t.getRootID();
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    if (tau.prefixSize != 0) {
        for (int i = 0; i < nodesAtStep[0].size(); ++i) {
            VertexID cID = nodesAtStep[0][i];
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                patternVs[cID][j] = patternV[prefixPos[cID][j]];
                dataVs[cID][j] = dataV[prefixPos[cID][j]];
                startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                poses[cID][j] = pos[prefixPos[cID][j]];
                visits[cID][dataVs[cID][j]] = true;
            }
            if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
            else {
                if (keyPosSizes[cID] < n / 8 + 1) {
                    for (int j = 0; j < keyPosSizes[cID]; ++j) {
                        H[cID][keyPoses[cID][j]] = 0;
                    }
                }
                else memset(H[cID], 0, sizeof(Count) * n);
            }
            multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                visits[cID][dataVs[cID][j]] = false;
            }
        }
    }
    int mappingSize = tau.prefixSize;
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[0] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            EdgeID e = 0;
            if (triEdgeType[mappingSize] != 0) {
                int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                if (triEdgeType[mappingSize] == 1)
                    e = startOffset[pos2] + pos[pos2] - 1;
                else if (triEdgeType[mappingSize] == 2)
                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                else if (triEdgeType[mappingSize] == 5)
                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                else {
                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                    ++gNumEdgeID;
#endif
                    e = dout.getEdgeID(key1, key2);
                }
            }
            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                               nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    pos[0] = 0;
    while (mappingSize >= (int)tau.prefixSize) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            for (int i = 0; i < nodesAtStep[mappingSize].size(); ++i) {
                VertexID cID = nodesAtStep[mappingSize][i];
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    patternVs[cID][j] = patternV[prefixPos[cID][j]];
                    dataVs[cID][j] = dataV[prefixPos[cID][j]];
                    startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                    poses[cID][j] = pos[prefixPos[cID][j]];
                    visits[cID][dataVs[cID][j]] = true;
                }
                if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
                else if (t.getNode(cID).keySize == 1) {
                    if (keyPosSizes[cID] < n / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * n);
                }
                else {
                    if (keyPosSizes[cID] < m / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * m);
                }
                multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    visits[cID][dataVs[cID][j]] = false;
                }
            }
            if (mappingSize == tau.numVertices - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    if (childKeyPos[j].empty()) cnt *= H[cID][0];
                    else cnt *= H[cID][dataV[childKeyPos[j][0]]];
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else {
                        for (int j = 0; j < aggreV.size(); ++j) {
                            VertexID key = dataV[aggrePos[j]];
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
                    if (tau.keySize == 0) h[0] += cnt;
                    else {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePos[0]]] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = dataV[aggrePos[0]];
                            ++keyPosSize;
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    EdgeID e = 0;
                    if (triEdgeType[mappingSize] != 0) {
                        int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                        if (triEdgeType[mappingSize] == 1)
                            e = startOffset[pos2] + pos[pos2] - 1;
                        else if (triEdgeType[mappingSize] == 2)
                            e = outID[startOffset[pos2] + pos[pos2] - 1];
                        else if (triEdgeType[mappingSize] == 5)
                            e = unID[startOffset[pos2] + pos[pos2] - 1];
                        else {
                            VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                            ++gNumEdgeID;
#endif
                            e = dout.getEdgeID(key1, key2);
                        }
                    }
                    generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                                       nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --mappingSize;
        if (mappingSize >= (int)tau.prefixSize) {
            visited[dataV[mappingSize]] = false;
        }
    }
}

void multiJoinE(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    const std::vector<std::vector<int>> &posChildEdge = t.getPosChildEdge(nID);
    const std::vector<std::vector<int>> &posAggreEdge = t.getPosAggreEdge(nID);
    const std::vector<int> &childEdgeType = t.getChildEdgeType(nID);
    const std::vector<int> &aggreEdgeType = t.getAggreEdgeType(nID);
    const std::vector<std::vector<VertexID>> &nodesAtStep = t.getNodesAtStep()[nID];
    const std::vector<std::vector<int>> &prefixPos = t.getPrefixPos();
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    VertexID **candidate = candidates[nID];
    ui *candCount = candCounts[nID];
    ui *pos = poses[nID];
    VertexID *patternV = patternVs[nID];
    VertexID *dataV = dataVs[nID];
    EdgeID *startOffset = startOffsets[nID];
    bool *visited = visits[nID];
    ui *keyPos = keyPoses[nID];
    ui &keyPosSize = keyPosSizes[nID];
    ui sizeBound = sizeBounds[nID];
    bool isRoot = nID == t.getRootID();
    std::vector<EdgeID> childKey(child.size());
    std::vector<EdgeID> aggreKey;
    if (isRoot && orbitType == 2) aggreKey = std::vector<EdgeID>(aggreWeight.size());
    else if (tau.keySize == 2) aggreKey.push_back(0);
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    if (tau.prefixSize != 0) {
        for (int i = 0; i < nodesAtStep[0].size(); ++i) {
            VertexID cID = nodesAtStep[0][i];
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                patternVs[cID][j] = patternV[prefixPos[cID][j]];
                dataVs[cID][j] = dataV[prefixPos[cID][j]];
                startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                poses[cID][j] = pos[prefixPos[cID][j]];
                visits[cID][dataVs[cID][j]] = true;
            }
            if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
            else if (t.getNode(cID).keySize == 1) {
                if (keyPosSizes[cID] < n / 8 + 1) {
                    for (int j = 0; j < keyPosSizes[cID]; ++j) {
                        H[cID][keyPoses[cID][j]] = 0;
                    }
                }
                else memset(H[cID], 0, sizeof(Count) * n);
            }
            else {
                if (keyPosSizes[cID] < m / 8 + 1) {
                    for (int j = 0; j < keyPosSizes[cID]; ++j) {
                        H[cID][keyPoses[cID][j]] = 0;
                    }
                }
                else memset(H[cID], 0, sizeof(Count) * m);
            }
            multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                visits[cID][dataVs[cID][j]] = false;
            }
        }
    }
    int mappingSize = tau.prefixSize;
    for (int i = 1; i < mappingSize; ++i) {
        // for child edge keys, compute the key value
        for (int j: posChildEdge[i]) {
            childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
        }
        // for aggregation edge keys, compute the key value
        for (int j: posAggreEdge[i]) {
            aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
        }
    }
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[0] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                              nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    pos[0] = 0;
    while (mappingSize >= (int)tau.prefixSize) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            // for child edge keys, compute the key value
            for (int j: posChildEdge[mappingSize]) {
                childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
            }
            // for aggregation edge keys, compute the key value
            for (int j: posAggreEdge[mappingSize]) {
                aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
            }
            for (int i = 0; i < nodesAtStep[mappingSize].size(); ++i) {
                VertexID cID = nodesAtStep[mappingSize][i];
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    patternVs[cID][j] = patternV[prefixPos[cID][j]];
                    dataVs[cID][j] = dataV[prefixPos[cID][j]];
                    startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                    poses[cID][j] = pos[prefixPos[cID][j]];
                    visits[cID][dataVs[cID][j]] = true;
                }
                if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
                else if (t.getNode(cID).keySize == 1) {
                    if (keyPosSizes[cID] < n / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * n);
                }
                else {
                    if (keyPosSizes[cID] < m / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * m);
                }
                multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    visits[cID][dataVs[cID][j]] = false;
                }
            }
            if (mappingSize == tau.numVertices - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    if (childKeyPos[j].empty()) cnt *= H[cID][0];
                    else if (childKeyPos[j].size() == 1) {
                        cnt *= H[cID][dataV[childKeyPos[j][0]]];
                    }
                    else {
                        cnt *= H[cID][childKey[j]];
                    }
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else if (orbitType == 1) {
                        for (int j = 0; j < aggreV.size(); ++j) {
                            VertexID key = dataV[aggrePos[j]];
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                    else {
                        for (int j = 0; j < aggreKey.size(); ++j) {
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[aggreKey[j]] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
                    if (tau.keySize == 0) h[0] += cnt;
                    else if (tau.keySize == 1) {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePos[0]]] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = dataV[aggrePos[0]];
                            ++keyPosSize;
                        }
                    }
                    else {
                        EdgeID e = aggreKey[0];
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[e] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = e;
                            ++keyPosSize;
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                                      nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --mappingSize;
        if (mappingSize >= (int)tau.prefixSize) {
            visited[dataV[mappingSize]] = false;
        }
    }
}

void multiJoinET(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = t.getOrbitType();
    HashTable h = H[nID];
    const Node &tau = t.getNode(nID);
    const std::vector<int> &aggrePos = t.getAggrePos(nID);
    const std::vector<bool> &nodeInterPos = t.getNodeInterPos(nID);
    const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
    const std::vector<std::vector<int>> &nodeInPos = t.getNodeInPos(nID);
    const std::vector<std::vector<int>> &nodeOutPos = t.getNodeOutPos(nID);
    const std::vector<std::vector<int>> &nodeUnPos = t.getNodeUnPos(nID);
    const std::vector<std::vector<int>> &greaterPos = t.getNodeGreaterPos(nID);
    const std::vector<std::vector<int>> &lessPos = t.getNodeLessPos(nID);
    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
    const std::vector<VertexID> &aggreV = t.getAggreV();
    const std::vector<int> &aggreWeight = t.getAggreWeight();
    const std::vector<std::vector<int>> &posChildEdge = t.getPosChildEdge(nID);
    const std::vector<std::vector<int>> &posAggreEdge = t.getPosAggreEdge(nID);
    const std::vector<int> &childEdgeType = t.getChildEdgeType(nID);
    const std::vector<int> &aggreEdgeType = t.getAggreEdgeType(nID);
    const std::vector<std::pair<int, int>> &nodeTriPos = t.getNodeTriPos(nID);
    const std::vector<int> &triEdgeType = t.getTriEdgeType(nID);
    const std::vector<int> &triEndType = t.getTriEndType(nID);
    const std::vector<std::vector<VertexID>> &nodesAtStep = t.getNodesAtStep()[nID];
    const std::vector<std::vector<int>> &prefixPos = t.getPrefixPos();
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    VertexID **candidate = candidates[nID];
    ui *candCount = candCounts[nID];
    ui *pos = poses[nID];
    VertexID *patternV = patternVs[nID];
    VertexID *dataV = dataVs[nID];
    EdgeID *startOffset = startOffsets[nID];
    bool *visited = visits[nID];
    ui *keyPos = keyPoses[nID];
    ui &keyPosSize = keyPosSizes[nID];
    ui sizeBound = sizeBounds[nID];
    bool isRoot = nID == t.getRootID();
    std::vector<EdgeID> childKey(child.size());
    std::vector<EdgeID> aggreKey;
    if (isRoot && orbitType == 2) aggreKey = std::vector<EdgeID>(aggreWeight.size());
    else if (tau.keySize == 2) aggreKey.push_back(0);
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    if (tau.prefixSize != 0) {
        for (int i = 0; i < nodesAtStep[0].size(); ++i) {
            VertexID cID = nodesAtStep[0][i];
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                patternVs[cID][j] = patternV[prefixPos[cID][j]];
                dataVs[cID][j] = dataV[prefixPos[cID][j]];
                startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                poses[cID][j] = pos[prefixPos[cID][j]];
                visits[cID][dataVs[cID][j]] = true;
            }
            if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
            else if (t.getNode(cID).keySize == 1) {
                if (keyPosSizes[cID] < n / 8 + 1) {
                    for (int j = 0; j < keyPosSizes[cID]; ++j) {
                        H[cID][keyPoses[cID][j]] = 0;
                    }
                }
                else memset(H[cID], 0, sizeof(Count) * n);
            }
            else {
                if (keyPosSizes[cID] < m / 8 + 1) {
                    for (int j = 0; j < keyPosSizes[cID]; ++j) {
                        H[cID][keyPoses[cID][j]] = 0;
                    }
                }
                else memset(H[cID], 0, sizeof(Count) * m);
            }
            multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
            for (int j = 0; j < prefixPos[cID].size(); ++j) {
                visits[cID][dataVs[cID][j]] = false;
            }
        }
    }
    int mappingSize = tau.prefixSize;
    for (int i = 1; i < mappingSize; ++i) {
        // for child edge keys, compute the key value
        for (int j: posChildEdge[i]) {
            childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
        }
        // for aggregation edge keys, compute the key value
        for (int j: posAggreEdge[i]) {
            aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, i, pos,
                                         aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
        }
    }
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[0] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            EdgeID e = 0;
            if (triEdgeType[mappingSize] != 0) {
                int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                if (triEdgeType[mappingSize] == 1)
                    e = startOffset[pos2] + pos[pos2] - 1;
                else if (triEdgeType[mappingSize] == 2)
                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                else if (triEdgeType[mappingSize] == 5)
                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                else {
                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                    ++gNumEdgeID;
#endif
                    e = dout.getEdgeID(key1, key2);
                }
            }
            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                               nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    pos[0] = 0;
    while (mappingSize >= (int)tau.prefixSize) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            // for child edge keys, compute the key value
            for (int j: posChildEdge[mappingSize]) {
                childKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             childEdgeType[j], dataV[childKeyPos[j][0]], dataV[childKeyPos[j][1]]);
            }
            // for aggregation edge keys, compute the key value
            for (int j: posAggreEdge[mappingSize]) {
                aggreKey[j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                             aggreEdgeType[j], dataV[aggrePos[2 * j]], dataV[aggrePos[2 * j + 1]]);
            }
            for (int i = 0; i < nodesAtStep[mappingSize].size(); ++i) {
                VertexID cID = nodesAtStep[mappingSize][i];
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    patternVs[cID][j] = patternV[prefixPos[cID][j]];
                    dataVs[cID][j] = dataV[prefixPos[cID][j]];
                    startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                    poses[cID][j] = pos[prefixPos[cID][j]];
                    visits[cID][dataVs[cID][j]] = true;
                }
                if (t.getNode(cID).keySize == 0) H[cID][0] = 0;
                else if (t.getNode(cID).keySize == 1) {
                    if (keyPosSizes[cID] < n / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * n);
                }
                else {
                    if (keyPosSizes[cID] < m / 8 + 1) {
                        for (int j = 0; j < keyPosSizes[cID]; ++j) {
                            H[cID][keyPoses[cID][j]] = 0;
                        }
                    }
                    else memset(H[cID], 0, sizeof(Count) * m);
                }
                multiJoinWrapper(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    visits[cID][dataVs[cID][j]] = false;
                }
            }
            if (mappingSize == tau.numVertices - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                Count cnt = 1;
                // multiply count in children
                for (int j = 0; j < child.size(); ++j) {
                    VertexID cID = child[j];
                    if (childKeyPos[j].empty()) cnt *= H[cID][0];
                    else if (childKeyPos[j].size() == 1) {
                        cnt *= H[cID][dataV[childKeyPos[j][0]]];
                    }
                    else {
                        cnt *= H[cID][childKey[j]];
                    }
                }
                if (isRoot) {
                    if (orbitType == 0) h[0] += cnt * aggreWeight[0];
                    else if (orbitType == 1) {
                        for (int j = 0; j < aggreV.size(); ++j) {
                            VertexID key = dataV[aggrePos[j]];
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[key] += cnt * aggreWeight[j];
                        }
                    }
                    else {
                        for (int j = 0; j < aggreKey.size(); ++j) {
#ifdef COLLECT_STATISTICS
                            ++gNumUpdate;
#endif
                            h[aggreKey[j]] += cnt * aggreWeight[j];
                        }
                    }
                }
                else {
                    if (tau.keySize == 0) h[0] += cnt;
                    else if (tau.keySize == 1) {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePos[0]]] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = dataV[aggrePos[0]];
                            ++keyPosSize;
                        }
                    }
                    else {
                        EdgeID e = aggreKey[0];
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[e] += cnt;
                        if (keyPosSize < sizeBound) {
                            keyPos[keyPosSize] = e;
                            ++keyPosSize;
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    EdgeID e = 0;
                    if (triEdgeType[mappingSize] != 0) {
                        int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                        if (triEdgeType[mappingSize] == 1)
                            e = startOffset[pos2] + pos[pos2] - 1;
                        else if (triEdgeType[mappingSize] == 2)
                            e = outID[startOffset[pos2] + pos[pos2] - 1];
                        else if (triEdgeType[mappingSize] == 5)
                            e = unID[startOffset[pos2] + pos[pos2] - 1];
                        else {
                            VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                            ++gNumEdgeID;
#endif
                            e = dout.getEdgeID(key1, key2);
                        }
                    }
                    generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                                       nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --mappingSize;
        if (mappingSize >= (int)tau.prefixSize) {
            visited[dataV[mappingSize]] = false;
        }
    }
}

void multiJoinWrapper(
        VertexID nID,
        const Tree &t,
        const std::vector<VertexID> &child,
        VertexID ***candidates,
        ui **candCounts,
        HashTable *H,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool useTriangle,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffsets,
        VertexID **patternVs,
        VertexID **dataVs,
        bool **visits,
        ui **poses,
        ui **keyPoses,
        ui *keyPosSizes,
        ui *sizeBounds,
        VertexID *&tmp,
        VertexID *allV
) {
    bool edgeKey = t.getNode(nID).edgeKey;
    if (!useTriangle) {
        if (!edgeKey) {
            multiJoin(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                      startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
        else {
            multiJoinE(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                       startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
    }
    else {
        if (!edgeKey) {
            multiJoinT(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                       startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
        else {
            multiJoinET(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                        startOffsets, patternVs, dataVs, visits, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
    }
}

void multiJoinTree(
        const Tree &t,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool useTriangle,
        const Triangle &tri,
        const Pattern &p,
        HashTable *H,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID **startOffset,
        VertexID **patternV,
        VertexID **dataV,
        bool **visited,
        VertexID *&tmp,
        VertexID *allV
) {
    int numNodes = (int)t.getNumNodes();
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    VertexID ***candidate = new VertexID **[numNodes];
    ui **candCount = new ui *[numNodes];
    ui **pos = new ui *[numNodes];
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        const Node &tau = t.getNode(nID);
        candCount[nID] = new ui[tau.numVertices];
        candidate[nID] = new VertexID *[tau.numVertices];
        pos[nID] = new ui[tau.numVertices];
        memset(pos[nID], 0, sizeof(ui) * tau.numVertices);
        const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
        for (ui i = tau.prefixSize; i < tau.numVertices; ++i) {
            if (nodeCandPos[i])
                candidate[nID][i] = new VertexID[dout.getNumVertices()];
        }
    }
    ui **keyPos = new ui *[t.getNumNodes()];
    ui *keyPosSize = new ui[t.getNumNodes()];
    memset(keyPosSize, 0, sizeof(ui) * (t.getNumNodes()));
    ui *sizeBounds = new ui[t.getNumNodes()];
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        if (t.getNode(nID).keySize == 0)
            sizeBounds[nID] = 1;
        else if (t.getNode(nID).keySize == 1)
            sizeBounds[nID] = n / 8 + 1;
        else if (t.getNode(nID).keySize == 2)
            sizeBounds[nID] = m / 8 + 1;
        keyPos[nID] = new ui[sizeBounds[nID]];
    }
    // call simple nodes
    for (VertexID nID: t.getPostOrder()) {
        if (t.getNode(nID).prefixSize == 0) {
            multiJoinWrapper(nID, t, t.getChild()[nID], candidate, candCount, H, din, dout, dun, useTriangle, tri, p,
                             outID, unID, reverseID, startOffset, patternV, dataV, visited,
                             pos, keyPos, keyPosSize, sizeBounds, tmp, allV);
        }
    }
    for (VertexID nID = 0; nID < numNodes; ++nID) {
        const Node &tau = t.getNode(nID);
        const std::vector<bool> &nodeCandPos = t.getNodeCandPos(nID);
        for (ui i = tau.prefixSize; i < tau.numVertices; ++i) {
            if (nodeCandPos[i])
                delete[] candidate[nID][i];
        }
        delete[] candidate[nID];
        delete[] candCount[nID];
        delete[] pos[nID];
        delete[] keyPos[nID];
    }
    delete[] candCount;
    delete[] candidate;
    delete[] pos;
    delete[] keyPos;
}

void executeSharedNode(
        const Node &rep,
        const std::vector<std::vector<VertexID>> &children,
        const std::vector<VertexID> &nIDs,
        const std::vector<Tree> &trees,
        std::vector<HashTable *> &hashTables,
        const std::vector<std::vector<int>> &aggreWeights,
        const std::vector<int> &id2AggreKey,
        const std::vector<int> &id2ChildKey,
        const std::vector<bool> &isRoot,
        const std::vector<std::vector<std::vector<int>>> &childKeyPoses,
        const std::vector<std::vector<int>> &aggrePoses,
        const std::vector<int> &tableIDs,
        std::vector<ui *> &keyPos,
        std::vector<ui> &keyPosSize,
        const std::vector<ui> &sizeBound,
        const std::vector<bool> &nodeInterPos,
        const std::vector<std::vector<int>> &nodeInPos,
        const std::vector<std::vector<int>> &nodeOutPos,
        const std::vector<std::vector<int>> &nodeUnPos,
        const std::vector<std::vector<int>> &greaterPos,
        const std::vector<std::vector<int>> &lessPos,
        const std::vector<int> &nonLeafIDs,
        VertexID **candidate,
        ui *candCount,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    int orbitType = trees[0].getOrbitType();
    int numSharedNodes = (int)nIDs.size();
    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                              nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = rep.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            if (depth == rep.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                for (int id = 0; id < numSharedNodes; ++id) {
                    Count cnt = 1;
                    HashTable *H = hashTables[id];
                    int ck = id2ChildKey[id];
                    int ak = id2AggreKey[id];
                    // multiply count in children
                    for (int j = 0; j < children[id].size(); ++j) {
                        VertexID cID = children[id][j];
                        cnt *= H[cID][dataV[childKeyPoses[ck][j][0]]];
                    }
                    // update hash table
                    HashTable h = H[nIDs[id]];
                    if (isRoot[id]) {
                        if (orbitType == 0) h[0] += cnt * aggreWeights[id][0];
                        else if (orbitType == 1) {
                            for (int j = 0; j < aggrePoses[ak].size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                h[dataV[aggrePoses[ak][j]]] += cnt * aggreWeights[id][j];
                            }
                        }
                    }
                    else {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePoses[ak][0]]] += cnt;
                        int tableID = tableIDs[id];
                        if (keyPosSize[tableID] < sizeBound[tableID]) {
                            keyPos[tableID][keyPosSize[tableID]] = dataV[aggrePoses[ak][0]];
                            ++keyPosSize[tableID];
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                                      nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
}

void executeSharedNodeT(
        const Node &rep,
        const std::vector<std::vector<VertexID>> &children,
        const std::vector<VertexID> &nIDs,
        const std::vector<Tree> &trees,
        std::vector<HashTable *> &hashTables,
        const std::vector<std::vector<int>> &aggreWeights,
        const std::vector<int> &id2AggreKey,
        const std::vector<int> &id2ChildKey,
        const std::vector<bool> &isRoot,
        const std::vector<std::vector<std::vector<int>>> &childKeyPoses,
        const std::vector<std::vector<int>> &aggrePoses,
        const std::vector<int> &tableIDs,
        std::vector<ui *> &keyPos,
        std::vector<ui> &keyPosSize,
        const std::vector<ui> &sizeBound,
        const std::vector<bool> &nodeInterPos,
        const std::vector<std::vector<int>> &nodeInPos,
        const std::vector<std::vector<int>> &nodeOutPos,
        const std::vector<std::vector<int>> &nodeUnPos,
        const std::vector<bool> &nodeCandPos,
        const std::vector<std::vector<int>> &greaterPos,
        const std::vector<std::vector<int>> &lessPos,
        const std::vector<std::pair<int, int>> &nodeTriPos,
        const std::vector<int> &triEdgeType,
        const std::vector<int> &triEndType,
        const std::vector<int> &nonLeafIDs,
        VertexID **candidate,
        ui *candCount,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    int orbitType = trees[0].getOrbitType();
    int numSharedNodes = (int)nIDs.size();
    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            EdgeID e = 0;
            if (triEdgeType[mappingSize] != 0) {
                int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                if (triEdgeType[mappingSize] == 1)
                    e = startOffset[pos2] + pos[pos2] - 1;
                else if (triEdgeType[mappingSize] == 2)
                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                else if (triEdgeType[mappingSize] == 5)
                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                else {
                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                    ++gNumEdgeID;
#endif
                    e = dout.getEdgeID(key1, key2);
                }
            }
            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                               nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = rep.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            if (depth == rep.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                for (int id = 0; id < numSharedNodes; ++id) {
                    Count cnt = 1;
                    HashTable *H = hashTables[id];
                    int ck = id2ChildKey[id];
                    int ak = id2AggreKey[id];
                    // multiply count in children
                    for (int j = 0; j < children[id].size(); ++j) {
                        VertexID cID = children[id][j];
                        cnt *= H[cID][dataV[childKeyPoses[ck][j][0]]];
                    }
                    // update hash table
                    HashTable h = H[nIDs[id]];
                    if (isRoot[id]) {
                        if (orbitType == 0) h[0] += cnt * aggreWeights[id][0];
                        else if (orbitType == 1) {
                            for (int j = 0; j < aggrePoses[ak].size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                h[dataV[aggrePoses[ak][j]]] += cnt * aggreWeights[id][j];
                            }
                        }
                    }
                    else {
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[dataV[aggrePoses[ak][0]]] += cnt;
                        int tableID = tableIDs[id];
                        if (keyPosSize[tableID] < sizeBound[tableID]) {
                            keyPos[tableID][keyPosSize[tableID]] = dataV[aggrePoses[ak][0]];
                            ++keyPosSize[tableID];
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    EdgeID e = 0;
                    if (triEdgeType[mappingSize] != 0) {
                        int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                        if (triEdgeType[mappingSize] == 1)
                            e = startOffset[pos2] + pos[pos2] - 1;
                        else if (triEdgeType[mappingSize] == 2)
                            e = outID[startOffset[pos2] + pos[pos2] - 1];
                        else if (triEdgeType[mappingSize] == 5)
                            e = unID[startOffset[pos2] + pos[pos2] - 1];
                        else{
                            VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                            ++gNumEdgeID;
#endif
                            e = dout.getEdgeID(key1, key2);
                        }
                    }
                    generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                                       nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
}

void executeSharedNodeEdgeKey(
        const Node &rep,
        const std::vector<std::vector<VertexID>> &children,
        const std::vector<VertexID> &nIDs,
        const std::vector<Tree> &trees,
        std::vector<HashTable *> &hashTables,
        const std::vector<std::vector<int>> &aggreWeights,
        const std::vector<int> &id2AggreKey,
        const std::vector<int> &id2ChildKey,
        const std::vector<bool> &isRoot,
        const std::vector<std::vector<std::vector<int>>> &childKeyPoses,
        const std::vector<std::vector<int>> &aggrePoses,
        const std::vector<int> &tableIDs,
        std::vector<ui *> &keyPos,
        std::vector<ui> &keyPosSize,
        const std::vector<ui> &sizeBound,
        const std::vector<bool> &nodeInterPos,
        const std::vector<std::vector<int>> &nodeInPos,
        const std::vector<std::vector<int>> &nodeOutPos,
        const std::vector<std::vector<int>> &nodeUnPos,
        const std::vector<std::vector<int>> &greaterPos,
        const std::vector<std::vector<int>> &lessPos,
        const std::vector<std::vector<std::pair<int, int>>> &posChildEdge,
        const std::vector<std::vector<std::pair<int, int>>> &posAggreEdge,
        const std::vector<std::vector<int>> &childEdgeType,
        const std::vector<std::vector<int>> &aggreEdgeType,
        const std::vector<int> &nonLeafIDs,
        VertexID **candidate,
        ui *candCount,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    int orbitType = trees[0].getOrbitType();
    ui **childKey = new ui *[childKeyPoses.size()];
    for (int group = 0; group < childKeyPoses.size(); ++group) {
        childKey[group] = new ui[childKeyPoses[group].size()];
    }
    ui **aggreKey = new ui *[aggrePoses.size()];
    for (int group = 0; group < aggrePoses.size(); ++group) {
        if (!isRoot[group]) aggreKey[group] = new ui[1];
        else if (orbitType == 1) aggreKey[group] = new ui[aggrePoses[group].size()];
        else aggreKey[group] = new ui[aggrePoses[group].size() / 2];
    }
    int numSharedNodes = (int)nIDs.size();
    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                              nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    for (int k = 1; k < mappingSize; ++k) {
        // for child edge keys, compute the key value
        for (const std::pair<int, int> &pr: posChildEdge[k]) {
            int i = pr.first, j = pr.second;
            childKey[i][j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, k, pos,
                                            childEdgeType[i][j], dataV[childKeyPoses[i][j][0]], dataV[childKeyPoses[i][j][1]]);
        }
        // for aggregation edge keys, compute the key value
        for (const std::pair<int, int> &pr: posAggreEdge[k]) {
            int i = pr.first, j = pr.second;
            aggreKey[i][j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, k, pos,
                                            aggreEdgeType[i][j], dataV[aggrePoses[i][2 * j]], dataV[aggrePoses[i][2 * j + 1]]);
        }
    }
    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = rep.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            // for child edge keys, compute the key value
            for (const std::pair<int, int> &pr: posChildEdge[mappingSize]) {
                int i = pr.first, j = pr.second;
                childKey[i][j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                                childEdgeType[i][j], dataV[childKeyPoses[i][j][0]], dataV[childKeyPoses[i][j][1]]);
            }
            // for aggregation edge keys, compute the key value
            for (const std::pair<int, int> &pr: posAggreEdge[mappingSize]) {
                int i = pr.first, j = pr.second;
                aggreKey[i][j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, mappingSize, pos,
                                                aggreEdgeType[i][j], dataV[aggrePoses[i][2 * j]], dataV[aggrePoses[i][2 * j + 1]]);
            }

            if (depth == rep.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // multiply count in children
                for (int id = 0; id < numSharedNodes; ++id) {
                    Count cnt = 1;
                    HashTable *H = hashTables[id];
                    HashTable h = H[nIDs[id]];
                    int ck = id2ChildKey[id];
                    int ak = id2AggreKey[id];
                    for (int j = 0; j < children[id].size(); ++j) {
                        VertexID cID = children[id][j];
                        if (childKeyPoses[ck][j].size() == 1)
                            cnt *= H[cID][dataV[childKeyPoses[ck][j][0]]];
                        else
                            cnt *= H[cID][childKey[ck][j]];
                    }
                    if (isRoot[id]) {
                        if (orbitType == 0) h[0] += cnt * aggreWeights[id][0];
                        else if (orbitType == 1) {
                            for (int j = 0; j < aggrePoses[ak].size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                h[dataV[aggrePoses[ak][j]]] += cnt * aggreWeights[id][j];
                            }
                        }
                        else {
                            for (int j = 0; j < aggrePoses[ak].size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                h[aggreKey[ak][j]] += cnt * aggreWeights[id][j];
                            }
                        }
                    }
                    else {
                        ui key = 0;
                        if (aggrePoses[ak].size() == 1) key = dataV[aggrePoses[ak][0]];
                        else key = aggreKey[ak][0];
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[key] += cnt;
                        int tableID = tableIDs[id];
                        if (keyPosSize[tableID] < sizeBound[tableID]) {
                            keyPos[tableID][keyPosSize[tableID]] = key;
                            ++keyPosSize[tableID];
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    generateCandidate(din, dout, dun, dataV, mappingSize, nodeInPos[mappingSize],
                                      nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
    for (int group = 0; group < childKeyPoses.size(); ++group) {
        delete[] childKey[group];
    }
    for (int group = 0; group < aggrePoses.size(); ++group) {
        delete[] aggreKey[group];
    }
    delete[] childKey;
    delete[] aggreKey;
}

void executeSharedNodeEdgeKeyT(
        const Node &rep,
        const std::vector<std::vector<VertexID>> &children,
        const std::vector<VertexID> &nIDs,
        const std::vector<Tree> &trees,
        std::vector<HashTable *> &hashTables,
        const std::vector<std::vector<int>> &aggreWeights,
        const std::vector<int> &id2AggreKey,
        const std::vector<int> &id2ChildKey,
        const std::vector<bool> &isRoot,
        const std::vector<std::vector<std::vector<int>>> &childKeyPoses,
        const std::vector<std::vector<int>> &aggrePoses,
        const std::vector<int> &tableIDs,
        std::vector<ui *> &keyPos,
        std::vector<ui> &keyPosSize,
        const std::vector<ui> &sizeBound,
        const std::vector<bool> &nodeInterPos,
        const std::vector<std::vector<int>> &nodeInPos,
        const std::vector<std::vector<int>> &nodeOutPos,
        const std::vector<std::vector<int>> &nodeUnPos,
        const std::vector<bool> &nodeCandPos,
        const std::vector<std::vector<int>> &greaterPos,
        const std::vector<std::vector<int>> &lessPos,
        const std::vector<std::vector<std::pair<int, int>>> &posChildEdge,
        const std::vector<std::vector<std::pair<int, int>>> &posAggreEdge,
        const std::vector<std::vector<int>> &childEdgeType,
        const std::vector<std::vector<int>> &aggreEdgeType,
        const std::vector<std::pair<int, int>> &nodeTriPos,
        const std::vector<int> &nonLeafIDs,
        const std::vector<int> &triEdgeType,
        const std::vector<int> &triEndType,
        VertexID **candidate,
        ui *candCount,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    int orbitType = trees[0].getOrbitType();
    ui **childKey = new ui *[childKeyPoses.size()];
    for (int group = 0; group < childKeyPoses.size(); ++group) {
        childKey[group] = new ui[childKeyPoses[group].size()];
    }
    ui **aggreKey = new ui *[aggrePoses.size()];
    for (int group = 0; group < aggrePoses.size(); ++group) {
        if (!isRoot[group]) aggreKey[group] = new ui[1];
        else if (orbitType == 1) aggreKey[group] = new ui[aggrePoses[group].size()];
        else aggreKey[group] = new ui[aggrePoses[group].size() / 2];
    }
    int numSharedNodes = (int)nIDs.size();
    ui n = din.getNumVertices();
    if (mappingSize == 0) {
        candidate[0] = allV;
        candCount[0] = n;
        pos[mappingSize] = 0;
    }
    else {
        if (!nodeInterPos[mappingSize]) {
            if (!nodeOutPos[mappingSize].empty()) {
                VertexID v = dataV[nodeOutPos[mappingSize][0]];
                startOffset[mappingSize] = outOffset[v];
                candidate[mappingSize] = outNbors + outOffset[v];
                candCount[mappingSize] = outOffset[v + 1] - outOffset[v];
            }
            else if (!nodeInPos[mappingSize].empty()){
                VertexID v = dataV[nodeInPos[mappingSize][0]];
                startOffset[mappingSize] = inOffset[v];
                candidate[mappingSize] = inNbors + inOffset[v];
                candCount[mappingSize] = inOffset[v + 1] - inOffset[v];
            }
            else {
                VertexID v = dataV[nodeUnPos[mappingSize][0]];
                startOffset[mappingSize] = unOffset[v];
                candidate[mappingSize] = unNbors + unOffset[v];
                candCount[mappingSize] = unOffset[v + 1] - unOffset[v];
            }
        }
        else {
            EdgeID e = 0;
            if (triEdgeType[mappingSize] != 0) {
                int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                if (triEdgeType[mappingSize] == 1)
                    e = startOffset[pos2] + pos[pos2] - 1;
                else if (triEdgeType[mappingSize] == 2)
                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                else if (triEdgeType[mappingSize] == 5)
                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                else {
                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                    ++gNumEdgeID;
#endif
                    e = dout.getEdgeID(key1, key2);
                }
            }
            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                               nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
        }
        pos[mappingSize] = 0;
        // apply the symmetry breaking rules
        if (!greaterPos[mappingSize].empty()) {
            ui maxTarget = dataV[greaterPos[mappingSize][0]];
            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                    maxTarget = dataV[greaterPos[mappingSize][i]];
            }
            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
        }
        if (!lessPos[mappingSize].empty()) {
            ui minTarget = dataV[lessPos[mappingSize][0]];
            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                if (dataV[lessPos[mappingSize][i]] < minTarget)
                    minTarget = dataV[lessPos[mappingSize][i]];
            }
            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
        }
    }
    for (int k = 1; k < mappingSize; ++k) {
        // for child edge keys, compute the key value
        for (const std::pair<int, int> &pr: posChildEdge[k]) {
            int i = pr.first, j = pr.second;
            childKey[i][j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, k, pos,
                                            childEdgeType[i][j], dataV[childKeyPoses[i][j][0]], dataV[childKeyPoses[i][j][1]]);
        }
        // for aggregation edge keys, compute the key value
        for (const std::pair<int, int> &pr: posAggreEdge[k]) {
            int i = pr.first, j = pr.second;
            aggreKey[i][j] = computeEdgeKey(din, dout, dun, outID, unID, reverseID, startOffset, k, pos,
                                            aggreEdgeType[i][j], dataV[aggrePoses[i][2 * j]], dataV[aggrePoses[i][2 * j + 1]]);
        }
    }
    int depth = 0;
    while (depth >= 0) {
        while (pos[mappingSize] < candCount[mappingSize]) {
            VertexID v = candidate[mappingSize][pos[mappingSize]];
            ++pos[mappingSize];
            if (visited[v]) continue;
            visited[v] = true;
            patternV[mappingSize] = rep.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            // for child edge keys, compute the key value
            for (const std::pair<int, int> &pr: posChildEdge[mappingSize]) {
                int i = pr.first, j = pr.second;
                switch (childEdgeType[i][j]) {
                    case 1:
                        childKey[i][j] = startOffset[mappingSize] + pos[mappingSize] - 1;
                        break;
                    case 2:
                        childKey[i][j] = outID[startOffset[mappingSize] + pos[mappingSize] - 1];
                        break;
                    case 3:
#ifdef COLLECT_STATISTICS
                        ++gNumEdgeID;
#endif
                        childKey[i][j] = dout.getEdgeID(dataV[childKeyPoses[i][j][0]], dataV[childKeyPoses[i][j][1]]);
                        break;
                    case 4:
#ifdef COLLECT_STATISTICS
                        ++gNumEdgeID;
#endif
                        childKey[i][j] = dun.getUndirectedEID(dataV[childKeyPoses[i][j][0]], dataV[childKeyPoses[i][j][1]]);
                        break;
                    case 6:
                        childKey[i][j] = reverseID[startOffset[mappingSize] + pos[mappingSize] - 1];
                        break;
                    case 7:
#ifdef COLLECT_STATISTICS
                        ++gNumEdgeID;
#endif
                        childKey[i][j] = dun.getUndirectedEID(dataV[childKeyPoses[i][j][1]], dataV[childKeyPoses[i][j][0]]);
                        break;
                }
            }
            // for aggregation edge keys, compute the key value
            for (const std::pair<int, int> &pr: posAggreEdge[mappingSize]) {
                int i = pr.first, j = pr.second;
                switch (aggreEdgeType[i][j]) {
                    case 1:
                        aggreKey[i][j] = startOffset[mappingSize] + pos[mappingSize] - 1;
                        break;
                    case 2:
                        aggreKey[i][j] = outID[startOffset[mappingSize] + pos[mappingSize] - 1];
                        break;
                    case 3:
#ifdef COLLECT_STATISTICS
                        ++gNumEdgeID;
#endif
                        aggreKey[i][j] = dout.getEdgeID(dataV[aggrePoses[i][2 * j]], dataV[aggrePoses[i][2 * j + 1]]);
                        break;
                    case 4:
#ifdef COLLECT_STATISTICS
                        ++gNumEdgeID;
#endif
                        aggreKey[i][j] = dun.getUndirectedEID(dataV[aggrePoses[i][2 * j]], dataV[aggrePoses[i][2 * j + 1]]);
                        break;
                    case 5:
                        aggreKey[i][j] = unID[startOffset[mappingSize] + pos[mappingSize] - 1];
                        break;
                    case 6:
                        aggreKey[i][j] = reverseID[startOffset[mappingSize] + pos[mappingSize] - 1];
                        break;
                    case 7:
#ifdef COLLECT_STATISTICS
                        ++gNumEdgeID;
#endif
                        aggreKey[i][j] = dun.getUndirectedEID(dataV[aggrePoses[i][2 * j + 1]], dataV[aggrePoses[i][2 * j]]);
                        break;
                }
            }

            if (depth == rep.localOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // multiply count in children
                for (int id = 0; id < numSharedNodes; ++id) {
                    Count cnt = 1;
                    HashTable *H = hashTables[id];
                    HashTable h = H[nIDs[id]];
                    int ck = id2ChildKey[id];
                    int ak = id2AggreKey[id];
                    for (int j = 0; j < children[id].size(); ++j) {
                        VertexID cID = children[id][j];
                        if (childKeyPoses[ck][j].size() == 1)
                            cnt *= H[cID][dataV[childKeyPoses[ck][j][0]]];
                        else
                            cnt *= H[cID][childKey[ck][j]];
                    }
                    if (isRoot[id]) {
                        if (orbitType == 0) h[0] += cnt * aggreWeights[id][0];
                        else if (orbitType == 1) {
                            for (int j = 0; j < aggrePoses[ak].size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                h[dataV[aggrePoses[ak][j]]] += cnt * aggreWeights[id][j];
                            }
                        }
                        else {
                            for (int j = 0; j < aggrePoses[ak].size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                h[aggreKey[ak][j]] += cnt * aggreWeights[id][j];
                            }
                        }
                    }
                    else {
                        ui key = 0;
                        if (aggrePoses[ak].size() == 1) key = dataV[aggrePoses[ak][0]];
                        else key = aggreKey[ak][0];
#ifdef COLLECT_STATISTICS
                        ++gNumUpdate;
#endif
                        h[key] += cnt;
                        int tableID = tableIDs[id];
                        if (keyPosSize[tableID] < sizeBound[tableID]) {
                            keyPos[tableID][keyPosSize[tableID]] = key;
                            ++keyPosSize[tableID];
                        }
                    }
                }
                visited[dataV[mappingSize]] = false;
            }
            else {
#ifdef COLLECT_STATISTICS
                ++gNumIntermediate;
#endif
                ++mappingSize;
                ++depth;
                pos[mappingSize] = 0;
                if (!nodeInterPos[mappingSize]) {
                    if (!nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    EdgeID e = 0;
                    if (triEdgeType[mappingSize] != 0) {
                        int pos1 = nodeTriPos[mappingSize].first, pos2 = nodeTriPos[mappingSize].second;
                        if (triEdgeType[mappingSize] == 1)
                            e = startOffset[pos2] + pos[pos2] - 1;
                        else if (triEdgeType[mappingSize] == 2)
                            e = outID[startOffset[pos2] + pos[pos2] - 1];
                        else if (triEdgeType[mappingSize] == 5)
                            e = unID[startOffset[pos2] + pos[pos2] - 1];
                        else {
                            VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                            ++gNumEdgeID;
#endif
                            e = dout.getEdgeID(key1, key2);
                        }
                    }
                    generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !nodeCandPos[mappingSize], e, triEndType[mappingSize],
                                       nodeInPos[mappingSize], nodeOutPos[mappingSize], nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[greaterPos[mappingSize][0]];
                    for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                        if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[lessPos[mappingSize][0]];
                    for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                        if (dataV[lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
            }
        }
        --depth;
        if (depth >= 0) {
            --mappingSize;
            visited[dataV[mappingSize]] = false;
        }
    }
    for (int group = 0; group < childKeyPoses.size(); ++group) {
        delete[] childKey[group];
    }
    for (int group = 0; group < aggrePoses.size(); ++group) {
        delete[] aggreKey[group];
    }
    delete[] childKey;
    delete[] aggreKey;
//    delete[] cnt;
}

void executeConNode(
        const ConNode &cn,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool useTriangle,
        const Triangle &tri,
        const Pattern &p,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    int orbitType = cn.orbitType;
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    Count *cnt, *cntCopy;
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    if (cn.childKeyPos.size() == 1) {
        cnt = new Count[n];
        cntCopy = new Count[n];
        memset(cnt, 0, sizeof(Count) * n);
        memset(cntCopy, 0, sizeof(Count) * n);
    }
    else {
        if (p.u.isEOrbitDir()) {
            cnt = new Count[m];
            cntCopy = new Count[m];
            memset(cnt, 0, sizeof(Count) * m);
            memset(cntCopy, 0, sizeof(Count) * m);
        }
        else {
            cnt = new Count[m * 2];
            cntCopy = new Count[m * 2];
            memset(cnt, 0, sizeof(Count) * m * 2);
            memset(cntCopy, 0, sizeof(Count) * m * 2);
        }
    }
    int mappingSize = 0;
    int prefixSize = (int)cn.prefixOrder.size();
    VertexID **candidate = new VertexID *[cn.nodeOrder.size()];
    ui *candCount = new ui[cn.nodeOrder.size()];
    for (int i = 0; i < cn.nodeOrder.size(); ++i) {
        if (cn.candPos[i])
            candidate[i] = new VertexID[n];
    }
    candidate[0] = allV;
    candCount[0] = n;
    pos[0] = 0;
    if (!cn.edgeKey) {
        while (mappingSize >= 0) {
            // matching the prefix
            if (mappingSize < prefixSize) {
                while (pos[mappingSize] < candCount[mappingSize]) {
                    VertexID v = candidate[mappingSize][pos[mappingSize]];
                    ++pos[mappingSize];
                    if (visited[v]) continue;
                    visited[v] = true;
                    patternV[mappingSize] = cn.nodeOrder[mappingSize];
                    dataV[mappingSize] = v;
                    ++mappingSize;
                    pos[mappingSize] = 0;
                    if (!cn.interPos[mappingSize]) {
                        if (!cn.nodeOutPos[mappingSize].empty()) {
                            VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                            startOffset[mappingSize] = outOffset[w];
                            candidate[mappingSize] = outNbors + outOffset[w];
                            candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                        }
                        else if (!cn.nodeInPos[mappingSize].empty()){
                            VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                            startOffset[mappingSize] = inOffset[w];
                            candidate[mappingSize] = inNbors + inOffset[w];
                            candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                        }
                        else {
                            VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                            startOffset[mappingSize] = unOffset[w];
                            candidate[mappingSize] = unNbors + unOffset[w];
                            candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                        }
                    }
                    else {
                        EdgeID e = 0;
                        if (cn.triEdgeType[mappingSize] != 0) {
                            int pos1 = cn.triPos[mappingSize].first, pos2 = cn.triPos[mappingSize].second;
                            if (cn.triEdgeType[mappingSize] == 1)
                                e = startOffset[pos2] + pos[pos2] - 1;
                            else if (cn.triEdgeType[mappingSize] == 2)
                                e = outID[startOffset[pos2] + pos[pos2] - 1];
                            else if (cn.triEdgeType[mappingSize] == 5)
                                e = unID[startOffset[pos2] + pos[pos2] - 1];
                            else {
                                VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                                ++gNumEdgeID;
#endif
                                e = dout.getEdgeID(key1, key2);
                            }
                        }
                        generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !cn.candPos[mappingSize], e, cn.triEndType[mappingSize],
                                           cn.nodeInPos[mappingSize], cn.nodeOutPos[mappingSize], cn.nodeUnPos[mappingSize], candidate, candCount, tmp);
                    }
                    // apply the symmetry breaking rules
                    if (!cn.greaterPos[mappingSize].empty()) {
                        ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                        for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                            if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                        }
                        pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                    }
                    if (!cn.lessPos[mappingSize].empty()) {
                        ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                        for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                            if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                minTarget = dataV[cn.lessPos[mappingSize][i]];
                        }
                        candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                    }
                    if (mappingSize == prefixSize) break;
                }
                if (mappingSize != prefixSize) {
                    --mappingSize;
                    if (mappingSize >= 0)
                        visited[dataV[mappingSize]] = false;
                }
            }
            else {
                // matching the remaining vertices to compute the child table
                while (mappingSize >= prefixSize) {
                    while (pos[mappingSize] < candCount[mappingSize]) {
                        VertexID v = candidate[mappingSize][pos[mappingSize]];
                        ++pos[mappingSize];
                        if (visited[v]) continue;
                        visited[v] = true;
                        patternV[mappingSize] = cn.nodeOrder[mappingSize];
                        dataV[mappingSize] = v;
                        if (mappingSize == cn.nodeOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                            ++gNumMatch;
#endif
                            ++cnt[dataV[cn.childKeyPos[0]]];
                            cntCopy[dataV[cn.childKeyPos[0]]] = cnt[dataV[cn.childKeyPos[0]]];
                            visited[dataV[mappingSize]] = false;
                        }
                        else {
#ifdef COLLECT_STATISTICS
                            ++gNumIntermediate;
#endif
                            ++mappingSize;
                            pos[mappingSize] = 0;
                            if (!cn.interPos[mappingSize]) {
                                if (!cn.nodeOutPos[mappingSize].empty()) {
                                    VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                    startOffset[mappingSize] = outOffset[w];
                                    candidate[mappingSize] = outNbors + outOffset[w];
                                    candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                                }
                                else if (!cn.nodeInPos[mappingSize].empty()){
                                    VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                    startOffset[mappingSize] = inOffset[w];
                                    candidate[mappingSize] = inNbors + inOffset[w];
                                    candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                                }
                                else {
                                    VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                    startOffset[mappingSize] = unOffset[w];
                                    candidate[mappingSize] = unNbors + unOffset[w];
                                    candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                                }
                            }
                            else {
                                EdgeID e = 0;
                                if (cn.triEdgeType[mappingSize] != 0) {
                                    int pos1 = cn.triPos[mappingSize].first, pos2 = cn.triPos[mappingSize].second;
                                    if (cn.triEdgeType[mappingSize] == 1)
                                        e = startOffset[pos2] + pos[pos2] - 1;
                                    else if (cn.triEdgeType[mappingSize] == 2)
                                        e = outID[startOffset[pos2] + pos[pos2] - 1];
                                    else if (cn.triEdgeType[mappingSize] == 5)
                                        e = unID[startOffset[pos2] + pos[pos2] - 1];
                                    else {
                                        VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                                        ++gNumEdgeID;
#endif
                                        e = dout.getEdgeID(key1, key2);
                                    }
                                }
                                generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !cn.candPos[mappingSize], e, cn.triEndType[mappingSize],
                                                   cn.nodeInPos[mappingSize], cn.nodeOutPos[mappingSize], cn.nodeUnPos[mappingSize], candidate, candCount, tmp);
                            }
                            // apply the symmetry breaking rules
                            if (!cn.greaterPos[mappingSize].empty()) {
                                ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                                for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                        maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                                }
                                pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                            }
                            if (!cn.lessPos[mappingSize].empty()) {
                                ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                                for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                        minTarget = dataV[cn.lessPos[mappingSize][i]];
                                }
                                candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                            }
                        }
                    }
                    --mappingSize;
                    if (mappingSize >= prefixSize)
                        visited[dataV[mappingSize]] = false;
                }
                // matching the remaining vertices to query the child table and compute orbit count
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!cn.greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                    for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                        if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!cn.lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                    for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                        if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[cn.lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
                while (mappingSize >= prefixSize) {
                    while (pos[mappingSize] < candCount[mappingSize]) {
                        VertexID v = candidate[mappingSize][pos[mappingSize]];
                        ++pos[mappingSize];
                        if (visited[v]) continue;
                        visited[v] = true;
                        patternV[mappingSize] = cn.nodeOrder[mappingSize];
                        dataV[mappingSize] = v;
                        if (mappingSize == cn.nodeOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                            ++gNumMatch;
#endif
                            ui childKey = dataV[cn.childKeyPos[0]];
                            if (orbitType == 0) {
                                cn.hashTables[0][0] += choosec(cnt[childKey], cn.num);
                            }
                            else {
                                for (int j = 0; j < cn.aggrePos.size(); ++j) {
#ifdef COLLECT_STATISTICS
                                    ++gNumUpdate;
#endif
                                    ui aggreKey = dataV[cn.aggrePos[j]];
                                    if (cn.aggreType[j])
                                        cn.hashTables[j][aggreKey] += choosec(cnt[childKey], cn.num);
                                    else
                                        cn.hashTables[j][aggreKey] += choosec(cntCopy[childKey] - 1, cn.num - 1);
                                }
                            }
                            cnt[childKey] = 0;
                            visited[dataV[mappingSize]] = false;
                        }
                        else {
#ifdef COLLECT_STATISTICS
                            ++gNumIntermediate;
#endif
                            ++mappingSize;
                            pos[mappingSize] = 0;
                            if (!cn.interPos[mappingSize]) {
                                if (!cn.nodeOutPos[mappingSize].empty()) {
                                    VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                    startOffset[mappingSize] = outOffset[w];
                                    candidate[mappingSize] = outNbors + outOffset[w];
                                    candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                                }
                                else if (!cn.nodeInPos[mappingSize].empty()){
                                    VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                    startOffset[mappingSize] = inOffset[w];
                                    candidate[mappingSize] = inNbors + inOffset[w];
                                    candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                                }
                                else {
                                    VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                    startOffset[mappingSize] = unOffset[w];
                                    candidate[mappingSize] = unNbors + unOffset[w];
                                    candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                                }
                            }
                            else {
                                EdgeID e = 0;
                                if (cn.triEdgeType[mappingSize] != 0) {
                                    int pos1 = cn.triPos[mappingSize].first, pos2 = cn.triPos[mappingSize].second;
                                    if (cn.triEdgeType[mappingSize] == 1)
                                        e = startOffset[pos2] + pos[pos2] - 1;
                                    else if (cn.triEdgeType[mappingSize] == 2)
                                        e = outID[startOffset[pos2] + pos[pos2] - 1];
                                    else if (cn.triEdgeType[mappingSize] == 5)
                                        e = unID[startOffset[pos2] + pos[pos2] - 1];
                                    else {
                                        VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                                        ++gNumEdgeID;
#endif
                                        e = dout.getEdgeID(key1, key2);
                                    }
                                }
                                generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !cn.candPos[mappingSize], e, cn.triEndType[mappingSize],
                                                   cn.nodeInPos[mappingSize], cn.nodeOutPos[mappingSize], cn.nodeUnPos[mappingSize], candidate, candCount, tmp);

                            }
                            // apply the symmetry breaking rules
                            if (!cn.greaterPos[mappingSize].empty()) {
                                ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                                for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                        maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                                }
                                pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                            }
                            if (!cn.lessPos[mappingSize].empty()) {
                                ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                                for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                    if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                        minTarget = dataV[cn.lessPos[mappingSize][i]];
                                }
                                candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                            }
                        }
                    }
                    --mappingSize;
                    if (mappingSize >= 0)
                        visited[dataV[mappingSize]] = false;
                }
            }
        }
        delete[] cnt;
        delete[] cntCopy;
        return;
    }
    // the edge key case
    EdgeID childKey = 0;
    std::vector<EdgeID> aggreKey;
    if (orbitType == 2) aggreKey = std::vector<EdgeID>(cn.aggreEdgeType.size());
    while (mappingSize >= 0) {
        // matching the prefix
        if (mappingSize < prefixSize) {
            while (pos[mappingSize] < candCount[mappingSize]) {
                VertexID v = candidate[mappingSize][pos[mappingSize]];
                ++pos[mappingSize];
                if (visited[v]) continue;
                visited[v] = true;
                patternV[mappingSize] = cn.nodeOrder[mappingSize];
                dataV[mappingSize] = v;
                ++mappingSize;
                pos[mappingSize] = 0;
                if (!cn.interPos[mappingSize]) {
                    if (!cn.nodeOutPos[mappingSize].empty()) {
                        VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                        startOffset[mappingSize] = outOffset[w];
                        candidate[mappingSize] = outNbors + outOffset[w];
                        candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                    }
                    else if (!cn.nodeInPos[mappingSize].empty()){
                        VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                        startOffset[mappingSize] = inOffset[w];
                        candidate[mappingSize] = inNbors + inOffset[w];
                        candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                    }
                    else {
                        VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                        startOffset[mappingSize] = unOffset[w];
                        candidate[mappingSize] = unNbors + unOffset[w];
                        candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                    }
                }
                else {
                    EdgeID e = 0;
                    if (cn.triEdgeType[mappingSize] != 0) {
                        int pos1 = cn.triPos[mappingSize].first, pos2 = cn.triPos[mappingSize].second;
                        if (cn.triEdgeType[mappingSize] == 1)
                            e = startOffset[pos2] + pos[pos2] - 1;
                        else if (cn.triEdgeType[mappingSize] == 2)
                            e = outID[startOffset[pos2] + pos[pos2] - 1];
                        else if (cn.triEdgeType[mappingSize] == 5)
                            e = unID[startOffset[pos2] + pos[pos2] - 1];
                        else {
                            VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                            ++gNumEdgeID;
#endif
                            e = dout.getEdgeID(key1, key2);
                        }
                    }
                    generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !cn.candPos[mappingSize], e, cn.triEndType[mappingSize],
                                       cn.nodeInPos[mappingSize], cn.nodeOutPos[mappingSize], cn.nodeUnPos[mappingSize], candidate, candCount, tmp);
                }
                // apply the symmetry breaking rules
                if (!cn.greaterPos[mappingSize].empty()) {
                    ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                    for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                        if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                            maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                    }
                    pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                }
                if (!cn.lessPos[mappingSize].empty()) {
                    ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                    for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                        if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                            minTarget = dataV[cn.lessPos[mappingSize][i]];
                    }
                    candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                }
                if (mappingSize == prefixSize) break;
            }
            if (mappingSize != prefixSize) {
                --mappingSize;
                if (mappingSize >= 0)
                    visited[dataV[mappingSize]] = false;
            }
        }
        else {
            // matching the remaining vertices to compute the child table
            while (mappingSize >= prefixSize) {
                while (pos[mappingSize] < candCount[mappingSize]) {
                    VertexID v = candidate[mappingSize][pos[mappingSize]];
                    ++pos[mappingSize];
                    if (visited[v]) continue;
                    visited[v] = true;
                    patternV[mappingSize] = cn.nodeOrder[mappingSize];
                    dataV[mappingSize] = v;
                    if (mappingSize == cn.nodeOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                        ++gNumMatch;
#endif
                        if (cn.childEdgeType == 0) childKey = dataV[cn.childKeyPos[0]];
                        else {
                            if (cn.childEdgeType == 1) childKey = startOffset[cn.posChildEdge] + pos[cn.posChildEdge] - 1;
                            else if (cn.childEdgeType == 2) childKey = outID[startOffset[cn.posChildEdge] + pos[cn.posChildEdge] - 1];
                            else {
                                VertexID key1 = dataV[cn.childKeyPos[0]], key2 = dataV[cn.childKeyPos[1]];
                                if (cn.childEdgeType == 3) childKey = dout.getEdgeID(key1, key2);
                                else childKey = dun.getUndirectedEID(key1, key2);
#ifdef COLLECT_STATISTICS
                                ++gNumEdgeID;
#endif
                            }
                        }

                        ++cnt[childKey];
                        cntCopy[childKey] = cnt[childKey];
                        visited[dataV[mappingSize]] = false;
                    }
                    else {
#ifdef COLLECT_STATISTICS
                        ++gNumIntermediate;
#endif
                        ++mappingSize;
                        pos[mappingSize] = 0;
                        if (!cn.interPos[mappingSize]) {
                            if (!cn.nodeOutPos[mappingSize].empty()) {
                                VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                startOffset[mappingSize] = outOffset[w];
                                candidate[mappingSize] = outNbors + outOffset[w];
                                candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                            }
                            else if (!cn.nodeInPos[mappingSize].empty()){
                                VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                startOffset[mappingSize] = inOffset[w];
                                candidate[mappingSize] = inNbors + inOffset[w];
                                candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                            }
                            else {
                                VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                startOffset[mappingSize] = unOffset[w];
                                candidate[mappingSize] = unNbors + unOffset[w];
                                candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                            }
                        }
                        else {
                            EdgeID e = 0;
                            if (cn.triEdgeType[mappingSize] != 0) {
                                int pos1 = cn.triPos[mappingSize].first, pos2 = cn.triPos[mappingSize].second;
                                if (cn.triEdgeType[mappingSize] == 1)
                                    e = startOffset[pos2] + pos[pos2] - 1;
                                else if (cn.triEdgeType[mappingSize] == 2)
                                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                                else if (cn.triEdgeType[mappingSize] == 5)
                                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                                else {
                                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                                    ++gNumEdgeID;
#endif
                                    e = dout.getEdgeID(key1, key2);
                                }
                            }
                            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !cn.candPos[mappingSize], e, cn.triEndType[mappingSize],
                                               cn.nodeInPos[mappingSize], cn.nodeOutPos[mappingSize], cn.nodeUnPos[mappingSize], candidate, candCount, tmp);
                        }
                        // apply the symmetry breaking rules
                        if (!cn.greaterPos[mappingSize].empty()) {
                            ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                            for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                    maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                            }
                            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                        }
                        if (!cn.lessPos[mappingSize].empty()) {
                            ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                            for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                    minTarget = dataV[cn.lessPos[mappingSize][i]];
                            }
                            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                        }
                    }
                }
                --mappingSize;
                if (mappingSize >= prefixSize)
                    visited[dataV[mappingSize]] = false;
            }
            // matching the remaining vertices to query the child table and compute orbit count
            ++mappingSize;
            pos[mappingSize] = 0;
            if (!cn.greaterPos[mappingSize].empty()) {
                ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                    if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                        maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                }
                pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
            }
            if (!cn.lessPos[mappingSize].empty()) {
                ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                    if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                        minTarget = dataV[cn.lessPos[mappingSize][i]];
                }
                candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
            }
            while (mappingSize >= prefixSize) {
                for (int j: cn.posAggreEdge[mappingSize]) {
                    if (cn.aggreEdgeType[j] == 1) aggreKey[j] = startOffset[mappingSize] + pos[mappingSize] - 1;
                    else if (cn.aggreEdgeType[j] == 2) aggreKey[j] = outID[startOffset[mappingSize] + pos[mappingSize] - 1];
                    else if (cn.aggreEdgeType[j] != 5){
                        VertexID key1 = dataV[cn.aggrePos[2 * j]], key2 = dataV[cn.aggrePos[2 * j + 1]];
                        if (cn.aggreEdgeType[j] == 3) aggreKey[j] = dout.getEdgeID(key1, key2);
                        else aggreKey[j] = dun.getUndirectedEID(key1, key2);
#ifdef COLLECT_STATISTICS
                        ++gNumEdgeID;
#endif
                    }
                    else aggreKey[j] = unID[startOffset[mappingSize] + pos[mappingSize] - 1];
                }
                while (pos[mappingSize] < candCount[mappingSize]) {
                    VertexID v = candidate[mappingSize][pos[mappingSize]];
                    ++pos[mappingSize];
                    if (visited[v]) continue;
                    visited[v] = true;
                    patternV[mappingSize] = cn.nodeOrder[mappingSize];
                    dataV[mappingSize] = v;
                    if (mappingSize == cn.nodeOrder.size() - 1) {
#ifdef COLLECT_STATISTICS
                        ++gNumMatch;
#endif
                        if (cn.childEdgeType == 0) childKey = cn.childKeyPos[0];
                        else {
                            if (cn.childEdgeType == 1) childKey = startOffset[cn.posChildEdge] + pos[cn.posChildEdge] - 1;
                            else if (cn.childEdgeType == 2) childKey = outID[startOffset[cn.posChildEdge] + pos[cn.posChildEdge] - 1];
                            else {
                                VertexID key1 = dataV[cn.childKeyPos[0]], key2 = dataV[cn.childKeyPos[1]];
                                if (cn.childEdgeType == 3) childKey = dout.getEdgeID(key1, key2);
                                else childKey = dun.getUndirectedEID(key1, key2);
#ifdef COLLECT_STATISTICS
                                ++gNumEdgeID;
#endif
                            }
                        }
                        if (orbitType == 0) {
                            cn.hashTables[0][0] += choosec(cnt[childKey], cn.num);
                        }
                        else if (orbitType == 1) {
                            for (int j = 0; j < cn.aggrePos.size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                VertexID ak = dataV[cn.aggrePos[j]];
                                if (cn.aggreType[j])
                                    cn.hashTables[j][ak] += choosec(cnt[childKey], cn.num);
                                else
                                    cn.hashTables[j][ak] += choosec(cntCopy[childKey] - 1, cn.num - 1);
                            }
                        }
                        else {
                            for (int j = 0; j < aggreKey.size(); ++j) {
#ifdef COLLECT_STATISTICS
                                ++gNumUpdate;
#endif
                                if (cn.aggreType[j])
                                    cn.hashTables[j][aggreKey[j]] += choosec(cnt[childKey], cn.num);
                                else
                                    cn.hashTables[j][aggreKey[j]] += choosec(cntCopy[childKey] - 1, cn.num - 1);
                            }
                        }
                        cnt[childKey] = 0;
                        visited[dataV[mappingSize]] = false;
                    }
                    else {
#ifdef COLLECT_STATISTICS
                        ++gNumIntermediate;
#endif
                        ++mappingSize;
                        pos[mappingSize] = 0;
                        if (!cn.interPos[mappingSize]) {
                            if (!cn.nodeOutPos[mappingSize].empty()) {
                                VertexID w = dataV[cn.nodeOutPos[mappingSize][0]];
                                startOffset[mappingSize] = outOffset[w];
                                candidate[mappingSize] = outNbors + outOffset[w];
                                candCount[mappingSize] = outOffset[w + 1] - outOffset[w];
                            }
                            else if (!cn.nodeInPos[mappingSize].empty()){
                                VertexID w = dataV[cn.nodeInPos[mappingSize][0]];
                                startOffset[mappingSize] = inOffset[w];
                                candidate[mappingSize] = inNbors + inOffset[w];
                                candCount[mappingSize] = inOffset[w + 1] - inOffset[w];
                            }
                            else {
                                VertexID w = dataV[cn.nodeUnPos[mappingSize][0]];
                                startOffset[mappingSize] = unOffset[w];
                                candidate[mappingSize] = unNbors + unOffset[w];
                                candCount[mappingSize] = unOffset[w + 1] - unOffset[w];
                            }
                        }
                        else {
                            EdgeID e = 0;
                            if (cn.triEdgeType[mappingSize] != 0) {
                                int pos1 = cn.triPos[mappingSize].first, pos2 = cn.triPos[mappingSize].second;
                                if (cn.triEdgeType[mappingSize] == 1)
                                    e = startOffset[pos2] + pos[pos2] - 1;
                                else if (cn.triEdgeType[mappingSize] == 2)
                                    e = outID[startOffset[pos2] + pos[pos2] - 1];
                                else if (cn.triEdgeType[mappingSize] == 5)
                                    e = unID[startOffset[pos2] + pos[pos2] - 1];
                                else {
                                    VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                                    ++gNumEdgeID;
#endif
                                    e = dout.getEdgeID(key1, key2);
                                }
                            }
                            generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !cn.candPos[mappingSize], e, cn.triEndType[mappingSize],
                                               cn.nodeInPos[mappingSize], cn.nodeOutPos[mappingSize], cn.nodeUnPos[mappingSize], candidate, candCount, tmp);

                        }
                        // apply the symmetry breaking rules
                        if (!cn.greaterPos[mappingSize].empty()) {
                            ui maxTarget = dataV[cn.greaterPos[mappingSize][0]];
                            for (int i = 1; i < cn.greaterPos[mappingSize].size(); ++i) {
                                if (dataV[cn.greaterPos[mappingSize][i]] > maxTarget)
                                    maxTarget = dataV[cn.greaterPos[mappingSize][i]];
                            }
                            pos[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], maxTarget);
                        }
                        if (!cn.lessPos[mappingSize].empty()) {
                            ui minTarget = dataV[cn.lessPos[mappingSize][0]];
                            for (int i = 1; i < cn.lessPos[mappingSize].size(); ++i) {
                                if (dataV[cn.lessPos[mappingSize][i]] < minTarget)
                                    minTarget = dataV[cn.lessPos[mappingSize][i]];
                            }
                            candCount[mappingSize] = firstPosGreaterThan(candidate[mappingSize], 0, candCount[mappingSize], minTarget);
                        }
                    }
                }
                --mappingSize;
                if (mappingSize >= 0)
                    visited[dataV[mappingSize]] = false;
            }
        }
    }
    delete[] cnt;
    delete[] cntCopy;
}

void executeShared(
        ExeParam &param,
        std::vector<ui *> &keyPoses,
        std::vector<ui> &keyPosSize,
        const std::vector<ui> &sizeBound,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        bool useTriangle,
        const Triangle &tri,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
        int mappingSize,
        bool *visited,
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    if (!useTriangle) {
        if (!param.rep.edgeKey) {
            executeSharedNode(param.rep, param.children, param.nIDs, param.trees, param.hashTables,
                              param.aggreWeights, param.id2AggreKey, param.id2ChildKey, param.isRoot,
                              param.childKeyPoses, param.aggrePoses, param.tableIDs, keyPoses, keyPosSize, sizeBound,
                              param.nodeInterPos, param.nodeInPos, param.nodeOutPos, param.nodeUnPos, param.greaterPos,
                              param.lessPos, param.nonLeafIDs, param.candidate,param.candCount, din, dout, dun,
                              param.p, outID, unID, reverseID, startOffset, patternV, dataV, mappingSize, visited, pos, tmp, allV);
        }
        else {
            executeSharedNodeEdgeKey(param.rep, param.children, param.nIDs, param.trees, param.hashTables,
                                     param.aggreWeights, param.id2AggreKey, param.id2ChildKey, param.isRoot,
                                     param.childKeyPoses, param.aggrePoses, param.tableIDs, keyPoses, keyPosSize, sizeBound,
                                     param.nodeInterPos, param.nodeInPos, param.nodeOutPos, param.nodeUnPos, param.greaterPos,
                                     param.lessPos, param.posChildEdge, param.posAggreEdge, param.childEdgeType, param.aggreEdgeType,
                                     param.nonLeafIDs, param.candidate, param.candCount, din, dout, dun, param.p,
                                     outID, unID, reverseID, startOffset, patternV, dataV, mappingSize, visited, pos, tmp, allV);

        }
    }
    else {
        if (!param.rep.edgeKey) {
            executeSharedNodeT(param.rep, param.children, param.nIDs, param.trees, param.hashTables,
                               param.aggreWeights, param.id2AggreKey, param.id2ChildKey, param.isRoot,
                               param.childKeyPoses, param.aggrePoses, param.tableIDs, keyPoses, keyPosSize, sizeBound,
                               param.nodeInterPos, param.nodeInPos, param.nodeOutPos, param.nodeUnPos, param.nodeCandPos,
                               param.greaterPos, param.lessPos, param.nodeTriPos, param.triEdgeType, param.triEndType,
                               param.nonLeafIDs, param.candidate, param.candCount, din, dout, dun, tri, param.p,
                               outID, unID, reverseID, startOffset, patternV, dataV, mappingSize, visited, pos, tmp, allV);
        }
        else {
            executeSharedNodeEdgeKeyT(param.rep, param.children, param.nIDs, param.trees, param.hashTables,
                                      param.aggreWeights, param.id2AggreKey, param.id2ChildKey, param.isRoot,
                                      param.childKeyPoses, param.aggrePoses, param.tableIDs, keyPoses, keyPosSize, sizeBound,
                                      param.nodeInterPos, param.nodeInPos, param.nodeOutPos, param.nodeUnPos, param.nodeCandPos,
                                      param.greaterPos, param.lessPos, param.posChildEdge, param.posAggreEdge, param.childEdgeType,
                                      param.aggreEdgeType, param.nodeTriPos, param.nonLeafIDs, param.triEdgeType,
                                      param.triEndType, param.candidate, param.candCount, din, dout, dun, tri, param.p,
                                      outID, unID, reverseID, startOffset, patternV, dataV, mappingSize, visited, pos, tmp, allV);
        }
    }
}

void executeForest(Forest &f, std::vector<HashTable> &result, const DataGraph &din,
                   const DataGraph &dout, const DataGraph &dun, bool useTriangle, const Triangle &tri,
                   std::vector<HashTable> &allocatedHashTable, EdgeID *outID, EdgeID *unID, EdgeID *reverseID,
                   EdgeID *startOffset, VertexID *patternV, VertexID *dataV, bool *visited, ui *candPos,
                   VertexID *&tmp, VertexID *allV, specialsparse *sg, VertexID *cliqueVertices) {
    if (f.numQuery == 0) return;
    EdgeID *inOffset = din.getOffsets();
    VertexID *inNbors = din.getNbors();
    EdgeID *outOffset = dout.getOffsets();
    VertexID *outNbors = dout.getNbors();
    EdgeID *unOffset = dun.getOffsets();
    VertexID *unNbors = dun.getNbors();
    int orbitType = f.orbitType;
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    // allocate memory. 1. for each query, one hash table for the overall result and factor 1, one for each other
    // factor. 2. for each non-root. root node directly update hash table in queryFactorH
    int nTable = 0, mTable = 0;
    HashTable **queryFactorH = new HashTable *[result.size()];
    for (int i = 0; i < result.size(); ++i) {
        queryFactorH[i] = new HashTable[f.queryFactor[i].size()];
    }
    HashTable **treeNodeH = new HashTable *[f.allTree.size()];
    for (int tID = 0; tID < f.allTree.size(); ++tID)
        treeNodeH[tID] = new HashTable[f.allTree[tID].getNumNodes()];

    for (int qID = 0; qID < result.size(); ++qID) {
        for (int fID = 0; fID < f.queryFactor[qID].size(); ++fID) {
            int divideFactor = f.queryFactor[qID][fID];
            int firstTree = f.queryFactorFirstTree[qID][fID];
            int numTree = f.queryFactorNumTree[qID][fID];
            if (fID == 0) queryFactorH[qID][fID] = result[qID];
            else {
                if (orbitType == 0) {
                    queryFactorH[qID][fID] = new Count[1];
                    queryFactorH[qID][fID][0] = 0;
                }
                else if (orbitType == 1) {
                    queryFactorH[qID][fID] = new Count[n];
                    memset(queryFactorH[qID][fID], 0, sizeof(Count) * n);
                    ++nTable;
                }
                else {
                    if (f.allPattern[firstTree].u.isEOrbitDir()) {
                        queryFactorH[qID][fID] = new Count[m / 2];
                        memset(queryFactorH[qID][fID], 0, sizeof(Count) * m / 2);
                        ++mTable;
                    }
                    else {
                        queryFactorH[qID][fID] = new Count[m];
                        memset(queryFactorH[qID][fID], 0, sizeof(Count) * m);
                        mTable += 2;
                    }
                }
                allocatedHashTable.push_back(queryFactorH[qID][fID]);
            }

            for (int tID = firstTree; tID < numTree + firstTree; ++tID) {
                const Tree &t = f.allTree[tID];
                VertexID nID = t.getRootID();
                treeNodeH[tID][nID] = queryFactorH[qID][fID];
#ifdef DEBUG
                treeNodeH[tID][nID] = new Count[n];
                memset(treeNodeH[tID][nID], 0, sizeof(Count) * n);
#endif
            }
        }
    }
    f.allocateTables(treeNodeH, allocatedHashTable, n, m, nTable, mTable);
    double tableSize = (double)(nTable * n + mTable * m) * sizeof(Count) / 1e9;
//    std::cout << "allocated memory for hash tables in executeForest: size = " <<  tableSize << "GB" << std::endl;
//    std::cout << "number of vertex tables: " << nTable << ", number of edge tables: " << mTable << std::endl;
    nTable = 0;
    mTable = 0;
    // allocate keyPos tables
    std::vector<ui *> keyPoses(f.childTables.size(), nullptr);
    std::vector<ui> keyPosSize(f.childTables.size(), 0);
    std::vector<ui> sizeBound(f.childTables.size(), 0);
    for (int i = 0; i < f.childTables.size(); ++i) {
        const Node &rep = f.allTree[f.childTID[i][0]].getNode(f.childNID[i][0]);
        if (rep.prefixSize != 0 && rep.keySize == 0) {
            keyPoses[i] = new ui[1];
            memset(keyPoses[i], 0, sizeof(ui));
            sizeBound[i] = 1;
        }
        else if (rep.prefixSize != 0 && rep.keySize == 1) {
            keyPoses[i] = new ui[n / 8 + 1];
            memset(keyPoses[i], 0, sizeof(ui) * (n / 8 + 1));
            sizeBound[i] = n / 8 + 1;
            ++nTable;
        }
        else if (rep.prefixSize != 0 && rep.keySize == 2){
            keyPoses[i] = new ui[m / 8 + 1];
            memset(keyPoses[i], 0, sizeof(ui) * (m / 8 + 1));
            sizeBound[i] = m / 8 + 1;
            ++mTable;
        }
    }
    tableSize = (double)(nTable * n + mTable * m) * sizeof(ui) / 8000000000;
//    std::cout << "allocated memory for key poses: size = " <<  tableSize << "GB" << std::endl;
    std::vector<std::vector<VertexID **>> treeNodeCand(f.allTree.size());
    std::vector<std::vector<ui *>> treeNodeCandCount(f.allTree.size());
    for (int i = 0; i < f.allTree.size(); ++i) {
        treeNodeCand[i] = std::vector<VertexID **>(f.allTree[i].getNumNodes(), nullptr);
        treeNodeCandCount[i] = std::vector<ui *>(f.allTree[i].getNumNodes(), nullptr);
    }
    // group the partitions that can be executed together. Then execute these
    // partitions together using exeSharedNode. group ID: 0 for partition size = 0, 1 for size = 1, 3 for size = 2,
    // out edge, 4 for size = 2, out edge. for types of partition order, generate group ID according to canonValue
    int maxGroupID = 5;
    std::map<CanonType, int> partCanon2GroupID;
    std::vector<std::vector<VertexID>> groupPID(5);
    std::vector<std::vector<int>> groupTID(5);
    std::vector<std::vector<int>> groupDepthNodeNum(5);
    groupDepthNodeNum[1] = std::vector<int>(1, 0);
    groupDepthNodeNum[2] = std::vector<int>(2, 0);
    groupDepthNodeNum[3] = std::vector<int>(2, 0);
    groupDepthNodeNum[4] = std::vector<int>(2, 0);
    for (int qID = 0; qID < result.size(); ++qID) {
        for (int fID = 0; fID < f.queryFactor[qID].size(); ++fID) {
            int firstTree = f.queryFactorFirstTree[qID][fID];
            int numTree = f.queryFactorNumTree[qID][fID];
            for (int tID = firstTree; tID < numTree + firstTree; ++tID) {
                const Tree &t = f.allTree[tID];
                const std::vector<std::vector<VertexID>> &globalOrder = t.getGlobalOrder();
                for (VertexID pID = 0; pID < t.getGlobalOrder().size(); ++pID) {
                    const std::vector<VertexID> &partitionOrder = globalOrder[pID];
                    const std::vector<std::vector<int>> &partitionInPos = t.getPartitionInPos(pID);
                    const std::vector<std::vector<int>> &partitionOutPos = t.getPartitionOutPos(pID);
                    const std::vector<std::vector<int>> &partitionUnPos = t.getPartitionUnPos(pID);
                    const std::vector<std::vector<std::vector<VertexID>>> &nodesAtStep = t.getNodesAtStep();
                    const Pattern &p = f.allPattern[tID];
                    if (partitionOrder.empty()) {
                        groupPID[0].push_back(pID);
                        groupTID[0].push_back(tID);
                    }
                    else if (partitionOrder.size() == 1) {
                        groupPID[1].push_back(pID);
                        groupTID[1].push_back(tID);
                        for (VertexID nID: nodesAtStep[pID][0]) {
                            if (f.nodeNeedExecute[tID][nID])
                                ++groupDepthNodeNum[1][0];
                        }
                    }
                    else if (partitionOrder.size() == 2) {
                        int groupID;
                        if (!partitionUnPos[1].empty()) groupID = 2;
                        else if (!partitionOutPos[1].empty()) groupID = 3;
                        else groupID = 4;
                        groupPID[groupID].push_back(pID);
                        groupTID[groupID].push_back(tID);
                        for (int i = 0; i < 2; ++i) {
                            for (VertexID nID: nodesAtStep[pID][i]) {
                                if (f.nodeNeedExecute[tID][nID])
                                    ++groupDepthNodeNum[groupID][i];
                            }
                        }
                    }
                    else {
                        CanonType canonValue;
                        if (p.useDAG())
                            canonValue = subgraphCanonValue(p.out, partitionOrder);
                        else
                            canonValue = subgraphCanonValue(p.u, partitionOrder);
                        if (partCanon2GroupID.find(canonValue) != partCanon2GroupID.end()) {
                            int gID = partCanon2GroupID[canonValue];
                            groupPID[gID].push_back(pID);
                            groupTID[gID].push_back(tID);
                            for (int i = 0; i < partitionOrder.size(); ++i) {
                                for (VertexID nID: nodesAtStep[pID][i]) {
                                    if (f.nodeNeedExecute[tID][nID])
                                        ++groupDepthNodeNum[gID][i];
                                }
                            }
                        }
                        else {
                            partCanon2GroupID[canonValue] = maxGroupID;
                            groupPID.emplace_back();
                            groupTID.emplace_back();
                            groupDepthNodeNum.emplace_back();
                            groupPID[maxGroupID].push_back(pID);
                            groupTID[maxGroupID].push_back(tID);
                            groupDepthNodeNum[maxGroupID] = std::vector<int>(partitionOrder.size(), 0);
                            for (int i = 0; i < partitionOrder.size(); ++i) {
                                for (VertexID nID: nodesAtStep[pID][i]) {
                                    if (f.nodeNeedExecute[tID][nID])
                                        ++groupDepthNodeNum[maxGroupID][i];
                                }
                            }
                            ++maxGroupID;
                        }
                    }
                }
            }
        }
    }
    // execute cliques first
    if (sg != nullptr) {
        for (int qID = 0; qID < result.size(); ++qID) {
            int tID = f.queryFactorFirstTree[qID][0];
            if (f.allPattern[tID].u.isClique() && sg != nullptr) {
                int k = f.allPattern[tID].u.getNumVertices();
                if (k < 4) continue;
                int aggreWeight = f.allTree[tID].getAggreWeight()[0];
                mkspecial(sg, k);
                kclique(k, k, sg, cliqueVertices, result[qID], orbitType);
                freesub(sg, k);
                f.nodeExecuted[tID][0] = true;
                if (aggreWeight != 1) {
                    if (orbitType == 0) result[qID][0] *= aggreWeight;
                    else if (orbitType == 1) {
                        for (VertexID v = 0; v < n; ++v)
                            result[qID][v] *= aggreWeight;
                    }
                    else {
                        for (EdgeID e = 0; e < m; ++e)
                            result[qID][e] *= aggreWeight;
                    }
                }
            }
        }
    }

    // execute some partitions with no prefix.
    std::vector<ExeParam> noPrefixParams;
    bool noNewNode = false;
    while (!noNewNode) {
        noNewNode = true;
        for (int i = 0; i < groupTID[0].size(); ++i) {
            int tID = groupTID[0][i];
            VertexID pID = groupPID[0][i];
            const Tree &t = f.allTree[tID];
            const std::vector<VertexID> &postOrder = t.getPostOrder();
            const std::vector<int> &partitionPos = t.getPartitionPos();
            VertexID nID;
            if (pID == 0) nID = postOrder[0];
            else nID = postOrder[partitionPos[pID - 1] + 1];
            if (!f.nodeExecutable[tID][nID]) continue;
            ExeParam param;
            int pos = f.posInMergedNode[tID][nID];
            param.p = f.mergedNodePattern(pos);
            for (const NodeInfo& inf: f.toTreeNode[pos]) {
                if (f.executable(inf)) {
                    f.nodeExecuted[inf.tID][inf.nID] = true;
                    const Tree &t2 = f.allTree[inf.tID];
                    // add this node for execution
                    param.nIDs.push_back(inf.nID);
                    param.tIDs.push_back(inf.tID);
                    int tableID = f.posInChildTables[inf.tID][inf.nID];
                    param.tableIDs.push_back(tableID);
                    param.trees.push_back(t2);
                    bool isRt = t2.getRootID() == inf.nID;
                    param.aggreWeights.push_back(inf.aggreWeight);
                    const std::vector<std::vector<VertexID>> &child = t2.getChild();
                    param.children.push_back(child[inf.nID]);
                    param.assignGroup(inf.aggrePos, inf.childKeyPos, isRt, 0, n, m);
                    param.hashTables.push_back(treeNodeH[inf.tID]);
                }
            }
            if (param.nIDs.empty()) continue;
            else noNewNode = false;
            // allocate candidate and cand count before execution
            VertexID **candidate;
            ui *candCount;
            int exeTID = f.mergedNodeTID[pos];
            VertexID exeNID = f.mergedNodeNID[pos];
            const Tree &exeTree = f.allTree[exeTID];
            param.rep = exeTree.getNode(exeNID);
            if (treeNodeCand[exeTID][exeNID] == nullptr) {
                candidate = new VertexID *[param.rep.localOrder.size()];
                const std::vector<bool> &nodeCandPos = f.allTree[exeTID].getNodeCandPos(exeNID);
                for (int j = 0; j < nodeCandPos.size(); ++j) {
                    if (nodeCandPos[j])
                        candidate[j] = new VertexID[dout.getNumVertices()];
                }
                candCount = new ui[param.rep.localOrder.size()];
                treeNodeCand[exeTID][exeNID] = candidate;
                treeNodeCandCount[exeTID][exeNID] = candCount;
            }
            else {
                candidate = treeNodeCand[exeTID][exeNID];
                candCount = treeNodeCandCount[exeTID][exeNID];
            }
            param.candidate = candidate;
            param.candCount = candCount;
            param.nodeInterPos = exeTree.getNodeInterPos(exeNID);
            param.nodeInPos = exeTree.getNodeInPos(exeNID);
            param.nodeOutPos = exeTree.getNodeOutPos(exeNID);
            param.nodeUnPos = exeTree.getNodeUnPos(exeNID);
            param.nodeCandPos = exeTree.getNodeCandPos(exeNID);
            param.greaterPos = exeTree.getNodeGreaterPos(exeNID);
            param.lessPos = exeTree.getNodeLessPos(exeNID);
            param.nodeTriPos = exeTree.getNodeTriPos(exeNID);
            param.triEdgeType = exeTree.getTriEdgeType(exeNID);
            param.triEndType = exeTree.getTriEndType(exeNID);
            param.buildPosesAndTypes(f);
            for (int tableID : param.tableIDs) {
                // update executable nodes
                if (tableID != -1) {
                    ++f.numExecuted[tableID];
                    if (f.numExecuted[tableID] == f.numTables[tableID])
                        f.updateNodeExecutable(tableID);
                }
            }
            noPrefixParams.push_back(param);
        }
    }
    for (ExeParam &param : noPrefixParams) {
        executeShared(param, keyPoses, keyPosSize, sizeBound, din, dout, dun, useTriangle, tri, outID, unID, reverseID,
                      startOffset, patternV, dataV, 0, visited, candPos, tmp, allV);
    }
    bool noNewNodeOuter = false;
    while (!noNewNodeOuter) {
        noNewNodeOuter = true;
        for (int groupID = 1; groupID < maxGroupID; ++groupID) {
            nTable = 0;
            mTable = 0;
            if (groupTID[groupID].empty()) continue;
            const std::vector<int> &gTID = groupTID[groupID];
            const std::vector<VertexID> &gPID = groupPID[groupID];
            const Tree &firstTree = f.allTree[gTID[0]];
            const std::vector<VertexID> &partitionOrder = firstTree.getPartitionOrder(gPID[0]);
            const std::vector<std::vector<int>> &partitionInPos = firstTree.getPartitionInPos(gPID[0]);
            const std::vector<std::vector<int>> &partitionOutPos = firstTree.getPartitionOutPos(gPID[0]);
            const std::vector<std::vector<int>> &partitionUnPos = firstTree.getPartitionUnPos(gPID[0]);
            const std::vector<bool> &partitionInterPos = firstTree.getPartitionInterPos(gPID[0]);
            const std::vector<std::vector<int>> &greaterPos = firstTree.getPartitionGreaterPos(gPID[0]);
            const std::vector<std::vector<int>> &lessPos = firstTree.getPartitionLessPos(gPID[0]);
            const std::vector<bool> &partitionCandPos = firstTree.getPartitionCandPos(gPID[0]);
            const std::vector<std::pair<int, int>> &partitionTriPos = firstTree.getPartitionTriPos(gPID[0]);
            const std::vector<int> &triEdgeType = firstTree.getPartitionEdgeType(gPID[0]);
            const std::vector<int> &triEndType = firstTree.getPartitionEndType(gPID[0]);
            // first traverse the partitions to decide at each step which nodes are executed
            int length = (int)partitionOrder.size();
            std::vector<std::vector<ExeParam>> depthParams(length);
            for (int depth = 0; depth < length; ++depth) {
                bool noNewNodeInner = false;
                std::vector<ExeParam> params;
                while (!noNewNodeInner) {
                    noNewNodeInner = true;
                    for (int i = 0; i < gTID.size(); ++i) {
                        int tID = gTID[i];
                        VertexID pID = gPID[i];
                        const Tree &t = f.allTree[tID];
                        const std::vector<std::vector<std::vector<VertexID>>> &nodesAtStep = t.getNodesAtStep();
                        // make sure that child partitions are executed
                        VertexID nID = MAX_PATTERN_SIZE;
                        for (VertexID j : nodesAtStep[pID][depth]) {
                            if (f.nodeExecuted[tID][j]) continue;
                            else {
                                nID = j;
                                break;
                            }
                        }
                        if (nID == MAX_PATTERN_SIZE) continue;
                        // execute the first unexecuted node
                        int pos = f.posInMergedNode[tID][nID];
                        ExeParam param;
                        // collect other nodes to be executed together
                        for (const NodeInfo& inf: f.toTreeNode[pos]) {
                            if (f.executable(inf)) {
                                f.nodeExecuted[inf.tID][inf.nID] = true;
                                const Tree &t2 = f.allTree[inf.tID];
                                VertexID partID = f.nodePID[inf.tID][inf.nID];
                                if (!t2.getPartitionOrder(partID).empty())
                                    if (param.nIDs.empty()) {
                                        const Pattern &pt = f.allPattern[inf.tID];
                                        param.p = pt;
                                        param.rep = t2.getNode(inf.nID);
                                        param.nodeInterPos = t2.getNodeInterPos(inf.nID);
                                        param.nodeInPos = t2.getNodeInPos(inf.nID);
                                        param.nodeOutPos = t2.getNodeOutPos(inf.nID);
                                        param.nodeUnPos = t2.getNodeUnPos(inf.nID);
                                        param.nodeCandPos = t2.getNodeCandPos(inf.nID);
                                        param.greaterPos = t2.getNodeGreaterPos(inf.nID);
                                        param.lessPos = t2.getNodeLessPos(inf.nID);
                                        param.nodeTriPos = t2.getNodeTriPos(inf.nID);
                                        param.triEdgeType = t2.getTriEdgeType(inf.nID);
                                        param.triEndType = t2.getTriEndType(inf.nID);
                                    }
                                // add this node for execution
                                param.nIDs.push_back(inf.nID);
                                param.tIDs.push_back(inf.tID);
                                int tableID = f.posInChildTables[inf.tID][inf.nID];
                                param.tableIDs.push_back(tableID);
                                param.trees.push_back(t2);
                                bool isRt = t2.getRootID() == inf.nID;
                                param.aggreWeights.push_back(inf.aggreWeight);
                                const std::vector<std::vector<VertexID>> &child = t2.getChild();
                                param.children.push_back(child[inf.nID]);
                                ui prefixSize = t2.getNode(inf.nID).prefixSize;
                                ui keySize = t2.getNode(inf.nID).keySize;
                                if (prefixSize == 0) keySize = 0;
                                param.assignGroup(inf.aggrePos, inf.childKeyPos, isRt, keySize, n, m);
                                param.hashTables.push_back(treeNodeH[inf.tID]);
                            }
                        }
                        if (param.nIDs.empty()) continue;
                        else noNewNodeInner = false;
                        // allocate candidate and cand count before execution
                        VertexID **candidate;
                        ui *candCount;
                        int exeTID = f.mergedNodeTID[pos];
                        VertexID exeNID = f.mergedNodeNID[pos];
                        const Node &tau = f.allTree[exeTID].getNode(exeNID);
                        if (treeNodeCand[exeTID][exeNID] == nullptr) {
                            candidate = new VertexID *[tau.nodeOrder.size()];
                            const std::vector<bool> &nodeCandPos = f.allTree[exeTID].getNodeCandPos(exeNID);
                            for (int j = 0; j < nodeCandPos.size(); ++j) {
                                if (nodeCandPos[j])
                                    candidate[j] = new VertexID[dout.getNumVertices()];
                            }
                            candCount = new ui[tau.nodeOrder.size()];
                            treeNodeCand[exeTID][exeNID] = candidate;
                            treeNodeCandCount[exeTID][exeNID] = candCount;
                        }
                        else {
                            candidate = treeNodeCand[exeTID][exeNID];
                            candCount = treeNodeCandCount[exeTID][exeNID];
                        }
                        param.candidate = candidate;
                        param.candCount = candCount;
                        param.buildPosesAndTypes(f);
                        for (int tableID : param.tableIDs) {
                            // update executable nodes
                            if (tableID != -1 && !f.isPartition[tableID]) {
                                ++f.numExecuted[tableID];
                                if (f.numExecuted[tableID] == f.numTables[tableID])
                                    f.updateNodeExecutable(tableID);
                            }
                        }
                        params.push_back(param);
                    }
                }
                depthParams[depth] = params;
            }
            VertexID **partitionCandidate = new VertexID *[partitionOrder.size()];
            for (int i = 0; i < partitionOrder.size(); ++i) {
                if (partitionInterPos[i])
                    partitionCandidate[i] = new VertexID[dout.getNumVertices()];
            }
            ui *partitionCandCount = new ui[partitionOrder.size()];
            partitionCandidate[0] = allV;
            partitionCandCount[0] = n;
            candPos[0] = 0;
            int depth = 0;
            int mappingSize = 0;
            while (depth >= 0) {
                while (candPos[depth] < partitionCandCount[depth]) {
                    VertexID v = partitionCandidate[depth][candPos[depth]];
                    ++candPos[depth];
                    if (visited[v]) continue;
                    visited[v] = true;
                    patternV[depth] = partitionOrder[depth];
                    dataV[depth] = v;
                    std::vector<ExeParam> &params = depthParams[depth];
                    // for each shared node, clear count
                    for (auto & param : params) {
                        for (int tableID: param.tableIDs) {
                            if (tableID == -1) continue;
                            if (sizeBound[tableID] == n / 8 + 1) {
                                if (keyPosSize[tableID] < n / 8 + 1) {
                                    for (int j = 0; j < keyPosSize[tableID]; ++j)
                                        f.childTables[tableID][keyPoses[tableID][j]] = 0;
                                } else {
                                    memset(f.childTables[tableID], 0, sizeof(Count) * n);
                                }
                            }
                            else if (sizeBound[tableID] == 1) {
                                f.childTables[tableID][keyPoses[tableID][0]] = 0;
                            }
                            else {
                                if (keyPosSize[tableID] < m / 8 + 1) {
                                    for (int j = 0; j < keyPosSize[tableID]; ++j)
                                        f.childTables[tableID][keyPoses[tableID][j]] = 0;
                                } else {
                                    memset(f.childTables[tableID], 0, sizeof(Count) * m);
                                }
                            }
                            keyPosSize[tableID] = 0;
                        }
                    }

                    for (ExeParam &param : params) {
                        executeShared(param, keyPoses, keyPosSize, sizeBound, din, dout, dun, useTriangle, tri, outID, unID, reverseID,
                                      startOffset, patternV, dataV, mappingSize + 1, visited, candPos, tmp, allV);
                    }
                    if (depth == length - 1) {
                        visited[dataV[mappingSize]] = false;
                    }
                    else {
                        ++mappingSize;
                        ++depth;
                        candPos[depth] = 0;
                        if (!partitionInterPos[depth]) {
                            if (!partitionOutPos[depth].empty()) {
                                VertexID w = dataV[partitionOutPos[depth][0]];
                                startOffset[mappingSize] = outOffset[w];
                                partitionCandidate[depth] = outNbors + outOffset[w];
                                partitionCandCount[depth] = outOffset[w + 1] - outOffset[w];
                            }
                            else if (!partitionInPos[depth].empty()){
                                VertexID w = dataV[partitionInPos[depth][0]];
                                startOffset[mappingSize] = inOffset[w];
                                partitionCandidate[depth] = inNbors + inOffset[w];
                                partitionCandCount[depth] = inOffset[w + 1] - inOffset[w];
                            }
                            else {
                                VertexID w = dataV[partitionUnPos[depth][0]];
                                startOffset[mappingSize] = unOffset[w];
                                partitionCandidate[depth] = unNbors + unOffset[w];
                                partitionCandCount[depth] = unOffset[w + 1] - unOffset[w];
                            }
                        }
                        else {
                            if (useTriangle) {
                                EdgeID e = 0;
                                if (triEdgeType[mappingSize] != 0) {
                                    int pos1 = partitionTriPos[mappingSize].first, pos2 = partitionTriPos[mappingSize].second;
                                    if (triEdgeType[mappingSize] == 1)
                                        e = startOffset[pos2] + candPos[pos2] - 1;
                                    else if (triEdgeType[mappingSize] == 2)
                                        e = outID[startOffset[pos2] + candPos[pos2] - 1];
                                    else if (triEdgeType[mappingSize] == 5)
                                        e = unID[startOffset[pos2] + candPos[pos2] - 1];
                                    else{
                                        VertexID key1 = dataV[pos1], key2 = dataV[pos2];
#ifdef COLLECT_STATISTICS
                                        ++gNumEdgeID;
#endif
                                        e = dout.getEdgeID(key1, key2);
                                    }
                                }
                                generateCandidateT(din, dout, dun, tri, dataV, mappingSize, !partitionCandPos[mappingSize], e, triEndType[mappingSize],
                                                   partitionInPos[mappingSize], partitionOutPos[mappingSize], partitionUnPos[mappingSize], partitionCandidate, partitionCandCount, tmp);

                            }
                            else {
                                generateCandidate(din, dout, dun, dataV, depth, partitionInPos[depth], partitionOutPos[depth],
                                                  partitionUnPos[depth], partitionCandidate, partitionCandCount, tmp);
                            }
                        }
                        // apply the symmetry breaking rules
                        if (!greaterPos[mappingSize].empty()) {
                            ui maxTarget = dataV[greaterPos[mappingSize][0]];
                            for (int i = 1; i < greaterPos[mappingSize].size(); ++i) {
                                if (dataV[greaterPos[mappingSize][i]] > maxTarget)
                                    maxTarget = dataV[greaterPos[mappingSize][i]];
                            }
                            candPos[mappingSize] = firstPosGreaterThan(partitionCandidate[mappingSize], 0, partitionCandCount[mappingSize], maxTarget);
                        }
                        if (!lessPos[mappingSize].empty()) {
                            ui minTarget = dataV[lessPos[mappingSize][0]];
                            for (int i = 1; i < lessPos[mappingSize].size(); ++i) {
                                if (dataV[lessPos[mappingSize][i]] < minTarget)
                                    minTarget = dataV[lessPos[mappingSize][i]];
                            }
                            partitionCandCount[mappingSize] = firstPosGreaterThan(partitionCandidate[mappingSize], 0, partitionCandCount[mappingSize], minTarget);
                        }
                    }
                }
                --depth;
                --mappingSize;
                if (depth >= 0)
                    visited[dataV[depth]] = false;
            }
            for (int j = 0; j < partitionOrder.size(); ++j) {
                if (partitionCandPos[j])
                    delete[] partitionCandidate[j];
            }
            delete[] partitionCandidate;
            delete[] partitionCandCount;
            // set the parent nodes of these partitions as executable
            const std::vector<ExeParam> &lastParams = depthParams[length - 1];
            if (!lastParams.empty()) noNewNodeOuter = false;
            for (int i = 0; i < lastParams.size(); ++i) {
                for (int tableID : lastParams[i].tableIDs) {
                    if (tableID != -1 && f.isPartition[tableID]) {
                        ++f.numExecuted[tableID];
                        if (f.numExecuted[tableID] == f.numTables[tableID])
                            f.updateNodeExecutable(tableID);
                    }
                }
            }
        }
    }

//    std::cout << "time for nodes that need prefix: " << exeTime - firstTime << "s" << std::endl;
    // execute the remaining partitions
    noNewNode = false;
    noPrefixParams.clear();
    while (!noNewNode) {
        noNewNode = true;
        for (int i = 0; i < groupTID[0].size(); ++i) {
            int tID = groupTID[0][i];
            VertexID pID = groupPID[0][i];
            const Tree &t = f.allTree[tID];
            const std::vector<VertexID> &postOrder = t.getPostOrder();
            const std::vector<int> &partitionPos = t.getPartitionPos();
            VertexID nID;
            if (pID == 0) nID = postOrder[0];
            else nID = postOrder[partitionPos[pID - 1] + 1];
            if (!f.nodeExecutable[tID][nID]) continue;
            ExeParam param;
            int pos = f.posInMergedNode[tID][nID];
            param.p = f.mergedNodePattern(pos);
            for (const NodeInfo& inf: f.toTreeNode[pos]) {
                if (f.executable(inf)) {
                    f.nodeExecuted[inf.tID][inf.nID] = true;
                    const Tree &t2 = f.allTree[inf.tID];
                    // add this node for execution
                    param.nIDs.push_back(inf.nID);
                    param.tIDs.push_back(inf.tID);
                    int tableID = f.posInChildTables[inf.tID][inf.nID];
                    param.tableIDs.push_back(tableID);
                    param.trees.push_back(t2);
                    bool isRt = t2.getRootID() == inf.nID;
                    param.aggreWeights.push_back(inf.aggreWeight);
                    const std::vector<std::vector<VertexID>> &child = t2.getChild();
                    param.children.push_back(child[inf.nID]);
                    param.assignGroup(inf.aggrePos, inf.childKeyPos, isRt, 0, n, m);
                    param.hashTables.push_back(treeNodeH[inf.tID]);
                }
            }
            if (param.nIDs.empty()) continue;
            else noNewNode = false;
            // allocate candidate and cand count before execution
            VertexID **candidate;
            ui *candCount;
            int exeTID = f.mergedNodeTID[pos];
            VertexID exeNID = f.mergedNodeNID[pos];
            const Tree &exeTree = f.allTree[exeTID];
            param.rep = exeTree.getNode(exeNID);
            if (treeNodeCand[exeTID][exeNID] == nullptr) {
                candidate = new VertexID *[param.rep.localOrder.size()];
                const std::vector<bool> &nodeCandPos = f.allTree[exeTID].getNodeCandPos(exeNID);
                for (int j = 0; j < nodeCandPos.size(); ++j) {
                    if (nodeCandPos[j])
                        candidate[j] = new VertexID[dout.getNumVertices()];
                }
                candCount = new ui[param.rep.localOrder.size()];
                treeNodeCand[exeTID][exeNID] = candidate;
                treeNodeCandCount[exeTID][exeNID] = candCount;
            }
            else {
                candidate = treeNodeCand[exeTID][exeNID];
                candCount = treeNodeCandCount[exeTID][exeNID];
            }
            param.candidate = candidate;
            param.candCount = candCount;
            param.nodeInterPos = exeTree.getNodeInterPos(exeNID);
            param.nodeInPos = exeTree.getNodeInPos(exeNID);
            param.nodeOutPos = exeTree.getNodeOutPos(exeNID);
            param.nodeUnPos = exeTree.getNodeUnPos(exeNID);
            param.nodeCandPos = exeTree.getNodeCandPos(exeNID);
            param.greaterPos = exeTree.getNodeGreaterPos(exeNID);
            param.lessPos = exeTree.getNodeLessPos(exeNID);
            param.nodeTriPos = exeTree.getNodeTriPos(exeNID);
            param.triEdgeType = exeTree.getTriEdgeType(exeNID);
            param.triEndType = exeTree.getTriEndType(exeNID);
            param.buildPosesAndTypes(f);
            for (int tableID : param.tableIDs) {
                // update executable nodes
                if (tableID != -1) {
                    ++f.numExecuted[tableID];
                    if (f.numExecuted[tableID] == f.numTables[tableID])
                        f.updateNodeExecutable(tableID);
                }
            }
            noPrefixParams.push_back(param);
        }
    }
    for (ExeParam &param : noPrefixParams) {
        executeShared(param, keyPoses, keyPosSize, sizeBound, din, dout, dun, useTriangle, tri, outID, unID, reverseID,
                      startOffset, patternV, dataV, 0, visited, candPos, tmp, allV);
    }
//    std::cout << "time for other remaining nodes: " << elapsedSeconds.count() << std::endl;
    // release memory and prepare count for output
    for (int tID = 0; tID < treeNodeCand.size(); ++tID) {
        for (VertexID nID = 0; nID < treeNodeCand[tID].size(); ++nID) {
            const std::vector<bool> &nodeCandPos = f.allTree[tID].getNodeCandPos(nID);
            for (int i = 0; i < nodeCandPos.size(); ++i) {
                if (nodeCandPos[i] && treeNodeCand[tID][nID] != nullptr)
                    delete[] treeNodeCand[tID][nID][i];
            }
            delete[] treeNodeCand[tID][nID];
            delete[] treeNodeCandCount[tID][nID];
        }
    }
#ifdef DEBUG
    HashTable hPattern = new Count[n];
    memset(hPattern, 0, sizeof(Count) * n);
#endif
    for (int qID = 0; qID < result.size(); ++qID) {
        for (int fID = 0; fID < f.queryFactor[qID].size(); ++fID) {
            int divideFactor = f.queryFactor[qID][fID];
            int firstTree = f.queryFactorFirstTree[qID][fID];
            int multiFactor = f.allTree[firstTree].getMultiFactor();
#ifdef DEBUG
            if (divideFactor != 2) continue;
            int firstID = f.queryFactorFirstTree[qID][fID];
            int lastID = firstID + f.queryFactorNumTree[qID][fID];
            for (int tID = firstID; tID < lastID; ++tID) {
                if (f.allPattern[tID].u.getCanonValue() == 282161) {
                    VertexID rootID = f.allTree[tID].getRootID();
                    for (int i = 0; i < n; ++i)
                        hPattern[i] += treeNodeH[tID][rootID][i];
                }
            }
#endif
            if (f.allPattern[firstTree].useDAG()) {
                if (orbitType == 0) {
                    if (fID == 0) {
                        result[qID][0] = result[qID][0]  * multiFactor / divideFactor;
                    }
                    else result[qID][0] += queryFactorH[qID][fID][0] * multiFactor / divideFactor;
                } else if (orbitType == 1) {
                    if (fID == 0) {
                        for (int i = 0; i < n; ++i)
                            result[qID][i] = result[qID][i]  * multiFactor / divideFactor;
                    }
                    else {
                        for (int i = 0; i < n; ++i)
                            result[qID][i] += queryFactorH[qID][fID][i]  * multiFactor / divideFactor;
                    }
                } else {
                    if (f.allPattern[firstTree].u.isEOrbitDir()) {
                        if (fID == 0) {
                            for (int i = 0; i < m / 2; ++i)
                                result[qID][i] = result[qID][i]  * multiFactor / divideFactor;
                        }
                        else {
                            for (int i = 0; i < m / 2; ++i)
                                result[qID][i] += queryFactorH[qID][fID][i]  * multiFactor / divideFactor;
                        }
                    }
                    else {
                        if (fID == 0) {
                            for (int i = 0; i < m; ++i)
                                result[qID][i] = result[qID][i]  * multiFactor / divideFactor;
                        }
                        else {
                            for (int i = 0; i < m; ++i)
                                result[qID][i] += queryFactorH[qID][fID][i]  * multiFactor / divideFactor;
                        }
                    }
                }
            }
            else {
                if (orbitType == 0) {
                    result[qID][0] = result[qID][0]  * multiFactor;
                } else if (orbitType == 1) {
                    for (int i = 0; i < n; ++i)
                        result[qID][i] = result[qID][i]  * multiFactor;
                } else {
                    if (f.allPattern[firstTree].u.isEOrbitDir()) {
                        for (int i = 0; i < m / 2; ++i)
                            result[qID][i] = result[qID][i]  * multiFactor;
                    }
                    else {
                        for (int i = 0; i < m; ++i)
                            result[qID][i] = result[qID][i]  * multiFactor;
                    }
                }
            }
        }
    }
    for (int tableID = 0; tableID < keyPoses.size(); ++tableID) {
        if (keyPoses[tableID] != nullptr)
            delete[] keyPoses[tableID];
    }
}