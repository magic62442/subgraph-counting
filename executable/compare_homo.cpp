// baseline type: 0: homo + prefix, 1: homo + prefix, total homo count 2: tiso, 3: scope, 5: homo, no prefix

#include "command.h"
#include "execution.h"

void executeNodeTHomo(
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
                    if (orbitType == 0) h[0] += cnt;
                    else if (orbitType == 1) {
                        for (int j = 0; j < aggreV.size(); ++j) {
//                            assert(patternV[aggrePos[j]] == aggreV[j]);
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
//                        assert(patternV[aggrePos[0]] == tau.key[0]);
                    h[dataV[aggrePos[0]]] += cnt;
                    if (keyPosSize < sizeBound) {
                        keyPos[keyPosSize] = dataV[aggrePos[0]];
                        ++keyPosSize;
                    }
                }
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
        }
    }
}

void executeNodeEdgeKeyTHomo(
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
                    if (orbitType == 0) h[0] += cnt;
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
        }
    }
}

void executePartitionHomo(
        VertexID pID,
        const Tree &t,
        VertexID ***candidate,
        ui **candCount,
        HashTable *H,
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
    ui n = dout.getNumVertices(), m = dout.getNumEdges();
    int mappingSize = 0;
    const std::vector<std::vector<VertexID>> &allChild = t.getChild();
    ui argument = 0;
    if (partitionOrder.empty()) {
        for (int i = startPos; i < endPos; ++i) {
            VertexID nID = postOrder[i];
            isRoot = endPos == postOrder.size() && i == endPos - 1;
            if (!t.nodeEdgeKey(nID)) {
                executeNodeTHomo(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                             startOffset, patternV, dataV, 0, pos, nullptr, argument, argument, tmp, allV);

            }
            else {
                executeNodeEdgeKeyTHomo(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                                    startOffset, patternV, dataV, 0, pos, nullptr, argument, argument, tmp, allV);
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
            patternV[depth] = partitionOrder[depth];
            dataV[depth] = v;
            for (VertexID nID: nodesAtStep[pID][depth]) {
                if (nID == postOrder[endPos - 1]) {
                    if (!t.nodeEdgeKey(nID)) {
                        executeNodeTHomo(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                                     startOffset, patternV, dataV, mappingSize + 1, pos, nullptr, argument, argument, tmp, allV);
                    }

                    else {
                        executeNodeEdgeKeyTHomo(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, isRoot, outID, unID, reverseID,
                                            startOffset, patternV, dataV, mappingSize + 1, pos, nullptr, argument, argument, tmp, allV);
                    }
                }
                else if (t.getNode(nID).keySize < 2) {
                    if (keyPosSize[nID] < n / 8) {
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
                    if (!t.nodeEdgeKey(nID)) {
                        executeNodeTHomo(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, false, outID, unID, reverseID,
                                     startOffset, patternV, dataV, mappingSize + 1, pos, keyPos[nID], keyPosSize[nID], sizeBound, tmp, allV);
                    }
                    else {
                        executeNodeEdgeKeyTHomo(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, false, outID, unID, reverseID,
                                            startOffset, patternV, dataV, mappingSize + 1,  pos, keyPos[nID], keyPosSize[nID], sizeBound, tmp, allV);
                    }
                }
                else {
                    if (keyPosSize[nID] < m / 8) {
                        for (int j = 0; j < keyPosSize[nID]; ++j) {
                            H[nID][keyPos[nID][j]] = 0;
                        }
                    }
                    else
                        memset(H[nID], 0, sizeof(Count) * 2 * m);
                    keyPosSize[nID] = 0;
                    executeNodeEdgeKeyTHomo(nID, t, allChild[nID], candidate[nID], candCount[nID], H, din, dout, dun, tri, p, false, outID, unID, reverseID, startOffset,
                                        patternV, dataV, mappingSize + 1, pos, keyPos[nID], keyPosSize[nID], m / 8, tmp, allV);
                }
            }
            if (depth != partitionOrder.size() - 1) {
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
            }
        }
        --depth;
        --mappingSize;
    }
    for (int j = 0; j < partitionOrder.size(); ++j) {
        if (partitionCandPos[j])
            delete[] partitionCandidate[j];
    }
    delete[] partitionCandidate;
    delete[] partitionCandCount;
    for (int j = startPos; j < endPos; ++j) {
        VertexID nID = postOrder[j];
        if (t.getNode(nID).keySize != 0)
            delete[] keyPos[postOrder[j]];
    }
    delete[] keyPos;
    delete[] keyPosSize;
}

void executeTreeHomo (
        const Tree &t,
        const DataGraph &din,
        const DataGraph &dout,
        const DataGraph &dun,
        const Triangle &tri,
        const Pattern &p,
        HashTable *H,
        EdgeID *outID,
        EdgeID *unID,
        EdgeID *reverseID,
        EdgeID *startOffset,
        VertexID *patternV,
        VertexID *dataV,
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
        executePartitionHomo(pID, t, candidate, candCount, H, din, dout, dun, tri,
                         p, outID, unID, reverseID, startOffset, patternV, dataV, pos, tmp, allV);
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

void multiJoinWrapperHomo(VertexID nID, const Tree &t, const std::vector<VertexID> &child, VertexID ***candidates,
        ui **candCounts, HashTable *H, const DataGraph &din, const DataGraph &dout, const DataGraph &dun,
        bool useTriangle, const Triangle &tri, const Pattern &p, EdgeID *outID, EdgeID *unID, EdgeID *reverseID,
        EdgeID **startOffsets, VertexID **patternVs, VertexID **dataVs, ui **poses, ui **keyPoses, ui *keyPosSizes,
        ui *sizeBounds, VertexID *&tmp, VertexID *allV
);

void multiJoinHomo(
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
            multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            for (int i = 0; i < nodesAtStep[mappingSize].size(); ++i) {
                VertexID cID = nodesAtStep[mappingSize][i];
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    patternVs[cID][j] = patternV[prefixPos[cID][j]];
                    dataVs[cID][j] = dataV[prefixPos[cID][j]];
                    startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                    poses[cID][j] = pos[prefixPos[cID][j]];
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
                multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
    }
}

void multiJoinTHomo(
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
            multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
            patternV[mappingSize] = tau.nodeOrder[mappingSize];
            dataV[mappingSize] = v;
            for (int i = 0; i < nodesAtStep[mappingSize].size(); ++i) {
                VertexID cID = nodesAtStep[mappingSize][i];
                for (int j = 0; j < prefixPos[cID].size(); ++j) {
                    patternVs[cID][j] = patternV[prefixPos[cID][j]];
                    dataVs[cID][j] = dataV[prefixPos[cID][j]];
                    startOffsets[cID][j] = startOffset[prefixPos[cID][j]];
                    poses[cID][j] = pos[prefixPos[cID][j]];
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
                multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
    }
}

void multiJoinEHomo(
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
            multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
                multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
    }
}

void multiJoinETHomo(
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
            multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                             reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
                multiJoinWrapperHomo(cID, t, t.getChild()[cID], candidates, candCounts, H, din, dout, dun, true, tri, p, outID, unID,
                                 reverseID, startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
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
    }
}

void multiJoinWrapperHomo(
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
            multiJoinHomo(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                      startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
        else {
            multiJoinEHomo(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                       startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
    }
    else {
        if (!edgeKey) {
            multiJoinTHomo(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                       startOffsets, patternVs, dataVs, poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
        else {
            multiJoinETHomo(nID, t, child, candidates, candCounts, H, din, dout, dun, tri, p, outID, unID, reverseID,
                        startOffsets, patternVs, dataVs,  poses, keyPoses, keyPosSizes, sizeBounds, tmp, allV);
        }
    }
}

void multiJoinTreeHomo(
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
            multiJoinWrapperHomo(nID, t, t.getChild()[nID], candidate, candCount, H, din, dout, dun, useTriangle, tri, p,
                             outID, unID, reverseID, startOffset, patternV, dataV, pos, keyPos, keyPosSize, sizeBounds, tmp, allV);
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

void executeSharedNodeHomo(
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
                        if (orbitType == 0) h[0] += cnt;
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
        }
    }
}

void executeSharedNodeEdgeKeyHomo(
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
                        if (orbitType == 0) h[0] += cnt;
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
        --depth;
        if (depth >= 0) {
            --mappingSize;
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

void executeSharedHomo(
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
        ui *pos,
        VertexID *&tmp,
        VertexID *allV
) {
    if (!param.rep.edgeKey) {
        executeSharedNodeHomo(param.rep, param.children, param.nIDs, param.trees, param.hashTables,
                           param.aggreWeights, param.id2AggreKey, param.id2ChildKey, param.isRoot,
                           param.childKeyPoses, param.aggrePoses, param.tableIDs, keyPoses, keyPosSize, sizeBound,
                           param.nodeInterPos, param.nodeInPos, param.nodeOutPos, param.nodeUnPos, param.nodeCandPos,
                           param.greaterPos, param.lessPos, param.nodeTriPos, param.triEdgeType, param.triEndType,
                           param.nonLeafIDs, param.candidate, param.candCount, din, dout, dun, tri, param.p,
                           outID, unID, reverseID, startOffset, patternV, dataV, mappingSize, pos, tmp, allV);
    }
    else {
        executeSharedNodeEdgeKeyHomo(param.rep, param.children, param.nIDs, param.trees, param.hashTables,
                                  param.aggreWeights, param.id2AggreKey, param.id2ChildKey, param.isRoot,
                                  param.childKeyPoses, param.aggrePoses, param.tableIDs, keyPoses, keyPosSize, sizeBound,
                                  param.nodeInterPos, param.nodeInPos, param.nodeOutPos, param.nodeUnPos, param.nodeCandPos,
                                  param.greaterPos, param.lessPos, param.posChildEdge, param.posAggreEdge, param.childEdgeType,
                                  param.aggreEdgeType, param.nodeTriPos, param.nonLeafIDs, param.triEdgeType,
                                  param.triEndType, param.candidate, param.candCount, din, dout, dun, tri, param.p,
                                  outID, unID, reverseID, startOffset, patternV, dataV, mappingSize, pos, tmp, allV);
    }
}
void executeForestHomo(Forest &f, std::vector<HashTable> &result, const DataGraph &din,
                   const DataGraph &dout, const DataGraph &dun, bool useTriangle, const Triangle &tri,
                   std::vector<HashTable> &allocatedHashTable, EdgeID *outID, EdgeID *unID, EdgeID *reverseID,
                   EdgeID *startOffset, VertexID *patternV, VertexID *dataV, ui *candPos,
                   VertexID *&tmp, VertexID *allV) {
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
        if (rep.prefixSize != 0 && rep.keySize == 1) {
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
        executeSharedHomo(param, keyPoses, keyPosSize, sizeBound, din, dout, dun, useTriangle, tri, outID, unID, reverseID,
                      startOffset, patternV, dataV, 0, candPos, tmp, allV);
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
            const std::vector<bool> &partitionCandPos = firstTree.getPartitionCandPos(gPID[0]);
            const std::vector<std::vector<int>> &partitionInPos = firstTree.getPartitionInPos(gPID[0]);
            const std::vector<std::vector<int>> &partitionOutPos = firstTree.getPartitionOutPos(gPID[0]);
            const std::vector<std::vector<int>> &partitionUnPos = firstTree.getPartitionUnPos(gPID[0]);
            const std::vector<bool> &partitionInterPos = firstTree.getPartitionInterPos(gPID[0]);
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
                            else if (sizeBound[tableID] == m / 8 + 1) {
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
                        executeSharedHomo(param, keyPoses, keyPosSize, sizeBound, din, dout, dun, useTriangle, tri, outID, unID, reverseID,
                                      startOffset, patternV, dataV, mappingSize + 1, candPos, tmp, allV);
                    }
                    if (depth != length - 1) {
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
                    }
                }
                --depth;
                --mappingSize;
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
        executeSharedHomo(param, keyPoses, keyPosSize, sizeBound, din, dout, dun, useTriangle, tri, outID, unID, reverseID,
                      startOffset, patternV, dataV, 0, candPos, tmp, allV);
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

std::vector<HashTable>
baseline(const std::vector<PatternGraph> &patternGraphs, const DataGraph &din, const DataGraph &dout,
         const DataGraph &dun, const Triangle &triangle, double &totalPlanTime, double &totalExeTime,
         EdgeID *outID, EdgeID *unID, EdgeID *reverseID, EdgeID **startOffset, VertexID **patternV, VertexID **dataV,
         bool **visited, ui *candPos, VertexID *&tmp, VertexID *allV, int type, specialsparse *sg,
         VertexID *cliqueVertices, const std::vector<std::string> &files, Count &numMatch) {
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    HashTable factorSum = new Count[m + 1];
    HashTable ht[MAX_NUM_NODE];
    for (auto & h: ht)
        h = new Count[m + 1];
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::vector<HashTable> result(patternGraphs.size());
    int patternSize = patternGraphs[0].getNumVertices();
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
        }
        else {
            if (pg.isEOrbitDir()) {
                result[i] = new Count[m / 2];
                memset(result[i], 0, sizeof(Count) * m / 2);
            }
            else {
                result[i] = new Count[m];
                memset(result[i], 0, sizeof(Count) * m);
            }
        }
    }
    std::vector<HashTable> mathCalH(patternGraphs.size());
    for (int i = 0; i < patternGraphs.size(); ++i) {
        numMatch += gNumMatch;
//        gNumMatch = 0;
        const PatternGraph &pg = patternGraphs[i];
        int orbitType = pg.getOrbitType();
        HashTable H;
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
        std::map<int, std::vector<Pattern>> patterns;
        std::map<int, std::vector<std::vector<Tree>>> trees;
        ConNode cn;
        start = std::chrono::steady_clock::now();
        if (type == 1)
            homoEquation(pg, patterns, trees, true);
        else if (type == 2)
            genEquation(pg, patterns, trees, cn, true, false, true);
        else
            genEquation(pg, patterns, trees, cn, true, true, true);
        end = std::chrono::steady_clock::now();
        elapsedSeconds = end - start;
        totalPlanTime += elapsedSeconds.count();
        start = std::chrono::steady_clock::now();
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
            executeConNode(cn, din, dout, dun, true, triangle, p, outID, unID, startOffset[0], patternV[0],
                           dataV[0], visited[0], candPos, tmp, allV);
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
        for (auto it = patterns.begin(); it != patterns.end(); ++it) {
            int divideFactor = it->first;
            memset(factorSum, 0, sizeof(Count) * (m + 1));
            for (int j = 0; j < trees[divideFactor].size(); ++j) {
                for (int j2 = 0; j2 < trees[divideFactor][j].size(); ++j2) {
                    int multiFactor = trees[divideFactor][j][j2].getMultiFactor();
                    if (multiFactor == 0) continue;
                    for (int l = 0; l < trees[divideFactor][j][0].getNumNodes(); ++l) {
                        memset(ht[l], 0, sizeof(Count) * m);
                    }
                    const Tree &t = trees[divideFactor][j][j2];
                    if (type == 1) {
                        if (t.getExecuteMode()) {
                            executeTreeHomo(t, din, dout, dun, triangle, patterns[divideFactor][j],
                                            ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], candPos, tmp, allV);
                        }
                        else {
                            multiJoinTreeHomo(t, din, dout, dun, true, triangle, patterns[divideFactor][j],
                                              ht, outID, unID, reverseID, startOffset, patternV, dataV, tmp, allV);
                        }
                    }
                    else {
                        if (t.getExecuteMode()) {
                            executeTree(t, din, dout, dun, true, triangle, patterns[divideFactor][j],
                                        ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
                        }
                        else {
                            multiJoinTree(t, din, dout, dun, true, triangle, patterns[divideFactor][j],
                                          ht, outID, unID, reverseID, startOffset, patternV, dataV, visited, tmp, allV);
                        }
                    }
                    HashTable h = ht[trees[divideFactor][j][j2].getRootID()];
                    if (orbitType == 0) factorSum[0] += h[0];
                    else if (orbitType == 1)
                        for (VertexID l = 0; l < n; ++l) {
                            factorSum[l] += h[l] * multiFactor;
                        }
                    else
                        for (EdgeID l = 0; l < m + 1; ++l) {
                            factorSum[l] += h[l] * multiFactor;
                        }
                    if (j == 0 && type == 1) {
                        for (VertexID l = 0; l < n; ++l) {
                            H[l] = factorSum[l];
                        }
                    }
                }
            }
            if (type == 2) {
                if (orbitType == 0) H[0] += factorSum[0];
                else if (orbitType == 1)
                    for (VertexID l = 0; l < n; ++l) {
                        H[l] += factorSum[l];
                    }
                else
                    for (EdgeID l = 0; l < m + 1; ++l) {
                        H[l] += factorSum[l];
                    }
            }

        }
        end = std::chrono::steady_clock::now();
        elapsedSeconds = end - start;
        const std::vector<std::vector<Tree>> &allTree = trees.begin()->second;
        ui maxTW = 0, length = 0;
        double maxFHW = 0.0;
        int maxTWCnt = 0, maxFHWCnt = 0;
//        std::cout << files[i] << " " << elapsedSeconds.count() << " ";
        for (const auto & ts : allTree) {
            if (ts[0].getTreeWidth() > maxTW) maxTW = ts[0].getTreeWidth();
            if (ts[0].getFhw() > maxFHW) maxFHW = ts[0].getFhw();
        }
        for (const auto & ts : allTree) {
            if (ts[0].getMultiFactor() != 0) ++length;
            if (ts[0].getTreeWidth() == maxTW) ++maxTWCnt;
            if (ts[0].getFhw() == maxFHW) ++maxFHWCnt;
        }
//        std::cout << length << " " << maxFHW << " " << maxFHWCnt << " " << maxTW << " " << maxTWCnt << " " << gNumMatch << std::endl;
        totalExeTime += elapsedSeconds.count();
        result[i] = H;
    }

//    for (int i = 0; i < conIDs.size(); ++i) {
//        int divideFactor = conFactors[i];
//        int orbitType = patternGraphs[i].getOrbitType();
//        bool eOrbitDir = patternGraphs[i].isEOrbitDir();
//        int id = conIDs[i];
//        if (orbitType == 0) result[id][0] /= divideFactor;
//        else if (orbitType == 1) {
//            for (int k = 0; k < n; ++k)
//                result[id][k] /= divideFactor;
//        }
//        else {
//            if (eOrbitDir) {
//                for (int k = 0; k < m / 2; ++k)
//                    result[id][k] /= divideFactor;
//            }
//            else {
//                for (int k = 0; k < m; ++k)
//                    result[id][k] /= divideFactor;
//            }
//        }
//    }

    return result;
}

// collect statistics: total number of iso-matches, homo-matches, tiso-matches
// homo: num trees, num tree nodes, tree node size distribution

void
collect(const std::vector<PatternGraph> &patternGraphs, const DataGraph &din, const DataGraph &dout,
         const DataGraph &dun, const Triangle &triangle, double &totalPlanTime, double &totalExeTime,
         EdgeID *outID, EdgeID *unID, EdgeID *reverseID, EdgeID **startOffset, VertexID **patternV, VertexID **dataV,
         bool **visited, ui *candPos, VertexID *&tmp, VertexID *allV, int type, specialsparse *sg,
         VertexID *cliqueVertices, const std::vector<std::string> &files, Count &numMatch) {
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    const PatternGraph &pg = patternGraphs[0];
    int orbitType = pg.getOrbitType();
    HashTable H = new Count[dun.getNumVertices()];
        memset(H, 0, sizeof(Count) * n);
    HashTable ht[MAX_NUM_NODE];
    for (auto & h: ht)
        h = new Count[m + 1];
    Count isoMatch = 0, tisoMatch = 0, homoMatch = 0;
    ui homoNumTree = 0, tisoNumTree = 0;
    std::map<ui, int> homoNodeSize, tisoNodeSize;
    std::map<int, std::vector<Pattern>> patterns;
    std::map<int, std::vector<std::vector<Tree>>> trees;
    ConNode cn;
    homoEquation(pg, patterns, trees, true);
    int divideFactor = patterns.begin()->first;
    Tree t = trees[divideFactor][0][0];
    int multiFactor = trees[divideFactor][0][0].getMultiFactor();
    for (int l = 0; l < t.getNumNodes(); ++l) {
        memset(ht[l], 0, sizeof(Count) * m);
    }
    if (t.getExecuteMode()) {
        executeTreeHomo(t, din, dout, dun, triangle, patterns[divideFactor][0],
                        ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], candPos, tmp, allV);
    }
    else {
        multiJoinTreeHomo(t, din, dout, dun, true, triangle, patterns[divideFactor][0],
                          ht, outID, unID, reverseID, startOffset, patternV, dataV, tmp, allV);
    }
    HashTable h = ht[t.getRootID()];
    for (int i = 0; i < n; ++i) {
        homoMatch += h[i] * multiFactor ;
    }
    homoNumTree = trees.begin()->second.size();
    for (int j = 0; j < trees[divideFactor].size(); ++j) {
        for (int j2 = 0; j2 < trees[divideFactor][j].size(); ++j2) {
            t = trees[divideFactor][j][j2];
            for (int i = 0; i < t.getNumNodes(); ++i) {
                const Node &tau = t.getNode(i);
                if (homoNodeSize.find(tau.numVertices) == homoNodeSize.end())
                    homoNodeSize[tau.numVertices] = 1;
                else homoNodeSize[tau.numVertices] += 1;
            }
        }
    }
    cn.num = 0;
    patterns.clear();
    trees.clear();
    genEquation(pg, patterns, trees, cn, true, true, true);
    divideFactor = patterns.begin()->first;
    t = trees[divideFactor][0][0];
    multiFactor = trees[divideFactor][0][0].getMultiFactor();
    for (int l = 0; l < t.getNumNodes(); ++l) {
        memset(ht[l], 0, sizeof(Count) * m);
    }
    if (t.getExecuteMode()) {
        executeTreeHomo(t, din, dout, dun, triangle, patterns[divideFactor][0],
                        ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], candPos, tmp, allV);
    }
    else {
        multiJoinTreeHomo(t, din, dout, dun, true, triangle, patterns[divideFactor][0],
                          ht, outID, unID, reverseID, startOffset, patternV, dataV, tmp, allV);
    }
    h = ht[t.getRootID()];
    for (int i = 0; i < n; ++i) {
        tisoMatch += h[i] * multiFactor ;
    }
    tisoNumTree = trees.begin()->second.size();
    HashTable factorSum = new Count[m + 1];
    memset(factorSum, 0, sizeof(Count) * (m + 1));
    for (int j = 0; j < trees[divideFactor].size(); ++j) {
        for (int j2 = 0; j2 < trees[divideFactor][j].size(); ++j2) {
            t = trees[divideFactor][j][j2];
            for (int i = 0; i < t.getNumNodes(); ++i) {
                const Node &tau = t.getNode(i);
                if (tisoNodeSize.find(tau.numVertices) == tisoNodeSize.end())
                    tisoNodeSize[tau.numVertices] = 1;
                else tisoNodeSize[tau.numVertices] += 1;
            }
            for (int l = 0; l < t.getNumNodes(); ++l) {
                memset(ht[l], 0, sizeof(Count) * m);
            }
            if (t.getExecuteMode()) {
                executeTree(t, din, dout, dun, true, triangle, patterns[divideFactor][j],
                            ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], visited[0], candPos, tmp, allV);
            }
            else {
                multiJoinTree(t, din, dout, dun, true, triangle, patterns[divideFactor][j],
                              ht, outID, unID, reverseID, startOffset, patternV, dataV, visited, tmp, allV);
            }
            h = ht[trees[divideFactor][j][j2].getRootID()];
            multiFactor = t.getMultiFactor();
            for (VertexID l = 0; l < n; ++l) {
                factorSum[l] += h[l] * multiFactor;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        isoMatch += factorSum[i] ;
    }
    std::cout << isoMatch << " " << tisoMatch << " " << homoMatch << std::endl;
    std::cout << homoNumTree << " " << tisoNumTree << std::endl;
    for (auto it = homoNodeSize.begin(); it != homoNodeSize.end(); ++it) {
        std::cout << it->first << " " << it->second << std::endl;
    }
    std::cout << "------" << std::endl;
    for (auto it = tisoNodeSize.begin(); it != tisoNodeSize.end(); ++it) {
        std::cout << it->first << " " << it->second << std::endl;
    }
}

int main(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string queryGraphPath = cmd.getQueryGraphPath();
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    std::string trianglePath = cmd.getTrianglePath();
    int baselineType = cmd.getBaselineType();
    bool batchQuery = cmd.getBatchQuery();
    std::vector<std::string> files;
    std::vector<PatternGraph> patternGraphs = loadPatternGraph(queryGraphPath, batchQuery, files);
    DataGraph dun = DataGraph();
    dun.loadDataGraph(dataGraphPath);
    const DataGraph din = constructDirectedDataGraph(dun, false);
    const DataGraph dout = constructDirectedDataGraph(dun, true);
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    Triangle triangle;
    triangle.load(trianglePath, m / 2);
    HashTable ht[MAX_NUM_NODE];
    for (auto & h: ht)
        h = new Count[m + 1];
    EdgeID *outID = buildInID2OutID(din, dout);
    EdgeID *unID = buildUnID2OutID(dun, dout);
    EdgeID *reverseID = buildReverseUnID(dun);
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
    specialsparse *sg = (specialsparse *)malloc(sizeof(specialsparse));
    dout.initSpecialSparse(sg);
    VertexID *cliqueVertices = new VertexID[MAX_PATTERN_SIZE + 1];
    double totalPlanTime = 0.0, totalExeTime = 0.0;
    std::vector<HashTable> result = std::vector<HashTable>(patternGraphs.size(), nullptr);
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
        }
        else {
            if (pg.isEOrbitDir()) {
                result[i] = new Count[m / 2];
                memset(result[i], 0, sizeof(Count) * m / 2);
            }
            else {
                result[i] = new Count[m];
                memset(result[i], 0, sizeof(Count) * m);
            }
        }
    }
    Count numMatch = 0;
//    bool prefix = true;
//    if (baselineType > 3) prefix = false;
//    if (baselineType == 0 || baselineType == 5) {
//        int patternSize = patternGraphs[0].getNumVertices();
//        bool forestShare = patternSize <= 5;
//        auto start = std::chrono::steady_clock::now();
//        auto end = std::chrono::steady_clock::now();
//        std::chrono::duration<double> elapsedSeconds = end - start;
//        Forest fu, fuShrink;
//        std::vector<CanonType> shrinkCanon;
//        std::vector<int> shrinkOrbit;
//        std::vector<int> shrinkRootOrbit;
//        std::vector<HashTable> shrinkageH;
//        std::vector<HashTable> uH;
//        std::vector<std::vector<HashTable>> dH;
//        // for the ith Pattern, jth shrinkage, its position in the shrinkageH and the factor
//        std::vector<std::vector<int>> shrinkPos(patternGraphs.size());
//        // for the ith pattern, position j, the multiply factor of that shrinkage
//        std::vector<std::vector<int>> shrinkMultiFactors(patternGraphs.size());
//        std::vector<int> divideFactors(patternGraphs.size());
//        std::vector<bool> isUndirected(patternGraphs.size(), false);
//        std::vector<HashTable> allocatedHashTable;
//        allocatedHashTable.reserve(10000);
//        for (int i = 0; i < patternGraphs.size(); ++i) {
//            const PatternGraph &pg = patternGraphs[i];
//            int orbitType = pg.getOrbitType();
//            start = std::chrono::steady_clock::now();
//            std::map<int, std::vector<Pattern>> patterns;
//            std::map<int, std::vector<std::vector<Tree>>> trees;
//            ConNode cn;
//            start = std::chrono::steady_clock::now();
//            homoEquation(pg, patterns, trees, prefix);
//            isUndirected[i] = true;
//            // put the first pattern to fu, the remaining patterns to fShrinkage
//            int divideFactor = patterns.begin() -> first;
//            divideFactors[i] = divideFactor;
//            const std::vector<Pattern> &allPattern = patterns.begin()->second;
//            std::vector<std::vector<Tree>> &allTree = trees.begin()->second;
//            std::map<int, std::vector<Pattern>> rootPattern;
//            std::map<int, std::vector<std::vector<Tree>>> rootTree;
//            rootPattern[divideFactor].push_back(allPattern[0]);
//            rootTree[divideFactor].push_back(allTree[0]);
//            shrinkPos[i] = std::vector<int>(allPattern.size() - 1);
//            shrinkMultiFactors[i] = std::vector<int>(allPattern.size() - 1 + shrinkCanon.size());
//            for (int j = 1; j < allPattern.size(); ++j) {
//                bool exists = false;
//                int rootOrbit = allTree[j][0].getNode(allTree[j][0].getRootID()).v2o[0];
//                for (int k = 0; k < shrinkCanon.size(); ++k) {
//                    if (allPattern[j].u.getCanonValue() == shrinkCanon[k] && allPattern[j].u.getOrbit(0) == shrinkOrbit[k]
//                        && rootOrbit == shrinkRootOrbit[k]) {
//                        shrinkPos[i][j - 1] = k;
//                        shrinkMultiFactors[i][k] = allTree[j][0].getMultiFactor();
//                        exists = true;
//                        break;
//                    }
//                }
//                if (!exists) {
//                    HashTable h;
//                    if (orbitType == 0) {
//                        h = new Count[1];
//                        h[0] = 0;
//                    }
//                    else if (orbitType == 1) {
//                        h = new Count[n];
//                        memset(h, 0, sizeof(Count) * n);
//                    }
//                    else {
//                        if (pg.isEOrbitDir()) {
//                            h = new Count[m / 2];
//                            memset(h, 0, sizeof(Count) * m / 2);
//                        }
//                        else {
//                            h = new Count[m];
//                            memset(h, 0, sizeof(Count) * m);
//                        }
//                    }
//                    shrinkPos[i][j - 1] = shrinkCanon.size();
//                    shrinkMultiFactors[i][shrinkCanon.size()] = allTree[j][0].getMultiFactor();
//                    allTree[j][0].setMultiFactor(1);
//                    shrinkCanon.push_back(allPattern[j].u.getCanonValue());
//                    shrinkOrbit.push_back(allPattern[j].u.getOrbit(0));
//                    shrinkRootOrbit.push_back(rootOrbit);
//                    std::map<int, std::vector<Pattern>> shrinkPattern;
//                    std::map<int, std::vector<std::vector<Tree>>> shrinkTree;
//                    shrinkPattern[1].push_back(allPattern[j]);
//                    shrinkTree[1].push_back(allTree[j]);
//                    fuShrink.loadQuery(shrinkPattern, shrinkTree, forestShare);
//                    shrinkageH.push_back(h);
//                }
//            }
//            uH.push_back(result[i]);
//            fu.loadQuery(rootPattern, rootTree, forestShare);
//            end = std::chrono::steady_clock::now();
//            std::chrono::duration<double> elapsedSeconds = end - start;
//            totalPlanTime += elapsedSeconds.count();
//        }
//#ifndef ONLY_PLAN
//        // 2. execute undirected shrinkCanon, undirected root patterns and directed patterns
//        if (patternSize <= 5) {
//            start = std::chrono::steady_clock::now();
//            executeForestHomo(fuShrink, shrinkageH, din, dout, dun, true, triangle, allocatedHashTable,
//                          outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], candPos, tmp, allV);
//            for (HashTable h: allocatedHashTable) delete[] h;
//            allocatedHashTable.clear();
//            end = std::chrono::steady_clock::now();
//            elapsedSeconds = end - start;
//            totalExeTime += elapsedSeconds.count();
//            std::cout << "finished undirected shrinkages. time: " << totalExeTime << std::endl;
//            start = std::chrono::steady_clock::now();
//            executeForestHomo(fu, uH,  din, dout, dun, true, triangle, allocatedHashTable,
//                          outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], candPos, tmp, allV);
//            for (HashTable h: allocatedHashTable) delete[] h;
//            allocatedHashTable.clear();
//            end = std::chrono::steady_clock::now();
//            elapsedSeconds = end - start;
//            totalExeTime += elapsedSeconds.count();
//            std::cout << "finished undirected patterns. time: " << totalExeTime << std::endl;
//        }
//        else {
//            start = std::chrono::steady_clock::now();
//            for (int i = 0; i < fuShrink.numQuery; ++i) {
//                const Tree &t = fuShrink.allTree[i];
//                const Pattern &p = fuShrink.allPattern[i];
//                for (int l = 0; l < t.getNumNodes(); ++l) memset(ht[l], 0, sizeof(Count) * m);
//                if (t.getExecuteMode())
//                    executeTreeHomo(t, din, dout, dun, triangle, p, ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], candPos, tmp, allV);
//                else
//                    multiJoinTreeHomo(t, din, dout, dun, true, triangle, p, ht, outID, unID, reverseID, startOffset, patternV, dataV, tmp, allV);
//                HashTable h = ht[t.getRootID()];
//                for (int j = 0; j < n; ++j) {
//                    shrinkageH[i][j] = h[j];
//                }
//            }
//            end = std::chrono::steady_clock::now();
//            elapsedSeconds = end - start;
//            totalExeTime += elapsedSeconds.count();
//            std::cout << "finished undirected shrinkages. time: " << totalExeTime << std::endl;
//            start = std::chrono::steady_clock::now();
//            for (int i = 0; i < fu.numQuery; ++i) {
//                const Tree &t = fu.allTree[i];
//                const Pattern &p = fu.allPattern[i];
//                for (int l = 0; l < t.getNumNodes(); ++l) memset(ht[l], 0, sizeof(Count) * m);
//                if (t.getExecuteMode())
//                    executeTreeHomo(t, din, dout, dun, triangle, p, ht, outID, unID, reverseID, startOffset[0], patternV[0], dataV[0], candPos, tmp, allV);
//                else
//                    multiJoinTreeHomo(t, din, dout, dun, true, triangle, p, ht, outID, unID, reverseID, startOffset, patternV, dataV, tmp, allV);
//                HashTable h = ht[t.getRootID()];
//                int multiFactor = t.getMultiFactor();
//                for (int j = 0; j < n; ++j) {
//                    uH[i][j] = h[j] * multiFactor;
//                }
//            }
//            end = std::chrono::steady_clock::now();
//            elapsedSeconds = end - start;
//            totalExeTime += elapsedSeconds.count();
//            std::cout << "finished undirected patterns. time: " << totalExeTime << std::endl;
//        }
//
//        // 3. subtract shrinkage count for undirected patterns
//        for (int i = 0; i < patternGraphs.size(); ++i) {
//            if (!isUndirected[i]) continue;
//            const PatternGraph &pg = patternGraphs[i];
//            int orbitType = pg.getOrbitType();
//            int divideFactor = divideFactors[i];
//            for (int j = 0; j < shrinkPos[i].size(); ++j) {
//                HashTable h = shrinkageH[shrinkPos[i][j]];
//                if (orbitType == 0)
//                    result[i][0] += h[0] * shrinkMultiFactors[i][shrinkPos[i][j]];
//                else if (orbitType == 1) {
//                    for (int k = 0; k < n; ++k)
//                        result[i][k] += h[k] * shrinkMultiFactors[i][shrinkPos[i][j]];
//                }
//                else if (orbitType == 2 && pg.isEOrbitDir()) {
//                    for (int k = 0; k < m / 2; ++k)
//                        result[i][k] += h[k] * shrinkMultiFactors[i][shrinkPos[i][j]];
//                }
//                else {
//                    for (int k = 0; k < m; ++k)
//                        result[i][k] += h[k] * shrinkMultiFactors[i][shrinkPos[i][j]];
//                }
//            }
//            if (orbitType == 0)
//                result[i][0] /= divideFactor;
//            else if (orbitType == 1) {
//                for (int k = 0; k < n; ++k)
//                    result[i][k] /= divideFactor;
//            }
//            else if (orbitType == 2 && pg.isEOrbitDir()) {
//                for (int k = 0; k < m / 2; ++k)
//                    result[i][k] /= divideFactor;
//            }
//            else {
//                for (int k = 0; k < m; ++k)
//                    result[i][k] /= divideFactor;
//            }
//        }
//#endif
//    }
//    else  {
//        result = baseline(patternGraphs, din, dout, dun, triangle, totalPlanTime, totalExeTime, outID,
//                          unID, reverseID, startOffset, patternV, dataV, visited, candPos, tmp, allV, baselineType,
//                          sg, cliqueVertices, files, numMatch);
//    }
////    std::cout << "total planning time: " << totalPlanTime << ", total execution time: " << totalExeTime
////              << ", total time: " << totalPlanTime + totalExeTime << std::endl;
////    std::cout << "number of match: " << gNumMatch << ", number of intermediate: " << gNumIntermediate << std::endl;
//    std::cout << gNumIntermediate << " " << gNumMatch << " " << totalExeTime << " " << totalPlanTime << std::endl;
//    int orbitType = patternGraphs[0].getOrbitType();
//    std::vector<int> orbitTypes(files.size(), orbitType);
//    if (!resultPath.empty()) saveCount(resultPath, result, dun, batchQuery, files, orbitTypes);
    collect(patternGraphs, din, dout, dun, triangle, totalPlanTime, totalExeTime, outID,
             unID, reverseID, startOffset, patternV, dataV, visited, candPos, tmp, allV, baselineType,
             sg, cliqueVertices, files, numMatch);
    return 0;
}
