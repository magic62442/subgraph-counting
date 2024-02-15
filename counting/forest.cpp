//
// Created by anonymous author on 2022/9/15.
//

#include "forest.h"

void
Forest::loadQuery(std::map<int, std::vector<Pattern>> &patterns, std::map<int, std::vector<std::vector<Tree>>> &trees,
                  bool share) {
    if (!share) {
        for (auto it = patterns.begin(); it != patterns.end(); ++it) {
            int divideFactor = it->first;
            const std::vector<Pattern> &pts = it->second;
            for (int i = 0; i < it->second.size(); ++i) {
                allPattern.push_back(patterns[divideFactor][i]);
                allTree.push_back(trees[divideFactor][i][0]);
                ++numQuery;
            }
        }
        return;
    }
    queryFactor.emplace_back();
    queryFactorFirstTree.emplace_back(patterns.size());
    queryFactorNumTree.emplace_back(patterns.size());
    int fID = 0;
    for (auto it = patterns.begin(); it != patterns.end(); ++it) {
        int divideFactor = it -> first;
        const std::vector<Pattern> &pts = it -> second;
        queryFactor[numQuery].push_back(divideFactor);
        queryFactorFirstTree[numQuery][fID] = (int)allTree.size();
        queryFactorNumTree[numQuery][fID] = 0;
        // sort trees for the same factor.
        std::vector<TreeAndPattern> tps;
        for (int i = 0; i < it->second.size(); ++i) {
            for (int j = 0; j < trees[divideFactor][i].size(); ++j)
                tps.emplace_back(trees[divideFactor][i][j], i);
        }
        std::sort(tps.begin(), tps.end());
        queryFactorNumTree[numQuery][fID] = (int)tps.size();
        for (const TreeAndPattern& tp : tps) {
            int tID = (int)allTree.size();
            const Tree &t = tp.t;
            allPattern.push_back(pts[tp.patternPos]);
            if (allTree.empty()) orbitType = t.getOrbitType();
            allTree.push_back(t);
            nodeExecutable.emplace_back(t.getNumNodes(), false);
            nodeExecuted.emplace_back(t.getNumNodes());
            nodeNeedExecute.emplace_back(t.getNumNodes());
            posInChildTables.emplace_back(t.getNumNodes());
            nodeParent.emplace_back(t.getNumNodes());
            posInMergedNode.emplace_back(t.getNumNodes());
            // deal with tree nodes
            const std::vector<std::vector<VertexID>> &children = t.getChild();
            for (VertexID nID : t.getPostOrder()) {
                const Node &tau = t.getNode(nID);
                std::vector<int> newAggrePos;
                std::vector<int> aggreWeight;
                std::vector<int> newAggreWeight;
                std::vector<std::vector<int>> newChildKeyPos;
                bool nodeExists = false;
                int mergedNodePos = (int)mergedNode.size();
                bool isRoot = false;
                if (t.getRootID() == nID) {
                    isRoot = true;
                    aggreWeight = t.getAggreWeight();
                }
                for (int k = 0; k < mergedNode.size(); ++k) {
                    if (mergeable(mergedNodeNID[k], nID, allTree[mergedNodeTID[k]], t, allPattern[mergedNodeTID[k]],
                                  pts[tp.patternPos], newAggrePos, newChildKeyPos, aggreWeight, newAggreWeight,
                                  isRoot, orbitType))  {
                        toTreeNode[k].emplace_back(numQuery, tID, nID, newAggrePos, newAggreWeight, newChildKeyPos);
                        if (mergedNode[k].localOrder.size() == tau.localOrder.size()) {
                            nodeExists = true;
                            mergedNodePos = k;
                        }
                    }
                }
                posInMergedNode[tID][nID] = mergedNodePos;
                if (!nodeExists) {
                    mergedNode.push_back(tau);
                    mergedNodeTID.push_back(tID);
                    mergedNodeNID.push_back(nID);
                    toTreeNode.emplace_back();
                    const std::vector<int> &aggrePos = t.getAggrePos(nID);
                    const std::vector<std::vector<int>> &childKeyPos = t.getChildKeyPos(nID);
                    toTreeNode.back().emplace_back(numQuery, tID, nID, aggrePos, aggreWeight, childKeyPos);
                }
//                if (t.isLeaf(nID)) nodeExecutable[tID][nID] = true;
//                else nodeExecutable[tID][nID] = false;
                nodeExecuted[tID][nID] = false;
                nodeNeedExecute[tID][nID] = false;
                posInChildTables[tID][nID] = -1;
                nodeParent[tID][nID] = t.getParent(nID);
            }
            // deal with tree partitions
            const std::vector<int> &partitionPos = t.getPartitionPos();
            const std::vector<VertexID> &postOrder = t.getPostOrder();
            nodePID.emplace_back(t.getNumNodes());
            for (VertexID pID = 0; pID < t.getGlobalOrder().size(); ++pID) {
                int startPos;
                if (pID == 0) startPos = 0;
                else startPos = partitionPos[pID - 1] + 1;
                int endPos = partitionPos[pID] + 1;
                for (int k = startPos; k < endPos; ++k) {
                    VertexID nID = postOrder[k];
                    nodePID[tID][nID] = pID;
                }
            }
            // leaf nodes of a partition is executable if for all partition nodes their children are still in this partition
            for (VertexID pID = 0; pID < t.getGlobalOrder().size(); ++pID) {
                bool partExecutable = true;
                int startPos;
                if (pID == 0) startPos = 0;
                else startPos = partitionPos[pID - 1] + 1;
                int endPos = partitionPos[pID] + 1;
                for (int k = startPos; k < endPos; ++k) {
                    VertexID nID = postOrder[k];
                    for (VertexID cID : children[nID]) {
                        if (nodePID[tID][cID] != pID) {
                            partExecutable = false;
                            break;
                        }
                    }
                    if (!partExecutable) break;
                }
                if (partExecutable) {
                    for (int k = startPos; k < endPos; ++k) {
                        VertexID nID = postOrder[k];
                        if (t.isLeaf(nID)) nodeExecutable[tID][nID] = true;
                    }
                }
            }
        }
        ++fID;
    }
    ++numQuery;
    numMergedNodes = (int)mergedNode.size();
//    for (int i = 0; i < allTree.size(); ++i) {
//        printf("tree %d:\n", i);
//        allTree[i].print();
//    }
}

bool Forest::executable(const NodeInfo &inf) const {
    return nodeNeedExecute[inf.tID][inf.nID] && !nodeExecuted[inf.tID][inf.nID] && nodeExecutable[inf.tID][inf.nID];
}

const Pattern &Forest::mergedNodePattern(int pos) const {
    int tID = mergedNodeTID[pos];
    return allPattern[tID];
}

void Forest::updateNodeExecutable(int tableID) {
    for (int i = 0; i < childTID[tableID].size(); ++i) {
        int tID = childTID[tableID][i];
        VertexID nID = childNID[tableID][i];
        nodeExecuted[tID][nID] = true;
        VertexID parentID = allTree[tID].getParent(nID);
        bool flag = true;
        if (nodeNeedExecute[tID][parentID]) {
            for (VertexID cID: allTree[tID].getChild()[parentID]) {
                if (allTree[tID].isLeaf(cID)) nodeExecutable[tID][cID] = true;
                if (!nodeExecuted[tID][cID]) {
                    nodeExecutable[tID][parentID] = false;
                    flag = false;
                    break;
                }
            }
            if (flag) nodeExecutable[tID][parentID] = true;
        }
    }
}

void Forest::allocateTables(HashTable **treeNodeH, std::vector<HashTable> &allocatedHashTable, ui n, ui m, int &nTable,
                            int &mTable) {
    for (int qID = 0; qID < numQuery; ++qID) {
        for (int fID = 0; fID < queryFactor[qID].size(); ++fID) {
            int firstTree = queryFactorFirstTree[qID][fID];
            int numTree = queryFactorNumTree[qID][fID];
            std::map<CanonType, std::vector<std::pair<int, VertexID>>> canon2RootNode;
            for (int tID = firstTree; tID < numTree + firstTree; ++tID) {
                const Tree &t = allTree[tID];
                VertexID rootID = t.getRootID();
                const Pattern &p = allPattern[tID];
                CanonType canon = p.u.getCanonValue();
                canon2RootNode[canon].emplace_back(tID, rootID);
            }
            std::queue<std::vector<std::pair<int, VertexID>>> q;
            std::queue<int> siblingQ;
            for (auto it = canon2RootNode.begin(); it != canon2RootNode.end(); ++it)
#ifdef DEBUG
                if (fID == 1 && it -> first == 282161) {
                    q.push(it->second);
                    siblingQ.push(1);
                    for (int i = 0; i < it->second.size(); ++i) {
                        int tID = it->second[i].first;
                        printf("tID = %d\n", tID);
                        printf("pattern out canon = %lld\n", allPattern[tID].out.getCanonValue());
                        allPattern[tID].out.printGraph();
                        allTree[tID].print();
                    }
                }
//                    q.push(it->second);
#else
            {
                q.push(it->second);
                siblingQ.push(1);
            }
#endif
            while (!q.empty()) {
                std::vector<std::pair<int, VertexID>> treeNodeID = q.front();
                int siblingNum = siblingQ.front();
                int setTRep = treeNodeID[0].first;
                VertexID setNRep = treeNodeID[0].second;
                bool isRoot = allTree[setTRep].getRootID() == setNRep;
                bool isPart = !allTree[setTRep].needPrefix(setNRep);
                q.pop();
                siblingQ.pop();
                std::set<PatternInfo> setPatterns;
                // group nodes in the same set, collect the child of each group.
                // assume that if the subtree canonical value is equal, then the node order is isomorphic
                // nodes are grouped together if the followings are the same:
                // 1. directed canon value of the node 2. child key pos 3. aggre pos 4. aggre weight 5. node order isomorphic
                std::vector<CanonType> groupNodeCanon;
                std::vector<std::vector<std::vector<int>>> groupChildKeyPos;
                std::vector<std::vector<int>> groupAggrePos;
                std::vector<std::vector<int>> groupAggreWeight;
                std::vector<std::vector<int>> groupTIDs;
                std::vector<std::vector<VertexID>> groupNIDs;
                for (int i = 0; i < treeNodeID.size(); ++i) {
                    int tID = treeNodeID[i].first;
                    const Tree &t = allTree[tID];
                    const Pattern &p = allPattern[tID];
                    VertexID nID = treeNodeID[i].second;
                    int pID = nodePID[tID][nID];
                    int pIDCopy = pID;
                    for (int j = 0; j < pIDCopy; ++j) {
                        if (t.getGlobalOrder()[j].empty()) --pID;
                    }
                    const Node &tau = t.getNode(nID);
                    std::vector<std::vector<int>> childKeyPos = t.getChildKeyPos(nID);
                    std::vector<int> aggrePos = t.getAggrePos(nID);
                    std::vector<int> aggreWeight;
                    if (t.getRootID() == nID) {
                        aggreWeight = t.getAggreWeight();
                        // sort aggrePos and aggreWeight
                        const std::vector<int> aggrePosCopy = aggrePos;
                        const std::vector<int> aggreWeightCopy = aggreWeight;
                        sortAggrePosWeight(orbitType, aggrePos, aggreWeight, aggrePosCopy, aggreWeightCopy);
                    }
                    int pos = posInMergedNode[tID][nID];
                    setPatterns.emplace(tau.subPatternCanon, tau.keyOrbit, tau.prefixOrbit, tau.subPatternRules, tau.unEdgeKey, pID, pos);
                    bool groupExists = false;
                    for (int j = 0; j < groupNodeCanon.size(); j++) {
                        const Pattern &pRep = allPattern[groupTIDs[j][0]];
                        const Node &nRep = allTree[groupTIDs[j][0]].getNode(groupNIDs[j][0]);
                        if (groupNodeCanon[j] == tau.canonValue && groupChildKeyPos[j] == childKeyPos &&
                            groupAggrePos[j] == aggrePos && groupAggreWeight[j] == aggreWeight &&
                            checkMapping(p, pRep, tau.nodeOrder, nRep.nodeOrder)) {
                            groupTIDs[j].push_back(tID);
                            groupNIDs[j].push_back(nID);
                            groupExists = true;
                            break;
                        }
                    }
                    if (!groupExists) {
                        // build a new group if it is leaf node or can not assign to existing groups
                        groupNodeCanon.push_back(tau.canonValue);
                        groupChildKeyPos.push_back(childKeyPos);
                        groupAggrePos.push_back(aggrePos);
                        std::vector<int> tIDs = {tID};
                        std::vector<VertexID> nIDs = {nID};
                        groupTIDs.push_back(tIDs);
                        groupNIDs.push_back(nIDs);
                        groupAggreWeight.push_back(aggreWeight);
                    }
                }
                // build the node set for each child for each group representative.
                bool valid = true;
                std::vector<std::vector<int>> groupChildNum(groupNodeCanon.size());
                for (int group = 0; group < groupNodeCanon.size(); ++group) {
                    const std::vector<std::vector<int>> &childKeyPos = groupChildKeyPos[group];
                    const std::vector<int> &tIDs = groupTIDs[group];
                    const std::vector<VertexID> &nIDs = groupNIDs[group];
//                    if (childKeyPos.size() <= 1) continue;
                    int childNum = 1;
                    for (int j = 0; j < childKeyPos.size(); ++j) {
                        std::set<PatternInfo> childPattern;
                        for (int k = 0; k < tIDs.size(); ++k) {
                            int tID = tIDs[k];
                            const Tree &t = allTree[tID];
                            VertexID cID = t.getChild()[nIDs[k]][j];
                            const Node &child = t.getNode(cID);
                            int cPID = nodePID[tID][cID];
                            for (int l = 0; l < cPID; ++l) {
                                if (t.getGlobalOrder()[l].empty()) --cPID;
                            }
                            int cPos = posInMergedNode[tID][cID];
                            childPattern.emplace(child.subPatternCanon, child.keyOrbit, child.prefixOrbit,
                                                 child.subPatternRules, child.unEdgeKey, cPID, cPos);
                        }
                        groupChildNum[group].push_back(childPattern.size());
                        childNum *= (int)childPattern.size();
                    }
                    if (childNum != tIDs.size() / siblingNum) {
                        // this group is invalid. all rooted subtrees must be executed separately.
                        valid = false;
                        break;
                    }
                }
                bool vertexKey = allTree[setTRep].getNode(setNRep).keySize < 2;
                bool tableExists = false;
                if (!isRoot) {
                    for (int tableID = 0; tableID < childTables.size(); ++tableID) {
                        if (childPatterns[tableID] == setPatterns) {
                            for (auto & id : treeNodeID) {
                                int tID = id.first;
                                VertexID nID = id.second;
                                treeNodeH[tID][nID] = childTables[tableID];
                                childTID[tableID].push_back(tID);
                                childNID[tableID].push_back(nID);
                                posInChildTables[tID][nID] = tableID;
                            }
                            tableExists = true;
                            break;
                        }
                    }
                }
                if (!tableExists) {
                    if (!isRoot) {
                        HashTable h;
                        if (vertexKey) {
                            h = new Count[n];
                            memset(h, 0, sizeof(Count) * n);
                            allocatedHashTable.push_back(h);
                            ++nTable;
                        } else {
                            h = new Count[m];
                            memset(h, 0, sizeof(Count) * m);
                            allocatedHashTable.push_back(h);
                            mTable += 2;
                        }
                        for (auto &id: treeNodeID) {
                            treeNodeH[id.first][id.second] = h;
                        }
                        childTables.push_back(h);
                        if (valid) numTables.push_back(groupTIDs.size());
                        else numTables.push_back(treeNodeID.size());
                        childTID.emplace_back();
                        childNID.emplace_back();
                        childPatterns.push_back(setPatterns);
                        isPartition.push_back(isPart);
                        for (int i = 0; i < treeNodeID.size(); ++i) {
                            childTID.back().push_back(treeNodeID[i].first);
                            childNID.back().push_back(treeNodeID[i].second);
                            posInChildTables[treeNodeID[i].first][treeNodeID[i].second] = (int)childTables.size() - 1;
                        }
                    }
                    if (valid) {
                        for (int group = 0; group < groupNodeCanon.size(); ++group) {
                            const std::vector<std::vector<int>> &childKeyPos = groupChildKeyPos[group];
                            const std::vector<int> &tIDs = groupTIDs[group];
                            const std::vector<VertexID> &nIDs = groupNIDs[group];
                            nodeNeedExecute[tIDs[0]][nIDs[0]] = true;
                            for (int j = 0; j < childKeyPos.size(); ++j) {
                                std::vector<std::pair<int, VertexID>> next;
                                for (int k = 0; k < tIDs.size(); ++k) {
                                    const Tree &t = allTree[tIDs[k]];
                                    VertexID cID = t.getChild()[nIDs[k]][j];
                                    next.emplace_back(tIDs[k], cID);
                                }
                                q.push(next);
                                siblingQ.push(groupTIDs[group].size() / groupChildNum[group][j]);
                            }
                        }
                    }
                    else {
                        for (auto &idPair: treeNodeID) {
                            int tID = idPair.first;
                            VertexID nID = idPair.second;
                            nodeNeedExecute[tID][nID] = true;
                            const Tree &t = allTree[tID];
                            const std::vector<std::vector<VertexID>> &children = t.getChild();
                            for (int j = 0; j < children[nID].size(); ++j) {
                                VertexID cID = children[nID][j];
                                std::vector<std::pair<int, VertexID>> next;
                                next.emplace_back(tID, cID);
                                q.push(next);
                                siblingQ.push(1);
                            }
                        }
                    }
                }
            }
        }
    }
    numExecuted = std::vector<int>(childTables.size(), 0);
}

void sortAggrePosWeight(int orbitType, std::vector<int> &sortedAggrePos, std::vector<int> &sortedAggreWeight,
                        const std::vector<int> &oldAggrePos, const std::vector<int> &oldAggreWeight) {
    if (sortedAggrePos.empty()) {
        sortedAggreWeight = oldAggreWeight;
        return;
    }
    std::vector<int> pos(sortedAggrePos.size());
    if (orbitType == 2) {
        sortEOrbitAggrePos(sortedAggrePos);
        for (int i = 0; i < sortedAggrePos.size(); i = i + 2) {
            for (int j = 0; j < oldAggrePos.size(); j = j + 2) {
                if (sortedAggrePos[i] == oldAggrePos[j] && sortedAggrePos[i + 1] == oldAggrePos[j + 1]) {
                    pos[i] = j;
                    break;
                }
            }
        }
    }
    else {
        std::sort(sortedAggrePos.begin(), sortedAggrePos.end());
        for (int i = 0; i < sortedAggrePos.size(); ++i) {
            for (int j = 0; j < oldAggrePos.size(); ++j) {
                if (sortedAggrePos[i] == oldAggrePos[j]) {
                    pos[i] = j;
                    break;
                }
            }
        }
    }
    for (int i = 0; i < sortedAggreWeight.size(); ++i) {
        sortedAggreWeight[i] = oldAggreWeight[pos[i]];
    }
}

// check whether we can merge tau2 to tau1
bool mergeable(VertexID nID1, VertexID nID2, const Tree &t1, const Tree &t2, const Pattern &p1, const Pattern &p2,
               std::vector<int> &sortedAggrePos, std::vector<std::vector<int>> &sortedChildKeyPos,
               const std::vector<int> &aggreWeight, std::vector<int> &sortedAggreWeight, bool isRoot, int orbitType) {
    const Node &tau1 = t1.getNode(nID1);
    const Node &tau2 = t2.getNode(nID2);
    if (tau1.canonValue != tau2.canonValue) return false;
    if (tau1.edgeKey != tau2.edgeKey) return false;
    if (tau1.numRules != tau2.numRules) return false;
    const std::vector<int> &oldAggrePos = t2.getAggrePos(nID2);
    const std::vector<std::vector<int>> &oldChildKeyPos = t2.getChildKeyPos(nID2);
    const std::vector<VertexID> &nodeOrder1 = tau1.nodeOrder;
    const std::vector<VertexID> &nodeOrder2 = tau2.nodeOrder;
    int prefixSize1 = tau1.nodeOrder.size() - tau1.localOrder.size();
    int prefixSize2 = tau2.nodeOrder.size() - tau2.localOrder.size();
    if (prefixSize1 != prefixSize2) return false;
    const std::vector<std::vector<int>> &nodeInPos1 = t1.getNodeInPos(nID1);
    const std::vector<std::vector<int>> &nodeOutPos1 = t1.getNodeOutPos(nID1);
    const std::vector<std::vector<int>> &nodeUnPos1 = t1.getNodeUnPos(nID1);
    const std::vector<std::vector<int>> &nodeInPos2 = t2.getNodeInPos(nID2);
    const std::vector<std::vector<int>> &nodeOutPos2 = t2.getNodeOutPos(nID2);
    const std::vector<std::vector<int>> &nodeUnPos2 = t2.getNodeUnPos(nID2);
    const std::vector<std::vector<int>> &nodeGreaterPos1 = t1.getNodeGreaterPos(nID1);
    const std::vector<std::vector<int>> &nodeGreaterPos2 = t2.getNodeGreaterPos(nID2);
    if (!(nodeInPos1 == nodeInPos2 && nodeOutPos1 == nodeOutPos2 && nodeUnPos1 == nodeUnPos2 && nodeGreaterPos1 == nodeGreaterPos2)) return false;
    VertexID *mapping = new VertexID[MAX_PATTERN_SIZE];
    for (int i = 0; i < tau1.numVertices; ++i) {
        mapping[nodeOrder2[i]] = nodeOrder1[i];
    }
    if(!checkMapping(p1, p2, nodeOrder2, mapping)) {
        delete[] mapping;
        return false;
    }

    sortedAggrePos = oldAggrePos;
    sortedChildKeyPos = oldChildKeyPos;
    // build the new aggrePos and childKeyPos
    for (int i = 0; i < oldAggrePos.size(); ++i) {
        VertexID oldAggreV = nodeOrder2[oldAggrePos[i]];
        for (int j = 0; j < nodeOrder1.size(); ++j) {
            if (mapping[oldAggreV] == nodeOrder1[j]) {
                sortedAggrePos[i] = j;
                break;
            }
        }
    }
    for (int i = 0; i < oldChildKeyPos.size(); ++i) {
        for (int j = 0; j < oldChildKeyPos[i].size(); j++) {
            VertexID oldChildKey = nodeOrder2[oldChildKeyPos[i][j]];
            for (int k = 0; k < nodeOrder1.size(); ++k) {
                if (mapping[oldChildKey] == nodeOrder1[k]) {
                    sortedChildKeyPos[i][j] = k;
                    break;
                }
            }
        }
    }

    // sort aggregation vertices in ascending order and adjust aggreWeight accordingly
    if (isRoot) {
        const std::vector<int> aggrePosCopy = sortedAggrePos;
        sortedAggreWeight = aggreWeight;
        sortAggrePosWeight(orbitType, sortedAggrePos, sortedAggreWeight, aggrePosCopy, aggreWeight);
        std::vector<int> pos(sortedAggrePos.size());
        if (orbitType == 2) {
            sortEOrbitAggrePos(sortedAggrePos);
            for (int i = 0; i < sortedAggrePos.size(); i = i + 2) {
                for (int j = 0; j < aggrePosCopy.size(); j = j + 2) {
                    if (sortedAggrePos[i] == aggrePosCopy[j] && sortedAggrePos[i + 1] == aggrePosCopy[j + 1]) {
                        pos[i] = j;
                        break;
                    }
                }
            }
        }
        else {
            std::sort(sortedAggrePos.begin(), sortedAggrePos.end());
            for (int i = 0; i < sortedAggrePos.size(); ++i) {
                for (int j = 0; j < aggrePosCopy.size(); ++j) {
                    if (sortedAggrePos[i] == aggrePosCopy[j]) {
                        pos[i] = j;
                        break;
                    }
                }
            }
        }
        sortedAggreWeight = aggreWeight;
        if (!pos.empty()) {
            for (int i = 0; i < sortedAggreWeight.size(); ++i) {
                sortedAggreWeight[i] = aggreWeight[pos[i]];
            }
        }
    }

    delete[] mapping;
    return true;
}

void ExeParam::assignGroup(const std::vector<int> &newAggrePos, const std::vector<std::vector<int>> &newChildKeyPos,
                           bool isRt, ui keySize, ui n, ui m) {
    int shareID = (int)hashTables.size();
    bool aggreExists = false;
    for (int i = 0; i < aggrePoses.size(); ++i) {
        if (aggrePoses[i] == newAggrePos && (isRt == isRoot[i])) {
            aggreShareIDs[i].push_back(shareID);
            id2AggreKey.push_back(i);
            aggreExists = true;
        }
    }
    if (!aggreExists) {
        std::vector<int> id;
        id.push_back(shareID);
        id2AggreKey.push_back(aggrePoses.size());
        aggrePoses.push_back(newAggrePos);
        aggreShareIDs.push_back(id);
        isRoot.push_back(isRt);
    }
    bool childKeyExists = false;
    for (int i = 0; i < childKeyPoses.size(); ++i) {
        if (childKeyPoses[i] == newChildKeyPos) {
            id2ChildKey.push_back(i);
            childKeyExists = true;
        }
    }
    if (!childKeyExists) {
        id2ChildKey.push_back(childKeyPoses.size());
        childKeyPoses.push_back(newChildKeyPos);
    }
}

void ExeParam::buildPosesAndTypes(Forest &f) {
    for (int i = 0; i < nIDs.size(); ++i) {
        if (!trees[i].isLeaf(nIDs[i]))
            nonLeafIDs.push_back(i);
    }
    if (!rep.edgeKey) {
        isRoot = std::vector<bool>(nIDs.size());
        for (int shareID = 0 ; shareID < nIDs.size(); ++shareID) {
            if (trees[shareID].getRootID() == nIDs[shareID])
                isRoot[shareID] = true;
            else
                isRoot[shareID] = false;
        }
        return;
    }
    int orbitType = trees[0].getOrbitType();
    const std::vector<VertexID> &nodeOrder = rep.nodeOrder;
    // build posChildEdge and childEdgeType
    posChildEdge = std::vector<std::vector<std::pair<int, int>>>(nodeOrder.size());
    childEdgeType = std::vector<std::vector<int>>(childKeyPoses.size());
    for (int i = 0; i < childKeyPoses.size(); ++i) {
        childEdgeType[i] = std::vector<int>(childKeyPoses[i].size(), 0);
        for (int j = 0; j < childKeyPoses[i].size(); ++j) {
            if (childKeyPoses[i][j].size() == 2) {
                int pos1 = childKeyPoses[i][j][0], pos2 = childKeyPoses[i][j][1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                int smallerPos = pos1 < pos2 ? pos1 : pos2;
                posChildEdge[pos].emplace_back(i, j);
                if (!nodeInterPos[pos]) {
                    if (!nodeInPos[pos].empty()) childEdgeType[i][j] = 2;
                    else if (!nodeOutPos[pos].empty()) childEdgeType[i][j]= 1;
                    else {
                        if (pos == pos2) childEdgeType[i][j] = 1;
                        else childEdgeType[i][j] = 6;
                    }
                }
                else {
                    if (p.useDAG()) childEdgeType[i][j] = 3;
                    else {
                        if (std::find(nodeInPos[pos].begin(), nodeInPos[pos].end(), smallerPos) != nodeInPos[pos].end() ||
                            std::find(nodeOutPos[pos].begin(), nodeOutPos[pos].end(), smallerPos) != nodeOutPos[pos].end() ||
                            std::find(greaterPos[pos].begin(), greaterPos[pos].end(), smallerPos) != greaterPos[pos].end() ||
                            std::find(lessPos[pos].begin(), lessPos[pos].end(), smallerPos) != lessPos[pos].end())
                            childEdgeType[i][j] = 3;
                        else {
                            if (pos == pos2) childEdgeType[i][j] = 4;
                            else childEdgeType[i][j] = 7;
                        }
                    }
                }
            }
        }
    }
    // build posAggreEdge and aggreEdgeType
    posAggreEdge = std::vector<std::vector<std::pair<int, int>>>(nodeOrder.size());
    aggreEdgeType = std::vector<std::vector<int>>(aggrePoses.size());
    for (int i = 0; i < aggrePoses.size(); ++i) {
        if (!isRoot[i]) {
            if (aggrePoses[i].size() == 2) {
                int pos1 = aggrePoses[i][0], pos2 = aggrePoses[i][1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                int smallerPos = pos1 < pos2 ? pos1 : pos2;
                posAggreEdge[pos].emplace_back(i, 0);
                if (!nodeInterPos[pos]) {
                    if (!nodeInPos[pos].empty()) aggreEdgeType[i].push_back(2);
                    else if (!nodeOutPos[pos].empty()) aggreEdgeType[i].push_back(1);
                    else {
                        if (pos == pos2) aggreEdgeType[i].push_back(1);
                        else aggreEdgeType[i].push_back(6);
                    }
                }
                else {
                    if (p.useDAG()) aggreEdgeType[i].push_back(3);
                    else {
                        if (std::find(nodeInPos[pos].begin(), nodeInPos[pos].end(), smallerPos) != nodeInPos[pos].end() ||
                            std::find(nodeOutPos[pos].begin(), nodeOutPos[pos].end(), smallerPos) != nodeOutPos[pos].end() ||
                            std::find(greaterPos[pos].begin(), greaterPos[pos].end(), smallerPos) != greaterPos[pos].end() ||
                            std::find(lessPos[pos].begin(), lessPos[pos].end(), smallerPos) != lessPos[pos].end())
                            aggreEdgeType[i].push_back(3);
                        else {
                            if (pos == pos2) aggreEdgeType[i].push_back(4);
                            else aggreEdgeType[i].push_back(7);
                        }
                    }
                }
            }
        }
        else if (orbitType == 2) {
            aggreEdgeType[i] = std::vector<int>(aggrePoses[i].size() / 2);
            for (int j = 0; j < aggrePoses[j].size(); j = j + 2) {
                int pos1 = aggrePoses[i][j], pos2 = aggrePoses[i][j + 1];
                int pos = pos1 > pos2 ? pos1 : pos2;
                posAggreEdge[pos].emplace_back(i, j / 2);
                if (!nodeInterPos[pos]) {
                    if (!nodeOutPos[pos].empty()) aggreEdgeType[i][j / 2] = 1;
                    else if (!nodeInPos[pos].empty()) aggreEdgeType[i][j / 2] = 2;
                    else aggreEdgeType[i][j / 2] = 5;
                }
                else aggreEdgeType[i][j / 2] = 3;
            }
        }
    }
    isRoot = std::vector<bool>(nIDs.size());
    for (int shareID = 0 ; shareID < nIDs.size(); ++shareID) {
        if (trees[shareID].getRootID() == nIDs[shareID])
            isRoot[shareID] = true;
        else
            isRoot[shareID] = false;
    }
}