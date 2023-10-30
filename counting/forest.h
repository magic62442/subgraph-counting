//
// Created by anonymous author on 2022/9/15.
//

#ifndef SCOPE_FOREST_H
#define SCOPE_FOREST_H

#include "equation.h"
#include "aggregation.h"
#include <stack>

struct TreeAndPattern {
    Tree t;
    int patternPos;

    TreeAndPattern(const Tree &t, int patternPos) : t(t), patternPos(patternPos) {}

    bool operator<(const TreeAndPattern &rhs) const {
        return t < rhs.t;
    }
};

struct NodeInfo {
    int qID;
    int tID;
    VertexID nID;
    std::vector<int> aggrePos;
    std::vector<int> aggreWeight;
    std::vector<std::vector<int>> childKeyPos;

    NodeInfo(int qId, int tId, VertexID nId, const std::vector<int> &aggrePos, const std::vector<int> &aggreWeight,
             const std::vector<std::vector<int>> &childKeyPos) : qID(qId), tID(tId), nID(nId), aggrePos(aggrePos),
                                                                 aggreWeight(aggreWeight), childKeyPos(childKeyPos) {}
};

struct PatternInfo {
    CanonType canon;
    int keyOrbit;
    std::vector<int> prefixOrbit;
    int subPatternRules;
    bool unEdgeKey;
    int pID;
    int pos;

    PatternInfo(CanonType canon, int keyOrbit, const std::vector<int> &prefixOrbit, int subPatternRules, bool unEdgeKey,
                int pId, int pos) : canon(canon), keyOrbit(keyOrbit), prefixOrbit(prefixOrbit),
                                    subPatternRules(subPatternRules), unEdgeKey(unEdgeKey), pID(pId), pos(pos) {}

    bool operator==(const PatternInfo &rhs) const {
        return std::tie(canon, keyOrbit, prefixOrbit, subPatternRules, unEdgeKey, pID, pos) ==
               std::tie(rhs.canon, rhs.keyOrbit, rhs.prefixOrbit, rhs.subPatternRules, rhs.unEdgeKey, rhs.pID, rhs.pos);
    }

    bool operator<(const PatternInfo &rhs) const {
        return std::tie(canon, keyOrbit, prefixOrbit, subPatternRules, unEdgeKey, pID, pos) <
               std::tie(rhs.canon, rhs.keyOrbit, rhs.prefixOrbit, rhs.subPatternRules, rhs.unEdgeKey, rhs.pID, rhs.pos);
    }
};

// after the equation is generated, we organize trees into a forest to allow the computation of nodes to be shared
struct Forest {
    int numQuery;
    int numMergedNodes;
    int orbitType;
    std::vector<Node> mergedNode;
    // we merge node i, j if: 1. canonValue is equal 2. prefix is empty or isomorphic
    std::vector<int> mergedNodeTID;
    std::vector<VertexID> mergedNodeNID;
    std::vector<std::vector<NodeInfo>> toTreeNode;
    std::vector<Tree> allTree;
    // for the nID-th node in the tID-th tree, store its pID
    std::vector<std::vector<int>> nodePID;
    // for the nID-th node in the tID-th tree, store the nID of its parent
    std::vector<std::vector<VertexID>> nodeParent;
    std::vector<Pattern> allPattern;
    // for the nID-th node in the tID-th tree, its position in mergedNode is posInMergedNode[i][nID]
    std::vector<std::vector<int>> posInMergedNode;
    // whether a node is finished. node executable means that prefixed hash table is computed.
    std::vector<std::vector<bool>> nodeExecutable;
    std::vector<std::vector<bool>> nodeExecuted;
    // for the nID-th node in the tID-th tree, whether it should be executed
    std::vector<std::vector<bool>> nodeNeedExecute;
    std::vector<std::vector<int>> queryFactor;
    std::vector<std::vector<int>> queryFactorFirstTree;
    std::vector<std::vector<int>> queryFactorNumTree;
    /******************* structures used for merging children *******************/
    std::vector<HashTable> childTables;
    std::vector<std::set<PatternInfo>> childPatterns;
    // whether this set of tables represents a complete partition or prefixed nodes
    std::vector<bool> isPartition;
    // for each set of tables, the number of tables needs to cover the set
    std::vector<int> numTables;
    std::vector<int> numExecuted;   // set to 0 initially
    // for each set of tables, store nodes to cover the set
    std::vector<std::vector<int>> childTID;
    std::vector<std::vector<VertexID>> childNID;
    // for the nID-th node in the tID-th tree, the position of its corresponding child table
    std::vector<std::vector<int>> posInChildTables;

    Forest() {
        numQuery = 0;
        numMergedNodes = 0;
        orbitType = 0;
    }

    void loadQuery(std::map<int, std::vector<Pattern>> &patterns, std::map<int, std::vector<std::vector<Tree>>> &trees,
                   bool share);
    const Pattern &mergedNodePattern(int pos) const;
    void updateNodeExecutable(int tableID);
    bool executable(const NodeInfo &inf) const;
    void allocateTables(HashTable **treeNodeH, std::vector<HashTable> &allocatedHashTable, ui n, ui m, int &nTable,
                        int &mTable);
};

// parameters needed for calling exeSharedNode
struct ExeParam {
    Node rep;
    /*********************** indexes of these containers are shareID ***********************/
    std::vector<std::vector<VertexID>> children;
    std::vector<VertexID> nIDs;
    std::vector<int> tIDs;
    std::vector<Tree> trees;
    std::vector<HashTable *> hashTables;
    std::vector<std::vector<int>> aggreWeights;
    std::vector<int> id2AggreKey;
    std::vector<int> id2ChildKey;
    std::vector<bool> isRoot;
    std::vector<int> tableIDs;
    /***************************************************************************************/
    /*********************** indexes of these containers are groupID ***********************/
    // group nodes by: 1. childKeyPose 2. aggrePos and isRoot
    std::vector<std::vector<std::vector<int>>> childKeyPoses;
    std::vector<std::vector<int>> aggrePoses;
    std::vector<std::vector<int>> aggreShareIDs;
    /***************************************************************************************/
    /*********************** indexes of these containers are mappingSize *******************/
    std::vector<bool> nodeInterPos;
    std::vector<std::vector<int>> nodeInPos;
    std::vector<std::vector<int>> nodeOutPos;
    std::vector<std::vector<int>> nodeUnPos;
    std::vector<bool> nodeCandPos;
    std::vector<std::vector<int>> greaterPos;
    std::vector<std::vector<int>> lessPos;
    // when the mapping size is i, store the position of child edge keys that are matched
    std::vector<std::vector<std::pair<int, int>>> posChildEdge;
    // when the mapping size is i, store the position of aggregation edge keys that are matched
    std::vector<std::vector<std::pair<int, int>>> posAggreEdge;
    std::vector<std::vector<int>> childEdgeType;
    std::vector<std::vector<int>> aggreEdgeType;
    // when the mappingSize is i, the positions of triangle neighbors
    std::vector<std::pair<int, int>> nodeTriPos;
    std::vector<int> triEdgeType;   // 1: edge ID at p1 2: outID at p1 3: compute
    std::vector<int> triEndType;    // 1: in-in 2: out-out 3: in-out
    /***************************************************************************************/
    std::vector<int> nonLeafIDs;
    VertexID **candidate;
    ui *candCount;
    Pattern p;

    void assignGroup(const std::vector<int> &newAggrePos, const std::vector<std::vector<int>> &newChildKeyPos,
                     bool isRt, ui keySize, ui n, ui m);
    void buildPosesAndTypes(Forest &f);
};

void sortAggrePosWeight(int orbitType, std::vector<int> &sortedAggrePos, std::vector<int> &sortedAggreWeight,
                        const std::vector<int> &oldAggrePos, const std::vector<int> &oldAggreWeight);
bool mergeable(VertexID nID1, VertexID nID2, const Tree &t1, const Tree &t2, const Pattern &p1, const Pattern &p2,
               std::vector<int> &sortedAggrePos, std::vector<std::vector<int>> &sortedChildKeyPos,
               const std::vector<int> &aggreWeight, std::vector<int> &sortedAggreWeight, bool isRoot, int orbitType);

#endif //SCOPE_FOREST_H