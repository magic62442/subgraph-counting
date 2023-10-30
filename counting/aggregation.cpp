//
// Created by Qiyan LI on 2022/10/7.
//

#include "aggregation.h"

KeyStruct::KeyStruct() {
    bound = 0;
    size = 0;
    groupKey = nullptr;
    aggreSize = 0;
    childSize = 0;
    keyType = false;
    keyPos = 0;
    exceedsBound = false;
}

KeyStruct::KeyStruct(ui bound, int aggreSize, int childSize) :
        bound(bound), aggreSize(aggreSize), childSize(childSize) {
    groupKey = new GroupKey[bound + 1];
    for (int i = 0; i <= bound; ++i) {
        groupKey[i].aggreKey = new ui[aggreSize];
        groupKey[i].childKey = new ui[childSize];
    }
    size = 0;
    keyType = false;
    keyPos = 0;
    exceedsBound = false;
}

KeyStruct::~KeyStruct() {
    for (int i = 0; i < bound; ++i) {
        delete[] groupKey[i].aggreKey;
        delete[] groupKey[i].childKey;
    }
    delete[] groupKey;
}

void KeyStruct::sort() {
    if (keyType) {
        std::sort(groupKey, groupKey + size, [this](const GroupKey &lhs, const GroupKey &rhs) {
            return lhs.aggreKey[keyPos] < rhs.aggreKey[keyPos];
        });
    }
    else {
        std::sort(groupKey, groupKey + size, [this](const GroupKey &lhs, const GroupKey &rhs) {
            return lhs.childKey[keyPos] < rhs.childKey[keyPos];
        });
    }
}

bool KeyStruct::addKey(ui *aggreKey, ui *childKey) {
    memcpy(groupKey[size].aggreKey, aggreKey, sizeof(ui) * aggreSize);
    memcpy(groupKey[size].childKey, childKey, sizeof(ui) * childSize);
    ++size;
    return size < bound;
}

KeyStruct::KeyStruct(const KeyStruct &rhs) {
    bound = rhs.bound;
    size = rhs.size;
    aggreSize = rhs.aggreSize;
    childSize = rhs.childSize;
    keyType = rhs.keyType;
    keyPos = rhs.keyPos;
    exceedsBound = rhs.exceedsBound;
    groupKey = new GroupKey[bound + 1];
    for (int i = 0; i <= bound; ++i) {
        groupKey[i].aggreKey = new ui[aggreSize];
        groupKey[i].childKey = new ui[childSize];
        memcpy(groupKey[i].aggreKey, rhs.groupKey[i].aggreKey, sizeof(ui) * aggreSize);
        memcpy(groupKey[i].childKey, rhs.groupKey[i].childKey, sizeof(ui) * childSize);
    }
}

KeyStruct &KeyStruct::operator=(const KeyStruct &rhs) {
    if (this == &rhs) return *this;
    bound = rhs.bound;
    size = rhs.size;
    aggreSize = rhs.aggreSize;
    childSize = rhs.childSize;
    keyType = rhs.keyType;
    keyPos = rhs.keyPos;
    exceedsBound = rhs.exceedsBound;
    for (int i = 0; i < bound; ++i) {
        delete[] groupKey[i].aggreKey;
        delete[] groupKey[i].childKey;
    }
    delete[] groupKey;
    groupKey = new GroupKey[bound + 1];
    for (int i = 0; i <= bound; ++i) {
        groupKey[i].aggreKey = new ui[aggreSize];
        groupKey[i].childKey = new ui[childSize];
        memcpy(groupKey[i].aggreKey, rhs.groupKey[i].aggreKey, sizeof(ui) * aggreSize);
        memcpy(groupKey[i].childKey, rhs.groupKey[i].childKey, sizeof(ui) * childSize);
    }

    return *this;
}

void aggregation(KeyStruct &keys, const std::vector<int> &shareIDs, bool isRoot, std::vector<HashTable *> &hashTables,
                 const std::vector<VertexID> &nIDs, const std::vector<std::vector<VertexID>> &children,
                 const std::vector<std::vector<int>> &aggreWeights, int orbitType) {
    Count **allCount = new Count *[hashTables.size()];
    for (int shareID : shareIDs) {
        allCount[shareID] = new Count[keys.size];
        for (int i = 0; i < keys.size; ++i)
            allCount[shareID][i] = 1;
    }

    // multiply count in children
    keys.keyType = false;
    for (int i = 0; i < keys.childSize; ++i) {
        keys.keyPos = i;
        keys.sort();
        for (int shareID: shareIDs) {
            VertexID cID = children[shareID][i];
            HashTable h = hashTables[shareID][cID];
            Count *cnt = allCount[shareID];
            for (int j = 0; j < keys.size; ++j) {
                cnt[j] *= h[keys.groupKey[j].childKey[i]];
            }
        }
    }
    // update hash table
    keys.keyType = true;
    if (isRoot) {
        for (int i = 0; i < keys.aggreSize; ++i) {
            keys.keyPos = i;
            keys.sort();
            if (orbitType == 0) {
                for (int shareID: shareIDs) {
                    HashTable h = hashTables[shareID][nIDs[shareID]];
                    Count *cnt = allCount[shareID];
                    for (int j = 0; j < keys.size; ++j) {
                        h[0] += cnt[j];
                    }
                }
            }
            else {
                for (int shareID: shareIDs) {
                    HashTable h = hashTables[shareID][nIDs[shareID]];
                    Count *cnt = allCount[shareID];
                    int weight = aggreWeights[shareID][i];
                    if (weight == 1) {
                        for (int j = 0; j < keys.size; ++j)
                            h[keys.groupKey[j].aggreKey[i]] += cnt[j];
                    }
                    else if (weight == -1) {
                        for (int j = 0; j < keys.size; ++j)
                            h[keys.groupKey[j].aggreKey[i]] -= cnt[j];
                    }
                    else {
                        for (int j = 0; j < keys.size; ++j)
                            h[keys.groupKey[j].aggreKey[i]] += cnt[j] * weight;
                    }
                }
            }
        }
    }
    else {
        keys.keyPos = 0;
        keys.sort();
        for (int shareID: shareIDs) {
            HashTable h = hashTables[shareID][nIDs[shareID]];
            Count *cnt = allCount[shareID];
            for (int j = 0; j < keys.size; ++j) {
                h[keys.groupKey[j].aggreKey[0]] += cnt[j];
            }
        }
    }
    // check whether size exceeds bound and reset size
    if (keys.size >= keys.bound) {
        keys.size = 0;
        keys.exceedsBound = true;
    }
    for (int shareID : shareIDs) {
        delete[] allCount[shareID];
    }
    delete[] allCount;
}