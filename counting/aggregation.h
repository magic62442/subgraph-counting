//
// Created by ??¨¬¡Â?? on 2022/10/7.
//
// This is to enforce sequential visit of hashtable. So we need to sort the key.

#ifndef SCOPE_AGGREGATION_H
#define SCOPE_AGGREGATION_H

#include "equation.h"

struct GroupKey {
    ui *aggreKey;
    ui *childKey;
};

struct KeyStruct {
    ui bound;
    ui size;
    GroupKey *groupKey;
    int aggreSize;
    int childSize;
    bool keyType;   // true: aggreKey, false: childKey
    int keyPos;
    bool exceedsBound;

    KeyStruct();

    KeyStruct(ui bound, int aggreSize, int childSize);

    KeyStruct(const KeyStruct &rhs);

    KeyStruct &operator=(const KeyStruct &rhs);

    virtual ~KeyStruct();

    void sort();
    bool addKey(ui *aggreKey, ui *childKey);
};

void aggregation(KeyStruct &keys, const std::vector<int> &shareIDs, bool isRoot, std::vector<HashTable *> &hashTables,
                 const std::vector<VertexID> &nIDs, const std::vector<std::vector<VertexID>> &children,
                 const std::vector<std::vector<int>> &aggreWeights, int orbitType);

#endif //SCOPE_AGGREGATION_H