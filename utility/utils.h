//
// Created by anonymous author on 2022/11/21.
//

#ifndef SCOPE_UTILS_H
#define SCOPE_UTILS_H

#include "config.h"
#include <vector>
#include <algorithm>
#ifndef __APPLE__
#include <immintrin.h>
#include <x86intrin.h>
#endif

int quickPow10(int n);
Count choosec(Count n, int k);
void leftShit(std::vector<bool> &arr, int k);
std::vector<std::vector<bool>> chooseK(ui n, int k);
void sortEOrbitAggrePos(std::vector<int> &aggrePos);
ui firstPosGreaterThan(const ui *candidate, ui begin, ui end, ui target);
void generatePermutations(ui n, ui** automorphisms, ui* current, bool* used, size_t& index, size_t currentSize);
void generateCliqueRulesHelper(const std::vector<VertexID> &vertices,
                               std::vector<std::vector<std::vector<VertexID>>> &rules,
                               std::vector<std::vector<VertexID>> current);
std::vector<std::vector<std::vector<VertexID>>> generateCliqueRules(const std::vector<VertexID> &vertices,
                                                                    const std::vector<std::vector<std::vector<VertexID>>> &candidateRules);

#endif //SCOPE_UTILS_H
