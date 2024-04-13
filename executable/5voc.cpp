//
// the ad-hoc version for 5-vertex orbit counting. some functions are adapted from EVOKE
//

#include "execution.h"
#include "command.h"

void fourCycle(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
               std::vector<HashTable> &eTables, EdgeID *outID) {
    ui n = dun.getNumVertices();
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    HashTable wedgeCount = new Count[n];
    memset(wedgeCount, 0, sizeof(Count) * n);
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
                ++wedgeCount[k];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
            }
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
                ++wedgeCount[k];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            EdgeID ji = outID[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
                // i <- j -> k
                vTables[8][j] += std::max((Count)0, wedgeCount[k] - 1);
                eTables[5][ji] += std::max((Count)0, wedgeCount[k] - 1);
                eTables[5][jk] += std::max((Count)0, wedgeCount[k] - 1);
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
            }
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                EdgeID kj = outID[jk];
                // i <- j <- k
                vTables[8][j] += std::max((Count)0, wedgeCount[k] - 1);
                eTables[5][ji] += std::max((Count)0, wedgeCount[k] - 1);
                eTables[5][kj] += std::max((Count)0, wedgeCount[k] - 1);
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // i <- j -> k
                vTables[8][i] += choosec(wedgeCount[k], 2);
                vTables[8][k] += choosec(wedgeCount[k], 2);
                wedgeCount[k] = 0;
            }
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // i <- j <- k
                vTables[8][i] += choosec(wedgeCount[k], 2);
                vTables[8][k] += choosec(wedgeCount[k], 2);
                wedgeCount[k] = 0;
            }
        }
    }

    delete[] wedgeCount;
}

void otherFourOrbit(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
                    std::vector<HashTable> &eTables, std::vector<HashTable> &revETables) {
    ui n = dun.getNumVertices();
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    for (VertexID i = 0; i < n; ++i) {
        vTables[5][i] = vTables[1][i] * (vTables[0][i] - 1) - 2 * vTables[3][i];
        vTables[7][i] = choosec(vTables[0][i], 3);
        vTables[11][i] = vTables[3][i] * (vTables[0][i] - 2);
        for (EdgeID ij = outOffset[i]; ij < outOffset[i + 1]; ++ij) {
            VertexID j = outNbors[ij];
#ifdef COLLECT_STATISTICS
            ++gNumMatch;
#endif
            vTables[9][i] += vTables[3][j] - eTables[1][ij];
            vTables[10][i] += eTables[1][ij] * (vTables[0][j] - 2);
            vTables[13][j] += choosec(eTables[1][ij], 2);
            ui degJ = vTables[0][j];
            vTables[4][i] += vTables[1][j];
            vTables[6][i] += choosec(degJ - 1, 2);
            for (EdgeID ik = ij + 1; ik < outOffset[i + 1]; ++ik) {
                VertexID k = outNbors[ik];
                int jk = dout.getEdgeID(j, k);
                if (jk != -1) {
#ifdef COLLECT_STATISTICS
                    ++gNumMatch;
#endif
                    vTables[12][i] += eTables[1][jk] - 1;
                    vTables[12][j] += eTables[1][ik] - 1;
                    vTables[12][k] += eTables[1][ij] - 1;
                    eTables[7][ij] += vTables[0][k]-2;
                    eTables[7][ik] += vTables[0][j]-2;
                    eTables[7][jk] += vTables[0][i]-2;
                    eTables[9][ij] += eTables[1][ik] - 1;
                    revETables[9][ij] += eTables[1][jk] - 1;
                    eTables[9][ik] += eTables[1][ij] - 1;
                    revETables[9][ik] += eTables[1][jk] - 1;
                    eTables[9][jk] += eTables[1][ij] - 1;
                    revETables[9][jk] += eTables[1][ik] - 1;
                }
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            EdgeID ji = dout.getEdgeID(i, j);
            vTables[9][i] += vTables[3][j] - eTables[1][ji];
            vTables[10][i] += eTables[1][ji] * (vTables[0][j] - 2);
            vTables[13][j] += (eTables[1][ji] * (eTables[1][ji] - 1)) / 2;
            vTables[4][i] += vTables[1][j];
            ui degJ = vTables[0][j];
            vTables[6][i] += choosec(degJ - 1, 2);
        }
        vTables[4][i] -= 2 * vTables[2][i] + 2 * vTables[3][i];
    }
}

void fourCliqueOrbit(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
                     std::vector<HashTable> &eTables) {
    ui n = dun.getNumVertices();
    VertexID *triends = new VertexID[n]; // array to store triangle ends
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    EdgeID *edge_pos_i_triend = new EdgeID[n]; // array to store triangle ends
    EdgeID *edge_pos_j_triend = new EdgeID[n]; // array to store triangle ends
    EdgeID edge_pos_ij, edge_pos_ik, edge_pos_il, edge_pos_jk, edge_pos_jl, edge_pos_kl;

    for (VertexID i = 0; i < n; ++i) // loop over vertices
    {
        for (EdgeID posj = outOffset[i]; posj < outOffset[i+1]; ++posj) // loop over out-neighbors of i
        {
            VertexID j = outNbors[posj]; // j is current out-neighbor
            edge_pos_ij = posj;

            VertexID count = 0;
            for (EdgeID posk = posj+1; posk < outOffset[i+1]; ++posk) // loop over another out-neighbor of i, that is "ahead" of j in list of out-neighbors
            {
                VertexID k = outNbors[posk]; // k is next out-neighbor
                edge_pos_jk = dout.getEdgeID(j, k);
                if (edge_pos_jk != -1) // check if edge (j,k) is present
                {
                    triends[count] = k;  // so (i,j,k) form a triangle. we store the fact that k forms a triangle with edge (i,j) in digraph gout
                    edge_pos_i_triend[count] = posk;
                    edge_pos_j_triend[count] = edge_pos_jk;

                    ++count;
                }
            }
            for (EdgeID posk = 0; posk < count; ++posk) // loop over all pairs of triangles formed by (i,j)
            {
                VertexID k = triends[posk]; // k is vertex as index posk in triends
                edge_pos_ik = edge_pos_i_triend[posk];
                edge_pos_jk = edge_pos_j_triend[posk];
                for (VertexID posell = posk+1; posell < count; ++posell)
                {
                    VertexID ell = triends[posell];

                    edge_pos_il = edge_pos_i_triend[posell];
                    edge_pos_jl = edge_pos_j_triend[posell];
                    edge_pos_kl = dout.getEdgeID(k, ell);
                    if (edge_pos_kl != -1) // (k,ell) is an end, thus (i,j,k,ell) form a 4-clique
                    {
#ifdef COLLECT_STATISTICS
                        ++gNumMatch;
#endif
                        vTables[14][i]++;
                        vTables[14][j]++;
                        vTables[14][k]++;
                        vTables[14][ell]++;
                        eTables[11][edge_pos_ij]++;
                        eTables[11][edge_pos_ik]++;
                        eTables[11][edge_pos_il]++;
                        eTables[11][edge_pos_jk]++;
                        eTables[11][edge_pos_jl]++;
                        eTables[11][edge_pos_kl]++;
                    }
                }
            }
        }
    }
    delete[] triends;
    delete[] edge_pos_i_triend;
    delete[] edge_pos_j_triend;
}

void fiveTree(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables) {
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    for (VertexID i=0; i < dun.getNumVertices(); i++) // loop over vertices
    {
#ifdef COLLECT_STATISTICS
        ++gNumMatch;
#endif
        vTables[16][i] = vTables[4][i] * (vTables[0][i]-1) - vTables[10][i] - 2 * vTables[8][i];
        vTables[17][i] = choosec(vTables[1][i], 2) - vTables[3][i] - vTables[6][i] - vTables[8][i] - vTables[10][i];
        vTables[23][i] = choosec(vTables[0][i], 4);
        for (EdgeID posj = outOffset[i]; posj < outOffset[i+1]; posj++) // loop over out-neighbors of i
        {
#ifdef COLLECT_STATISTICS
            ++gNumMatch;
#endif
            VertexID j = outNbors[posj];
            vTables[15][i] += vTables[4][j];
            vTables[18][i] += vTables[6][j];
            vTables[19][i] += vTables[5][j];
            vTables[20][i] += (vTables[0][i] - 1) * choosec(vTables[0][j] - 1, 2);
            vTables[21][i] += (vTables[0][j] - 1) * choosec(vTables[0][i] - 1, 2);
            vTables[22][i] += choosec(vTables[0][j] - 1, 3);
        }
        for (EdgeID posj = inOffset[i]; posj < inOffset[i+1]; posj++)   // loop over in-neighbor of i
        {
#ifdef COLLECT_STATISTICS
            ++gNumMatch;
#endif
            VertexID j = inNbors[posj];
            vTables[15][i] += vTables[4][j];
            vTables[18][i] += vTables[6][j];
            vTables[19][i] += vTables[5][j];
            vTables[20][i] += (vTables[0][i] - 1) * choosec(vTables[0][j] - 1, 2);
            vTables[21][i] += (vTables[0][j] - 1) * choosec(vTables[0][i] - 1, 2);
            vTables[22][i] += choosec(vTables[0][j] - 1, 3);
        }
        vTables[15][i] -= vTables[5][i] + 2 * vTables[8][i] + 2 * vTables[11][i];
        vTables[18][i] -= 3 * vTables[7][i] + vTables[10][i];
        vTables[19][i] -= vTables[4][i] + vTables[5][i] + vTables[10][i];
        vTables[20][i] -= vTables[10][i];
        vTables[21][i] -= 2*vTables[11][i];
    }
}

void otherFiveOrbits(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
                     std::vector<HashTable> &eTables, std::vector<HashTable> &revETables, EdgeID *outID) {
    ui n = dun.getNumVertices();
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    for (VertexID i = 0; i < n; ++i) {
#ifdef COLLECT_STATISTICS
        ++gNumMatch;
#endif
        vTables[33][i] = vTables[3][i] * choosec(vTables[0][i]-2, 2);
        vTables[38][i] = vTables[8][i] * (vTables[0][i] - 2) - vTables[13][i];
        vTables[44][i] = choosec(vTables[3][i], 2) - vTables[13][i];
        vTables[47][i] = vTables[12][i] * (vTables[0][i] - 2) - 3 * vTables[14][i];
        vTables[58][i] = vTables[14][i] * (vTables[0][i] - 3);
        for (EdgeID posj = outOffset[i]; posj < outOffset[i+1]; posj++) // loop over out-neighbors of i
        {
#ifdef COLLECT_STATISTICS
            ++gNumMatch;
#endif
            VertexID j = outNbors[posj];
            vTables[24][i] += vTables[10][j];
            vTables[27][i] += vTables[9][j];
            vTables[28][i] += (vTables[0][i] - 1) * (vTables[3][j]);
            vTables[30][i] += (vTables[0][j] - 1) * (vTables[3][i]);
            vTables[31][i] += (vTables[0][j] - 1) * vTables[3][j];
            vTables[35][i] += vTables[8][j];
            vTables[37][i] += eTables[5][posj] * (vTables[0][j] - 2);
            vTables[37][j] += eTables[5][posj] * (vTables[0][i] - 2);
            vTables[39][i] += vTables[13][j];
            vTables[40][i] += revETables[9][posj] * (vTables[0][j] - 3);
            vTables[40][j] += eTables[9][posj] * (vTables[0][i] - 3);
            vTables[41][i] += choosec(eTables[1][posj], 2) * (vTables[0][j] - 3);
            vTables[41][j] += choosec(eTables[1][posj], 2) * (vTables[0][i] - 3);
            vTables[42][i] += choosec(eTables[1][posj], 2) * (vTables[0][i] - 3);
            vTables[42][j] += choosec(eTables[1][posj], 2) * (vTables[0][j] - 3);
            vTables[45][i] += vTables[12][j];
            vTables[55][i] += choosec(eTables[1][posj], 3);
            vTables[55][j] += choosec(eTables[1][posj], 3);
            vTables[56][i] += vTables[14][j];
            vTables[57][i] += eTables[11][posj] * (vTables[0][j] - 3);
            vTables[57][j] += eTables[11][posj] * (vTables[0][i] - 3);
            vTables[60][i] += revETables[9][posj] * (eTables[1][posj] - 1);
            vTables[60][j] += eTables[9][posj] * (eTables[1][posj]  - 1);
            vTables[67][i] += eTables[11][posj] * (eTables[1][posj] - 2);
            vTables[67][j] += eTables[11][posj] * (eTables[1][posj] - 2);
            for (EdgeID posk = posj+1; posk < outOffset[i+1]; posk++) // loop over another out-neighbor of i, that is "ahead" of j in list of out-neighbors
            {
                VertexID k = outNbors[posk];
                EdgeID jk = dout.getEdgeID(j,k);
                if (jk != -1)
                {
#ifdef COLLECT_STATISTICS
                    ++gNumMatch;
#endif
                    vTables[25][i] += (vTables[0][j] - 2) * (vTables[0][k] - 2);
                    vTables[25][j] += (vTables[0][i] - 2) * (vTables[0][k] - 2);
                    vTables[25][k] += (vTables[0][i] - 2) * (vTables[0][j] - 2);
                    vTables[26][i] += (vTables[0][i] - 2) * ((vTables[0][j] - 2) + (vTables[0][k] - 2));
                    vTables[26][j] += (vTables[0][j] - 2) * ((vTables[0][i] - 2) + (vTables[0][k] - 2));
                    vTables[26][k] += (vTables[0][k] - 2) * ((vTables[0][i] - 2) + (vTables[0][j] - 2));
                    vTables[29][i] += vTables[1][j] + vTables[1][k];
                    vTables[29][j] += vTables[1][i] + vTables[1][k];
                    vTables[29][k] += vTables[1][i] + vTables[1][j];
                    vTables[32][i] += choosec(vTables[0][j]-2, 2) + choosec(vTables[0][k]-2, 2);
                    vTables[32][j] += choosec(vTables[0][i]-2, 2) + choosec(vTables[0][k]-2, 2);
                    vTables[32][k] += choosec(vTables[0][i]-2, 2) + choosec(vTables[0][j]-2, 2);
                    vTables[43][i] += vTables[3][j] - 1 + vTables[3][k] - 1;
                    vTables[43][j] += vTables[3][i] - 1 + vTables[3][k] - 1;
                    vTables[43][k] += vTables[3][i] - 1 + vTables[3][j] - 1;
                    vTables[46][i] += eTables[7][jk];
                    vTables[46][j] += eTables[7][posk];
                    vTables[46][k] += eTables[7][posj];
                    vTables[48][i] += (eTables[1][posj] - 1) * (vTables[0][k] - 2) + (eTables[1][posk] - 1) * (vTables[0][j] - 2);
                    vTables[48][j] += (eTables[1][posj] - 1) * (vTables[0][k] - 2) + (eTables[1][jk] - 1) * (vTables[0][i] - 2);
                    vTables[48][k] += (eTables[1][posk] - 1) * (vTables[0][j] - 2) + (eTables[1][jk] - 1) * (vTables[0][i] - 2);
                    vTables[52][i] += eTables[5][jk];
                    vTables[52][j] += eTables[5][posk];
                    vTables[52][k] += eTables[5][posj];
                    vTables[53][i] += eTables[5][posj] + eTables[5][posk];
                    vTables[53][j] += eTables[5][posj] + eTables[5][jk];
                    vTables[53][k] += eTables[5][posk] + eTables[5][jk];
                    vTables[54][i] += choosec(eTables[1][jk] - 1, 2);
                    vTables[54][j] += choosec(eTables[1][posk] - 1, 2);
                    vTables[54][k] += choosec(eTables[1][posj] - 1, 2);
                    vTables[59][i] += eTables[9][jk] + revETables[9][jk];
                    vTables[59][j] += eTables[9][posk] + revETables[9][posk];
                    vTables[59][k] += eTables[9][posj] + revETables[9][posj];
                    vTables[61][i] += (eTables[1][posj] - 1) * (eTables[1][posk] - 1);
                    vTables[61][j] += (eTables[1][posj] - 1) * (eTables[1][jk] - 1);
                    vTables[61][k] += (eTables[1][posk] - 1) * (eTables[1][jk] - 1);
                    vTables[65][i] += eTables[11][jk];
                    vTables[65][j] += eTables[11][posk];
                    vTables[65][k] += eTables[11][posj];
                }
            }
        }
        for (EdgeID posj = inOffset[i]; posj < inOffset[i+1]; posj++)   // loop over in-neighbor of i
        {
#ifdef COLLECT_STATISTICS
            ++gNumMatch;
#endif
            VertexID j = inNbors[posj];
            vTables[24][i] += vTables[10][j];
            vTables[27][i] += vTables[9][j];
            vTables[28][i] += (vTables[0][i] - 1) * (vTables[3][j]);
            vTables[30][i] += (vTables[0][j] - 1) * (vTables[3][i]);
            vTables[31][i] += (vTables[0][j] - 1) * vTables[3][j];
            vTables[35][i] += vTables[8][j];
            vTables[39][i] += vTables[13][j];
            vTables[45][i] += vTables[12][j];
            vTables[56][i] += vTables[14][j];
        }
        vTables[24][i] -= vTables[10][i] + 2 * vTables[11][i] + 2 * vTables[12][i];
        vTables[25][i] -= vTables[12][i];
        vTables[27][i] -= vTables[11][i] + 2 * vTables[13][i];
        vTables[26][i] -= 2 * vTables[13][i];
        vTables[28][i] -= 2 * vTables[3][i] + 2 * vTables[11][i] + 2 * vTables[12][i];
        vTables[29][i] -= 4 * vTables[3][i] + vTables[10][i] + 2 * vTables[11][i] + 2 * vTables[12][i] + 2 * vTables[13][i];
        vTables[30][i] -= 2 * vTables[3][i] + vTables[10][i] + 2 * vTables[13][i];
        vTables[31][i] -= 2 * vTables[3][i] + 2 * vTables[9][i] + vTables[10][i];
        vTables[35][i] -= 2 * vTables[8][i] + vTables[13][i];
        vTables[37][i] -= 2 * vTables[12][i];
        vTables[39][i] -= 2 * vTables[12][i] + vTables[13][i];
        vTables[43][i] -= 2 * vTables[12][i] + 2 * vTables[13][i];
        vTables[46][i] -= vTables[11][i] + 3 * vTables[14][i];
        vTables[48][i] -= 6 * vTables[14][i];
        vTables[45][i] -= 2 * vTables[13][i] + 3 * vTables[14][i];
        vTables[52][i] -= 2 * vTables[13][i];
        vTables[53][i] -= 2 * vTables[13][i] + 2 * vTables[12][i];
        vTables[56][i] -= 3 * vTables[14][i];
        vTables[59][i] -= 2 * vTables[13][i] + 6 * vTables[14][i];
        vTables[60][i] -= 6 * vTables[14][i];
        vTables[61][i] -= 3 * vTables[14][i];
        vTables[65][i] -= 3 * vTables[14][i];
    }
}

void fiveCycle(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
               const Triangle &tri, EdgeID *outID) {
    ui n = dun.getNumVertices();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    HashTable wedgeCount = new Count[n];
    HashTable pathCount = new Count[n];
    memset(wedgeCount, 0, sizeof(Count) * n);
    memset(pathCount, 0, sizeof(Count) * n);
    ui *keyPos = new ui[n / 8 + 1];
    ui keyPosSize = 0;
    // DP1, DP2, DP3
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                ++wedgeCount[k];
            }
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k != i) {
#ifdef COLLECT_STATISTICS
                    ++gNumMatch;
#endif
                    ++wedgeCount[k];
                }
            }
        }
        for (EdgeID ij = outOffset[i]; ij < outOffset[i + 1]; ++ij) {
            VertexID j = outNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                ++wedgeCount[k];
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k == i) continue;
                for (EdgeID kl = outOffset[k]; kl < outOffset[k + 1]; ++kl) {
                    VertexID l = outNbors[kl];
                    if (l == i) continue;
#ifdef COLLECT_STATISTICS
                    ++gNumMatch;
#endif
                    ++pathCount[l];
                    if (keyPosSize < n / 8) {
                        keyPos[keyPosSize] = l;
                        ++keyPosSize;
                    }
                    vTables[34][j] += wedgeCount[l];
                    vTables[34][k] += wedgeCount[l];
                    vTables[34][l] += wedgeCount[l];
                    vTables[34][i] += wedgeCount[l];
                }
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                vTables[34][j] += pathCount[k];
                wedgeCount[k] = 0;
            }
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (i == k) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                vTables[34][j] += pathCount[k];
                wedgeCount[k] = 0;
            }
        }
        for (EdgeID ij = outOffset[i]; ij < outOffset[i + 1]; ++ij) {
            VertexID j = outNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                vTables[34][j] += pathCount[k];
                wedgeCount[k] = 0;
            }
        }
        if (keyPosSize == n / 8) {
            memset(pathCount, 0, sizeof(Count) * n);
        }
        else {
            for (int pos = 0; pos < keyPosSize; ++pos)
                pathCount[keyPos[pos]] = 0;
        }
    }
    Count *vertexTriCnt = new Count[n];
    memset(vertexTriCnt, 0, sizeof(Count) * n);

    // shr4
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = outOffset[i]; ij < outOffset[i + 1]; ++ij) {
            VertexID j = outNbors[ij];
            VertexID *candidateK;
            ui cntK;
            candidateK = tri.getTriend(ij, 1, cntK);
            for (int pos = 0; pos < cntK; ++pos) {
                VertexID k = candidateK[pos];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                vTables[34][k] += 5;
                vTables[34][i] += 5;
                vTables[34][j] += 5;
                ++vertexTriCnt[k];
                ++vertexTriCnt[i];
                ++vertexTriCnt[j];
            }
        }
    }
    // shr1, shr2, shr3
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = outOffset[i]; ij < outOffset[i + 1]; ++ij) {
            VertexID j = outNbors[ij];
#ifdef COLLECT_STATISTICS
            ++gNumMatch;
#endif
            Count cnt = vertexTriCnt[i];
            vTables[34][i] -= 2 * cnt;
            vTables[34][j] -= cnt;
        }
    }
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = outOffset[i]; ij < outOffset[i + 1]; ++ij) {
            VertexID j = outNbors[ij];
            VertexID *candidateK;
            ui cntK;
            candidateK = tri.getTriend(ij, 1, cntK);
            Count cnt = outOffset[i + 1] - outOffset[i];
            for (int pos = 0; pos < cntK; ++pos) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                VertexID k = candidateK[pos];
                vTables[34][k] -= cnt;
                vTables[34][j] -= cnt;
            }
            candidateK = tri.getTriend(ij, 3, cntK);
            for (int pos = 0; pos < cntK; ++pos) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                VertexID k = candidateK[pos];
                vTables[34][j] -= cnt;
                vTables[34][k] -= cnt;
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            EdgeID ji = outID[ij];
            VertexID *candidateK;
            ui cntK;
            candidateK = tri.getTriend(ji, 2, cntK);
            Count cnt = outOffset[i + 1] - outOffset[i];
            for (int pos = 0; pos < cntK; ++pos) {
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                VertexID k = candidateK[pos];
                vTables[34][j] -= cnt;
                vTables[34][k] -= cnt;
            }
        }
    }
    delete[] wedgeCount;
    delete[] pathCount;
    delete[] keyPos;
    delete[] vertexTriCnt;
}

void wedgeCollision(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables) {
    ui n = dun.getNumVertices();
    const EdgeID *unOffset = dun.getOffsets();
    const VertexID *unNbors = dun.getNbors();
    Count *wedgeCount = new Count[n];
    Count *wedgeCountCopy = new Count[n];
    memset(wedgeCount, 0, sizeof(Count) * n);
    memset(wedgeCountCopy, 0, sizeof(Count) * n);
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = unOffset[i]; ij < unOffset[i + 1]; ++ij) {
            VertexID j = unNbors[ij];
            for (EdgeID jk = unOffset[j]; jk < unOffset[j + 1]; ++jk) {
                VertexID k = unNbors[jk];
                if (k <= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                ++wedgeCount[k];
                wedgeCountCopy[k] = wedgeCount[k];
            }
        }
        for (EdgeID ij = unOffset[i]; ij < unOffset[i + 1]; ++ij) {
            VertexID j = unNbors[ij];
            for (EdgeID jk = unOffset[j]; jk < unOffset[j + 1]; ++jk) {
                VertexID k = unNbors[jk];
                if (k <= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                vTables[49][j] += choosec(wedgeCountCopy[k]-1, 2);
                vTables[50][i] += choosec(wedgeCount[k], 3);
                vTables[50][k] += choosec(wedgeCount[k], 3);
                wedgeCount[k] = 0;
            }
        }
    }
    delete[] wedgeCount;
    delete[] wedgeCountCopy;
}

void orbit51(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
             std::vector<HashTable> &eTables, EdgeID *outID) {
    ui n = dun.getNumVertices();
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    HashTable ijTriangle = new Count[n];
    HashTable jkTriangle = new Count[n];
    memset(ijTriangle, 0, sizeof(Count) * n);
    memset(jkTriangle, 0, sizeof(Count) * n);
    ui *keyPos = new ui[n / 8 + 1];
    ui keyPosSize = 0;
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            EdgeID ji = outID[ij];
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                EdgeID kj = outID[jk];
                ijTriangle[k] += eTables[1][ji];
                jkTriangle[k] += eTables[1][kj];
                if (keyPosSize < n / 8) {
                    keyPos[keyPosSize] = k;
                    ++keyPosSize;
                }
            }
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                ijTriangle[k] += eTables[1][ji];
                jkTriangle[k] += eTables[1][jk];
                if (keyPosSize < n / 8) {
                    keyPos[keyPosSize] = k;
                    ++keyPosSize;
                }
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            EdgeID ji = outID[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                vTables[51][i] += jkTriangle[k] - eTables[1][jk];
                vTables[51][j] += jkTriangle[k] - eTables[1][jk];
                vTables[51][j] += ijTriangle[k] - eTables[1][ji];
                vTables[51][k] += ijTriangle[k] - eTables[1][ji];
            }
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                EdgeID kj = outID[jk];
                vTables[51][i] += jkTriangle[k] - eTables[1][kj];
                vTables[51][j] += jkTriangle[k] - eTables[1][kj];
                vTables[51][j] += ijTriangle[k] - eTables[1][ji];
                vTables[51][k] += ijTriangle[k] - eTables[1][ji];
            }
        }
        if (keyPosSize == n / 8) {
            memset(ijTriangle, 0, sizeof(Count) * n);
            memset(jkTriangle, 0, sizeof(Count) * n);
        }
        else {
            for (int pos = 0; pos < keyPosSize; ++pos) {
                ijTriangle[keyPos[pos]] = 0;
                jkTriangle[keyPos[pos]] = 0;
            }
        }
        keyPosSize = 0;
    }
    for (VertexID i = 0; i < n; ++i)
        vTables[51][i] -= 2 * vTables[12][i] + 2 * vTables[13][i];
    delete[] ijTriangle;
    delete[] jkTriangle;
    delete[] keyPos;
}

void orbit36(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables) {
    ui n = dun.getNumVertices();
    const EdgeID *inOffset = din.getOffsets();
    const VertexID *inNbors = din.getNbors();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    HashTable nborEdgeOutCount = new Count[n];
    HashTable wedgeCount = new Count[n];
    memset(nborEdgeOutCount, 0, sizeof(Count) * n);
    memset(wedgeCount, 0, sizeof(Count) * n);
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                ++wedgeCount[k];
                nborEdgeOutCount[k] += std::max((Count)0, vTables[0][j] - 2);
            }
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                ++wedgeCount[k];
                nborEdgeOutCount[k] += std::max((Count)0, vTables[0][j] - 2);
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // i <- j -> k
                vTables[36][j] += std::max(Count(0), nborEdgeOutCount[k] - vTables[0][j] + 2);
            }
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // i <- j <- k
                vTables[36][j] += std::max(Count(0), nborEdgeOutCount[k] - vTables[0][j] + 2);
            }
        }
        for (EdgeID ij = inOffset[i]; ij < inOffset[i + 1]; ++ij) {
            VertexID j = inNbors[ij];
            for (EdgeID jk = outOffset[j]; jk < outOffset[j + 1]; ++jk) {
                VertexID k = outNbors[jk];
                if (k >= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // i <- j -> k
                if (wedgeCount[k] > 1) {
                    vTables[36][i] += choosec(wedgeCount[k], 2) * std::max(Count(0), vTables[0][k] - 2);
                    vTables[36][k] += choosec(wedgeCount[k], 2) * std::max(Count(0), vTables[0][i] - 2);
                }
                nborEdgeOutCount[k] = 0;
                wedgeCount[k] = 0;
            }
            for (EdgeID jk = inOffset[j]; jk < inOffset[j + 1]; ++jk) {
                VertexID k = inNbors[jk];
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                // i <- j <- k
                if (wedgeCount[k] > 1) {
                    vTables[36][i] += choosec(wedgeCount[k], 2) * std::max(Count(0), vTables[0][k] - 2);
                    vTables[36][k] += choosec(wedgeCount[k], 2) * std::max(Count(0), vTables[0][i] - 2);
                }
                nborEdgeOutCount[k] = 0;
                wedgeCount[k] = 0;
            }
        }
    }
    for (VertexID i = 0; i < n; ++i)
        vTables[36][i] -= vTables[13][i];
    delete[] nborEdgeOutCount;
    delete[] wedgeCount;
}

void orbit66(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
             std::vector<HashTable> &eTables) {
    ui n = dun.getNumVertices();
    VertexID *triends = new VertexID[n]; // array to store triangle ends
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    EdgeID edge_pos_ij, edge_pos_ik, edge_pos_il, edge_pos_jk, edge_pos_jl, edge_pos_kl;

    VertexID *edge_pos_i_triend = new VertexID[n]; // array to store triangle ends
    VertexID *edge_pos_j_triend = new VertexID[n]; // array to store triangle ends

    for (VertexID i=0; i < n; ++i) // loop over vertices
    {
        for (VertexID posj = outOffset[i]; posj < outOffset[i+1]; ++posj) // loop over out-neighbors of i
        {
            VertexID j = outNbors[posj]; // j is current out-neighbor
            edge_pos_ij = posj;
            VertexID count = 0;
            for (VertexID posk = posj+1; posk < outOffset[i+1]; ++posk) // loop over another out-neighbor of i, that is "ahead" of j in list of out-neighbors
            {
                VertexID k = outNbors[posk]; // k is next out-neighbor
                edge_pos_jk = dout.getEdgeID(j, k);
                if (edge_pos_jk != -1) // check if edge (j,k) is present
                {
#ifdef COLLECT_STATISTICS
                    ++gNumMatch;
#endif
                    triends[count] = k;  // so (i,j,k) form a triangle. we store the fact that k forms a triangle with edge (i,j) in digraph gout
                    edge_pos_i_triend[count] = posk;
                    edge_pos_j_triend[count] = edge_pos_jk;
                    ++count;
                }
            }

            for (VertexID posk = 0; posk < count; ++posk) // loop over all pairs of triangles formed by (i,j)
            {
                VertexID k = triends[posk]; // k is vertex as index posk in triends
                edge_pos_ik = edge_pos_i_triend[posk];
                edge_pos_jk = edge_pos_j_triend[posk];
                // We will search all other vertices in triends in k's adj list
                for (VertexID posell = posk+1; posell < count; ++posell)
                {
                    VertexID ell = triends[posell];
                    edge_pos_il = edge_pos_i_triend[posell];
                    edge_pos_jl = edge_pos_j_triend[posell];
                    edge_pos_kl = dout.getEdgeID(k, ell);
                    if (edge_pos_kl != -1) // (k,ell) is present, thus (i,j,k,ell) form a 4-clique
                    {
#ifdef COLLECT_STATISTICS
                        ++gNumMatch;
#endif
                        vTables[66][i] += eTables[1][edge_pos_jk]-2;
                        vTables[66][i] += eTables[1][edge_pos_jl]-2;
                        vTables[66][i] += eTables[1][edge_pos_kl]-2;
                        vTables[66][j] += eTables[1][edge_pos_ik]-2;
                        vTables[66][j] += eTables[1][edge_pos_il]-2;
                        vTables[66][j] += eTables[1][edge_pos_kl]-2;
                        vTables[66][k] += eTables[1][edge_pos_ij]-2;
                        vTables[66][k] += eTables[1][edge_pos_il]-2;
                        vTables[66][k] += eTables[1][edge_pos_jl]-2;
                        vTables[66][ell] += eTables[1][edge_pos_ij]-2;
                        vTables[66][ell] += eTables[1][edge_pos_ik]-2;
                        vTables[66][ell] += eTables[1][edge_pos_jk]-2;
                    }
                }
            }
        }
    }

    delete[] triends;
    delete[] edge_pos_i_triend;
    delete[] edge_pos_j_triend;
}

void chordalWedge(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
                  const Triangle &tri, EdgeID *unID) {
    ui n = dun.getNumVertices();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    const EdgeID *unOffset = dun.getOffsets();
    const VertexID *unNbors = dun.getNbors();
    HashTable wedgeCount = new Count[n];
    HashTable diamondCount = new Count[n];
    HashTable diamondCountCopy = new Count[n];
    HashTable diamondCountI = new Count[n];
    HashTable midDiamondCount = new Count[n];
    memset(wedgeCount, 0, sizeof(Count) * n);
    memset(diamondCount, 0, sizeof(Count) * n);
    memset(diamondCountCopy, 0, sizeof(Count) * n);
    memset(diamondCountI, 0, sizeof(Count) * n);
    memset(midDiamondCount, 0, sizeof(Count) * n);
    ui *keyPos1 = new ui[n / 8 + 1];
    ui *keyPos2 = new ui[n / 8 + 1];
    ui keyPosSize1 = 0, keyPosSize2 = 0;
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = unOffset[i]; ij < unOffset[i + 1]; ++ij) {
            VertexID j = unNbors[ij];
            for (EdgeID jk = unOffset[j]; jk < unOffset[j + 1]; ++jk) {
                VertexID k = unNbors[jk];
                if (k <= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                ++wedgeCount[k];
                if (keyPosSize1 < n / 8) {
                    keyPos1[keyPosSize1] = k;
                    ++keyPosSize1;
                }
            }
        }
        for (EdgeID ij = unOffset[i]; ij < unOffset[i + 1]; ++ij) {
            VertexID j = unNbors[ij];
            VertexID *candidateK;
            ui cntK = 0;
            candidateK = tri.getTriend(unID[ij], 4, cntK);
            for (ui posK = 0; posK < cntK; ++posK) {
                VertexID k = candidateK[posK];
                EdgeID jk = dout.getEdgeID(j, k);
                VertexID *candidateL;
                ui cntL = 0;
                candidateL = tri.getTriend(jk, 4, cntL);
                // i < l
                ui startPosL = firstPosGreaterThan(candidateL, 0, cntL, i);
                for (ui posL = startPosL; posL < cntL; ++posL) {
                    VertexID l = candidateL[posL];
#ifdef COLLECT_STATISTICS
                    ++gNumMatch;
#endif
                    ++diamondCount[l];
                    diamondCountCopy[l] = diamondCount[l];
                    if (j < k) {
                        ++diamondCountI[l];
                        if (keyPosSize2 < n / 8) {
                            keyPos2[keyPosSize2] = l;
                            ++keyPosSize2;
                        }
                        ++midDiamondCount[j];
                        ++midDiamondCount[k];
                        vTables[63][i] += wedgeCount[l] - 2;
                        vTables[63][l] += wedgeCount[l] - 2;
                        vTables[64][j] += wedgeCount[l] - 2;
                        vTables[64][k] += wedgeCount[l] - 2;
                    }
                }
            }
            for (ui posK = 0; posK < cntK; ++posK) {
                VertexID k = candidateK[posK];
                EdgeID jk = dout.getEdgeID(j, k);
                VertexID *candidateL;
                ui cntL = 0;
                candidateL = tri.getTriend(jk, 4, cntL);
                // i < l
                ui startPosL = firstPosGreaterThan(candidateL, 0, cntL, i);
                for (ui posL = startPosL; posL < cntL; ++posL) {
                    VertexID l = candidateL[posL];
#ifdef COLLECT_STATISTICS
                    ++gNumMatch;
#endif
                    vTables[69][j] += choosec(diamondCount[l], 2);
                    vTables[68][i] += choosec(diamondCount[l], 2);
                    vTables[68][l] += choosec(diamondCount[l], 2);
                    vTables[68][k] += diamondCountCopy[l] - 1;
                    diamondCount[l] = 0;
                }
            }
        }
        for (EdgeID ij = unOffset[i]; ij < unOffset[i + 1]; ++ij) {
            VertexID j = unNbors[ij];
            for (EdgeID jk = unOffset[j]; jk < unOffset[j + 1]; ++jk) {
                VertexID k = unNbors[jk];
                if (k <= i) continue;
#ifdef COLLECT_STATISTICS
                ++gNumMatch;
#endif
                vTables[62][j] += diamondCountI[k];
            }
        }
        for (EdgeID ij = unOffset[i]; ij < unOffset[i + 1]; ++ij) {
            VertexID j = unNbors[ij];
            vTables[62][j] -= midDiamondCount[j];
            midDiamondCount[j] = 0;
        }
        if (keyPosSize1 == n / 8) memset(wedgeCount, 0, sizeof(Count) * n);
        else {
            for (int pos = 0; pos < keyPosSize1; ++pos)
                wedgeCount[keyPos1[pos]] = 0;
        }
        if (keyPosSize2 == n / 8) memset(diamondCountI, 0, sizeof(Count) * n);
        else {
            for (int pos = 0; pos < keyPosSize2; ++pos)
                diamondCountI[keyPos2[pos]] = 0;
        }
        keyPosSize1 = 0;
        keyPosSize2 = 0;
    }

    for (VertexID v = 0; v < n; ++v) {
        vTables[68][v] /= 2;
        vTables[69][v] /= 2;
    }
    delete[] wedgeCount;
    delete[] diamondCount;
    delete[] diamondCountCopy;
    delete[] diamondCountI;
    delete[] midDiamondCount;
    delete[] keyPos1;
    delete[] keyPos2;
}

void almostFiveClique(const DataGraph &dun, const DataGraph &din, const DataGraph &dout, std::vector<HashTable> &vTables,
                      const Triangle &tri) {
    ui n = dun.getNumVertices();
    const EdgeID *outOffset = dout.getOffsets();
    const VertexID *outNbors = dout.getNbors();
    VertexID *candidateL = new VertexID[n];
    // get a directed triangle
    for (VertexID i = 0; i < n; ++i) {
        for (EdgeID ij = outOffset[i]; ij < outOffset[i + 1]; ++ij) {
            VertexID j = outNbors[ij];
            ui candCount;
            VertexID *candidateK = tri.getTriend(ij, 3, candCount);
            for (int pos = 0; pos < candCount; ++pos) {
                VertexID k = candidateK[pos];
                EdgeID ik = dout.getEdgeID(i, k);
                EdgeID jk = dout.getEdgeID(j, k);
                // get undirected triangle of ik and jk
                ui cnt1, cnt2, fourClique;
                VertexID *triangle1 = tri.getTriend(ik, 4, cnt1);
                VertexID *triangle2 = tri.getTriend(jk, 4, cnt2);
                ComputeSetIntersection::ComputeCandidates(triangle1, cnt1, triangle2, cnt2, candidateL, fourClique);
#ifdef COLLECT_STATISTICS
                gNumMatch += fourClique;
#endif
                vTables[71][i] += choosec(fourClique, 2);
                vTables[71][j] += choosec(fourClique, 2);
                vTables[71][k] += choosec(fourClique, 2);
                for (int pos2 = 0; pos2 < fourClique; ++pos2)
                    vTables[70][candidateL[pos2]] += fourClique - 1;
            }
        }
    }

    delete[] candidateL;
}

int main(int argc, char **argv) {
    Command cmd(argc, argv);
    std::string dataGraphPath = cmd.getDataGraphPath();
    std::string resultPath = cmd.getResultPath();
    std::string trianglePath = cmd.getTrianglePath();
    DataGraph dun = DataGraph();
    dun.loadDataGraph(dataGraphPath);
    const DataGraph din = constructDirectedDataGraph(dun, false);
    const DataGraph dout = constructDirectedDataGraph(dun, true);
    specialsparse *sg = (specialsparse *)malloc(sizeof(specialsparse));
    dout.initSpecialSparse(sg);
    ui n = dun.getNumVertices(), m = dun.getNumEdges();
    Triangle triangle;
    triangle.load(trianglePath, m / 2);
    EdgeID *outID = buildInID2OutID(din, dout);
    EdgeID *unID = buildUnID2OutID(dun, dout);
    std::vector<HashTable> vTables(73);
    std::vector<HashTable> eTables(12);
    std::vector<int> eIndexes = {1, 5, 7, 9, 11};
    std::vector<HashTable> revETables(12);
    for (HashTable &h: vTables) {
        h = new Count[n];
        memset(h, 0, sizeof(Count) * n);
    }
    for (int i : eIndexes) {
        eTables[i] = new Count[m / 2];
        memset(eTables[i], 0, sizeof(Count) * m / 2);
    }
    revETables[9] = new Count[m / 2];
    memset(eTables[9], 0, sizeof(Count) * m / 2);
    const EdgeID *unOffset = dun.getOffsets();
    const VertexID *unNbors = dun.getNbors();
    for (VertexID i = 0; i < n; ++i) {
        Count degree = unOffset[i + 1] - unOffset[i];
        vTables[0][i] = degree;
        for (EdgeID ij = unOffset[i]; ij < unOffset[i + 1]; ++ij) {
            VertexID j = unNbors[ij];
            Count degJ = unOffset[j + 1] - unOffset[j];
            vTables[1][i] += degJ - 1;
            VertexID *triend;
            ui cnt;
            triend = triangle.getTriend(unID[ij], 4, cnt);
            ui startPos = firstPosGreaterThan(triend, 0, cnt, j);
            vTables[3][i] += cnt - startPos;
            eTables[1][unID[ij]] = cnt;
        }
        vTables[2][i] = choosec(degree, 2);
    }

    double totalTime = 0.0;
    auto start = std::chrono::steady_clock::now();
    otherFourOrbit(dun, din, dout, vTables, eTables, revETables);
    fourCycle(dun, din, dout, vTables, eTables, outID);
    fourCliqueOrbit(dun, din, dout, vTables, eTables);
    fiveTree(dun, din, dout, vTables);
    otherFiveOrbits(dun, din, dout, vTables, eTables, revETables, outID);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    start = std::chrono::steady_clock::now();
    fiveCycle(dun, din, dout, vTables, triangle, outID);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    std::cout << "Done: counting five cycles, time: " << elapsedSeconds.count() << std::endl;
    start = std::chrono::steady_clock::now();
    wedgeCollision(dun, din, dout, vTables);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    std::cout << "Done: counting orbit 49 and 50, time: " << elapsedSeconds.count() << std::endl;
    start = std::chrono::steady_clock::now();
    orbit51(dun, din, dout, vTables, eTables, outID);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    std::cout << "Done: counting orbit 51, time: " << elapsedSeconds.count() << std::endl;
    start = std::chrono::steady_clock::now();
    orbit36(dun, din, dout, vTables);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    std::cout << "Done: counting orbit 36, time: " << elapsedSeconds.count() << std::endl;
    start = std::chrono::steady_clock::now();
    chordalWedge(dun, din, dout, vTables, triangle, unID);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    std::cout << "Done: counting 62,63,64,68,69, time: " << elapsedSeconds.count() << std::endl;
    start = std::chrono::steady_clock::now();
    orbit66(dun, din, dout, vTables, eTables);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    std::cout << "Done: counting orbit 66, time: " << elapsedSeconds.count() << std::endl;
    start = std::chrono::steady_clock::now();
    almostFiveClique(dun, din, dout, vTables, triangle);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    totalTime += elapsedSeconds.count();
    std::cout << "Done: counting 70, 71, time: " << elapsedSeconds.count() << std::endl;
    start = std::chrono::steady_clock::now();
    VertexID *cliqueVertices = new VertexID[5];
    mkspecial(sg, 5);
    kclique(5, 5, sg, cliqueVertices, vTables[72], 1);
    freesub(sg, 5);
    end = std::chrono::steady_clock::now();
    elapsedSeconds = end - start;
    std::cout << "Done: counting 72 time: " << elapsedSeconds.count() << std::endl;
    totalTime += elapsedSeconds.count();
    std::cout << "total time: " << totalTime << std::endl;
    std::cout << "number of matches enumerated: " << gNumMatch;
    std::vector<int> orbitTypes(73, 1);
    std::vector<std::string> files(73);
    for (int i = 0; i < 73; ++i) {
        files[i] = std::to_string(i) + ".txt";
    }
    if (!resultPath.empty()) saveCount(resultPath, vTables, dun, true, files, orbitTypes);
    return 0;
}