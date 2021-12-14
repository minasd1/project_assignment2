#ifndef CUBE_H
#define CUBE_H

#include <iostream>
#include <utility>

#include "vector_ops.h"
#include "hypercube.h"
#include "hash_functions.h"
#include "knn_table_functions.h"

using namespace std;

vector<int> cube_range_search(int g, int radius, int probes, int dimensions, pair<pair<string, int>, vector<double>>& query_curve);
vector<dist_id_pair> cube_find_approximate_knn(pair<pair<string, int>, vector<double>>& query_curve, int k,  G_Hypercube& g, int probes, int dimensions, int& count_nn, int max_candidates=-1);












#endif
