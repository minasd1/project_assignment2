#ifndef FRECHET_FUNCTIONS_H
#define FRECHET_FUNCTIONS_H

#include <iostream>
#include <vector>
#include <algorithm>

#include "vector_ops.h"
#include "hash_functions.h"
#include "hashTable.h"
#include "lsh.h"
#include "knn_table_functions.h"

using namespace std;


vector<dist_id_pair> frechet_find_approximate_knn_discrete(pair<pair<string, int>, vector<double>>& query_curve, int k,  G_Frechet& g, int max_candidates=-1);
vector<dist_id_pair> frechet_find_approximate_knn_continuous(pair<pair<string, int>, vector<double>>& query_curve, int k,  G_Frechet& g, int frechet_grid, int max_candidates=-1);
vector<dist_id_pair> frechet_find_exact_knn_discrete(pair<pair<string, int>, vector<double>>& query_curve, int k, int num_of_points);







#endif