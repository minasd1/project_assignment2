#ifndef LSH_H
#define LSH_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

#include "vector_ops.h"
#include "hash_functions.h"
#include "hashTable.h"
#include "knn_table_functions.h"

using namespace std;

vector<dist_id_pair> lsh_find_approximate_knn(pair<pair<string, int>, vector<double>>& query_curve, int k,  G_Lsh& g, int frechet_grid, int max_candidates=-1);
vector<dist_id_pair> find_exact_knn(pair<pair<string, int>, vector<double>>& query_curve, int k, int num_of_points);
vector<int> lsh_range_search(vector<int>& g, int radius, pair<pair<string, int>, vector<double>>& query_curve);
bool is_there_someone_with_same_id(vector<int>& query_curve_ids, vector<int>& candidate_curves, /*G_Lsh& g,*/
                                   int& same_id_counter, vector<int>& curves_with_same_id);








#endif
