#ifndef LLOYDS_AUXILIARY_H
#define LLOYDS_AUXILIARY_H
#include <vector>

using namespace std;

double max_double(double& a, double& b);
bool already_in_that_cluster(vector<vector<int>> cluster_table, int index, int id);
int find_nearest_centroid(pair<pair<string, int>, vector<double>>& current_curve);
int find_second_nearest_centroid(pair<pair<string, int>, vector<double>>& current_curve);
void print_cluster_table (vector<vector<int>>& table);
#endif
