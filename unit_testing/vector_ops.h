#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <utility>

using namespace std;

/*--------------------OPERATIONS BETWEEN VECTORS FUNCTIONS-----------------*/
int calculate_dot_product(const pair<pair<string, int>, vector<double>>& curve, vector <int>& d_vector);
vector<double> add_vectors(const vector<double>& curve1, const vector<double>& curve2);
double calculate_distance(vector<double>& curve1, const vector<double>& curve2, int k=2);
bool non_zero_coordinates(vector<double>& coordinates);
pair<pair<string, int>, vector<double>> get_mean_curve_vector(vector<double> vector_of_sums, int num_of_vectors, int& last_known_id);

#endif
