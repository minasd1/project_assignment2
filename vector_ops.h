#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <iostream>
#include <vector>
#include <string>
#include <utility>

using namespace std;


static vector<pair<pair<string, int>, vector<double>>> curve_vector;


/*-------------------------CURVE VECTOR FUNCTIONS------------------------*/
void curve_vector_insert_curve(pair<pair<string, int>, vector<double>>& curve);
pair<pair<string, int>, vector<double>> curve_vector_get_curve(int curve_id);
int curve_vector_get_size();
void curve_vector_print_values();





#endif