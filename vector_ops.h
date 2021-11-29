#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <iostream>
#include <vector>
#include <string>
#include <utility>

using namespace std;

//A VECTOR OF THE CURVES READ FROM THE INPUT FILE
static vector<pair<pair<string, int>, vector<double>>> curve_vector;

//EVERY CURVE HAS A COLLECTION OF L ID VALUES
static vector<vector<int>> curves_ID_vector;


/*-------------------------CURVE VECTOR FUNCTIONS------------------------*/
void curve_vector_insert_curve(pair<pair<string, int>, vector<double>>& curve);
pair<pair<string, int>, vector<double>> curve_vector_get_curve(int curve_id);
int curve_vector_get_size();
void curve_vector_print_values();

/*--------------------POINTS ID VECTOR FUNCTIONS--------------------------*/
void curves_ID_vector_initialize(int num_of_curves, int L);
void curves_ID_vector_insert(int index_value, vector<int>& curve_id_values);
int curves_ID_vector_get_curve_value(int index_value, int k);





#endif