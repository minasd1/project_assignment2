#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <cmath>
#include <chrono>
#include <limits>
#include <algorithm>

#include "hash_functions.h"

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
void curves_ID_vector_insert_lsh(int index_value, vector<int>& curve_id_values);
void curves_ID_vector_insert_frechet(int index_value, int curve_id_value, int grid);
int curves_ID_vector_get_curve_value(int index_value, int k);

/*-------------------------V_VECTOR FUNCTIONS----------------------------*/
void v_vectors_initialization(vector<vector<int>>& v_vectors, int num_of_v_vectors, int dimensions);
void v_vectors_assign_coordinances(vector<vector<int>>& v_vectors, int num_of_vectors, int dimensions,
                                   std::default_random_engine& generator);
void v_vectors_printdata(vector<vector<int>>& v_vectors);

/*------------------------T-VECTOR FUNCTIONS------------------------------*/
void create_vector_t(vector<float>& t, int k, int w, std::default_random_engine& generator);
void print_vector_t(vector<float>& t);
void create_vector_int(vector<int>& ints, int k, int w, std::default_random_engine& generator);

/*--------------------OPERATIONS BETWEEN VECTORS FUNCTIONS-----------------*/
int calculate_dot_product(const pair<pair<string, int>, vector<double>>& curve, vector <int>& d_vector);
vector<double> add_vectors(const pair<pair<string, int>, vector<double>>& curve1, const pair<pair<string, int>, vector<double>>& curve2);
double calculate_distance(vector<double>& point1, const vector<double>& point2, int k=2);





#endif
