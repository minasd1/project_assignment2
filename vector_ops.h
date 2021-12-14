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
#include <initializer_list>

#include "hash_functions.h"
#include "frechet.hpp"
#include "curve.hpp"
#include "point.hpp"

using namespace std;

//A VECTOR OF THE CURVES READ FROM THE INPUT FILE
static vector<pair<pair<string, int>, vector<double>>> curve_vector;

//EVERY CURVE HAS A COLLECTION OF L ID VALUES
static vector<vector<int>> curves_ID_vector;

//VECTOR OF CENTROID POINTS - USED IN CLUSTERING
static vector<int> centroids;


/*-------------------------CURVE VECTOR FUNCTIONS------------------------*/
void curve_vector_insert_curve(pair<pair<string, int>, vector<double>>& curve);
pair<pair<string, int>, vector<double>> curve_vector_get_curve(int curve_id);
int curve_vector_get_size();
void curve_vector_print_values();

/*--------------------CURVES ID VECTOR FUNCTIONS--------------------------*/
void curves_ID_vector_initialize(int num_of_curves, int L);
void curves_ID_vector_insert_lsh(int index_value, vector<int>& curve_id_values);
void curves_ID_vector_insert_frechet(int index_value, int curve_id_value, int grid);
int curves_ID_vector_get_curve_value(int index_value, int k);
vector<int> curves_ID_vector_get_curve_all_ids(int curve_id);
void curves_ID_vector_print();

/*-----------------------CENTROID FUNCTIONS-------------------------------*/
void centroids_insert_curve(int id);
int centroids_get_centroid(int index);
int centroids_get_size();
vector<int> centroids_get_table();
void centroids_pick_first_centroid();
void centroids_pick_next_centroid(vector<float>& partial_sums);
float centroids_calculate_min_distance_curve(vector<double>& curve);
void centroids_calculate_min_distance_input(vector<float>& curves_min_distances);
void centroids_print_data();

/*-------------------------V_VECTOR FUNCTIONS----------------------------*/
void v_vectors_initialization(vector<vector<int>>& v_vectors, int num_of_v_vectors, int dimensions);
void v_vectors_assign_coordinances(vector<vector<int>>& v_vectors, int num_of_vectors, int dimensions,
                                   std::default_random_engine& generator);
void v_vectors_printdata(vector<vector<int>>& v_vectors);

/*------------------------T-VECTOR FUNCTIONS------------------------------*/
void create_vector_t(vector<float>& t, int k, int w, std::default_random_engine& generator);
void print_vector_t(vector<float>& t);
void create_vector_int(vector<int>& ints, int k, int w, std::default_random_engine& generator);

/*--------------------OPERATIONS BETWEEN CURVES FUNCTIONS-----------------*/
double curve_calculate_dfd(const pair<pair<string, int>, vector<double>>& curve1, const pair<pair<string, int>, vector<double>>& curve2);
Curve convert_for_continuous_frechet(const pair<pair<string, int>, vector<double>>& curve, const unsigned long dimensions);

/*--------------------OPERATIONS BETWEEN VECTORS FUNCTIONS-----------------*/
int calculate_dot_product(const pair<pair<string, int>, vector<double>>& curve, vector <int>& d_vector);
vector<double> add_vectors(const pair<pair<string, int>, vector<double>>& curve1, const pair<pair<string, int>, vector<double>>& curve2);
double calculate_distance(vector<double>& point1, const vector<double>& point2, int k=2);

/*-------------------------OTHER FUNCTIONS---------------------------------*/
float calculate_partial_sums(vector<float>& min_distances, vector<float>& partial_sums);



#endif
