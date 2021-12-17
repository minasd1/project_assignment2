#ifndef HASH_FUNCTIONS_H
#define HASH_FUNCTIONS_H

#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <cmath>
#include <algorithm>
#include <utility>

using namespace std;

typedef std::default_random_engine engine;

class G_Lsh {
private:
    int k; // THE NUMBER OF H FUNCTIONS USED IN G FUNCTIONS
    int h_num; // THE TOTAL NUMBER OF H FUNCTIONS TO BE CREATED, MULTIPLE OF k
    int dimension; // THE NUMBER OF COORDINATES EACH V_VECTOR AND POINT IN THE DATASET HAS
    engine generator; // THE GENERATOR OF PSEUDO-RANDOM NUMBERS
    vector<vector<int>> v_vectors; // A STRUCT TO STORE THE RANDOMLY GENERATED V_VECTORS (USED IN H FUNCTION)
    int num_of_vectors; // THE SAME WITH h_num
    vector<float> t; // A STRUCT TO STORE THE RANDOM FLOATS t (USED IN H FUNCTION)
    int w; // THE WINDOW VARIABLE (USED IN H FUNCTION)
    int m; // THE FIRST MODULO IN G FUNCTIONS
    int table_size; // THE SECOND MODULO IN G FUNCTION
    int l; // THE NUMBER OF G FUNCTIONS (OR THE NUMBER OF HASH TABLES)
    vector<vector<int>> h_functions; //A STRUCT TO STORE WHICH H FUNCTIONS ARE USED IN EACH G FUNCTION
    vector<vector<int>> r; //A STRUCT TO STORE THE r RANDOM INTEGERS OF EACH G FUNCTION
public:
    G_Lsh(int k_num, int dim, engine gen, int win, int m_mod, int tab_s, int l_num);
    void id(const pair<pair<string, int>, vector<double>>& curve, vector<int>& id_vector, bool is_query, int frechet_grid);
    void hash(const pair<pair<string, int>, vector<double>>& curve, vector<int>& hash_vector, bool is_query, int frechet_grid);
    void print_hash_functions_data(void);
};

class G_Hypercube {
private:
    vector<vector<int>> v_vectors; // A STRUCT TO STORE THE RANDOMLY GENERATED V_VECTORS (USED IN H FUNCTION)
    vector<float> t; // A STRUCT TO STORE THE RANDOM FLOATS t (USED IN H FUNCTION)
    int dimension; // THE NUMBER OF COORDINATES EACH V_VECTOR AND POINT IN THE DATASET HAS
    engine generator; // THE GENERATOR OF PSEUDO-RANDOM NUMBERS
    int string_dimensions; // d'   THE NUMBER OF BITS OF THE STRING PRODUCED BY THE HASHFUNCTION
    int w; // THE WINDOW VARIABLE (USED IN H FUNCTION)
    vector<int> rand_ints; // A TABLE OF d' RANDOM INTEGER VALUES, ONE FOR EACH f FUNCTION
public:
    G_Hypercube(int dim, engine gen, int win, int str_dim);
    void hash(const pair<pair<string, int>, vector<double>>& curve, unsigned int& hash_value, bool is_query);
};

class G_Frechet{
    private:
        G_Lsh g;
        vector<vector<float>> t;        //A STRUCT TO STORE THE RANDOM DOUBLE T (USED IN GRID CURVES)
        engine generator;               //THE GENERATOR OF PSEUDO RANDOM NUMBERS
        int L;                          //NUMBER OF GRID CURVES WHEN LSH IS USED
        double delta;                   //SPACE BETWEEN GRID CURVE VALUES
        int num_of_curve_values;        //MAXIMUM NUMBER OF VALUES USED IN GIRD CURVE
        double max_coordinate_value;    //MAXIMUM VALUE THAT IS RECORDED FOR ALL COORDINATES
        int extra_values_factor;        //NEEDED WHEN FRECHET CLUSTERING IS USED
        bool flag_frechet_cluster;  
        int max_length;                 //THE MAXIMUM SIZE OF A MEAN CURVE - CLUSTERING       
    
    public:
        G_Frechet(G_Lsh g, engine gen, int L_num, double delta_value, double max_value, int num_of_curve_values, bool flag_fc);
        void hash(const pair<pair<string, int>, vector<double>>& curve, vector<int>& hash_vector, vector<int>& id_vector, bool is_query, bool is_mean, int grid_dimensions);
        void filter(const vector<double>& curve, vector<double>& filtered_curve, double epsilon, bool is_mean);
        void snap_to_grid(const pair<pair<string, int>, vector<double>>& curve, vector<double>& snapped_curve, int grid_dimensions, int grid_num);
        void padding(vector<double>& snapped_curve, int grid_dimensions, int extra_values_factor);
        void minima_maxima(vector<double>& snapped_curve);

};



#endif
