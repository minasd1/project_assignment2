#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <random>

#include "vector_ops.h"
#include "file_functions.h"
#include "hashTable.h"
#include "hash_functions.h"
#include "frechet.hpp"
#include "curve.hpp"
#include "lsh.h"
#include "hypercube.h"
#include "cube.h"
#include "frechet_functions.h"
#include "read_file_data.h"


int main(int argc, char* argv[]){

    int k = 4;                      //NUMBER OF H FUNCTIONS USED IN FUNCTION G
    int L = 5;                      //NUMBER OF HASHTABLES(LSH)
    int N=1 ;                         //NUMBER OF NEAREST MEIGHBORS
    float R ;                       //SEARCH RADIUS
    int probes=2 ;                    //MAX NUMBER OF HYPERCUBE VERTICES TO BE CHECKED
    int k_cube=3;                     //D'
    int M_cube=1000 ;                    //MAX NUMBER OF CANDIDATE POINTS TO BE CHECKED
    int window= 50;                
    int k_cluster ;                 //NUMBER OF CENTROIDS - CLUSTERING
    double delta = 2.5;
    double epsilon = 0.1;
    bool complete= false;

    fstream input_file;             //FILE WE READ INPUT FROM
    fstream query_file;             //FILE WE READ QUERIES FROM
    fstream config_file;            //CLUSTER CONFIGURATION FILE
    fstream output_file;            //FILE TO WRITE OUTPUT TO
    string line;
    string input_file_name, query_file_name, output_file_name, config_file_name;
    int start;
    int finish = 0;
    bool first_iteration = true;
    string token;
    unsigned int number_of_curves = 0;
    int buckets;
    int num_of_curve_values = 0;
    int query_num_of_curve_values = 0;
    int curves_divider = 16;        //USED TO GET TOTAL CURVES IN EACH HASHTABLE
    int continue_execution = 1;
    string name;
    int id;
    unsigned int hash_value;
    vector<double> values;
    vector<int> hash_vector;
    vector<int> id_vector;
    int M = pow(2, 31) - 5;
    int count = 0;
    double max_value = 0.0;
    double max_coordinate_value;
    double query_max_value = 0.0;
    double query_max_coordinate_value;
    bool is_query_curve = true;
    string algorithm = "Frechet"; 
    string metric = "continuous";

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    input_file_name = argv[1];      //NEEDS TO BE CHANGED(BY USING CMD_LINE_ARGS)

    open_file(&input_file, input_file_name, fstream::in);   //OPEN INPUT FILE

    while(getline(input_file, line)){                       //READ FILE LINE BY LINE
       
        pair<pair<string, int>, vector<double>> curve;      //CURRENT CURVE
        pair<string, int> name_and_id;          
        vector<double> curve_values;

        start = 0;
        while(start < line.size()){                         //TOKENIZE EVERY LINE IN IT'S SEPERATED STRINGS
            finish = line.find_first_of('\t', start);

            if(finish == string::npos){

                finish = line.size();
            }

            if(start < line.size() - 1){
                token = line.substr(start, finish - start);
                if(start == 0){

                    name_and_id = make_pair(token, count);  //ASSIGN NAME TO AN INCREASING ID
                }
                else{
                    //KEEP TRACK OF THE MAX VALUE RECORDED IN FILE
                    if((stod(token)) > max_value){

                        max_value = stod(token);
                    }
                    curve_values.push_back(stod(token));    //CONVERT THE STRING TO DOUBLE AND PASS THEM TO CURVE VALUES
                }
                
            }

            start = finish + 1;
            if(first_iteration == true){
                num_of_curve_values++;                      //KEEP TRACK OF THE NUMBER OF CURVE VALUES OF GIVEN INSTANCE
            }
        }

        //INSERT CURRENT CURVE TO CURVE VECTOR
        curve = make_pair(name_and_id, curve_values);
        curve_vector_insert_curve(curve);

        number_of_curves++;
        first_iteration = false;

        count++;
    }

    max_coordinate_value = max(max_value, (double)num_of_curve_values);
    
    curves_ID_vector_initialize(number_of_curves, L);

    //NUMBER OF BUCKETS IN EACH HASHTABLE
    buckets = number_of_curves/curves_divider;

    //INITIALIZE G FUNCTION THAT LEADS US TO LSH HASHTABLE BUCKETS
    G_Lsh g_lsh(k, num_of_curve_values, generator, window, M, buckets, L);

    //INITIALIZE G FUNCTION THAT LEADS US TO HYPERCUBE BUCKETS
    G_Hypercube g_cube(num_of_curve_values, generator, window, k_cube);

    query_file_name = argv[2];

    if(strcmp(argv[0], "./search") == 0){

        query_max_coordinate_value = file_get_max_value(query_file, query_file_name); 
    }
    
    max_coordinate_value = max(max_coordinate_value, query_max_coordinate_value);

    //INITIALIZE G FRECHET FUNCTION THAT USES g_lsh
    G_Frechet g_frechet(g_lsh, generator, L, delta, max_coordinate_value);

    if(strcmp(argv[0], "./search") == 0){

        //INITIALIZE L HASHTABLES WITH HASHTABLESIZE BUCKETS AND ZERO POINTS IN EACH BUCKET
        hashTable_initialization(L, buckets);

        if((algorithm == "Frechet")){

            //INSERT ALL INPUT CURVES TO THE HASHTABLES
            for(int i = 0; i < number_of_curves; i++){
                if(metric == "discrete"){

                    g_frechet.hash(curve_vector_get_curve(i), hash_vector, id_vector, !is_query_curve, 2);
                }
                else if(metric == "continuous"){

                    g_frechet.hash(curve_vector_get_curve(i), hash_vector, id_vector,  !is_query_curve, 1);
                }    
                
            }
            
            
            hash_vector.clear();
        }
        else if (algorithm == "LSH") {
            //INSERT ALL INPUT CURVES TO THE HASHTABLES
            for(int i = 0; i < number_of_curves; i++){

                g_lsh.hash(curve_vector_get_curve(i), hash_vector, 0, 0);
            }
            hash_vector.clear();
        }
        else if (algorithm == "Hypercube") {
            //INITIALIZE A HYPERCUBE WITH 2^D' BUCKETS AND ZERO POINTS IN EACH BUCKET
            hyperCube_initialization(pow(2, k_cube));

            //INSERT ALL THE POINTS FROM INPUT FILE TO HYPERCUBE BUCKETS
            for(int i = 0; i < number_of_curves; i++){
                g_cube.hash(curve_vector_get_curve(i), hash_value, 0);
            }
        }

    }

    output_file_name = argv[3];

    query_file.clear();
    query_file.seekg(0, query_file.beg);

    while(continue_execution == 1){

        if(strcmp(argv[0], "./search") == 0){

            //OPEN FILE TO WRITE RESULTS TO
            open_file(&output_file, output_file_name, fstream::out);

            finish = 0;

            //FOR EVERY QUERY CURVE
            while(getline(query_file, line)){                       //READ QUERY FILE LINE BY LINE

                pair<pair<string, int>, vector<double>> query_curve;
                pair<string, int> name_and_id;          
                vector<double> query_curve_values;

                start = 0;

                while(start < line.size()){                         //TOKENIZE EVERY LINE IN IT'S SEPERATED STRINGS
                    finish = line.find_first_of('\t', start);

                    if(finish == string::npos){

                        finish = line.size();
                    }

                    if(start < line.size() - 1){
                        token = line.substr(start, finish - start);
                        if(start == 0){

                            name_and_id = make_pair(token, -1);         //ASSIGN NAME TO AND ID
                        }
                        else{
                            query_curve_values.push_back(stod(token));    //CONVERT THE STRING TO DOUBLE AND PASS THEM TO CURVE VALUES
                        }
                        
                    }

                    start = finish + 1;
                }

                query_curve = make_pair(name_and_id, query_curve_values);

                if(algorithm == "LSH"){
                    g_lsh.hash(query_curve, hash_vector, 1, 0);
                    vector<dist_id_pair> curves_lsh, curves_brute;

                    //LSH NEAREST NEIGHBORS
                    //FIND TIME LSH - APPROXIMATE NEIGHBORS
                    
                    auto start_time = std::chrono::high_resolution_clock::now();
                    curves_lsh= lsh_find_approximate_knn(query_curve, N, g_lsh, 0);
                    auto stop_time = std::chrono::high_resolution_clock::now();
                    auto time_lsh = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
                    //FIND TIME BRUTE FORCE - EXACT NEIGHBORS
                   
                    auto start_time2 = std::chrono::high_resolution_clock::now();
                    curves_brute= find_exact_knn(query_curve, N, number_of_curves);
                    auto stop_time2 = std::chrono::high_resolution_clock::now();
                    auto time_brute = std::chrono::duration_cast<std::chrono::microseconds>(stop_time2 - start_time2);
                    //PRINTING IN OUTPUT FILE
                    output_file << "Query: " << query_curve.first.first << endl;
                    output_file << "Algorithm: " << algorithm << endl;
                    for (int i= 1; i <= N ; i++) {
                        if (i > curves_lsh.size()) {
                            output_file << "Nearest neighbor-" << i<< ": Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            output_file << "distanceLSH: Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            continue;
                        }
                        output_file << "Nearest neighbor-" << i<< ": " << curves_lsh[i-1].id << endl;
                        output_file << "distanceLSH: " << curves_lsh[i-1].dist << endl;
                        output_file << "distanceTrue: " << curves_brute[i-1].dist << endl;
                        
                    }
                    output_file <<  "tLSH: " << time_lsh.count() << " microseconds" << endl;
                    output_file << "tTrue: " << time_brute.count() << " microseconds" << endl;
                    output_file << endl;
                }
                else if(algorithm == "Hypercube"){
                    g_cube.hash(query_curve, hash_value, 1);
                    vector<dist_id_pair> curves_cube, curves_brute;
                    int count_nn;
                    output_file << "Query: " << query_curve.first.first << endl;
                    output_file << "Algorithm: " << algorithm << endl;
                    //--------------
                    //CUBE NEAREST NEIGHBORS
                    //FIND TIME CUBE - APPROXIMATE NEIGHBORS
                    auto start_time = std::chrono::high_resolution_clock::now();
                    curves_cube= cube_find_approximate_knn(query_curve, N, g_cube, probes, k_cube, count_nn, M_cube);
                    auto stop_time = std::chrono::high_resolution_clock::now();
                    auto time_cube = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
                    //FIND TIME BRUTE FORCE - EXACT NEIGHBORS
                    auto start_time2 = std::chrono::high_resolution_clock::now();
                    curves_brute= find_exact_knn(query_curve, N, number_of_curves);
                    auto stop_time2 = std::chrono::high_resolution_clock::now();
                    auto time_brute = std::chrono::duration_cast<std::chrono::microseconds>(stop_time2 - start_time2);
                    
                    // //PRINTING IN OUTPUT FILE
                    for (int i= 1; i <= N ; i++) {
                        //IF NO POINTS FOUND
                        if(curves_cube[i-1].id == -1){
                            output_file << "No points found in query's bucket or it's relative buckets" << endl;
                        }
                        //IF AT LEAST i POINTS HAVE BEEN FOUND
                        else if (i <= count_nn){
                            output_file << "Nearest neighbor-" << i<< ": " << curves_cube[i-1].id << endl;
                            output_file << "distanceCUBE: " << curves_cube[i-1].dist << endl;
                        }
                        //IF LESS THAN N POINTS HAVE BEEN FOUND, FOR THE REMAINING (NOT FOUND) POINTS
                        else {
                            output_file << "Nearest neighbor-" << i<< ": No more points found in query's bucket or (relative) buckets" << endl;
                            output_file << "distanceCUBE: - " << endl;
                        }

                        output_file << "distanceTrue: " << curves_brute[i-1].dist << endl;
                        
                    }
                    output_file <<  "tCUBE: " << time_cube.count() << " microseconds" << endl;
                    output_file << "tTrue: " << time_brute.count() << " microseconds" << endl;
                    output_file << endl;
                    //-----------------
                }
                else if(algorithm == "Frechet"){
                        
                    //auto tApproximateAverage = 0;
                    vector<dist_id_pair> curves_discrete_frechet, curves_brute;

                    //DISCRETE FRECHET NEAREST NEIGHBORS
                    //FIND TIME DISCRETE FRECHET - APPROXIMATE NEIGHBORS

                    auto start_time = std::chrono::high_resolution_clock::now();
                    curves_discrete_frechet= frechet_find_approximate_knn(query_curve, N, g_frechet, metric);
                    auto stop_time = std::chrono::high_resolution_clock::now();
                    auto time_discrete_frechet = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
                    //tApproximateAverage = (tApproximateAverage + time_discrete_frechet);

                    //FIND TIME BRUTE FORCE - EXACT NEIGHBORS
                    auto start_time2 = std::chrono::high_resolution_clock::now();
                    curves_brute= frechet_find_exact_knn(query_curve, N, number_of_curves, metric);
                    auto stop_time2 = std::chrono::high_resolution_clock::now();
                    auto time_brute = std::chrono::duration_cast<std::chrono::microseconds>(stop_time2 - start_time2);
                    //PRINTING IN OUTPUT FILE
                    output_file << "Query: " << query_curve.first.second << endl;
                    if(metric == "discrete"){
                        output_file << "Algorithm: LSH_Frechet_Discrete" << endl;
                    }
                    else if(metric == "continuous"){
                        output_file << "Algorithm: LSH_Frechet_Continuous" << endl;
                    }
                    
                    for (int i= 1; i <= N ; i++) {
                        if (i > curves_discrete_frechet.size()) {
                            output_file << "Nearest neighbor-" << i<< ": Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            output_file << "distanceApproximate: Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            continue;
                        }
                        output_file << "Approximate Nearest neighbor-" << i<< ": " << curves_discrete_frechet[i-1].id << endl;
                        output_file << "True Nearest neighbor: " << curves_brute[i-1].id << endl;
                        output_file << "distanceApproximate: " << curves_discrete_frechet[i-1].dist << endl;
                        output_file << "distanceTrue: " << curves_brute[i-1].dist << endl;
                    }
                    output_file <<  "tApproximate: " << time_discrete_frechet.count() << " microseconds" << endl;
                    output_file << "tTrue: " << time_brute.count() << " microseconds" << endl;
                    output_file << endl;
                    
                
                }


                hash_vector.clear();

            }

            close_file(&query_file);
            close_file(&output_file);

            continue_execution = 0;//only for now
        }
    }

    close_file(&input_file);

    return 0;
}