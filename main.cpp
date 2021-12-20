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
#include "cmd_line_args.h"
#include "user_input.h"
#include "cluster.h"
#include "conf_file.h"
#include "mean_curve.h"



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
    double delta = 0.7;             //experiments done
    double epsilon = 0.1;

    fstream input_file;             //FILE WE READ INPUT FROM
    fstream query_file;             //FILE WE READ QUERIES FROM
    fstream config_file;            //CLUSTER CONFIGURATION FILE
    fstream output_file;            //FILE TO WRITE OUTPUT TO
    
    string line, token, algorithm, metric, name, assignment, update;
    string input_file_name, query_file_name, output_file_name, config_file_name;
   
    int start, finish, buckets, id, count, M, num_of_curve_values, error,
        query_num_of_curve_values, continue_execution, sum_approximate_time, 
        sum_exact_time, query_curves_counted, factor, max_mean_curve_length; 
    
    int curves_divider = 16;        //USED TO GET TOTAL CURVES IN EACH HASHTABLE
    
    bool is_query_curve, is_mean_curve, first_iteration, complete_flag, new_query_file, 
         flag_frechet_cluster, datapath_given, query_given, output_given, algorithm_given,
         metric_given, silhouette_flag;

    unsigned int hash_value, number_of_curves ;
 
    vector<double> values;
    vector<int> hash_vector, id_vector;
    pair<pair<string, int>, vector<double>> current_curve;

    double max_coordinate_value, max_value, query_max_value, query_max_coordinate_value;
    double maf; //MAXIMUM APPROXIMATION FACTOR

    finish = 0;
    count = 0;
    M = pow(2, 31) - 5;
    num_of_curve_values = 0;
    query_num_of_curve_values = 0;
    number_of_curves = 0;
    continue_execution = 1;
    factor = 1;
    max_value = 0.0;
    query_max_value = 0.0;
    is_query_curve = true;
    first_iteration = true;
    new_query_file = false;
    flag_frechet_cluster = false;
    is_mean_curve = true;

    error= read_cmd_args(argc, argv, input_file_name, query_file_name, k, k_cube, L, output_file_name, 
                  N, R, M_cube, probes, config_file_name, assignment, complete_flag, algorithm, 
                  metric, delta, update, datapath_given, query_given, output_given, algorithm_given,
                  metric_given, silhouette_flag);

    if (error) {
        return -1;
    }
                
    read_path(input_file_name, query_file_name, output_file_name, algorithm, metric, 
                   datapath_given, query_given, output_given, algorithm_given, metric_given);

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
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

    //MARK ALL THE INPUT POINTS AS UNASSIGNED
    is_assigned_initialize();

    if(strcmp(argv[0], "./cluster") == 0){

        read_configuration_file(config_file, config_file_name, k_cluster, L, k, M_cube, k_cube, probes);
        if((assignment == "Classic" && update == "Mean_Frechet") || (assignment == "LSH_Frechet")){
            //V VECTORS IN LSH MUBT HAVE THE SAME SIZE AS THE CONCATENATED CURVE VECTOR
            factor = 8;//experiments done

            flag_frechet_cluster = true;
            max_mean_curve_length= num_of_curve_values; //experiments done
        }
    }

    max_coordinate_value = max(max_value, (double)num_of_curve_values);

    if(strcmp(argv[0], "./search") == 0){
   
        query_max_coordinate_value = file_get_max_value(query_file, query_file_name); 
        max_coordinate_value = max(max_coordinate_value, query_max_coordinate_value);

        if(algorithm == "Frechet" && metric == "discrete"){
            factor = 2;
        }
    }
    
    curves_ID_vector_initialize(number_of_curves, L);

    //NUMBER OF BUCKETS IN EACH HASHTABLE
    buckets = number_of_curves/curves_divider;

    //INITIALIZE G FUNCTION THAT LEADS US TO LSH HASHTABLE BUCKETS
    G_Lsh g_lsh(k, num_of_curve_values*factor, generator, window, M, buckets, L);

    //INITIALIZE G FUNCTION THAT LEADS US TO HYPERCUBE BUCKETS
    G_Hypercube g_cube(num_of_curve_values, generator, window, k_cube);

    //INITIALIZE G FRECHET FUNCTION THAT USES g_lsh
    G_Frechet g_frechet(g_lsh, generator, L, delta, max_coordinate_value, num_of_curve_values, flag_frechet_cluster);

    if(strcmp(argv[0], "./search") == 0 || strcmp(argv[0], "./cluster") == 0){

        //INITIALIZE L HASHTABLES WITH HASHTABLESIZE BUCKETS AND ZERO POINTS IN EACH BUCKET
        hashTable_initialization(L, buckets);

        if((algorithm == "Frechet" || assignment == "LSH_Frechet")){

            //INSERT ALL INPUT CURVES TO THE HASHTABLES
            for(int i = 0; i < number_of_curves; i++){
                if(metric == "discrete" || assignment == "LSH_Frechet"){

                    g_frechet.hash(curve_vector_get_curve(i), hash_vector, id_vector, !is_query_curve, !is_mean_curve, 2);
                }
                else if(metric == "continuous"){

                    g_frechet.hash(curve_vector_get_curve(i), hash_vector, id_vector,  !is_query_curve, !is_mean_curve, 1);
                }    
                
            }
            
            
            hash_vector.clear();
        }
        else if (algorithm == "LSH" || assignment == "LSH") {
            //INSERT ALL INPUT CURVES TO THE HASHTABLES
            for(int i = 0; i < number_of_curves; i++){

                g_lsh.hash(curve_vector_get_curve(i), hash_vector, 0, 0);
            }
            hash_vector.clear();
        }
        else if (algorithm == "Hypercube" || assignment == "Hypercube") {
            //INITIALIZE A HYPERCUBE WITH 2^D' BUCKETS AND ZERO POINTS IN EACH BUCKET
            hyperCube_initialization(pow(2, k_cube));

            //INSERT ALL THE POINTS FROM INPUT FILE TO HYPERCUBE BUCKETS
            for(int i = 0; i < number_of_curves; i++){
                g_cube.hash(curve_vector_get_curve(i), hash_value, 0);
            }
        }

    }

    

    while(continue_execution == 1){

        if(strcmp(argv[0], "./search") == 0){

            if(new_query_file == true){

                //OPEN FILE TO READ QUERY FILES FROM
                open_file(&query_file, query_file_name, fstream::in);
            }

            //OPEN FILE TO WRITE RESULTS TO
            open_file(&output_file, output_file_name, fstream::out);
            finish = 0;

            maf= 0.0;
            sum_approximate_time= 0;
            sum_exact_time= 0;
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
                    output_file << "Algorithm: LSH_Vector" << endl;
                    for (int i= 1; i <= N ; i++) {
                        if (i > curves_lsh.size()) {
                            output_file << "Approximate Nearest neighbor-" << i<< ": Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            output_file << "distanceLSH: Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            continue;
                        }
                        output_file << "Approximate Nearest neighbor-" << i<< ": ";
                        current_curve= curve_vector_get_curve(curves_lsh[i-1].id);
                        output_file << current_curve.first.first << endl;
                        output_file << "True Nearest neighbor-" << i<< ": ";
                        current_curve= curve_vector_get_curve(curves_brute[i-1].id);
                        output_file << current_curve.first.first << endl;
                        output_file << "distanceApproximate: " << curves_lsh[i-1].dist << endl;
                        output_file << "distanceTrue: " << curves_brute[i-1].dist << endl;
                        output_file << endl;
                    }
                    if ((curves_lsh[0].dist / (double) curves_brute[0].dist) > maf) {
                        maf= curves_lsh[0].dist / (double) curves_brute[0].dist;
                    }
                    query_curves_counted++;
                    sum_approximate_time+= time_lsh.count();
                    sum_exact_time+= time_brute.count();
                }
                else if(algorithm == "Hypercube"){
                    g_cube.hash(query_curve, hash_value, 1);
                    vector<dist_id_pair> curves_cube, curves_brute;
                    int count_nn;
                    output_file << "Query: " << query_curve.first.first << endl;
                    output_file << "Algorithm: " << "Hypercube" << endl;
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
                            output_file << "Approximate Neighrest neighbor - " << i << ": ";
                            output_file << "No points found in query's bucket or it's relative buckets" << endl;
                        }
                        //IF AT LEAST i POINTS HAVE BEEN FOUND
                        else if (i <= count_nn){
                            output_file << "Approximate Nearest Neighbor - " << i<< ": ";
                            current_curve= curve_vector_get_curve(curves_cube[i-1].id);
                            output_file << current_curve.first.first << endl;
                            output_file << "True Nearest Neighbor - " << i<< ": ";
                            current_curve= curve_vector_get_curve(curves_brute[i-1].id);
                            output_file << current_curve.first.first << endl;
                            output_file << "distanceApproximate: " << curves_cube[i-1].dist << endl;
                            }
                        //IF LESS THAN N POINTS HAVE BEEN FOUND, FOR THE REMAINING (NOT FOUND) POINTS
                        else {
                            output_file << "Nearest neighbor - " << i<< ": No more points found in query's bucket or (relative) buckets" << endl;
                            output_file << "distanceCUBE: - " << endl;
                        }

                        output_file << "distanceTrue: " << curves_brute[i-1].dist << endl;
                        output_file << endl;
                        
                    }
                    if ((curves_cube[0].dist / curves_brute[0].dist) > maf) {
                        maf= curves_cube[0].dist /(double) curves_brute[0].dist;
                    }
                    query_curves_counted++;
                    sum_approximate_time+= time_cube.count();
                    sum_exact_time+= time_brute.count();
                    
                    //-----------------
                }
                else if(algorithm == "Frechet"){
                        
                    //auto tApproximateAverage = 0;
                    vector<dist_id_pair> curves_frechet, curves_brute;

                    //DISCRETE FRECHET NEAREST NEIGHBORS
                    //FIND TIME DISCRETE FRECHET - APPROXIMATE NEIGHBORS

                    auto start_time = std::chrono::high_resolution_clock::now();
                    curves_frechet= frechet_find_approximate_knn(query_curve, N, g_frechet, metric);
                    auto stop_time = std::chrono::high_resolution_clock::now();
                    auto time_frechet = std::chrono::duration_cast<std::chrono::microseconds>(stop_time - start_time);
                    //tApproximateAverage = (tApproximateAverage + time_discrete_frechet);

                    //FIND TIME BRUTE FORCE - EXACT NEIGHBORS
                    auto start_time2 = std::chrono::high_resolution_clock::now();
                    curves_brute= frechet_find_exact_knn(query_curve, N, number_of_curves, metric);
                    auto stop_time2 = std::chrono::high_resolution_clock::now();
                    auto time_brute = std::chrono::duration_cast<std::chrono::microseconds>(stop_time2 - start_time2);
                    //PRINTING IN OUTPUT FILE
                    output_file << "Query: " << query_curve.first.first << endl;
                    if(metric == "discrete"){
                        output_file << "Algorithm: LSH_Frechet_Discrete" << endl;
                    }
                    else if(metric == "continuous"){
                        output_file << "Algorithm: LSH_Frechet_Continuous" << endl;
                    }
                    
                    for (int i= 1; i <= N ; i++) {
                        if (i > curves_frechet.size()) {
                            output_file << "Nearest neighbor-" << i<< ": Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            output_file << "distanceApproximate: Not enough points in buckets (Consider to decrease hash table size or window)" << endl;
                            continue;
                        }
                        output_file << "Approximate Nearest neighbor-" << i<< ": ";
                        current_curve= curve_vector_get_curve(curves_frechet[i-1].id);
                        output_file << current_curve.first.first << endl;
                        output_file << "True Nearest neighbor-" << i<< ": ";
                        current_curve= curve_vector_get_curve(curves_brute[i-1].id);
                        output_file << current_curve.first.first << endl;
                        output_file << "distanceApproximate: " << curves_frechet[i-1].dist << endl;
                        output_file << "distanceTrue: " << curves_brute[i-1].dist << endl;
                        output_file << endl;
                    }
                    if ((curves_frechet[0].dist / (double) curves_brute[0].dist) > maf) {
                        maf= curves_frechet[0].dist / (double) curves_brute[0].dist;
                    }
                    query_curves_counted++;
                    sum_approximate_time+= time_frechet.count();
                    sum_exact_time+= time_brute.count();
                    
                
                }


                hash_vector.clear();

            }
            output_file << "tApproximateAverage: " << sum_approximate_time / (double) query_curves_counted << endl;
            output_file << "tTrueAverage: " << sum_exact_time / (double) query_curves_counted << endl;
            output_file << "MAF: " << maf << endl;

            read_user_input(query_file_name, &continue_execution);
            new_query_file = true;
            close_file(&query_file);
            close_file(&output_file);
        }
        else if(strcmp(argv[0], "./cluster") == 0){
            
            //OPEN FILE TO WRITE RESULTS TO
            open_file(&output_file, output_file_name, fstream::out);
            if (assignment == "Classic") {

                lloyds(k_cluster, output_file, assignment, update, 
                epsilon, max_mean_curve_length, generator, silhouette_flag, complete_flag);
            }
            else if(assignment == "LSH"){
                reverse_assignment_lsh(g_lsh, output_file, k_cluster, assignment, update, silhouette_flag, complete_flag);
            }
            else if(assignment == "Hypercube"){

                reverse_assignment_cube(g_cube, output_file, k_cluster, probes, assignment, update, silhouette_flag, complete_flag);
            }
            else if(assignment == "LSH_Frechet"){

                reverse_assignment_frechet(g_frechet, output_file, k_cluster, assignment, epsilon, 
                                            max_mean_curve_length, generator, silhouette_flag, complete_flag);
            }
            

            continue_execution = 0;
            close_file(&output_file);
        }
    }

    close_file(&input_file);
    

    return 0;
}