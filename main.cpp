#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <cmath>
#include <chrono>
#include <random>

#include "vector_ops.h"
#include "file_functions.h"
#include "hashTable.h"
#include "hash_functions.h"


int main(int argc, char* argv[]){

    int k = 4;                         //NUMBER OF H FUNCTIONS USED IN FUNCTION G
    int L = 5;                      //NUMBER OF HASHTABLES(LSH)
    int N ;                         //NUMBER OF NEAREST MEIGHBORS
    float R ;                       //SEARCH RADIUS
    int probes ;                    //MAX NUMBER OF HYPERCUBE VERTICES TO BE CHECKED
    int k_cube ;                    //D'
    int M_cube ;                    //MAX NUMBER OF CANDIDATE POINTS TO BE CHECKED
    int window= 100;                
    int k_cluster ;                 //NUMBER OF CENTROIDS - CLUSTERING
    double delta = 2.5;
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
    int curves_divider = 16;        //USED TO GET TOTAL CURVES IN EACH HASHTABLE
    int continue_execution = 1;
    string name;
    int id;
    vector<double> values;
    vector<int> hash_vector;
    int M = pow(2, 31) - 5;
    int count = 0;
    string algorithm = "Frechet", metric = "discrete";

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

    curves_ID_vector_initialize(number_of_curves, L);

    //NUMBER OF BUCKETS IN EACH HASHTABLE
    buckets = number_of_curves/curves_divider;

    //INITIALIZE G FUNCTION THAT LEADS US TO LSH HASHTABLE BUCKETS
    G_Lsh g_lsh(k, num_of_curve_values, generator, window, M, buckets, L);

    if((strcmp(argv[0], "./search") == 0) && ((algorithm == "LSH") || ((algorithm == "Frechet") && (metric == "discrete")))){
        //INITIALIZE L HASHTABLES WITH HASHTABLESIZE BUCKETS AND ZERO POINTS IN EACH BUCKET
        hashTable_initialization(L, buckets);
        
        if(algorithm == "Frechet"){
            //INITIALIZE G FRECHET FUNCTION THAT USES g_lsh
            G_Frechet g_frechet(g_lsh, generator, L, delta, num_of_curve_values);

            //INSERT ALL INPUT CURVES TO THE HASHTABLES
            for(int i = 0; i < number_of_curves; i++){

                g_lsh.hash(curve_vector_get_curve(i), hash_vector, 0, 0);
            }
        }
    }

    
    

    return 0;
}