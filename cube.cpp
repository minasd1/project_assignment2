#include <algorithm>
#include "cube.h"
#include "knn_table_functions.h"

//RETURN ALL POINTS THAT ARE INSIDE THE GIVEN RANGE OF QUERY
vector<int> cube_range_search(int g, int radius, int probes, int dimensions, pair<pair<string, int>, vector<double>>& query_curve){

    int retrieved_items = 0;
    int max_retrieved_items = 100;

    //A VECTOR WITH ALL THE RELATIVE BUCKETS OF QUERY CURVE IN HYPERCUBE
    vector <unsigned int> relative_buckets_indexes;
    relative_buckets_indexes = get_relative_buckets(g, probes, dimensions);

    //CURVES THAT ARE IN THE GIVEN RADIUS
    vector<int> curves_in_range;

    //FOR ALL THE ELEMENTS IN QUERY CURVE'S BUCKET (INDICATED BY G)
    for(int i = 0; i < hyperCube_get_bucket_size(g); i++){
        //IF INPUT CURVE IS NOT THE SAME AS THE QUERY CURVE AND DOES NOT ALREADY EXIST IN CURVES IN RANGE
        if((query_curve.first.second != hyperCube_get_curve(g, i))
        && (!already_exists(curves_in_range, hyperCube_get_curve(g, i)))
        && (!already_assigned(hyperCube_get_curve(g, i)))){
            if(calculate_distance(query_curve.second, curve_vector_get_curve(hyperCube_get_curve(g, i)).second, 2) < radius){
                //ADD THE POINT TO POINTS IN RANGE
                curves_in_range.push_back(hyperCube_get_curve(g, i));
                retrieved_items++;

                if(retrieved_items == max_retrieved_items){

                    return curves_in_range;
                }
            }
        }

    }

    //FOR ALL THE RELATIVE BUCKETS OF QUERY CURVE
    for(int i = 0; i < relative_buckets_indexes.size(); i++){
        //FOR ALL THEIR ELEMENTS
        for(int j = 0; j < hyperCube_get_bucket_size(relative_buckets_indexes[i]); j++){
            //IF INPUT POINT IS NOT THE SAME AS THE QUERY POINT AND DOES NOT ALREADY EXIST IN POINTS IN RANGE
            if((query_curve.first.second != hyperCube_get_curve(relative_buckets_indexes[i], j))
            && (!already_exists(curves_in_range, hyperCube_get_curve(relative_buckets_indexes[i], j)))
            && (!already_assigned(hyperCube_get_curve(relative_buckets_indexes[i], j)))){
                //IF DISTANCE BETWEEN QUERY AND A CURVE IN RELATIVE BUCKET IS IN THE RADIUS
                if(calculate_distance(query_curve.second, curve_vector_get_curve(hyperCube_get_curve(relative_buckets_indexes[i], j)).second, 2) < radius){
                    //ADD IT TO POINTS IN RANGE
                    curves_in_range.push_back(hyperCube_get_curve(relative_buckets_indexes[i], j));
                    retrieved_items++;

                    if(retrieved_items == max_retrieved_items){

                        return curves_in_range;
                    }
                }
            }
            
        }
    }

    return curves_in_range;
} 

//RECEIVES A QUERY CURVE AND RETURNS THE FIRST k NEAREST NEIGHBORS IN ASCENDING DISTANCE ORDER
//THE LAST ARGUMENT IS OPTIONAL AND SHOWS THE MAXIMUM NUMBER OF CURVES TO BE EXAMINED AS POSSIBLE NEAREST NEIGHBORS
//IF NO FORTH ARGUMENT IS GIVEN THEN ALL THE CURVES IN THE SAME BUCKETS WITH query_curve ARE EXAMINED
vector<dist_id_pair> cube_find_approximate_knn(pair<pair<string, int>, vector<double>>& query_curve, int k,  G_Hypercube& g, int probes, int dimensions, int& count_nn, int max_candidates)
{
    unsigned int query_hash;
    int curves_in_table_counter= 0;//A COUNTER OF THE ELEMENTS INSIDE THE nn_table
    vector<dist_id_pair> nn_table; //A TABLE IN WHICH THE PAIRS OF {DISTANCE, ID} OF THE NEAREST NEIGHBORING CURVES ARE STORED IN ASCENDIND DISTANCE ORDER
    vector <unsigned int> relative_buckets_indexes; //A VECTOR WITH ALL THE RELATIVE BUCKETS OF QUERY CURVE IN HYPERCUBE
    vector<int> candidate_curves;  //A TABLE OF INPUT CURVES TO BE CHECKED AS POSSIBLE NEAREST NEIGHBORS
    double max_distance;            //THE DISTANCE BETWEEN THE QUERY CURVE AND IT'S k-th NEAREST NEIGHBORS
    double distance;
    int candidates_counter = 0;        //A COUNTER OF THE CURVES THAT HAVE BEEN COMPARED
    pair<pair<string, int>, vector<double>> current_candidate;
    dist_id_pair current_pair;
    int i, j;

    //GET THE BUCKET OF THE HASH TABLE IN WHICH THE QUERY CURVE WOULD BELONG
    //IN OTHER WORDS GET QUERY CURVE'S HASH VALUE
    g.hash(query_curve, query_hash, true);
    //GET THE CURVES THAT BELONG IN THE SAME BUCKET WITH THE QUERY
    candidate_curves= hyperCube_get_curves_in_bucket(query_hash);

    //FOR ALL THE ELEMENTS IN QUERY CURVE'S BUCKET (INDICATED BY query_hash)
    for (i= 0; i < hyperCube_get_bucket_size(query_hash); i++) {
        if (candidates_counter == max_candidates) { //DONT COMPARE MORE THAN max_candidates CURVES
            break;
        }
        //GET THE CURRENT CANDIDATE CURVE'S COORDINATES AND IT's DISTANCE FROM QUERY CURVE
        current_candidate= curve_vector_get_curve(hyperCube_get_curve(query_hash, i)); 
        distance= calculate_distance(query_curve.second, current_candidate.second);
        //CREATE A PAIR WITH THESE TWO VALUES
        current_pair.dist= distance;
        current_pair.id= current_candidate.first.second;

        if (curves_in_table_counter < k) { //THE FIRST k CANDIDATES TO COME ARE DIRECTLY PUSHED INTO THE NEAREST NEIGHBORS VECTOR
            if (!already_exist(nn_table, current_pair.id)) {
                nn_table.push_back(current_pair);
                curves_in_table_counter++;
            }
            candidates_counter++;
        }
        else if (curves_in_table_counter == k) { //WHEN THE nn_table HAS BEEN FILLED FULLY FOR THE FIRST TIME
            sort(nn_table.begin(), nn_table.end(), compare_distance); //SORT THE FIRST k PAIRS IN ASCENDING DISTANCE ORDER
            max_distance= nn_table[k-1].dist; //THE LAST ELEMENT OF THE nn_table HAS THE MAXIMUM DISTANCE FROM QUERY CURVE
            if (current_pair.dist < max_distance) {
                insert_at_correct_place(nn_table, current_pair);
                max_distance= nn_table[k-1].dist; 
            }
            candidates_counter++;
        }
        else { //IF THE nn_table ALREADY HAS k ID-DISTANCE PAIRS
            if (current_pair.dist < max_distance) {
                insert_at_correct_place(nn_table, current_pair);
                max_distance= nn_table[k-1].dist; 
            }
            candidates_counter++;
        }
        count_nn= nn_table.size(); //STORE THE NUMBER OF NEIGHBORS THAT ARE CONTAINED IN THE NN_TABLE
    }

    //IF THE CURVES IN THE QUERY'S BUCKET ARE LESS THAN K
    //FILL THE REMAINING NN_TABLE CELLS WITH CURVES OF THE RELATIVE BUCKETS
    if (curves_in_table_counter < k) {
        relative_buckets_indexes = get_relative_buckets(query_hash, probes, dimensions);
        //FOR ALL THE RELATIVE BUCKETS OF QUERY CURVE
        for(i = 0; i < relative_buckets_indexes.size(); i++){
            if (candidates_counter == max_candidates) { //DONT COMPARE MORE THAN max_candidates CURVES
                break;
            }
            //GET THE CURVESS OF THE CURRENT BUCKET
            candidate_curves= hyperCube_get_curves_in_bucket(relative_buckets_indexes[i]);

            //FOR EVERY CURVE IN THIS BUCKET
            for (j= 0; j < hyperCube_get_bucket_size(relative_buckets_indexes[i]); j++) {
                if (candidates_counter == max_candidates) { //DONT COMPARE MORE THAN max_candidates CURVES
                    break;
                }
                //GET THE CURRENT CANDIDATE CURVE'S COORDINATES AND IT's DISTANCE FROM QUERY CURVE
                current_candidate= curve_vector_get_curve(hyperCube_get_curve(relative_buckets_indexes[i], j)); 
                distance= calculate_distance(query_curve.second, current_candidate.second);
                //CREATE A PAIR WITH THESE TWO VALUES
                current_pair.dist= distance;
                current_pair.id= current_candidate.first.second;

                if (curves_in_table_counter < k) { //THE FIRST k CANDIDATES TO COME ARE DIRECTLY PUSHED INTO THE NEAREST NEIGHBORS VECTOR
                    if (!already_exist(nn_table, current_pair.id)) {
                        nn_table.push_back(current_pair);
                        curves_in_table_counter++;
                    }
                    candidates_counter++;
                }
                else if (curves_in_table_counter == k) { //WHEN THE nn_table HAS BEEN FILLED FULLY FOR THE FIRST TIME
                    sort(nn_table.begin(), nn_table.end(), compare_distance); //SORT THE FIRST k PAIRS IN ASCENDING DISTANCE ORDER
                    max_distance= nn_table[k-1].dist; //THE LAST ELEMENT OF THE nn_table HAS THE MAXIMUM DISTANCE FROM QUERY CURVE
                    if (current_pair.dist < max_distance) {
                        insert_at_correct_place(nn_table, current_pair);
                        max_distance= nn_table[k-1].dist;
                    }
                    candidates_counter++;
                }
                else { //IF THE nn_table ALREADY HAS k ID-DISTANCE PAIRS
                    if (current_pair.dist < max_distance) {
                        insert_at_correct_place(nn_table, current_pair);
                    }
                    candidates_counter++;
                }
            }
        }
        count_nn= nn_table.size(); //STORE THE NUMBER OF NEIGHBORS THAT ARE CONTAINED IN THE NN_TABLE
    }

    if(nn_table.size() == 0){
    //IF THERE ARE NO CURVES IN QUERY'S BUCKET OR QUERY'S ADJACENT BUCKETS
        count_nn= nn_table.size(); //STORE THE NUMBER OF NEIGHBORS THAT ARE CONTAINED IN THE NN_TABLE
        current_pair.id= -1;
        current_pair.dist = -1;

        nn_table.push_back(current_pair);
    }

    //JUST IN CASE THERE ARE LESS THAN k CANDIDATE CURVES IN ALL THE BUCKETS
    //(BOTH QUERY'S AND RELATIVE BUCKETS)
    if (nn_table.size() < k) {
        count_nn= nn_table.size(); //STORE THE NUMBER OF NEIGHBORS THAT ARE CONTAINED IN THE NN_TABLE
        sort(nn_table.begin(), nn_table.end(), compare_distance);
    }

    return nn_table;
}
