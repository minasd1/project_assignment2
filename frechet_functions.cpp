#include "frechet_functions.h"

//RECEIVES A QUERY CURVE AND RETURNS THE FIRST k NEAREST NEIGHBORS IN ASCENDING DISTANCE ORDER
//THE LAST ARGUMENT IS OPTIONAL AND SHOWS THE MAXIMUM NUMBER OF CURVES TO BE EXAMINED AS POSSIBLE NEAREST NEIGHBORS
//IF NO FORTH ARGUMENT IS GIVEN THEN ALL THE CURVES IN THE SAME BUCKETS WITH query_curve ARE EXAMINED
//DISTANCE BETWEEN CURVES IS CALCULATED USING DISCRETE FRECHET METRIC
vector<dist_id_pair> frechet_find_approximate_knn(pair<pair<string, int>, vector<double>>& query_curve, int k,  G_Frechet& g, string metric, int max_candidates){
    int curves_in_table_counter= 0; //A COUNTER OF THE ELEMENTS INSIDE THE nn_table
    vector<dist_id_pair> nn_table; //A TABLE IN WHICH THE PAIRS OF {DISTANCE, ID} OF THE NEAREST NEIGHBORING CURVES ARE STORED IN ASCENDIND DISTANCE ORDER
    vector<int> candidate_curves; //A TABLE OF ALL THE {HASHTABLE_INDEX, POINT_ID} PAIRS, WHERE POINT_ID BELONGS TO A CURVE THAT BELONG TO THE SAME BUCKET WITH THE query_curve IN ANY OF THE L HASHTABLES
                                  //PERHAPS ONE CURVE EXISTS MORE THAN ONCE IN THE candidate_curves TABLE BECAUSE IT APPEARS IN MORE THAN ONE HASHTABLES
    vector<int> curves_with_same_id; //A TABLE OF PAIRS WHERE THE CANDIDATES WITH THE SAME ID WITH QUERY_CURVE'S ID ARE STORED
    vector<int> buckets_indexes; //THE L BUCKETS' INDEXES WHERE THE QUERY CURVE WOULD BELONG IF IT WAS HASHED BY THE L G_HASH_FUNCTIONS
    vector<int> query_curves_ids;//A TABLE OF L QUERY_CURVE'S IDs, ONE FOR EACH HASH TABLE
    double max_distance; //THE DISTANCE BETWEEN THE QUERY CURVE AND IT'S k-th NEAREST NEIGHBORS
    double distance;
    int candidates_counter; //A COUNTER OF THE POINTS THAT HAVE BEEN COMPARED
    pair<pair<string, int>, vector<double>> current_candidate;
    vector<int> current_candidate_ids; //THE IDs OF A CANDIDATE POINT
    dist_id_pair current_pair;
    int i;
    bool at_least_one_candidate_has_same_id;
    int same_id_counter;
    bool is_mean= false;

    Curve query_curve_continuous = convert_for_continuous_frechet(query_curve, 1);

    if(metric == "discrete"){

        //HASH THE QUERY CURVE TO FIND THE BUCKETS WHERE IT WOULD BELONG
        g.hash(query_curve, buckets_indexes, query_curves_ids, true, is_mean, 2);
    }
    else if(metric == "continuous"){

        g.hash(query_curve, buckets_indexes, query_curves_ids, true, is_mean, 1);
    }
    
    
    //GET ALL THE CANDIDATE POINTS IDS AND THE HASHTABLE WHERE EACH OF THE WAS FOUND
    candidate_curves= hashTable_get_curves_in_buckets(buckets_indexes);
   
    //SEARCH IF THERE IS AT LEAST ONE POINT IN CANDIDATES SET WITH SAME ID WITH THE QUERY POINT
    at_least_one_candidate_has_same_id= is_there_someone_with_same_id(query_curves_ids,
    candidate_curves, same_id_counter, curves_with_same_id);

    //IF THE BOOLEAN ABOVE IS TRUE ALL THE CANDIDATE POINTS WITH DIFFRENT ID THAN THE QUERY WILL BE SKIPED
    //IF ITS FALSE THEN ALL THE CANDIDATES WILL BE EXAMINED AS POSSIBLE NEAREST NEIGHBORS

    candidates_counter= 0;
    //FIRST CHECK THE CANDIDATE_POINTS WITH THE SAME ID WITH QUERY POINT
    for (i= 0; i < curves_with_same_id.size() ; i+=2) {
        if (candidates_counter == max_candidates) { //DONT COMPARE MORE THAN max_candidates POINTS
            break;
        }
        //GET THE CURRENT CANDIDATE POINT'S COORDINATES AND IT's DISTANCE FROM QUERY POINT
        current_candidate= curve_vector_get_curve(curves_with_same_id[i+1]); //+1 BECAUSE candidate_curves IS [HASHTABLE, CURVE_INDEX, HASHTABLE, ...]

        if(metric == "discrete"){
            //CALCULATE DISCRETE FRECHET DISTANCE
            distance= curve_calculate_dfd(query_curve.second, current_candidate.second);
        }
        else if(metric == "continuous"){
            //CONVERT CANDIDATE CURVE TO PROPER FORM FOR CONTINUOUS FRECHET CALCULATIONS
            Curve current_candidate_continuous = convert_for_continuous_frechet(current_candidate, 1);
            //CALCULATE CONTINUOUS FRECHET DISTANCE
            Frechet::Continuous::Distance cfd = Frechet::Continuous::distance(query_curve_continuous, current_candidate_continuous);
            distance = cfd.value;
        }
        
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
            max_distance= nn_table[k-1].dist; //THE LAST ELEMENT OF THE nn_table HAS THE MAXIMUM DISTANCE FROM QUERY POINT
            if (current_pair.dist < max_distance && !already_exist(nn_table, current_pair.id)) {
                insert_at_correct_place(nn_table, current_pair);
                max_distance= nn_table[k-1].dist;
            }
            candidates_counter++;
        }
        else { //IF THE nn_table ALREADY HAS k ID-DISTANCE PAIRS
            if (current_pair.dist < max_distance && !already_exist(nn_table, current_pair.id)) {
                insert_at_correct_place(nn_table, current_pair);
                max_distance= nn_table[k-1].dist;
            }
            candidates_counter++;
        }
    }

    //IF CANDIDATE CURVES WITH SAME ID WITH QUERY_CURVE ARE NOT ENOUGH
    if (curves_in_table_counter < k) {
        //SEARCH FOR THE REST BEST NEIGHBORS AMONG THE WHOLE CANDIDATE_CURVES DATASET

        //FOR ALL THE CANDIDATE CURVES
        for (i=0; i < candidate_curves.size() ; i+=2) { //FOR EVERY CANDIDATE CURVE
            if (candidates_counter == max_candidates) { //DONT COMPARE MORE THAN max_candidates CURVES
                break;
            }
            //GET THE CURRENT CANDIDATE CURVE'S COORDINATES AND IT'S DISTANCE FROM QUERY CURVE
            current_candidate= curve_vector_get_curve(candidate_curves[i+1]);//+1 BECAUSE candidate_points IS [HASHTABLE, POINT_INDEX, HASHTABLE, ...]
            if(metric == "discrete"){
                //CALCULATE DISCRETE FRECHET DISTANCE
                distance= curve_calculate_dfd(query_curve.second, current_candidate.second);
            }
            else if(metric == "continuous"){
                //CONVERT CANDIDATE CURVE TO PROPER FORM FOR CONTINUOUS FRECHET CALCULATIONS
                Curve current_candidate_continuous = convert_for_continuous_frechet(current_candidate, 1);
                //CALCULATE CONTINUOUS FRECHET DISTANCE
                Frechet::Continuous::Distance cfd = Frechet::Continuous::distance(query_curve_continuous, current_candidate_continuous);
                distance = cfd.value;
            }
            
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
                max_distance= nn_table[k-1].dist; //THE LAST ELEMENT OF THE nn_table HAS THE MAXIMUM DISTANCE FROM QUERY POINT
                if (current_pair.dist < max_distance && !already_exist(nn_table, current_pair.id) ) {
                    insert_at_correct_place(nn_table, current_pair);
                    max_distance= nn_table[k-1].dist;
                }
                candidates_counter++;
            }
            else { //IF THE nn_table ALREADY HAS k ID-DISTANCE PAIRS
                if (current_pair.dist < max_distance && !already_exist(nn_table, current_pair.id)) {

                    insert_at_correct_place(nn_table, current_pair);
                    max_distance= nn_table[k-1].dist;
                }
                candidates_counter++;
            }
        }
        if (nn_table.size() < k) { //JUST IN CASE THERE ARE LESS THAN k CANDIDATE POINTS IN THE BUCKETS
            sort(nn_table.begin(), nn_table.end(), compare_distance);
        }
    }
    return nn_table;


}

vector<dist_id_pair> frechet_find_exact_knn(pair<pair<string, int>, vector<double>>& query_curve, int k, int num_of_curves, string metric){
    int i;
    vector<dist_id_pair> nn_table; //A TABLE IN WHICH THE PAIRS OF {DISTANCE, ID} OF THE NEAREST //NEIGHBORING POINTS ARE STORED IN ASCENDIND DISTANCE ORDER
    pair<pair<string, int>, vector<double>> current_candidate;
    dist_id_pair current_pair;
    int curves_in_table_counter= 0;
    double max_distance, distance;

    Curve query_curve_continuous = convert_for_continuous_frechet(query_curve, 1);

    //FOR ALL THE CANDIDATE CURVES
    for (i=0; i < num_of_curves ; i++) {
        current_candidate= curve_vector_get_curve(i);
        if(metric == "discrete"){
            //CALCULATE DISCRETE FRECHET DISTANCE
            distance= curve_calculate_dfd(query_curve.second, current_candidate.second);
        }
        else if(metric == "continuous"){
            //CONVERT CANDIDATE CURVE TO PROPER FORM FOR CONTINUOUS FRECHET CALCULATIONS
            Curve current_candidate_continuous = convert_for_continuous_frechet(current_candidate, 1);
            //CALCULATE CONTINUOUS FRECHET DISTANCE
            Frechet::Continuous::Distance cfd = Frechet::Continuous::distance(query_curve_continuous, current_candidate_continuous);
            distance = cfd.value;
        }
        
        //CREATE A PAIR WITH THESE TWO VALUES
        current_pair.dist= distance;
        current_pair.id= current_candidate.first.second;

        if (curves_in_table_counter < k) { //THE FIRST k CANDIDATES TO COME ARE DIRECTLY PUSHED INTO THE NEAREST NEIGHBORS VECTOR
            if (!already_exist(nn_table, current_pair.id)) {
                nn_table.push_back(current_pair);
                curves_in_table_counter++;
            }
        }
        else if (curves_in_table_counter == k) { //WHEN THE nn_table HAS BEEN FILLED FULLY FOR THE FIRST TIME
            sort(nn_table.begin(), nn_table.end(), compare_distance); //SORT THE FIRST k PAIRS IN ASCENDING DISTANCE ORDER
            max_distance= nn_table[k-1].dist; //THE LAST ELEMENT OF THE nn_table HAS THE MAXIMUM DISTANCE FROM QUERY POINT
            if (current_pair.dist < max_distance) {
                insert_at_correct_place(nn_table, current_pair);
                max_distance= nn_table[k-1].dist;
            }
        }
        else { //IF THE nn_table ALREADY HAS k ID-DISTANCE PAIRS
            if (current_pair.dist < max_distance) {
                insert_at_correct_place(nn_table, current_pair);
                max_distance= nn_table[k-1].dist;
            }
        }
    }
    if (nn_table.size() < k) { //JUST IN CASE THERE ARE LESS THAN k CANDIDATE POINTS IN THE BUCKETS
        sort(nn_table.begin(), nn_table.end(), compare_distance);
    }

    return nn_table;
}

vector<int> frechet_range_search(vector<int>& g, int radius, pair<pair<string, int>, vector<double>>& query_curve){
    
    int retrieved_items = 0;
    int count = 0;
    int max_retrieved_items = 20* hashTable_get_num_of_htables();

    vector<int> curves_in_range;

    //ACCESS TO I-TH HASHTABLE
    for(int i = 0; i < hashTable_get_num_of_htables(); i++){
        
        //ACCESS TO THE SPECIFIC BUCKET THAT I-TH G FUNCTION INDICATES
        for(int j = 0; j < hashTable_get_bucket_size(i, g[i]); j++){
            //IF INPUT POINT IS NOT THE SAME AS THE QUERY POINT AND DOES NOT ALREADY EXIST IN POINTS IN RANGE
            if((query_curve.first.second != hashTable_get_curve(i, g[i], j))
        && (!already_exists(curves_in_range, hashTable_get_curve(i, g[i], j))) 
        && (!already_assigned(hashTable_get_curve(i, g[i], j)))){   //ALSO CHECK THAT IT IS NOT ASSIGNED ALREADY - CLUSTERING
                //IF DISTANCE OF J-TH POINT IN THIS BUCKET IS IN THE GIVEN RADIUS
                if((curve_calculate_dfd(query_curve.second, curve_vector_get_curve(hashTable_get_curve(i, g[i], j)).second)) < radius){
                    //THEN ADD IT'S ID TO POINTS_IN_RANGE VECTOR
                    curves_in_range.push_back(hashTable_get_curve(i, g[i], j));
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