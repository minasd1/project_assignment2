#include "lsh.h"

//RECEIVES THE QUERY_CURVE'S IDS AND THE CANDIDATE_CURVES' IDS
//RETURNS IF THERE IS AT LEAST ONE CURVE IN THE CANDIDATES SET WITH THE SAME ID WITH THE QUERY_CURVE
bool is_there_someone_with_same_id(vector<int>& query_curve_ids, vector<int>& candidate_curves, /*G_Lsh& g,*/
                                   int& same_id_counter, vector<int>& curves_with_same_id) {
    int i, j;
    vector<int> current_candidate_ids; //THE IDs OF A CANDIDATE CURVE
    pair<pair<string, int>, vector<double>> current_candidate;
    bool is_there_at_least_one= false;

    same_id_counter= 0;
    for (i=0; i < candidate_curves.size() ; i+=2) {
        current_candidate= curve_vector_get_curve(candidate_curves[i+1]);
        current_candidate_ids= curves_ID_vector_get_curve_all_ids(current_candidate.first.second);
        
        if (query_curve_ids[candidate_curves[i]] == current_candidate_ids[candidate_curves[i]] ) {
            curves_with_same_id.push_back(candidate_curves[i]); //IN WHICH HASHTABLE WAS THE POINT FOUND
            curves_with_same_id.push_back(candidate_curves[i+1]); //WHICH POINT IS IT
            is_there_at_least_one= true;
            same_id_counter++;
        }
    }
    return is_there_at_least_one;
}

//RECEIVES A QUERY CURVE AND RETURNS THE FIRST k NEAREST NEIGHBORS IN ASCENDING DISTANCE ORDER
//THE LAST ARGUMENT IS OPTIONAL AND SHOWS THE MAXIMUM NUMBER OF CURVES TO BE EXAMINED AS POSSIBLE NEAREST NEIGHBORS
//IF NO FORTH ARGUMENT IS GIVEN THEN ALL THE CURVES IN THE SAME BUCKETS WITH query_curve ARE EXAMINED
vector<dist_id_pair> lsh_find_approximate_knn(pair<pair<string, int>, vector<double>>& query_curve, int k, G_Lsh& g, int frechet_grid, int max_candidates)// THIRD ARG IS THE OBJECT OF CLASS G DECLEARED IN MAIN FUNCTION
{
    int curves_in_table_counter= 0; //A COUNTER OF THE ELEMENTS INSIDE THE nn_table
    vector<dist_id_pair> nn_table; //A TABLE IN WHICH THE PAIRS OF {DISTANCE, ID} OF THE NEAREST NEIGHBORING CURVES ARE STORED IN ASCENDIND DISTANCE ORDER
    vector<int> candidate_curves; //A TABLE OF ALL THE {HASHTABLE_INDEX, CURVE_ID} PAIRS, WHERE CURVE_ID BELONGS TO A CURVE THAT BELONG TO THE SAME BUCKET WITH THE query_curve IN ANY OF THE L HASHTABLES
                                  //PERHAPS ONE CURVE EXISTS MORE THAN ONCE IN THE candidate_curves TABLE BECAUSE IT APPEARS IN MORE THAN ONE HASHTABLES
    vector<int> curves_with_same_id; //A TABLE OF PAIRS WHERE THE CANDIDATES WITH THE SAME ID WITH QUERY_CURVE'S ID ARE STORED
    vector<int> buckets_indexes; //THE L BUCKETS' INDEXES WHERE THE QUERY CURVE WOULD BELONG IF IT WAS HASHED BY THE L G_HASH_FUNCTIONS
    vector<int> query_curves_ids;//A TABLE OF L QUERY_CURVE'S IDs, ONE FOR EACH HASH TABLE
    double max_distance; //THE DISTANCE BETWEEN THE QUERY CURVE AND IT'S k-th NEAREST NEIGHBORS
    double distance;
    int candidates_counter; //A COUNTER OF THE CURVES THAT HAVE BEEN COMPARED
    pair<pair<string, int>, vector<double>> current_candidate;
    vector<int> current_candidate_ids; //THE IDs OF A CANDIDATE CURVE
    dist_id_pair current_pair;
    int i;
    bool at_least_one_candidate_has_same_id;
    int same_id_counter;

    //HASH THE QUERY CURVE TO FIND THE BUCKETS WHERE IT WOULD BELONG
    g.hash(query_curve, buckets_indexes, true, frechet_grid);
    //FIND QUERY CURVE'S ID FOR EACH HASH TABLE
    g.id(query_curve, query_curves_ids, true, frechet_grid);
    //GET ALL THE CANDIDATE CURVES IDS AND THE HASHTABLE WHERE EACH OF THE WAS FOUND
    candidate_curves= hashTable_get_curves_in_buckets(buckets_indexes);
    //SEARCH IF THERE IS AT LEAST ONE CURVE IN CANDIDATES SET WITH SAME ID WITH THE QUERY CURVE
    at_least_one_candidate_has_same_id= is_there_someone_with_same_id(query_curves_ids,
    candidate_curves, same_id_counter, curves_with_same_id);
    //IF THE BOOLEAN ABOVE IS TRUE ALL THE CANDIDATE CURVES WITH DIFFRENT ID THAN THE QUERY WILL BE SKIPED
    //IF ITS FALSE THEN ALL THE CANDIDATES WILL BE EXAMINED AS POSSIBLE NEAREST NEIGHBORS

    candidates_counter= 0;
    //FIRST CHECK THE CANDIDATE_CURVES WITH THE SAME ID WITH QUERY CURVE
    for (i= 0; i < curves_with_same_id.size() ; i+=2) {
        if (candidates_counter == max_candidates) { //DONT COMPARE MORE THAN max_candidates CURVES
            break;
        }
        //GET THE CURRENT CANDIDATE CURVE'S COORDINATES AND IT's DISTANCE FROM QUERY CURVE
        current_candidate= curve_vector_get_curve(curves_with_same_id[i+1]); //+1 BECAUSE candidate_curves IS [HASHTABLE, CURVE_INDEX, HASHTABLE, ...]
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
            current_candidate= curve_vector_get_curve(candidate_curves[i+1]);//+1 BECAUSE candidate_curves IS [HASHTABLE, CURVE_INDEX, HASHTABLE, ...]
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
        if (nn_table.size() < k) { //JUST IN CASE THERE ARE LESS THAN k CANDIDATE CURVES IN THE BUCKETS
            sort(nn_table.begin(), nn_table.end(), compare_distance);
        }
    }
    return nn_table;
}

//RECEIVES A QUERY CURVE AND RETURNS THE FIRST k EXACT NEAREST NEIGHBORS
//IN ASCENDING DISTANCE ORDER USING BRUTE FORCE
vector<dist_id_pair> find_exact_knn(pair<pair<string, int>, vector<double>>& query_curve, int k, int num_of_curves)
{
    int i;
    vector<dist_id_pair> nn_table; //A TABLE IN WHICH THE PAIRS OF {DISTANCE, ID} OF THE NEAREST //NEIGHBORING CURVES ARE STORED IN ASCENDIND DISTANCE ORDER
    pair<pair<string, int>, vector<double>> current_candidate;
    dist_id_pair current_pair;
    int curves_in_table_counter= 0;
    double max_distance, distance;

    //FOR ALL THE CANDIDATE CURVES
    for (i=0; i < num_of_curves ; i++) {
        current_candidate= curve_vector_get_curve(i);
        distance= calculate_distance(query_curve.second, current_candidate.second);
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
            max_distance= nn_table[k-1].dist; //THE LAST ELEMENT OF THE nn_table HAS THE MAXIMUM DISTANCE FROM QUERY CURVE
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
    if (nn_table.size() < k) { //JUST IN CASE THERE ARE LESS THAN k CANDIDATE CURVES IN THE BUCKETS
        sort(nn_table.begin(), nn_table.end(), compare_distance);
    }

    return nn_table;
}



//IMPLEMENT APPROXIMATE RANGE SEARCH ALGORITHM
vector<int> lsh_range_search(vector<int>& g, int radius, pair<pair<string, int>, vector<double>>& query_curve){

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
                if((calculate_distance(query_curve.second, curve_vector_get_curve(hashTable_get_curve(i, g[i], j)).second, 2)) < radius){
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