#include "cluster.h"

void k_means_plus_plus(int k){

    int t = 1;                              //CURRENT NUMBER OF CENTROIDS

    vector<float> curves_min_distances;     //EVERY NON CENTROID CURVE MIN DISTANCE TO CENTROIDS
    vector <float> partial_sums;

    centroids_pick_first_centroid();

    //CALCULATE DISTANCE OF ALL THE NON INPUT CURVES TO THE FIRST CENTROID
    centroids_calculate_min_distance_input(curves_min_distances);

    //CALCULATE ALL THE PARTIAL SUMS NEEDED TO CHOOSE THE NEXT CENTROID
    calculate_partial_sums(curves_min_distances, partial_sums);

    while(t < k){

        centroids_pick_next_centroid(partial_sums);
        t++;
    }

}

void reverse_assignment_lsh(G_Lsh g, fstream& output_file, int k, bool complete_flag){

    pair<pair<string, int>, vector<double>> centroid; //HERE CENTROIDS ARE THE QUERY POINTS
    vector<int> appending_curves;        //NEW CURVES THAT WILL BE ADDED TO CLUSTERS   
    vector<pair<vector<int>,int>> curves_in_range; //ALL THE CURVES THAT ARE INSIDE THE GIVEN RADIUS FROM EVERY CENTROID
    vector<vector<int>> cluster_table;   //TABLE WITH ALL THE CURVES IN EACH CLUSTER
    vector<vector<int>> hashes;          //INDEXES THAT LEAD CENTROIDS TO HASHTABLE BUCKETS
    vector<int> current_point;
    int radius = centroids_get_radii();  //SET MINIMUM DISTANCE BETWEEN CENTROIDS DIVIDED BY 2 AS FIRST RADIUS
    int last_id = curve_vector_get_size();
    bool first_iteration = true;
    int new_curves_assigned = 0;
    int previous_curves_assigned = 0;
    int num_of_curves = curve_vector_get_size();
    int complexity = curve_vector_get_curve(0).second.size(); //COMPLEXITY OF CURVES //dimensions previously

    //START COUNTING TIME
    auto start_time = std::chrono::high_resolution_clock::now();
    k_means_plus_plus(k);                //INITIALIZE K CENTROIDS USING K-MEANS++ ALGORITHM

    do{

        //GET EVERY CENTROID'S HASHTABLE BUCKET HASHES
        centroids_get_hashtable_hashes(g, hashes);
        //FOR EVERY CENTROID
        for(int i = 0; i < centroids_get_size(); i++){
            //GET CENTROID'S COORDINATES BY ACCESSING THE POINT VECTOR DATA
            centroid = curve_vector_get_curve(centroids_get_centroid(i));
            if(first_iteration){ //0 IS THE FIRST INDEX OF AN UNASSIGNED CURVE AT START - ALL CURVES UNASSIGNED
                curves_in_range.push_back(make_pair(lsh_range_search(hashes[i], radius, centroid), 0));
            }
            else{
                appending_curves = lsh_range_search(hashes[i], radius, centroid);
                curves_in_range[i].first.insert(curves_in_range[i].first.end(), appending_curves.begin(),
                                                                                appending_curves.end());
                appending_curves.clear();
            }

            //PARTITION POINTS IN RANGE OF CENTROID i TO ALREADY ASSIGNED AND NOT ASSIGNED
            partition_assigned_unassigned(curves_in_range[i]);

            //hash_vector.clear();    //WE MUST CLEAR THE HASH VECTOR AFTER EVERY ITERATION
        }

        //IF A POINT HAS BEEN ASSIGNED TO MORE THAN ONE CENTROIDS, ASSIGN IT TO THE NEAREST CENTROID
        centroids_duplicates_assign_to_nearest_centroid(curves_in_range);

        radius = radius*2;

        for(int i = 0; i < centroids_get_size(); i++){

            //GET CENTROID'S COORDINANCES BY ACCESSING THE POINT VECTOR DATA
            centroid = curve_vector_get_curve(centroids_get_centroid(i));
            //AND PERFORM RANGE SEARCH - 0 IS THE FIRST INDEX OF AN UNASSIGNED CURVE AT START
            appending_curves = lsh_range_search(hashes[i], radius, centroid);

            curves_in_range[i].first.insert(curves_in_range[i].first.end(), appending_curves.begin(),
                                                                                    appending_curves.end());
            //PARTITION CURVES IN RANGE OF CENTROID i TO ALREADY ASSIGNED AND NOT ASSIGNED
            partition_assigned_unassigned(curves_in_range[i]);


            appending_curves.clear();
            //hash_vector.clear();    //WE MUST CLEAR THE HASH VECTOR AFTER EVERY ITERATION
        }

        //IF A CURVE HAS BEEN ASSIGNED TO MORE THAN ONE CENTROIDS, ASSIGN IT TO THE NEAREST CENTROID
        centroids_duplicates_assign_to_nearest_centroid(curves_in_range);
        new_curves_assigned = is_assigned_count_assigned();

        get_cluster_table(curves_in_range, cluster_table);

        //IF THERE ARE NO ANY NEW ASSIGNMENTS TO THE CLUSTERS
        if(new_curves_assigned == previous_curves_assigned){
            break;
        }

        previous_curves_assigned = new_curves_assigned;

        first_iteration = false;

    }while(1);

    //WHEN THE CLUSTERS HAVE BEEN DEFINITIVELY FORMED STOP COUNTING TIME
    auto stop_time = std::chrono::high_resolution_clock::now();

}