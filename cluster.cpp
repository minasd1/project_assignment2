#include <fstream>
#include <limits>
#include <chrono>
#include "cluster.h"
#include "lloyds_auxiliary.h"
#include "silhouette.h"
#include "vector_ops.h"

void k_means_plus_plus(int k, string assignment){

    int t = 1;                              //CURRENT NUMBER OF CENTROIDS

    vector<float> curves_min_distances;     //EVERY NON CENTROID CURVE MIN DISTANCE TO CENTROIDS
    vector <float> partial_sums;

    centroids_pick_first_centroid();

    //CALCULATE DISTANCE OF ALL THE NON INPUT CURVES TO THE FIRST CENTROID
    centroids_calculate_min_distance_input(curves_min_distances, assignment);

    //CALCULATE ALL THE PARTIAL SUMS NEEDED TO CHOOSE THE NEXT CENTROID
    calculate_partial_sums(curves_min_distances, partial_sums);

    while(t < k){

        centroids_pick_next_centroid(partial_sums);
        t++;
    }

}

//IMPLEMENTATION OF THE LLOYDS ALGORITHM
void lloyds(int number_of_clusters, fstream& output_file, string assignment, bool complete_flag)
{
    int i, dimensions, j;
    int nearest_centroid; //NEAREST CENTROID'S INDEX IN THE centroid TABLE
    int changes_made= 0; //HOW MANY CURVES CHANGED CLUSTER IN A NEW ASSIGNMENT
    pair<pair<string, int>, vector<double>> current_curve;
    vector<vector<int>> new_cluster_table, previous_cluster_table;
    float change_rate; //THE RATE OF CURVES THAT CHANGED CLUSTER TO THE TOTAL NUMBER OF CURVES
    int num_of_curves = curve_vector_get_size();
    int last_known_id= num_of_curves-1;


    dimensions= curve_vector_get_curve(1).second.size();
    previous_cluster_table.resize(number_of_clusters);
    //START COUNTING TIME
    auto start_time = std::chrono::high_resolution_clock::now();
    //INITIALIZE THE CENDROID CURVES
    k_means_plus_plus(number_of_clusters, assignment);
    //ASSIGN CURVESS IN CLUSTERS FOR THE FIRST TIME
    for (i=0 ; i < num_of_curves ; i++) { //FOR EVERY POINT
        current_curve= curve_vector_get_curve(i);
        nearest_centroid= find_nearest_centroid(current_curve);
        previous_cluster_table[nearest_centroid].push_back(current_curve.first.second);
        changes_made++;
    }
    change_rate= float(changes_made)/float(num_of_curves); //INITIALLY change_rate WILL BE 1 (100%)

    //LOOP UNTIL A SMALL PERCENTAGE OF CURVES CHANGE CLUSTER
    //OR THE MAXIMUM NUMBER OF ITERATIONS HAS BEEN REACHED
    while (change_rate > 0.01) {
        //UPDATE THE CENTROIDS
        update_as_vector(previous_cluster_table, last_known_id);
        //PREPARE THE NEW CLUSTER TABLE FOR THE NEW CENTROIDS ASSIGNMENT
        changes_made= 0;
        new_cluster_table.clear();
        new_cluster_table.resize(number_of_clusters);

        //MAKE A NEW ASSIGNMENT FOR ALL THE CURVES
        for (i=0 ; i < num_of_curves ; i++) { //FOR EVERY CURVE
            current_curve= curve_vector_get_curve(i);
            nearest_centroid= find_nearest_centroid(current_curve);
            
            //IF A CURVE IS BEING ASSIGNED IN A DIFFERENT CLUSTER THAN THE ONE IT WAS ASSIGNED IN THE PREVIOUS ASSIGNMENT
            if (!already_in_that_cluster(previous_cluster_table, nearest_centroid, current_curve.first.second)) {
                changes_made++;
            }
            new_cluster_table[nearest_centroid].push_back(current_curve.first.second);
        }
        previous_cluster_table= new_cluster_table;
        change_rate= float(changes_made)/float(num_of_curves);
        
    }
    //WHEN THE CLUSTERS HAVE BEEN DEFINITIVELY FORMED STOP COUNTING TIME
    auto stop_time = std::chrono::high_resolution_clock::now();

    //PRINT THE RESULTS
    output_file << "Algorithm: Assignment Classic Update Mean Frechet/Vector" << endl;
    for (i= 0 ; i < number_of_clusters ; i++) {
        output_file << "CLUSTER-" << i+1 << " {size: " << previous_cluster_table[i].size();
        output_file << " centroid: ";
        current_curve= curve_vector_get_curve(centroids_get_centroid(i));
        for (j=0 ; j < dimensions ; j++) {
            output_file <<  current_curve.second[j]<< " ";
        }
        output_file << "}" << endl;
    }
    auto time_passed = std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time);
    output_file << "clustering_time: " << time_passed.count() << " seconds" << endl;
    output_file << "Silhouette: ";
    print_silhouette(previous_cluster_table, output_file);
    if (complete_flag) {
        for(i=0 ; i < number_of_clusters ; i++) {
            output_file << "CLUSTER-" << i+1 << " {size: " << previous_cluster_table[i].size();
            output_file << " centroid: [";
            current_curve= curve_vector_get_curve(centroids_get_centroid(i));
            for (j=0 ; j < dimensions ; j++) {
                output_file << current_curve.second[j] << " ";
            }
            output_file << "]";
            output_file << ", ";
            for (j= 0; j < previous_cluster_table[i].size(); j++) {
                output_file << previous_cluster_table[i][j] << " ";
            }
            output_file << "}" << endl;
        }
    }
    
    
}

void reverse_assignment_lsh(G_Lsh g, fstream& output_file, int k, string assignment, bool complete_flag){

    pair<pair<string, int>, vector<double>> centroid; //HERE CENTROIDS ARE THE QUERY POINTS
    vector<int> appending_curves;        //NEW CURVES THAT WILL BE ADDED TO CLUSTERS   
    vector<pair<vector<int>,int>> curves_in_range; //ALL THE CURVES THAT ARE INSIDE THE GIVEN RADIUS FROM EVERY CENTROID
    vector<vector<int>> cluster_table;   //TABLE WITH ALL THE CURVES IN EACH CLUSTER
    vector<vector<int>> hashes;          //INDEXES THAT LEAD CENTROIDS TO HASHTABLE BUCKETS
    vector<int> current_point;
    int radius; 
    int last_id = curve_vector_get_size()-1;
    bool first_iteration = true;
    int new_curves_assigned = 0;
    int previous_curves_assigned = 0;
    int num_of_curves = curve_vector_get_size();
    int complexity = curve_vector_get_curve(0).second.size(); //COMPLEXITY OF CURVES //dimensions previously

    //START COUNTING TIME
    auto start_time = std::chrono::high_resolution_clock::now();
    k_means_plus_plus(k, assignment);                //INITIALIZE K CENTROIDS USING K-MEANS++ ALGORITHM
    radius = centroids_get_radii();      //SET MINIMUM DISTANCE BETWEEN CENTROIDS DIVIDED BY 2 AS FIRST RADIUS

    do{

        //GET EVERY CENTROID'S HASHTABLE BUCKET HASHES
        centroids_get_hashtable_hashes_lsh(g, hashes);
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

            //GET CENTROID'S COORDINATES BY ACCESSING THE POINT VECTOR DATA
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

        //UPDATE THE CENTROIDS
        update_as_vector(cluster_table, last_id);

        first_iteration = false;

    }while(1);

    //WHEN THE CLUSTERS HAVE BEEN DEFINITIVELY FORMED STOP COUNTING TIME
    auto stop_time = std::chrono::high_resolution_clock::now();

}

void reverse_assignment_cube(G_Hypercube g, fstream& output_file, int k, int probes, string assignment, bool complete_flag){

    pair<pair<string, int>, vector<double>> centroid; //HERE CENTROIDS ARE THE QUERY POINTS
    vector<int> appending_curves;        //NEW CURVES THAT WILL BE ADDED TO CLUSTERS   
    vector<pair<vector<int>,int>> curves_in_range; //ALL THE CURVES THAT ARE INSIDE THE GIVEN RADIUS FROM EVERY CENTROID
    vector<vector<int>> cluster_table;   //TABLE WITH ALL THE CURVES IN EACH CLUSTER
    vector<int> hashes;          //INDEXES THAT LEAD CENTROIDS TO HASHTABLE BUCKETS
    vector<int> current_point;
    int radius; 
    int last_id = curve_vector_get_size()-1;
    bool first_iteration = true;
    int new_curves_assigned = 0;
    int previous_curves_assigned = 0;
    int num_of_curves = curve_vector_get_size();
    int complexity = curve_vector_get_curve(0).second.size(); //COMPLEXITY OF CURVES //dimensions previously

    //START COUNTING TIME
    auto start_time = std::chrono::high_resolution_clock::now();
    k_means_plus_plus(k, assignment);                //INITIALIZE K CENTROIDS USING K-MEANS++ ALGORITHM
    radius = centroids_get_radii();      //SET MINIMUM DISTANCE BETWEEN CENTROIDS DIVIDED BY 2 AS FIRST RADIUS

    do{

        //GET EVERY CENTROID'S HYPERCUBE BUCKET HASH
        centroids_get_hypercube_hashes(g, hashes);
        //FOR EVERY CENTROID
        for(int i = 0; i < centroids_get_size(); i++){
            //GET CENTROID'S COORDINATES BY ACCESSING THE POINT VECTOR DATA
            centroid = curve_vector_get_curve(centroids_get_centroid(i));
            if(first_iteration){ //0 IS THE FIRST INDEX OF AN UNASSIGNED CURVE AT START - ALL CURVES UNASSIGNED
                curves_in_range.push_back(make_pair(cube_range_search(hashes[i], radius, probes, complexity, centroid), 0));
            }
            else{
                appending_curves = cube_range_search(hashes[i], radius, probes, complexity, centroid);
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

            //GET CENTROID'S COORDINATES BY ACCESSING THE POINT VECTOR DATA
            centroid = curve_vector_get_curve(centroids_get_centroid(i));
            //AND PERFORM RANGE SEARCH - 0 IS THE FIRST INDEX OF AN UNASSIGNED CURVE AT START
            appending_curves = cube_range_search(hashes[i], radius, probes, complexity, centroid);

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

        //UPDATE THE CENTROIDS
        update_as_vector(cluster_table, last_id);

        first_iteration = false;

    }while(1);

    //WHEN THE CLUSTERS HAVE BEEN DEFINITIVELY FORMED STOP COUNTING TIME
    auto stop_time = std::chrono::high_resolution_clock::now();

}

void reverse_assignment_frechet(G_Frechet g, fstream& output_file, int k, string assignment, bool complete_flag){

    pair<pair<string, int>, vector<double>> centroid; //HERE CENTROIDS ARE THE QUERY POINTS
    vector<int> appending_curves;        //NEW CURVES THAT WILL BE ADDED TO CLUSTERS   
    vector<pair<vector<int>,int>> curves_in_range; //ALL THE CURVES THAT ARE INSIDE THE GIVEN RADIUS FROM EVERY CENTROID
    vector<vector<int>> cluster_table;   //TABLE WITH ALL THE CURVES IN EACH CLUSTER
    vector<vector<int>> hashes;          //INDEXES THAT LEAD CENTROIDS TO HASHTABLE BUCKETS
    vector<int> current_point;
    int radius; 
    int last_id = curve_vector_get_size()-1;
    bool first_iteration = true;
    int new_curves_assigned = 0;
    int previous_curves_assigned = 0;
    int num_of_curves = curve_vector_get_size();
    int complexity = curve_vector_get_curve(0).second.size(); //COMPLEXITY OF CURVES //dimensions previously

    //START COUNTING TIME
    auto start_time = std::chrono::high_resolution_clock::now();
    k_means_plus_plus(k, assignment);    //INITIALIZE K CENTROIDS USING K-MEANS++ ALGORITHM
    radius = centroids_get_radii();      //SET MINIMUM DISTANCE BETWEEN CENTROIDS DIVIDED BY 2 AS FIRST RADIUS
    
    do{

        //GET EVERY CENTROID'S HASHTABLE BUCKET HASHES
        centroids_get_hashtable_hashes_frechet(g, hashes);
        
        //FOR EVERY CENTROID
        for(int i = 0; i < centroids_get_size(); i++){
            //GET CENTROID'S COORDINATES BY ACCESSING THE POINT VECTOR DATA
            centroid = curve_vector_get_curve(centroids_get_centroid(i));
            if(first_iteration){ //0 IS THE FIRST INDEX OF AN UNASSIGNED CURVE AT START - ALL CURVES UNASSIGNED
                
                curves_in_range.push_back(make_pair(frechet_range_search(hashes[i], radius, centroid), 0));
            }
            else{
                appending_curves = frechet_range_search(hashes[i], radius, centroid);
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

            //GET CENTROID'S COORDINATES BY ACCESSING THE POINT VECTOR DATA
            centroid = curve_vector_get_curve(centroids_get_centroid(i));
            //AND PERFORM RANGE SEARCH - 0 IS THE FIRST INDEX OF AN UNASSIGNED CURVE AT START
            appending_curves = frechet_range_search(hashes[i], radius, centroid);

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

//RECEIVES A TABLE OF THE CLUSTERS. EACH ROW CORRESPENDS TO A CLUSTER
//IN EACH ROW ARE STORED THE IDS OF THE INPUT CURVES THAT BELONG TO THAT CLUSTER
//RETURNS A TABLE WITH EACH CLUSTER'S NEW CENTROID
void update_as_vector(vector<vector<int>>& cluster_table, int& last_known_id){

    int row, column;  //ITERATORS
    int dimensions; //THE NUMBER OF COORDINATES AN INPUT CURVE HAS
    vector<double> coordinates_sum;  //A TABLE OF THE SUMS OF THE 1ST, 2ND, ..., NTH COORDINATE OF THE CURVES IN THE SAME CLUSTER
    pair<pair<string, int>, vector<double>> current_curve, mean_curve;
    vector<int> centroids_cp;

    dimensions = curve_vector_get_curve(1).second.size();
    centroids_cp = centroids_get_table();
    centroids_clear();

    for (row=0 ; row < cluster_table.size() ; row++) { //FOR EACH CLUSTER
        coordinates_sum.assign(dimensions, 0); //INITIALIZE ALL SUMS (ONE SUM FOR EACH COORDINATE) WITH 0
        for (column=0 ; column < cluster_table[row].size(); column++) { //FOR EVERY POINT IN THE CLUSTER
            current_curve= curve_vector_get_curve((cluster_table[row][column]));
            coordinates_sum= add_vectors(current_curve.second, coordinates_sum);
        }
        if(non_zero_coordinates(coordinates_sum)){
            //mean_vector.push_back(++last_known_id);
            mean_curve= get_mean_curve_vector(coordinates_sum, cluster_table[row].size(),last_known_id);
            curve_vector_insert_curve(mean_curve);
            //centroids_ids.push_back(mean_vector[0]);
            centroids_insert_curve(mean_curve.first.second);
        }
        else{

            centroids_insert_curve(centroids_cp[row]);
        }

    }
    
}