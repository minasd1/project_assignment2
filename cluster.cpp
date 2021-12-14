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

    centroids_print_data();
}