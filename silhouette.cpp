#include <fstream>
#include <iostream>
#include <vector>
#include "vector_ops.h"
#include "lloyds_auxiliary.h"

using namespace std;

//PRINTS THE SILHOUETTE'S DATA
void print_silhouette(vector<vector<int>> cluster_table, fstream& output_file, string update)
{
    int i, j, k, second_nearest;
    double sum_a, sum_b, a, b, silhouette, cluster_mean_silhouette, dataset_mean_silhouette, cluster_sum, dataset_sum;
    pair<pair<string, int>, vector<double>> current_curve, other_curve;

    dataset_sum= 0.0;
    output_file << "[ ";

    for (i= 0; i < cluster_table.size() ; i++) { //FOR EVERY CLUSTER
        cluster_sum= 0.0;
        for (j= 0; j < cluster_table[i].size() ; j++) { //FOR EVERY POINT IN THIS CLUSTER
            sum_a= 0.0;
            current_curve= curve_vector_get_curve(cluster_table[i][j]);
            for (k= 0; k < cluster_table[i].size(); k++) { //FOR EVERY CURVE IN THE SAME CLUSTER
                other_curve= curve_vector_get_curve(cluster_table[i][k]);
                sum_a+= calculate_distance(current_curve.second, other_curve.second);
            }
            if (cluster_table[i].size()==1) {
                a= 0; //IF THERE IS ONLY ONE ELEMENT IN THE CLUSTER THE DISTANCE WITH ITSELF IS ZERO
            }
            else {
                a= sum_a/double(cluster_table[i].size()-1);
            }
            //cout << "Sum_a = " << sum_a << " and a= "<< a<<  endl;
            second_nearest= find_second_nearest_centroid(current_curve, update);
            //cout << "second_nearest= " << second_nearest << endl;
            sum_b= 0.0;
            for (k= 0; k < cluster_table[second_nearest].size(); k++) { //FOR EVERY CURVE IN THE SECOND NEAREST CLUSTER
                other_curve= curve_vector_get_curve(cluster_table[second_nearest][k]);
                sum_b+= calculate_distance(current_curve.second, other_curve.second);
            }
            b= sum_b/double(cluster_table[second_nearest].size());
            //cout << "Sum_b= " << sum_b << " and b= " << b << endl;
            silhouette= (b-a)/max_double(a,b);
            cluster_sum+= silhouette;
        }
        cluster_mean_silhouette= cluster_sum/double(cluster_table[i].size());
        output_file << cluster_mean_silhouette << ", ";
        dataset_sum+= cluster_mean_silhouette;
    }
    dataset_mean_silhouette= dataset_sum/double(cluster_table.size());
    output_file << dataset_mean_silhouette << "]" << endl;
}