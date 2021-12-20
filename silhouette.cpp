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
    double sum_a, sum_b, a, b, silhouette, cluster_mean_silhouette, 
           dataset_mean_silhouette, cluster_sum, dataset_sum, distance;
    pair<pair<string, int>, vector<double>> current_curve, other_curve;
    dataset_sum= 0.0;
    output_file << "[ ";

    for (i= 0; i < cluster_table.size() ; i++) { //FOR EVERY CLUSTER
       
        cluster_sum= 0.0;
     
        for (j= 0; j < cluster_table[i].size() ; j++) { //FOR EVERY CURVE IN THIS CLUSTER
           
            sum_a= 0.0;
            current_curve= curve_vector_get_curve(cluster_table[i][j]);
           
            for (k= 0; k < cluster_table[i].size() ; k++) { //FOR EVERY CURVE IN THE SAME CLUSTER

                other_curve= curve_vector_get_curve(cluster_table[i][k]);
               
                if(update == "Mean_Vector"){
                    distance= calculate_distance(current_curve.second, other_curve.second);
                   
                    sum_a+= distance;
                   
                }
                else if(update == "Mean_Frechet"){
                    distance= curve_calculate_dfd(current_curve.second, other_curve.second);
                   
                    sum_a+= distance;
                }
                
                
            }
          
            if (cluster_table[i].size()==1 ) {
                //IF THERE IS ONLY ONE ELEMENT IN THE CLUSTER THE DISTANCE WITH ITSELF IS ZERO
                a = 0;
            }
            else {//IF THERE ARE MORE THAN ONE CURVES IN THE CLUSTER
                a= sum_a/double(cluster_table[i].size()-1);
            }
            
            second_nearest= find_second_nearest_centroid(current_curve, update);
            sum_b= 0.0;
         
            for (k= 0; k < cluster_table[second_nearest].size(); k++) { //FOR EVERY CURVE IN THE SECOND NEAREST CLUSTER
               
                other_curve= curve_vector_get_curve(cluster_table[second_nearest][k]);
                if(update == "Mean_Vector"){
                    distance= calculate_distance(current_curve.second, other_curve.second);
                   
                    sum_b+= distance;
                }
                else if(update == "Mean_Frechet"){
                    distance+= curve_calculate_dfd(current_curve.second, other_curve.second);
                   
                    sum_b+= distance;
                }
               
            }
            if (cluster_table[second_nearest].size() == 0 ) {
               
                // silhouette= 0.0;
                b = 0;
                if(a == 0){
                   
                   silhouette= 0.0; 
                }
                else{
                   
                    silhouette = -1.0;
                }
            }
            else {
                b= sum_b/(double)(cluster_table[second_nearest].size());
                
                silhouette= (b-a)/max_double(a,b);
            }
           
            cluster_sum+= silhouette;

        }
     
        if (cluster_table[i].size() == 0) {

            cluster_mean_silhouette= 0.0;
        }
        else {

            cluster_mean_silhouette= cluster_sum/double(cluster_table[i].size());
        }
        output_file << cluster_mean_silhouette << ", ";
        dataset_sum+= cluster_mean_silhouette;
    }
   
    dataset_mean_silhouette= dataset_sum/double(cluster_table.size());
    output_file << dataset_mean_silhouette << "]" << endl;
}