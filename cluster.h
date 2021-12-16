#ifndef CLUSTER_H
#define CLUSTER_H

#include <iostream>
#include <vector>

#include "vector_ops.h"
#include "lsh.h"
#include "cube.h"
#include "frechet_functions.h"

using namespace std;

void k_means_plus_plus(int k, string assignment);
void lloyds(int number_of_clusters, fstream& output_file, string assignment, bool complete_flag);
void reverse_assignment_lloyds(vector<vector<int>>& cluster_table, int number_of_clusters, int num_of_curves, 
                                                                                            int last_id);
void reverse_assignment_lsh(G_Lsh g, fstream& output_file, int k, string assignment, bool complete_flag);
void reverse_assignment_cube(G_Hypercube g, fstream& output_file, int k, int probes, string assignment, bool complete_flag);
void reverse_assignment_frechet(G_Frechet g, fstream& output_file, int k, string assignment, bool complete_flag);
void update_as_vector(vector<vector<int>>& cluster_table, int& last_known_id);







#endif