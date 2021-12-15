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
void reverse_assignment_lsh(G_Lsh g, fstream& output_file, int k, string assignment, bool complete_flag);
void reverse_assignment_cube(G_Hypercube g, fstream& output_file, int k, int probes, string assignment, bool complete_flag);
void reverse_assignment_frechet(G_Frechet g, fstream& output_file, int k, string assignment, bool complete_flag);







#endif