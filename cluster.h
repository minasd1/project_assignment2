#ifndef CLUSTER_H
#define CLUSTER_H

#include <iostream>
#include <vector>

#include "vector_ops.h"
#include "lsh.h"
#include "cube.h"

using namespace std;

void k_means_plus_plus(int k);
void reverse_assignment_lsh(G_Lsh g, fstream& output_file, int k, bool complete_flag);
void reverse_assignment_cube(G_Hypercube g, fstream& output_file, int k, int probes, bool complete_flag);












#endif