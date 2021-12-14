#ifndef CLUSTER_H
#define CLUSTER_H

#include <iostream>
#include <vector>

#include "vector_ops.h"
#include "lsh.h"

using namespace std;

void k_means_plus_plus(int k);
void reverse_assignment_lsh(G_Lsh g, fstream& output_file, int k, bool complete_flag);












#endif