#ifndef USER_INPUT_H
#define USER_INPUT_H

#include <iostream>
#include <string>

using namespace std;



void read_user_input(string& query_file, int* continue_execution);
void read_path(string& dataset_path, string& queryset_path, string& output_file, string& algorithm, 
               string& metric, bool datapath_given, bool query_given, bool output_given, 
               bool algorithm_given, bool metric_given);





#endif