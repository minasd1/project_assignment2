#ifndef CMD_LINE_ARGS_H
#define CMD_LINE_ARGS_H
#include <string>

using namespace std;

int read_cmd_args(int argc, char** argv, string& input_file, string& query_file,
                      int& k_lsh, int& k_cube, int& l, string& output_file, int& n,
                      float& r, int& m, int& probes, string& config_file, string& assignment, 
                      bool& complete_flag, string& algorithm, string& metric, double& delta, 
                      string& update, bool& datapath_given, bool& query_given, bool& output_given, 
                      bool& algorithm_given, bool& metric_given);
#endif
