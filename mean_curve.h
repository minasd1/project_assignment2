#include <utility>
#include <vector>
#include <string>

typedef std::default_random_engine engine;

using namespace std;
pair<pair<string,int>, vector<double>> find_mean_curve_Macchu_Picchu(const int& cluster_size, 
         const vector<int>& cluster, engine gen, double e, int max_length, int& last_known_id);

int create_Macchu_Picchu (const int& cluster_size, 
                         vector <vector <pair<pair<string,int>,vector<double>>>>& curves);

void create_base_layer (engine gen, const int& floors_to_the_peak, const vector<int>& cluster, 
                        vector <vector<pair<pair<string,int>,vector<double>>>>& curves);

pair<pair<string,int>, vector<double>> reach_the_peak (const int& floors_to_the_peak, double e, 
               int max_length, vector <vector<pair<pair<string,int>, vector<double>>>>& curves);

bool compare_distance(const pair<int,int>& pair1, const pair<int, int>& pair2);