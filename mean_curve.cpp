#include <cmath>
#include <random>
#include <algorithm>
#include "vector_ops.h"
#include "mean_curve.h"

using namespace std;

//TAKES A CLUSTER AND RETURNS ITS MEAN CURVE USING THE TREE TECHNIQUE BUT IN A DIFFERENT WAY
pair<pair<string,int>, vector<double>> find_mean_curve_Macchu_Picchu(const int& cluster_size, const vector<int>& cluster, engine gen, double e, int max_length, int& last_known_id) 
{
    pair<pair<string,int>, vector<double>> current_curve, mean_curve;
   
    //BEGGINING FROM THE BOTTOM LAYER, HOW MANY FLOORS ABOVE US IS THE PYRAMID'S PEAK
    int floors_to_the_peak; 
    //THE PYRAMID'S STRUCT FOR STORING THE MEAN CURVES PRODUCED BY THE ORIGINAL CURVES OF THE CLUSTER
    vector <vector <pair<pair<string,int>, vector<double>>>> curves; 

    //CREATE THE PYRAMID STRUCT AND RETURN IT'S HIGHEST FLOOR'S NUMBER (BASE LAYER IS FLOOR 0)
    floors_to_the_peak= create_Macchu_Picchu(cluster_size, curves);
    
    //INITIALIZE RANDOMLY THE BASE LAYER WITH CURVES FROM THE CLUSTER
    create_base_layer(gen, floors_to_the_peak, cluster, curves);
    
    //FILL THE PYRAMID LAYER BY LAYER UNTIL THE TOP WHERE THE MEAN CURVE OF THE CLUSTER WILL BE FOUND
    mean_curve= reach_the_peak(floors_to_the_peak, e, max_length, curves);
    mean_curve.first.first= "MEAN";
    mean_curve.first.second= ++last_known_id;

    return mean_curve;
} 


int create_Macchu_Picchu (const int& cluster_size, vector <vector<pair<pair<string,int>, vector<double>>>>& curves)
{
    int floors_to_the_peak, i;

    //FIND THE PYRAMID'S HIGH
    floors_to_the_peak = 0;
    while (pow(2.0, floors_to_the_peak) < cluster_size) {
        floors_to_the_peak++;
    }

    //CONSTRUCT THE STRUCT
    curves.resize(floors_to_the_peak+1); //THE PYRAMID HAS floors_to_the_peak+1 LAYERS
    for (i= 0; i <= floors_to_the_peak ; i ++) { //EACH OF THE LAYER HAS 
        curves[i].resize(pow(2.0 , i));    //2^i CURVES
    }

    return floors_to_the_peak;
}
 
void create_base_layer (engine gen, const int& floors_to_the_peak, const vector<int>& cluster, vector <vector<pair<pair<string,int>,vector<double>>>>& curves)
{
    int i,j;
    pair<int,int> current_pair, temp;
    vector <pair<int,int>> pairs;
    pair<pair<string,int>, vector<double>> current_curve;
    int limit= pow(cluster.size(), 3);
    uniform_int_distribution<int> distribution(0,limit);

    //CREATE PAIRS OF INTEGERS {TABLE_INDEX, VALUE} , WHERE TABLE_INDEX=[0, 1, ...], VALUE=[RANDOM(0,limit)]
    //AND SORT THEM ACCORDING TO THE SECOND OF THE PAIR (THE RANDOM VALUE)
    //GOAL IS TO MIX UP THE TABLE_INDEXES (E.G. INSTEAD OF [0,1,2,3,4,5] ---> [3,5,1,4,2,0])
  
    for (i= 0 ; i < cluster.size() ; i++) {
      
        current_pair.first= i;
        current_pair.second= distribution(gen);

        pairs.push_back(current_pair);
    }
    //SORTING THE PAIRS BY SECOND INTEGER ASCENDING ORDER
    sort(pairs.begin(), pairs.end(), compare_distance);

    //INITIALIZE THE BASE LAYER WITH EMPTY POSITIONS
    for (i= 0;  i < curves[floors_to_the_peak].size() ; i++) {
        current_curve.first.first= "NULL";
        current_curve.first.second= -1;
        curves[floors_to_the_peak][i]= current_curve;
    }

    //NOW OVERWRITE THE BASE LAYER WITH REAL CURVES FROM THE CLUSTER
    //TAKE EACH CURVE FROM THE CLUSTER TABLE  
    //AND PLACE IT WHERE THE FIRST INTEGER OF THE PAIR IN THE SAME INDEX IN THE PAIRS TABLE "POINTS" 
    //E.G. CLUSTER_TABLE:[curve1, curve2, ...,curveN]
    //     PAIRS_TABLE: [{4,120}, {18, 140}, {3, 149}, ...]
    //     THEN curve1 GOES TO curves[base_layer][4], curve2 TO curves[curve's_cluster][18] AND SO ON
    for (i= 0 ; i < cluster.size() ; i++) {
       
        current_curve= curve_vector_get_curve(cluster[i]); //TAKE THE iTH CURVE OF THE CLUSTER
      
        curves[floors_to_the_peak][pairs[i].first]= current_curve;
       
    }
} 

pair<pair<string,int>, vector<double>> reach_the_peak (const int& floors_to_the_peak, double e, 
                int max_length, vector <vector<pair<pair<string,int>, vector<double>>>>& curves)
{
    int i, j;
    pair <pair<string, int>, vector<double>> mean_curve;

    for (i= floors_to_the_peak ; i > 0 ; i--) {
        
        for (j= 0 ; j < curves[i].size() ; j+=2) {
            //IF BOTH CURVES EXIST
            if (curves[i][j].first.first != "NULL" && curves[i][j+1].first.first != "NULL") {
                
                mean_curve= get_mean_curve(curves[i][j].second,  curves[i][j+1].second);
                if (mean_curve.second.size() > max_length) {
                    filter(mean_curve.second, e,  max_length);
                } 
                curves[i-1][j/2]= mean_curve;
              
                curves[i-1][j/2].first.first= "MEAN";
                
            }
            //IF BOTH ARE EMPTY
            else if (curves[i][j].first.first == "NULL" && curves[i][j+1].first.first == "NULL") {
             
                curves[i-1][j/2].first.first= "NULL";
             
            }
            //IF THE LEFT CURVE EXISTS AND THE RIGHT IS EMPTY
            else if (curves[i][j].first.first != "NULL" && curves[i][j+1].first.first == "NULL"){ 
             
                curves[i-1][j/2]= curves[i][j];
             
            }
            
        }
        curves[i].clear();
    }
    return curves[0][0];
}


bool compare_distance(const pair<int,int>& pair1, const pair<int, int>& pair2)
{
    return (pair1.second < pair2.second);
}