#include "vector_ops.h"

//CALCULATE DOT PRODUCT OF TWO VECTORS
int calculate_dot_product(const pair<pair<string, int>, vector<double>>& curve, vector <int>& d_vector){
    double product = 0.0;
    for(int i = 0; i < curve.second.size(); i++){
        product = product + curve.second[i] * d_vector[i];
    }
    return (int)product;
}

//CALCULATES THE ADDITION OF TWO DOUBLE VECTORS OF THE SAME SIZE
//EG. v1= [0.0, 2.0, 4.0, 8.0]  v2= [0.0, 1.0, -1.0, -5.0]  v1+v2= v3= [0.0, 3.0, 3.0, 3.0]
vector<double> add_vectors(const vector<double>& curve1, const vector<double>& curve2)
{
    int i;
    vector<double> sum_curve;

    if (curve1.size() != curve2.size()) {
        cerr << "Error in add_vectors: Can not add vectors of different size" << endl;
        sum_curve.assign(curve1.size(), -666);
    }
    else {
        sum_curve.assign(curve1.size(), 0);
        for (i=0 ; i < curve1.size(); i++) {
            sum_curve[i]= curve1[i] + curve2[i];
        }
    }
    return sum_curve;
}

//COMPUTES THE DISTANCE BETWEEN 2 VECTORS USING THE k-NORM
double calculate_distance(vector<double>& curve1, const vector<double>& curve2, int k)
{
    double distance = 0.0;
    double sum = 0;

    for (int i=0 ; i < curve1.size() ; i++) {
        sum+= pow(abs(curve1[i]-curve2[i]), k);
    }
    distance = pow(sum, 1.0/(double)k);

    return distance;
}

//CHECK IF A VECTOR HAS NON ZERO COORDINATES OR NOT
bool non_zero_coordinates(vector<double>& coordinates){

    for(int i = 0; i < coordinates.size(); i++){
        if(coordinates[i] != 0){
            return true;
        }
    }

    return false;
}

//GET MEAN CURVE BETWEEN TWO CURVES AS A VECTOR
pair<pair<string, int>, vector<double>> get_mean_curve_vector(vector<double> vector_of_sums, int num_of_vectors, int& last_known_id)
{

    double mean_value;
    pair<pair<string, int>, vector<double>> mean_curve;

    mean_curve.first.second = ++last_known_id;
    for(int i = 0; i < vector_of_sums.size(); i++){
        mean_value = vector_of_sums[i]/num_of_vectors;
        mean_curve.second.push_back(mean_value);
    }

    return mean_curve;
}

