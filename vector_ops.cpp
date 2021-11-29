#include "vector_ops.h"

//INSERT A CURVE TO CURVE VECTOR
void curve_vector_insert_curve(pair<pair<string, int>, vector<double>>& curve){

    curve_vector.push_back(curve);
}

//GET A CURVE FROM CURVE VECTOR
pair<pair<string, int>, vector<double>> curve_vector_get_curve(int curve_id){

    return curve_vector[curve_id];
}

//GET NUMBER OF CURVES IN CURVE VECTOR
int curve_vector_get_size(){

    return curve_vector.size();
}

//PRINT CURVE VECTOR VALUES - USED FOR CHECKING PURPOSES
void curve_vector_print_values(){

    for(int i = 0; i < curve_vector.size(); i++){

        cout << curve_vector[i].first.second << " " << curve_vector[i].first.first;
        for(int j = 0; j < curve_vector[i].second.size(); j++){

            cout << " " << curve_vector[i].second[j];
        }

        cout << endl;
    }
    
}

