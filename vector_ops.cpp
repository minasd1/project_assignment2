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

//INITIALIZE THE CURVES ID VECTOR
void curves_ID_initialize(int num_of_curves, int L){

    curves_ID_vector.resize(num_of_curves, vector<int>(L));
}

//INSERT CURVE'S ID VALUES TO ID VECTOR
//INDEX VALUE IS THE ID VALUE OF THE CURVE (KEY)
void curves_ID_insert(int index_value, vector<int>& curve_id_values){

    for(int i = 0; i < curves_ID_vector[index_value].size(); i++){

        curves_ID_vector[index_value][i] = curve_id_values[i];
    }
}

//GET ID VALUE OF A CURVE CORRESPONDING TO A SPECIFIC HASHTABLE
//K IS THE HASHTABLE THAT CONCERNS US - (NUMBER OF G-FUNCTION USED IN THIS HASHTABLE ({0,1,...,L}))
//INDEX VALUE IS THE KEY VALUE OF THE CURVE
int curves_ID_get_curve_value(int index_value, int k){

    int id_value;
    id_value = curves_ID_vector[index_value][k];

    return id_value;
}

