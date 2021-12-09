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
void curves_ID_vector_initialize(int num_of_curves, int L){

    curves_ID_vector.resize(num_of_curves, vector<int>(L));
}

//USED FOR CLASSIC LSH ALGORITHM
//INSERT CURVE'S ID VALUES (A VECTOR OF IDS) TO ID VECTOR
//INDEX VALUE IS THE ID VALUE OF THE CURVE (KEY)
void curves_ID_vector_insert_lsh(int index_value, vector<int>& curve_id_values){

    for(int i = 0; i < curves_ID_vector[index_value].size(); i++){

        curves_ID_vector[index_value][i] = curve_id_values[i];
    }
}

//USED FOR FRECHET ALGORITHM
//INSERT A CURVE'S ID VALUE TO ID VECTOR
//INDEX VALUE IS THE ID VALUE OF THE CURVE (KEY)
void curves_ID_vector_insert_frechet(int index_value, int curve_id_value, int grid){

        curves_ID_vector[index_value][grid-1] = curve_id_value;
}

//GET ID VALUE OF A CURVE CORRESPONDING TO A SPECIFIC HASHTABLE
//K IS THE HASHTABLE THAT CONCERNS US - (NUMBER OF G-FUNCTION USED IN THIS HASHTABLE ({0,1,...,L}))
//INDEX VALUE IS THE KEY VALUE OF THE CURVE
int curves_ID_vector_get_curve_value(int index_value, int k){

    int id_value;
    id_value = curves_ID_vector[index_value][k];

    return id_value;
}

//RETURNS ALL THE IDS OF THE GIVEN CURVE
//THE FIRST ID CORRESPONDS TO THE FIRST GRID'S HASH TABLE, 
//THE SECOND ID CORRESPONDS TO THE SECOND GRID'S HASH TABLE AND SO ON 
//INDEX VALUE IS THE KEY VALUE OF THE CURVE
vector<int> curves_ID_vector_get_curve_all_ids(int index_value) 
{    
    int i;
    vector<int> curve_ids;

    for (i= 0; i < curves_ID_vector[index_value].size() ; i++) {
        curve_ids.push_back(curves_ID_vector[index_value][i]);
    }
    return curve_ids;
}

void curves_ID_vector_print()
{
    int i, j;

    cout << "ID VECTOR" << endl;
    for (i= 0; i < curves_ID_vector.size() ; i++) {
        cout << "i :   ";
        for (j= 0; j < curves_ID_vector[0].size() ; j++) {
            cout << curves_ID_vector[i][j] << " ";
        }
        cout << endl;
    }
}

//SET NUMBER AND SIZE OF V-VECTORS
void v_vectors_initialization(vector<vector<int>>& v_vectors, int num_of_v_vectors, int dimensions){
    v_vectors.resize(num_of_v_vectors, vector<int>(dimensions));
}

//ASSIGN COORDINANCES TO EACH VECTOR V
//THE COORDINATES COME FROM THE GAUSSIAN DISTRIBUTION ~N(0,1)
void v_vectors_assign_coordinances(vector<vector<int>>& v_vectors, int num_of_vectors, int dimensions,
                                   std::default_random_engine& generator)
{
    int rand_int;
    normal_distribution<double> distribution(0,1);
    //INITIALIZE THE VECTORS
    v_vectors_initialization(v_vectors, num_of_vectors, dimensions);

    //ASSING VALUES USING RANDOM - INSTEAD POISSON DISTRIBUTION SHOULD BE USED
    for (int i = 0; i < v_vectors.size(); i++) {
        for (int j = 0; j < v_vectors[i].size(); j++){
            //EVERY COORDINANCE OF EACH VECTOR IS ASSIGNED TO BE BETWEEN 0 AND 101
            rand_int = distribution(generator)/1;
            v_vectors[i][j] = rand_int;
        }
    }
}

//PRINT CONTENTS OF V-VECTORS - USED FOR CHECKING PURPOSES
void v_vectors_printdata(vector<vector<int>>& v_vectors){
    for (int i = 0; i < v_vectors.size(); i++) {
        for (int j = 0; j < v_vectors[i].size(); j++){
            cout<< v_vectors[i][j] << " ";
        }
        cout << endl;
    }
}

//CREATE A VECTOR t WHICH CONTAINS k RANDOM FLOATS IN RANGE [0,w)
//THE RANDOM FLOATS COME FROM THE UNIFORM DISTRIBUTION ~Unif[0,w)
void create_vector_t(vector<float>& t, int k, int w, std::default_random_engine& generator)
{
    int i;
    float rand_float;
    uniform_real_distribution<float> distribution(0.0,(float)w);
    t.resize(k);
    for (i=0; i < t.size() ; i++) {
        rand_float= distribution(generator);
        t[i]= rand_float;
    }
}

//CREATE A VECTOR OF k RANDOM INTEGERS  IN RANGE [0,w)
//THE RANDOM INEGERS COME FROM THE UNIFORM DISTRIBUTION ~Unif[0,w)
void create_vector_int(vector<int>& ints, int k, int w, std::default_random_engine& generator)
{
    int i;
    int rand_number;
    uniform_int_distribution<> distribution(0, w);

    ints.resize(k);
    for (i=0; i < ints.size() ; i++) {
        rand_number= distribution(generator);
        ints[i]= rand_number;
    }
}

//PRINT CONTENTS OF t-VECTORS - USED FOR CHECKING PURPOSES
void print_vector_t(vector<float>& t){
    for (int i = 0; i < t.size(); i++) {
        cout << t[i] << " ";
    }
    cout << endl;
}

//CALCULATE DISCRETE FRECHET DISTANCE BETWEEN TWO CURVES
double curve_calculate_dfd(const pair<pair<string, int>, vector<double>>& curve1, const pair<pair<string, int>, vector<double>>& curve2){
    vector<vector<double>> DFD;     //ARRAY THAT KEEPS DISCRETE FRECHET DISTANCE FOR EVERY POSSIBLE TRAVERSAL
    vector<double> point1, point2;  //CURVE POINTS
    int num_of_points = curve1.second.size();

    DFD.resize(num_of_points, vector<double>(num_of_points));
    
    //FOR ALL THE VERTICES OF CURVE 1
    for(int i = 0; i < num_of_points; i++){
        //PAIR THEIR X AND Y VALUES
        point1.push_back(i);
        point1.push_back(curve1.second[i]);

        //FOR ALL THE VERTICES IN CURVE 2
        for(int j = 0; j < num_of_points; j++){
            //DO THE SAME
            point2.push_back(j);
            point2.push_back(curve2.second[j]);

            if(i == 0 && j == 0){   //FIRST ELEMENT OF THE ARRAY
                DFD[i][j] = calculate_distance(point1, point2);
            }
            else if(i == 0){        //FIRST ROW OF THE ARRAY
                DFD[i][j] = max(DFD[i][j-1], calculate_distance(point1, point2));
            }
            else if(j == 0){        //FIRST COLUMN OF THE ARRAY
                DFD[i][j] = max(DFD[i-1][j], calculate_distance(point1,point2));
            }
            else{                   //RANDOM ARRAY INDEX
                DFD[i][j] = max(min({DFD[i-1][j], DFD[i-1][j-1], DFD[i][j-1]}), calculate_distance(point1, point2));
            }

            point2.clear();
        }

        point1.clear();
    }

    return DFD[num_of_points-1][num_of_points-1];

}

//CONVERT A CURVE TO USE IT FOR CONTINUOUS FRECHET
Curve convert_for_continuous_frechet(const pair<pair<string, int>, vector<double>>& curve, const unsigned long dimensions){

    Points points(dimensions);

    for(int i = 0; i < curve.second.size(); i++){

        Point point(dimensions);
        point.set(dimensions-1, curve.second[i]);

        points.add(point);
    }

    Curve modified_curve(points, curve.first.first);

    return modified_curve;
}

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
vector<double> add_vectors(const pair<pair<string, int>, vector<double>>& curve1, const pair<pair<string, int>, vector<double>>& curve2)
{
    int i;
    vector<double> sum_curve;

    if (curve1.second.size() != curve2.second.size()) {
        cerr << "Error in add_vectors: Can not add vectors of different size" << endl;
        sum_curve.assign(curve1.second.size(), -666);
    }
    else {
        sum_curve.assign(curve1.second.size(), 0);
        for (i=0 ; i < curve1.second.size(); i++) {
            sum_curve[i]= curve1.second[i] + curve2.second[i];
        }
    }
    return sum_curve;
}

//COMPUTES THE DISTANCE BETWEEN 2 VECTORS USING THE k-NORM
double calculate_distance(vector<double>& point1, const vector<double>& point2, int k)
{
    double distance = 0.0;
    double sum = 0;

    for (int i=0 ; i < point1.size() ; i++) {
        sum+= pow(abs(point1[i]-point2[i]), k);
    }
    distance = pow(sum, 1.0/(double)k);

    return distance;
}
