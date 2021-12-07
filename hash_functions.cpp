#include <iostream>
#include <random>
#include <cmath>
#include "hash_functions.h"
#include "vector_ops.h"
#include "mod.h"
#include "hashTable.h"

//CONSTRUCTOR OF ALL THE G FUNCTIONS
//GENERATES THE RANDOM V_VECTOR AND t FLOATS USED IN H FUNCTION
//GENERATES THE r INTEGERS
//ALL THESE 3 RANDOMLY GENERATED ELEMENTS COMBINED PRODUCES THE G HASH FUNCTIONS
G_Lsh::G_Lsh(int k_num, int dim, engine gen, int win, int m_mod, int tab_s, int l_num)
    : k(k_num), h_num(4*k_num), dimension(dim), generator(gen), w(win), m(m_mod),
      table_size(tab_s), l(l_num)
{
    int i, j, limit, rand_number;
    limit= pow(2,30); //AN INTEGER <= 32 BITS
    //THE k H FUNCTIONS AND THE k r INTEGERS WILL BE CHOSEN BY THE UNIFORM DISTRIBUTION
    uniform_real_distribution<float> h_distribution(0.0,(float)h_num); // FOR H FUNCTIONS
    uniform_real_distribution<float> r_distribution(0.0,(float)limit); // FOR r VALUES

    //INITIALIZE THE V_VECTORS AND t FLOATS
    v_vectors_assign_coordinances(v_vectors, h_num, dimension, generator);
    create_vector_t(t, h_num, w, generator);
    v_vectors_initialization(h_functions, l, k);
    v_vectors_initialization(r, l, k);

    for (i=0; i < l ; i++){  //FOR EVERY G FUNCTION
        for (j=0; j < k ; j++) { //CHOOSE RANDOMLY k OUT OF 4*k H FUNCTIONS AND k r INTEGERS TO BE USED
            rand_number= h_distribution(generator)/1;
            h_functions[i][j]= rand_number;
            rand_number= r_distribution(generator)/1;
            r[i][j]= rand_number;
        }
    }
}

//RETURNS A SINGLE ID (IN FRECHET'S ALGORITHM CASE) OR A VECTOR THAT CONTAINS
//THE ID OF A GIVEN CURVE FOR EVERY HASH FUNCTION (FOR CLASSIC LSH ALGORITHM)
//THE LAST ARGUMENT IS USED TO DEFINE THE ALGORITHM USED (POSITIVE VALUE FOR FRECHET, ZERO FOR CLASSIC LSH)
//IN FRECHET'S ALGORITHM THE LAST ARGUMENTS VALUE IS THE HASH_GRID OF WHICH THE INSERTED CURVE WAS PRODUCED (NUMBERING BEGINS FROM 1)
void G_Lsh::id(const pair<pair<string, int>, vector<double>>& curve, vector<int>& id_vector, bool is_query, int frechet_grid)
{
    vector<int> curve_id_values;
    int i, j, sum, h, id;

    for (i=0; i < l ; i++) {  //FOR EVERY G FUNCTION
        sum= 0;
        for (j=0 ; j < k ; j++) { //ADD k r*h(p) PRODUCTS
            h= (calculate_dot_product(curve, v_vectors[h_functions[i][j]]) + t[h_functions[i][j]])/w;
            sum+= mod(r[i][j], m) * mod(h, m);
        }
        id= mod(sum, m);
        //cout << "MOD (" << sum << " " << m << ")= " << id << endl; //CHECKING PURPOSE, TO BE REMOVED
        curve_id_values.push_back(id);

		//IF THE FUNCTION IS USED FOR THE FRECHET ALGORITHM THEN ONLY 1 ID IS ENOUGH
        if (frechet_grid != 0) {
        	break;
		}
    }

    if(!is_query){
		if (frechet_grid != 0) {
			curves_ID_vector_insert_frechet(curve.first.second, curve_id_values[0], frechet_grid);
		}
		else {
			curves_ID_vector_insert_lsh(curve.first.second, curve_id_values);
		}
    }
    else {
		id_vector = curve_id_values;
    }

}

//RETURNS A SINGLE HASH VALUE (IN FRECHET'S ALGORITHM CASE) OR A VECTOR THAT CONTAINS
//THE HASH VALUE OF A GIVEN CURVE FOR EVERY HASH FUNCTION (FOR CLASSIC LSH ALGORITHM)
//THE LAST ARGUMENT IS USED TO DEFINE THE ALGORITHM USED (POSITIVE VALUE FOR FRECHET, ZERO FOR CLASSIC LSH)
//IN FRECHET'S ALGORITHM THE LAST ARGUMENTS VALUE IS THE HASH_GRID OF WHICH THE INSERTED CURVE WAS PRODUCED (NUMBERING BEGINS FROM 1)
void G_Lsh::hash(const pair<pair<string, int>, vector<double>>& curve, vector<int>& hash_vector, bool is_query, int frechet_grid)
{
    int i;
    vector<int> curve_hash_values;
    vector<int> curve_id_values;


    this->id(curve, curve_id_values, true, frechet_grid);
    for (i=0; i< curve_id_values.size() ; i++) {
        curve_hash_values.push_back(mod(curve_id_values[i], table_size));
    }
    if (!is_query && !frechet_grid){
    	
		hashTable_push_back(curve_hash_values, curve.first.second);

    }
    else {
    	hash_vector = curve_hash_values;
	}
}

//PRINTS THE V_VECTORS, t, r (THE DATA OF EACH G HASH FUNCTION)
//FOR DEBUGGING PURPOSES
void G_Lsh::print_hash_functions_data(void)
{
    int i, j;

    cout << "Total number of g hash functions: " << l << endl;
    cout << "Each g hash function uses " << k << " h functions" << endl;
    cout << "    v                                  t" << endl;
    cout << "h";
    for (i=0 ; i < h_num ; i++) {
        for (j=0; j < dimension ; j++ ){
            cout << v_vectors[i][j] << " " ;
        }
        cout << "t= " << t[i] << endl;
    }
}

/*
//CONSTRUCTOR OF THE HYPERCUBE'S HASH FUNCTION
//GENERATES d' RANDOM h HASH FUNCTIONS
//GENERATES d' RANDOM f FUNCTIONS
//(ACTUALLY ONE f FUNCTION THAT RECEIVE d' DIFFERENT VALUES IN ONE OF ITS ARGUMENTS)
G_Hypercube::G_Hypercube(int dim, engine gen, int win, int str_dim)
    : dimension(dim), generator(gen), w(win), string_dimensions(str_dim)
{
    int limit, rand_number;

    limit= pow(2, 30); //A BIG INTEGER, THE MAX NUMBER THE UNIFORM DISTRIBUTION CAN PRODUCE (OR DI'S SECRETARY PHONE NUMBER)

    //INITIALIZE THE V_VECTORS AND t FLOATS (USED IN h FUNCTIONS)
    v_vectors_assign_coordinances(v_vectors, string_dimensions, dimension, generator);
    create_vector_t(t, string_dimensions, w, generator);
    //INITIALIZE d' RANDOM INTEGERS THAT WILL BE PASSED IN f FUNCTION
    create_vector_int(rand_ints, string_dimensions, limit, generator);
}

//MAPS A GIVEN NUMBER (h(p)) IN {0,1} UNIFORMLY
int f(int& h, int& random_int_2)
{
    int seed, random_int_1, xor_int, result;

    //SET THE UNIFORM DISTRIBUTION PARAMETERS
    seed = h;
    default_random_engine generator(seed);
    uniform_int_distribution<> distribution(0, 2107275161);

    random_int_1= distribution(generator);
    xor_int= random_int_1 ^ random_int_2;
    result= mod(xor_int, 2);
    return result;
}

//RECEIVES AND MODIFIES ARGUMENT hash_value SO THAT IT SHOWS IN WHICH VERTEX
//OF THE HYPERCUBE THE GIVEN POINT CORRESPONDS TO
//IF is_query IS false THE RECEIVED POINT point IS STORED IN THE HYPERCUBE'S HASH TABLE
void G_Hypercube::hash(const vector<int>& point, unsigned int& hash_value, bool is_query)
{
    int i, h, result;

    hash_value= 0; // 0000...000
    for (i=0; i < string_dimensions ; i++){
        h= (calculate_dot_product(point, v_vectors[i]) + t[i])/w;
        result= f(h, rand_ints[i]); //0 or 1
        hash_value <<= 1;
        hash_value+= result; // Last bit gets the value of result
    }
    if (!is_query) {
        hyperCube_push_back(hash_value, point[0]);
    }
}
*/

//GENERATES VECTOR OF FLOATS USED IN GRID CURVE
G_Frechet::G_Frechet(G_Lsh g_func, engine gen, int L_num, double delta_value, int num_grid_values)
    :g(g_func), generator(gen), L(L_num), delta(delta_value), num_of_grid_values(num_grid_values)
{
    //INITIALIZE A FLOAT VECTOR TO DISTORT THE GRID CURVE
    //create_vector_t(t, num_of_grid_values, delta, generator);

} 

void G_Frechet::hash(const pair<pair<string, int>, vector<double>>& curve, vector<int>& hash_vector, bool is_query)
{

    pair<pair<string, int>, vector<double>> snapped_curve;
    vector<int> hash_values;

    //SNAPPED CURVE MUST HAVE THE SAME STRING AND ID VALUE WITH THE INITIAL CURVE
    snapped_curve.first.first = curve.first.first;
    snapped_curve.first.second = curve.first.second;

    //FOR EVERY HASHTABLE
    for(int i = 1; i <= L; i++){
        //PRODUCE A GRID
        create_vector_t(t, num_of_grid_values, delta, generator);

        //SNAP CURVE TO CURRENT GRID
        snap_to_grid(curve, snapped_curve.second); //snapped_curve
        padding(snapped_curve.second);

        //GET THE BUCKET THAT CURVE MUST BE INSERTED
        g.hash(snapped_curve, hash_values, 0, i);       //i CAN POSSIBLY BE A  BOOLEAN FLAG
        //AND ADD IT TO CURVE'S HASH VALUES
        hash_vector.push_back(hash_values[0]);

        snapped_curve.second.clear();
    }
    
    //INSERT CURVE'S ID TO EACH ONE OF THE L HASHTABLES INDICATED BY THE HASH_VECTOR
    hashTable_push_back(hash_vector, curve.first.second);

    hash_vector.clear();
}

//GET A CURVE FROM INPUT AND SNAP IT TO GRID
//GET THE SNAPPED CURVE PRODUCED IN THE PROCESS
void G_Frechet::snap_to_grid(const pair<pair<string, int>, vector<double>>& curve, vector<double>& snapped_curve)
{
    double delta_multiple, previous_delta_multiple;                     
    double previous_distance, current_distance;
    double grid_previous_x_value, grid_previous_y_value;
    double grid_current_x_value = 0, grid_current_y_value = 0;
    double curve_current_x_value;
    int count;

    //FOR ALL THE CURVE'S VERTICES
    for(int i = 0; i < curve.second.size(); i++){
        
        count = 0;
        delta_multiple = t[count];  //GRID COORDINATES ARE DISTORTED BY VECTOR t SO INSTEAD
                                    //OF ZERO THE FIRST VALUE OF delta_multiple IS t[0]
        previous_delta_multiple = delta_multiple;
        curve_current_x_value = i;

        //SNAP CURVE'S CURRENT X_VALUE TO AN X_VALUE OF THE GRID
        current_distance = abs(curve_current_x_value - delta_multiple);
        previous_distance = current_distance;

        //WHILE DISTANCE BETWEEN CURVE'S AND GRID'S CURRENT X_VALUE KEEPS DECREASING
        //OR IS NOT ZERO
        while((current_distance <= previous_distance) && (current_distance != 0)){
            count++;
            previous_delta_multiple = delta_multiple;
            //GO TO THE NEXT X_VALUE OF GRID CURVE
            delta_multiple += (delta+t[count]);       

            previous_distance = current_distance;
            current_distance = abs(curve_current_x_value - delta_multiple); 
        }
        
        //KEEP TRACK OF THE PREVIOUS GRID X VALUE - NEEDED TO AVOID DUPLICATES
        grid_previous_x_value = grid_current_x_value;
        //SNAP CURVE'S CURRENT X COORDINATE TO CLOSEST GRID X VALUE
        grid_current_x_value = previous_delta_multiple;

        count = 0;
        delta_multiple = t[count];
        previous_delta_multiple = delta_multiple;

        //DO THE SAME THING FOR CURVE'S CURRENT Y COORDINATE
        current_distance = abs(curve.second[i] - delta_multiple);
        previous_distance = current_distance;

        while((current_distance <= previous_distance) && (current_distance != 0)){
            count++;

            previous_delta_multiple = delta_multiple;
            delta_multiple += (delta+t[count]);                      

            previous_distance = current_distance;
            current_distance = abs(curve.second[i] - delta_multiple);
        }

        //KEEP TRACK OF THE PREVIOUS GRID Y VALUE - NEEDED TO AVOID DUPLICATES
        grid_previous_y_value = grid_current_y_value;
        //SNAP CURVE'S CURRENT Y COORDINATE TO CLOSEST GRID Y VALUE
        grid_current_y_value = previous_delta_multiple;
        
        //IF AT LEAST ONE CURVE VERTICE HAS BEEN ASSIGNED TO THE PROPER GRID VERTICE
        if(snapped_curve.size() >= 2){

            //AVOID CONSECUTIVE DUPLICATES OF GRID VERTICES WHILE SNAPPING CURVE
            if((grid_previous_x_value != grid_current_x_value) || (grid_previous_y_value != grid_current_y_value)){
                //APPEND GRID'S CLOSEST (X,Y) COORDINATES TO CURRENT CURVE'S COORDINATES TO SNAPPED_CURVE
                snapped_curve.push_back(grid_current_x_value);
                snapped_curve.push_back(grid_current_y_value);
            }
        }
        else{   //THE FIRST GRID VERTICE TO BE INSERTED TO GRID CURVE VALUES

            snapped_curve.push_back(grid_current_x_value);
            snapped_curve.push_back(grid_current_y_value);
        }

    }
}

//GET CURVE AFTER SNAPPING TO GRID
//ADD WITH A LARGE VALUE UNTIL MAXIMUM SIZE IS REACHED
void G_Frechet::padding(vector<double>& snapped_curve){

    double padding_value = 500.0;

    //IF SNAPPED CURVE HAS LESS THAN LESS THAN THE MAXIMUM VALUES
    if(snapped_curve.size() < num_of_grid_values){

        for(int i = snapped_curve.size(); i < num_of_grid_values; i++){
            //FILL WITH A VERY LARGE VALUE UNTIL MAXIMUM SIZE IS REACHED
            snapped_curve.push_back(padding_value);
        }
    }
}
