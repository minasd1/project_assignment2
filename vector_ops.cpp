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

void centroids_insert_curve(int id){

    centroids.push_back(id);
}

int centroids_get_centroid(int index){

    return centroids[index];
}

//GET NUMBER OF CENTROID POINTS
int centroids_get_size(){

    return centroids.size();
}

vector<int> centroids_get_table(){

    vector<int> centroid_ids;
    for(int i = 0; i < centroids.size(); i++){

        centroid_ids.push_back(centroids_get_centroid(i));
    }

    return centroid_ids;
}

//GET THE MINIMUM DISTANCE BETWEEN CENTROIDS
//USED AS A STARTING RANGE IN REVERSE ASSIGNMENT
double centroids_get_radii(){

    double min_distance = numeric_limits<double>::max();
    double current_distance;
   
    for (int i = 0; i < centroids.size(); i++){

        for(int j = 0; j < centroids.size(); j++){

            if(i != j){

                current_distance = calculate_distance(curve_vector[centroids[i]].second, curve_vector[centroids[j]].second, 2);
                if(current_distance < min_distance){

                    min_distance = current_distance;
                }
            }
        }
    }

    return min_distance/2;
}

void centroids_pick_first_centroid(){

    int first_centroid_id;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    //WE USE UNIFORM DISTRIBUTION TO GET A POINT FROM INPUT POINTS
    uniform_int_distribution<int> p_distribution(0, curve_vector.size()-1);

    first_centroid_id = p_distribution(generator);

    //INSERT THE FIRST CENTROID ID TO CENTROIDS VECTOR
    centroids.push_back(first_centroid_id);
}

void centroids_pick_next_centroid(vector<float>& partial_sums){

    double x;
    int r, i;

    vector<int>::iterator it;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);

    uniform_real_distribution<float> distribution(0.0, partial_sums[partial_sums.size() - 1]);

    do{

        i = 0;
        r = 1;

        //ASSIGN TO X A VALUE IN [0, P(n-t)] USING UNIFORM DISTRIBUTION
        x = distribution(generator);

        //GET INDEX OF THE NEXT CENTROID - WE NEED AN r THAT SATISFIES THE RELATION P(r-1) < x <= P(r)
        while(x > partial_sums[i]){

            r++;
            i++;
        }

        it = find(centroids.begin(), centroids.end(), r);

    }while(it != centroids.end());  //WHILE r ALREADY EXISTS IN CENTROIDS - CALCULATE A NEW r


    //PUSH THE NEW CENTROID ID TO CENTROIDS VECTOR
    centroids.push_back(r);
}

//FIND THE NEAREST CENTROID OF A GIVEN CURVE
int centroids_find_nearest_centroid(vector<int>& centroids_with_same_id, int id){

    //GET CURVE'S COORDINATES
    pair<pair<string, int>, vector<double>> curve = curve_vector_get_curve(id);

    //INITIALIZE THE MIN DISTANCE OF POINT FROM THE CENTROIDS WITH THE MAX VALUE
    double min_distance = numeric_limits<double>::max();
    double current_distance;
    int centroid_with_min_distance; //INDEX OF CENTROID WITH MIN DISTANCE
    int size = centroids_with_same_id.size();

    for(int i = 0; i < size; i++){

        current_distance = calculate_distance(curve_vector[centroids_with_same_id[i]].second, curve.second);

        if(current_distance < min_distance){

            min_distance = current_distance;
            centroid_with_min_distance = i;
        }
    }

    return centroid_with_min_distance;
}

void centroids_duplicates_assign_to_nearest_centroid(vector<pair<vector<int>,int>>& curves_in_range){

    vector<int> centroids_with_same_id;
    int nearest_centroid;
    //FOR ALL THE CENTROIDS
    for(int i = 0; i < curves_in_range.size(); i++){
        //FOR ALL THE POINTS ASSIGNED TO THEM
        for(int k = 0; k < curves_in_range[i].first.size(); k++){

            centroids_with_same_id.push_back(centroids_get_centroid(i));
            //FOR ALL THE OTHER CENTROIDS
            for(int j = i+1; j < curves_in_range.size(); j++){
                //SEARCH IF THERE IS ANOTHER CENTROID THAT HAS THE SAME ID
                search_if_in_range(curves_in_range[j], centroids_with_same_id, curves_in_range[i].first[k], j);

            }   //IF THERE ARE MORE THAN ONE CENTROIDS THAT HAVE THE SAME ID
            if(centroids_with_same_id.size() >= 2){
                //FIND THE NEAREST CENTROID TO THE POINT WITH THIS ID
                nearest_centroid = centroids_find_nearest_centroid(centroids_with_same_id, curves_in_range[i].first[k]);

                for(int centroid = 0; centroid < centroids.size(); centroid++){
                    
                    if(centroid != nearest_centroid){
                        //REMOVE CURVE ID FROM ALL THE CENTROID THAT DO NOT HAVE MINIMUM DISTANCE WITH IT
                        update_curves_in_range(curves_in_range[centroid], curves_in_range[i].first[k]);
                    }
                    
    
                }
            }
            centroids_with_same_id.clear();
        }
    }
 
    //LABEL ALL THE ASSIGNED CURVES AS ASSIGNED
    label_assigned_curves(curves_in_range);

}

//CALCULATE MIN DISTANCE OF A CURVE FROM ALL THE CENTROIDS
float centroids_calculate_min_distance_curve(vector<double>& curve, string assignment){

    float min_distance = numeric_limits<float>::max();
    float current_distance;
    int centroid_with_min_distance;

    for(int i = 0; i < centroids_get_size(); i++){

        if(assignment == "LSH" || assignment == "Hypercube" || assignment == "Classic"){

            //CALCULATE EUCLIDEAN DISTANCE
            current_distance = calculate_distance(curve_vector[centroids[i]].second, curve);
        }
        else if(assignment == "LSH_Frechet"){

            //CALCULATE DISCRETE FRECHET DISTANCE
            current_distance = curve_calculate_dfd(curve_vector[centroids[i]].second, curve);
        }

        if(current_distance < min_distance){

            min_distance = current_distance;
            centroid_with_min_distance = i;
        }


    }

    return min_distance;
}

//CALCULATE MIN DISTANCE OF EACH INPUT CURVE FROM ALL THE CENTROIDS
void centroids_calculate_min_distance_input(vector<float>& curves_min_distances, string assignment){

    float curve_min_distance;

    for(int i = 0; i < curve_vector.size(); i++){

        curve_min_distance = centroids_calculate_min_distance_curve(curve_vector[i].second, assignment);

        //IF CURRENT CURVE IS NOT A CENTROID
        if(curve_min_distance != 0){

            curves_min_distances.push_back(curve_min_distance);
        }
    }
}

//GET EVERY CENTROID'S HASHTABLE BUCKET HASHES - USING FRECHET HASHING
void centroids_get_hashtable_hashes_frechet(G_Frechet g, vector<vector<int>>& hashes, bool is_mean){

    vector<int> hash_vector;
    vector<int> id_vector;

    for(int i = 0; i < centroids.size(); i++){
       
        g.hash(curve_vector_get_curve(centroids[i]), hash_vector, id_vector, true, is_mean, 2);
        hashes.push_back(hash_vector);

        hash_vector.clear();
    }
}

//GET EVERY CENTROID'S HASHTABLE BUCKET HASHES - USING LSH HASHING
void centroids_get_hashtable_hashes_lsh(G_Lsh g, vector<vector<int>>& hashes){

    vector<int> hash_vector;

    for(int i = 0; i < centroids.size(); i++){
        g.hash(curve_vector_get_curve(centroids[i]), hash_vector, 1, 0);
        hashes.push_back(hash_vector);
    }
}

//GET EVERY CENTROID'S HYPERCUBE BUCKET HASH
void centroids_get_hypercube_hashes(G_Hypercube g, vector<int>& hashes){

    unsigned int hash_value;

    for(int i = 0; i < centroids.size(); i++){
        g.hash(curve_vector_get_curve(centroids[i]), hash_value, 1);
        hashes.push_back(hash_value);
    }
}

//CLEAR CENTROIDS DATA STRUCTURE
void centroids_clear(){

    centroids.clear();
}

//PRINT CENTROIDS IDS - USED FOR CHECKING PURPOSES
void centroids_print_data(){

    for(int i = 0; i < centroids.size(); i++){
        cout << centroids[i] << endl;
    }
}

//INITIALIZE THE VECTOR WITH POINT_VECTOR SIZE
//AND ALL IT'S VALUES TO FALSE
void is_assigned_initialize(){

    int size = curve_vector_get_size();
    is_assigned.resize(size);

}

//GET SIZE OF IS_ASSIGNED DATA STRUCTURE (IT MUST BE EQUAL WITH NUM_OF_POINTS)
int is_assigned_get_size(){

    return is_assigned.size();
}

//COUNT NUMBER OF POINTS THAT ARE ASSIGNED
int is_assigned_count_assigned(){
    int assigned_num = 0;

    for(int i = 0; i < is_assigned.size(); i++){
        if(is_assigned[i] == true){
            assigned_num++;
        }
    }

    return assigned_num;
}

//COUNT NUMBER OF POINTS THAT ARE NOT ASSIGNED
int is_assigned_count_unassigned(){
    int assigned_num = 0;

    for(int i = 0; i < is_assigned.size(); i++){
        if(is_assigned[i] == false){
            assigned_num++;
        }
    }

    return assigned_num;
}

//REVERSE ASSIGNMENT - CLUSTERING: IF A POINT IS ASSIGNED TO A CENTROID
//THEN MARK IT AS ASSIGNED
void mark_as_assigned(int index){

    is_assigned[index] = true;
}

//CHECK IF AN INPUT POINT IS ALREADY ASSIGNED TO A CLUSTER
bool already_assigned(int index){

    if(is_assigned[index] == true){

        return true;
    }
    else{

        return false;
    }
}

//PARTITION POINTS IN RANGE OF CENTROID TO ASSIGNED AND UNASSIGNED
void partition_assigned_unassigned(pair<vector<int>,int>& curves_in_range){

    int temp;
    int left = 0;
    int right = curves_in_range.first.size() - 1;
    int partition_pointer = 0;

    while(left < right){

        //WHILE CURVE IN LEFT INDEX IS ALREADY ASSIGNED
        while(already_assigned(curves_in_range.first[left])){
            left++;
            if(left == right){
                curves_in_range.second = left + 1;
                return;
            }
        }
        while(!already_assigned(curves_in_range.first[right])){
            right--;
            if(right == left){
                curves_in_range.second = left + 1;
                return;
            }
        }
        //SWAP THE UNASSIGNED CURVE FROM THE LEFT WITH THE ASSIGNED CURVE FROM THE RIGHT
        temp = curves_in_range.first[left];
        curves_in_range.first[left] = curves_in_range.first[right];
        curves_in_range.first[right] = temp;
        left++;
        right--;

    }

    //KEEP TRACK OF THE FIRST INDEX WITH UNASSIGNED POINT
    curves_in_range.second = left;
}

//CHECK IF A POINT IS ASSIGNED OR NOT
bool is_assigned_get_value(int index){

    return is_assigned[index];
}

//LABEL ALL THE POINTS THAT HAVE BEEN ASSIGNED TO A CLUSTER AS ASSIGNED
void label_assigned_curves(vector<pair<vector<int>,int>>& curves_in_range){
    //FOR ALL THE CLUSTERS
    for(int i = 0; i < curves_in_range.size(); i++){
        //FOR ALL THE POINTS ASSIGNED TO THEM
        for(int j = 0; j < curves_in_range[i].first.size(); j++){
            mark_as_assigned(curves_in_range[i].first[j]);
            //KEEP TRACK OF THE FIRST INDEX OF EVERY CLUSTER THAT HAS NOT BEEN ASSIGNED
            curves_in_range[i].second++;
        }
    }
}

//PRINT POINTS DATA REGARDING TO WHETHER THEY ARE ASSIGNED OR NOT
//USED FOR CHECKING PURPOSES
void assigned_print_assigned(){
    int count = 0;
    for(int i = 0; i < is_assigned.size(); i++){
        if(is_assigned[i] == true){
            cout << i+1 << " point is assigned " << endl;
            count++;
        }
    }
    cout << "number of assigned points is " << count << endl;
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

vector<vector<double>> get_dfd_array(const vector<double>& curve1, const vector<double>& curve2){
    vector<vector<double>> DFD;     //ARRAY THAT KEEPS DISCRETE FRECHET DISTANCE FOR EVERY POSSIBLE TRAVERSAL
    vector<double> point1, point2;  //CURVE POINTS
    int num_of_curve1_points = curve1.size();
    int num_of_curve2_points= curve2.size();

    DFD.resize(num_of_curve1_points, vector<double>(num_of_curve2_points));
    
    //FOR ALL THE VERTICES OF CURVE 1
    for(int i = 0; i < num_of_curve1_points; i++){
        //PAIR THEIR X AND Y VALUES
        point1.push_back(i);
        point1.push_back(curve1[i]);

        //FOR ALL THE VERTICES IN CURVE 2
        for(int j = 0; j < num_of_curve2_points; j++){
            //DO THE SAME
            point2.push_back(j);
            point2.push_back(curve2[j]);

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

    return DFD;
}

//CALCULATE DISCRETE FRECHET DISTANCE BETWEEN TWO CURVES
double curve_calculate_dfd(const vector<double>& curve1, const vector<double>& curve2){

    vector<vector<double>> DFD;     //ARRAY THAT KEEPS DISCRETE FRECHET DISTANCE FOR EVERY POSSIBLE TRAVERSAL
    vector<double> point1, point2;  //CURVE POINTS
    int num_of_curve1_points = curve1.size();
    int num_of_curve2_points= curve2.size();

    DFD.resize(num_of_curve1_points, vector<double>(num_of_curve2_points));
    
    //FOR ALL THE VERTICES OF CURVE 1
    for(int i = 0; i < num_of_curve1_points; i++){
        //PAIR THEIR X AND Y VALUES
        point1.push_back(i);
        point1.push_back(curve1[i]);

        //FOR ALL THE VERTICES IN CURVE 2
        for(int j = 0; j < num_of_curve2_points; j++){
            //DO THE SAME
            point2.push_back(j);
            point2.push_back(curve2[j]);

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

    return DFD[num_of_curve1_points-1][num_of_curve2_points-1];

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

//FIND OPTIMAL TRAVERSAL BETWEEN TWO CURVES
vector<pair<int,int>> find_optimal_traversal(const vector<double>& curve1, const vector<double>& curve2){

    vector<vector<double>> DFD = get_dfd_array(curve1, curve2);

    int min_index;
    vector<double> neighbor_cells;
    vector<pair<int,int>> optimal_traversal;
    int i = curve1.size()-1;
    int j = curve2.size()-1;

    pair<int,int> last_step = make_pair(i, j);

    optimal_traversal.push_back(last_step);

    while(i > 0 && j > 0){
        
        neighbor_cells.push_back(DFD[i-1][j]);
        neighbor_cells.push_back(DFD[i][j-1]);
        neighbor_cells.push_back(DFD[i-1][j-1]);

        min_index = std::min_element(neighbor_cells.begin(), neighbor_cells.end()) - neighbor_cells.begin();
        
        if(min_index == 0){
            optimal_traversal.push_back(make_pair(--i, j));
        }
        else if(min_index == 1){
            optimal_traversal.push_back(make_pair(i,--j));
        }
        else{
            optimal_traversal.push_back(make_pair(--i,--j));
        }

        neighbor_cells.clear();

    }
    
    if(i == 0 && j > 0){
        
        while(j > 0){
            optimal_traversal.push_back(make_pair(i, --j));
        }
        
    }
    else if(i > 0 && j == 0){
        
        while(i > 0){
           optimal_traversal.push_back(make_pair(--i, j)); 
        }
    }

    return optimal_traversal;
}

//GET MEAN CURVE BETWEEN TWO CURVES - FRECHET
pair<pair<string, int>, vector<double>> get_mean_curve(const vector<double>& curve1, const vector<double>& curve2){

    pair<pair<string, int>, vector<double>> mean_curve;
    double mean_value;
    int count = 0;

    vector<pair<int,int>> optimal_traversal = find_optimal_traversal(curve1, curve2);
    mean_curve.second.resize(optimal_traversal.size());
    
    for(int i = optimal_traversal.size()-1; i >=0; i--){

        mean_value = (curve1[optimal_traversal[i].first] + curve2[optimal_traversal[i].second])/2;
        mean_curve.second[count] = mean_value;

        count++;
    }

    return mean_curve;
}

//GET MEAN CURVE BETWEEN TWO CURVES AS A VECTOR
pair<pair<string, int>, vector<double>> get_mean_curve_vector(vector<double> vector_of_sums, int num_of_vectors, int& last_known_id){

    double mean_value;
    pair<pair<string, int>, vector<double>> mean_curve;

    mean_curve.first.second = ++last_known_id;
    for(int i = 0; i < vector_of_sums.size(); i++){
        mean_value = vector_of_sums[i]/num_of_vectors;
        mean_curve.second.push_back(mean_value);
    }

    return mean_curve;
}


//FILTER A CURVE TO HAVE A PREDEFINED SIZE
void filter(vector<double>& curve, double epsilon, int max_length){

    while(curve.size() > max_length){

        for(int i = 1; i < curve.size() - 1; i++){

            if(curve.size() > max_length){
                if((abs(curve[i-1] - curve[i]) <= epsilon) && (abs(curve[i] - curve[i+1]) <= epsilon)){
                    curve.erase(curve.begin() + i);
                }
            }
            else{
                break;
            }
            
        }

        if(curve.size() > max_length){
            epsilon = epsilon*2;
        }
        else{

            break;
        }
    }

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
double calculate_distance(const vector<double>& curve1, const vector<double>& curve2, int k)
{
    double distance = 0.0;
    double sum = 0;

    for (int i=0 ; i < curve1.size() ; i++) {
        sum+= pow(abs(curve1[i]-curve2[i]), k);
    }
    distance = pow(sum, 1.0/(double)k);

    return distance;
}

float calculate_partial_sums(vector<float>& min_distances, vector<float>& partial_sums){

    float sums = 0;
    float last_partial_sum;

    for (int i = 0; i < min_distances.size(); i++){

        sums = sums + pow(min_distances[i], 2);
        partial_sums.push_back(sums);

        if(i == min_distances.size()){
            last_partial_sum = sums;
        }
    }

    return last_partial_sum;
}

//CHECK IF A GIVEN ID ALREADY EXISTS IN ID VECTOR
bool already_exists(vector<int>& ids, int id){

    if(find(ids.begin(), ids.end(), id) != ids.end() )
        return true;
    else
        return false;
}

//REMOVE POINT ID FROM ALL THE CLUSTERS OF CENTROID'S THAT DO NOT HAVE MINIMUM DISTANCE WITH IT
void update_curves_in_range(pair<vector<int>,int>& curves_in_range, int id){

    curves_in_range.first.erase(std::remove(curves_in_range.first.begin(),
    curves_in_range.first.end(), id), curves_in_range.first.end());
}

//SEARCH IF A CURVE ID IS IN ASSIGNED TO A CENTROID - CLUSTERING
void search_if_in_range(pair<vector<int>,int>& curves_in_range, vector<int>& centroid, int id, int num){
    //FOR ALL THE CURVES IN RANGE THAT HAVE NOT ALREADY BEEN ASSIGNED
    for(int i = 0; i < curves_in_range.first.size(); i++){
        //IF GIVEN POINT ID IS ASSIGNED TO THIS CENTROID
        if(curves_in_range.first[i] == id){
            //ADD THE CENTROID TO CENTROIDS WITH SAME ID
            centroid.push_back(centroids_get_centroid(num));

        }
    }
}

//UPDATE CLUSTER TABLE WITH THE IDS OF THE POINTS IN EACH CLUSTER
void get_cluster_table(vector<pair<vector<int>,int>>& points_in_range, vector<vector<int>>& cluster_table){

    cluster_table.clear();

    for(int i = 0; i < points_in_range.size(); i++){

        cluster_table.push_back(points_in_range[i].first);
    }
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