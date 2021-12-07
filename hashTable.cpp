#include "hashTable.h"

void hashTable_initialization(int num_of_hashTables, int num_of_buckets){

    //INITIALIZE L HASHTABLES WITH HASHTABLESIZE BUCKETS AND 1 CURVE (WITH ID 0) IN EACH BUCKET
    HashTables.resize(num_of_hashTables, vector<vector<int> >(num_of_buckets, vector<int>(0)));
    //IT CAN ALSO BE INITIALIZED AS ZERO
}

//PUSH A CURVE'S ID VALUES TO A BUCKET OF EACH OF THE HASHTABLES USING G FUNCTION
//EACH BUCKET IS INDICATED BY THE INT VALUES PRODUCED BY FUNCTION G
void hashTable_push_back(vector<int>& g, int key_val){

    for(int i = 0; i < HashTables.size(); i++){
        HashTables[i][g[i]].push_back(key_val);
        //HashTables[i][g[i]].end() = key_val;
    }
}

//PUSH A CURVE'S ID VALUE TO A BUCKET OF THE HASHTABLES THAT THE HASH_GRID CORRESPONDS
//USING G FUNCTION
//THE BUCKET IS INDICATED BY THE INT VALUES PRODUCED BY FUNCTION G
// void hashTable_push_back_frechet(int& g, int key_val, int grid){

//         HashTables[grid-1][g].push_back(key_val);
//         //HashTables[i][g[i]].end() = key_val;
// }

//GET KEY VALUES OF ALL THE CURVES OF A BUCKET OF EACH HASH TABLE
//EACH BUCKET IS INDICATED BY THE INT VALUES PRODUCED BY FUNCTION G
vector<int> hashTable_get_curves_in_buckets(vector<int>& g){

    vector<int> points_in_buckets;

    for(int i = 0; i < HashTables.size(); i++){
        for(int j = 0; j < HashTables[i][g[i]].size(); j++){
            //INSERT THE NUMBER OF HASHTABLE OF THE CURVE (USE TO GET APPROPRIATE ID VALUE)
            points_in_buckets.push_back(i);
            //PUSH THE ID OF THE CURVE TO POINTS IN HASHTABLE BUCKETS
            points_in_buckets.push_back(HashTables[i][g[i]][j]);
        }
    }

    return points_in_buckets;
}

//GET NUMBER OF HASHTABLES - L
int hashTable_get_num_of_htables(){

    return HashTables.size();
}

//GET NUMBER OF CURVES IN A HASHTABLE BUCKET
int hashTable_get_bucket_size(int htable_num, int bucket){

    return HashTables[htable_num][bucket].size();
}

//GET A SPECIFIC CURVE FROM A HASHTABLE BUCKET
int hashTable_get_curve(int htable_num, int bucket, int place){

    return HashTables[htable_num][bucket][place];
}

//PRINT HASHTABLES DATA - USED FOR CHECKING PURPOSES
void hashTable_print_data(){
    for (int i = 0; i < HashTables.size(); i++) {
        for (int j = 0; j < HashTables[i].size(); j++){
            for(int k = 0; k < HashTables[i][j].size(); k++){
                cout<< HashTables[i][j][k] << " ";
            }
        }
        cout << endl;
    }
    cout << "Hashtable size is: " << HashTables.size() << endl;
    cout << "buckets of each hashtable are: " << HashTables[0].size() << endl;
    cout << "elements in each bucket are: " << HashTables[0][0].size() << endl;
}

//PRINTS HASHTABLE - USED FOR CHECKING PORPUSES
void hashTable_print(){
    int sum= 0;
    for (int i = 0; i < HashTables.size(); i++) {
        cout << "HASH TABLE NUMBER " << i+1 << endl;
        cout << "(size: " << HashTables[i].size() << ")" << endl;
        cout << "---------------------" << endl << endl;
        for (int j = 0; j < HashTables[i].size(); j++){
            cout << "Bucket " << j << endl;
            cout << "(size: " << HashTables[i][j].size() << ")" << endl;
            sum+= HashTables[i][j].size();
            cout << "Contains (curve ids): ";
            for(int k = 0; k < HashTables[i][j].size(); k++){
                cout<< HashTables[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << "Total number of hash tables: " << HashTables.size() << endl;
    cout << "Each hash table has: " << HashTables[0].size() << " buckets" << endl;
    cout << "Mean bucket size is: " << sum/(HashTables.size() * HashTables[0].size()) << endl;
}
