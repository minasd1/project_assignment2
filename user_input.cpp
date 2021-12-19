#include "user_input.h"


//READ INPUT FROM USER AFTER FIRST EXECUTION
//ACT ACCORDINGLY IF USER WISHES TO CONTINUE OR NOT
//INFORM USER IF INPUT IS WRONG
void read_user_input(string& query_file, int* continue_execution){
    
    string user_input;
    string query_file_name;

    int valid_input = 0;

    cout << endl;
    cout << "Do you wish to continue execution for another query file?" << endl;
    cout << "Press \"N\" (and then Enter) if not , or press \"Y\" (and then Enter) if yes: ";

    while(valid_input == 0){

        
        cin >> user_input;
        getchar();
    
        if(user_input.compare("N") == 0){
            
            valid_input = 1;
            *continue_execution = 0;                    //STOP PROGRAM EXECUTION

        }
        else if(user_input == "Y"){

            cout << "Enter query file name: ";
            cin >> query_file_name;
            getchar();
            
            // if(query_file_name.compare("query_small_id") == 0){
            //     cout << "it is the same!" << endl;
            // }

            valid_input = 1;

            query_file = query_file_name;
            
            *continue_execution = 1;                    //CONTINUE PROGRAM EXECUTION WITH ANOTHER QUERY FILE

        }
        else{

            cout << "Wrong input, please insert a valid command: ";

            valid_input = 0;
        }
    }

}

void read_path(string& dataset_path, string& queryset_path, string& output_file, string& algorithm, 
               string& metric, bool datapath_given, bool query_given, bool output_given, 
               bool algorithm_given, bool metric_given)
{
    if (!datapath_given) {
        cout << endl << "Please insert the path to the dataset: ";
        cin >> dataset_path;
        getchar();
    }
    if (!query_given) {
        cout << endl << "Please insert the path to the queryset: ";
        cin >> queryset_path;
        getchar();
    }
    if (!output_given) {
        cout << endl << "Please insert the output file: ";
        cin >> output_file;
        getchar();
    }
    if (!algorithm_given) {
        cout << endl << "Please insert the algorithm (LSH, Hypercube, Frechet): ";
        cin >> algorithm;
        getchar();
        while (algorithm != "LSH" && algorithm != "Hypercube" && algorithm != "Frechet") {
            cout << "Please insert a valid algorithm (LSH, Hypercube, Frechet): ";
            cin >> algorithm;
            getchar();
        }
    }
    if (!metric_given && (algorithm == "Frechet" || algorithm == "frechet")) {
        cout << endl << "Please insert the metric (discrete, continuous): ";
        cin >> metric;
        getchar();
        while (metric != "discrete" && metric != "continuous") {
            cout << "Please insert a valid metric (discrete, continuous): ";
            cin >> metric;
            getchar();
        }
    }
}