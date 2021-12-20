#include "read_file_data.h"


double file_get_max_value(fstream& file, string file_name){

    int start, finish;
    int num_of_values = 0;
    bool first_iteration = true;
    string line, token;
    double max_value = 0.0;
    double max_coordinate_value;

    open_file(&file, file_name, fstream::in);

    finish = 0;

    //FOR EVERY QUERY CURVE
    while(getline(file, line)){                       //READ QUERY FILE LINE BY LINE

        start = 0;

        while(start < line.size()){                         //TOKENIZE EVERY LINE IN IT'S SEPERATED STRINGS
            finish = line.find_first_of('\t', start);

            if(finish == string::npos){

                finish = line.size();
            }

            if(start < line.size() - 1){
                token = line.substr(start, finish - start);
                if(start == 0){

                }
                else{
                    //KEEP TRACK OF THE MAX VALUE RECORDED IN FILE
                    if((stod(token)) > max_value){

                        max_value = stod(token);
                    }
        
                }
                
            }

            if(first_iteration == true && start != 0){

                num_of_values++;
            }
            start = finish + 1;
        
        }

        first_iteration = false;
    }

    max_coordinate_value = max(max_value, (double)num_of_values);

    file.clear();
    file.seekg(0, file.beg);

    return max_coordinate_value;

}