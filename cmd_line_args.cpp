#include <string>
#include <cstdio>
#include <iostream>

using namespace std;

//Initializes program's variables with command line arguments
int read_cmd_args(int argc, char** argv, string& input_file, string& query_file,
                      int& k_lsh, int& k_cube, int& l, string& output_file, int& n,
                      float& r, int& m, int& probes, string& config_file, 
                      bool& complete_flag, string& algorithm, string& metric, double& delta)
{
    int i;
    bool input_file_flag, k_lsh_flag, k_cube_flag, l_flag, output_file_flag, query_file_flag, 
         n_flag, r_flag, probes_flag, algorithm_flag, metric_flag, delta_flag, m_flag, c_flag, 
         method_flag;

    //Flags for given arguments (false for missing args)
    input_file_flag= false;
    k_lsh_flag= false;
    k_cube_flag= false;
    l_flag= false;
    output_file_flag= false;
    query_file_flag= false;
    n_flag= false;
    r_flag= false;
    probes_flag= false;
    m_flag= false;
    c_flag= false;
    method_flag= false;
    complete_flag= false;
    algorithm_flag= false;
    metric_flag= false;
    delta_flag= false;

    for (i=1; i < argc ; i+=2) { //For every other argument
        if((string)argv[i] == "-i"){
            input_file_flag= true;
            input_file= argv[i+1];
        }
        else if((string)argv[i] == "-q" ) {
            query_file_flag= true;
            query_file= argv[i+1];
        }
        else if((string)argv[i] == "-k") {
            k_lsh_flag= true;
            k_cube_flag= true;
            k_lsh= stoi(argv[i+1]);
            k_cube= stoi(argv[i+1]);
        }
        else if((string)argv[i] == "-L") {
            l_flag= true;
            l= stoi(argv[i+1]);
        }
        else if((string)argv[i] == "-o") {
            output_file_flag= true;
            output_file= argv[i+1];
        }
        else if((string)argv[i] == "-N") {
            n_flag= true;
            n= stoi(argv[i+1]);
        }
        else if((string)argv[i] == "-R") {
            r_flag= true;
            r= stof(argv[i+1]);
        }
        else if((string)argv[i] == "-M") {
            m_flag= true;
            m= stoi(argv[i+1]);
        }
        else if((string)argv[i] == "-probes") {
            probes_flag= true;
            probes= stoi(argv[i+1]);
        }
        else if((string)argv[i] == "-algorithm") {
            algorithm_flag= true;
            algorithm= argv[i+1];
        }
        else if((string)argv[i] == "-metric") {
            metric_flag= true;
            metric= argv[i+1];
        }
        else if((string)argv[i] == "-delta") {
            delta_flag= true;
            delta= stod(argv[i+1]);
        }
        /*else if ((string)argv[i] == "-c" && (string)argv[0] == "./cluster") {
            c_flag= true;
            config_file= argv[i+1];
        }
        else if ((string)argv[i] == "-complete" && (string)argv[0] == "./cluster") {
            complete_flag= true;
        }
        else if ((string)argv[i] == "-m" && (string)argv[0] == "./cluster") {
            if ((string)argv[i+1] == "Classic" || (string)argv[i+1] == "classic"
                || (string)argv[i+1] == "CLASSIC"){
                    method_flag= true;
                    method= "classic";
            }
            else if ((string)argv[i+1] == "Lsh" || (string)argv[i+1] == "lsh"
                || (string)argv[i+1] == "LSH") {
                    method_flag= true;
                    method= "lsh";
            }
            else if ((string)argv[i+1] == "Hypercube" || (string)argv[i+1] == "hypercube"
                || (string)argv[i+1] == "HYPERCUBE") {
                    method_flag= true;
                    method= "hypercube";
            }
            else {
                cerr << "Unknown method \"" << argv[i+1] << "\". " << "Please run the program with one of the acceptable methods:" << endl;
                cerr << "(Classic, Lsh, Hypercube)" << endl;
                return -1;
            }
        }*/
        else {
            cerr << "Wrong input arguent: " << argv[i] << endl;
            return -1;
        }
    }

    if (!algorithm_flag && (string)argv[0] == "./search") {
        cerr << "Algorithm not defined!" << endl;
    }
    
    if ((string)argv[0] == "search") {
        if (input_file_flag && query_file_flag) {
            //Initialize missing arguments with default values
            if (!k_lsh_flag)
                k_lsh= 4;
            if  (!k_cube_flag)
                k_cube= 14;
            if (!l_flag)
                l= 5;
            if(!n_flag)
                n= 1;
            if(!r_flag)
                r= 10000.0;
            if (!m_flag && (string)argv[0]=="./cube")
                m= 10;
            if (!probes_flag && (string)argv[0]=="./cube")
                probes= 2;
            if (!output_file_flag)
                output_file= "output_file";
            if (!input_file_flag)
                input_file= "NULL";
            if (!query_file_flag)
                query_file= "NULL";
            return 0;
        }
    }
    /*else if ((string)argv[0] == "./cluster"){
        if (!i_flag || !method_flag || !c_flag) {
            cerr << "At least one of the following is missing: ";
            cerr << "(Input file, Configuration file, Method)" << endl;
            cerr << "Please run the program as below:\n" << endl;
            cerr << "./cluster -i <input_file> -c <configuration_file> -m <Classic OR LSH or Hypercube>" << endl;
            return -1;
        }
        else {
            //INITIALIZE LSH AND HYPERCUBE ARGUMENTS WITH THE DEFAULT VALUES
            l= 3;
            k_lsh= 4;
            k_cube= 3;
            m= 10;
            probes= 2;
            n= 1;
            if (!o_flag) {
                output_file= "output_file";
            }
        }
    }*/
    return 0;
}
