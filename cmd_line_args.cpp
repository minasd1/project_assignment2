#include <string>
#include <cstdio>
#include <iostream>

using namespace std;

//Initializes program's variables with command line arguments
int read_cmd_args(int argc, char** argv, string& input_file, string& query_file,
                      int& k_lsh, int& k_cube, int& l, string& output_file, int& n,
                      float& r, int& m, int& probes, string& config_file, string& assignment,
                      bool& complete_flag, string& algorithm, string& metric, double& delta, 
                      string& update, bool& datapath_given, bool& query_given, bool& output_given, 
                      bool& algorithm_given, bool& metric_given, bool& silhouette_flag)
{
    if (argc < 3) {
        cerr << "Arguments are missing!" << endl;
        return -1;
    }

    int i;
    bool k_lsh_flag, k_cube_flag, l_flag, n_flag, r_flag, probes_flag, delta_flag, m_flag, c_flag, 
         method_flag, assignment_flag, update_flag;

    //Flags for given arguments (false for missing args)
    datapath_given= false;
    k_lsh_flag= false;
    k_cube_flag= false;
    l_flag= false;
    output_given= false;
    query_given= false;
    n_flag= false;
    r_flag= false;
    probes_flag= false;
    m_flag= false;
    c_flag= false;
    assignment_flag= false;
    complete_flag= false;
    algorithm_given= false;
    metric_given= false;
    delta_flag= false;
    update_flag= false;
    silhouette_flag= false;


    for (i=1; i < argc ; i+=2) { //For every other argument
        if((string)argv[i] == "-i"){
            datapath_given= true;
            input_file= argv[i+1];
        }
        else if((string)argv[i] == "-q" ) {
            query_given= true;
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
            output_given= true;
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
            if ((string) argv[i+1] != "Classic" && (string) argv[i+1] != "classic" 
              && (string) argv[i+1] != "LSH"    && (string) argv[i+1] != "Lsh"
              && (string) argv[i+1] != "lsh"  && (string) argv[i+1] != "Hypercube"
              && (string) argv[i+1] != "hypercube" && (string) argv[i+1] != "Frechet"
              && (string) argv[i+1] != "frechet")
            {
                cerr << "Invalid algorithm given! {LSH, Hypercube, Frechet}" << endl;
                return -1;
            }
            algorithm_given= true;
            algorithm= argv[i+1];
        }
        else if((string)argv[i] == "-metric") {
            metric_given= true;
            metric= argv[i+1];
        }
        else if((string)argv[i] == "-delta") {
            delta_flag= true;
            delta= stod(argv[i+1]);
        }
        else if ((string)argv[i] == "-c") {
            c_flag= true;
            config_file= argv[i+1];
        }
        else if (((string)argv[i] == "-silhouette"  || (string)argv[i] == "-complete") && (i != argc-1   && i != argc-2) ) {
            cerr << "Please give the \"-silhouette\"  and \"-complete\" args in the end" << endl;
            return -1;
        }
        else if ((string) argv[i] == "-silhouette") {
            continue;
        }
        else if ((string) argv[i] == "-complete") {
            continue;
        }

        else if ((string)argv[i] == "-assignment") {
            if ((string)argv[i+1] == "Classic" || (string)argv[i+1] == "classic"
                || (string)argv[i+1] == "CLASSIC"){
                    assignment_flag= true;
                    assignment= "Classic";
            }
            else if ((string)argv[i+1] == "Lsh" || (string)argv[i+1] == "lsh"
                || (string)argv[i+1] == "LSH") {
                    assignment_flag= true;
                    assignment= "LSH";
            }
            else if ((string)argv[i+1] == "Hypercube" || (string)argv[i+1] == "hypercube"
                || (string)argv[i+1] == "HYPERCUBE") {
                    assignment_flag= true;
                    assignment= "Hypercube";
            }
            else if ((string)argv[i+1] == "LSH_Frechet" || (string)argv[i+1] == "lsh_frechet"
                || (string)argv[i+1] == "Lsh_Frechet"   || (string)argv[i+1] == "Lsh_frechet"
                || (string)argv[i+1] == "lsh_Frechet"   || (string)argv[i+1] == "LSH_FRECHET") {
                    assignment_flag= true;
                    assignment= "LSH_Frechet";
            }
            else {
                cerr << "Unknown assignment method \"" << argv[i+1] << "\". " << "Please run the program with one of the acceptable methods:" << endl;
                cerr << "(Classic, LSH, Hypercube, LSH_Frechet)" << endl;
                return -1;
            }
        }
        else if((string)argv[i] == "-update" && ((string) argv[i+1] != "Mean_Frechet" && (string) argv[i+1] != "Mean_Vector" )) {
            cerr << "Please give the update method using one word (Mean_Frechet or Mean_Vector)" << endl;
            return -1;
        }
        else if ((string) argv[i] == "-update" ) {
            update_flag= true;
            update= argv[i+1];
        }
        else {
            cerr << "Wrong input arguent: " << argv[i] << endl;
            return -1;
        }
    }

    /*if (!algorithm_flag && (string)argv[0] == "./search") {
        cerr << "Algorithm not defined!" << endl;
        return -1;
    }*/
    if ((string)argv[argc-1] == "-silhouette"  || (string) argv[argc-2] == "-silhouette") {
        silhouette_flag= true;
    }
    if ((string)argv[argc-1] == "-complete"  || (string) argv[argc-2] == "-complete") {
        complete_flag= true;
    }
    if (!assignment_flag && (string)argv[0] == "./cluster"){
        cerr << "Assignment method not defined!" << endl;
        return -1;
    }
    if (!update_flag && (string)argv[0] == "./cluster"){
        cerr << "Update method not defined!" << endl;
        return -1;
    }
    if (!c_flag && (string) argv[0] == "./cluster") {
        cerr << "Configuration file is missing!" << endl;
        return -1;
    }    
    if ((string)argv[0] == "search") {
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
            if (!m_flag )
                m= 10;
            if (!probes_flag )
                probes= 2;
            return 0;
    }
    else if ((string)argv[0] == "./cluster"){
            metric_given= true;
            query_given= true;
            algorithm_given= true;
            //INITIALIZE LSH AND HYPERCUBE ARGUMENTS WITH THE DEFAULT VALUES
            l= 3;
            k_lsh= 4;
            k_cube= 3;
            m= 10;
            probes= 2;
            n= 1;
            
    }
    return 0;
}
