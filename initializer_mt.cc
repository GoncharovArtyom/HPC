//
// Created by artyom on 29.03.18.
//
#include "parallel_mt.hh"
#include <iostream>
#include <fstream>

using namespace std;
using namespace autoreg;

int main(int argc, char ** argv){
    if (not(argc == 3 or argc == 4)){
        cerr<<"Incorrect parameter. Usage: initialize_mt [n_configs] [file_name] [seed]"<<endl;
    }

    int n_configs = stoi(argv[1]);
    string file_name = argv[2];
    int seed = 0;
    if (argc == 4){
        seed = stoi(argv[3]);
    }

    ofstream file(file_name);
    parallel_mt_seq<> initializer(seed);

    for (int i=0; i<n_configs; ++i){
        mt_config configuration = initializer();
        file << configuration;
    }

    cout<<"Done."<<endl;

}

