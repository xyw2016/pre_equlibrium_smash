//     g++ convert_to_binary_SMASH.cpp -lz -o convert_to_binary_SMASH.e

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>


using namespace std;

int main(int argc, char *argv[]) {
    string input_filename;
    string output_filename;

    if (argc < 3) {
        cout << "Usage: " << argv[0] << " inputfile outputfile" << endl;
        exit(0);
    }else if(argc ==3){
        input_filename = argv[1];
        output_filename = argv[2];    
    }
    else{
        cout<<" Error: wrong number of input file "<<endl;
        exit(0);
    }

    ifstream SMASH_file(input_filename.c_str());
    if (!SMASH_file.is_open()) {
        cout << "can not open " << input_filename << endl;
        cout << "exit" << endl;
        exit(1);
    }

    ofstream outbin;
    //outbin.open(output_filename.c_str(), ios::app| ios::out | ios::binary);
    outbin.open(output_filename.c_str(), ios::app| ios::out);

    string temp_string;
    double dummy;
    string str_dummy;
    int pdg, charge, ncoll, baryon_number;
    double mass, t, x, y, z, E, px, py, pz;
    

    // start to read in SMASH outputs
    // skip the header
    for (int i = 0; i < 3; i++) {
        getline(SMASH_file, temp_string);
    }



    
    std::istringstream iss;
    
    int ievent=0;
    
    while(SMASH_file.good()){
        int nline=0;
        
        while(getline(SMASH_file,temp_string))
        {
           iss.str(temp_string);
           iss >> str_dummy >> str_dummy >> str_dummy >> str_dummy;
           if (str_dummy == "in"){
            iss.clear();
            float n_particles[ ]={static_cast<float>(888888)};
            outbin.write((char*) &(n_particles[0]), sizeof(float));
            continue;
           }
           else if(str_dummy == "end"){
            iss.clear();
            std::cout<<"Finish read event " << ievent <<" "<<"particle number: "<<nline<<std::endl;
            ievent++;
            float n_particles[ ]={static_cast<float>(999999)};
            outbin.write((char*) &(n_particles[0]), sizeof(float));
            break;
           }
           else{
            iss.seekg(0);
           }
           

           

            iss >>  t  >> x >> y >> z 
                >>  mass >> E >> px >> py
                >>  pz   >> pdg >> str_dummy >> charge >> ncoll
                >> str_dummy>> str_dummy >> str_dummy >> str_dummy >> str_dummy>> str_dummy>> str_dummy>> baryon_number;
            iss.clear();

            float particle_array [] = {static_cast<float>(t), 
                                       static_cast<float>(x),
                                       static_cast<float>(y),
                                       static_cast<float>(z),
                                       static_cast<float>(mass),
                                       static_cast<float>(E),
                                       static_cast<float>(px),
                                       static_cast<float>(py),
                                       static_cast<float>(pz),
                                       static_cast<float>(pdg),
                                       static_cast<float>(charge),
                                       static_cast<float>(ncoll),
                                       static_cast<float>(baryon_number)};
            for (int ii = 0; ii < 13; ii++) {
                outbin.write((char*) &(particle_array[ii]), sizeof(float));
            }
           
            nline++;
        }    
        if (SMASH_file.eof()) break;
    
    }

    outbin.close();
    SMASH_file.close();   
    return 0;
}
