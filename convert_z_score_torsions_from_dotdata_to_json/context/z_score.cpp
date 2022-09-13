#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <cstring>
#include <sstream>      

struct Bin { 
    float lower_phi; 
    float higher_phi;
    
    float lower_psi; 
    float higher_psi;

    int count; 
};

struct Data {  
    std::string linkage;
    std::vector<Bin> bins;
    std::string file_name;

    float count_mean, count_stddev;

    int get_bin_size();
    void print_bin(const Bin& bin);

    void load_data_from_file();

    float calculate_z_score(float phi, float psi);

    private:
        float count(float phi, float psi );
};


float Data::count(float phi, float psi) { 

    std::cout << bins.size() << std::endl;

    // If efficiency is an issue can get index from (rounded_bin_angle + 180)/binspacing but a loop is more readable.
    for(const auto& bin: bins) { 
            if ((bin.lower_phi <= phi)) {
                if (phi < bin.higher_phi) {
                    if (bin.lower_psi <= psi){
                        if  (psi < bin.higher_psi) { 
                            return bin.count;
                        }
                    }
                } 
            }    
        }
    
    return 0.0f;
}

float Data::calculate_z_score(const float phi, const float psi) {

    float count = Data::count(phi, psi);
    float z_score = (count - count_mean) / count_stddev;
    return z_score;
}

void Data::load_data_from_file() {
    std::ifstream input_file(file_name);
    std::string entry;

    if(!input_file.is_open()) {
        std::runtime_error("Could not open file!");
        exit(EXIT_FAILURE);
    }

    while (std::getline(input_file,entry)) { 
        float column_1, column_2, column_3, column_4, column_5; 
        
        std::istringstream string_stream(entry);

        if (string_stream >> column_1 >> column_2) { 
            if (string_stream >> column_3 >> column_4 >> column_5) { 
                Bin new_bin; 
                new_bin.lower_phi = column_1;
                new_bin.higher_phi = column_2;
                new_bin.lower_psi = column_3;
                new_bin.higher_psi = column_4;
                new_bin.count = column_5;

                bins.push_back(new_bin);
            }
            else { 
                count_mean = column_1; 
                count_stddev = column_2;
            }
        }
    }
}

int Data::get_bin_size() { 

    if (bins.size() > 0) { 
        int phi_bin_size = bins[0].higher_phi-bins[0].lower_phi;
        int psi_bin_size = bins[0].higher_psi-bins[0].lower_psi;

        if (phi_bin_size != psi_bin_size) { 
            throw std::runtime_error("Bin sizes from data file are not consistent.");
        }
        return phi_bin_size;
    } else { 
        throw std::runtime_error("Data files are not loaded.");
    }
}

int main() { 
    
    std::cout << "BEGINNING SCRIPT" << std::endl;

    Data dat; 
    dat.file_name = "./data/bins/NAG-1-4-NAG.data"; 
    dat.linkage = "NAG-NAG";
    dat.load_data_from_file();


    // To calculate a z_score pass in phi and psi into this function calculate_z_score(phi,psi)
    float z_score = dat.calculate_z_score(-79.9, -129.9);
    std::cout << z_score << std::endl;
    

    // Potential implementation: 
    // Input PDB file
    // Call Privateer functions to analyse glycans
    // Collect a list of linkages with all the torsions
    // sort by linkage 
    // load the correct linkages data from the file
    // establish the z scores for each and store
    // calculate average z score simply by taking sum(z_score)/number of z_scores

    // To obtain a measure of quality of protein, we need to get the z scores for every protein in the dataset



    return 0; 
}