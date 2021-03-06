/*
 * hybridzone_simulator is a program in the project assortative_crows
 *
 * Copyright (C) 2020 Dirk Metzler, http://evol.bio.lmu.de/_statgen/
 *
 * This file is part of hybridzone_simulator.
 *
 * hybridzone_simulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "hybridzone_genotypes.hh"

int main(int argc, char * argv[]) {
    if(argc<5) {
	cerr << "missing line arguments\n";
	exit(1);
    }
    //double assort_strength=atof(argv[1]), assort_threshold=0.0;
    string outputname(argv[4]);
    initialize_w(argv[1],argv[2]);
    // initialize_w_speed(200,0.05); 
    //initialize_w_neutral();
    
    // initialize_disp_kernel(5);
    initialize_disp_kernel_mix2norm_siefke();
    initialize_disp_and_dxxp();
    initialize_freq(atof(argv[3]));
    // initialize_freq_speed();
    // initialize_freq_speed_eq();
    
    print_allele_frequencies();
    print_w();
    cout << "Dispersal probabilities:\n";
    for(unsigned i=0; i<21; ++i) cout << disp_kernel[i] << ',';
    cout << endl;

    ofstream of(outputname+"simres.csv");
    // for(unsigned t=0; t<2000; ++t) {
    for(unsigned t=0; t<=6000; ++t) {
	cout << "t = " << t << '\n';
	if(t % 10 == 0) {
	    of << t;
	    for(unsigned x=0; x<200; ++x) {
		double s[3] = {0.0,0.0,0.0};
		for(unsigned i=0; i<3; ++i)
		    for(unsigned j=0; j<3; ++j) 
			for(unsigned k=0; k<3; ++k) {
			    s[k] += freq[x][k][i][j];
			}
		    
		of << "," << (0.5*s[1]+s[2])/(s[0]+s[1]+s[2]);
	    }
	    for(unsigned x=0; x<200; ++x) {
		double s[3] = {0.0,0.0,0.0};
		for(unsigned i=0; i<3; ++i)
		    for(unsigned j=0; j<3; ++j) 
			for(unsigned k=0; k<3; ++k) {
			    s[k] += freq[x][i][k][j];
			}
		    
		of << "," << (0.5*s[1]+s[2])/(s[0]+s[1]+s[2]);
	    }
	    for(unsigned x=0; x<200; ++x) {
		double s[3] = {0.0,0.0,0.0};
		for(unsigned i=0; i<3; ++i)
		    for(unsigned j=0; j<3; ++j) 
			for(unsigned k=0; k<3; ++k) {
			    s[k] += freq[x][i][j][k];
			}
		    
		of << "," << (0.5*s[1]+s[2])/(s[0]+s[1]+s[2]);
	    }
	    of << endl;
	}
	if(t<6000) simulate_mating(true);
    }
    of.close();
    return 0;
}
