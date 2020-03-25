
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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "hybridzone_genotypes.hh"

int main(int argc, char * argv[]) {
    if(argc<5) {
	cerr << "missing line arguments\n";
	exit(1);
    }
    //double assort_strength=atof(argv[1]), assort_threshold=0.0;
    string outputname(argv[4]);
    if(argc>6 && string(argv[6])=="additive")
	initialize_w_additive(argv[1],argv[2]);
    else
      	initialize_w(argv[1],argv[2]);
	// initialize_w_pure(atof(argv[1]),atof(argv[2]));
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

    if(argc>5 && string(argv[5])=="selection") {
	initialize_v(1,1);
	simulate_mating(true);
    } else {
	initialize_v();
	gsl_vector *vv = gsl_vector_alloc (1800), *ff = gsl_vector_alloc (1800);
	for(unsigned x=0; x<200; ++x) for(unsigned g1=0; g1<3; ++g1) for(unsigned g2=0; g2<3; ++g2)
									 gsl_vector_set(vv,9*x+3*g1+g2,v[x][g1][g2]);     
	weight_v_deviation_f(vv,NULL,ff);
	double sum=0;
	for(unsigned i=0; i<1800; ++i) sum+=abs(gsl_vector_get(ff,i));
	cout << "Summe :" << sum << endl;
	
	optimize_v_with_J(v);
	simulate_mating(false);
    }

    ofstream of(outputname+"simres.csv");
    for(unsigned t=0; t<2000; ++t) {
	cout << "t = " << t << '\n';
	if(argc>5 && string(argv[5])=="selection") {
	    simulate_mating(true);
	} else {
	    optimize_v_with_J(v);
	    simulate_mating(false);
	}
	if((t+1) % 10 == 0) {
	    of << t+1;
	    for(unsigned x=0; x<200; ++x) {
		double s = freq[x][0][0]+freq[x][0][1]+freq[x][0][2] +
		    freq[x][1][0]+freq[x][1][1]+freq[x][1][2] +
		    freq[x][2][0]+freq[x][2][1]+freq[x][2][2];
		of << "," << 0.5*(freq[x][1][0]+freq[x][1][1]+freq[x][1][2])/s+(freq[x][2][0]+freq[x][2][1]+freq[x][2][2])/s;
	    }
	    for(unsigned x=0; x<200; ++x) {
		double s = freq[x][0][0]+freq[x][0][1]+freq[x][0][2] +
		    freq[x][1][0]+freq[x][1][1]+freq[x][1][2] +
		    freq[x][2][0]+freq[x][2][1]+freq[x][2][2];
		of << "," << 0.5*(freq[x][0][1]+freq[x][1][1]+freq[x][2][1])/s+(freq[x][0][2]+freq[x][1][2]+freq[x][2][2])/s;
	    }
	    of << endl;
	}
    }
    of.close();
    ofstream offreq(outputname+"finfreqs.csv");
    for(unsigned x=0; x<200; ++x) {
	offreq << x;
	for(unsigned i=0; i<3; ++i) {
	    for(unsigned j=0; j<3; ++j) {
		offreq << ',' << freq[x][i][j]; 
	    }
	}
       offreq << endl;
    }
    offreq.close();
    return 0;
}
