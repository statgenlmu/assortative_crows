/*
 * genotype_llh is a program in the project assortative_crows
 *
 * Copyright (C) 2020 Dirk Metzler, http://evol.bio.lmu.de/_statgen/
 *
 * This file is part of genotype_llh.
 *
 * genotype_llh is free software: you can redistribute it and/or modify
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
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "genotype_data.hh"

using namespace std;

#include "hybridzone_genotypes.hh"

double pairprob[200][3][3][3][3]; 

void calculate_pairprobs() {
    for(unsigned z=0; z<200; ++z) {
	for(unsigned g1=0; g1<3; ++g1) { // mother genes
	    for(unsigned g2=0; g2<3; ++g2) {
		for(unsigned gp1=0; gp1<3; ++gp1) { // father genes
		    for(unsigned gp2=0; gp2<3; ++gp2) {
			pairprob[z][g1][g2][gp1][gp2]=0;
		    }
		}
	    }
	}
	for(unsigned x= (z>10 ? z-10 : 0) ; x<= (z<190 ? z+10 : 199); ++x) { // origin of mother
	    for(unsigned xp= (z>10 ? z-10 : 0) ; xp <= (z<190 ? z+10 : 199); ++xp) { // origin of father
		for(unsigned g1=0; g1<3; ++g1) { // mother genes
		    for(unsigned g2=0; g2<3; ++g2) {
			for(unsigned gp1=0; gp1<3; ++gp1) { // father genes
			    for(unsigned gp2=0; gp2<3; ++gp2) {
				pairprob[z][g1][g2][gp1][gp2] += w[g1][g2][gp1][gp2] * freq[x][g1][g2] *
				    disp[x][z] * v[x][g1][g2] * freq[xp][gp1][gp2] * disp[xp][z] * v[xp][gp1][gp2];
			    }
			}
		    }
		}
	    }
	}
    }
}

int main(int argc, char * argv[]) {
    if(argc < 5) {
	cerr << "please give the following command line arguments:\n"
	     << "input file name for genotype frequencies\n"
	     << "assortativity parameter value\n"
	     << "assortativity mininmum value\n"
	     << "data file name\n";
    }
    initialize_w(argv[2],argv[3]); 
    //  initialize_w_pure(atof(argv[2]),atof(argv[3])); 
    string ffn(argv[1]), dfn(argv[4]);
    read_freq(ffn);
    initialize_disp_kernel_mix2norm_siefke();
    initialize_disp_and_dxxp();
    if(argc>5 && string(argv[5])=="selection") {
	initialize_v(1,1);
    } else {
	initialize_v();
	optimize_v_with_J(v);
    }
    calculate_pairprobs();
    dataset dat(dfn);
    for(double h=300; h<=700; ++h) {
	// std::cout << "Calculate likelihood for the case that hybrid zone is at km " << h << std::endl;
	std::cout << h << '\t' << std::setprecision(9) 
		  << dat.logprob(pairprob,h) << std::endl;
    }
    return 0;
}
