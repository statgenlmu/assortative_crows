/*
 * hybridzone_simulator is a program in the project assortative_crows
 *
 * Copyright (C) 2021 Dirk Metzler, http://evol.bio.lmu.de/_statgen/
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

#include <iostream>
#include "mt_random.hh"
#include "crows.hh"

preference_matrix::preference_matrix(double zz, double mmp) {
    z=zz;
    min_mate_prob=mmp;
    for(unsigned i=0; i<3; ++i) {
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		for(unsigned l=0; l<3; ++l) {
		    w[i][j][k][l] = std::exp(-z*pow(pheno[i][j]-pheno[k][l],2));
		    if(w[i][j][k][l] < min_mate_prob)
			w[i][j][k][l] = min_mate_prob; 
		}
	    }
	}
    }
}

bool preference_matrix::decide(const crow & a, const crow & b) const {
    // random decision according to mating probability
    std::bernoulli_distribution
	bern(w[a.get_loc1()][a.get_loc2()][b.get_loc1()][b.get_loc2()]);
    // std::cout << "mating probability: " <<
    //    w[a.get_loc1()][a.get_loc2()][b.get_loc1()][b.get_loc2()] << std::endl;
    return bern(random_generator);
}
    
bool preference_matrix::decide(double mppa, const crow & b) const {
    // random decision according to mating probability
    std::bernoulli_distribution
	bern(wi(mppa, b.get_loc1(), b.get_loc2()));
    return bern(random_generator);
}
    
