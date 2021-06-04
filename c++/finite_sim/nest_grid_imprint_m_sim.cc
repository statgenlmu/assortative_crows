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
#include <iomanip>  
#include "crows.hh"
#include "nest.hh"
#include "mt_random.hh"
#include "nest_grid.hh"

std::random_device rdev;
std::mt19937 random_generator(rdev());

int main() {
    nest_grid ng(0.7776,2.171,0.6173);
    for(int g=0; g<5000; ++g) {
	// std::cout << "generation " << g << std::endl; 
	ng.simulate_generation("mother"); // here imprint is set to true !
	if(g==0 || g==99 || g==999 || g==1999 || g==4999) {
	    std::cout << "##\n## Phenotypes after " << g+1 << " generations:\n";
	    std::vector<std::vector<double>> v=ng.get_phenos();
	    for(unsigned i=0; i<v.size(); ++i) {
		for(unsigned j=0; j<v[i].size()-1; ++j) {
		     std::cout << std::setprecision(2) << v[i][j] << ',';
		 }
		std::cout << std::setprecision(2) << v[i][v[i].size()-1] << std::endl;
	    }
	}
	if(g==0 || g==1999 || g==4999 ) {
	    std::cout << "##\n## Genotypes after " << g+1 << " generations:\n";
	    std::vector<std::vector<unsigned>> v=ng.get_genos();
	    for(unsigned i=0; i<v.size(); ++i) {
		for(unsigned j=0; j<v[i].size()-1; ++j) {
		    std::cout << v[i][j] << ',';
		}
		std::cout << v[i][v[i].size()-1] << std::endl;
	    }
	}
    }
    return 0;
}
