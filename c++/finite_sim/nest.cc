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

#include "mt_random.hh"
#include "crows.hh"
#include "nest.hh"

crow nest::offspring() const {
    // simulates and returns the genotype of a chick of this nest
    std::bernoulli_distribution bern(0.5);
    
    unsigned l1=0, l2=0;
    if(father.get_loc1()>0) {
	if(father.get_loc1()>1) {
	    ++l1;
	} else {
	    l1 += bern(random_generator);
	}
    }
    if(mother.get_loc1()>0) {
	if(mother.get_loc1()>1) {
	    ++l1;
	} else {
	      l1 += bern(random_generator);
	}
    }
    if(father.get_loc2()>0) {
	if(father.get_loc2()>1) {
	    ++l2;
	} else {
	        l2 += bern(random_generator);
	}
    }
    if(mother.get_loc2()>0) {
	if(mother.get_loc2()>1) {
	    ++l2;
	} else {
	   	    l2 += bern(random_generator);
	}
    }
    return crow(l1,l2);
}
