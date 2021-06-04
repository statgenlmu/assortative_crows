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
