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
