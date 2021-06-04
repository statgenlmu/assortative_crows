#include <iostream>
#include "../crows.hh"
#include "../nest.hh"
#include "../mt_random.hh"
#include "../nest_grid.hh"

std::random_device rdev;
std::mt19937 random_generator(rdev());

int main() {
    nest_grid ng(0.7776,2.171,0.6173);
    for(int i=0; i<2000; ++i) {
	std::cout << "generation " << i << std::endl; 
	ng.simulate_generation();
    }
    return 0;   
}
