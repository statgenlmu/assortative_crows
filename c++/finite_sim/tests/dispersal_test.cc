#include <iostream>
#include <utility>
#include "../crows.hh"
#include "../nest.hh"
#include "../mt_random.hh"
#include "../nest_grid.hh"

std::random_device rdev;
std::mt19937 random_generator(rdev());

int main() {
    nest_grid ng(0.6,2.8,0.6);
    std::pair<unsigned, unsigned> p;
    for(int i=0; i<10000; ++i) {
	p = ng.dispers(400,450);
	std::cout << p.first << "," << p.second << std::endl;
    }
    return 0;   
}
