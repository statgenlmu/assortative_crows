all: hybridzone_simulator

CXXFLAGS = -O2 -std=c++11 -Wall

hybridzone_simulator: hybridzone_simulator.cc  hybridzone_genotypes.o hybridzone_genotypes.hh
	g++ -o hybridzone_simulator hybridzone_simulator.cc hybridzone_genotypes.o $(CXXFLAGS) -lgsl -lgslcblas

hybridzone_genotypes.o: hybridzone_genotypes.cc hybridzone_genotypes.hh
	g++ -c hybridzone_genotypes.cc $(CXXFLAGS) -lgsl -lgslcblas
