all: hybridzone_simulator genotype_llh

CXXFLAGS = -O2 -std=c++11 -Wall

genotype_llh: genotype_llh.cc hybridzone_genotypes.o hybridzone_genotypes.hh genotype_data.o
	g++ -o genotype_llh genotype_llh.cc hybridzone_genotypes.o genotype_data.o $(CXXFLAGS) -lgsl -lgslcblas

hybridzone_simulator: hybridzone_simulator.cc  hybridzone_genotypes.o hybridzone_genotypes.hh
	g++ -o hybridzone_simulator hybridzone_simulator.cc hybridzone_genotypes.o $(CXXFLAGS) -lgsl -lgslcblas

hybridzone_genotypes.o: hybridzone_genotypes.cc hybridzone_genotypes.hh
	g++ -c hybridzone_genotypes.cc $(CXXFLAGS) -lgsl -lgslcblas

genotype_data.o: genotype_data.cc genotype_data.hh
	g++ -c genotype_data.cc $(CXXFLAGS)
