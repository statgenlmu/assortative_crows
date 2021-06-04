all: nest_grid_sim

CXXFLAGS = -g -std=c++11 -Wall

nest_grid_sim: ../nest_grid_sim.cc nest.o crows.o nest_grid.o
	g++ -o nest_grid_sim ../nest_grid_sim.cc nest.o crows.o nest_grid.o $(CXXFLAGS)

nest.o : ../nest.cc ../mt_random.hh ../crows.hh ../nest.hh
	g++ -c ../nest.cc $(CXXFLAGS)

crows.o : ../crows.cc ../mt_random.hh ../crows.hh
	g++ -c ../crows.cc $(CXXFLAGS)

nest_grid.o : ../nest_grid.cc ../mt_random.hh ../crows.hh ../nest.hh ../nest_grid.hh
	g++ -c ../nest_grid.cc $(CXXFLAGS)
