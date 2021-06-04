all: nest_grid_sim    nest_grid_imprint_sim   nest_grid_imprint_f_sim    nest_grid_imprint_m_sim

CXXFLAGS = -O2 -std=c++11 -Wall

nest_grid_sim: nest_grid_sim.cc nest.o crows.o nest_grid.o

nest_grid_imprint_sim: nest_grid_imprint_sim.cc nest.o crows.o nest_grid.o

nest_grid_imprint_f_sim: nest_grid_imprint_f_sim.cc nest.o crows.o nest_grid.o

nest_grid_imprint_m_sim: nest_grid_imprint_m_sim.cc nest.o crows.o nest_grid.o

nest.o : nest.cc mt_random.hh crows.hh nest.hh

crows.o : crows.cc mt_random.hh crows.hh 

nest_grid.o : nest_grid.cc mt_random.hh crows.hh nest.hh nest_grid.hh
