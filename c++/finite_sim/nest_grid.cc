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
#include "mt_random.hh"
#include "crows.hh"
#include "nest.hh"
#include "nest_grid.hh"

nest_grid::nest_grid(double daf, double z, double min_mate_prob) :
    A(height,std::vector<nest>(width)), B(height,std::vector<nest>(width)),
    preference(z, min_mate_prob) {
    // daf is the locus-2 dark allele frequency in the west.
    crow m(0,0), f(0,0);
    std::bernoulli_distribution bern(daf);
    for(unsigned i=0; i<height; ++i) {
	for(unsigned j=width/2; j<width; ++j) {
	    A[i][j].set_couple(f,m);
	}
    }
    f.set1(2);
    m.set1(2);
    for(unsigned i=0; i<height; ++i) {
	for(unsigned j=0; j<width/2; ++j) {
	    f.set2(bern(random_generator)+bern(random_generator));
	    m.set2(bern(random_generator)+bern(random_generator));
	    A[i][j].set_couple(f,m);
	}
    }
    uptodate='A';
};

std::pair<unsigned,unsigned> nest_grid::dispers(unsigned i, unsigned j) const {
// propose coordinates of the nest from which a visitor comes to nest at i, j.
    std::bernoulli_distribution bern(0.767);
    int x, y;
    if(bern(random_generator)) {
	std::normal_distribution<double> gauss(0, 6.07);
	y= i + (int) round(gauss(random_generator));
	x= j + (int) round(gauss(random_generator));
    } else {
	std::normal_distribution<double> gauss(0, 28.12);
	y= i + (int) round(gauss(random_generator));
	x= j + (int) round(gauss(random_generator));
    }
    while( x<0 || x>=(int)width ) {
	if(x<0) x = -x;
	if(x>=(int)width) x = 2*(width-1) - x;
    }
    while( y<0 || y>=(int)height ) {
	if(y<0) y = -y;
	if(y>=(int)height) y = 2*(height-1) - y;
    }
    return std::pair<unsigned,unsigned>(unsigned(y),unsigned(x));
};

void nest_grid::simulate_generation(std::string imprint) {
    std::vector<std::vector<nest>> *kids, *parents;
    if(uptodate=='A') {
	parents = &A;
	kids = &B;
    } else {
	parents = &B;
	kids = &A;
    }
    for(unsigned i=0; i<height; ++i) {
        // std::cout << "line " << i << std::endl;
	for(unsigned j=0; j<width; ++j) {
	    std::pair<unsigned,unsigned> fyx, myx;
	    crow fc, mc;
	    bool mate = false;
	    while(mate == false) {
		fyx = dispers(i,j);
		myx = dispers(i,j);
		fc = (*parents)[fyx.first][fyx.second].offspring();
		mc = (*parents)[myx.first][myx.second].offspring();
		if(imprint=="mean") {
		    // use mid-parental phenotype of crow mc for decision
		    mate = preference.decide((*parents)[myx.first][myx.second].get_pheno(),fc);
		} else if (imprint=="none") {
		    mate = preference.decide(fc,mc);
		} else if (imprint=="father") {
		    mate = preference.decide(fc,(*parents)[myx.first][myx.second].get_father());
		} else if (imprint=="mother") {
		    mate = preference.decide(fc,(*parents)[myx.first][myx.second].get_mother());
		} else {
		    std::cerr << "unknown imprinting model\n";
		    exit(1);
		}
	    }
	    (*kids)[i][j].set_couple(fc,mc);
	}
    }
    if(uptodate=='A') {
	uptodate='B';
    } else {
	uptodate='A';
    }
};

std::vector<std::vector<double>> nest_grid::get_phenos() const{
    std::vector<std::vector<double>> v;
    const std::vector<std::vector<nest>> *p;
    v.resize(height);
    if(uptodate=='A') p=&A; else p=&B;
    for (unsigned i=0; i<height; ++i) {
	v[i].resize(width);
	for(unsigned j=0; j<width; ++j)
	    v[i][j]=(*p)[i][j].get_pheno();
    }
    return v;
}

std::vector<std::vector<unsigned>> nest_grid::get_genos() const{
    std::vector<std::vector<unsigned>> v;
    const std::vector<std::vector<nest>> *p;
    v.resize(height);
    if(uptodate=='A') p=&A; else p=&B;
    for (unsigned i=0; i<height; ++i) {
	v[i].resize(width);
	for(unsigned j=0; j<width; ++j)
	    v[i][j]=(*p)[i][j].get_geno();
    }
    return v;
}
