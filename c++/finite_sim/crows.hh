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

#include <cmath>

const double pheno[3][3] = {0.0, 0.3743380, 0.7339315,
			    0.2748254, 0.6409001, 0.9094524,
			    1.0000000, 1.0000000, 1.0000000};
class crow {
private:
    unsigned loc1, loc2;
public:
    crow(unsigned l1, unsigned l2) {
	loc1=l1; loc2=l2;
    }
    crow(){loc1=0; loc2=0;};
    inline unsigned get_loc1() const {return loc1;}
    inline unsigned get_loc2() const {return loc2;}
    void set1(unsigned a) {loc1=a;}
    void set2(unsigned b) {loc2=b;}
    inline double get_pheno() const {return pheno[loc1][loc2];}
    inline unsigned get_geno() const {return loc1*10+loc2;}
};

class preference_matrix {
private:
    double w[3][3][3][3], z, min_mate_prob;
public:
    preference_matrix(double zz, double mmp);
    bool decide(const crow & a, const crow & b) const;
    // random decision according to mating probability
    bool decide(double mppa, const crow & b) const;
    // mppa is mid-parental phenotype of crow a, which is used in imprinting model
    inline double wi(double mppa, unsigned b_loc1, unsigned b_loc2) const {
	// used instead of preference matrix in the case of the imprinting model
	double r=std::exp(-z*pow(mppa-pheno[b_loc1][b_loc2],2));
	if(r>min_mate_prob) return r;
	return min_mate_prob;
    }
};
