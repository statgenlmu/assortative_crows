/*
 * genotype_llh is a program in the project assortative_crows
 *
 * Copyright (C) 2020 Dirk Metzler, http://evol.bio.lmu.de/_statgen/
 *
 * This file is part of genotype_llh.
 *
 * genotype_llh is free software: you can redistribute it and/or modify
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


class chick {
private :
    unsigned loc1, loc2;
public :
    chick(unsigned l1, unsigned l2) {loc1=l1; loc2=l2;}
    double prob(unsigned gm1, unsigned gm2, unsigned gf1, unsigned gf2);
    // returns probability of chick genotype if parents have
    // genotype (loc1=gm1,loc2=gm2) and (loc1=gf1,loc2=gf2)
    
};

class nest {
private:
    double dtrans;
    std::vector<chick> ch;
public:
    nest(double dt, chick chi) : ch(1,chi) {
	dtrans=dt;
    }
    void add(chick chi) {ch.push_back(chi);}
    double prob(unsigned gm1, unsigned gm2, unsigned gf1, unsigned gf2);
    // returns probability of all chicks genotypes if parents have
    // genotype (loc1=gm1,loc2=gm2) and (loc1=gf1,loc2=gf2)
    double logprob(const double pp[3][3][3][3]);
    // for frequencies of pairs of parental genotypes return
    // log total probability of chick genotypes
    double logprob(const double pp[200][3][3][3][3], double hzlon);
    // for frequencies of pairs of parental genotypes return
    // log total probability of chick genotypes, assuming that hzlon
    // is the location of the hybrid zone origin.
    void print_summary() const {
	std::cout << dtrans << '\t' << ch.size() ;
    }
};

class dataset {
private:
    std::map<std::string, nest> nestm;
public:
    dataset(std::string filename);
    double logprob(const double pp[200][3][3][3][3],double hzlon);
    // for frequencies of pairs of parental genotypes return
    // log total probability all genotypes, assuming that hzlon
    // is the location of the hybrid zone origin.
};
