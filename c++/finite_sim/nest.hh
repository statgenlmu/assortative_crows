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

class nest {
private:
    crow father, mother;
public:
    inline void set_couple(const crow & f, const crow & m) {
	father=f, mother=m;
    }
    inline void set_father(const crow & f) {
	father=f;
    }
    inline void set_mother(const crow & m) {
	mother=m;
    }
    inline crow get_father() {
	return father;
    }
    inline crow get_mother() {
	return mother;
    }
    crow offspring() const;
    // simulates and returns the genotype of a chick of this nest
    
    inline double get_pheno() const {return (father.get_pheno()+mother.get_pheno())/2.0;}
    inline unsigned get_geno() const {return father.get_geno()*100+mother.get_geno();}
};
