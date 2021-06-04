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

#include <utility>
#include <vector>
#include <string>
#include <cmath>

class nest_grid {
private:
    std::vector<std::vector<nest>> A, B;
    char uptodate;
    preference_matrix preference;
public:
    const static unsigned width=1000, height=500;
    nest_grid(double daf, double z, double min_mate_prob);
    // daf is the locus-2 dark allele frequency in the west.
    std::pair<unsigned,unsigned> dispers(unsigned i, unsigned j) const;
    // propose coordinates of the nest from which a visitor comes to nest at i, j.
    void simulate_generation(std::string imprint="none");
    std::vector<std::vector<double>> get_phenos() const;
    std::vector<std::vector<unsigned>> get_genos() const;
};
