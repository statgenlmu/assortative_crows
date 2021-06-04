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
