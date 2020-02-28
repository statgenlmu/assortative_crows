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
