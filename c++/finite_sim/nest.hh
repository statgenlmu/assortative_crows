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
