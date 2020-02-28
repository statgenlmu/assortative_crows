#include<map>
#include<vector>
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<cmath>
#include "genotype_data.hh"

double genoprob(unsigned m, unsigned f, unsigned k) {
    const double g[27]={1 , 0.5 , 0 , 0.5 , 0.25, 0  , 0 , 0 , 0 ,
			0 , 0.5 , 1 , 0.5 , 0.5 , 0.5, 1, 0.5, 0 ,
			0 , 0   , 0 , 0   , 0.25, 0.5, 0, 0.5, 1};
    if(k==999) return 1.0;
    return g[m+3*f+9*k];
}

double chick::prob(unsigned gm1, unsigned gm2, unsigned gf1, unsigned gf2) {
    // returns probability of chick genotype if parents have
    // genotype (loc1=gm1,loc2=gm2) and (loc1=gf1,loc2=gf2)
    return genoprob(gm1,gf1,loc1)*genoprob(gm2,gf2,loc2);
}

double nest::prob(unsigned gm1, unsigned gm2, unsigned gf1, unsigned gf2) {
    // returns probability of all chicks genotypes if parents have
    // genotype (loc1=gm1,loc2=gm2) and (loc1=gf1,loc2=gf2)

    double prod=1.0;
    for(auto it=ch.begin(); it!=ch.end(); ++it) {
	prod *= it->prob(gm1,gm2,gf1,gf2);
    }
    return prod;
}

double nest::logprob(const double pp[3][3][3][3]) {
    // for frequencies of pairs of parental genotypes return
    // log total probability of chick genotypes
    double sum=0.0;
    for(unsigned gm1=0; gm1<3; ++gm1)
	for(unsigned gm2=0; gm2<3; ++gm2)
	    for(unsigned gf1=0; gf1<3; ++gf1)
		for(unsigned gf2=0; gf2<3; ++gf2) {
		    sum += pp[gm1][gm2][gf1][gf2] * prob(gm1,gm2,gf1,gf2);
		}
    return log(sum);
}

double nest::logprob(const double pp[200][3][3][3][3], double hzlon) {
    // for frequencies of pairs of parental genotypes return
    // log total probability of chick genotypes, assuming that hzlon
    // is the longitude location of the hybrid zone origin.
    int li=int(std::round((dtrans-hzlon-2.5)/5.0)+100.0);
    if(li < 0) {
	li = 0;
    } else {
	if (li > 199) {
	    li = 199;
	}
    }
    // std::cout << longi << '\t' << hzlon << '\t' << unsigned(li) << std::endl;
    return logprob(pp[unsigned(li)]);
}

dataset::dataset(std::string filename) : nestm() {
    std::ifstream f(filename);
    std::string s;
    std::getline(f,s);
    unsigned nestcount=0, chickcount=0;
    while(!f.eof()) {
	std::string nestid;
	double lat,lon,dt;
	unsigned lo1,lo2;
	std::getline(f,s);
	std::stringstream ss(s);
	ss >> lat >> lon >> dt >> nestid;
	ss >> lo1;
	if(ss.fail()) {
	    lo1=999;
	}
	ss >> lo2;
	if(ss.fail()) {
	    lo2=999;
	}
	//std::cout << "***********************\n";
	//std::cout << dt << " " << nestid << " " << lo1 << " " << lo2 << std::endl;
	auto it=nestm.find(nestid);
	if(it==nestm.end()) {
	    nestm.insert(std::pair<std::string,nest> (nestid,nest(dt,chick(lo1,lo2))));
	    ++chickcount;
	    ++nestcount;
	    
	} else {
	    it -> second.add(chick(lo1,lo2));
	    ++chickcount;
	}
    }
    std::cout << "read data of " << chickcount << " chicks in " << nestcount << "nests.\n";
    f.close();
}

double dataset::logprob(const double pp[200][3][3][3][3], double hzlon) {
    // for frequencies of pairs of parental genotypes return
    // log total probability all genotypes, assuming that hzlon
    // is the longitude location of the hybrid zone origin.
    double sum = 0.0;
    for(auto it=nestm.begin(); it!=nestm.end(); ++it) {
	// it->second.print_summary();
	// std::cout << '\t' << it->second.logprob(pp, hzlon) << std::endl;
	sum += it->second.logprob(pp, hzlon);
    }
    return sum;
}
