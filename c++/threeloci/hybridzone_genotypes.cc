/*
 * hybridzone_simulator is a program in the project assortative_crows
 *
 * Copyright (C) 2020 Dirk Metzler, http://evol.bio.lmu.de/_statgen/
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
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "hybridzone_genotypes.hh"

double freq[200][3][3][3];
// freq[x][i][j][k] is frequency of genotype (i,j,k) in bin x

double disp_kernel[21]; // disp_kernel[y] is dispersal
// probability for a difference of y-10

double disp[200][200]; // disp[x][y] is dispersal
// probability from x to y

double dxxp[200][41]; // dxxp[x][dx] is sum over all suitable z of d[x][z]*d[x+dx-20][z]

double w[3][3][3][3]; // w[a][b][c][d] is the mating preference
// value for birds with genotypes (a,b) and (c,d)


void initialize_disp_kernel(double sigma) {
    const double invsqrt2pis=1/sqrt(2*pi)/sigma;
    double s=0;
    for(unsigned y=0; y<21; ++y) {
	disp_kernel[y] = invsqrt2pis*exp(-pow(5*(int(y)-10),2)/(2*sigma*sigma));
	s+=disp_kernel[y];
    }
    for(unsigned y=0; y<21; ++y) {
	disp_kernel[y]/=s;
    }
}

void initialize_disp_kernel_mix2norm_siefke() {
    double dk[21]={0.010601454, 0.004593536, 0.006005553, 0.007607918, 0.009340766, 0.011184832,
		   0.014224777, 0.027894418, 0.083391491, 0.194232956, 0.261844598, 0.194232956,
		   0.083391491, 0.027894418, 0.014224777, 0.011184832, 0.009340766, 0.007607918,
		   0.006005553, 0.004593536, 0.010601454};
    for(unsigned i=0; i<21; ++i) disp_kernel[i]=dk[i];
}

void initialize_disp_and_dxxp() {
    for(unsigned x=0; x<200; ++x) {
	for(unsigned y=0; y<200; ++y) {
	    // here y is the other index
	    disp[x][y]=0;
	}
	disp[x][x]=disp_kernel[10];
	for(unsigned y = 1; y<11; ++y) {
	    // now y is the difference between the bin indices
	    if(y>x) { // reflect at left boundary
		disp[x][y-x-1]+=disp_kernel[10-y];
	    } else {
		disp[x][x-y]=disp_kernel[10-y];
	    }
	    if(x+y>199) { // reflect at right boundary
		disp[x][399-(x+y)]+=disp_kernel[10+y];
	    } else {
		disp[x][x+y]=disp_kernel[10+y];
	    }
	}
	if(x<10) {
	    // here y is the other index
	    double s=0.0;
	    for(unsigned y=0; y<=x+10; ++y) {
		s+=disp[x][y];
	    }
	    for(unsigned y=0; y<=x+10; ++y) {
		disp[x][y]/=s;
	    }
	}
	if(x>189) {
	    // here y is the other index
	    double s=0.0;
	    for(unsigned y=x-10; y<=199; ++y) {
		s+=disp[x][y];
	    }
	    for(unsigned y=x-10; y<=199; ++y) {
		disp[x][y]/=s;
	    }
	}
    }
    
    for(unsigned x=0; x<200; ++x) {
	for(unsigned dx=0; dx<41; ++dx) {
	    dxxp[x][dx]=0; // remember: xp=x+dx-20, xp being the other index
	    if(x+dx >= 20 && x+dx < 220) {   // this means xp >= 0 && xp xp <200
		if(dx<20) { 
		    // xp = x-20+dx < x  =>  z goes from x-10 to xp+10=x+dx-10
		    for(unsigned z=(x > 10 ? x-10 : 0); z+10 <= ((x+dx) < 209 ? x+dx : 209)  ; ++z) {
			dxxp[x][dx]+=disp[x][z]*disp[x+dx-20][z];
		    }
		} else {
		    // xp = x-20+dx >= x  =>  z goes from xp-10 = x-30+dx to x+10
		    for(unsigned z=((x+dx>30) ? x-30+dx : 0) ; z <= ((x+10) < 200 ? x+10 : 199); ++z) {
			dxxp[x][dx]+=disp[x][z]*disp[x+dx-20][z];
		    }
		}
	    }
	}
    }
}

void initialize_freq(double f) {
    // f is the frequency of the second locus light allele in the west.
    if (f<0) f=0.0;
    if (f>1) f=1.0;
    for(unsigned x=0; x<100; ++x) {
	for(unsigned i=0; i<3; ++i) {
	    for(unsigned j=0; j<3; ++j) {
		for(unsigned k=0; k<3; ++k) {
		    freq[x][i][j][k]=(1e-12)*x;
		}
	    }
	}
	freq[x][2][0][2]=(1-f)*(1-f);
	freq[x][2][1][2]=2*f*(1-f);
	freq[x][2][2][2]=f*f-(1e-12)*6*x;
    }
    for(unsigned x=100; x<200; ++x) {
	for(unsigned i=0; i<3; ++i) {
	    for(unsigned j=0; j<3; ++j) {
		for(unsigned k=0; k<3; ++k) {
		    freq[x][i][j][k]=1e-12*(200-x);
		}
	    }
	}
	freq[x][0][0][0]=1-1e-12*8*(200-x);
    }
}

void initialize_freq() {
    initialize_freq(0);
}

void read_freq(string filename) {
    ifstream fin(filename);
    for(unsigned x=0; x<200; ++x) {
	unsigned y;
	fin >> y;
	if (x!=y) cerr << "unexpected beginning of line " << y << " instead of " << x << endl; 
	for(unsigned i=0; i<3; ++i) {
	    for(unsigned j=0; j<3; ++j) {
		for(unsigned k=0; k<3; ++k) {
		    char comma;
		    fin >> comma >> freq[x][i][j][k]; 
		}
	    }
	}
    }
    fin.close();
}

void initialize_w(char s1[], char s2[]) {
    if(string(s1)=="sbd") {
	initialize_w_same_but_different(atof(s2));
    } else {
	initialize_w(atof(s1),atof(s2));
    }
}

void initialize_w(double z, double min_mate_prob) {
    double pheno[3][3] = {0.0, 0.3743380, 0.7339315, 0.2748254, 0.6409001, 0.9094524,
			  1.0000000, 1.0000000, 1.0000000};
    for(unsigned i=0; i<3; ++i) {
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		for(unsigned l=0; l<3; ++l) {
		    w[i][j][k][l] = exp(-z*pow(pheno[i][j]-pheno[k][l],2));
		    if(w[i][j][k][l] < min_mate_prob)
			w[i][j][k][l] = min_mate_prob; 
		}
	    }
	}
    }
}

void initialize_w_speed(double z, double min_mate_prob) {
    double pheno[3][3] = {0.0, 0.0, 0.0,
			  0.5, 0.5, 0.5,
			  1.0000000, 1.0000000, 1.0000000};
    for(unsigned i=0; i<3; ++i) {
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		for(unsigned l=0; l<3; ++l) {
		    w[i][j][k][l] = exp(-z*pow(pheno[i][j]-pheno[k][l],2));
		    if(w[i][j][k][l] < min_mate_prob)
			w[i][j][k][l] = min_mate_prob; 
		}
	    }
	}
    }
}

void initialize_w_neutral() {
    for(unsigned i=0; i<3; ++i) {
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		for(unsigned l=0; l<3; ++l) {
		    w[i][j][k][l] = 1;
		}
	    }
	}
    }
}

void initialize_w_same_but_different(double thr) {
    for(unsigned i=0; i<3; ++i) {
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		for(unsigned l=0; l<3; ++l) {
		    if(i==0 && j==0) {
			if (k==0 && l==0)
			    w[i][j][k][l] = 1.0;
			else
			    w[i][j][k][l] = thr;
		    } else {
			if(i==2) {
			    if(k==2)
				w[i][j][k][l] = 1.0;
			    else
				w[i][j][k][l] = thr;
			} else {
			    if((k==0 && l==0) || k==2)
				w[i][j][k][l] = thr;
			    else
				w[i][j][k][l] = 1.0;
			}
		    }
		}
	    }
	}
    }
    
}



void print_w() {
    for(unsigned i=0; i<3; ++i) {
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		for(unsigned l=0; l<3; ++l) {
		    cout << w[i][j][k][l] << " ";
		}
	    }
	    cout << "\n";
	}
    }
}

void print_allele_frequencies() {
    for(unsigned i=0; i<200; ++i) {
	double s1=0.0, s2=0.0;
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		s1+=freq[i][1][j][k];
		s1+=freq[i][2][j][k];
	    }
	}
	cout << 0.5*s1+s2 << (i < 199 ? "," : "\n");
    }
    for(unsigned i=0; i<200; ++i) {
	double s1=0.0, s2=0.0;
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		s1+=freq[i][j][1][k];
		s1+=freq[i][j][2][k];
	    }
	}
	cout << 0.5*s1+s2 << (i < 199 ? "," : "\n");
    }
    for(unsigned i=0; i<200; ++i) {
	double s1=0.0, s2=0.0;
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		s1+=freq[i][j][k][1];
		s1+=freq[i][j][k][2];
	    }
	}
	cout << 0.5*s1+s2 << (i < 199 ? "," : "\n");
    }
}

void simulate_mating(bool rescale=true) {
    
    double freq_new[200][3][3][3];

    for(unsigned z=0; z<200; ++z) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		for(unsigned g3=0; g3<3; ++g3) {
		    freq_new[z][g1][g2][g3] = 0.0;
		}
	    }
	}
	for(unsigned x= (z>10 ? z-10 : 0) ; x<= (z<190 ? z+10 : 199); ++x) { // origin of mother
	    for(unsigned xp= (z>10 ? z-10 : 0) ; xp <= (z<190 ? z+10 : 199); ++xp) { // origin of father
		for(unsigned g1=0; g1<3; ++g1) { // mother genes
		    for(unsigned g2=0; g2<3; ++g2) {
			for(unsigned g3=0; g3<3; ++g3) {
			    for(unsigned gp1=0; gp1<3; ++gp1) { // father genes
				for(unsigned gp2=0; gp2<3; ++gp2) {
				    for(unsigned gp3=0; gp3<3; ++gp3) {
					// frequency of this pairing
					double fp = w[g1][g2][gp1][gp2] * freq[x][g1][g2][g3] * disp[x][z] * 
					    freq[xp][gp1][gp2][gp3] * disp[xp][z];
					
					// prob to get allele 1 from mother/father at locus 1 or 2:
					double  m1=(g1==0 ? 0 : (g1==1 ? 0.5 : 1)),
					    m2=(g2==0 ? 0 : (g2==1 ? 0.5 : 1)),
					    m3=(g3==0 ? 0 : (g3==1 ? 0.5 : 1)),
					    f1=(gp1==0 ? 0 : (gp1==1 ? 0.5 : 1)),
					    f2=(gp2==0 ? 0 : (gp2==1 ? 0.5 : 1)),
					    f3=(gp3==0 ? 0 : (gp3==1 ? 0.5 : 1)); 
					freq_new[z][0][0][0] += fp * ((1-m1)*(1-f1)) * ((1-m2)*(1-f2)) * (1-m3)*(1-f3);
					freq_new[z][0][1][0] += fp * ((1-m1)*(1-f1)) * (m2*(1-f2)+(1-m2)*f2) * (1-m3)*(1-f3);
					freq_new[z][0][2][0] += fp * ((1-m1)*(1-f1)) * (m2*f2) * (1-m3)*(1-f3);
					freq_new[z][1][0][0] += fp * (m1*(1-f1)+(1-m1)*f1) * ((1-m2)*(1-f2)) * (1-m3)*(1-f3);
					freq_new[z][1][1][0] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*(1-f2)+(1-m2)*f2) * (1-m3)*(1-f3);
					freq_new[z][1][2][0] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*f2) * (1-m3)*(1-f3);
					freq_new[z][2][0][0] += fp * (m1*f1) * ((1-m2)*(1-f2)) * (1-m3)*(1-f3);
					freq_new[z][2][1][0] += fp * (m1*f1) * (m2*(1-f2)+(1-m2)*f2) * (1-m3)*(1-f3);
					freq_new[z][2][2][0] += fp * (m1*f1) * (m2*f2) * (1-m3)*(1-f3);
					
					freq_new[z][0][0][1] += fp * ((1-m1)*(1-f1)) * ((1-m2)*(1-f2)) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][0][1][1] += fp * ((1-m1)*(1-f1)) * (m2*(1-f2)+(1-m2)*f2) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][0][2][1] += fp * ((1-m1)*(1-f1)) * (m2*f2) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][1][0][1] += fp * (m1*(1-f1)+(1-m1)*f1) * ((1-m2)*(1-f2)) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][1][1][1] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*(1-f2)+(1-m2)*f2) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][1][2][1] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*f2) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][2][0][1] += fp * (m1*f1) * ((1-m2)*(1-f2)) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][2][1][1] += fp * (m1*f1) * (m2*(1-f2)+(1-m2)*f2) * (m3*(1-f3)+(1-m3)*f3);
					freq_new[z][2][2][1] += fp * (m1*f1) * (m2*f2) * (m3*(1-f3)+(1-m3)*f3);
					
					freq_new[z][0][0][2] += fp * ((1-m1)*(1-f1)) * ((1-m2)*(1-f2)) * m3*f3;
					freq_new[z][0][1][2] += fp * ((1-m1)*(1-f1)) * (m2*(1-f2)+(1-m2)*f2) * m3*f3;
					freq_new[z][0][2][2] += fp * ((1-m1)*(1-f1)) * (m2*f2) * m3*f3;
					freq_new[z][1][0][2] += fp * (m1*(1-f1)+(1-m1)*f1) * ((1-m2)*(1-f2)) * m3*f3;
					freq_new[z][1][1][2] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*(1-f2)+(1-m2)*f2) * m3*f3;
					freq_new[z][1][2][2] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*f2) * m3*f3;
					freq_new[z][2][0][2] += fp * (m1*f1) * ((1-m2)*(1-f2)) * m3*f3;
					freq_new[z][2][1][2] += fp * (m1*f1) * (m2*(1-f2)+(1-m2)*f2) * m3*f3;
					freq_new[z][2][2][2] += fp * (m1*f1) * (m2*f2) * m3*f3;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	double sum=0.0;
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		for(unsigned g3=0; g3<3; ++g3) {
		    sum+=freq_new[z][g1][g2][g3];
		}
	    }
	}
	cout << "\t" << sum;
	if(rescale) {
	    for(unsigned g1=0; g1<3; ++g1) {
		for(unsigned g2=0; g2<3; ++g2) {
		    for(unsigned g3=0; g3<3; ++g3) {
			freq_new[z][g1][g2][g3]/=sum;
		    }
		}
	    }
	}
    }
    for(unsigned z=0; z<200; ++z) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		for(unsigned g3=0; g3<3; ++g3) {
		    freq[z][g1][g2][g3]=freq_new[z][g1][g2][g3];
		}
	    }
        }
    }
    cout << endl;
}
    
