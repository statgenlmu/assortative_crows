#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "hybridzone_genotypes.hh"

double freq[200][3][3];
// freq[x][i][j] is frequency of genotype (i,j) in bin x

double disp_kernel[21]; // disp_kernel[y] is dispersal
// probability for a difference of y-10

double disp[200][200]; // disp[x][y] is dispersal
// probability from x to y

double dxxp[200][41]; // dxxp[x][dx] is sum over all suitable z of d[x][z]*d[x+dx-20][z]

double v[200][3][3];  // crow partner search intensities to
// make sure that each bird has the same fitness

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

void initialize_freq_speed() {
    for(unsigned x=0; x<200; ++x) {
	for(unsigned i=0; i<3; ++i) {
	    for(unsigned j=0; j<3; ++j) {
		freq[x][i][j]=1e-12*x;
	    }
	}
    }
    for(unsigned x=0; x<100; ++x) {
	freq[x][2][0]=1-1e-12*8*x;
    }
    for(unsigned x=100; x<150; ++x) {
	freq[x][0][0]=1-1e-12*8*x;
    }
    for(unsigned x=150; x<200; ++x) {
	freq[x][0][2]=1-1e-12*8*x;
    }
}

void initialize_freq_speed_eq() {              
    for(unsigned x=0; x<200; ++x) {
        for(unsigned i=0; i<3; ++i) {
            for(unsigned j=0; j<3; ++j) {
                freq[x][i][j]=1e-12*x;
            }
        }
    }
    for(unsigned x=0; x<100; ++x) {
        freq[x][2][0]=1-1e-12*8*x;
    }
    for(unsigned x=100; x<200; ++x) {                                                                  
        freq[x][0][2]=1-1e-12*8*x;                                                                     
    }                                                                                                 
}

void initialize_freq(double f) {
    // f is the frequency of the second locus light allele in the west.
    if (f<0) f=0.0;
    if (f>1) f=1.0;
    for(unsigned x=0; x<100; ++x) {
	for(unsigned i=0; i<3; ++i) {
	    for(unsigned j=0; j<3; ++j) {
		freq[x][i][j]=(1e-12)*x;
	    }
	}
	freq[x][2][0]=(1-f)*(1-f);
	freq[x][2][1]=2*f*(1-f);
	freq[x][2][2]=f*f-(1e-12)*6*x;
    }
    for(unsigned x=100; x<200; ++x) {
	for(unsigned i=0; i<3; ++i) {
	    for(unsigned j=0; j<3; ++j) {
		freq[x][i][j]=1e-12*(200-x);
	    }
	}
	freq[x][0][0]=1-1e-12*8*(200-x);
	// alternative: second locus heterozygous: freq[x][0][1]=1-1e-12*8*(200-x);
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
		char comma;
		fin >> comma >> freq[x][i][j]; 
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

void initialize_w_additive(char s1[], char s2[]) {
    initialize_w_additive(atof(s1),atof(s2));
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

void initialize_w_additive(double z, double min_mate_prob) {
    double pheno[3][3] = {0.0,0.2954393,0.5792417,0.2103791,0.5058184,0.7896209,0.4207583,
			  0.7161976,1.0};
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

void initialize_w_pure(double pL, double pD) {
    // hooded crows have a preference value pL to mate with other than hooded crows
    // black crows have preference value pD to make with other than dark crows.
    for(unsigned i=0; i<3; ++i) {
	for(unsigned j=0; j<3; ++j) {
	    for(unsigned k=0; k<3; ++k) {
		for(unsigned l=0; l<3; ++l) {
		    double x=1.0;
		    if(i==0 && j==0 && (k!=0 || l!=0)) x *= pL;
		    else if(k==0 && l==0 && (i!=0 || j!=0)) x *= pL;
		    if(i==2 && k!=2 || i!=2 && k==2) x *= pD;
		    w[i][j][k][l]=x;
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

void initialize_v(double a, double b) {
    for(unsigned x=0; x<100; ++x) {
	for(unsigned g1=0; g1<2; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		v[x][g1][g2]=a;
	    }
	}
	for(unsigned g2=0; g2<3; ++g2) {
	    v[x][2][g2]=b;
	}
    }

    for(unsigned x=100; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		v[x][g1][g2]=a;
	    }
	}
	v[x][0][0]=b;
    }
}

void print_allele_frequencies() {
    for(unsigned i=0; i<199; ++i) {
	cout << 0.5*(freq[i][1][0] + freq[i][1][1] + freq[i][1][2]) + freq[i][2][0] + freq[i][2][1] + freq[i][2][2] << ",";
    }
    cout << 0.5*(freq[199][1][0] + freq[199][1][1] + freq[199][1][2]) + freq[199][2][0] + freq[199][2][1] + freq[199][2][2] << "\n";
    for(unsigned i=0; i<199; ++i) {
	cout << 0.5*(freq[i][0][1] + freq[i][1][1] + freq[i][2][1]) + freq[i][0][2] + freq[i][1][2] + freq[i][2][2] << ",";
    }
    cout << 0.5*(freq[199][0][1] + freq[199][1][1] + freq[199][2][1]) + freq[199][0][2] + freq[199][1][2] + freq[199][2][2] << "\n";
}

int weight_v_deviation_f(const gsl_vector * v, void *params,
              gsl_vector * f) {
// optimal values of v are the root of the functions in f

    for(unsigned x=0; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		    double sum_xp=0.0;
		    unsigned uxp = (x > 20 ? x-20 : 0) , oxp = ((x+20) < 199 ? x+20 : 199);
		    for(unsigned xp=uxp; xp <= oxp; ++xp) {
			double  sum_gp=0.0;
			for(unsigned gp1=0; gp1<3; ++gp1) {
			    for(unsigned gp2=0; gp2<3; ++gp2) {
				double summand = freq[xp][gp1][gp2] * w[g1][g2][gp1][gp2] * gsl_vector_get(v,xp*9+gp1*3+gp2);
				sum_gp += summand;
			    }   
			} // end of loops over gp1 and gp2 but still within xp loop
			sum_xp += sum_gp * dxxp[x][xp+20-x];
		    }
		    gsl_vector_set(f,x*9+g1*3+g2,
				   sum_xp * gsl_vector_get(v,x*9+g1*3+g2) - 1);
	    } // end of g2 loop
	}
    }
    return GSL_SUCCESS;
}

int weight_v_deviation_df(const gsl_vector * v, void *params,
               gsl_matrix * J) {
// J will be derivative of weight_v_deviation_f(...)

    for(unsigned x=0; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		    double sum_xp_J=0.0;
		    unsigned uxp = (x > 20 ? x-20 : 0) , oxp = ((x+20) < 199 ? x+20 : 199);
		    for(unsigned xp=uxp; xp <= oxp; ++xp) {
			double sum_gp_J_wfv=0.0;
			for(unsigned gp1=0; gp1<3; ++gp1) {
			    for(unsigned gp2=0; gp2<3; ++gp2) {
				double xp_J_wf = freq[xp][gp1][gp2] * w[g1][g2][gp1][gp2],
				    summand = xp_J_wf * gsl_vector_get(v,xp*9+gp1*3+gp2);
				if(x==xp && g1==gp1 && g2==gp2) {
				    summand *= 2;
				} else {
				    // now comes derivative of f(x,g) wrt v_(x',g')
				    gsl_matrix_set(J,x*9+g1*3+g2,xp*9+gp1*3+gp2,
						   gsl_vector_get(v,x*9+g1*3+g2) * dxxp[x][xp+20-x] * xp_J_wf);
				}
				// next line will be used for derivative of f(x,g) wrt v_(x,g)
				sum_gp_J_wfv+=summand; 
			    }   
			} // end of loops over gp1 and gp2 but still within xp loop
			sum_xp_J += sum_gp_J_wfv * dxxp[x][xp+20-x];
		    }
		    gsl_matrix_set(J,x*9+g1*3+g2,x*9+g1*3+g2,sum_xp_J);
	    } // end of g2 loop
	}
    }
    return GSL_SUCCESS;
}

int weight_v_deviation_fdf(const gsl_vector * v, void *params,
                gsl_vector * f, gsl_matrix * J){
// both value and derivative of weight_v_deviation_f(...)
    
    for(unsigned x=0; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		    double sum_xp=0.0, sum_xp_J=0.0;
		    unsigned uxp = (x > 20 ? x-20 : 0) , oxp = ((x+20) < 199 ? x+20 : 199);
		    for(unsigned xp=uxp; xp <= oxp; ++xp) {
			double  sum_gp=0.0, sum_gp_J_wfv=0.0;
			for(unsigned gp1=0; gp1<3; ++gp1) {
			    for(unsigned gp2=0; gp2<3; ++gp2) {
				double xp_J_wf = freq[xp][gp1][gp2] * w[g1][g2][gp1][gp2],
				    summand = xp_J_wf * gsl_vector_get(v,xp*9+gp1*3+gp2);
				sum_gp += summand;
				if(x==xp && g1==gp1 && g2==gp2) {
				    summand *= 2;
				} else {
				    // now comes derivative of f(x,g) wrt v_(x',g')
				    gsl_matrix_set(J,x*9+g1*3+g2,xp*9+gp1*3+gp2,
						   gsl_vector_get(v,x*9+g1*3+g2) * dxxp[x][xp+20-x] * xp_J_wf);
				}
				// next line will be used for derivative of f(x,g) wrt v_(x,g)
				sum_gp_J_wfv+=summand; 
			    }   
			} // end of loops over gp1 and gp2 but still within xp loop
			sum_xp += sum_gp * dxxp[x][xp+20-x];
			sum_xp_J += sum_gp_J_wfv * dxxp[x][xp+20-x];
		    }
		    gsl_matrix_set(J,x*9+g1*3+g2,x*9+g1*3+g2,sum_xp_J);
		    gsl_vector_set(f,x*9+g1*3+g2,
				   sum_xp * gsl_vector_get(v,x*9+g1*3+g2) - 1);
	    } // end of g2 loop
	}
    }
    return GSL_SUCCESS;
}

void optimize_v_with_J(double v[200][3][3]) {
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;
    
    int status;
    size_t iter = 0;
    
    const size_t n = 1800;
    gsl_multiroot_function_fdf f = {&weight_v_deviation_f,&weight_v_deviation_df,&weight_v_deviation_fdf, 1800, NULL};
    
    gsl_vector *vv = gsl_vector_alloc (n);

    for(unsigned x=0; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		gsl_vector_set (vv, x*9+g1*3+g2, v[x][g1][g2]);
	    }
	}
    }
    
    T = gsl_multiroot_fdfsolver_gnewton;
    s = gsl_multiroot_fdfsolver_alloc (T, n);
    gsl_matrix_set_zero (s -> J);
    
    gsl_multiroot_fdfsolver_set (s, &f, vv);
   
    // for(unsigned x=0; x<200; ++x) {
    // 	cout << x;
    // 	for(unsigned g1=0; g1<3; ++g1)
    // 	    for(unsigned g2=0; g2<3; ++g2)
    // 			cout << '\t' << gsl_vector_get(s->f,x*9+g1*3+g2);	
    //     cout << endl;
    // }
    
    do {
	iter++;
	status = gsl_multiroot_fdfsolver_iterate (s);
	
	double sum=0;
	for(unsigned i=0; i<1800; ++i) sum+=abs(gsl_vector_get(s->f,i));
	cout << "Summe :" << sum << endl;
	// print_state (iter, s);
	
	if (status)
	    break;
	
	status = gsl_multiroot_test_residual (s->f, 1e-10);
    } while (status == GSL_CONTINUE && iter < 1000);
    
    cerr << "status = " << gsl_strerror (status) << "\n";
    

    for(unsigned x=0; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		v[x][g1][g2] = gsl_vector_get (s->x, x*9+g1*3+g2);
	    }
	}
    }
    
    gsl_multiroot_fdfsolver_free (s);
    gsl_vector_free (vv);
}

void optimize_v_without_J(double v[200][3][3]) {
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;
    
    int status;
    size_t iter = 0;
    
    const size_t n = 1800;
    gsl_multiroot_function f = {&weight_v_deviation_f, n, NULL};
    
    gsl_vector *vv = gsl_vector_alloc (n);

    for(unsigned x=0; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		gsl_vector_set (vv, x*9+g1*3+g2, v[x][g1][g2]);
	    }
	}
    }
    
    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, n);
    gsl_multiroot_fsolver_set (s, &f, vv);
    
    do {
	iter++;
	status = gsl_multiroot_fsolver_iterate (s);
	
	gsl_vector * dx=gsl_multiroot_fsolver_dx(s);
	for(unsigned i=0; i<n; ++i) {
	    if(gsl_vector_get(dx,i)!=0) {
		cout << i << " ***** " << gsl_vector_get(dx,i) << endl;
	    }
	}
	    
	if (status)
	    break;
	
	status = gsl_multiroot_test_residual (s->f, 1e-10);
    } while (status == GSL_CONTINUE && iter < 10000);
    
    cerr << "status = " << gsl_strerror (status) << "\n";
    
    for(unsigned x=0; x<200; ++x) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		v[x][g1][g2] = gsl_vector_get (s -> x, x*9+g1*3+g2);
	    }
	}
    }
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (vv);
}

void simulate_mating(bool rescale=true) {
    
    double freq_new[200][3][3];

    for(unsigned z=0; z<200; ++z) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		freq_new[z][g1][g2] = 0.0;
	    }
	}
	for(unsigned x= (z>10 ? z-10 : 0) ; x<= (z<190 ? z+10 : 199); ++x) { // origin of mother
	    for(unsigned xp= (z>10 ? z-10 : 0) ; xp <= (z<190 ? z+10 : 199); ++xp) { // origin of father
		for(unsigned g1=0; g1<3; ++g1) { // mother genes
		    for(unsigned g2=0; g2<3; ++g2) {
			for(unsigned gp1=0; gp1<3; ++gp1) { // father genes
			    for(unsigned gp2=0; gp2<3; ++gp2) {
				// frequency of this pairing
				double fp = w[g1][g2][gp1][gp2] * freq[x][g1][g2] * disp[x][z] * v[x][g1][g2] *
				    freq[xp][gp1][gp2] * disp[xp][z] * v[xp][gp1][gp2];
				// cout << w[g1][g2][gp1][gp2] << '\t' <<  freq[x][g1][g2]<< '\t' << disp[x][z] << '\t' << v[x][g1][g2] 
				//    << '\t' << freq[xp][gp1][gp2] << '\t' << disp[xp][z] << '\t' << v[xp][gp1][gp2]<< endl;
				
				// prob to get allele 1 from mother/father at locus 1 or 2:
				double  m1=(g1==0 ? 0 : (g1==1 ? 0.5 : 1)),
				    m2=(g2==0 ? 0 : (g2==1 ? 0.5 : 1)),
				    f1=(gp1==0 ? 0 : (gp1==1 ? 0.5 : 1)),
				    f2=(gp2==0 ? 0 : (gp2==1 ? 0.5 : 1)); 
				freq_new[z][0][0] += fp * ((1-m1)*(1-f1)) * ((1-m2)*(1-f2));
				freq_new[z][0][1] += fp * ((1-m1)*(1-f1)) * (m2*(1-f2)+(1-m2)*f2);
				freq_new[z][0][2] += fp * ((1-m1)*(1-f1)) * (m2*f2);
				freq_new[z][1][0] += fp * (m1*(1-f1)+(1-m1)*f1) * ((1-m2)*(1-f2));
				freq_new[z][1][1] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*(1-f2)+(1-m2)*f2);
				freq_new[z][1][2] += fp * (m1*(1-f1)+(1-m1)*f1) * (m2*f2);
				freq_new[z][2][0] += fp * (m1*f1) * ((1-m2)*(1-f2));
				freq_new[z][2][1] += fp * (m1*f1) * (m2*(1-f2)+(1-m2)*f2);
				freq_new[z][2][2] += fp * (m1*f1) * (m2*f2);
			    }
			}
		    }
		}
	    }
	}
	double sum=0.0;
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		//cout << z << '\t'  << g1 << '\t'  << g2 << '\t' <<  freq_new[z][g1][g2] << endl;
		sum+=freq_new[z][g1][g2];
	    }
	}
	cout << "\t" << sum;
	if(rescale) {
	    for(unsigned g1=0; g1<3; ++g1) {
		for(unsigned g2=0; g2<3; ++g2) {
		    freq_new[z][g1][g2]/=sum;
		}
	    }
	}
    }
    for(unsigned z=0; z<200; ++z) {
	for(unsigned g1=0; g1<3; ++g1) {
	    for(unsigned g2=0; g2<3; ++g2) {
		freq[z][g1][g2]=freq_new[z][g1][g2];
	    }
        }
    }
    cout << endl;
}
    
