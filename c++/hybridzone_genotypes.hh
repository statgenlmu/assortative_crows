#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

static const double pi=3.141593;

extern double freq[200][3][3];
// freq[x][i][j] is frequency of genotype (i,j) in bin x

extern double disp_kernel[21]; // disp_kernel[y] is dispersal
// probability for a difference of y-10

extern double disp[200][200]; // disp[x][y] is dispersal
// probability from x to y

extern double dxxp[200][41]; // dxxp[x][dx] is sum over all suitable z of d[x][z]*d[x+dx-20][z]

extern double v[200][3][3];  // crow partner search intensities to
// make sure that each bird has the same fitness

extern double w[3][3][3][3]; // w[a][b][c][d] is the mating preference
// value for birds with genotypes (a,b) and (c,d)

void initialize_disp_kernel(double sigma);
void initialize_disp_kernel_mix2norm_siefke();
void initialize_disp_and_dxxp();
void initialize_freq_speed();
void initialize_freq_speed_eq();
void initialize_freq();
void initialize_freq(double f); // f is the initial frequency of locus 2 light allele in the west
void initialize_v(double a=3, double b=1);
void read_freq(string filename);
void initialize_w(char s1[], char s2[]);
void initialize_w(double z, double min_mate_prob);
void initialize_w_additive(char s1[], char s2[]);
void initialize_w_additive(double z, double min_mate_prob);
void initialize_w_speed(double z, double min_mate_prob);
void initialize_w_neutral();
void initialize_w_same_but_different(double thr);
void initialize_w_pure(double thr1,double thr2);
void print_w();
void print_allele_frequencies();
int weight_v_deviation_f(const gsl_vector * v, void *params,
			 gsl_vector * f);
int weight_v_deviation_df(const gsl_vector * v, void *params,
			  gsl_matrix * J);
int weight_v_deviation_fdf(const gsl_vector * v, void *params,
			   gsl_vector * f, gsl_matrix * J);
void optimize_v_with_J(double v[200][3][3]);
void optimize_v_without_J(double v[200][3][3]);
void simulate_mating(bool rescale);
