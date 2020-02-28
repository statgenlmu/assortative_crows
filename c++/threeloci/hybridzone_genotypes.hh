#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

static const double pi=3.141593;

extern double freq[200][3][3][3];
// freq[x][i][j] is frequency of genotype (i,j,k) in bin x

extern double disp_kernel[21]; // disp_kernel[y] is dispersal
// probability for a difference of y-10

extern double disp[200][200]; // disp[x][y] is dispersal
// probability from x to y

extern double dxxp[200][41]; // dxxp[x][dx] is sum over all suitable z of d[x][z]*d[x+dx-20][z]

extern double w[3][3][3][3]; // w[a][b][c][d] is the mating preference
// value for birds with genotypes (a,b) and (c,d)

void initialize_disp_kernel(double sigma);
void initialize_disp_kernel_mix2norm_siefke();
void initialize_disp_and_dxxp();
void initialize_freq_speed();
void initialize_freq_speed_eq();
void initialize_freq();
void initialize_freq(double f); // f is the initial frequency of locus 2 light allele in the west
void read_freq(string filename);
void initialize_w(char s1[], char s2[]);
void initialize_w(double z, double min_mate_prob);
void initialize_w_speed(double z, double min_mate_prob);
void initialize_w_neutral();
void initialize_w_same_but_different(double thr);
void print_w();
void print_allele_frequencies();
void simulate_mating(bool rescale);
