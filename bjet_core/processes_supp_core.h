// exact copy of the processes_supp from bjet
#ifndef __PROCESSES_SUPP_V02_H_
#define __PROCESSES_SUPP_V02_H_
/** @file */

extern double Simpson(double func[], double ics[], int res, int start, int end);
extern double linint(double x, double xvec[], int xdim, double x_min, double x_max);

extern double j_syn(double (*elec_spec)(double), double gamma_min, double gamma_max, double nu, double B, int prec1, int prec2);
extern double k_esa(double (*elec_spec)(double), double gamma_min, double gamma_max, double nu, double B, int prec1, int prec2);
extern double CylTransfEquat(double I_inp, double jj, double kk, double ll);
extern double SphTransfEquat(double jj, double kk, double ll);
extern double Intens2Flux(double Intens, double Radius, double Doppler, double z, double Hubble);
extern double CylIntens2Flux(double Intens1, double Intens2, double Radius, double Length, double Doppler, double z, double Hubble, double Theta);
extern double RingIntens2Flux(double Intens, double InnRadius, double OutRadius,  double Doppler, double z, double Hubble, int check);
extern double FreqTransS2O(double nu, double Doppler, double z);
extern double FreqTransO2S(double nu, double Doppler, double z);
extern double j_com(double (*elec_spec)(double), double gamma_min, double gamma_max, 
                    double I_rad[], double nu_rad_min, double nu_rad_max, int    nu_rad_dim,
	            double nu, int prec1, int prec2); 
extern double tau_IRA_Kneiske(double nu, double zz, int range);
extern double tau_IRA_Franceschini(double nu, double zz);
extern double tau_IRA_Finke(double nu, double zz);
extern double tau_IRA_Franceschini17(double nu, double zz);
extern double tau_IIR(double nu, double z, int level);		    
//extern double gg_abs_ssc(double (*elec_spec)(double), double gamma_min, double gamma_max, double nu, double B, double r, int sph_cyl, int prec1, int prec2);
extern double gg_abs(double nu_c, double I_syn[], int nu_dim, double nu_min, double nu_max, int prec1, int prec2); 
//extern double gg_abs_ssc2nd(double (*elec_spec)(double), double gamma_min, double gamma_max, double nu, double I_rad2nd[], double nu_rad_min, double nu_rad_max, int nu_rad_dim, double r, int sph_cyl, int prec1, int prec2);
//extern double gg_abs_eic(double (*elec_spec)(double), double gamma_min, double gamma_max, double I_rad[], double nu_rad_min, double nu_rad_max, int    nu_rad_dim, double nu_c, double r, int sph_cyl, int prec1, int prec2);
//extern double sig_gg(double,double, double NU[]);
extern double sigma_gg(double nu_s, double nu_c);
extern double intgl(double (*func)(double), double a, double b, int n, int m);
extern double Planck(double nu_BB, double T_BB);

//extern double LuminDist(const double redshift);
extern double Distance_Luminosity(double z, double H0, double WM);
extern double piondecay(double energy);// E^2 dN/dE  eV/(cm^2.s)
extern double piondecaytest(double energy);// E^2 dN/dE  eV/(cm^2.s)

extern void setHadronicParameters();
extern void setHadronicParametersTest();

extern double piondecay_ahaprecise_electron(double energy); // production of secondary electrons from pion decay

extern void energetics();

#endif
