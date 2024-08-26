// Slightly modified version of processes_supp_core.cpp for use in MCMC code; modified by Sarah Youngquist, Feb 2022 (see below)
// Only change is switching to using namespace bj_core02 and including bj_core.h instead of
// bj02 and bj02.h
/** @file */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>

//#include "catalog.h"
#include "bj_core.h"
#include "processes_supp_core.h"


using namespace std;
//using namespace catalog;
using namespace bj_core02;
//



// APPROXIMATION CONSTANTS

const double C1     = 0.78;
const double C2     = 0.25;
const double C3     = 2.175;

// OTHER CONSTANTS

const double MIN_FREQ = 1.0e+4;
const double MAX_FREQ = 1.0e+40;//





double fDistance; 		// in pc
double fBtesla;		// in Tesla
double Ksync;         //Facteur Ksync*E*E pour le calcul des pertes d'energie synchrotron
double Kcompt;       //Facteur Kcompt*E*E pour le calcul des pertes d'energie ic
double Kloss1;           //Facteur Kloss1*E*E pour le calcul des pertes totales en E*E  
double Kbremss;      //Facteur Kbremss*E pour le calcul des pertes d'energie bremsstrahlung

//Characteristic quantities
double Rcar;
double Tcar;
double Vcar;
double fVelocity;

// electron spectra
double fEmin;           // in eV
double fEmax; 		// in eV
double fEcut; 		// exp. cutoff in eV
double fSlope;		// slope of the electronic distribution
double fEnergy50;		// in 10^50 ergs
double fNorm;			// spectrum normalization
double fEnergy;		// average electron energy
double fNbelec;			// no. of electrons

// proton spectra
double fEnergy_prot;		// av. proton energy
double fNb_prot;			// no. of protons

double fFluxfact;		// per m^2 1/(4*PI*fDistance^2)

int NBSECONDARYELECTRONS; //fElectrons;//nb of secondary electrons for the synchrotron process
double fEnergy_electrons[10000];//energy of the secondary electrons produced by pp interaction
double fFlux_electrons[10000];//Flux of the secondary electrons produced by pp interaction


double fPerDecade;              //Numbers of steps per decade used for the integration in energy for radiative process

//

// PROTON SPECTRUM

double Kp;       // Normalization of the proton spectrum, in terms of Lorentz factor [cm^-3]
double ENERGY50; // total energy in protons [10^50 erg]
 double EPMIN;    // minimum proton energy [eV]
 double EPMAX;    // maximum proton energy [eV]
 double EPCUTOFF; // minimum proton energy [eV]
 double GAMMA_P_MAX; // maximum proton energy in terms of Lorentz factor [E/mp c^2]
 double SLOPEP;   // slope of proton energy spectrum
 double NP;       // slope of proton energy spectrum
 double ATOMICDENSITY; // atomic density of target protons [cm^-3]
 double SIZE; // radius of emitting zone [pc]
 double DISTANCE; // distance of the source [pc]
 double AGE; // Age of the source [s]
 double fNorm_prot;
 double RATIO; // ratio proton over electron number (ratio of total number of particles)
  
 int          CASE_PION;

/**
 * Calculates the value of N_p_PowLawExpCutoff function for a given energy.
 *
 * This function takes the energy in electron volts (eV) as input and calculates
 * the value of N_p_PowLawExpCutoff function. The N_p_PowLawExpCutoff function is
 * defined as fNorm_prot multiplied by energy raised to the power of -NP and
 * multiplied by the exponential of -energy divided by EPCUTOFF.
 *
 * @param energy The energy in eV for which to calculate the N_p_PowLawExpCutoff.
 * @return The calculated value of N_p_PowLawExpCutoff for the given energy.
 */
  double N_p_PowLawExpCutoff(double energy) { // energy in eV
    //return(Kp * pow(energy,-NP) * exp(-energy/EPCUTOFF));
    return(fNorm_prot * pow(energy,-NP) * exp(-energy/EPCUTOFF));
  }
  
  // Switch between different proton spectrum shapes
/**
 * Calculates the value of N_p function for a given energy.
 *
 * This function takes the energy in electron volts (eV) as input and calculates
 * the value of the N_p function. The N_p function is defined as the result of
 * the N_p_PowLawExpCutoff function.
 *
 * @param energy The energy in eV for which to calculate the N_p value.
 * @return The calculated value of N_p for the given energy.
 */
  double N_p(double energy){ // energy in eV
    double foo=0.;
    foo=N_p_PowLawExpCutoff(energy);
    return foo;
  }

/**
 * @brief Calculates the value of J_p_PowLawExpCutoff for a given energy.
 *
 * This function calculates the value of J_p_PowLawExpCutoff using the provided energy value. The energy is expected to be in electron volts (eV).
 *
 * @param energy The energy value for which to calculate the J_p_PowLawExpCutoff.
 * @return The calculated value of J_p_PowLawExpCutoff.
 *
 * @pre The variables Kp, NP, and EPCUTOFF must be defined and have valid values.
 */
  double J_p_PowLawExpCutoff(double energy) { // energy in eV
    //return(Kp * pow(energy,-NP) * exp(-energy/EPCUTOFF));
    return(Kp * pow(energy,-NP) * exp(-energy/EPCUTOFF));
  }
  
  // Switch between different proton spectrum shapes
/**
 * @brief Calculates the value of J_p for a given energy.
 *
 * This function calculates the value of J_p using the provided energy value. The energy is expected to be in electron volts (eV). The function internally calls J_p_PowLawExpCutoff to perform the calculation.
 *
 * @param energy The energy value for which to calculate the J_p.
 * @return The calculated value of J_p.
 *
 * @see J_p_PowLawExpCutoff
 *
 * @pre The variables Kp, NP, and EPCUTOFF must be defined and have valid values.
 */
  double J_p(double energy){ // energy in eV
    double foo=0.;
    foo=J_p_PowLawExpCutoff(energy);
    return foo;
  }


/*
********************************************************************************
*
* MODIFIED SIMPSON INTEGRATION ROUTINE
*
********************************************************************************
*/


/**
 * @brief Modified Simpson Integration Routine
 *
 * This function calculates the integral of a given function using the Modified
 * Simpson integration method.
 *
 * @param func Pointer to the array storing the function values
 * @param ics Pointer to the array storing the x-values of the function
 * @param res Length of the arrays func[] and ics[]
 * @param start Index of the starting point in the arrays func[] and ics[]
 * @param end Index of the ending point in the arrays func[] and ics[]
 *
 * @return The integral value calculated by the Modified Simpson integration method
 */
double Simpson(double func[], double ics[], int res, int start, int end){
  //ics[]: x, res: length x

  double partial = 0.0;
  double stepp[res+1];
  
   
  for (int i=1;i<=res-1;i++){
   stepp[i] = ics[i+1]-ics[i];
  // printf(" %e %e %e \n",ics[i],stepp[i],func[i]);
   }
  
     
  for (int i = start; i<end-1; i++){
       
    partial+= (func[i+2]/3.+ func[i+1]*4./3. + func[i]/3.) * (stepp[i+1]+stepp[i])/2.;
    i+=1;
    //printf("%e %d %d\n",partial, start, end);
  }
  
  // printf(" PAIR? %d\n ",(end-start+1)%2 );
  
  if( (end-start)%2 != 0 ) {
  
    partial+= (func[end] + func[end-1])*stepp[end-1]/2.;
   // printf(" %e %e %e \n ", func[end],func[end-1],stepp[end-1]);
  }
    
  
  //printf(" TOTAL : %e %d %d\n",partial, start, end);
 return partial; 

}

/*
********************************************************************************
*
* WEIGHTS FOR GAUSS-LEGENDRE INTEGRATING ROUTINE
*
********************************************************************************
*/
/**
 * \brief This function calculates the weights for a Gauss-Legendre integrating routine.
 *
 * The gauleg function calculates the weights for a Gauss-Legendre integrating routine
 * using the specified lower and upper bounds, number of points, and arrays to store
 * the calculated values. The weights are calculated using the Gauss-Legendre method.
 *
 * \param x1     The lower bound of the interval.
 * \param x2     The upper bound of the interval.
 * \param x      Pointer to the array to store the calculated points.
 * \param w      Pointer to the array to store the calculated weights.
 * \param n      The number of points to calculate.
 *
 * \return void
 *
 * \remark This function assumes that the arrays x and w have already been allocated
 *         with enough memory to store the calculated values.
 *
 * \remark The calculated points and weights are stored in the arrays x and w respectively,
 *         in the range from index 1 to n. The values at index 0 and index n+1 are not used.
 *
 * \remark The tolerance for convergence is set to EPS = 3.0e-11.
 *         If the maximum absolute difference between z and z1 is smaller than EPS,
 *         the iteration is considered to have converged and the loop is exited.
 */
void gauleg(double x1, double x2, double x[], double w[], int n) {
    int    m, j, i;
    double z1, z, xm, xl, pp, p3, p2, p1;
    double EPS = 3.0e-11;

    m  = (n + 1) / 2;
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);
    for (i = 1;i <= m;i++) {
       z = cos(3.141592654 * (i - 0.25) / (n + 0.5));
       do {
	 p1 = 1.0;
	 p2 = 0.0;
	 for (j = 1;j <= n;j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 -(j - 1.0) * p3) / j;
	 }
	 pp = n * (z * p1 - p2) / (z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1 / pp;
       } while (fabs(z - z1) > EPS);
       x[i]     = xm - xl * z;
       x[n+1-i] = xm + xl * z;
       w[i]     = 2.0 * xl / ((1.0 - z * z) * pp * pp);
       w[n+1-i] = w[i];
    }
}

/*
********************************************************************************
*
* COMBINED INTEGRATING ROUTINE, TRAPEZOID PLUS GAUSS-LEGENDRE
*
********************************************************************************
*/
/**
 * @brief This function calculates the integral of a given function using a combined integrating routine.
 * The routine uses both the trapezoid rule and Gauss-Legendre integration.
 *
 * @param func A pointer to the function to be integrated.
 * @param a The lower limit of integration.
 * @param b The upper limit of integration.
 * @param n The number of steps in the trapezoid rule.
 * @param m The number of points in the Gauss-Legendre quadrature.
 * @return The value of the integral.
 *
 * @details The function first checks if the lower limit 'a' is greater than or equal to the upper limit 'b'.
 * If so, it returns 0 and prints an error message. Otherwise, it calculates the step size based on the number of steps 'n' and
 * the logarithms of 'a' and 'b'.
 *
 * Then, it iterates from 1 to 'n-1' and calculates the lower and upper limits for each subinterval.
 * The gauleg function is called to calculate the Gauss-Legendre points and weights for each subinterval.
 * The function 'func' is evaluated at each point and multiplied by the corresponding weight.
 * The evaluations are summed up and multiplied by the step size to calculate the integral for each subinterval.
 * Finally, all the subinterval integrals are summed up to obtain the total integral, which is returned as the result.
 *
 * @note The function assumes that the provided function 'func' is continuous within the given limits 'a' and 'b'.
 * The function 'gauleg' is used internally to calculate the Gauss-Legendre points and weights for each subinterval.
 *
 * @see gauleg
 */
double intgl(double (*func)(double), double a, double b, int n, int m) {
  //a: start, b: end, n: nb steps
      int i, j;
      double x1, x2, stp, val, prt, sum;
      const int x_dim = m;
      double x[x_dim+1], w[x_dim+1];

      if (a >= b) { 
	std::cout << "ERROR: a >= b !!!\n" << std::endl;
        return 0.0;
      }	

      stp = (log10(b) - log10(a)) / (n - 1);
      val = log10(a);
      
      sum = 0.0;
      for (j = 1; j <= n-1; j++) {
         x1  = pow(10.0, val);
         x2  = pow(10.0, val + stp);
	 
	 gauleg(x1, x2, x, w, x_dim);
	 for (prt = 0.0, i = 1; i <= x_dim; i++) prt += w[i] * func(x[i]);
	 
	 sum = sum + prt;
	 
         val += stp;
      }

      return sum;
}
/*
********************************************************************************
*
* LINEARLY INTERPOLATION ON A LOG-LOG ARRAY
*
********************************************************************************
*/
/**
 * @brief Perform linear interpolation on a log-log array.
 *
 * @param x The value to interpolate at.
 * @param xvec The array of x values to interpolate between.
 * @param xdim The size of the xvec array.
 * @param x_min The minimum x value in the xvec array.
 * @param x_max The maximum x value in the xvec array.
 *
 * @return The interpolated value at x.
 *
 * This function performs linear interpolation on a log-log array. The x value provided is
 * used to determine the interpolated value between two adjacent values in the xvec array.
 * The x values in the xvec array must be sorted in ascending order. If the given x value is
 * less than or equal to zero, it is replaced with the x_min value. If any input parameters
 * are negative or if x_min is equal to x_max, an error message is displayed and the program exits.
 * The interpolated value is calculated based on the provided x and xvec values by determining the
 * position of x in the xvec array and using the surrounding values for interpolation. The interpolation
 * is performed in logarithmic space if the loglog flag is set, otherwise it is performed in linear space.
 * The resulting interpolated value is returned.
 */
double linint(double x,
	     double xvec[],
	     int xdim,
	     double x_min,
	     double x_max
	     ){

  double logx, logxmin, logxmax, dfound, out, step;

  int found; 

  if(x<=0.){
    x=x_min;
  }
  if(x_min<=0.){
    std::cout<<"error in linint, x_min="<<x_min<<std::endl;
    exit(1);
  }
  if(x_max<=0.){
    std::cout<<"error in linint, x_max="<<x_max<<std::endl;
    exit(1);
  } 
  if(x_min==x_max){
    std::cout<<"error in linint, x_max==x_min"<<std::endl;
    exit(1);
  }

    logx = log10(x);
    logxmin = log10(x_min);
    logxmax = log10(x_max);

     
  step = (logxmax-logxmin)/(xdim - 1);

  dfound = (logx-logxmin) / (logxmax-logxmin) * (xdim-1) + 1;
  // std::cout<<dfound<<std::endl;

  dfound = floor(dfound+0.5);
  found = int(dfound);
  if(found>xdim){
    found=xdim;
  }
  else if(found<1){
    found=1;
  }


  double lvec=-200., lvecp1=-200., lvecm1=-200.;
  const int loglog=1; // 1  interp. in log log   0  interp. in log lin

  if(loglog){
    if(xvec[found+1]>0.){lvecp1=log10(xvec[found+1]); }
    if(xvec[found]>0.){lvec=log10(xvec[found]); }
    if(xvec[found-1]>0.){lvecm1=log10(xvec[found-1]); }
    if(xvec[found-1]==0 || xvec[found]==0 || xvec[found+1]==0)return 0;
  }else{
    lvecp1=xvec[found+1];
    lvec=xvec[found];
    lvecm1=xvec[found-1];
  }

  if( (logxmin + (found-1)*step < logx) && found != xdim ){
    
    out = lvecp1*(logx - logxmin - (found-1)*step)/step - lvec*(logx - logxmin - (found)*step)/step;
    // out = pow(10.,out);

  }else if( (logxmin + (found-1)*step > logx) && found != 1){

    out = lvec*(logx - logxmin - (found-2)*step)/step - lvecm1*(logx - logxmin - (found-1)*step)/step;
    //out = pow(10., out);

  }else if( (logxmin + (found-1)*step < logx) && found == xdim ){

    out = lvec*(logx - logxmin - (found-2)*step)/step - lvecm1*(logx - logxmin - (found-1)*step)/step;
    //out = pow(10.,out);

  }else if( (logxmin + (found-1)*step > logx) && found == 1){

    out =  -lvecp1*(logx - logxmin - (found-1)*step)/step + lvec*(logx - 2.*logxmin - 2.*(found-1)*step + logxmin + (found)*step )/step;
    //out = pow(10., out);
  }
  
  if(loglog){
    out = pow(10.,out);
  }


  return out;

}

/*
********************************************************************************
*
* SYNCHROTRON EMISSION COEFFICIENT
*
* elec_spec   - function which defines electrons energy spectrum
* gamma_min   - minimal energy of electrons (Lorentz factor)
* gamma_max   - maximal energy of electrons (Lorentz factor)
* nu          - frequency for which spectrum will be calculated
* B           - value of magnetic field
* prec1       - integration precision for trapezoid integration
* prec2       - integration precision for Gauss-Legendre
*
********************************************************************************
*/
/**
 * @brief Calculates the synchrotron emission coefficient.
 *
 * This function calculates the synchrotron emission coefficient based on the given parameters.
 * It integrates over the electron energy spectrum defined by elec_spec and returns the result.
 *
 * @param elec_spec Function pointer to the electron energy spectrum function.
 * @param gamma_min The minimal energy of electrons (Lorentz factor).
 * @param gamma_max The maximal energy of electrons (Lorentz factor).
 * @param nu The frequency for which the spectrum will be calculated.
 * @param B The value of the magnetic field.
 * @param prec1 The integration precision for trapezoid integration.
 * @param prec2 The integration precision for Gauss-Legendre.
 *
 * @return The synchrotron emission coefficient.
 *
 * @note The function prints error messages and returns 0 if any of the input values are out of range.
 *
 * @see elec_spec
 */
double j_syn(double (*elec_spec)(double), 
             double gamma_min, 
	     double gamma_max, 
	     double nu, 
	     double B, 
	     int prec1,
	     int prec2
	    ) {

      char name[32] = "j_syn";

      if ((gamma_min < 1.0) || 
          (gamma_min > 1.0e+10) ||
	  (gamma_min > gamma_max)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_min'\n", name);
	return 0.0;
      }	
      if ((gamma_max < 1.0) ||
          (gamma_max > 1.0e+10)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_max'\n", name);
	return 0.0;
      }	
      if ((nu < MIN_FREQ) ||
          (nu > MAX_FREQ)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu'\n", name);
	return 0.0;
      }	
      if ((B < 1.0e-10) ||
          (B > 1.0e+5)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'B'\n", name);	
	return 0.0;
      }	
      if ((prec1 < 3) ||
          (prec1 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec1'\n", name);	
	return 0.0;
      }	
      if ((prec2 < 3) ||
          (prec2 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec2'\n", name);	
	return 0.0;
      }	
	    
      double t_min, t_max, U_B, NU_B;      
      
      U_B       = (B * B) / (8.0 * M_PI);   
      NU_B      = (e * B) / (2.0 * M_PI * m_e * c);

      t_min     = nu / (3.0 * gamma_min * gamma_min * NU_B); 
      t_max     = nu / (3.0 * gamma_max * gamma_max * NU_B); 

      int    i, j, n;
      double a, b;
      double stp, val, sum, prt, x1, x2;
      const  int x_dim = prec2;
      double x[x_dim+1], w[x_dim+1];

      a = t_max;
      b = t_min;
      n = prec1;

      // integration over electrons spectrum

      if (a >= b) { 
        fprintf(stderr,"ERROR: a >= b !!!\n");
        return 0.0;
      }	
      
      double tmp;

      tmp = (sqrt(nu/NU_B) * (9.0 * sig_T * c * U_B * C1)) / 
            (24.0 * M_PI * M_PI * NU_B);
	     
      stp = (log10(b) - log10(a)) / (n - 1);
      val = log10(a);
      sum = 0.0;
      
      for (j = 1; j <= n-1; j++) {
         x1  = pow(10.0, val);
         x2  = pow(10.0, val + stp);
	 
	 gauleg(x1, x2, x, w, x_dim);
	 
	 for (prt = 0.0, i = 1; i <= x_dim; i++) 
	    prt += w[i] * tmp * elec_spec(sqrt(nu/(3.0*x[i]*NU_B))) * 
	           pow(x[i], C2-1.5) * exp(-C3*x[i]);
	 
	 sum += prt;
	 
         val += stp;
      }

      return sum;
}

/*
********************************************************************************
*
* ELECTRONS SELF-ABSORPTION COEFFICIENT
*
* elec_spec   - function which defines electrons energy spectrum
* gamma_min   - minimal energy of electrons (Lorentz factor)
* gamma_max   - maximal energy of electrons (Lorentz factor)
* nu          - frequency for which spectrum will be calculated
* B           - value of magnetic field
* prec1       - integration precision for trapezoid integration
* prec2       - integration precision for Gauss-Legendre
*
********************************************************************************
*/
/**
 * @brief Calculate the electrons self-absorption coefficient.
 *
 * This function calculates the self-absorption coefficient of electrons for a given energy spectrum, frequency, and magnetic field.
 * It performs numerical integration to compute the coefficient.
 *
 * @param elec_spec Pointer to the function that defines the electrons energy spectrum
 * @param gamma_min Minimal energy of electrons (Lorentz factor)
 * @param gamma_max Maximal energy of electrons (Lorentz factor)
 * @param nu Frequency for which the spectrum will be calculated
 * @param B Value of the magnetic field
 * @param prec1 Integration precision for the trapezoid integration
 * @param prec2 Integration precision for the Gauss-Legendre method
 * @return The self-absorption coefficient
 *
 * @note Make sure to provide accurate and valid values for the parameters as specified.
 * @note The integration precisions (prec1, prec2) should be between 3 and 1.0e+4.
 * @note The gamma_min parameter should be between 1.0 and 1.0e+10.
 * @note The gamma_max parameter should be between 1.0 and 1.0e+10 and should be greater than gamma_min.
 * @note The nu parameter should be between MIN_FREQ and MAX_FREQ.
 * @note The B parameter should be between 1.0e-10 and 1.0e+5.
 * @note The elec_spec function should be defined to provide the electrons energy spectrum.
 * @note The function uses several constants: MIN_FREQ, MAX_FREQ, e, m_e, c, C2, sig_T, C1, C3.
 *        Make sure these constants are defined in the code.
 */
double k_esa(double (*elec_spec)(double), 
             double gamma_min, 
	     double gamma_max, 
	     double nu, 
	     double B, 
	     int prec1,
	     int prec2
	    ) {

      char name[32] = "k_esa";

      if ((gamma_min < 1.0) || 
          (gamma_min > 1.0e+10) ||
	  (gamma_min > gamma_max)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_min'\n", name);
	return 0.0;
      }	
      if ((gamma_max < 1.0) ||
          (gamma_max > 1.0e+10)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_max'\n", name);
	return 0.0;
      }	
      if ((nu < MIN_FREQ) ||
          (nu > MAX_FREQ)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu'\n", name);
	return 0.0;
      }	
      if ((B < 1.0e-10) ||
          (B > 1.0e+5)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'B'\n", name);	
	return 0.0;
      }	
      if ((prec1 < 3) ||
          (prec1 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec1'\n", name);	
	return 0.0;
      }	
      if ((prec2 < 3) ||
          (prec2 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec2'\n", name);	
	return 0.0;
      }	
	    
      double U_B, NU_B;      
      
      U_B       = (B * B) / (8.0 * M_PI);   
      NU_B      = (e * B) / (2.0 * M_PI * m_e * c);

      int    i, j, n;
      double a, b;
      double stp, val, prt, sum, x1, x2, ff;
      const  int x_dim = prec2;
      double x[x_dim+1], w[x_dim+1];


      a = gamma_min;
      b = gamma_max;
      n = prec1;

      // integration over electrons spectrum

      if (a >= b) { 
        fprintf(stderr,"ERROR: a >= b !!!\n");
        return 0.0;
      }	
      
      double dtmp;
      
      dtmp = - (pow(3.0, 0.5-C2)*sig_T*c*U_B*C1) / 
               (8.0*M_PI*M_PI*m_e*nu*nu*NU_B*NU_B);

      stp = (log10(b) - log10(a)) / (n - 1);
      val = log10(a);
      sum = 0.0;
      
      for (j = 1; j <= n-1; j++) {
         x1  = pow(10.0, val);
         x2  = pow(10.0, val + stp);
	 
	 gauleg(x1, x2, x, w, x_dim);
	 
	 for (prt = 0.0, i = 1; i <= x_dim; i++) {
	    if (x[i] == 1.0) { 
	      ff = 0.0;
  	    } else {
	      ff =  dtmp * elec_spec(x[i]) * pow(nu/(x[i]*x[i]*NU_B), C2) * exp(-C3*nu/(3.0*x[i]*x[i]*NU_B)) *
                                     (
	                              (
	                               -6.0*pow(x[i], 4.0)*NU_B+3.0*NU_B*x[i]*x[i]+6.0*C2*NU_B*pow(x[i], 4.0)
	                               -6.0*C2*NU_B*x[i]*x[i]-2.0*C3*nu*x[i]*x[i]+2.0*C3*nu
             	                      ) / (pow(x[i], 3.0)*(x[i]*x[i]-1.0))
	                             );
	    }  
	    prt += w[i] * ff;
	 }
	 sum += prt;
	 
         val += stp;	 
      }

      return sum;
}


/*
********************************************************************************
*
* SUBROUTINE NECESSARY FOR CALCULATIONS OF THE INVERSE-COMPTON EMISSIVITY
*
********************************************************************************
*/

/**
 * @brief Calculates the inverse-Compton emissivity.
 *
 * This function calculates the inverse-Compton emissivity based on the given parameters.
 *
 * @param ec Energy of the electron
 * @param gg Lorentz factor of the electron
 * @param es Energy of the soft photon
 *
 * @return Computed inverse-Compton emissivity
 */
double Compton_kernel(double ec, double gg, double es) {
      double r_e;
      double k_p;
      
      r_e = e * e / (m_e * c * c);
      k_p = ec / (4.0 * es * gg * (gg - ec));
      
      return ((2.0 * M_PI * r_e * r_e * c) / (gg * gg * es)) *
              (2.0 * k_p * log(k_p) + (1.0 + 2.0 * k_p) * (1.0 - k_p) +
	      (pow(4.0 * es * gg * k_p, 2.0) / (2.0*(1.0 + 4.0 * es * gg * k_p))) *
	      (1.0 - k_p)	     
	     );
}

/*
********************************************************************************
*
* INVERSE-COMPTON EMISSION COEFFICIENT (Inoue & Takahara 1996)
*
* elec_spec   - function which defines electrons energy spectrum
* gamma_min   - minimal energy of electrons (Lorentz factor)
* gamma_max   - maximal energy of electrons (Lorentz factor)
* I_rad       - matrix which contains intensity of radiation field
* nu_rad_min  - minimum frequency of radiation field photons
* nu_rad_max  - maximum frequency of radiation field photons
* nu_rad_dim  - dimension of matrix I_rad
* nu          - frequency for which spectrum will be calculated
* prec1       - integration precision for trapezoid integration
* prec2       - integration precision for Gauss-Legendre
*
********************************************************************************
*/
/**
 * @brief Calculates the inverse Compton emission coefficient.
 *
 * This function calculates the inverse Compton emission coefficient using the Inoue & Takahara 1996 formula. Given the input parameters and the electron energy spectrum, it calculates the emission coefficient for a given frequency.
 *
 * @param elec_spec Function pointer to the electron energy spectrum function. The function should take a double as an argument and return a double.
 * @param gamma_min Minimal energy of electrons (Lorentz factor).
 * @param gamma_max Maximal energy of electrons (Lorentz factor).
 * @param I_rad Array containing the intensity of the radiation field.
 * @param nu_rad_min Minimum frequency of the radiation field photons.
 * @param nu_rad_max Maximum frequency of the radiation field photons.
 * @param nu_rad_dim Dimension of the matrix I_rad.
 * @param nu Frequency for which the emission spectrum will be calculated.
 * @param prec1 Integration precision for the trapezoid integration.
 * @param prec2 Integration precision for the Gauss-Legendre.
 * @return The inverse Compton emission coefficient as a double value.
 *
 * @note Before calling this function, make sure to define the following constants:
 * - MIN_FREQ: Minimum frequency value (double)
 * - MAX_FREQ: Maximum frequency value (double)
 * - h: Planck's constant (double)
 * - m_e: Electron mass (double)
 * - c: Speed of light (double)
 *
 * @note The function uses Gauss-Legendre quadrature for numerical integration.
 *
 * @see The Inoue & Takahara 1996 paper for detailed information on the formula.
 */
double j_com(double (*elec_spec)(double), 
             double gamma_min, 
	     double gamma_max, 
	     double I_rad[], 
             double nu_rad_min, 
	     double nu_rad_max,
	     int    nu_rad_dim,
	     double nu, 
	     int prec1,
	     int prec2
	    ) {
	
      char name[32] = "j_com";

      if ((gamma_min < 1.0) || 
          (gamma_min > 1.0e+10) ||
	  (gamma_min > gamma_max)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_min'\n", name);
	return 0.0;
      }	
      if ((gamma_max < 1.0) ||
          (gamma_max > 1.0e+10)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_max'\n", name);
	return 0.0;
      }	
      if (nu_rad_min < MIN_FREQ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_rad_min'\n", name);
	return 0.0;
      }	
      if (nu_rad_max > MAX_FREQ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_rad_max'\n", name);
	return 0.0;
      }	
      if ((nu_rad_dim < 2.0) ||
          (nu_rad_dim > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_rad_dim'\n", name);
	return 0.0;
      }	      
      if ((nu < MIN_FREQ) ||
          (nu > MAX_FREQ)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu'\n", name);
	return 0.0;
      }	
      if ((prec1 < 3) ||
          (prec1 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec1'\n", name);	
	return 0.0;
      }		
      if ((prec2 < 3) ||
          (prec2 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec2'\n", name);	
	return 0.0;
      }	
	
      int    k;    
      double epss[nu_rad_dim+1];
      double n_p[nu_rad_dim+1];
      double tmp_min, tmp_max, tmp_stp, tmp_val, tmp_cur;
      
      tmp_min = log10(nu_rad_min);
      tmp_max = log10(nu_rad_max);
      tmp_stp = (tmp_max - tmp_min) / (nu_rad_dim); // changed (nu_rad_dim - 1) in (nu_rad_dim) for consistency with the main method
      tmp_val = tmp_min;        
   
      for (k = 1; k <= nu_rad_dim; k++) {
         tmp_cur  = pow(10.0, tmp_val);              
         epss[k]  = tmp_cur * h / (m_e * c * c);
         n_p[k]   = ((4.0 * M_PI) / (h * c * epss[k])) * I_rad[k];
         tmp_val  = tmp_val + tmp_stp;
      }     
 
      int   i, j;
      double epss_cur, epsc_cur, xf, xs, xd, yf, ys, sum;
      double stp, val, prt, x1, x2, ff;
      const  int x_dim = prec2;
      double x[x_dim+1], w[x_dim+1];
      
      
      epsc_cur = nu * h / (m_e * c *c); 
      
      sum = 0.0;
         
      for (k = 1; k <= nu_rad_dim-1; k++) {
         if (epss[k] <= epsc_cur) {
      
   	   xf       = epss[k];
   	   xs       = epss[k+1];
   	   xd       = xs - xf; 
	 	   	 
	   epss_cur = epss[k];	 	   	   
	   
     //discretization of e- lorentz factors
	   stp = (log10(gamma_max) - log10(gamma_min)) / (prec1 - 1);
           val = log10(gamma_min);
           yf  = 0.0;
      
     for (j = 1; j <= prec1 - 1; j++) {
        x1  = pow(10.0, val);
        x2  = pow(10.0, val + stp);	      	      
	      	      
	      gauleg(x1, x2, x, w, x_dim);
	 
	      for (prt = 0.0, i = 1; i <= x_dim; i++) {
	         if (epsc_cur <= ((x[i] * 4.0 * x[i] * epss_cur) / (1.0 + 4.0 * epss_cur * x[i]))) { 
	           ff = elec_spec(x[i]) * Compton_kernel(epsc_cur, x[i], epss_cur);
  	         } else {
	           ff = 0.0; 
	         }  
	         prt += w[i] * ff;
	      }	      
	      yf  += prt;	 
              val += stp;
      }	   
   	   yf  = yf * n_p[k];		 	 
	 
	   epss_cur = epss[k+1];   	
	   
	   stp = (log10(gamma_max) - log10(gamma_min)) / (prec1 - 1);
           val = log10(gamma_min);
           ys  = 0.0;
      
           for (j = 1; j <= prec1-1; j++) {
              x1  = pow(10.0, val);
              x2  = pow(10.0, val + stp);
	      
	      gauleg(x1, x2, x, w, x_dim);
	 
	  for (prt = 0.0, i = 1; i <= x_dim; i++) {
		if (epsc_cur <= ((x[i] * 4.0 * x[i] * epss_cur) / (1.0 + 4.0 * epss_cur * x[i]))) {
		  ff = elec_spec(x[i]) * Compton_kernel(epsc_cur, x[i], epss_cur);
		} else {
		  ff = 0.0;
		}
		prt += w[i] * ff;
	      }
	      ys  += prt;	 
              val += stp;	      
           }	   
   	   ys       = ys * n_p[k+1];         	   
	   
   	   sum      = sum + (xd * (yf + ys) / 2.0);   	
	 }
     }      

     return (sum * h * epsc_cur) / (4.0 * M_PI); 
}


// IR ABSORPTION ACCORDING TO KNEISKE et al. 2002, 2004

const int LDIM = 51;
const int MDIM = 2601; // 51x51;

double Z1[LDIM] = { 0.000,  0.002,  0.004,  0.006,  0.008,  0.010,  0.012, 
                    0.014,  0.016,  0.018,  0.020,  0.022,  0.024,  0.026,
                    0.028,  0.030,  0.032,  0.034,  0.036,  0.038,  0.040,
                    0.042,  0.044,  0.046,  0.048,  0.050,  0.052,  0.054,
                    0.056,  0.058,  0.060,  0.062,  0.064,  0.066,  0.068,
                    0.070,  0.072,  0.074,  0.076,  0.078,  0.080,  0.082,
                    0.084,  0.086,  0.088,  0.090,  0.092,  0.094,  0.096,
                    0.098,  0.100};
                    
double E1[LDIM] = {-1.000, -0.820, -0.640, -0.460, -0.280, -0.100,  0.080,
                    0.260,  0.440,  0.620,  0.800,  0.980,  1.160,  1.340,
                    1.520,  1.700,  1.880,  2.060,  2.240,  2.420,  2.600,
                    2.780,  2.960,  3.140,  3.320,  3.500,  3.680,  3.860,
                    4.040,  4.220,  4.400,  4.580,  4.760,  4.940,  5.120,
                    5.300,  5.480,  5.660,  5.840,  6.020,  6.200,  6.380,
                    6.560,  6.740,  6.920,  7.100,  7.280,  7.460,  7.640,
                    7.820,  8.000};
double T1[MDIM] = {
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.001,   0.002,   0.002,   0.002,   0.002,   0.002,   0.002,   0.002,   0.002,   0.002,   0.002,   0.002,   0.003,   0.003,   0.003,   0.003,   0.003,   0.003,   0.003,   0.003,   0.003,   0.004,   0.004,   0.004,   0.004,   0.004,   0.004, 
  0.000,   0.000,   0.001,   0.001,   0.001,   0.002,   0.002,   0.002,   0.002,   0.003,   0.003,   0.003,   0.004,   0.004,   0.005,   0.005,   0.005,   0.006,   0.006,   0.006,   0.007,   0.007,   0.008,   0.008,   0.008,   0.009,   0.009,   0.010,   0.010,   0.010,   0.011,   0.011,   0.012,   0.012,   0.013,   0.013,   0.013,   0.014,   0.014,   0.015,   0.015,   0.016,   0.016,   0.017,   0.017,   0.018,   0.018,   0.019,   0.019,   0.020,   0.020, 
  0.000,   0.001,   0.002,   0.003,   0.003,   0.004,   0.005,   0.006,   0.007,   0.008,   0.009,   0.010,   0.011,   0.012,   0.013,   0.014,   0.015,   0.015,   0.016,   0.017,   0.018,   0.019,   0.020,   0.021,   0.022,   0.024,   0.025,   0.026,   0.027,   0.028,   0.029,   0.030,   0.031,   0.032,   0.033,   0.034,   0.036,   0.037,   0.038,   0.039,   0.040,   0.041,   0.043,   0.044,   0.045,   0.046,   0.047,   0.049,   0.050,   0.051,   0.053, 
  0.000,   0.002,   0.004,   0.006,   0.008,   0.010,   0.012,   0.014,   0.017,   0.019,   0.021,   0.023,   0.025,   0.027,   0.030,   0.032,   0.034,   0.036,   0.039,   0.041,   0.043,   0.046,   0.048,   0.050,   0.053,   0.055,   0.058,   0.060,   0.063,   0.065,   0.068,   0.070,   0.073,   0.075,   0.078,   0.081,   0.083,   0.086,   0.089,   0.091,   0.094,   0.097,   0.100,   0.103,   0.105,   0.108,   0.111,   0.114,   0.117,   0.120,   0.123, 
  0.000,   0.005,   0.010,   0.014,   0.019,   0.024,   0.029,   0.034,   0.039,   0.044,   0.049,   0.054,   0.060,   0.065,   0.070,   0.075,   0.081,   0.086,   0.092,   0.097,   0.103,   0.108,   0.114,   0.119,   0.125,   0.131,   0.137,   0.142,   0.148,   0.154,   0.160,   0.166,   0.172,   0.178,   0.184,   0.190,   0.196,   0.202,   0.209,   0.215,   0.221,   0.227,   0.234,   0.240,   0.247,   0.254,   0.260,   0.267,   0.274,   0.280,   0.287, 
  0.000,   0.010,   0.021,   0.031,   0.042,   0.053,   0.063,   0.074,   0.085,   0.096,   0.107,   0.118,   0.129,   0.140,   0.152,   0.163,   0.175,   0.186,   0.198,   0.209,   0.221,   0.233,   0.245,   0.257,   0.269,   0.281,   0.293,   0.306,   0.318,   0.330,   0.343,   0.356,   0.368,   0.381,   0.394,   0.407,   0.420,   0.433,   0.446,   0.459,   0.472,   0.486,   0.499,   0.513,   0.526,   0.540,   0.554,   0.568,   0.582,   0.596,   0.610, 
  0.000,   0.020,   0.039,   0.059,   0.079,   0.099,   0.119,   0.139,   0.159,   0.179,   0.200,   0.221,   0.241,   0.262,   0.283,   0.304,   0.325,   0.346,   0.368,   0.389,   0.411,   0.432,   0.454,   0.476,   0.498,   0.520,   0.542,   0.564,   0.587,   0.609,   0.632,   0.655,   0.678,   0.700,   0.723,   0.746,   0.770,   0.793,   0.816,   0.840,   0.864,   0.887,   0.911,   0.935,   0.959,   0.984,   1.008,   1.033,   1.057,   1.082,   1.107, 
  0.000,   0.031,   0.061,   0.092,   0.123,   0.154,   0.185,   0.217,   0.248,   0.280,   0.312,   0.343,   0.375,   0.408,   0.440,   0.472,   0.505,   0.537,   0.570,   0.603,   0.636,   0.669,   0.702,   0.736,   0.769,   0.803,   0.836,   0.870,   0.904,   0.938,   0.972,   1.006,   1.041,   1.075,   1.110,   1.145,   1.180,   1.214,   1.250,   1.285,   1.320,   1.356,   1.391,   1.427,   1.463,   1.499,   1.535,   1.571,   1.607,   1.644,   1.680, 
  0.000,   0.040,   0.081,   0.122,   0.162,   0.203,   0.244,   0.285,   0.326,   0.368,   0.409,   0.451,   0.493,   0.534,   0.576,   0.618,   0.661,   0.703,   0.745,   0.788,   0.830,   0.873,   0.916,   0.959,   1.002,   1.045,   1.088,   1.132,   1.175,   1.219,   1.263,   1.306,   1.351,   1.395,   1.439,   1.483,   1.528,   1.572,   1.617,   1.662,   1.706,   1.752,   1.797,   1.842,   1.887,   1.933,   1.978,   2.024,   2.070,   2.115,   2.161, 
  0.000,   0.048,   0.096,   0.144,   0.192,   0.241,   0.289,   0.338,   0.387,   0.436,   0.485,   0.534,   0.583,   0.633,   0.683,   0.732,   0.782,   0.832,   0.882,   0.932,   0.983,   1.033,   1.084,   1.134,   1.185,   1.236,   1.287,   1.338,   1.389,   1.441,   1.492,   1.544,   1.596,   1.647,   1.699,   1.752,   1.804,   1.856,   1.909,   1.961,   2.014,   2.066,   2.119,   2.173,   2.226,   2.279,   2.332,   2.386,   2.439,   2.493,   2.547, 
  0.000,   0.055,   0.111,   0.166,   0.222,   0.278,   0.334,   0.390,   0.446,   0.503,   0.560,   0.617,   0.674,   0.731,   0.789,   0.846,   0.904,   0.962,   1.020,   1.078,   1.137,   1.195,   1.254,   1.313,   1.372,   1.431,   1.491,   1.550,   1.610,   1.670,   1.730,   1.791,   1.851,   1.912,   1.973,   2.034,   2.095,   2.156,   2.218,   2.280,   2.341,   2.403,   2.466,   2.528,   2.591,   2.654,   2.716,   2.779,   2.842,   2.906,   2.969, 
  0.000,   0.067,   0.134,   0.202,   0.269,   0.337,   0.406,   0.474,   0.543,   0.612,   0.681,   0.750,   0.820,   0.890,   0.960,   1.031,   1.102,   1.173,   1.244,   1.315,   1.387,   1.459,   1.531,   1.604,   1.677,   1.750,   1.823,   1.896,   1.970,   2.044,   2.118,   2.193,   2.267,   2.342,   2.418,   2.493,   2.569,   2.645,   2.721,   2.798,   2.874,   2.951,   3.029,   3.106,   3.184,   3.262,   3.340,   3.419,   3.497,   3.576,   3.655, 
  0.000,   0.084,   0.168,   0.252,   0.337,   0.422,   0.507,   0.593,   0.679,   0.765,   0.852,   0.939,   1.026,   1.114,   1.202,   1.290,   1.379,   1.468,   1.557,   1.647,   1.737,   1.828,   1.918,   2.009,   2.101,   2.193,   2.285,   2.377,   2.470,   2.563,   2.656,   2.750,   2.844,   2.938,   3.033,   3.128,   3.223,   3.319,   3.415,   3.512,   3.609,   3.706,   3.803,   3.901,   3.999,   4.098,   4.196,   4.295,   4.395,   4.495,   4.595, 
  0.000,   0.106,   0.213,   0.320,   0.427,   0.535,   0.644,   0.753,   0.862,   0.972,   1.082,   1.193,   1.305,   1.417,   1.529,   1.642,   1.755,   1.869,   1.984,   2.098,   2.214,   2.330,   2.446,   2.563,   2.680,   2.798,   2.917,   3.035,   3.155,   3.275,   3.395,   3.516,   3.637,   3.759,   3.882,   4.004,   4.128,   4.252,   4.376,   4.502,   4.627,   4.753,   4.880,   5.007,   5.135,   5.263,   5.392,   5.521,   5.651,   5.782,   5.913, 
  0.000,   0.152,   0.306,   0.460,   0.615,   0.772,   0.929,   1.087,   1.247,   1.407,   1.568,   1.731,   1.894,   2.059,   2.225,   2.392,   2.559,   2.729,   2.899,   3.070,   3.243,   3.416,   3.591,   3.767,   3.944,   4.123,   4.302,   4.483,   4.665,   4.849,   5.033,   5.219,   5.406,   5.595,   5.784,   5.975,   6.167,   6.361,   6.556,   6.753,   6.950,   7.150,   7.350,   7.552,   7.755,   7.960,   8.166,   8.373,   8.583,   8.794,   9.006, 
  0.000,   0.337,   0.678,   1.023,   1.371,   1.723,   2.078,   2.437,   2.800,   3.167,   3.538,   3.913,   4.292,   4.674,   5.061,   5.452,   5.847,   6.246,   6.650,   7.057,   7.469,   7.885,   8.306,   8.731,   9.160,   9.595,  10.033,  10.476,  10.924,  11.376,  11.834,  12.297,  12.763,  13.234,  13.711,  14.193,  14.679,  15.171,  15.667,  16.170,  16.677,  17.190,  17.707,  18.231,  18.759,  19.292,  19.832,  20.378,  20.928,  21.484,  22.047, 
  0.000,   0.934,   1.878,   2.831,   3.794,   4.766,   5.748,   6.740,   7.742,   8.754,   9.777,  10.809,  11.852,  12.905,  13.968,  15.042,  16.126,  17.221,  18.327,  19.443,  20.569,  21.707,  22.856,  24.017,  25.186,  26.370,  27.564,  28.765,  29.981,  31.207,  32.449,  33.702,  34.961,  36.229,  37.515,  38.816,  40.127,  41.444,  42.776,  44.125,  45.480,  46.849,  48.230,  49.629,  51.033,  52.445,  53.878,  55.321,  56.770,  58.234,  59.724, 
  0.000,   1.980,   3.975,   5.986,   8.011,  10.052,  12.109,  14.181,  16.269,  18.372,  20.491,  22.627,  24.778,  26.945,  29.127,  31.326,  33.542,  35.774,  38.021,  40.286,  42.564,  44.860,  47.174,  49.503,  51.842,  54.206,  56.585,  58.974,  61.383,  63.805,  66.250,  68.708,  71.185,  73.666,  76.172,  78.708,  81.257,  83.806,  86.378,  88.969,  91.571,  94.202,  96.844,  99.510, 102.179, 104.875, 107.593, 110.319, 113.058, 115.820, 118.613, 
  0.000,   3.144,   6.306,   9.485,  12.681,  15.894,  19.126,  22.377,  25.645,  28.932,  32.235,  35.557,  38.898,  42.257,  45.634,  49.029,  52.442,  55.874,  59.323,  62.792,  66.277,  69.780,  73.306,  76.851,  80.404,  83.987,  87.589,  91.198,  94.835,  98.485, 102.164, 105.856, 109.566, 113.281, 117.029, 120.808, 124.606, 128.391, 132.202, 136.047, 139.898, 143.771, 147.668, 151.591, 155.511, 159.456, 163.436, 167.416, 171.404, 175.425, 179.480, 
  0.000,   4.068,   8.152,  12.255,  16.374,  20.509,  24.663,  28.834,  33.022,  37.228,  41.450,  45.690,  49.948,  54.223,  58.514,  62.821,  67.146,  71.486,  75.844,  80.219,  84.608,  89.015,  93.441,  97.884, 102.332, 106.811, 111.306, 115.805, 120.330, 124.863, 129.426, 133.997, 138.588, 143.178, 147.798, 152.458, 157.130, 161.781, 166.459, 171.168, 175.887, 180.628, 185.380, 190.161, 194.938, 199.748, 204.581, 209.410, 214.253, 219.129, 224.032, 
  0.000,   4.499,   9.012,  13.538,  18.076,  22.627,  27.193,  31.771,  36.363,  40.968,  45.585,  50.215,  54.859,  59.516,  64.186,  68.868,  73.565,  78.273,  82.995,  87.732,  92.477,  97.237, 102.013, 106.803, 111.595, 116.411, 121.242, 126.072, 130.926, 135.783, 140.664, 145.551, 150.453, 155.350, 160.275, 165.231, 170.196, 175.137, 180.098, 185.091, 190.087, 195.095, 200.113, 205.155, 210.189, 215.249, 220.330, 225.399, 230.477, 235.586, 240.717, 
  0.000,   4.406,   8.821,  13.244,  17.677,  22.118,  26.568,  31.027,  35.495,  39.972,  44.456,  48.949,  53.450,  57.961,  62.481,  67.007,  71.544,  76.087,  80.640,  85.201,  89.767,  94.343,  98.929, 103.525, 108.119, 112.731, 117.353, 121.971, 126.608, 131.243, 135.897, 140.553, 145.219, 149.879, 154.559, 159.265, 163.976, 168.660, 173.359, 178.085, 182.813, 187.545, 192.280, 197.038, 201.786, 206.555, 211.334, 216.098, 220.873, 225.674, 230.487, 
  0.000,   3.939,   7.882,  11.831,  15.785,  19.743,  23.707,  27.676,  31.650,  35.630,  39.613,  43.601,  47.595,  51.594,  55.597,  59.606,  63.620,  67.637,  71.660,  75.688,  79.718,  83.754,  87.796,  91.843,  95.888,  99.945, 104.008, 108.066, 112.137, 116.204, 120.285, 124.367, 128.453, 132.533, 136.627, 140.741, 144.856, 148.947, 153.049, 157.171, 161.292, 165.414, 169.536, 173.675, 177.806, 181.951, 186.102, 190.240, 194.383, 198.545, 202.716, 
  0.000,   3.296,   6.595,   9.896,  13.200,  16.506,  19.816,  23.127,  26.442,  29.760,  33.078,  36.400,  39.724,  43.052,  46.382,  49.714,  53.049,  56.385,  59.725,  63.067,  66.409,  69.755,  73.103,  76.456,  79.805,  83.162,  86.523,  89.878,  93.242,  96.601,  99.970, 103.337, 106.707, 110.070, 113.443, 116.829, 120.216, 123.582, 126.955, 130.343, 133.729, 137.115, 140.498, 143.895, 147.284, 150.684, 154.085, 157.473, 160.868, 164.277, 167.691, 
  0.000,   2.635,   5.270,   7.907,  10.545,  13.184,  15.824,  18.465,  21.108,  23.752,  26.395,  29.040,  31.686,  34.334,  36.983,  39.632,  42.283,  44.934,  47.587,  50.241,  52.894,  55.549,  58.205,  60.864,  63.519,  66.180,  68.842,  71.500,  74.164,  76.823,  79.489,  82.153,  84.819,  87.479,  90.146,  92.823,  95.498,  98.157, 100.822, 103.495, 106.168, 108.839, 111.507, 114.185, 116.857, 119.536, 122.214, 124.882, 127.555, 130.237, 132.923, 
  0.000,   2.038,   4.076,   6.115,   8.153,  10.192,  12.232,  14.272,  16.312,  18.353,  20.394,  22.435,  24.476,  26.518,  28.560,  30.603,  32.646,  34.688,  36.732,  38.776,  40.818,  42.862,  44.906,  46.952,  48.994,  51.041,  53.088,  55.131,  57.179,  59.223,  61.271,  63.318,  65.366,  67.409,  69.456,  71.511,  73.564,  75.605,  77.649,  79.700,  81.751,  83.799,  85.845,  87.898,  89.947,  92.000,  94.052,  96.097,  98.144, 100.199, 102.257, 
  0.000,   1.538,   3.076,   4.615,   6.153,   7.691,   9.229,  10.768,  12.306,  13.845,  15.383,  16.921,  18.459,  19.997,  21.536,  23.074,  24.612,  26.150,  27.688,  29.227,  30.764,  32.302,  33.839,  35.378,  36.914,  38.453,  39.991,  41.527,  43.066,  44.602,  46.140,  47.678,  49.215,  50.749,  52.286,  53.828,  55.369,  56.900,  58.434,  59.972,  61.510,  63.046,  64.580,  66.119,  67.654,  69.192,  70.729,  72.260,  73.793,  75.331,  76.870, 
  0.000,   1.132,   2.263,   3.394,   4.524,   5.655,   6.785,   7.915,   9.044,  10.173,  11.302,  12.430,  13.558,  14.686,  15.813,  16.940,  18.067,  19.193,  20.320,  21.445,  22.570,  23.695,  24.819,  25.944,  27.066,  28.190,  29.314,  30.435,  31.558,  32.678,  33.799,  34.920,  36.040,  37.157,  38.275,  39.397,  40.518,  41.631,  42.745,  43.863,  44.979,  46.094,  47.206,  48.322,  49.434,  50.549,  51.662,  52.770,  53.879,  54.992,  56.105, 
  0.000,   0.796,   1.591,   2.385,   3.179,   3.972,   4.765,   5.556,   6.348,   7.138,   7.928,   8.716,   9.505,  10.292,  11.079,  11.865,  12.650,  13.435,  14.219,  15.002,  15.784,  16.565,  17.346,  18.127,  18.905,  19.684,  20.462,  21.239,  22.016,  22.790,  23.565,  24.339,  25.112,  25.884,  26.655,  27.428,  28.200,  28.967,  29.734,  30.502,  31.270,  32.035,  32.799,  33.565,  34.328,  35.093,  35.855,  36.614,  37.373,  38.135,  38.896, 
  0.000,   0.520,   1.039,   1.558,   2.075,   2.592,   3.108,   3.623,   4.138,   4.651,   5.164,   5.675,   6.186,   6.696,   7.205,   7.714,   8.221,   8.728,   9.234,   9.739,  10.243,  10.746,  11.248,  11.750,  12.251,  12.751,  13.250,  13.748,  14.246,  14.741,  15.237,  15.732,  16.226,  16.718,  17.210,  17.703,  18.194,  18.682,  19.170,  19.658,  20.146,  20.631,  21.115,  21.600,  22.084,  22.567,  23.049,  23.528,  24.007,  24.487,  24.967, 
  0.000,   0.311,   0.621,   0.930,   1.239,   1.547,   1.854,   2.160,   2.466,   2.771,   3.075,   3.378,   3.681,   3.982,   4.283,   4.584,   4.883,   5.182,   5.480,   5.778,   6.074,   6.370,   6.665,   6.959,   7.252,   7.545,   7.837,   8.128,   8.419,   8.709,   8.998,   9.286,   9.574,   9.860,  10.146,  10.432,  10.717,  11.000,  11.283,  11.566,  11.848,  12.129,  12.408,  12.688,  12.967,  13.246,  13.523,  13.799,  14.074,  14.350,  14.625, 
  0.000,   0.169,   0.338,   0.507,   0.674,   0.841,   1.008,   1.174,   1.340,   1.505,   1.669,   1.833,   1.997,   2.159,   2.322,   2.484,   2.645,   2.806,   2.966,   3.125,   3.284,   3.443,   3.601,   3.759,   3.916,   4.072,   4.228,   4.384,   4.539,   4.693,   4.847,   5.000,   5.153,   5.305,   5.457,   5.609,   5.760,   5.910,   6.060,   6.209,   6.358,   6.506,   6.653,   6.801,   6.948,   7.094,   7.240,   7.385,   7.529,   7.674,   7.818, 
  0.000,   0.085,   0.170,   0.254,   0.338,   0.422,   0.506,   0.589,   0.672,   0.754,   0.836,   0.918,   1.000,   1.081,   1.162,   1.242,   1.322,   1.402,   1.482,   1.561,   1.640,   1.719,   1.797,   1.875,   1.953,   2.030,   2.107,   2.184,   2.260,   2.336,   2.412,   2.488,   2.563,   2.638,   2.713,   2.787,   2.861,   2.935,   3.008,   3.082,   3.155,   3.227,   3.299,   3.371,   3.443,   3.515,   3.586,   3.657,   3.727,   3.798,   3.868, 
  0.000,   0.040,   0.080,   0.120,   0.160,   0.200,   0.239,   0.278,   0.317,   0.356,   0.395,   0.433,   0.472,   0.510,   0.548,   0.586,   0.624,   0.661,   0.699,   0.736,   0.773,   0.810,   0.846,   0.883,   0.919,   0.956,   0.992,   1.028,   1.063,   1.099,   1.135,   1.170,   1.205,   1.240,   1.275,   1.310,   1.344,   1.379,   1.413,   1.447,   1.481,   1.515,   1.548,   1.582,   1.616,   1.649,   1.682,   1.715,   1.748,   1.780,   1.813, 
  0.000,   0.018,   0.037,   0.055,   0.073,   0.091,   0.109,   0.127,   0.145,   0.162,   0.180,   0.197,   0.215,   0.232,   0.250,   0.267,   0.284,   0.301,   0.318,   0.335,   0.352,   0.368,   0.385,   0.402,   0.418,   0.435,   0.451,   0.467,   0.483,   0.500,   0.516,   0.532,   0.548,   0.563,   0.579,   0.595,   0.611,   0.626,   0.642,   0.657,   0.672,   0.688,   0.703,   0.718,   0.733,   0.748,   0.763,   0.778,   0.793,   0.807,   0.822, 
  0.000,   0.008,   0.016,   0.024,   0.033,   0.041,   0.049,   0.057,   0.065,   0.072,   0.080,   0.088,   0.096,   0.104,   0.111,   0.119,   0.127,   0.134,   0.142,   0.149,   0.157,   0.164,   0.172,   0.179,   0.186,   0.194,   0.201,   0.208,   0.215,   0.223,   0.230,   0.237,   0.244,   0.251,   0.258,   0.265,   0.272,   0.279,   0.286,   0.293,   0.299,   0.306,   0.313,   0.320,   0.326,   0.333,   0.340,   0.346,   0.353,   0.359,   0.366};

double Z2[LDIM] = { 0.000,  0.100,  0.200,  0.300,  0.400,  0.500,  0.600,  
                    0.700,  0.800,  0.900,  1.000,  1.100,  1.200,  1.300,  
                    1.400,  1.500,  1.600,  1.700,  1.800,  1.900,  2.000,  
                    2.100,  2.200,  2.300,  2.400,  2.500,  2.600,  2.700,  
                    2.800,  2.900,  3.000,  3.100,  3.200,  3.300,  3.400,  
                    3.500,  3.600,  3.700,  3.800,  3.900,  4.000,  4.100,  
                    4.200,  4.300,  4.400,  4.500,  4.600,  4.700,  4.800,  
                    4.900,  5.000};

double E2[LDIM] = {-1.000, -0.820, -0.640, -0.460, -0.280, -0.100,  0.080,
                    0.260,  0.440,  0.620,  0.800,  0.980,  1.160,  1.340,
                    1.520,  1.700,  1.880,  2.060,  2.240,  2.420,  2.600,
                    2.780,  2.960,  3.140,  3.320,  3.500,  3.680,  3.860,
                    4.040,  4.220,  4.400,  4.580,  4.760,  4.940,  5.120,
                    5.300,  5.480,  5.660,  5.840,  6.020,  6.200,  6.380,
                    6.560,  6.740,  6.920,  7.100,  7.280,  7.460,  7.640,
                    7.820,  8.000};
                    
double T2[MDIM] = {                    
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.001,   0.002,   0.003,   0.005,   0.007,   0.009,   0.011,   0.014,   0.018,   0.020,   0.021,   0.026,   0.028,   0.030,   0.033,   0.036,   0.038,   0.041,   0.044,   0.044, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.001,   0.003,   0.007,   0.012,   0.018,   0.026,   0.031,   0.041,   0.047,   0.056,   0.064,   0.075,   0.082,   0.093,   0.101,   0.110,   0.122,   0.131,   0.137,   0.142,   0.156,   0.161,   0.170,   0.176,   0.182,   0.188,   0.197,   0.198,   0.200,   0.209,   0.217,   0.221,   0.219,   0.224, 
  0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.002,   0.007,   0.016,   0.031,   0.047,   0.063,   0.082,   0.105,   0.122,   0.148,   0.165,   0.195,   0.215,   0.234,   0.263,   0.280,   0.301,   0.326,   0.346,   0.362,   0.381,   0.401,   0.423,   0.437,   0.446,   0.468,   0.489,   0.508,   0.515,   0.530,   0.543,   0.557,   0.573,   0.580,   0.599,   0.609,   0.613,   0.616,   0.619,   0.640,   0.653,   0.660,   0.662, 
  0.000,   0.000,   0.000,   0.002,   0.006,   0.016,   0.032,   0.057,   0.091,   0.131,   0.178,   0.229,   0.278,   0.325,   0.375,   0.416,   0.463,   0.509,   0.558,   0.598,   0.646,   0.687,   0.732,   0.778,   0.814,   0.854,   0.890,   0.930,   0.970,   0.999,   1.036,   1.067,   1.089,   1.119,   1.160,   1.180,   1.213,   1.234,   1.257,   1.278,   1.294,   1.311,   1.333,   1.340,   1.360,   1.377,   1.398,   1.415,   1.418,   1.426,   1.432, 
  0.000,   0.004,   0.014,   0.032,   0.063,   0.104,   0.158,   0.222,   0.297,   0.382,   0.475,   0.572,   0.670,   0.763,   0.854,   0.946,   1.036,   1.125,   1.212,   1.296,   1.377,   1.460,   1.540,   1.611,   1.687,   1.755,   1.828,   1.890,   1.958,   2.021,   2.069,   2.126,   2.178,   2.231,   2.286,   2.325,   2.363,   2.400,   2.436,   2.475,   2.503,   2.539,   2.569,   2.587,   2.602,   2.618,   2.658,   2.676,   2.698,   2.710,   2.725, 
  0.000,   0.021,   0.055,   0.107,   0.178,   0.267,   0.377,   0.509,   0.661,   0.826,   1.002,   1.184,   1.365,   1.545,   1.721,   1.894,   2.063,   2.229,   2.388,   2.542,   2.693,   2.838,   2.975,   3.107,   3.229,   3.348,   3.466,   3.576,   3.676,   3.772,   3.863,   3.948,   4.019,   4.107,   4.190,   4.266,   4.325,   4.376,   4.439,   4.491,   4.542,   4.590,   4.638,   4.688,   4.721,   4.751,   4.786,   4.817,   4.841,   4.865,   4.889, 
  0.000,   0.053,   0.132,   0.244,   0.394,   0.584,   0.815,   1.086,   1.388,   1.711,   2.048,   2.390,   2.729,   3.053,   3.371,   3.666,   3.958,   4.242,   4.513,   4.775,   5.020,   5.260,   5.492,   5.709,   5.916,   6.113,   6.308,   6.483,   6.652,   6.811,   6.969,   7.119,   7.250,   7.379,   7.505,   7.618,   7.720,   7.822,   7.921,   8.010,   8.098,   8.179,   8.251,   8.307,   8.363,   8.414,   8.468,   8.517,   8.564,   8.599,   8.624, 
  0.000,   0.124,   0.304,   0.554,   0.878,   1.273,   1.731,   2.242,   2.797,   3.379,   3.973,   4.569,   5.153,   5.716,   6.263,   6.776,   7.275,   7.753,   8.208,   8.641,   9.054,   9.448,   9.819,  10.173,  10.507,  10.822,  11.116,  11.394,  11.655,  11.897,  12.128,  12.349,  12.541,  12.730,  12.908,  13.072,  13.218,  13.351,  13.480,  13.601,  13.716,  13.829,  13.932,  14.005,  14.073,  14.136,  14.203,  14.270,  14.328,  14.373,  14.410, 
  0.000,   0.287,   0.686,   1.207,   1.855,   2.616,   3.469,   4.399,   5.380,   6.387,   7.397,   8.391,   9.346,  10.253,  11.114,  11.910,  12.671,  13.384,  14.054,  14.686,  15.273,  15.828,  16.346,  16.829,  17.289,  17.710,  18.106,  18.472,  18.817,  19.142,  19.445,  19.732,  19.979,  20.220,  20.459,  20.676,  20.873,  21.055,  21.238,  21.407,  21.565,  21.726,  21.876,  21.993,  22.095,  22.191,  22.298,  22.410,  22.508,  22.590,  22.657, 
  0.000,   0.610,   1.407,   2.396,   3.568,   4.886,   6.311,   7.798,   9.313,  10.818,  12.284,  13.689,  15.012,  16.247,  17.404,  18.463,  19.467,  20.394,  21.260,  22.069,  22.824,  23.533,  24.203,  24.832,  25.420,  25.968,  26.491,  26.988,  27.460,  27.900,  28.315,  28.705,  29.075,  29.446,  29.812,  30.152,  30.463,  30.746,  31.031,  31.311,  31.582,  31.852,  32.103,  32.325,  32.528,  32.719,  32.912,  33.100,  33.271,  33.411,  33.539, 
  0.000,   1.106,   2.469,   4.072,   5.875,   7.808,   9.812,  11.834,  13.836,  15.781,  17.653,  19.431,  21.106,  22.660,  24.122,  25.461,  26.731,  27.924,  29.039,  30.095,  31.097,  32.055,  32.969,  33.842,  34.678,  35.476,  36.251,  36.999,  37.715,  38.402,  39.069,  39.714,  40.331,  40.941,  41.539,  42.108,  42.637,  43.134,  43.620,  44.088,  44.538,  44.972,  45.374,  45.738,  46.074,  46.384,  46.686,  46.974,  47.237,  47.463,  47.659, 
  0.000,   1.677,   3.626,   5.798,   8.135,  10.574,  13.067,  15.567,  18.039,  20.457,  22.794,  25.026,  27.141,  29.132,  31.022,  32.793,  34.506,  36.139,  37.710,  39.222,  40.684,  42.100,  43.467,  44.797,  46.074,  47.301,  48.507,  49.676,  50.793,  51.857,  52.871,  53.862,  54.814,  55.743,  56.639,  57.499,  58.295,  59.034,  59.746,  60.433,  61.093,  61.724,  62.308,  62.849,  63.356,  63.819,  64.262,  64.668,  65.044,  65.375,  65.680, 
  0.000,   2.160,   4.576,   7.212,  10.022,  12.954,  15.980,  19.065,  22.176,  25.273,  28.315,  31.270,  34.112,  36.825,  39.451,  41.943,  44.380,  46.720,  48.989,  51.189,  53.317,  55.387,  57.373,  59.304,  61.169,  62.965,  64.688,  66.365,  67.965,  69.490,  70.962,  72.402,  73.779,  75.107,  76.389,  77.621,  78.783,  79.886,  80.958,  81.975,  82.949,  83.877,  84.742,  85.560,  86.326,  87.041,  87.720,  88.354,  88.947,  89.479,  89.954, 
  0.000,   2.543,   5.375,   8.487,  11.855,  15.458,  19.275,  23.262,  27.350,  31.481,  35.585,  39.599,  43.494,  47.235,  50.869,  54.325,  57.722,  60.995,  64.177,  67.267,  70.269,  73.198,  76.029,  78.771,  81.432,  84.017,  86.518,  88.942,  91.276,  93.529,  95.719,  97.856,  99.905, 101.899, 103.847, 105.738, 107.527, 109.228, 110.884, 112.481, 114.018, 115.507, 116.921, 118.279, 119.573, 120.795, 121.980, 123.119, 124.230, 125.281, 126.279, 
  0.000,   2.967,   6.367,  10.207,  14.469,  19.112,  24.092,  29.343,  34.783,  40.338,  45.898,  51.380,  56.728,  61.883,  66.922,  71.747,  76.510,  81.139,  85.672,  90.105,  94.451,  98.726, 102.899, 106.995, 111.007, 114.953, 118.829, 122.686, 126.504, 130.297, 134.099, 137.956, 141.807, 145.679, 149.638, 153.663, 157.724, 161.838, 166.035, 170.336, 174.695, 179.081, 183.519, 188.018, 192.487, 196.929, 201.357, 205.748, 210.043, 214.168, 218.122, 
  0.000,   3.652,   7.934,  12.838,  18.338,  24.394,  30.957,  37.973,  45.339,  52.962,  60.694,  68.412,  76.029,  83.493,  90.944,  98.285, 105.791, 113.399, 121.220, 129.278, 137.661, 146.390, 155.491, 164.999, 174.964, 185.306, 196.097, 207.316, 218.930, 231.036, 243.619, 256.731, 270.282, 284.332, 298.845, 313.819, 329.094, 344.730, 360.666, 376.692, 392.633, 408.478, 424.043, 439.447, 454.342, 468.691, 482.380, 495.620, 508.014, 519.537, 530.318, 
  0.000,   4.592,  10.043,  16.371,  23.585,  31.666,  40.671,  50.727,  61.968,  74.536,  88.417, 103.554, 119.844, 137.093, 155.535, 174.959, 195.973, 218.400, 242.629, 268.653, 296.553, 326.564, 358.354, 391.683, 426.357, 462.283, 499.028, 536.608, 574.631, 613.333, 652.155, 691.635, 730.796, 769.902, 809.203, 847.653, 885.245, 922.625, 958.747, 993.512, 1027.036, 1059.378, 1089.958, 1118.639, 1146.002, 1171.957, 1195.691, 1217.644, 1238.596, 1257.107, 1274.110, 
  0.000,   5.911,  13.202,  22.158,  33.239,  47.163,  65.028,  88.085, 117.450, 153.809, 196.898, 245.956, 299.586, 356.469, 416.804, 479.116, 545.189, 613.522, 685.139, 758.936, 835.181, 913.127, 992.404, 1071.295, 1150.584, 1229.543, 1306.549, 1382.978, 1458.044, 1530.251, 1601.541, 1671.298, 1738.967, 1803.260, 1865.427, 1927.034, 1984.922, 2040.286, 2093.451, 2143.816, 2191.080, 2235.514, 2278.534, 2318.042, 2354.268, 2387.500, 2418.925, 2448.089, 2475.510, 2499.783, 2521.147, 
  0.000,   9.009,  21.905,  40.924,  69.354, 110.389, 166.619, 239.642, 329.369, 434.502, 551.374, 676.077, 805.074, 935.130, 1068.127, 1198.761, 1331.921, 1464.838, 1597.544, 1729.686, 1861.429, 1991.163, 2117.745, 2242.157, 2363.210, 2480.272, 2592.943, 2700.797, 2805.663, 2905.764, 3002.492, 3095.326, 3181.724, 3266.074, 3346.168, 3422.500, 3494.678, 3562.859, 3626.573, 3686.188, 3743.668, 3796.101, 3843.547, 3887.201, 3929.151, 3967.500, 4002.189, 4035.082, 4062.839, 4087.740, 4110.684, 
  0.000,  22.030,  58.345, 114.949, 196.895, 305.804, 441.434, 601.890, 783.661, 981.775, 1189.916, 1402.180, 1613.917, 1821.036, 2026.358, 2223.118, 2418.201, 2607.655, 2792.789, 2972.401, 3145.981, 3315.884, 3478.184, 3632.535, 3782.126, 3923.559, 4057.134, 4186.818, 4307.988, 4421.621, 4529.106, 4635.345, 4731.974, 4821.854, 4909.586, 4991.282, 5066.075, 5135.250, 5202.233, 5262.064, 5317.433, 5370.475, 5420.160, 5463.793, 5503.402, 5539.215, 5571.581, 5603.335, 5631.482, 5655.009, 5674.906, 
  0.000,  59.652, 150.250, 276.051, 438.318, 634.789, 861.566, 1114.402, 1386.498, 1671.122, 1959.655, 2245.073, 2521.427, 2786.160, 3042.478, 3282.021, 3517.607, 3739.585, 3953.734, 4157.298, 4352.503, 4536.787, 4711.166, 4876.395, 5032.541, 5180.596, 5317.014, 5446.549, 5567.992, 5681.536, 5787.791, 5888.478, 5981.157, 6067.212, 6148.353, 6224.281, 6295.156, 6358.479, 6418.829, 6475.290, 6527.557, 6573.733, 6615.369, 6652.971, 6688.016, 6719.824, 6748.950, 6774.913, 6799.349, 6821.378, 6840.438, 
  0.000, 118.494, 278.927, 481.595, 724.009, 1000.130, 1303.285, 1626.507, 1961.350, 2300.040, 2633.008, 2953.707, 3258.802, 3544.348, 3817.979, 4069.297, 4310.603, 4537.821, 4752.204, 4951.775, 5139.664, 5319.047, 5484.217, 5638.147, 5784.602, 5920.299, 6045.535, 6163.309, 6272.021, 6371.124, 6465.534, 6555.511, 6635.531, 6710.660, 6781.522, 6847.047, 6906.352, 6960.062, 7010.595, 7056.398, 7099.781, 7138.903, 7174.945, 7206.892, 7235.127, 7259.714, 7283.960, 7306.937, 7326.664, 7343.925, 7359.449, 
  0.000, 179.343, 403.482, 668.360, 966.825, 1290.621, 1631.734, 1983.335, 2336.308, 2684.140, 3018.321, 3334.500, 3630.858, 3904.786, 4162.600, 4397.559, 4620.018, 4826.435, 5018.161, 5197.133, 5363.299, 5518.294, 5663.241, 5796.275, 5921.378, 6036.505, 6143.151, 6241.814, 6331.801, 6415.154, 6493.183, 6567.396, 6632.354, 6693.691, 6750.846, 6803.575, 6852.335, 6896.752, 6938.549, 6974.947, 7009.145, 7040.074, 7067.378, 7092.088, 7114.832, 7133.886, 7152.984, 7171.624, 7189.091, 7203.260, 7214.619, 
  0.000, 223.905, 486.639, 780.145, 1096.240, 1426.756, 1764.579, 2103.032, 2435.165, 2755.930, 3058.755, 3342.148, 3604.041, 3844.103, 4068.533, 4268.962, 4458.561, 4632.823, 4794.785, 4942.437, 5079.423, 5208.149, 5325.903, 5434.681, 5536.922, 5630.042, 5714.156, 5793.381, 5865.615, 5931.717, 5993.094, 6052.480, 6103.568, 6150.938, 6196.689, 6239.087, 6276.180, 6308.829, 6341.020, 6369.706, 6396.636, 6420.688, 6441.667, 6460.185, 6477.395, 6493.233, 6507.595, 6521.445, 6534.891, 6545.570, 6554.562, 
  0.000, 240.622, 508.819, 796.089, 1095.344, 1399.919, 1703.633, 2001.777, 2289.522, 2563.372, 2819.096, 3056.468, 3273.865, 3471.321, 3654.232, 3817.822, 3971.044, 4110.555, 4239.157, 4357.532, 4466.750, 4567.834, 4660.479, 4745.954, 4824.617, 4897.160, 4963.151, 5024.792, 5080.570, 5131.459, 5178.210, 5223.147, 5263.375, 5299.957, 5334.083, 5365.871, 5394.335, 5420.211, 5444.386, 5465.703, 5485.696, 5503.714, 5519.862, 5533.871, 5546.482, 5558.239, 5569.310, 5579.874, 5589.993, 5598.566, 5605.766, 
  0.000, 230.425, 477.497, 734.078, 995.078, 1255.276, 1510.209, 1757.013, 1992.392, 2214.309, 2420.055, 2609.483, 2782.140, 2938.057, 3081.465, 3209.623, 3328.595, 3437.413, 3537.528, 3628.418, 3711.829, 3789.643, 3860.361, 3925.063, 3985.562, 4041.236, 4091.024, 4137.058, 4179.067, 4217.457, 4252.827, 4286.727, 4316.502, 4343.948, 4369.604, 4393.746, 4415.326, 4434.192, 4451.867, 4467.756, 4482.504, 4495.454, 4507.267, 4517.729, 4527.243, 4535.838, 4544.152, 4552.424, 4559.841, 4565.887, 4570.890, 
  0.000, 202.681, 413.887, 628.576, 843.212, 1053.882, 1257.986, 1453.662, 1638.800, 1812.328, 1972.309, 2118.488, 2251.346, 2370.954, 2480.718, 2578.164, 2668.782, 2751.280, 2826.663, 2895.299, 2958.374, 3016.486, 3069.574, 3118.214, 3163.258, 3203.940, 3240.938, 3275.439, 3306.513, 3335.046, 3361.552, 3386.303, 3407.941, 3428.204, 3447.177, 3464.437, 3479.884, 3493.665, 3506.807, 3518.197, 3528.676, 3537.935, 3546.198, 3553.652, 3560.503, 3566.669, 3572.549, 3578.374, 3583.826, 3587.953, 3591.393, 
  0.000, 167.672, 339.009, 510.649, 679.928, 844.437, 1002.568, 1153.087, 1294.873, 1426.964, 1548.102, 1658.490, 1758.473, 1848.227, 1930.501, 2003.217, 2070.778, 2131.868, 2187.917, 2238.535, 2284.667, 2327.384, 2366.315, 2401.557, 2434.226, 2464.240, 2491.030, 2515.753, 2537.950, 2557.996, 2576.503, 2594.125, 2609.473, 2623.305, 2636.373, 2648.479, 2659.148, 2668.310, 2676.819, 2684.441, 2691.551, 2697.816, 2703.416, 2708.291, 2712.807, 2716.778, 2720.493, 2724.093, 2727.433, 2730.024, 2732.022, 
  0.000, 132.916, 267.007, 399.770, 529.464, 654.698, 774.278, 887.563, 993.712, 1092.009, 1181.721, 1263.300, 1336.798, 1402.501, 1462.408, 1515.257, 1563.981, 1607.879, 1647.603, 1683.403, 1716.156, 1746.100, 1772.958, 1797.455, 1819.896, 1839.977, 1857.870, 1874.344, 1888.938, 1902.108, 1914.228, 1925.560, 1935.252, 1944.066, 1952.271, 1959.681, 1966.154, 1971.705, 1976.940, 1981.504, 1985.686, 1989.340, 1992.627, 1995.382, 1997.820, 1999.950, 2002.034, 2003.986, 2005.820, 2007.309, 2008.568, 
  0.000, 102.254, 204.423, 304.611, 401.832, 495.029, 583.364, 666.410, 743.425, 814.174, 878.273, 935.997, 987.565, 1033.272, 1074.528, 1110.612, 1143.470, 1172.766, 1199.214, 1222.729, 1243.674, 1262.755, 1279.814, 1294.948, 1308.732, 1321.029, 1331.865, 1341.678, 1350.293, 1357.935, 1364.801, 1371.326, 1376.822, 1381.752, 1386.240, 1390.301, 1393.874, 1396.846, 1399.637, 1402.039, 1404.177, 1406.027, 1407.659, 1409.011, 1410.209, 1411.200, 1412.325, 1413.441, 1414.562, 1415.391, 1416.009, 
  0.000,  76.870, 152.822, 226.560, 297.319, 364.232, 426.760, 484.591, 537.404, 585.214, 627.839, 665.663, 698.986, 728.122, 754.067, 776.400, 796.516, 814.238, 829.861, 843.546, 855.744, 866.633, 876.154, 884.605, 892.252, 898.906, 904.695, 909.881, 914.314, 918.308, 921.931, 925.277, 927.992, 930.471, 932.790, 934.870, 936.627, 938.046, 939.402, 940.546, 941.605, 942.458, 943.183, 943.817, 944.361, 944.873, 945.458, 946.067, 946.674, 947.127, 947.480, 
  0.000,  56.108, 110.568, 162.481, 211.195, 256.230, 297.341, 334.503, 367.737, 397.178, 422.901, 445.333, 464.760, 481.476, 496.131, 508.586, 519.606, 529.164, 537.542, 544.773, 551.044, 556.629, 561.508, 565.781, 569.562, 572.858, 575.697, 578.274, 580.453, 582.388, 584.071, 585.670, 586.985, 588.180, 589.257, 590.213, 591.033, 591.687, 592.292, 592.793, 593.266, 593.638, 593.950, 594.243, 594.533, 594.782, 595.081, 595.379, 595.664, 595.878, 596.027, 
  0.000,  38.901,  75.612, 109.542, 140.419, 168.142, 192.747, 214.435, 233.370, 249.758, 263.800, 275.839, 286.096, 294.799, 302.318, 308.647, 314.168, 318.912, 322.993, 326.495, 329.535, 332.234, 334.519, 336.531, 338.334, 339.905, 341.227, 342.406, 343.374, 344.238, 345.026, 345.773, 346.373, 346.903, 347.393, 347.833, 348.205, 348.500, 348.760, 348.979, 349.194, 349.371, 349.524, 349.662, 349.779, 349.894, 350.031, 350.169, 350.299, 350.391, 350.464, 
  0.000,  24.974,  47.678,  67.898,  85.693, 101.188, 114.579, 126.096, 135.924, 144.276, 151.318, 157.274, 162.286, 166.503, 170.091, 173.104, 175.700, 177.911, 179.814, 181.445, 182.838, 184.058, 185.117, 186.041, 186.850, 187.544, 188.129, 188.657, 189.108, 189.515, 189.860, 190.189, 190.463, 190.709, 190.927, 191.117, 191.279, 191.409, 191.530, 191.628, 191.716, 191.786, 191.844, 191.900, 191.958, 192.009, 192.070, 192.132, 192.193, 192.238, 192.275, 
  0.000,  14.631,  27.400,  38.373,  47.738,  55.683,  62.394,  68.053,  72.808,  76.795,  80.122,  82.909,  85.236,  87.179,  88.824,  90.199,  91.378,  92.379,  93.237,  93.966,  94.593,  95.148,  95.615,  96.018,  96.379,  96.690,  96.956,  97.197,  97.393,  97.566,  97.717,  97.863,  97.984,  98.090,  98.185,  98.268,  98.339,  98.397,  98.450,  98.491,  98.530,  98.562,  98.591,  98.618,  98.643,  98.666,  98.694,  98.721,  98.746,  98.765,  98.781, 
  0.000,   7.822,  14.404,  19.899,  24.480,  28.292,  31.461,  34.100,  36.295,  38.122,  39.636,  40.897,  41.946,  42.818,  43.553,  44.170,  44.696,  45.137,  45.516,  45.845,  46.123,  46.363,  46.571,  46.751,  46.910,  47.047,  47.161,  47.263,  47.348,  47.429,  47.497,  47.562,  47.615,  47.663,  47.705,  47.742,  47.774,  47.800,  47.823,  47.841,  47.858,  47.872,  47.883,  47.893,  47.904,  47.913,  47.924,  47.936,  47.947,  47.956,  47.963, 
  0.000,   3.871,   7.043,   9.639,  11.769,  13.521,  14.963,  16.156,  17.143,  17.961,  18.636,  19.197,  19.661,  20.048,  20.373,  20.646,  20.877,  21.073,  21.240,  21.382,  21.504,  21.612,  21.703,  21.781,  21.850,  21.909,  21.959,  22.006,  22.045,  22.080,  22.109,  22.138,  22.161,  22.182,  22.200,  22.215,  22.229,  22.241,  22.251,  22.259,  22.265,  22.271,  22.276,  22.280,  22.285,  22.289,  22.295,  22.300,  22.305,  22.309,  22.313, 
  0.000,   1.814,   3.276,   4.459,   5.420,   6.205,   6.849,   7.379,   7.816,   8.178,   8.476,   8.723,   8.927,   9.097,   9.240,   9.360,   9.462,   9.546,   9.620,   9.684,   9.738,   9.784,   9.824,   9.858,   9.888,   9.915,   9.937,   9.957,   9.973,   9.989,  10.002,  10.015,  10.025,  10.034,  10.042,  10.049,  10.056,  10.061,  10.065,  10.069,  10.072,  10.074,  10.076,  10.078,  10.080,  10.082,  10.085,  10.088,  10.090,  10.092,  10.093, 
  0.000,   0.823,   1.479,   2.007,   2.433,   2.781,   3.065,   3.299,   3.491,   3.650,   3.781,   3.889,   3.979,   4.054,   4.117,   4.169,   4.214,   4.252,   4.284,   4.311,   4.335,   4.356,   4.373,   4.388,   4.401,   4.413,   4.422,   4.431,   4.439,   4.445,   4.451,   4.456,   4.461,   4.465,   4.468,   4.471,   4.474,   4.476,   4.478,   4.479,   4.480,   4.481,   4.482,   4.483,   4.484,   4.485,   4.486,   4.487,   4.488,   4.488,   4.489, 
  0.000,   0.366,   0.657,   0.890,   1.078,   1.231,   1.356,   1.458,   1.543,   1.613,   1.670,   1.718,   1.757,   1.790,   1.817,   1.840,   1.860,   1.876,   1.890,   1.902,   1.912,   1.921,   1.929,   1.935,   1.941,   1.946,   1.950,   1.953,   1.957,   1.960,   1.962,   1.964,   1.966,   1.968,   1.969,   1.971,   1.972,   1.973,   1.974,   1.974,   1.975,   1.976,   1.976,   1.977,   1.977,   1.977,   1.978,   1.978,   1.979,   1.979,   1.979};

double Z_MIN = 0.0001;
double Z_MAX = 0.0999;
double E_MIN = log10(500);    // 500 GeV
double E_MAX = log10(5.0e+4); //  50 TeV

const  int IDIM = 100;

/**
 * @brief Interpolates a value using bilinear interpolation.
 *
 * This function computes the interpolated value of a given point (ee, zz) using bilinear interpolation.
 * It checks if the given point falls within the range of the provided data arrays E1 and Z1.
 * If the point is within the range, it performs bilinear interpolation using the formula:
 * t = a*ee*zz + b*zz + c*ee + d
 *
 * @param ee The x-coordinate of the point to be interpolated.
 * @param zz The y-coordinate of the point to be interpolated.
 * @param print_flag A flag indicating whether to print debug information.
 *                   If set to 1, debug information will be printed to stderr.
 *                   If set to 0, no debug information will be printed.
 *
 * @return The interpolated value of the given point (ee, zz).
 *         If the point is outside the range of the data arrays, it returns 0.
 */
double Interpolate1(double ee, double zz, int print_flag) {
      int ei = 1, zi = 1;
      int col, row;
      
      double x1=0, y1=0, x2=0, y2=0;
      double z11=0, z12=0, z21=0, z22=0;
      double a, b, c, d, t;

      if ((ee >= E1[0]) && (ee <= E1[LDIM-1]) &&
          (zz >= Z1[0]) && (zz <= Z1[LDIM-1])) {

        if (print_flag) fprintf(stderr, "Input: E=%7.5f z=%7.5f\n\n", ee, zz);
      
        for (col=0; col <= LDIM-2; col++) {
           if ((zz >= Z1[col]) && (zz <= Z1[col+1])) {
    	     zi = col;
	     break;
	   }
        }
            
        for (row=0; row <= LDIM-2; row++) {
           if ((ee >= E1[row]) && (ee <= E1[row+1])) {
	     ei = row;
	     break;
	   }
        }
            
        if (print_flag) fprintf(stderr, "Found: E_min=%7.5f E_max=%7.5f\n", E1[ei], E1[ei+1]);
        if (print_flag) fprintf(stderr, "Found: z_min=%7.5f z_max=%7.5f\n\n", Z1[zi], Z1[zi+1]);
      
        x1  = Z1[zi];
        x2  = Z1[zi+1];
        y1  = E1[ei];
        y2  = E1[ei+1];
      
        z11 = T1[LDIM*ei+zi];           // [ei][zi];
        z12 = T1[LDIM*(ei+1)+zi];       // [ei+1][zi];
        z21 = T1[LDIM*ei+zi+1];         // [ei][zi+1];
        z22 = T1[LDIM*(ei+1)+zi+1];     // [ei+1][zi+1];
      
        a = (z11 - z12 - z21 + z22)             / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        b = (y1*z12 - y2*z11 - y1*z22 + y2*z21) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        c = (x1*z21 - x2*z11 - x1*z22 + x2*z12) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        d = (x1*y1*z22 - x1*y2*z21 - x2*y1*z12 + x2*y2*z11) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
      
        // the function
        t = a*ee*zz + b*zz + c*ee + d;
      } else {
        t = 0;
        //fprintf(stderr, "Incorrect range:\n");
        //fprintf(stderr, "E_min=%e  E=%e  E_max=%e\n", E1[0], ee, E1[LDIM-1]);
        //fprintf(stderr, "Z_min=%e  Z=%e  Z_max=%e\n", Z1[0], zz, Z1[LDIM-1]);
      }
      
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y1=%7.5f)                  (x2=%7.5f, y1=%7.5f)\n", x1, y1, x2, y1);
      if (print_flag) fprintf(stderr, "                    (x=%7.5f, y=%7.5f)\n", zz, ee);
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y2=%7.5f)                  (x2=%7.5f, y2=%7.5f)\n\n", x1, y2, x2, y2);
     
      if (print_flag) fprintf(stderr, "(z11=%7.3f)           (z21=%7.3f)\n", z11, z21);
      if (print_flag) fprintf(stderr, "             (z=%7.3f)\n", t);
      if (print_flag) fprintf(stderr, "(z12=%7.3f)           (z22=%7.3f)\n", z12, z22);
            
      return t;
}

/**
 * Interpolates a value using a bilinear interpolation method.
 *
 * This function takes an coordinate (ee, zz) and returns the interpolated value based on a bilinear interpolation method.
 * The function checks if the coordinate lies within the valid range of the given datasets. If the coordinate is within the range, the interpolation is performed. Otherwise, the function returns 0.
 *
 * @param ee The E coordinate.
 * @param zz The Z coordinate.
 * @param print_flag Flag to indicate whether to print debug information.
 * @return The interpolated value.
 */
double Interpolate2(double ee, double zz, int print_flag) {
      int ei = 1, zi = 1;
      int col, row;
      
      double x1=0, y1=0, x2=0, y2=0;
      double z11=0, z12=0, z21=0, z22=0;
      double a, b, c, d, t;

      if ((ee >= E2[0]) && (ee <= E2[LDIM-1]) &&
          (zz >= Z2[0]) && (zz <= Z2[LDIM-1])) {

        if (print_flag) fprintf(stderr, "Input: E=%7.5f z=%7.5f\n\n", ee, zz);
      
        for (col=0; col <= LDIM-2; col++) {
           if ((zz >= Z2[col]) && (zz <= Z2[col+1])) {
    	     zi = col;
	     break;
	   }
        }
            
        for (row=0; row <= LDIM-2; row++) {
           if ((ee >= E2[row]) && (ee <= E2[row+1])) {
	     ei = row;
	     break;
	   }
        }
            
        if (print_flag) fprintf(stderr, "Found: E_min=%7.5f E_max=%7.5f\n", E2[ei], E2[ei+1]);
        if (print_flag) fprintf(stderr, "Found: z_min=%7.5f z_max=%7.5f\n\n", Z2[zi], Z2[zi+1]);
      
        x1  = Z2[zi];
        x2  = Z2[zi+1];
        y1  = E2[ei];
        y2  = E2[ei+1];
      
        z11 = T2[LDIM*ei+zi];           // [ei][zi];
        z12 = T2[LDIM*(ei+1)+zi];       // [ei+1][zi];
        z21 = T2[LDIM*ei+zi+1];         // [ei][zi+1];
        z22 = T2[LDIM*(ei+1)+zi+1];     // [ei+1][zi+1];
      
        a = (z11 - z12 - z21 + z22)             / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        b = (y1*z12 - y2*z11 - y1*z22 + y2*z21) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        c = (x1*z21 - x2*z11 - x1*z22 + x2*z12) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        d = (x1*y1*z22 - x1*y2*z21 - x2*y1*z12 + x2*y2*z11) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
      
        // the function
        t = a*ee*zz + b*zz + c*ee + d;
      } else {
        t = 0;
        //fprintf(stderr, "Incorrect range:\n");
        //fprintf(stderr, "E_min=%e  E=%e  E_max=%e\n", E2[0], ee, E2[LDIM-1]);
        //fprintf(stderr, "Z_min=%e  Z=%e  Z_max=%e\n", Z2[0], zz, Z2[LDIM-1]);
      }
      
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y1=%7.5f)                  (x2=%7.5f, y1=%7.5f)\n", x1, y1, x2, y1);
      if (print_flag) fprintf(stderr, "                    (x=%7.5f, y=%7.5f)\n", zz, ee);
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y2=%7.5f)                  (x2=%7.5f, y2=%7.5f)\n\n", x1, y2, x2, y2);
     
      if (print_flag) fprintf(stderr, "(z11=%7.3f)           (z21=%7.3f)\n", z11, z21);
      if (print_flag) fprintf(stderr, "             (z=%7.3f)\n", t);
      if (print_flag) fprintf(stderr, "(z12=%7.3f)           (z22=%7.3f)\n", z12, z22);
            
      return t;
}

/**
 * Calculates the optical depth tau based on the frequency nu, redshift zz, and range parameter.
 *
 * @param nu - the frequency in nu
 * @param zz - the redshift in zz
 * @param range - the range parameter
 * @return the calculated optical depth tau
 */
double tau_IRA_Kneiske(double nu, double zz, int range) {
      double E_keV, E_GeV, E_TeV, tau;
      
      E_keV = nu / keV;
      E_GeV = E_keV / 1.0e+6;
      E_TeV = E_keV / 1.0e+9;
      
      if (range == 0) 
        tau = Interpolate1(log10(E_GeV), zz, 0);
      else 
        tau = Interpolate2(log10(E_GeV), zz, 0);
      
      return tau;
}




/**************************************************************************
*  IR absorbtion 
*
*  franceschini 2008
*
*
****************************************************************************/

const int ZDIM_Franceschini = 9;
const int EDIM_Franceschini = 50;
const int MDIM_Franceschini = 450;

double Z3[ZDIM_Franceschini] = {0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0} ;

double E3[EDIM_Franceschini] =
{0.0200,0.0240,0.0289,0.0347,0.0417,0.0502,0.0603,0.0726,0.0873,0.104 ,0.126 ,0.151,
0.182,0.219 ,0.263 ,0.316 ,0.381 ,0.458 ,0.550 ,0.662 ,0.796 ,0.957 ,1.15  ,1.38 
,1.66  ,2.000 ,2.40  ,2.89  ,3.47  ,4.17  ,5.02  ,6.03  ,7.26  ,8.73  ,10.5  ,12.6 
,15.2  ,18.2  ,21.9  ,26.4  ,31.7  ,38.1  ,45.8  ,55.1  ,66.2  ,79.6  ,95.7  ,115. 
,138.  ,166.   } ; 

double T3[MDIM_Franceschini] = {0.0000, 0.0000, 0.0000, 0.000000, 0.0000021,
0.004933, 0.0399, 0.1157, 0.2596,
0.0000, 0.0000, 0.0000, 0.000000, 0.000188, 0.01284, 0.0718, 0.1783, 0.3635,
0.0000, 0.0000, 0.0000, 0.000000, 0.001304, 0.0279, 0.1188, 0.2598, 0.4919,
0.0000, 0.0000, 0.0000, 0.000488, 0.004558, 0.0533, 0.1833, 0.3635, 0.6517,
0.0000, 0.0000, 5.254E-05, 0.002276, 0.01157, 0.0921, 0.2689, 0.4967, 0.8548,
0.0000, 9.445E-05, 5.408E-04, 0.006575, 0.02436, 0.1480, 0.3836, 0.6745, 1.118,
1.0976E-04, 4.241E-04, 1.915E-03, 0.014592, 0.04512, 0.2275, 0.5434, 0.9179, 1.465,
3.0882E-04, 1.103E-03, 4.548E-03, 0.02771, 0.07684, 0.3430, 0.7707, 1.251, 1.917,
6.5619E-04, 2.258E-03, 8.903E-03, 0.04808, 0.1248, 0.5137, 1.092, 1.703, 2.503,
1.2130E-03, 4.097E-03, 1.582E-02, 0.07958, 0.1984, 0.7640, 1.537, 2.302, 3.249,
2.1063E-03, 7.039E-03, 2.685E-02, 0.1284, 0.3109, 1.120, 2.133, 3.073, 4.181,
3.5291E-03, 1.167E-02, 4.406E-02, 0.2031, 0.4780, 1.607, 2.905, 4.042, 5.318,
5.7051E-03, 1.872E-02, 7.010E-02, 0.3134, 0.7163, 2.247, 3.875, 5.225, 6.673,
8.9183E-03, 2.907E-02, 0.1082, 0.4696, 1.040, 3.056, 5.055, 6.627, 8.241,
1.3517E-02, 4.378E-02, 0.1618, 0.6809, 1.461, 4.042, 6.438, 8.226, 9.997,
1.9793E-02, 6.367E-02, 0.2338, 0.9517, 1.981, 5.192, 7.989, 9.977, 11.89,
2.7938E-02, 8.935E-02, 0.3256, 1.281, 2.594, 6.474, 9.650, 11.81, 13.89,
3.7957E-02, 0.1205, 0.4356, 1.661, 3.284, 7.836, 11.34, 13.67, 15.93,
4.9558E-02, 0.1563, 0.5607, 2.082, 4.023, 9.214, 13.01, 15.51, 18.08,
6.2291E-02, 0.1953, 0.6961, 2.524, 4.779, 10.55, 14.63, 17.39, 20.45,
7.5753E-02, 0.2364, 0.8373, 2.967, 5.517, 11.82, 16.25, 19.49, 23.27,
8.9194E-02, 0.2768, 0.9750, 3.389, 6.210, 13.03, 18.04, 22.02, 26.81,
0.1019, 0.3152, 1.105, 3.779, 6.846, 14.29, 20.21, 25.22, 31.33,
0.1136, 0.3501, 1.223, 4.129, 7.432, 15.73, 22.98, 29.37, 37.23,
0.1240, 0.3810, 1.327, 4.444, 8.010, 17.54, 26.58, 34.78, 45.09,
0.1329, 0.4076, 1.419, 4.747, 8.652, 19.87, 31.31, 41.95, 55.80,
0.1409, 0.4318, 1.504, 5.079, 9.452, 22.96, 37.67, 51.72, 70.71,
0.1486, 0.4560, 1.596, 5.498, 10.52, 27.08, 46.30, 65.17, 92.14,
0.1579, 0.4863, 1.714, 6.075, 11.96, 32.66, 58.24, 84.48, 124.0,
0.1710, 0.5284, 1.879, 6.875, 13.92, 40.39, 75.45, 113.2, 172.1,
0.1896, 0.5887, 2.113, 7.952, 16.57, 51.39, 101.3, 157.7, 245.9,
0.2162, 0.6732, 2.431, 9.421, 20.24, 67.70, 141.4, 227.3, 357.0,
0.2512, 0.7847, 2.858, 11.45, 25.57, 92.73, 204.9, 335.3, 519.3,
0.3017, 0.9447, 3.464, 14.41, 33.52, 132.2, 304.3, 496.4, 747.0,
0.3732, 1.171, 4.334, 18.87, 45.81, 195.0, 454.1, 724.1, 1048.,
0.4795, 1.513, 5.663, 25.76, 65.21, 292.5, 666.6, 1027., 1426.,
0.6455, 2.048, 7.723, 36.60, 95.98, 435.4, 948.5, 1407., 1873.,
0.8984, 2.871, 10.93, 53.79, 143.7, 630.5, 1299., 1852., 2372.,
1.297, 4.162, 15.99, 80.41, 214.3, 878.2, 1705., 2339., 2897.,
1.917, 6.177, 23.86, 119.9, 311.5, 1172., 2145., 2843., 3412.,
2.856, 9.181, 35.47, 174.8, 435.3, 1494., 2587., 3323., 3887.,
4.211, 13.47, 51.78, 245.0, 580.9, 1824., 3004., 3752., 4303.,
6.038, 19.14, 72.76, 327.4, 739.9, 2134., 3362., 4102., 4618.,
8.285, 26.00, 97.51, 416.9, 899.0, 2403., 3638., 4352., 4835.,
10.82, 33.59, 124.3, 506.3, 1045., 2616., 3825., 4499., 4940.,
13.48, 41.42, 151.2, 587.7, 1169., 2756., 3915., 4540., 4945.,
16.04, 48.81, 175.8, 655.4, 1263., 2823., 3915., 4540., 4945.,
18.24, 54.98, 195.8, 705.7, 1320., 2823., 3915., 4540., 4945.,
20.01, 59.82, 210.8, 735.5, 1340., 2823., 3915., 4540., 4945.,
21.20, 63.03, 219.8, 744.0, 1340., 2823., 3915., 4540., 4945.};



/**
 * Interpolate3 function calculates the interpolated value for a given input E and z.
 *
 * @param ee - The input value for E.
 * @param zz - The input value for z.
 * @param print_flag - A flag indicating whether to print debug information.
 * @return The interpolated value.
 *
 * The Interpolate3 function takes in two parameters ee (the input value for E) and zz (the input value for z), and calculates the interpolated value using the function:
 * t = a * ee * zz + b * zz + c * ee + d
 * The coefficients a, b, c, and d are calculated based on the surrounding data points and used to determine the interpolated value.
 * If the input values ee and zz fall within the defined range of E and z values, the function will return the interpolated value. Otherwise, it returns 0.
 *
 * The function uses the following data arrays:
 *
 * E3 - An array of size EDIM_Franceschini containing the E values.
 * EDIM_Franceschini - The size of the E3 array.
 *
 * Z3 - An array of size ZDIM_Franceschini containing the z values.
 * ZDIM_Franceschini - The size of the Z3 array.
 *
 * Note: The values of E3, EDIM_Franceschini, Z3, and ZDIM_Franceschini are not shown in the generated documentation.
 */
double Interpolate3(double ee, double zz, int print_flag) {
      int ei = 1, zi = 1;
      int col, row;
      
      double x1, y1, x2, y2;
      double z11, z12, z21, z22;
      double a, b, c, d, t;

      if ((ee >= E3[0]) && (ee <= E3[EDIM_Franceschini-1]) &&
          (zz >= Z3[0]) && (zz <= Z3[ZDIM_Franceschini-1])) {

        if (print_flag) fprintf(stderr, "Input: E=%7.5f z=%7.5f\n\n", ee, zz);
      
        for (col=0; col <= ZDIM_Franceschini-2; col++) {
           if ((zz >= Z3[col]) && (zz <= Z3[col+1])) {
                 zi = col;
             break;
           }
        }
            
        for (row=0; row <= EDIM_Franceschini-2; row++) {
           if ((ee >= E3[row]) && (ee <= E3[row+1])) {
             ei = row;
             break;
           }
        }
            
        if (print_flag) fprintf(stderr, "Found: E_min=%7.5f E_max=%7.5f ei=%d\n",
E3[ei], E3[ei+1], ei);
        if (print_flag) fprintf(stderr, "Found: z_min=%7.5f z_max=%7.5f zi=%d\n\n",
Z3[zi], Z3[zi+1], zi);
      
        x1  = Z3[zi];
        x2  = Z3[zi+1];
        y1  = E3[ei];
        y2  = E3[ei+1];
      
        z11 = T3[ZDIM_Franceschini*ei+zi];           // [ei][zi];
        z12 = T3[ZDIM_Franceschini*(ei+1)+zi];       // [ei+1][zi];
        z21 = T3[ZDIM_Franceschini*ei+zi+1];         // [ei][zi+1];
        z22 = T3[ZDIM_Franceschini*(ei+1)+zi+1];     // [ei+1][zi+1];
      
        a = (z11 - z12 - z21 + z22)             / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        b = (y1*z12 - y2*z11 - y1*z22 + y2*z21) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        c = (x1*z21 - x2*z11 - x1*z22 + x2*z12) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        d = (x1*y1*z22 - x1*y2*z21 - x2*y1*z12 + x2*y2*z11) / (x1*y1 - x1*y2 - x2*y1
+ x2*y2);
      
        // the function
        t = a*ee*zz + b*zz + c*ee + d;
      } else {
        t = 0;


        //fprintf(stderr, "Incorrect range:\n");
        //fprintf(stderr, "E_min=%e  E=%e  E_max=%e\n", E1[0], ee, E1[LDIM-1]);
        //fprintf(stderr, "Z_min=%e  Z=%e  Z_max=%e\n", Z1[0], zz, Z1[LDIM-1]);
      }
      
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y1=%7.5f)    (x2=%7.5f, y1=%7.5f)\n", x1, y1, x2, y1);
      if (print_flag) fprintf(stderr, "(x=%7.5f, y=%7.5f)\n",
zz, ee);
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y2=%7.5f)    (x2=%7.5f, y2=%7.5f)\n\n", x1, y2, x2, y2);
     
    if (print_flag) fprintf(stderr, "(z11=%7.3f, z11_dim=%d)   (z21=%7.3f,z12_dim=%d)\n", z11, ZDIM_Franceschini*ei+zi,z21,ZDIM_Franceschini*ei+zi+1 );
      if (print_flag) fprintf(stderr, "             (z=%7.3f)\n", t);
      if (print_flag) fprintf(stderr, "(z12=%7.3f, z12_dim=%d)    (z22=%7.3f,z12_dim=%d)\n", z12, ZDIM_Franceschini*(ei+1)+zi,  z22,ZDIM_Franceschini*(ei+1)+zi+1);
            
      return t;
}

/**
 * Calculates the optical depth using the IRA Franseschini method.
 *
 * @param nu - The frequency in hertz.
 * @param zz - The redshift.
 * @return The optical depth.
 *
 * This function calculates the optical depth using the IRA Franseschini method. It takes in two parameters:
 * - nu: The frequency in hertz.
 * - zz: The redshift.
 *
 * It performs the following steps:
 * 1. Convert the frequency from hertz to keV.
 * 2. Convert the frequency from keV to GeV and TeV.
 * 3. Interpolate the value of the optical depth using the Interpolate3 function.
 *
 * The function returns the calculated optical depth.
 *
 * Note: The function relies on the Interpolate3 function to perform the interpolation. The Interpolate3 function takes in two parameters ee (the input value for E) and zz (the input value for z), and calculates the interpolated value using the function:
 * t = a * ee * zz + b * zz + c * ee + d
 * The coefficients a, b, c, and d are calculated based on the surrounding data points and used to determine the interpolated value.
 * If the input values ee and zz fall within the defined range of E and z values, the function will return the interpolated value. Otherwise, it returns 0.
 *
 * The Interpolate3 function uses the following data arrays:
 *
 * - E3: An array of size EDIM_Franceschini containing the E values.
 * - EDIM_Franceschini: The size of the E3 array.
 *
 * - Z3: An array of size ZDIM_Franceschini containing the z values.
 * - ZDIM_Franceschini: The size of the Z3 array.
 *
 * Note: The values of E3, EDIM_Franceschini, Z3, and ZDIM_Franceschini are not shown in the generated documentation.
 */
double tau_IRA_Franceschini(double nu, double zz) {
      double E_keV, E_GeV, E_TeV, tau;
      
      E_keV = nu / keV;
      E_GeV = E_keV / 1.0e+6;
      E_TeV = E_keV / 1.0e+9;
      
        tau = Interpolate3(E_TeV, zz, 0);
      
      return tau;
}


/**************************************************************************
*  IR absorbtion 
*
*  Finke 2010
*
*
****************************************************************************/

const int ZDIM_Finke = 500;
const int EDIM_Finke = 99;
const int MDIM_Finke = 49500;

double Z4[ZDIM_Finke] = 
{0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65, 1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.0, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.5, 2.51, 2.52, 2.53, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 2.77, 2.78, 2.79, 2.8, 2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.0, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21, 3.22, 3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.3, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49, 3.5, 3.51, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.8, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4.0, 4.01, 4.02, 4.03, 4.04, 4.05, 4.06, 4.07, 4.08, 4.09, 4.1, 4.11, 4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.24, 4.25, 4.26, 4.27, 4.28, 4.29, 4.3, 4.31, 4.32, 4.33, 4.34, 4.35, 4.36, 4.37, 4.38, 4.39, 4.4, 4.41, 4.42, 4.43, 4.44, 4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51, 4.52, 4.53, 4.54, 4.55, 4.56, 4.57, 4.58, 4.59, 4.6, 4.61, 4.62, 4.63, 4.64, 4.65, 4.66, 4.67, 4.68, 4.69, 4.7, 4.71, 4.72, 4.73, 4.74, 4.75, 4.76, 4.77, 4.78, 4.79, 4.8, 4.81, 4.82, 4.83, 4.84, 4.85, 4.86, 4.87, 4.88, 4.89, 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99
} ;

double E4[EDIM_Finke] =
{
0.001, 0.001122018, 0.001258925, 0.001412538, 0.001584893, 0.001778279, 0.001995262, 0.002238721, 0.002511886, 0.002818383, 0.003162278, 0.003548134, 0.003981072, 0.004466836, 0.005011872, 0.005623413, 0.006309573, 0.007079458, 0.007943282, 0.008912509, 0.01, 0.01122018, 0.01258925, 0.01412538, 0.01584893, 0.01778279, 0.01995262, 0.02238721, 0.02511886, 0.02818383, 0.03162278, 0.03548134, 0.03981072, 0.04466836, 0.05011872, 0.05623413, 0.06309573, 0.07079458, 0.07943282, 0.08912509, 0.1, 0.1122018, 0.1258925, 0.1412538, 0.1584893, 0.1778279, 0.1995262, 0.2238721, 0.2511886, 0.2818383, 0.3162278, 0.3548134, 0.3981072, 0.4466836, 0.5011872, 0.5623413, 0.6309573, 0.7079458, 0.7943282, 0.8912509, 1.0, 1.122018, 1.258925, 1.412538, 1.584893, 1.778279, 1.995262, 2.238721, 2.511886, 2.818383, 3.162278, 3.548134, 3.981072, 4.466836, 5.011872, 5.623413, 6.309573, 7.079458, 7.943282, 8.912509, 10.0, 11.22018, 12.58925, 14.12538, 15.84893, 17.78279, 19.95262, 22.38721, 25.11886, 28.18383, 31.62278, 35.48134, 39.81072, 44.66836, 50.11872, 56.23413, 63.09573, 70.79458, 79.43282
} ;

double T4[MDIM_Finke] =
{1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1.546542e-60, 1.719893e-60, 4.880344e-60, 5.405836e-60, 5.9832e-60, 9.774622e-60, 1.076661e-59, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 2.210451e-60, 2.49603e-60, 2.823011e-60, 8.07931e-60, 9.097583e-60, 1.503252e-59, 1.68214e-59, 2.345338e-59, 2.606411e-59, 3.343169e-59, 3.689987e-59, 4.054556e-59, 4.938142e-59, 5.39379e-59, 6.575495e-59, 7.151923e-59, 8.64877e-59, 9.377525e-59, 1.059449e-58, 1.189789e-58, 1.285418e-58, 1.475882e-58, 1.585457e-58, 1.789735e-58, 1.91646e-58, 2.089287e-58, 2.279883e-58, 2.427891e-58, 2.694399e-58, 2.863069e-58, 3.121229e-58, 3.358977e-58, 3.554421e-58, 3.895887e-58, 4.11688e-58, 4.426574e-58, 4.716791e-58, 5.004019e-58, 5.358685e-58, 5.631839e-58, 5.998419e-58, 6.340696e-58, 6.703727e-58, 7.113759e-58, 7.439294e-58, 7.907026e-58, 8.307456e-58, 8.744875e-58, 9.169569e-58, 9.643412e-58, 1.012754e-57, 1.058933e-57, 1.109153e-57, 1.158013e-57, 1.210945e-57, 1.265684e-57, 1.312457e-57, 1.374883e-57, 1.429699e-57, 1.491474e-57, 1.548221e-57, 1.605615e-57, 1.67184e-57, 1.731658e-57, 1.796205e-57, 1.866772e-57, 1.92994e-57, 1.996194e-57, 2.065924e-57, 2.135129e-57, 2.212625e-57, 2.28568e-57, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 2.962195e-60, 3.411003e-60, 9.968835e-60, 8.318019e-59, 1.918063e-59, 1.679771e-58, 3.087217e-59, 1.832545e-58, 4.523831e-59, 2.009294e-58, 6.238183e-59, 2.211216e-58, 8.531264e-59, 2.471402e-58, 1.150201e-58, 2.800923e-58, 1.49831e-58, 3.178819e-58, 1.899451e-58, 3.607432e-58, 2.356083e-58, 4.088423e-58, 2.869636e-58, 4.622972e-58, 3.466548e-58, 5.240889e-58, 4.155339e-58, 5.951441e-58, 4.868743e-58, 6.73543e-58, 5.711701e-58, 7.594718e-58, 6.575037e-58, 8.530907e-58, 9.145375e-58, 8.149428e-58, 1.020351e-57, 9.310425e-58, 1.13701e-57, 1.058124e-57, 1.265403e-57, 1.188988e-57, 1.25726e-57, 1.500738e-57, 1.409603e-57, 1.642266e-57, 1.57121e-57, 1.799249e-57, 1.737504e-57, 1.822446e-57, 2.088203e-57, 2.011752e-57, 2.273937e-57, 2.209164e-57, 2.464208e-57, 2.594616e-57, 2.526184e-57, 2.811731e-57, 2.756855e-57, 3.033064e-57, 2.991322e-57, 3.121699e-57, 3.411897e-57, 3.372456e-57, 3.661586e-57, 3.808974e-57, 3.786078e-57, 4.083484e-57, 4.063245e-57, 4.208895e-57, 4.534532e-57, 4.517536e-57, 4.826472e-57, 4.823649e-57, 4.982623e-57, 5.322083e-57, 5.308848e-57, 5.491275e-57, 5.834783e-57, 5.823057e-57, 6.230071e-57, 6.232376e-57, 6.523921e-57, 6.875181e-57, 6.890607e-57, 7.173935e-57, 7.53481e-57, 7.661167e-57, 7.851102e-57, 8.323586e-57, 8.340363e-57, 8.551561e-57, 9.036201e-57, 9.046642e-57, 9.615403e-57, 9.835277e-57, 9.859482e-57, 1.052754e-56, 1.074638e-56, 1.105995e-56, 1.143993e-56, 1.185066e-56, 1.197882e-56, 1.236852e-56, 1.288805e-56, 1.289658e-56, 1.348753e-56, 1.381147e-56, 1.383518e-56, 1.443288e-56, 1.483519e-56, 1.521405e-56, 1.554137e-56, 1.593304e-56, 1.641253e-56, 1.663657e-56, 1.738414e-56, 1.751179e-56, 1.77235e-56, 1.847783e-56, 1.859975e-56, 1.898956e-56, 1.967843e-56, 1.989363e-56, 2.032775e-56, 2.064053e-56, 2.133546e-56, 2.161515e-56, 2.192384e-56, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 3.554395e-60, 1.061848e-59, 1.239937e-59, 2.127957e-59, 2.462939e-59, 3.557655e-59, 4.727276e-59, 5.382419e-59, 6.778773e-59, 7.633834e-59, 9.269381e-59, 1.13499e-58, 1.327677e-58, 1.541972e-58, 1.774042e-58, 2.028562e-58, 2.303519e-58, 2.664784e-58, 2.91822e-58, 3.328393e-58, 3.693544e-58, 4.081122e-58, 4.525884e-58, 5.029185e-58, 5.567279e-58, 6.067734e-58, 6.668353e-58, 7.230905e-58, 7.89611e-58, 8.539368e-58, 9.333896e-58, 9.965438e-58, 1.082943e-57, 1.156077e-57, 1.249771e-57, 1.333237e-57, 1.435653e-57, 1.527385e-57, 1.643536e-57, 1.737045e-57, 1.862039e-57, 1.964505e-57, 2.097585e-57, 2.207891e-57, 2.349404e-57, 2.473765e-57, 2.62985e-57, 2.758287e-57, 2.924005e-57, 3.062836e-57, 3.237951e-57, 3.386291e-57, 3.571459e-57, 3.727848e-57, 6.670379e-57, 4.091896e-57, 9.844683e-57, 4.476686e-57, 1.028552e-56, 4.8856e-57, 1.07466e-56, 5.313456e-57, 1.121722e-56, 5.756407e-57, 1.171199e-56, 6.227593e-57, 1.222251e-56, 6.721441e-57, 1.274242e-56, 7.236702e-57, 1.328404e-56, 7.833209e-57, 1.391007e-56, 8.522822e-57, 1.462428e-56, 9.224536e-57, 1.533407e-56, 9.935592e-57, 1.605956e-56, 1.06695e-56, 1.679511e-56, 1.141689e-56, 1.756223e-56, 1.219742e-56, 1.846904e-56, 1.913855e-56, 1.369987e-56, 2.004832e-56, 1.475803e-56, 2.109382e-56, 1.583201e-56, 2.215706e-56, 1.692674e-56, 1.726707e-56, 2.417383e-56, 1.837993e-56, 2.526929e-56, 1.962746e-56, 2.60867e-56, 2.099706e-56, 2.136778e-56, 2.861439e-56, 2.273733e-56, 2.997819e-56, 2.412224e-56, 3.131316e-56, 3.188231e-56, 2.641246e-56, 3.324337e-56, 2.732801e-56, 3.470245e-56, 2.884899e-56, 2.999638e-56, 3.688919e-56, 3.165718e-56, 3.846022e-56, 3.931714e-56, 3.385548e-56, 4.063654e-56, 3.537465e-56, 3.589235e-56, 4.329696e-56, 3.757978e-56, 4.489208e-56, 3.920391e-56, 3.975893e-56, 4.738605e-56, 4.146938e-56, 4.278516e-56, 5.006605e-56, 4.438751e-56, 5.157894e-56, 4.570244e-56, 4.692013e-56, 5.423685e-56, 2.016715e-11, 6.626979e-11, 1.462541e-10, 3.032946e-10, 5.00459e-10, 7.876226e-10, 1.160267e-09, 1.596836e-09, 2.128997e-09, 2.750923e-09, 3.48638e-09, 4.28739e-09, 5.173187e-09, 6.25702e-09, 8.697354e-09, 1.011945e-08, 1.443103e-08, 1.645694e-08, 1.874978e-08, 2.419341e-08, 2.714962e-08, 3.32182e-08, 3.722207e-08, 4.119824e-08, 4.853915e-08, 5.369908e-08, 6.158945e-08, 6.753739e-08, 7.375358e-08, 8.58197e-08, 9.347997e-08, 1.094386e-07, 1.183699e-07, 1.335793e-07, 1.465882e-07, 1.581501e-07, 1.782651e-07, 1.910767e-07, 2.106707e-07, 2.277029e-07, 2.436118e-07, 2.677954e-07, 2.850277e-07, 3.083188e-07, 3.361099e-07, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 3.649026e-60, 4.368227e-60, 1.340265e-59, 2.352651e-59, 3.488634e-59, 4.083731e-59, 5.526313e-59, 7.105325e-59, 8.15091e-59, 1.007488e-58, 1.253176e-58, 1.485607e-58, 1.75126e-58, 2.11507e-58, 2.435697e-58, 2.792384e-58, 3.24979e-58, 3.667554e-58, 4.121063e-58, 4.678679e-58, 5.236661e-58, 5.871284e-58, 6.479485e-58, 7.279517e-58, 8.05204e-58, 8.831126e-58, 9.775695e-58, 1.069601e-57, 1.159725e-57, 1.261717e-57, 1.372344e-57, 1.492712e-57, 1.604572e-57, 1.737471e-57, 1.875487e-57, 2.014004e-57, 2.164899e-57, 2.308794e-57, 2.480093e-57, 2.635102e-57, 2.816552e-57, 2.992435e-57, 3.191498e-57, 3.389986e-57, 3.588791e-57, 3.808308e-57, 4.021385e-57, 4.257946e-57, 4.47847e-57, 4.728763e-57, 4.966884e-57, 5.236097e-57, 5.500384e-57, 5.785555e-57, 6.066469e-57, 6.358518e-57, 6.666157e-57, 6.973437e-57, 7.294029e-57, 7.609988e-57, 7.95017e-57, 8.285838e-57, 8.646769e-57, 9.001904e-57, 9.457614e-57, 9.980817e-57, 1.037688e-56, 1.091108e-56, 1.132585e-56, 1.187989e-56, 1.229551e-56, 1.286996e-56, 1.346038e-56, 1.390051e-56, 1.450071e-56, 1.51264e-56, 1.574235e-56, 1.651017e-56, 1.714338e-56, 1.792278e-56, 1.886133e-56, 1.9358e-56, 2.031746e-56, 2.08158e-56, 2.179002e-56, 2.247331e-56, 2.344408e-56, 2.413714e-56, 2.526795e-56, 2.596232e-56, 2.710067e-56, 2.783691e-56, 2.896145e-56, 3.00157e-56, 3.084347e-56, 3.191636e-56, 3.275342e-56, 3.382047e-56, 3.483612e-56, 3.571946e-56, 5.420882e-56, 3.7971e-56, 7.377786e-56, 3.991511e-56, 7.625088e-56, 4.218573e-56, 7.875833e-56, 4.446865e-56, 8.12625e-56, 4.678632e-56, 8.3902e-56, 4.925931e-56, 8.671033e-56, 5.188443e-56, 8.934792e-56, 5.450411e-56, 9.217346e-56, 5.716188e-56, 9.464805e-56, 5.984727e-56, 2.443722e-12, 6.039346e-11, 1.586089e-10, 3.900598e-10, 7.130955e-10, 1.226157e-09, 1.84274e-09, 2.693082e-09, 3.650971e-09, 4.851941e-09, 6.17749e-09, 7.744465e-09, 9.526162e-09, 1.362017e-08, 1.623663e-08, 2.326866e-08, 2.710394e-08, 3.540292e-08, 4.060257e-08, 5.052061e-08, 5.723238e-08, 6.892088e-08, 7.705781e-08, 8.628245e-08, 1.003332e-07, 1.115429e-07, 1.312176e-07, 1.445972e-07, 1.711629e-07, 1.881928e-07, 2.174323e-07, 2.372558e-07, 2.622777e-07, 2.935801e-07, 3.187204e-07, 3.570496e-07, 3.860086e-07, 4.271249e-07, 4.602565e-07, 5.126254e-07, 5.486458e-07, 5.939728e-07, 6.522751e-07, 6.957475e-07, 7.673591e-07, 8.155733e-07, 8.914434e-07, 9.476424e-07, 1.025752e-06, 1.090603e-06, 1.153992e-06, 1.244801e-06, 1.313822e-06, 1.413738e-06, 1.496866e-06, 1.602818e-06, 1.699373e-06, 1.784615e-06, 1.910087e-06, 2.009827e-06, 2.144767e-06, 2.249712e-06, 2.382412e-06, 2.502144e-06, 2.622601e-06, 2.775537e-06, 2.89916e-06, 3.063972e-06, 3.202365e-06, 3.3709e-06, 3.523358e-06, 3.674429e-06, 3.871523e-06, 4.033768e-06, 4.225617e-06, 4.402848e-06, 4.582697e-06, 4.805135e-06, 4.97844e-06, 5.208347e-06, 5.409998e-06, 5.609165e-06, 5.847442e-06, 6.068949e-06, 6.331244e-06, 6.553121e-06, 6.796592e-06, 7.063419e-06, 7.304289e-06, 7.588457e-06, 7.837294e-06, 8.121238e-06, 8.413661e-06, 8.664194e-06, 8.985365e-06, 9.266905e-06, 9.580708e-06, 9.882209e-06, 1.018452e-05, 1.058094e-05, 1.089015e-05, 1.122441e-05, 1.168442e-05, 1.199991e-05, 1.252791e-05, 1.287674e-05, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 3.596806e-60, 4.400357e-60, 1.385575e-59, 2.501579e-59, 3.797456e-59, 5.284717e-59, 6.97362e-59, 8.168584e-59, 1.031201e-58, 1.307907e-58, 1.651171e-58, 1.95861e-58, 2.314027e-58, 2.783084e-58, 3.294362e-58, 3.849087e-58, 4.31246e-58, 4.966296e-58, 5.706076e-58, 6.461165e-58, 7.359486e-58, 8.179133e-58, 9.213801e-58, 1.023399e-57, 1.132654e-57, 1.249592e-57, 1.37365e-57, 1.515964e-57, 1.662171e-57, 1.809809e-57, 1.975137e-57, 2.148034e-57, 2.328905e-57, 2.519795e-57, 2.712378e-57, 2.922156e-57, 3.145244e-57, 3.373119e-57, 3.61393e-57, 3.871499e-57, 4.131063e-57, 4.405468e-57, 4.694706e-57, 4.986963e-57, 5.300066e-57, 5.620318e-57, 5.944946e-57, 6.299193e-57, 6.641248e-57, 7.015353e-57, 7.403757e-57, 7.791703e-57, 8.198471e-57, 8.60831e-57, 9.051201e-57, 9.472655e-57, 9.947684e-57, 1.049913e-56, 1.098108e-56, 1.165155e-56, 1.231718e-56, 1.284883e-56, 1.35288e-56, 1.425246e-56, 1.480475e-56, 1.552151e-56, 1.632299e-56, 1.706954e-56, 1.784917e-56, 1.895561e-56, 1.992223e-56, 2.071978e-56, 2.169495e-56, 2.28732e-56, 2.369361e-56, 2.471687e-56, 2.590061e-56, 2.694664e-56, 2.799552e-56, 2.921066e-56, 3.02656e-56, 3.185036e-56, 3.291161e-56, 3.418493e-56, 3.545871e-56, 3.688401e-56, 3.819593e-56, 3.929655e-56, 4.080111e-56, 4.191636e-56, 4.391222e-56, 4.524846e-56, 4.674964e-56, 4.844406e-56, 4.961194e-56, 5.132296e-56, 5.28153e-56, 5.459778e-56, 5.577143e-56, 5.768305e-56, 5.91876e-56, 6.130564e-56, 6.252377e-56, 6.46645e-56, 6.585598e-56, 6.801655e-56, 6.959878e-56, 7.176246e-56, 6.774458e-11, 2.172562e-10, 6.047296e-10, 1.193581e-09, 2.052321e-09, 3.203894e-09, 4.641204e-09, 6.443149e-09, 8.573706e-09, 1.108876e-08, 1.389943e-08, 2.017621e-08, 2.434105e-08, 3.573509e-08, 4.764291e-08, 5.650728e-08, 7.192754e-08, 8.278457e-08, 1.012021e-07, 1.149544e-07, 1.365338e-07, 1.531217e-07, 1.779158e-07, 2.044301e-07, 2.473701e-07, 2.732828e-07, 3.228933e-07, 3.548948e-07, 4.109874e-07, 4.507148e-07, 5.121515e-07, 5.589016e-07, 6.267438e-07, 6.910249e-07, 7.657819e-07, 8.472809e-07, 9.305225e-07, 1.023548e-06, 1.117358e-06, 1.221813e-06, 1.324754e-06, 1.442014e-06, 1.558151e-06, 1.68607e-06, 1.815174e-06, 1.954078e-06, 2.095999e-06, 2.256761e-06, 2.413168e-06, 2.599812e-06, 2.771403e-06, 2.972345e-06, 3.155516e-06, 3.37953e-06, 3.576464e-06, 3.815788e-06, 4.016455e-06, 4.282616e-06, 4.483536e-06, 4.794174e-06, 5.015489e-06, 5.347437e-06, 5.597603e-06, 5.953791e-06, 6.242953e-06, 6.568254e-06, 6.908177e-06, 7.260401e-06, 7.609296e-06, 7.964885e-06, 8.358838e-06, 8.684016e-06, 9.123493e-06, 9.524778e-06, 9.984623e-06, 1.040398e-05, 1.086587e-05, 1.134271e-05, 1.181421e-05, 1.228328e-05, 1.280586e-05, 1.329868e-05, 1.381259e-05, 1.433913e-05, 1.488334e-05, 1.549303e-05, 1.600953e-05, 1.665946e-05, 1.728458e-05, 1.794195e-05, 1.870962e-05, 1.932717e-05, 2.021907e-05, 2.081645e-05, 2.17147e-05, 2.2393e-05, 2.307671e-05, 2.397493e-05, 2.463841e-05, 2.556397e-05, 2.632049e-05, 2.74483e-05, 2.8193e-05, 2.918935e-05, 3.029573e-05, 3.113296e-05, 3.242139e-05, 3.327368e-05, 3.471555e-05, 3.55022e-05, 3.673803e-05, 3.771391e-05, 3.860027e-05, 4.008799e-05, 4.089348e-05, 4.262158e-05, 4.353248e-05, 4.48693e-05, 4.647084e-05, 4.734846e-05, 4.936982e-05, 5.031927e-05, 5.217264e-05, 5.326093e-05, 5.426276e-05, 5.628495e-05, 5.723598e-05, 5.916655e-05, 6.02735e-05, 6.217498e-05, 6.359491e-05, 6.469771e-05, 6.729444e-05, 6.849876e-05, 7.063353e-05, 7.200757e-05, 7.331349e-05, 7.543438e-05, 7.679045e-05, 7.922313e-05, 8.06714e-05, 8.185882e-05, 8.409653e-05, 8.548702e-05, 8.771248e-05, 8.972336e-05, 9.076994e-05, 9.38425e-05, 9.526489e-05, 9.768518e-05, 9.940866e-05, 0.0001009061, 0.0001035304, 0.0001050317, 0.0001077691, 0.0001095213, 0.0001104431, 0.0001137014, 0.00011544, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 3.835457e-60, 1.252463e-59, 2.33493e-59, 3.656254e-59, 5.235866e-59, 7.089807e-59, 9.231703e-59, 1.167331e-58, 1.484185e-58, 1.810408e-58, 2.265264e-58, 2.776548e-58, 3.348769e-58, 3.980233e-58, 4.675433e-58, 5.435387e-58, 6.301665e-58, 7.215081e-58, 8.21248e-58, 9.382736e-58, 1.064067e-57, 1.199263e-57, 1.344748e-57, 1.491031e-57, 1.644819e-57, 1.829216e-57, 2.024917e-57, 2.234364e-57, 2.446087e-57, 2.678805e-57, 2.90234e-57, 3.160394e-57, 3.428482e-57, 3.721051e-57, 4.020947e-57, 4.334611e-57, 4.653608e-57, 5.006607e-57, 5.3671e-57, 5.741561e-57, 6.137734e-57, 6.546953e-57, 6.969646e-57, 7.418573e-57, 7.884624e-57, 8.361661e-57, 8.854694e-57, 9.364659e-57, 9.900605e-57, 1.044918e-56, 1.102002e-56, 1.168492e-56, 1.248755e-56, 1.330105e-56, 1.410085e-56, 1.476547e-56, 1.561142e-56, 1.647258e-56, 1.738604e-56, 1.808766e-56, 1.920452e-56, 2.051364e-56, 2.181725e-56, 2.262276e-56, 2.395236e-56, 2.531936e-56, 2.652419e-56, 2.753944e-56, 2.897501e-56, 3.056783e-56, 3.217988e-56, 3.365137e-56, 3.51065e-56, 3.676168e-56, 3.827504e-56, 3.978281e-56, 4.14884e-56, 4.303e-56, 4.514993e-56, 4.654395e-56, 4.847264e-56, 5.061341e-56, 5.257789e-56, 5.421648e-56, 5.639815e-56, 5.84341e-56, 6.009845e-56, 6.232154e-56, 6.457296e-56, 6.64085e-56, 6.853368e-56, 7.09562e-56, 7.325973e-56, 7.540943e-56, 7.766079e-56, 8.022818e-56, 8.364214e-11, 3.192983e-10, 8.768276e-10, 1.798142e-09, 3.111419e-09, 4.915831e-09, 7.256416e-09, 1.006712e-08, 1.366724e-08, 1.802937e-08, 2.672937e-08, 4.064984e-08, 4.959955e-08, 6.737153e-08, 8.758631e-08, 1.033816e-07, 1.283436e-07, 1.487775e-07, 1.791496e-07, 2.121683e-07, 2.488239e-07, 2.89231e-07, 3.502895e-07, 4.101209e-07, 4.690628e-07, 5.390492e-07, 6.082437e-07, 6.927113e-07, 7.877996e-07, 8.775194e-07, 9.775707e-07, 1.092705e-06, 1.233414e-06, 1.366539e-06, 1.515141e-06, 1.668136e-06, 1.841465e-06, 2.03302e-06, 2.223815e-06, 2.422361e-06, 2.632535e-06, 2.853551e-06, 3.100292e-06, 3.383029e-06, 3.659474e-06, 3.946252e-06, 4.26185e-06, 4.579047e-06, 4.955385e-06, 5.286098e-06, 5.675134e-06, 6.045398e-06, 6.467413e-06, 6.906199e-06, 7.367843e-06, 7.829988e-06, 8.351202e-06, 8.882321e-06, 9.391604e-06, 9.956355e-06, 1.052328e-05, 1.112131e-05, 1.171917e-05, 1.237868e-05, 1.304087e-05, 1.376581e-05, 1.445688e-05, 1.522398e-05, 1.597692e-05, 1.67815e-05, 1.757921e-05, 1.840672e-05, 1.930693e-05, 2.021417e-05, 2.109243e-05, 2.205385e-05, 2.300361e-05, 2.402382e-05, 2.515793e-05, 2.62137e-05, 2.753644e-05, 2.923741e-05, 3.002314e-05, 3.240343e-05, 3.26173e-05, 3.505051e-05, 3.531637e-05, 3.777991e-05, 3.809612e-05, 4.063526e-05, 4.128134e-05, 4.386507e-05, 4.486657e-05, 4.747786e-05, 4.855237e-05, 5.11373e-05, 5.225797e-05, 5.488327e-05, 5.609724e-05, 5.876747e-05, 6.00588e-05, 6.2794e-05, 6.410581e-05, 6.724679e-05, 6.859044e-05, 7.212038e-05, 7.348795e-05, 7.705718e-05, 7.761295e-05, 8.209365e-05, 8.277471e-05, 8.729037e-05, 8.916056e-05, 9.137762e-05, 9.440846e-05, 9.723202e-05, 0.0001002787, 0.0001023797, 0.0001067497, 0.0001087786, 0.0001119446, 0.000115226, 0.0001184135, 0.0001217106, 0.0001247507, 0.0001282729, 0.0001314078, 0.0001336777, 0.0001393666, 0.00014099, 0.0001462589, 0.0001485382, 0.0001523022, 0.0001578551, 0.0001588413, 0.0001652657, 0.000166444, 0.0001732313, 0.0001744733, 0.0001791566, 0.0001835162, 0.0001849528, 0.0001920604, 0.0001949852, 0.0001994505, 0.0002036204, 0.0002084695, 0.0002116662, 0.0002159772, 0.000220034, 0.0002246457, 0.0002294301, 0.0002327188, 0.0002393258, 0.0002409876, 0.0002460202, 0.0002528283, 0.00025447, 0.0002619192, 0.0002642408, 0.0002697498, 0.0002766234, 0.0002792982, 0.0002841172, 0.0002894474, 0.0002947039, 0.0002982829, 0.0003048862, 0.0003075824, 0.0003113452, 0.0003186641, 0.0003212128, 0.0003298175, 0.0003347364, 0.0003389716, 0.0003447578, 0.0003487399, 0.0003543982, 0.000360134, 0.0003654548, 0.0003685693, 0.0003741534, 0.0003797215, 0.0003824815, 0.0003901171, 0.0003964997, 0.000397061, 0.0004049355, 0.0004118326, 0.0004158876, 0.0004200862, 0.0004278808, 0.000432013, 0.0004357385, 0.000443864, 0.0004467965, 0.0004509287, 0.0004586818, 0.0004614277, 0.0004674838, 0.0004751916, 0.0004791196, 0.0004827661, 0.0004887845, 0.0004947788, 0.0004994595, 0.0005058822, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 4.087752e-60, 1.382444e-59, 3.372333e-59, 4.991541e-59, 6.963129e-59, 9.309075e-59, 1.204844e-58, 1.563004e-58, 2.094416e-58, 2.627485e-58, 3.31383e-58, 4.003563e-58, 4.776079e-58, 5.63312e-58, 6.619506e-58, 7.829946e-58, 9.082365e-58, 1.045564e-57, 1.19591e-57, 1.358331e-57, 1.53448e-57, 1.727166e-57, 1.931633e-57, 2.15988e-57, 2.403608e-57, 2.665596e-57, 2.943563e-57, 3.239286e-57, 3.555956e-57, 3.895242e-57, 4.254529e-57, 4.624927e-57, 5.006925e-57, 5.425226e-57, 5.863834e-57, 6.322808e-57, 6.80652e-57, 7.308191e-57, 7.844833e-57, 8.397947e-57, 8.966782e-57, 9.566258e-57, 1.018994e-56, 1.084805e-56, 1.153044e-56, 1.23288e-56, 1.322294e-56, 1.415707e-56, 1.509406e-56, 1.607752e-56, 1.707981e-56, 1.81271e-56, 1.899295e-56, 2.025873e-56, 2.172542e-56, 2.320012e-56, 2.471982e-56, 2.626525e-56, 2.78208e-56, 2.885454e-56, 3.045174e-56, 3.225987e-56, 3.428931e-56, 3.633472e-56, 3.820382e-56, 3.954283e-56, 4.165926e-56, 4.374742e-56, 4.595636e-56, 4.814077e-56, 4.994488e-56, 5.216979e-56, 5.479423e-56, 5.722699e-56, 5.948875e-56, 6.18087e-56, 6.39252e-56, 6.668776e-56, 6.923084e-56, 7.237268e-56, 7.496923e-56, 7.739597e-56, 8.022267e-56, 8.283766e-56, 8.619269e-56, 5.332281e-11, 3.01421e-10, 9.625167e-10, 2.147664e-09, 3.904237e-09, 6.31308e-09, 9.419296e-09, 1.334809e-08, 1.806823e-08, 2.761956e-08, 4.26286e-08, 5.335101e-08, 7.433598e-08, 9.850131e-08, 1.270807e-07, 1.509424e-07, 1.880581e-07, 2.292225e-07, 2.831861e-07, 3.258335e-07, 4.023019e-07, 4.884037e-07, 5.723514e-07, 6.585077e-07, 7.705025e-07, 8.833632e-07, 9.942115e-07, 1.133721e-06, 1.292505e-06, 1.481863e-06, 1.672679e-06, 1.864754e-06, 2.103029e-06, 2.346289e-06, 2.587316e-06, 2.882239e-06, 3.180725e-06, 3.493243e-06, 3.832085e-06, 4.240332e-06, 4.64321e-06, 5.056119e-06, 5.536589e-06, 6.021123e-06, 6.500773e-06, 7.039367e-06, 7.614455e-06, 8.210153e-06, 8.829936e-06, 9.566799e-06, 1.028785e-05, 1.099643e-05, 1.180949e-05, 1.263788e-05, 1.349469e-05, 1.441211e-05, 1.526631e-05, 1.633747e-05, 1.730277e-05, 1.839725e-05, 1.947275e-05, 2.068821e-05, 2.189811e-05, 2.312162e-05, 2.443801e-05, 2.578281e-05, 2.719104e-05, 2.85786e-05, 3.0112e-05, 3.154215e-05, 3.345512e-05, 3.505609e-05, 3.710437e-05, 3.92077e-05, 4.100828e-05, 4.315697e-05, 4.501651e-05, 4.7359e-05, 4.959671e-05, 5.168762e-05, 5.411445e-05, 5.674396e-05, 5.996133e-05, 6.256787e-05, 6.565084e-05, 6.836005e-05, 7.152752e-05, 7.426224e-05, 7.754868e-05, 8.1146e-05, 8.385049e-05, 8.757127e-05, 9.122274e-05, 9.466557e-05, 9.894647e-05, 0.0001035714, 0.0001076207, 0.0001115348, 0.0001156973, 0.0001198956, 0.0001241514, 0.0001282387, 0.000133729, 0.0001375517, 0.0001431453, 0.0001476216, 0.0001534228, 0.000158, 0.0001637891, 0.0001684681, 0.0001756402, 0.0001791648, 0.0001863749, 0.0001899825, 0.0001972418, 0.0002020108, 0.0002091822, 0.0002145736, 0.0002235105, 0.0002274629, 0.0002364787, 0.0002403238, 0.0002492192, 0.0002531846, 0.0002674262, 0.0002665223, 0.0002857815, 0.000280912, 0.000300126, 0.0002962697, 0.0003148906, 0.0003111168, 0.0003321806, 0.0003266352, 0.0003478357, 0.0003423222, 0.0003623052, 0.0003581911, 0.0003777425, 0.0003740667, 0.0003933089, 0.0003920002, 0.0004093739, 0.0004099897, 0.0004256868, 0.0004266738, 0.0004431091, 0.0004442291, 0.0004605112, 0.0004616797, 0.000477315, 0.0004782125, 0.0004954728, 0.000497051, 0.0005160243, 0.0005259459, 0.0005250269, 0.0005431123, 0.000543463, 0.0005615758, 0.0005638887, 0.0005810978, 0.0005813205, 0.0005885757, 0.0006140198, 0.0006101997, 0.0006332849, 0.0006314628, 0.0006524221, 0.0006512857, 0.0006598755, 0.0006848375, 0.0006807927, 0.0007041908, 0.0007001773, 0.0007249315, 0.0007345886, 0.0007330003, 0.0007546118, 0.0007508704, 0.000776019, 0.0007734316, 0.0007877658, 0.0008052464, 0.0008107121, 0.0008296918, 0.0008424471, 0.0008412621, 0.0008585376, 0.0008615396, 0.0008715374, 0.0008946678, 0.0008928596, 0.0009188998, 0.0009150618, 0.0009251044, 0.0009494915, 0.0009449759, 0.0009646326, 0.0009843579, 0.0009818515, 0.001005693, 0.001004092, 0.001013816, 0.00103699, 0.001036258, 0.001044246, 0.001073736, 0.001072434, 0.001078391, 0.001107135, 0.0011033, 0.001113431, 0.001140606, 0.001136591, 0.001156974, 0.001169696, 0.00116798, 0.001188247, 0.001205156, 0.001203824, 0.001222164, 0.001239031, 0.001235816, 0.001257492, 0.001270938, 0.00126898, 0.001290111, 0.001302334, 0.001300039, 0.001308545, 0.001335793, 0.001333308, 0.001339932, 0.001365665, 0.001366402, 0.001374206, 0.001397326, 0.00139663, 0.001405604, 0.001425751, 0.001428369, 0.001438973, 0.001455323, 0.00147131, 0.001466215, 0.001474239, 0.001499162, 0.001499444, 0.001507427, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 3.70871e-60, 1.972477e-59, 3.323052e-59, 5.074794e-59, 7.987556e-59, 1.071673e-58, 1.431299e-58, 1.971326e-58, 2.604092e-58, 3.264764e-58, 4.098148e-58, 5.039429e-58, 6.056741e-58, 7.235535e-58, 8.721581e-58, 1.02783e-57, 1.192758e-57, 1.374015e-57, 1.572146e-57, 1.812707e-57, 2.056286e-57, 2.329066e-57, 2.61525e-57, 2.923253e-57, 3.267289e-57, 3.625014e-57, 4.010816e-57, 4.421095e-57, 4.866167e-57, 5.33717e-57, 5.841749e-57, 6.359237e-57, 6.910828e-57, 7.493916e-57, 8.111496e-57, 8.762241e-57, 9.426413e-57, 1.014271e-56, 1.088101e-56, 1.166547e-56, 1.256623e-56, 1.359852e-56, 1.465598e-56, 1.574085e-56, 1.686061e-56, 1.80124e-56, 1.919357e-56, 2.060788e-56, 2.221881e-56, 2.385293e-56, 2.553242e-56, 2.725204e-56, 2.899654e-56, 3.077103e-56, 3.259467e-56, 3.498585e-56, 3.722101e-56, 3.949316e-56, 4.180015e-56, 4.413226e-56, 4.650548e-56, 4.891913e-56, 5.114688e-56, 5.395703e-56, 5.67783e-56, 5.906496e-56, 6.199537e-56, 6.495031e-56, 6.795254e-56, 7.101573e-56, 7.425675e-56, 7.753903e-56, 8.103732e-56, 8.338968e-56, 8.69551e-56, 9.056468e-56, 6.535496e-11, 3.995353e-10, 1.232034e-09, 2.697002e-09, 4.89281e-09, 7.900173e-09, 1.180793e-08, 1.671045e-08, 2.705486e-08, 4.398096e-08, 6.464572e-08, 8.947565e-08, 1.174725e-07, 1.451825e-07, 1.843322e-07, 2.283011e-07, 2.863484e-07, 3.599866e-07, 4.417439e-07, 5.130846e-07, 6.198698e-07, 7.406013e-07, 8.702149e-07, 1.014444e-06, 1.166524e-06, 1.33944e-06, 1.561809e-06, 1.80385e-06, 2.055901e-06, 2.324584e-06, 2.608558e-06, 2.951648e-06, 3.306155e-06, 3.710353e-06, 4.157495e-06, 4.58828e-06, 5.116364e-06, 5.658967e-06, 6.271244e-06, 6.845176e-06, 7.53038e-06, 8.23656e-06, 8.979708e-06, 9.848045e-06, 1.064915e-05, 1.159594e-05, 1.25736e-05, 1.359291e-05, 1.465008e-05, 1.579402e-05, 1.702239e-05, 1.831712e-05, 1.964794e-05, 2.105689e-05, 2.257551e-05, 2.413563e-05, 2.574986e-05, 2.750211e-05, 2.929963e-05, 3.117498e-05, 3.311325e-05, 3.512228e-05, 3.752138e-05, 3.972104e-05, 4.234869e-05, 4.515108e-05, 4.758204e-05, 5.055735e-05, 5.348113e-05, 5.659595e-05, 5.94776e-05, 6.276311e-05, 6.663457e-05, 7.066221e-05, 7.42751e-05, 7.88291e-05, 8.310104e-05, 8.708321e-05, 9.170988e-05, 9.634839e-05, 0.0001005854, 0.00010534, 0.0001112974, 0.0001170635, 0.0001221911, 0.0001288086, 0.0001349348, 0.0001407281, 0.0001469057, 0.0001538218, 0.000159859, 0.000166517, 0.0001731068, 0.0001827634, 0.0001890294, 0.0001975396, 0.000204234, 0.0002138464, 0.0002216057, 0.0002293627, 0.0002383279, 0.0002469825, 0.0002567472, 0.0002645791, 0.0002754636, 0.0002856148, 0.0002962656, 0.0003057995, 0.0003154472, 0.0003267671, 0.0003371551, 0.0003484375, 0.0003571704, 0.000369862, 0.0003787719, 0.0003951118, 0.0004044025, 0.0004181507, 0.0004291588, 0.000441775, 0.0004547223, 0.0004676475, 0.000480878, 0.0004919855, 0.0005068312, 0.0005203729, 0.0005338796, 0.0005478525, 0.0005643311, 0.0005756254, 0.000591768, 0.0006029766, 0.0006223527, 0.0006311688, 0.000650609, 0.0006609607, 0.0006800149, 0.0006955967, 0.0007153609, 0.0007294416, 0.0007467096, 0.0007614151, 0.0007769169, 0.0007951032, 0.0008091845, 0.0008250855, 0.0008432896, 0.0008585833, 0.0008786648, 0.0008942851, 0.0009131897, 0.0009331019, 0.0009468795, 0.0009683226, 0.0009823612, 0.001004003, 0.001018522, 0.001036397, 0.001056621, 0.001073471, 0.001094763, 0.001113522, 0.001152055, 0.001150569, 0.00121259, 0.00118901, 0.001252226, 0.00122699, 0.001287778, 0.001257809, 0.001328216, 0.001299482, 0.001369995, 0.001340475, 0.001413517, 0.001380684, 0.001451013, 0.001420962, 0.001492222, 0.001464387, 0.001532819, 0.001505091, 0.001570735, 0.001542052, 0.001617737, 0.0015888, 0.00165537, 0.001631827, 0.001693892, 0.001671951, 0.001732603, 0.001711213, 0.001779006, 0.001804168, 0.001770833, 0.001845101, 0.00181389, 0.001890577, 0.001862976, 0.001930498, 0.001904192, 0.001923111, 0.001990613, 0.001968304, 0.002040134, 0.002005744, 0.002082061, 0.002048712, 0.002076932, 0.002146378, 0.0021192, 0.002187913, 0.002160287, 0.002227805, 0.002249939, 0.002228074, 0.002286756, 0.002271874, 0.002335809, 0.002304652, 0.002333406, 0.00239771, 0.002375334, 0.002443709, 0.00246945, 0.002441948, 0.002512446, 0.002479397, 0.002506291, 0.002580271, 0.002538638, 0.002620566, 0.002587469, 0.002605214, 0.002683352, 0.002652616, 0.002667405, 0.002740934, 0.002716156, 0.002776882, 0.002759944, 0.002776146, 0.002846317, 0.002815775, 0.002833698, 0.002906192, 0.002879134, 0.002898241, 0.00296189, 0.002945038, 0.002958724, 0.003023992, 0.003003759, 0.003069245, 0.003088505, 0.003060061, 0.003128091, 0.003147746, 0.003119436, 0.003190759, 0.003201474, 0.003180754, 0.003248657, 0.003258338, 0.003239287, 0.00330817, 0.003324087, 0.003290464, 0.003310132, 0.003380684, 0.003348262, 0.003365889, 0.003434327, 0.003402125, 0.003419216, 0.003487997, 0.003454613, 0.00347862, 0.003541703, 0.00350672, 0.003525212, 0.003599065, 0.003612405, 0.003577489, 0.003595013, 0.003654539, 0.003645693, 0.003645352, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 3.070657e-60, 1.133671e-59, 3.04146e-59, 5.482667e-59, 8.511622e-59, 1.156619e-58, 1.692136e-58, 2.305872e-58, 3.041818e-58, 3.967089e-58, 4.902379e-58, 6.138439e-58, 7.445776e-58, 9.136339e-58, 1.089394e-57, 1.287468e-57, 1.512293e-57, 1.756472e-57, 2.043004e-57, 2.342378e-57, 2.669684e-57, 3.026349e-57, 3.424562e-57, 3.847875e-57, 4.301077e-57, 4.782549e-57, 5.3278e-57, 5.887063e-57, 6.479366e-57, 7.112058e-57, 7.798563e-57, 8.523506e-57, 9.279985e-57, 1.006922e-56, 1.091816e-56, 1.180922e-56, 1.282605e-56, 1.396235e-56, 1.530839e-56, 1.654284e-56, 1.783935e-56, 1.931965e-56, 2.085274e-56, 2.261761e-56, 2.441502e-56, 2.624042e-56, 2.865669e-56, 3.058009e-56, 3.256239e-56, 3.475313e-56, 3.714812e-56, 4.030993e-56, 4.280195e-56, 4.552082e-56, 4.811079e-56, 5.073918e-56, 5.433452e-56, 5.738562e-56, 6.049555e-56, 6.362809e-56, 6.723628e-56, 7.047378e-56, 7.396869e-56, 7.805082e-56, 8.176848e-56, 8.561644e-56, 8.945133e-56, 9.333587e-56, 5.371123e-11, 4.202657e-10, 1.385376e-09, 3.118677e-09, 5.747419e-09, 9.371988e-09, 1.407457e-08, 2.293729e-08, 3.723838e-08, 5.541563e-08, 7.796867e-08, 1.052887e-07, 1.378579e-07, 1.761949e-07, 2.277004e-07, 2.969685e-07, 3.80045e-07, 4.761815e-07, 5.859198e-07, 7.094732e-07, 8.489718e-07, 1.005184e-06, 1.176496e-06, 1.396093e-06, 1.640504e-06, 1.904992e-06, 2.202846e-06, 2.52974e-06, 2.885786e-06, 3.270397e-06, 3.715346e-06, 4.216167e-06, 4.735749e-06, 5.323578e-06, 5.972198e-06, 6.662927e-06, 7.396218e-06, 8.181616e-06, 9.046018e-06, 9.956668e-06, 1.096199e-05, 1.205632e-05, 1.32261e-05, 1.447942e-05, 1.580678e-05, 1.72012e-05, 1.863416e-05, 2.027399e-05, 2.200817e-05, 2.384321e-05, 2.575021e-05, 2.775092e-05, 2.97514e-05, 3.198352e-05, 3.430777e-05, 3.683457e-05, 3.965148e-05, 4.25802e-05, 4.583311e-05, 4.916717e-05, 5.271347e-05, 5.591447e-05, 5.96074e-05, 6.340116e-05, 6.793477e-05, 7.308768e-05, 7.702144e-05, 8.230994e-05, 8.784025e-05, 9.355775e-05, 9.899426e-05, 0.0001041236, 0.000110344, 0.0001173295, 0.0001249721, 0.0001321019, 0.0001389155, 0.0001473235, 0.0001553851, 0.0001632849, 0.0001706769, 0.0001795252, 0.0001888101, 0.0002000739, 0.0002081849, 0.0002200094, 0.0002309047, 0.0002412622, 0.0002514368, 0.0002631644, 0.0002740579, 0.0002882209, 0.0002985264, 0.0003138649, 0.0003283826, 0.0003411437, 0.0003551923, 0.0003701007, 0.0003832776, 0.0003964134, 0.0004123707, 0.0004287742, 0.0004478804, 0.0004646016, 0.0004807542, 0.0004985369, 0.0005175098, 0.0005337384, 0.0005496723, 0.0005719689, 0.0005874321, 0.0006073066, 0.0006303471, 0.0006508029, 0.0006707186, 0.0006920362, 0.0007156946, 0.0007309697, 0.0007547821, 0.0007774564, 0.0007989362, 0.0008244339, 0.000848677, 0.0008762973, 0.0008935027, 0.0009212342, 0.0009477174, 0.0009707345, 0.0009967372, 0.001017052, 0.001045449, 0.00106979, 0.001102415, 0.001121683, 0.001154247, 0.001183912, 0.001210337, 0.001240097, 0.001262689, 0.001293504, 0.001320492, 0.001352803, 0.001376302, 0.001410062, 0.001439331, 0.001473081, 0.001503702, 0.001531752, 0.001566957, 0.001595952, 0.001631249, 0.001655417, 0.001693823, 0.001723527, 0.0017608, 0.001787585, 0.001826354, 0.001851603, 0.001892628, 0.001925478, 0.001962594, 0.001990293, 0.002027202, 0.002061052, 0.002091832, 0.00213179, 0.002163765, 0.00220777, 0.002235351, 0.002276053, 0.002305201, 0.002346333, 0.002379541, 0.002418162, 0.002450326, 0.002489516, 0.002516815, 0.002565064, 0.002594515, 0.002637578, 0.002671903, 0.002714582, 0.002748214, 0.002792689, 0.002819696, 0.00286074, 0.002887034, 0.002941259, 0.00296905, 0.00302076, 0.003049944, 0.00309633, 0.003125842, 0.003173093, 0.003197527, 0.003247815, 0.003330045, 0.003328686, 0.003468403, 0.003399296, 0.00354755, 0.003484725, 0.003632939, 0.003562915, 0.003702544, 0.003646956, 0.003787455, 0.003710058, 0.003866776, 0.00379554, 0.003945311, 0.003875006, 0.004021394, 0.003951629, 0.004108704, 0.0040161, 0.004185309, 0.00409876, 0.004268915, 0.004187639, 0.004340171, 0.004254486, 0.004416402, 0.004324658, 0.00450362, 0.004414772, 0.004579758, 0.004626457, 0.0045352, 0.004697437, 0.004617878, 0.004768563, 0.00468903, 0.004858139, 0.004757198, 0.004812435, 0.004961679, 0.004893055, 0.005051376, 0.004965392, 0.005129173, 0.005038966, 0.005075267, 0.005237594, 0.005154961, 0.005304144, 0.005236862, 0.005388911, 0.005424495, 0.0053433, 0.005502296, 0.005412384, 0.005578324, 0.005490195, 0.005530647, 0.005698819, 0.005600233, 0.005768286, 0.005811468, 0.005713983, 0.005873191, 0.005784326, 0.005817301, 0.005987414, 0.005891454, 0.006046966, 0.005971402, 0.006002421, 0.006154267, 0.006073744, 0.006112094, 0.006273024, 0.006173028, 0.006344919, 0.00623759, 0.006282781, 0.006450355, 0.006327697, 0.006395982, 0.006543566, 0.00643683, 0.006485801, 0.006646549, 0.006543083, 0.006575768, 0.006750295, 0.006642499, 0.006818918, 0.006847752, 0.006726343, 0.006920208, 0.006936268, 0.006825047, 0.00700675, 0.007038913, 0.006930579, 0.007096958, 0.007137144, 0.007022466, 0.007203198, 0.0072253, 0.007100856, 0.007149075, 0.007317302, 0.007191971, 0.007227906, 0.007404257, 0.007285241, 0.007317483, 0.007482837, 0.007372061, 0.007409181, 0.007569419, 0.007445977, 0.007496887, 0.007662874, 0.007684472, 0.007554872, 0.007592709, 0.007747921, 0.007653835, 0.007668855, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1.530511e-59, 3.381183e-59, 5.904343e-59, 9.172537e-59, 1.352581e-58, 1.983783e-58, 2.761292e-58, 3.696975e-58, 4.800052e-58, 6.107719e-58, 7.710627e-58, 9.496628e-58, 1.164364e-57, 1.401395e-57, 1.662788e-57, 1.963276e-57, 2.297927e-57, 2.661229e-57, 3.067723e-57, 3.510064e-57, 3.994813e-57, 4.530366e-57, 5.09898e-57, 5.722216e-57, 6.381596e-57, 7.100253e-57, 7.858104e-57, 8.67745e-57, 9.55097e-57, 1.046458e-56, 1.144971e-56, 1.257388e-56, 1.380962e-56, 1.52625e-56, 1.676571e-56, 1.819206e-56, 1.981605e-56, 2.151841e-56, 2.390279e-56, 2.586753e-56, 2.837994e-56, 3.044307e-56, 3.29086e-56, 3.576799e-56, 3.868249e-56, 4.178631e-56, 4.482014e-56, 4.756982e-56, 5.071971e-56, 5.461664e-56, 5.82101e-56, 6.152733e-56, 6.526043e-56, 6.959194e-56, 7.342518e-56, 7.749625e-56, 8.147596e-56, 8.619649e-56, 9.099838e-56, 9.515789e-56, 8.722522e-11, 6.148668e-10, 1.859407e-09, 3.971044e-09, 7.0878e-09, 1.133005e-08, 1.915657e-08, 3.20251e-08, 5.325095e-08, 7.545214e-08, 1.037161e-07, 1.381658e-07, 1.8042e-07, 2.534782e-07, 3.281309e-07, 4.183068e-07, 5.28605e-07, 6.509778e-07, 7.895429e-07, 9.668476e-07, 1.163834e-06, 1.392478e-06, 1.646246e-06, 1.933454e-06, 2.249401e-06, 2.619747e-06, 3.024004e-06, 3.490589e-06, 4.011278e-06, 4.584745e-06, 5.216649e-06, 5.903591e-06, 6.646426e-06, 7.462481e-06, 8.399013e-06, 9.382143e-06, 1.044512e-05, 1.159577e-05, 1.282809e-05, 1.415944e-05, 1.560574e-05, 1.716553e-05, 1.881615e-05, 2.059647e-05, 2.25015e-05, 2.456453e-05, 2.675697e-05, 2.907378e-05, 3.157057e-05, 3.439332e-05, 3.751856e-05, 4.07557e-05, 4.417796e-05, 4.772081e-05, 5.150643e-05, 5.546634e-05, 5.933374e-05, 6.398867e-05, 6.931442e-05, 7.479518e-05, 8.051159e-05, 8.632156e-05, 9.251134e-05, 9.899853e-05, 0.0001057235, 0.000112415, 0.0001204369, 0.0001289716, 0.0001376564, 0.0001467737, 0.0001561452, 0.000165832, 0.0001732594, 0.0001844186, 0.0001967126, 0.0002092894, 0.0002223492, 0.0002351366, 0.0002459424, 0.0002584451, 0.0002723835, 0.0002876366, 0.0003039988, 0.0003215638, 0.0003332834, 0.0003512483, 0.0003695707, 0.000387434, 0.0004048053, 0.0004236833, 0.0004403947, 0.0004631911, 0.0004848734, 0.0005087265, 0.0005298297, 0.0005486829, 0.0005719879, 0.0005935637, 0.0006213289, 0.0006471655, 0.0006687132, 0.0006975191, 0.0007241028, 0.0007535455, 0.0007759042, 0.0008050393, 0.0008340371, 0.0008638466, 0.0008978591, 0.0009244189, 0.000959793, 0.0009920747, 0.001024149, 0.001053346, 0.001089382, 0.001126266, 0.001159899, 0.001197274, 0.001232237, 0.001271876, 0.001309385, 0.001344067, 0.001386392, 0.001428466, 0.001465887, 0.00150674, 0.001544886, 0.001593196, 0.001630946, 0.001676055, 0.001724943, 0.001760563, 0.001809061, 0.001850157, 0.001896195, 0.001947499, 0.001993952, 0.002041292, 0.002088191, 0.002144346, 0.002191288, 0.002234642, 0.002293453, 0.002332023, 0.002388926, 0.002442317, 0.00249675, 0.00254731, 0.002596073, 0.002661998, 0.002698651, 0.00276926, 0.002822976, 0.002874569, 0.002928577, 0.002979884, 0.003039335, 0.003094252, 0.00316249, 0.00320479, 0.003279594, 0.003330082, 0.003391241, 0.003448894, 0.003506551, 0.003570412, 0.003622387, 0.003695493, 0.003736815, 0.00381898, 0.003876914, 0.003938253, 0.004003251, 0.004065019, 0.004126455, 0.004172604, 0.004251587, 0.004298802, 0.004385836, 0.004435243, 0.004509272, 0.004568592, 0.004631904, 0.004707378, 0.00476608, 0.004833712, 0.004895831, 0.004957867, 0.005035629, 0.005092069, 0.005162624, 0.005237718, 0.005287024, 0.00536714, 0.005416233, 0.005506055, 0.005551954, 0.005628798, 0.005690553, 0.005761633, 0.005832562, 0.005916187, 0.005975195, 0.006054496, 0.006106977, 0.00617952, 0.006246825, 0.006304596, 0.006386654, 0.006442319, 0.006522254, 0.006586897, 0.006661, 0.006736651, 0.006791738, 0.006871134, 0.006933412, 0.006990928, 0.007076261, 0.007132084, 0.007202066, 0.00727776, 0.007337033, 0.007414942, 0.007485816, 0.007549946, 0.007750219, 0.007674654, 0.008010439, 0.007820282, 0.00813259, 0.007964405, 0.008279533, 0.0080904, 0.008413441, 0.00822335, 0.008555814, 0.008370379, 0.008680456, 0.008506222, 0.00881159, 0.008639982, 0.00895849, 0.008761585, 0.009089389, 0.008900129, 0.009231284, 0.009036446, 0.009360149, 0.009188338, 0.009502162, 0.009272037, 0.009623784, 0.009411439, 0.009753879, 0.009544469, 0.009882038, 0.009954506, 0.009743377, 0.01006187, 0.00988727, 0.01021299, 0.009995461, 0.01033514, 0.01012739, 0.01018299, 0.01053051, 0.01030941, 0.01066366, 0.01043906, 0.01076851, 0.01057695, 0.01063295, 0.01097008, 0.0107462, 0.01110396, 0.01085955, 0.01122171, 0.01127772, 0.01103193, 0.01140932, 0.01116543, 0.01150399, 0.01129425, 0.01134031, 0.01168927, 0.01145193, 0.01180557, 0.01187178, 0.01162293, 0.01195226, 0.01176182, 0.01180285, 0.01214714, 0.01190458, 0.01225892, 0.01201206, 0.01207009, 0.01242934, 0.01217573, 0.0122616, 0.01260057, 0.01233553, 0.01271522, 0.01244453, 0.01250793, 0.01285691, 0.0126115, 0.0126687, 0.01300948, 0.01275709, 0.01281542, 0.01318024, 0.01290818, 0.01296187, 0.01333149, 0.01305682, 0.01341845, 0.01346906, 0.01321087, 0.01357814, 0.01360469, 0.01334448, 0.01371242, 0.01377134, 0.01349054, 0.0138209, 0.01391144, 0.0136177, 0.01398165, 0.01402503, 0.0137601, 0.01379724, 0.01417196, 0.01388486, 0.0139211, 0.01431282, 0.01401191, 0.01402966, 0.01444369, 0.01415013, 0.01417869, 0.01453956, 0.01426397, 0.01430249, 0.01467563, 0.01470097, 0.01440903, 0.01446625, 0.01480893, 0.01454562, 0.0145462, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 6.761826e-60, 2.50041e-59, 4.67105e-59, 8.159767e-59, 1.283952e-58, 1.95575e-58, 2.808545e-58, 3.854212e-58, 5.114007e-58, 6.670236e-58, 8.535484e-58, 1.07097e-57, 1.326117e-57, 1.60667e-57, 1.933111e-57, 2.311375e-57, 2.724447e-57, 3.174981e-57, 3.691767e-57, 4.256196e-57, 4.874601e-57, 5.555868e-57, 6.280253e-57, 7.07885e-57, 7.93637e-57, 8.85843e-57, 9.852402e-57, 1.089758e-56, 1.203476e-56, 1.356794e-56, 1.511399e-56, 1.670482e-56, 1.824216e-56, 2.006622e-56, 2.233392e-56, 2.480685e-56, 2.735677e-56, 2.999752e-56, 3.255796e-56, 3.580621e-56, 3.895174e-56, 4.235273e-56, 4.568674e-56, 4.892448e-56, 5.286244e-56, 5.735403e-56, 6.114394e-56, 6.562099e-56, 6.942531e-56, 7.424134e-56, 7.919514e-56, 8.403991e-56, 8.895855e-56, 9.378395e-56, 7.240348e-11, 6.062695e-10, 1.87975e-09, 4.095752e-09, 7.399505e-09, 1.347856e-08, 2.363177e-08, 4.307753e-08, 6.411336e-08, 9.548895e-08, 1.289034e-07, 1.780967e-07, 2.3787e-07, 3.210402e-07, 4.135589e-07, 5.275748e-07, 6.543005e-07, 8.094349e-07, 1.008144e-06, 1.243531e-06, 1.500516e-06, 1.799322e-06, 2.136604e-06, 2.515964e-06, 2.958731e-06, 3.45889e-06, 4.025022e-06, 4.625694e-06, 5.282247e-06, 6.058066e-06, 6.905427e-06, 7.814415e-06, 8.804318e-06, 9.915259e-06, 1.114269e-05, 1.241624e-05, 1.381013e-05, 1.536387e-05, 1.712845e-05, 1.896805e-05, 2.096649e-05, 2.306251e-05, 2.531817e-05, 2.77578e-05, 3.06255e-05, 3.368479e-05, 3.690299e-05, 4.059866e-05, 4.418738e-05, 4.800071e-05, 5.211829e-05, 5.671943e-05, 6.18471e-05, 6.785487e-05, 7.349571e-05, 7.940948e-05, 8.560978e-05, 9.238381e-05, 9.982433e-05, 0.0001079946, 0.0001164244, 0.0001262363, 0.0001352456, 0.0001443018, 0.0001540588, 0.0001648757, 0.0001769468, 0.000189298, 0.0002018563, 0.0002149639, 0.00022838, 0.0002423747, 0.0002579166, 0.0002749681, 0.0002909115, 0.0003081249, 0.0003263676, 0.0003450472, 0.0003639985, 0.0003835916, 0.0004049201, 0.0004284031, 0.0004492693, 0.0004725611, 0.0004955528, 0.000520341, 0.00054632, 0.0005737526, 0.0006009362, 0.0006305006, 0.0006624273, 0.0006875089, 0.0007193065, 0.0007519834, 0.0007843762, 0.0008202688, 0.0008531966, 0.0008904614, 0.0009232331, 0.0009628047, 0.001002144, 0.001044499, 0.001082346, 0.00112327, 0.001159691, 0.001207003, 0.00125445, 0.001303069, 0.001349187, 0.001397325, 0.001438267, 0.001486158, 0.001540419, 0.001597906, 0.001650734, 0.001704355, 0.00175209, 0.001807615, 0.001870348, 0.001922299, 0.00198188, 0.002037394, 0.002100044, 0.00216792, 0.002229485, 0.002289657, 0.002352616, 0.002419706, 0.002484892, 0.002552385, 0.002620237, 0.002687103, 0.002770199, 0.002833785, 0.00291257, 0.002980896, 0.003048378, 0.003133906, 0.003198172, 0.003284964, 0.003364244, 0.003433979, 0.003529696, 0.003602812, 0.003681833, 0.003767088, 0.003837389, 0.003933402, 0.004008908, 0.004097868, 0.004195371, 0.00427478, 0.004362447, 0.004457403, 0.004540723, 0.004629673, 0.004713162, 0.004814965, 0.004910578, 0.004993755, 0.005102444, 0.005179808, 0.005287167, 0.005380087, 0.005472824, 0.005579049, 0.005667033, 0.005772154, 0.005862298, 0.005980905, 0.006077001, 0.006171191, 0.006277639, 0.006371463, 0.006485868, 0.006569078, 0.006690217, 0.006779783, 0.00690596, 0.007013129, 0.00710663, 0.007225226, 0.007307775, 0.00743717, 0.0075281, 0.007644936, 0.007743824, 0.007855643, 0.007980352, 0.008069074, 0.008195291, 0.008301543, 0.008419197, 0.008522468, 0.00862429, 0.008739102, 0.008850138, 0.008970245, 0.009050066, 0.009185566, 0.009288441, 0.009414972, 0.009519822, 0.009639017, 0.009753594, 0.009847701, 0.009992516, 0.01007795, 0.01021179, 0.0103053, 0.01044012, 0.01053783, 0.01065807, 0.010765, 0.01089447, 0.0110022, 0.0111299, 0.01121932, 0.01134769, 0.01146453, 0.01157003, 0.01168443, 0.01178176, 0.01193688, 0.01202236, 0.01214956, 0.01225679, 0.01235806, 0.01248413, 0.01260581, 0.01271227, 0.01283232, 0.01291672, 0.0130615, 0.01314772, 0.01328998, 0.01338573, 0.01350634, 0.01361204, 0.01373271, 0.01383794, 0.0139636, 0.01403679, 0.0141806, 0.0142766, 0.01440337, 0.01450988, 0.0146202, 0.0147337, 0.01485535, 0.01492934, 0.01507279, 0.01514447, 0.01551898, 0.01537711, 0.01600476, 0.01560395, 0.01623364, 0.01583173, 0.01644367, 0.01600436, 0.01668324, 0.01624232, 0.01688798, 0.01646501, 0.01709077, 0.01668418, 0.01730656, 0.01686452, 0.01754226, 0.017105, 0.01771445, 0.01731664, 0.01792251, 0.01750621, 0.01811447, 0.01773531, 0.01831166, 0.01788083, 0.01848307, 0.01812067, 0.01870806, 0.01829005, 0.0189309, 0.01904918, 0.01860331, 0.01924501, 0.01879235, 0.01942383, 0.01898434, 0.01957667, 0.01919402, 0.01926313, 0.01991019, 0.01947266, 0.02009412, 0.01965272, 0.02027451, 0.01981376, 0.01994738, 0.0205964, 0.02005864, 0.02075416, 0.02026143, 0.02090857, 0.0210395, 0.0205559, 0.02119084, 0.02071986, 0.02135314, 0.02089054, 0.02098172, 0.02158624, 0.02117919, 0.02180409, 0.02190341, 0.02138457, 0.02202783, 0.02155973, 0.02164981, 0.02232289, 0.02178719, 0.02248432, 0.02198016, 0.02201759, 0.02270506, 0.02219359, 0.02230083, 0.02295721, 0.02243535, 0.02308332, 0.0225774, 0.02264628, 0.02333274, 0.02282215, 0.02288267, 0.02350093, 0.0230564, 0.02309389, 0.02375361, 0.02322141, 0.02331069, 0.02402477, 0.02345755, 0.02408963, 0.0242341, 0.02367343, 0.02431951, 0.0243895, 0.0238821, 0.02451123, 0.02462177, 0.02404094, 0.02470585, 0.02481517, 0.02423572, 0.02490186, 0.02500809, 0.02445077, 0.02448368, 0.02517214, 0.02460818, 0.0246519, 0.02532744, 0.02476862, 0.02484609, 0.02551306, 0.02496741, 0.0249857, 0.02566757, 0.02512383, 0.02518115, 0.02584167, 0.02590434, 0.0253333, 0.02537382, 0.0260256, 0.02549195, 0.02552717, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 9.095012e-60, 3.031298e-59, 5.789377e-59, 9.988601e-59, 1.643502e-58, 2.511711e-58, 3.574384e-58, 4.924648e-58, 6.61255e-58, 8.706968e-58, 1.114447e-57, 1.398484e-57, 1.730841e-57, 2.118454e-57, 2.551254e-57, 3.047301e-57, 3.597614e-57, 4.221818e-57, 4.896274e-57, 5.647367e-57, 6.484373e-57, 7.388686e-57, 8.366416e-57, 9.433982e-57, 1.057867e-56, 1.186565e-56, 1.355315e-56, 1.51809e-56, 1.702498e-56, 1.885077e-56, 2.126786e-56, 2.399207e-56, 2.657835e-56, 2.951836e-56, 3.258012e-56, 3.596326e-56, 3.973326e-56, 4.330871e-56, 4.699972e-56, 5.122712e-56, 5.552879e-56, 6.009141e-56, 6.475541e-56, 6.954476e-56, 7.471335e-56, 8.029326e-56, 8.569342e-56, 9.13496e-56, 8.38255e-11, 6.883267e-10, 2.148188e-09, 4.693413e-09, 8.46346e-09, 1.779229e-08, 3.233427e-08, 5.191913e-08, 7.465677e-08, 1.061059e-07, 1.475674e-07, 2.0907e-07, 2.819846e-07, 3.732092e-07, 4.755367e-07, 6.092214e-07, 7.884348e-07, 9.782089e-07, 1.223746e-06, 1.506091e-06, 1.815711e-06, 2.216264e-06, 2.647714e-06, 3.160257e-06, 3.728079e-06, 4.345408e-06, 5.048378e-06, 5.828422e-06, 6.733416e-06, 7.683434e-06, 8.76147e-06, 9.95037e-06, 1.127793e-05, 1.278572e-05, 1.438366e-05, 1.617381e-05, 1.805857e-05, 2.007763e-05, 2.231392e-05, 2.489529e-05, 2.775709e-05, 3.098501e-05, 3.420945e-05, 3.795991e-05, 4.169878e-05, 4.598046e-05, 5.086639e-05, 5.648782e-05, 6.177011e-05, 6.743477e-05, 7.374625e-05, 8.033546e-05, 8.778979e-05, 9.540604e-05, 0.0001049512, 0.0001131771, 0.000121886, 0.0001310397, 0.0001419921, 0.0001548548, 0.000166318, 0.0001784745, 0.0001911165, 0.0002054841, 0.0002208625, 0.0002365476, 0.0002530571, 0.0002697447, 0.0002869822, 0.0003065884, 0.0003261438, 0.0003462591, 0.0003681292, 0.0003901917, 0.0004131271, 0.0004391542, 0.0004633261, 0.0004883539, 0.0005155242, 0.0005466703, 0.0005769462, 0.0006079066, 0.0006423711, 0.0006746345, 0.0007070329, 0.000740965, 0.0007782243, 0.0008157311, 0.0008553172, 0.0008957758, 0.0009364339, 0.0009788551, 0.001024112, 0.001071714, 0.001121013, 0.001169485, 0.001218515, 0.001264698, 0.001316864, 0.00137048, 0.00142637, 0.001487528, 0.001545743, 0.001608373, 0.001668052, 0.001727467, 0.001784769, 0.001853302, 0.001925663, 0.001995684, 0.002072294, 0.002144025, 0.002212823, 0.002288088, 0.002355652, 0.002435013, 0.002520266, 0.002607323, 0.002689685, 0.002777268, 0.002859899, 0.00293341, 0.003024813, 0.003120092, 0.003213967, 0.003320406, 0.003415564, 0.003493857, 0.003596217, 0.003701106, 0.003796563, 0.003911319, 0.004021107, 0.004106696, 0.00422415, 0.004336416, 0.00444717, 0.004564481, 0.004668259, 0.004779936, 0.004906114, 0.005018417, 0.005151969, 0.005275548, 0.005376585, 0.005508174, 0.005641301, 0.005758828, 0.005897858, 0.006010626, 0.006146702, 0.006296811, 0.006426195, 0.006560037, 0.006697157, 0.006823318, 0.006972588, 0.007113558, 0.007241543, 0.007399796, 0.007541785, 0.007693984, 0.007842337, 0.007972673, 0.008124163, 0.008277701, 0.008413105, 0.008577193, 0.008720669, 0.008891205, 0.009050293, 0.009197018, 0.009369312, 0.009497398, 0.009672029, 0.009831094, 0.009982384, 0.01016713, 0.01031172, 0.01047362, 0.01066378, 0.01080812, 0.01098097, 0.01112827, 0.01131365, 0.01148325, 0.01163378, 0.01183038, 0.01197505, 0.01216211, 0.01233864, 0.01251268, 0.01267741, 0.01284961, 0.01303767, 0.01317083, 0.01340218, 0.01354637, 0.0137288, 0.01390769, 0.01407499, 0.01426245, 0.0144465, 0.01463966, 0.01478641, 0.01498874, 0.01514789, 0.01533668, 0.01553717, 0.01569038, 0.01588757, 0.0160593, 0.01625797, 0.0164212, 0.01660515, 0.01679223, 0.01695614, 0.01716104, 0.01729526, 0.01753674, 0.01769914, 0.01788458, 0.01806608, 0.0182563, 0.01841985, 0.01856382, 0.01879569, 0.01893749, 0.01915147, 0.01931651, 0.01952136, 0.01968704, 0.01985742, 0.02007484, 0.02024296, 0.02040518, 0.02059983, 0.02076831, 0.0209527, 0.02113296, 0.0213167, 0.02150723, 0.0216563, 0.02186721, 0.02200724, 0.02222381, 0.022367, 0.02256992, 0.02272798, 0.02289993, 0.02311142, 0.02327219, 0.02344973, 0.02365572, 0.02379468, 0.02397617, 0.0241505, 0.02432113, 0.02448392, 0.02464222, 0.02485115, 0.02501943, 0.02516554, 0.02537609, 0.02552046, 0.02569468, 0.02585566, 0.02601219, 0.0262169, 0.02632841, 0.02656803, 0.02669147, 0.0268385, 0.02702796, 0.02719585, 0.0273552, 0.02751479, 0.02808223, 0.02787438, 0.02887575, 0.02812888, 0.02920971, 0.02846692, 0.02950789, 0.02878027, 0.02982529, 0.02910651, 0.0301477, 0.02942193, 0.03050098, 0.02970668, 0.03078707, 0.03003013, 0.03109196, 0.03033068, 0.03139897, 0.03063744, 0.03169064, 0.03088886, 0.0320398, 0.03126806, 0.03226167, 0.0314505, 0.03258194, 0.03179465, 0.03285079, 0.03205893, 0.03314449, 0.03329793, 0.0324313, 0.03359689, 0.0327618, 0.03386641, 0.03302575, 0.03414808, 0.0333267, 0.03342784, 0.03452603, 0.03372342, 0.03482496, 0.03392258, 0.03509926, 0.03423178, 0.03435859, 0.03545354, 0.03460853, 0.03572938, 0.03484856, 0.0359476, 0.03611996, 0.03524199, 0.03632999, 0.03544003, 0.03659492, 0.0356906, 0.03584302, 0.03695087, 0.03607614, 0.03721649, 0.03730299, 0.03638258, 0.03755813, 0.03662444, 0.03672168, 0.03788105, 0.0369542, 0.0381077, 0.03718531, 0.0372982, 0.03846151, 0.03750948, 0.03759093, 0.03879409, 0.03784305, 0.03897091, 0.03799745, 0.03814362, 0.03926548, 0.03831073, 0.03841441, 0.03957766, 0.03863494, 0.03871438, 0.03981878, 0.03894178, 0.03903766, 0.04016595, 0.03920858, 0.04036777, 0.04043245, 0.03944401, 0.04063082, 0.0407142, 0.0397398, 0.04090596, 0.04094338, 0.03997838, 0.04120951, 0.04124815, 0.04020864, 0.04141601, 0.04151119, 0.04045631, 0.04052616, 0.04176543, 0.04074522, 0.04079115, 0.04196887, 0.04087799, 0.04102796, 0.04220329, 0.04115606, 0.04120842, 0.04243914, 0.04136197, 0.04142328, 0.04266164, 0.04271012, 0.04167533, 0.04171877, 0.04283293, 0.04189498, 0.04196528, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1.219706e-59, 3.744084e-59, 7.2429e-59, 1.250128e-58, 2.101972e-58, 3.196727e-58, 4.557609e-58, 6.346127e-58, 8.541541e-58, 1.121431e-57, 1.435796e-57, 1.812504e-57, 2.255656e-57, 2.754519e-57, 3.328432e-57, 3.976116e-57, 4.701696e-57, 5.522248e-57, 6.422411e-57, 7.416295e-57, 8.505097e-57, 9.694171e-57, 1.103561e-56, 1.279655e-56, 1.458084e-56, 1.646691e-56, 1.867191e-56, 2.138793e-56, 2.441686e-56, 2.726288e-56, 3.044748e-56, 3.42858e-56, 3.804134e-56, 4.204362e-56, 4.625696e-56, 5.065788e-56, 5.572677e-56, 6.072919e-56, 6.575886e-56, 7.143215e-56, 7.734556e-56, 8.316295e-56, 2.185465e-12, 2.048852e-10, 1.002232e-09, 2.627818e-09, 5.252212e-09, 1.151333e-08, 2.187393e-08, 3.818976e-08, 5.830973e-08, 8.420952e-08, 1.306461e-07, 1.869302e-07, 2.588855e-07, 3.49378e-07, 4.718829e-07, 6.356416e-07, 8.348753e-07, 1.068201e-06, 1.344717e-06, 1.660891e-06, 2.045534e-06, 2.485907e-06, 2.99475e-06, 3.567403e-06, 4.21067e-06, 4.960137e-06, 5.794921e-06, 6.724794e-06, 7.754751e-06, 8.885403e-06, 1.018323e-05, 1.160617e-05, 1.318641e-05, 1.49007e-05, 1.678601e-05, 1.895895e-05, 2.1423e-05, 2.418221e-05, 2.72127e-05, 3.026606e-05, 3.383208e-05, 3.818389e-05, 4.254655e-05, 4.775816e-05, 5.306781e-05, 5.845215e-05, 6.485055e-05, 7.142478e-05, 7.938122e-05, 8.670396e-05, 9.531018e-05, 0.0001033708, 0.0001133787, 0.0001242454, 0.000135463, 0.0001477484, 0.0001597515, 0.0001727816, 0.0001862284, 0.0002026767, 0.0002185645, 0.0002357469, 0.0002538484, 0.0002706332, 0.00029044, 0.0003125003, 0.0003364577, 0.0003582383, 0.0003822625, 0.0004052335, 0.0004326026, 0.0004619398, 0.0004918313, 0.0005202387, 0.0005528695, 0.000585315, 0.0006177929, 0.0006536719, 0.0006893525, 0.0007337072, 0.0007709153, 0.0008095205, 0.0008493675, 0.0008937993, 0.0009461762, 0.0009953691, 0.001043795, 0.001092007, 0.001142196, 0.001199424, 0.001255621, 0.001315928, 0.001373455, 0.001433727, 0.001494692, 0.001561967, 0.001627844, 0.001699499, 0.001771102, 0.001847334, 0.001923166, 0.0019985, 0.002079776, 0.002162863, 0.002247075, 0.002328184, 0.002413783, 0.002503894, 0.002597035, 0.002689017, 0.002786614, 0.002891098, 0.00299017, 0.003098932, 0.003202299, 0.003310721, 0.003417163, 0.003518583, 0.003635327, 0.003758286, 0.003882108, 0.004007797, 0.004128897, 0.004260919, 0.004384706, 0.004515174, 0.004646501, 0.004782641, 0.004924008, 0.005066602, 0.005218655, 0.005357815, 0.005505569, 0.005645422, 0.005800704, 0.005954277, 0.006120759, 0.006283608, 0.006439974, 0.006616354, 0.006778439, 0.006928798, 0.007102976, 0.007287638, 0.00745758, 0.00765226, 0.00782802, 0.008001671, 0.008171255, 0.008364529, 0.00855308, 0.008766837, 0.008962276, 0.009146706, 0.009337856, 0.009534993, 0.009733588, 0.009951933, 0.01015481, 0.01036481, 0.0105687, 0.0107726, 0.01100204, 0.01121782, 0.01142999, 0.01162927, 0.01187777, 0.01207945, 0.01233966, 0.01256603, 0.01276845, 0.01301665, 0.01322503, 0.01347138, 0.01371837, 0.01392294, 0.01417422, 0.01443587, 0.01465332, 0.01489462, 0.01514992, 0.01539058, 0.01563399, 0.01587383, 0.01612657, 0.01639705, 0.01663479, 0.01690658, 0.01712937, 0.01740093, 0.01767442, 0.01790422, 0.01816423, 0.01843558, 0.01870395, 0.01895627, 0.01922033, 0.01950337, 0.01975973, 0.02003644, 0.0202953, 0.02055099, 0.02083418, 0.02108953, 0.02135112, 0.02167067, 0.02188891, 0.02218998, 0.02242462, 0.02273855, 0.02299191, 0.02327877, 0.02354759, 0.02378385, 0.0241185, 0.02437842, 0.02465545, 0.02492356, 0.02521043, 0.0255109, 0.0257238, 0.02605465, 0.02628783, 0.02657594, 0.02689391, 0.02716592, 0.02745932, 0.02770247, 0.02801182, 0.02824099, 0.02856901, 0.02881436, 0.02909755, 0.02940688, 0.02964835, 0.02994485, 0.03017391, 0.03052155, 0.03075522, 0.03102888, 0.03130375, 0.03160399, 0.03188444, 0.03210853, 0.032422, 0.03266758, 0.03296224, 0.03321909, 0.03351185, 0.03379284, 0.03401609, 0.0343378, 0.0345687, 0.03489063, 0.0351117, 0.03540579, 0.03565285, 0.03591754, 0.0362102, 0.0364584, 0.0367489, 0.03702575, 0.03722684, 0.0375106, 0.03779088, 0.03803318, 0.03830963, 0.03851733, 0.03885327, 0.03905738, 0.03936359, 0.0395874, 0.03981815, 0.04005688, 0.04037169, 0.0406045, 0.04085995, 0.04107101, 0.04138092, 0.04157675, 0.04183442, 0.04208712, 0.04233191, 0.04257667, 0.04282241, 0.04305737, 0.04331145, 0.04345897, 0.04380376, 0.04399183, 0.04428426, 0.04449718, 0.04473339, 0.04496313, 0.04519971, 0.04544011, 0.04566253, 0.04583489, 0.04609906, 0.04705916, 0.046572, 0.04824492, 0.04701156, 0.04871429, 0.04744338, 0.04907479, 0.04792049, 0.0495822, 0.04831148, 0.05006036, 0.04874094, 0.05048619, 0.04916726, 0.05086307, 0.04963332, 0.05133105, 0.049992, 0.05179329, 0.05043495, 0.05215401, 0.05082105, 0.05256695, 0.05128652, 0.05295842, 0.05153617, 0.0533739, 0.05195194, 0.05377369, 0.05238944, 0.05414308, 0.05433617, 0.0529483, 0.05470906, 0.05332256, 0.05507102, 0.05371236, 0.05550298, 0.05401515, 0.05418915, 0.05603605, 0.05457197, 0.05638165, 0.05489867, 0.05675279, 0.05529996, 0.05541669, 0.05724515, 0.05581781, 0.05761896, 0.05607113, 0.05795288, 0.05810818, 0.056588, 0.05844551, 0.05689394, 0.05879521, 0.05719704, 0.05728409, 0.05927203, 0.05769608, 0.0595418, 0.05968425, 0.05816826, 0.06000958, 0.05841331, 0.05858674, 0.06052806, 0.05888634, 0.06073448, 0.05913541, 0.05933452, 0.06119154, 0.05956966, 0.05976905, 0.06162725, 0.05999371, 0.061871, 0.0602482, 0.060427, 0.06230894, 0.06058874, 0.06076608, 0.06271502, 0.06104497, 0.0611604, 0.06308378, 0.06138805, 0.06152404, 0.06341073, 0.0617939, 0.0637433, 0.06383053, 0.06207322, 0.06405011, 0.06418276, 0.06244156, 0.06443013, 0.06452651, 0.06280102, 0.06475594, 0.06478766, 0.06306175, 0.06509138, 0.06521918, 0.06340813, 0.06343364, 0.06552146, 0.06375631, 0.06383556, 0.06577645, 0.06397407, 0.06411216, 0.06607324, 0.06429189, 0.06437994, 0.06646041, 0.06457571, 0.06460806, 0.06665367, 0.06674593, 0.06491873, 0.06501023, 0.06699677, 0.06523856, 0.06527158, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 3.155122e-60, 2.251618e-59, 5.374217e-59, 1.023585e-58, 1.821036e-58, 2.909282e-58, 4.341688e-58, 6.24955e-58, 8.665239e-58, 1.158399e-57, 1.518751e-57, 1.944749e-57, 2.443454e-57, 3.026631e-57, 3.699613e-57, 4.460933e-57, 5.321788e-57, 6.294172e-57, 7.376265e-57, 8.564159e-57, 9.9156e-57, 1.168145e-56, 1.356941e-56, 1.554528e-56, 1.807937e-56, 2.105446e-56, 2.40322e-56, 2.717222e-56, 3.093784e-56, 3.503398e-56, 3.905817e-56, 4.350637e-56, 4.86189e-56, 5.383109e-56, 5.926439e-56, 6.484032e-56, 7.058685e-56, 7.717024e-56, 4.695902e-12, 2.834388e-10, 1.30715e-09, 3.350866e-09, 7.206203e-09, 1.674564e-08, 3.104171e-08, 5.094363e-08, 7.860349e-08, 1.227151e-07, 1.802419e-07, 2.498547e-07, 3.373599e-07, 4.58211e-07, 6.073209e-07, 7.793648e-07, 9.907748e-07, 1.247765e-06, 1.555498e-06, 1.92145e-06, 2.343934e-06, 2.828525e-06, 3.396849e-06, 4.080446e-06, 4.831039e-06, 5.683782e-06, 6.681496e-06, 7.809681e-06, 9.077443e-06, 1.052089e-05, 1.211833e-05, 1.391339e-05, 1.617095e-05, 1.862026e-05, 2.134021e-05, 2.43179e-05, 2.769481e-05, 3.170892e-05, 3.59838e-05, 4.056857e-05, 4.546593e-05, 5.098857e-05, 5.751824e-05, 6.431583e-05, 7.142926e-05, 7.878483e-05, 8.704225e-05, 9.640733e-05, 0.0001059004, 0.0001163667, 0.0001276191, 0.0001387518, 0.0001528023, 0.0001670505, 0.0001820006, 0.0001972509, 0.0002128777, 0.0002311982, 0.0002501838, 0.0002721875, 0.0002946109, 0.000315057, 0.0003393275, 0.0003639837, 0.0003930898, 0.0004222176, 0.0004493051, 0.0004803538, 0.0005102094, 0.0005469522, 0.0005804976, 0.000620397, 0.0006561628, 0.0006963872, 0.0007377681, 0.0007839628, 0.0008315179, 0.0008775566, 0.0009249296, 0.0009748863, 0.001031507, 0.001087432, 0.00114494, 0.001204743, 0.001265048, 0.001328706, 0.001394545, 0.001463792, 0.0015377, 0.001613315, 0.001683945, 0.001756043, 0.001844357, 0.001929921, 0.002013206, 0.002100269, 0.002195321, 0.002287325, 0.002379826, 0.002477505, 0.002577312, 0.002689238, 0.0028012, 0.002905863, 0.003013819, 0.003127515, 0.003257938, 0.003371859, 0.003500531, 0.003625816, 0.003752914, 0.003887413, 0.004031746, 0.004171904, 0.004317096, 0.00445587, 0.004613251, 0.004765915, 0.004918454, 0.00508995, 0.005253901, 0.005421531, 0.005577226, 0.005746911, 0.005922085, 0.006103878, 0.00628832, 0.00648446, 0.006678166, 0.006864938, 0.007066721, 0.007261311, 0.007461005, 0.007668621, 0.007870581, 0.008095824, 0.008323533, 0.008548593, 0.008770826, 0.008990571, 0.009225297, 0.009460893, 0.009703414, 0.0099492, 0.01017562, 0.01041406, 0.01067092, 0.0109361, 0.01117832, 0.01146365, 0.01172776, 0.01198693, 0.01224164, 0.01252419, 0.01280573, 0.01308902, 0.01339549, 0.01366673, 0.01394675, 0.0142291, 0.01453359, 0.01483307, 0.01514655, 0.01545251, 0.01576558, 0.01607316, 0.01636987, 0.01669273, 0.01703457, 0.01734841, 0.01766987, 0.01799906, 0.01829466, 0.01864832, 0.0190155, 0.01932017, 0.0196901, 0.02004424, 0.02033295, 0.02069812, 0.02106758, 0.0214139, 0.02179438, 0.02211862, 0.02248555, 0.02285504, 0.02322258, 0.02354644, 0.02392331, 0.02428627, 0.02468251, 0.02505768, 0.02541944, 0.02580243, 0.02620589, 0.02654623, 0.02695489, 0.027314, 0.02770008, 0.02812466, 0.02848255, 0.02889755, 0.02925976, 0.02965399, 0.03006334, 0.03042983, 0.03081052, 0.0312656, 0.03160142, 0.03206726, 0.03240339, 0.03282487, 0.03325681, 0.03362657, 0.03402627, 0.03441422, 0.03484784, 0.03523334, 0.03563201, 0.03604854, 0.03650101, 0.03687125, 0.03730765, 0.0376578, 0.03809043, 0.03850339, 0.03889629, 0.0393385, 0.03967467, 0.04012535, 0.04052709, 0.04094449, 0.04135778, 0.0417771, 0.04217627, 0.04252524, 0.04302465, 0.04338529, 0.04379624, 0.04419628, 0.04461958, 0.0450364, 0.04536684, 0.04591548, 0.04621272, 0.04668103, 0.04703925, 0.04744503, 0.04787888, 0.04825115, 0.04866593, 0.04902853, 0.04948961, 0.0498529, 0.05021968, 0.05063347, 0.05104135, 0.05144101, 0.05176866, 0.05225157, 0.05261818, 0.05300279, 0.05341818, 0.05378251, 0.05415224, 0.0544798, 0.05495373, 0.05525372, 0.05570971, 0.05604189, 0.05645067, 0.05682981, 0.05718569, 0.05761602, 0.05794137, 0.05836065, 0.05867878, 0.05904767, 0.05941526, 0.05981273, 0.06018985, 0.0605354, 0.06090484, 0.06127433, 0.06158502, 0.06200449, 0.06227935, 0.0626713, 0.06302358, 0.0633596, 0.06379301, 0.06407522, 0.0644179, 0.06483889, 0.06511219, 0.06545851, 0.06583994, 0.06616036, 0.06647679, 0.06676748, 0.06721837, 0.0674885, 0.06778457, 0.06814974, 0.06846514, 0.06878142, 0.06911977, 0.06943988, 0.06973287, 0.07003477, 0.07042788, 0.07069211, 0.07101584, 0.07136353, 0.0716057, 0.07196014, 0.07223591, 0.07250594, 0.07401193, 0.07311326, 0.07567477, 0.0738071, 0.07629228, 0.07430543, 0.07686839, 0.0749134, 0.0774587, 0.07545867, 0.07812179, 0.07605204, 0.07858078, 0.07661992, 0.07918728, 0.0771665, 0.07974165, 0.07772234, 0.08033006, 0.07821736, 0.08082162, 0.078843, 0.0813955, 0.0794119, 0.08196617, 0.0796276, 0.08242108, 0.08029192, 0.08293055, 0.08083819, 0.08347027, 0.08362419, 0.08156531, 0.08422214, 0.08201127, 0.08472851, 0.0824865, 0.08517079, 0.08297624, 0.08316811, 0.08597047, 0.08367269, 0.08626309, 0.08410813, 0.08683227, 0.08456635, 0.08474729, 0.08754721, 0.08523775, 0.08794779, 0.0856212, 0.088443, 0.08862522, 0.08620438, 0.08900896, 0.08666845, 0.08943439, 0.08702451, 0.08729437, 0.09005914, 0.0876659, 0.09046275, 0.09073943, 0.08832108, 0.09104871, 0.08857659, 0.08885726, 0.09167419, 0.08918749, 0.09199942, 0.0895389, 0.08975978, 0.0925389, 0.0900148, 0.0903636, 0.09321789, 0.09059133, 0.09342084, 0.09095198, 0.09110912, 0.09400961, 0.09139356, 0.09160702, 0.09453863, 0.09189307, 0.09199871, 0.0950538, 0.09239207, 0.09255431, 0.09542258, 0.09291864, 0.09572394, 0.09596093, 0.09328201, 0.09621761, 0.09638697, 0.09366916, 0.09656022, 0.09687307, 0.09414503, 0.0971438, 0.09722454, 0.09455623, 0.0975142, 0.09767184, 0.09495036, 0.09506527, 0.0981159, 0.09531316, 0.09532718, 0.09851444, 0.09573826, 0.09579925, 0.09883796, 0.09602962, 0.0961327, 0.09916011, 0.09634154, 0.0964898, 0.09957074, 0.09967064, 0.09678402, 0.09679175, 0.09995057, 0.0972113, 0.09725313, 
1e-100, 1e-100, 1e-100, 1e-100, 1e-100, 8.277813e-61, 1.955685e-59, 5.095862e-59, 1.028825e-58, 1.898539e-58, 3.098141e-58, 4.753624e-58, 6.955784e-58, 9.718737e-58, 1.320241e-57, 1.746326e-57, 2.252219e-57, 2.856123e-57, 3.561229e-57, 4.369766e-57, 5.30243e-57, 6.35378e-57, 7.538416e-57, 8.918344e-57, 1.070429e-56, 1.261433e-56, 1.468405e-56, 1.761987e-56, 2.057213e-56, 2.376233e-56, 2.758578e-56, 3.167153e-56, 3.582301e-56, 4.053704e-56, 4.582284e-56, 5.110077e-56, 5.651244e-56, 6.293261e-56, 6.94361e-56, 2.151733e-12, 1.848525e-10, 8.608051e-10, 2.21746e-09, 5.752734e-09, 1.309027e-08, 2.42148e-08, 4.097484e-08, 6.915643e-08, 1.08611e-07, 1.605313e-07, 2.287673e-07, 3.286225e-07, 4.567183e-07, 6.189609e-07, 8.340533e-07, 1.109657e-06, 1.451995e-06, 1.856549e-06, 2.333738e-06, 2.896161e-06, 3.546078e-06, 4.29298e-06, 5.14657e-06, 6.1527e-06, 7.258528e-06, 8.518743e-06, 9.890522e-06, 1.150399e-05, 1.343438e-05, 1.553163e-05, 1.787119e-05, 2.054067e-05, 2.360658e-05, 2.716038e-05, 3.090942e-05, 3.485875e-05, 3.98033e-05, 4.529959e-05, 5.083179e-05, 5.726849e-05, 6.394617e-05, 7.176816e-05, 8.062216e-05, 8.973109e-05, 9.91624e-05, 0.0001103012, 0.0001230234, 0.0001358553, 0.0001493132, 0.0001637287, 0.0001809693, 0.0001985033, 0.0002171908, 0.0002358359, 0.0002553205, 0.0002779722, 0.0003024573, 0.0003274865, 0.0003536593, 0.0003805237, 0.0004102102, 0.0004421566, 0.0004745683, 0.0005095271, 0.0005436378, 0.000582515, 0.0006231984, 0.000664345, 0.000707352, 0.0007492164, 0.000798668, 0.0008527645, 0.0009032827, 0.0009604633, 0.001011608, 0.001073358, 0.001136965, 0.001200912, 0.001273945, 0.001339319, 0.001409999, 0.001480471, 0.001567814, 0.001651815, 0.001731787, 0.001813874, 0.001903972, 0.002000461, 0.002093291, 0.002194331, 0.002297324, 0.00240279, 0.002507988, 0.002620816, 0.002745965, 0.00286654, 0.002979913, 0.003102157, 0.00323933, 0.003374822, 0.003511876, 0.003653865, 0.003807435, 0.003954143, 0.004107226, 0.004261391, 0.004429626, 0.004603075, 0.004766407, 0.00492944, 0.005105557, 0.005302551, 0.005491157, 0.005680236, 0.005867344, 0.006075651, 0.006268487, 0.006480798, 0.006696502, 0.006905238, 0.007138637, 0.007364639, 0.007594472, 0.007822611, 0.008063845, 0.008320871, 0.008577053, 0.008832842, 0.009088427, 0.009337921, 0.009607422, 0.009882467, 0.01018337, 0.01047284, 0.01074293, 0.0110267, 0.01132112, 0.01161115, 0.01193928, 0.01224665, 0.0125698, 0.01290046, 0.01320634, 0.01354763, 0.01388901, 0.01422412, 0.01456951, 0.01494636, 0.01528636, 0.01563037, 0.01599448, 0.01637456, 0.01673811, 0.01713802, 0.01753795, 0.01791071, 0.01829544, 0.01869201, 0.0190734, 0.01947686, 0.01989244, 0.0203262, 0.02072806, 0.02117075, 0.02158456, 0.0220196, 0.02243525, 0.02287579, 0.02330873, 0.02378892, 0.02426045, 0.0246909, 0.02518145, 0.02562587, 0.02604043, 0.02652886, 0.02704774, 0.02750332, 0.02800563, 0.02849416, 0.02894206, 0.02944305, 0.02995148, 0.03042465, 0.03097499, 0.0314764, 0.03193734, 0.03246506, 0.03299463, 0.03346533, 0.03405207, 0.03455993, 0.03504228, 0.03557384, 0.03617296, 0.03664324, 0.03723682, 0.03777067, 0.03826703, 0.03881714, 0.03938949, 0.03992124, 0.04048775, 0.04099398, 0.04155824, 0.04215843, 0.04268208, 0.04325657, 0.04377992, 0.04438168, 0.04493168, 0.0454958, 0.04603447, 0.04664449, 0.04716776, 0.04778788, 0.04831824, 0.0488735, 0.04946389, 0.05004956, 0.05062853, 0.05117396, 0.0517498, 0.05236439, 0.05293813, 0.05344883, 0.05406169, 0.05459824, 0.05525426, 0.05581771, 0.05637296, 0.05699215, 0.05751481, 0.05812418, 0.05867511, 0.05926205, 0.05986877, 0.06037511, 0.06096654, 0.06162034, 0.06212246, 0.06275838, 0.06327152, 0.06390658, 0.0644686, 0.06500423, 0.0655841, 0.06610327, 0.06675073, 0.06729298, 0.06788924, 0.0684595, 0.06902144, 0.06959623, 0.07012606, 0.07072264, 0.07122095, 0.07181285, 0.07228482, 0.07303198, 0.07354112, 0.07407519, 0.07463947, 0.07512829, 0.07575607, 0.07625307, 0.0767864, 0.07737168, 0.07791408, 0.07849061, 0.0789163, 0.07961394, 0.08005162, 0.08057328, 0.08113932, 0.08165574, 0.08221139, 0.08265689, 0.0832696, 0.08371938, 0.0842802, 0.08475816, 0.08531532, 0.08582764, 0.08632393, 0.08685683, 0.08729238, 0.08793576, 0.08834902, 0.08887366, 0.08936704, 0.08985994, 0.09032548, 0.09083868, 0.09135616, 0.09186821, 0.09226506, 0.09276164, 0.09324799, 0.09374782, 0.0942153, 0.09472576, 0.09518586, 0.09557018, 0.09617676, 0.09655638, 0.09705676, 0.09742687, 0.09789264, 0.09841924, 0.09883042, 0.09927472, 0.09979721, 0.1001492, 0.100627, 0.1010435, 0.1015334, 0.101913, 0.1023782, 0.1027738, 0.1032345, 0.1035449, 0.1040717, 0.1044323, 0.1048842, 0.1053138, 0.1056763, 0.1060852, 0.1064953, 0.1069825, 0.1073439, 0.1076406, 0.1081476, 0.10847, 0.1106424, 0.1092464, 0.1131218, 0.110087, 0.1138517, 0.1107798, 0.1147109, 0.1115775, 0.1153843, 0.1123115, 0.1161422, 0.1130158, 0.1169141, 0.1137571, 0.1176724, 0.114386, 0.1184135, 0.1152226, 0.1190563, 0.1158426, 0.1197785, 0.1165387, 0.1204201, 0.1173539, 0.1211776, 0.1176284, 0.1216924, 0.118484, 0.1224764, 0.1191357, 0.123083, 0.1234123, 0.120076, 0.1239988, 0.1205999, 0.1247054, 0.1213099, 0.1252645, 0.1218335, 0.1221377, 0.1262346, 0.1227622, 0.1268683, 0.1233242, 0.127422, 0.1237353, 0.1242564, 0.12834, 0.1247159, 0.1287923, 0.1252558, 0.1293199, 0.1296852, 0.1260079, 0.1302117, 0.1266001, 0.1306862, 0.1270103, 0.1273707, 0.1315322, 0.1278145, 0.1320164, 0.1322892, 0.1285254, 0.1327417, 0.1290498, 0.1293069, 0.1335384, 0.1296046, 0.1340581, 0.1302318, 0.1304315, 0.1346932, 0.1308097, 0.1311128, 0.135315, 0.1314942, 0.1358077, 0.1320151, 0.132141, 0.1363499, 0.1325009, 0.1327795, 0.1371205, 0.1330974, 0.1333367, 0.1377104, 0.133643, 0.1339286, 0.1382977, 0.1344676, 0.1387478, 0.1387793, 0.1348378, 0.1393521, 0.1395185, 0.1353271, 0.1398398, 0.1400374, 0.1358591, 0.140386, 0.1405752, 0.1364574, 0.1409041, 0.1410132, 0.1368252, 0.1370354, 0.1416039, 0.1373761, 0.1375143, 0.1420655, 0.1378001, 0.1378464, 0.1424957, 0.1382569, 0.1384012, 0.1429881, 0.1385256, 0.1387752, 0.1434764, 0.1436309, 0.1391633, 0.1393159, 0.1438895, 0.1396415, 0.1397422, 
1e-100, 2.859336e-58, 6.565779e-58, 1.112035e-57, 1.665879e-57, 2.331619e-57, 3.111468e-57, 4.024947e-57, 5.076125e-57, 6.27145e-57, 7.801827e-57, 9.633722e-57, 1.16422e-56, 1.430727e-56, 1.725843e-56, 2.041689e-56, 2.428395e-56, 2.838647e-56, 3.271535e-56, 3.777103e-56, 4.295676e-56, 4.849254e-56, 5.458473e-56, 6.105968e-56, 7.952365e-12, 2.563691e-10, 1.048133e-09, 2.963313e-09, 8.197567e-09, 1.652555e-08, 2.885057e-08, 5.082715e-08, 8.222571e-08, 1.245617e-07, 1.839797e-07, 2.655391e-07, 3.694208e-07, 4.994841e-07, 6.749773e-07, 8.888254e-07, 1.144064e-06, 1.454925e-06, 1.838327e-06, 2.291434e-06, 2.824449e-06, 3.481381e-06, 4.248833e-06, 5.148894e-06, 6.205685e-06, 7.504153e-06, 9.051346e-06, 1.082858e-05, 1.285674e-05, 1.538423e-05, 1.839077e-05, 2.164605e-05, 2.523478e-05, 2.983138e-05, 3.469036e-05, 4.010062e-05, 4.566614e-05, 5.246069e-05, 5.982784e-05, 6.764627e-05, 7.573587e-05, 8.522569e-05, 9.577424e-05, 0.0001065465, 0.0001182156, 0.000131376, 0.0001448501, 0.0001604398, 0.0001763045, 0.0001930602, 0.0002118739, 0.0002320245, 0.0002535116, 0.0002757718, 0.000299814, 0.000326071, 0.0003549206, 0.0003832736, 0.0004148466, 0.0004485057, 0.0004834901, 0.0005205645, 0.000559077, 0.0005996061, 0.0006456395, 0.0006940726, 0.0007427838, 0.0007932294, 0.0008487588, 0.0009066717, 0.0009666292, 0.001026958, 0.001090484, 0.001159434, 0.001231456, 0.001308137, 0.001379314, 0.001459016, 0.001538708, 0.001628154, 0.001714981, 0.001810171, 0.001903939, 0.002002155, 0.002112132, 0.002219145, 0.002333165, 0.002445234, 0.002564218, 0.002693626, 0.002819245, 0.002957132, 0.003085353, 0.003228952, 0.003378025, 0.003529235, 0.00369201, 0.003845091, 0.004011802, 0.004172473, 0.004356735, 0.004549705, 0.004727989, 0.004910066, 0.005104756, 0.005309071, 0.005511901, 0.005725892, 0.00594113, 0.006163094, 0.00638531, 0.006616404, 0.006869761, 0.007111048, 0.007359132, 0.007619551, 0.007874344, 0.008148102, 0.008427977, 0.008705347, 0.008991565, 0.009301836, 0.009591315, 0.009888667, 0.0101991, 0.01052583, 0.01086836, 0.0111875, 0.01151157, 0.01183198, 0.01221529, 0.0125694, 0.01291184, 0.0132787, 0.01365881, 0.01404163, 0.01442574, 0.01483409, 0.01523523, 0.01564674, 0.01607328, 0.01648976, 0.01689651, 0.01732918, 0.0177785, 0.01823169, 0.01868, 0.01914949, 0.01960823, 0.02005364, 0.02052261, 0.02103945, 0.02153118, 0.02202416, 0.02253692, 0.02304973, 0.02353913, 0.02408669, 0.02460619, 0.025149, 0.02570948, 0.0262671, 0.02677913, 0.02734912, 0.02792381, 0.02850586, 0.02907295, 0.02968412, 0.03028585, 0.03085916, 0.03148163, 0.03206245, 0.0326698, 0.03330131, 0.03394129, 0.03456148, 0.03517373, 0.03583646, 0.03647327, 0.03709937, 0.03777399, 0.03840443, 0.03907831, 0.03974362, 0.04044289, 0.04108567, 0.04182922, 0.04250047, 0.04314932, 0.04383058, 0.04459298, 0.04526133, 0.04600827, 0.04670908, 0.0474077, 0.04809366, 0.04886563, 0.04955405, 0.05033652, 0.05107214, 0.05181841, 0.05254209, 0.05327872, 0.05397473, 0.05477025, 0.05553517, 0.05629657, 0.05703806, 0.05779888, 0.05852367, 0.05932552, 0.06009608, 0.06084821, 0.06163413, 0.0623902, 0.06318512, 0.06399485, 0.06471052, 0.06546525, 0.06631581, 0.06699081, 0.06789567, 0.06865865, 0.06940529, 0.0702106, 0.07108179, 0.07181628, 0.07259739, 0.07332917, 0.0741647, 0.07496958, 0.07574952, 0.076551, 0.07729263, 0.07815818, 0.0789114, 0.0796715, 0.08046049, 0.08135378, 0.08207497, 0.08290438, 0.08363714, 0.08444209, 0.08524792, 0.08607514, 0.08685502, 0.08755036, 0.08845134, 0.08918766, 0.08996574, 0.09077224, 0.09158094, 0.09232554, 0.09318328, 0.09387191, 0.09467282, 0.09541753, 0.09616082, 0.09702986, 0.09768821, 0.09859079, 0.09925777, 0.1000625, 0.1008301, 0.1016524, 0.1023644, 0.1030671, 0.1038805, 0.1046217, 0.105414, 0.106141, 0.1069067, 0.1076343, 0.1083692, 0.1090331, 0.109911, 0.1107012, 0.1113442, 0.112092, 0.1128343, 0.113526, 0.1142644, 0.1148608, 0.1157754, 0.1164058, 0.1170996, 0.1178214, 0.1185698, 0.1192414, 0.1198708, 0.1206524, 0.1212623, 0.122036, 0.1227405, 0.1233839, 0.1240627, 0.124732, 0.1254511, 0.1260133, 0.1268059, 0.1274201, 0.1280989, 0.1287685, 0.1294383, 0.1301177, 0.1306431, 0.1314418, 0.1320321, 0.1326112, 0.1333264, 0.1339264, 0.1345527, 0.1351846, 0.1358496, 0.1364089, 0.1369623, 0.1376989, 0.1382481, 0.1388787, 0.1394392, 0.1400772, 0.1406958, 0.14119, 0.1419131, 0.1424092, 0.1429986, 0.1435976, 0.1441919, 0.1447714, 0.1453398, 0.1458247, 0.1464879, 0.1469764, 0.1475644, 0.1481953, 0.1486321, 0.1491634, 0.1497025, 0.1502968, 0.1507407, 0.1512607, 0.1519874, 0.1523611, 0.1528458, 0.1534565, 0.1539081, 0.1544384, 0.1549478, 0.1554388, 0.1559463, 0.1587625, 0.1570497, 0.1623087, 0.1579052, 0.1632095, 0.1589029, 0.1642544, 0.1597828, 0.1652409, 0.1607774, 0.1662044, 0.1615009, 0.1671749, 0.1625695, 0.1680708, 0.1634557, 0.1689769, 0.164346, 0.1697988, 0.165208, 0.170841, 0.1660921, 0.1715348, 0.1671897, 0.1724482, 0.1674366, 0.173274, 0.1685126, 0.1741637, 0.169221, 0.1749153, 0.1754378, 0.1704708, 0.1761324, 0.1712653, 0.1769074, 0.1719828, 0.1777146, 0.1727612, 0.1731885, 0.1788216, 0.1738313, 0.1796645, 0.1745354, 0.1803505, 0.1752939, 0.1756584, 0.1813924, 0.1763296, 0.1821325, 0.1770948, 0.1828469, 0.1830252, 0.1779309, 0.1838575, 0.1786441, 0.1844843, 0.1792718, 0.179596, 0.1853572, 0.1801246, 0.1861472, 0.1865077, 0.1810598, 0.1869394, 0.1816523, 0.1819666, 0.1880401, 0.1824997, 0.1885685, 0.1831028, 0.1832569, 0.189349, 0.1839225, 0.1842771, 0.1902171, 0.1846576, 0.1907281, 0.1852515, 0.1855619, 0.1916886, 0.1859676, 0.1862529, 0.1922873, 0.1866975, 0.1870464, 0.1932301, 0.1874169, 0.1876188, 0.1938277, 0.1882832, 0.1944397, 0.1946789, 0.1888478, 0.1951178, 0.1952401, 0.189471, 0.1958135, 0.1960506, 0.1901228, 0.1963933, 0.1966447, 0.1906724, 0.1970798, 0.1973515, 0.1913143, 0.1915257, 0.1979217, 0.1917555, 0.1921248, 0.1985232, 0.1924586, 0.1925209, 0.1990436, 0.1928559, 0.1931085, 0.1997102, 0.1934889, 0.1936827, 0.2001826, 0.200215, 0.1941849, 0.1943629, 0.2007122, 0.1947612, 0.1946787, 
1e-100, 2.766916e-57, 5.699466e-57, 9.091814e-57, 1.303998e-56, 1.715639e-56, 2.198163e-56, 2.714058e-56, 3.254761e-56, 3.878946e-56, 4.515771e-56, 7.640273e-11, 5.556417e-10, 1.661489e-09, 5.606746e-09, 1.252856e-08, 2.313425e-08, 4.381313e-08, 7.325682e-08, 1.143579e-07, 1.759107e-07, 2.574442e-07, 3.627392e-07, 5.054077e-07, 6.848435e-07, 9.045656e-07, 1.184035e-06, 1.525297e-06, 1.932068e-06, 2.419665e-06, 2.999351e-06, 3.672926e-06, 4.462394e-06, 5.463264e-06, 6.611318e-06, 7.941208e-06, 9.521646e-06, 1.146083e-05, 1.362617e-05, 1.604763e-05, 1.908635e-05, 2.246553e-05, 2.621772e-05, 3.036695e-05, 3.557632e-05, 4.105964e-05, 4.714179e-05, 5.426594e-05, 6.258207e-05, 7.148102e-05, 8.04911e-05, 9.169099e-05, 0.0001044779, 0.0001176944, 0.0001315355, 0.0001483172, 0.0001658403, 0.0001848321, 0.0002041318, 0.0002268726, 0.0002498782, 0.0002740921, 0.0002993749, 0.0003277546, 0.0003583285, 0.000390135, 0.0004219795, 0.0004597694, 0.0004968084, 0.000536469, 0.0005779202, 0.0006243558, 0.0006717336, 0.0007220068, 0.0007737603, 0.0008271524, 0.0008890095, 0.0009503708, 0.001013447, 0.00108191, 0.001153055, 0.001228051, 0.001305944, 0.001388273, 0.001471991, 0.001564179, 0.001660297, 0.001760162, 0.001861298, 0.001965641, 0.002081723, 0.002196569, 0.002317453, 0.002439433, 0.002571973, 0.00270678, 0.002852465, 0.002993215, 0.003145866, 0.003299012, 0.003461674, 0.003629086, 0.003799738, 0.003974627, 0.004152277, 0.004346455, 0.004542384, 0.004741059, 0.004951465, 0.0051532, 0.005384273, 0.005606069, 0.00584743, 0.006083087, 0.006319923, 0.00658395, 0.006840765, 0.007115746, 0.007386249, 0.007675547, 0.00796945, 0.00825687, 0.008571355, 0.008877714, 0.009205783, 0.009527707, 0.009850566, 0.01020974, 0.01055285, 0.01090785, 0.01126283, 0.01165427, 0.01202439, 0.01243335, 0.01281684, 0.01324361, 0.0136497, 0.0140662, 0.01448455, 0.01496408, 0.01541291, 0.01586338, 0.01631715, 0.01680108, 0.01729227, 0.01779733, 0.01827873, 0.01878224, 0.01931197, 0.01983315, 0.02035087, 0.02089014, 0.02145578, 0.02202332, 0.02259267, 0.02315269, 0.02373278, 0.02435547, 0.0249631, 0.02558028, 0.02620452, 0.0268149, 0.02748674, 0.02812796, 0.02876588, 0.02944983, 0.03010659, 0.03078799, 0.03146188, 0.03217375, 0.03286183, 0.03356229, 0.03429131, 0.03506133, 0.03576304, 0.03651515, 0.0372756, 0.03802864, 0.03880215, 0.03962122, 0.04040039, 0.04118, 0.04199103, 0.04278636, 0.04355182, 0.04438107, 0.04524671, 0.04608783, 0.04693637, 0.04779182, 0.04864915, 0.04948044, 0.050393, 0.05128156, 0.05212674, 0.05299619, 0.0539537, 0.05480696, 0.05574908, 0.05667673, 0.05755147, 0.05848948, 0.05941409, 0.06033945, 0.06127531, 0.06222831, 0.06320322, 0.06413566, 0.06514191, 0.06611261, 0.06705469, 0.06801154, 0.06903687, 0.0699715, 0.07100106, 0.07200727, 0.07301919, 0.07401972, 0.07497801, 0.07598533, 0.07700783, 0.07803142, 0.07901579, 0.08012196, 0.08107809, 0.08210283, 0.08313116, 0.08422128, 0.08515131, 0.08630324, 0.08728984, 0.08832171, 0.08935629, 0.09046118, 0.09146704, 0.09257173, 0.09356369, 0.09462121, 0.09573381, 0.09669537, 0.09782017, 0.0988405, 0.09988358, 0.1009963, 0.1020222, 0.1031052, 0.1041603, 0.105155, 0.1062697, 0.1073631, 0.1083985, 0.1093846, 0.1105946, 0.1115094, 0.1126373, 0.1136675, 0.1147322, 0.1158235, 0.1168312, 0.1179319, 0.1189654, 0.1200167, 0.1211293, 0.1221097, 0.1231168, 0.1242973, 0.1251403, 0.1263344, 0.1273581, 0.1284022, 0.1294753, 0.1304799, 0.1315526, 0.1325252, 0.1336102, 0.1345763, 0.1355847, 0.1366565, 0.1377277, 0.1386144, 0.139801, 0.1406543, 0.1418195, 0.1427299, 0.1437417, 0.1447166, 0.1456527, 0.1468123, 0.1477342, 0.1487284, 0.1497715, 0.1507417, 0.1517305, 0.1526871, 0.1537007, 0.1545256, 0.1556246, 0.1565448, 0.1573814, 0.1586696, 0.1594893, 0.1604233, 0.1611863, 0.1623314, 0.1631984, 0.1640921, 0.1650455, 0.166011, 0.166951, 0.1677688, 0.1688499, 0.169673, 0.1705158, 0.1714587, 0.172328, 0.1732324, 0.1741272, 0.1750257, 0.1758044, 0.1767711, 0.1775826, 0.1784763, 0.1793223, 0.1801859, 0.1810825, 0.1817389, 0.1828265, 0.1835444, 0.1844064, 0.1852655, 0.1860766, 0.1868559, 0.1875638, 0.1885095, 0.1893094, 0.1900023, 0.1909243, 0.1917064, 0.1924699, 0.1932369, 0.1940642, 0.1948876, 0.1953956, 0.1964233, 0.1970656, 0.197825, 0.1985975, 0.1993563, 0.199916, 0.2008379, 0.2015864, 0.202289, 0.2029402, 0.203703, 0.204338, 0.2050993, 0.2057894, 0.2065206, 0.2073458, 0.2078206, 0.2085153, 0.209374, 0.2098745, 0.210547, 0.2112443, 0.2119782, 0.2125282, 0.2132051, 0.213867, 0.2145588, 0.2149198, 0.2157953, 0.2163524, 0.2170215, 0.2210874, 0.2182618, 0.2258219, 0.21942, 0.227113, 0.2207671, 0.2282262, 0.2217293, 0.229509, 0.2229524, 0.2307299, 0.2241955, 0.2319499, 0.2252984, 0.2328047, 0.2265048, 0.2341867, 0.2274904, 0.2352415, 0.2286005, 0.236335, 0.2296428, 0.2374723, 0.2310333, 0.2385471, 0.2310923, 0.2395978, 0.2325979, 0.2406279, 0.2335727, 0.2416302, 0.2421358, 0.2349733, 0.2430497, 0.2361419, 0.2440496, 0.2368305, 0.2450603, 0.2378694, 0.2384109, 0.2464334, 0.2392183, 0.2474814, 0.2399944, 0.2482333, 0.2410337, 0.2415259, 0.2496228, 0.2421875, 0.2504588, 0.24308, 0.2514185, 0.2518847, 0.2443471, 0.2526241, 0.2449181, 0.2535425, 0.2458627, 0.2463106, 0.2546039, 0.2470039, 0.2553498, 0.2558257, 0.2481455, 0.2565073, 0.2489021, 0.2491912, 0.2575334, 0.2499221, 0.2583979, 0.2504957, 0.2508462, 0.2595105, 0.2514769, 0.2518641, 0.2605214, 0.252602, 0.2612575, 0.2530868, 0.25328, 0.2622742, 0.2540564, 0.2543508, 0.2631628, 0.2549794, 0.2552273, 0.2639993, 0.2558664, 0.2562418, 0.2650146, 0.2567669, 0.2653708, 0.2659674, 0.2575468, 0.2664519, 0.2666937, 0.2583338, 0.2672792, 0.2674566, 0.2590217, 0.2681089, 0.2684227, 0.2598219, 0.2686649, 0.2690932, 0.2605767, 0.2607881, 0.2698814, 0.2611742, 0.2613589, 0.2704773, 0.2617968, 0.262059, 0.2713581, 0.2624245, 0.262556, 0.2719329, 0.2631553, 0.2633762, 0.2726566, 0.272857, 0.2639197, 0.2640845, 0.2731243, 0.2645054, 0.2647127, 
1e-100, 1.049242e-09, 4.593549e-09, 1.113679e-08, 2.355887e-08, 4.347501e-08, 7.231234e-08, 1.174435e-07, 1.799782e-07, 2.65741e-07, 3.838315e-07, 5.363015e-07, 7.320564e-07, 9.817266e-07, 1.288242e-06, 1.669745e-06, 2.140307e-06, 2.704858e-06, 3.393844e-06, 4.277165e-06, 5.29738e-06, 6.482076e-06, 8.041969e-06, 9.791158e-06, 1.171226e-05, 1.421635e-05, 1.698763e-05, 1.999599e-05, 2.378213e-05, 2.797988e-05, 3.240981e-05, 3.765498e-05, 4.362634e-05, 5.007075e-05, 5.712876e-05, 6.572666e-05, 7.466851e-05, 8.426614e-05, 9.567349e-05, 0.0001081624, 0.0001214784, 0.0001354449, 0.0001526017, 0.0001701741, 0.0001891008, 0.0002103585, 0.0002336904, 0.0002578066, 0.0002844209, 0.0003149403, 0.0003458891, 0.000379187, 0.0004147707, 0.0004550139, 0.0004975493, 0.0005413843, 0.0005888238, 0.0006408291, 0.0006952915, 0.0007500293, 0.0008094259, 0.0008730766, 0.0009392388, 0.001007822, 0.001079092, 0.001156837, 0.001239267, 0.001320839, 0.001409258, 0.001503793, 0.001598707, 0.001702192, 0.001807809, 0.001918424, 0.002035448, 0.002155063, 0.002282463, 0.002412732, 0.002551785, 0.002693555, 0.00284579, 0.002997893, 0.003156648, 0.003332093, 0.003505748, 0.003677383, 0.003870071, 0.00406516, 0.004270692, 0.004487162, 0.004703622, 0.004922917, 0.00516628, 0.005401255, 0.005653203, 0.005899552, 0.006164916, 0.006431908, 0.006723214, 0.007001892, 0.007291892, 0.007587887, 0.007908583, 0.008228733, 0.008564492, 0.008912622, 0.009240439, 0.009612088, 0.009977257, 0.01036222, 0.01074752, 0.01112913, 0.0115392, 0.0119591, 0.01238901, 0.01283875, 0.01326491, 0.01373144, 0.01419207, 0.01468721, 0.01518138, 0.01566786, 0.01619488, 0.01669928, 0.01723261, 0.01775706, 0.01833466, 0.01889906, 0.01945445, 0.02004664, 0.02065343, 0.02126779, 0.02188496, 0.02250893, 0.02314708, 0.02381133, 0.02446115, 0.02516433, 0.02586703, 0.0265332, 0.02724332, 0.02795224, 0.02871705, 0.02945976, 0.03021913, 0.03094637, 0.0317419, 0.03250631, 0.03331976, 0.0341184, 0.03495324, 0.035795, 0.03661201, 0.03748391, 0.03833935, 0.03925008, 0.04013644, 0.0410216, 0.04188891, 0.0428332, 0.04379241, 0.04474155, 0.04568582, 0.04663968, 0.04759427, 0.04856179, 0.04957367, 0.05057749, 0.05154548, 0.05257737, 0.05360946, 0.05465297, 0.05569232, 0.05675923, 0.0578422, 0.05888279, 0.06001271, 0.06111246, 0.06220653, 0.06331937, 0.06445022, 0.06556442, 0.06669523, 0.06785635, 0.06905941, 0.07016331, 0.07133922, 0.07256254, 0.0737118, 0.07494458, 0.07614532, 0.07736901, 0.07855814, 0.07976184, 0.08096984, 0.08223413, 0.08351528, 0.08478813, 0.08602259, 0.08732783, 0.08849708, 0.08971597, 0.09101334, 0.09236252, 0.09365823, 0.09499793, 0.09631448, 0.09759744, 0.09885654, 0.100219, 0.1014703, 0.1028594, 0.1042058, 0.1056048, 0.1068473, 0.1082567, 0.1095914, 0.1108832, 0.1122346, 0.1137102, 0.1149586, 0.116413, 0.1177818, 0.1190538, 0.1204644, 0.1218203, 0.1231585, 0.12466, 0.126, 0.1273488, 0.1287622, 0.1301118, 0.1314831, 0.132889, 0.1342835, 0.135625, 0.1370601, 0.1383209, 0.1399351, 0.1412042, 0.1425894, 0.1439623, 0.145457, 0.1467235, 0.1482049, 0.1495263, 0.1509553, 0.1524215, 0.1536993, 0.155238, 0.1565028, 0.1578777, 0.1593044, 0.1606957, 0.1620334, 0.1634504, 0.1646762, 0.1662724, 0.1675341, 0.1689261, 0.1702945, 0.1717167, 0.1730258, 0.1744547, 0.1757297, 0.1772, 0.1784444, 0.1798347, 0.1811443, 0.182297, 0.1839211, 0.1851698, 0.1865581, 0.1878991, 0.1892538, 0.1904912, 0.1919384, 0.1931467, 0.1944438, 0.1957406, 0.1970279, 0.1983515, 0.1995702, 0.201067, 0.2022529, 0.2035298, 0.2047694, 0.2060804, 0.2073567, 0.208579, 0.2098969, 0.211058, 0.2124044, 0.2136171, 0.2148742, 0.2160918, 0.2173504, 0.2185733, 0.2193695, 0.2212867, 0.2222372, 0.223419, 0.2246374, 0.2258292, 0.2270337, 0.2280835, 0.2295012, 0.2304982, 0.2316494, 0.2328523, 0.2339904, 0.2352085, 0.2363245, 0.2374323, 0.2384517, 0.2395767, 0.240756, 0.2418453, 0.2429747, 0.2440618, 0.2452282, 0.2460694, 0.2474172, 0.2483587, 0.2495324, 0.250551, 0.2515677, 0.2526773, 0.2536561, 0.2547833, 0.2557848, 0.2568155, 0.2577149, 0.2588098, 0.2597659, 0.2608668, 0.2618966, 0.2628476, 0.2635764, 0.2648324, 0.2656594, 0.266763, 0.2676548, 0.2685392, 0.2695409, 0.2704217, 0.2714697, 0.2723952, 0.2730381, 0.2743117, 0.2750023, 0.2760138, 0.2768857, 0.2778092, 0.2787054, 0.2793856, 0.2804849, 0.2813075, 0.2821851, 0.2829374, 0.2838412, 0.2846639, 0.285394, 0.2862331, 0.2872786, 0.2880124, 0.2886366, 0.2895371, 0.2904473, 0.2912317, 0.2919934, 0.2927954, 0.2934844, 0.2939448, 0.2999168, 0.2957403, 0.3060105, 0.2972771, 0.3075291, 0.2987547, 0.3090925, 0.3000843, 0.3107224, 0.3016335, 0.3118385, 0.3030943, 0.3135495, 0.3044399, 0.3149313, 0.3057805, 0.3162997, 0.3070991, 0.3177565, 0.3085288, 0.3189981, 0.3095649, 0.3203687, 0.3114726, 0.3217374, 0.3118292, 0.3228972, 0.3133764, 0.3239582, 0.314793, 0.325538, 0.3261176, 0.3162994, 0.3272629, 0.3175393, 0.32844, 0.3187514, 0.3298136, 0.3198733, 0.3202006, 0.3315187, 0.321638, 0.3325935, 0.3225188, 0.3336436, 0.3234998, 0.3240561, 0.3353488, 0.3251411, 0.3364652, 0.3260665, 0.3371333, 0.3380903, 0.3275999, 0.3389022, 0.3285696, 0.3398791, 0.3293529, 0.3299405, 0.3413797, 0.3309612, 0.342302, 0.3426796, 0.332196, 0.3437694, 0.3331342, 0.333429, 0.3451089, 0.3343325, 0.3458823, 0.3351085, 0.3356299, 0.3474002, 0.3363027, 0.3367223, 0.3485678, 0.3374854, 0.3493973, 0.3382481, 0.3386554, 0.3506614, 0.3392936, 0.3397139, 0.3517707, 0.3406257, 0.3408211, 0.3527335, 0.3414564, 0.3418555, 0.3539707, 0.3426936, 0.3546075, 0.3550629, 0.343398, 0.3555655, 0.3561569, 0.3444878, 0.3567337, 0.3568848, 0.3452704, 0.3573431, 0.3580321, 0.3463285, 0.3586497, 0.3590176, 0.347059, 0.3472008, 0.3598408, 0.3479567, 0.3482844, 0.3606586, 0.3485174, 0.3487736, 0.3613523, 0.3494788, 0.3497115, 0.3623984, 0.3502564, 0.3503542, 0.3630032, 0.3635335, 0.3510618, 0.3513897, 0.3638118, 0.3519057, 0.3517536, 
1e-100, 1.686167e-07, 4.019961e-07, 7.113353e-07, 1.094399e-06, 1.569403e-06, 2.209328e-06, 2.959692e-06, 3.901715e-06, 5.101775e-06, 6.471436e-06, 8.193562e-06, 1.024032e-05, 1.250632e-05, 1.533263e-05, 1.858401e-05, 2.214313e-05, 2.646978e-05, 3.135644e-05, 3.659394e-05, 4.294445e-05, 4.985005e-05, 5.727393e-05, 6.602572e-05, 7.559654e-05, 8.589415e-05, 9.765122e-05, 0.0001105072, 0.0001241742, 0.0001396419, 0.0001565058, 0.0001747572, 0.0001939436, 0.0002165774, 0.0002402424, 0.0002652335, 0.0002936188, 0.0003242111, 0.000356467, 0.0003906926, 0.0004291929, 0.0004703943, 0.0005129901, 0.000559065, 0.0006095547, 0.0006640673, 0.0007199375, 0.0007824565, 0.0008462797, 0.0009162462, 0.0009875789, 0.001068458, 0.001152341, 0.001238483, 0.00133132, 0.001431694, 0.001534072, 0.001643694, 0.001756565, 0.001880482, 0.002004284, 0.002135124, 0.002266709, 0.00241254, 0.002561042, 0.002716492, 0.002872728, 0.003049781, 0.003220203, 0.003401876, 0.00358646, 0.003793217, 0.003998396, 0.004209765, 0.004428677, 0.004656371, 0.004897933, 0.005146276, 0.005398713, 0.005663904, 0.005939035, 0.006219528, 0.006517516, 0.006817869, 0.007128898, 0.00746451, 0.007801817, 0.008140511, 0.008500577, 0.008870308, 0.00924403, 0.009639519, 0.01003962, 0.01044773, 0.01087531, 0.01131631, 0.01176584, 0.01222671, 0.01269676, 0.0131841, 0.01367567, 0.01418264, 0.01470601, 0.0152465, 0.01577816, 0.01635118, 0.01690309, 0.01750397, 0.01808338, 0.01870372, 0.019337, 0.01997085, 0.02064722, 0.02128046, 0.02196956, 0.02267192, 0.02339107, 0.02413867, 0.02484614, 0.02559655, 0.0263437, 0.02716268, 0.02796055, 0.02873652, 0.02957647, 0.03039218, 0.03125177, 0.03211143, 0.03301808, 0.03389943, 0.03480516, 0.03570083, 0.03666951, 0.03763734, 0.03858563, 0.03955309, 0.04053599, 0.04156882, 0.04258774, 0.04361534, 0.0446688, 0.04573351, 0.0468083, 0.0478651, 0.04898078, 0.0501122, 0.05127281, 0.05239299, 0.05354669, 0.05472285, 0.0559311, 0.05712619, 0.05831614, 0.05957099, 0.06084205, 0.06205168, 0.06328806, 0.0645873, 0.06590501, 0.06723564, 0.06853462, 0.06985165, 0.07114417, 0.072519, 0.07385762, 0.07524836, 0.07665015, 0.07802135, 0.07944991, 0.08084402, 0.08231547, 0.08375736, 0.08520203, 0.08665002, 0.08813298, 0.08960252, 0.09112486, 0.09264597, 0.09418075, 0.09563546, 0.09717126, 0.09873086, 0.1002652, 0.1018517, 0.1034471, 0.1050701, 0.1065676, 0.1082831, 0.1098331, 0.1114307, 0.113044, 0.1146671, 0.1162824, 0.1179458, 0.1196555, 0.1213528, 0.1229582, 0.1247011, 0.1264076, 0.1280029, 0.1295649, 0.1314049, 0.1330306, 0.1348041, 0.1365611, 0.1382105, 0.1400063, 0.1417146, 0.143368, 0.1452199, 0.1469408, 0.1487467, 0.1504576, 0.1522541, 0.1539479, 0.1556998, 0.1574955, 0.1592683, 0.1610538, 0.1629059, 0.1646145, 0.1663312, 0.1681424, 0.1699138, 0.1718257, 0.1735905, 0.1753338, 0.1770978, 0.1789217, 0.1806829, 0.1824949, 0.1844022, 0.1861099, 0.187922, 0.1897859, 0.1914538, 0.1933334, 0.1951257, 0.1969604, 0.1986517, 0.2004879, 0.2022906, 0.2041725, 0.2058571, 0.2076312, 0.2094413, 0.2110941, 0.2131114, 0.214795, 0.2166725, 0.2183907, 0.220276, 0.2219007, 0.2237525, 0.2253666, 0.2272959, 0.2290005, 0.230766, 0.232455, 0.2341887, 0.2360656, 0.2377225, 0.2393806, 0.2411736, 0.2430156, 0.2446084, 0.2464005, 0.2480403, 0.2498492, 0.2515295, 0.2532118, 0.254975, 0.2566276, 0.2583975, 0.2599937, 0.2615875, 0.2633474, 0.2650241, 0.2665963, 0.268424, 0.2698963, 0.2716293, 0.2731947, 0.2748928, 0.276476, 0.2779752, 0.2797764, 0.2811873, 0.2828902, 0.2844898, 0.2860453, 0.287679, 0.2892551, 0.290745, 0.2921551, 0.2939691, 0.2954521, 0.29695, 0.2982785, 0.3003154, 0.301585, 0.3028551, 0.3047073, 0.3060237, 0.3074685, 0.3090406, 0.3104933, 0.3118808, 0.3133335, 0.3149502, 0.3162216, 0.3177015, 0.3190924, 0.3203935, 0.321889, 0.3232202, 0.3247166, 0.3258834, 0.3275177, 0.3287292, 0.3301907, 0.331494, 0.3328203, 0.3341514, 0.3353212, 0.3369149, 0.3380673, 0.3394533, 0.3407027, 0.3419302, 0.3431905, 0.3445513, 0.3457036, 0.3470598, 0.3480842, 0.3496147, 0.3508107, 0.3520246, 0.3532422, 0.3543999, 0.3555921, 0.356634, 0.3580051, 0.3590602, 0.3602361, 0.3613525, 0.3625166, 0.3636043, 0.3646415, 0.3660967, 0.3670852, 0.367994, 0.3693816, 0.3703459, 0.3715446, 0.372479, 0.3736512, 0.3746695, 0.375639, 0.3768445, 0.3778883, 0.3787983, 0.3796703, 0.3808044, 0.3818111, 0.3829804, 0.3837908, 0.3848764, 0.3858218, 0.3865241, 0.387796, 0.3886297, 0.3895959, 0.3904696, 0.3980132, 0.3922237, 0.4063935, 0.3942951, 0.4084294, 0.3959384, 0.4101328, 0.3975936, 0.4121244, 0.3994943, 0.4138513, 0.4011954, 0.4154897, 0.4026703, 0.4173998, 0.4045573, 0.4191334, 0.4061946, 0.4206168, 0.4076491, 0.4223808, 0.4093351, 0.4241292, 0.4113475, 0.4252547, 0.4118433, 0.4270216, 0.4137378, 0.4286811, 0.4153126, 0.4300216, 0.4309445, 0.4172589, 0.4325067, 0.4188525, 0.4335466, 0.4201373, 0.435369, 0.4215918, 0.4221924, 0.437437, 0.4235774, 0.4386735, 0.4246769, 0.4403525, 0.4260229, 0.4265941, 0.4419942, 0.427824, 0.4434026, 0.4290457, 0.4445708, 0.4453761, 0.430673, 0.4462513, 0.4321143, 0.4477359, 0.4329797, 0.4334882, 0.4493544, 0.4346662, 0.4505577, 0.4511684, 0.4364929, 0.4522998, 0.4372787, 0.4378957, 0.4540585, 0.4388206, 0.4549625, 0.4398211, 0.4403028, 0.4564095, 0.4412642, 0.4420055, 0.4581745, 0.4427029, 0.4588006, 0.4436, 0.444286, 0.4604971, 0.4449864, 0.4452994, 0.4617641, 0.4460786, 0.4466396, 0.4631569, 0.4475217, 0.4479336, 0.4643027, 0.4488281, 0.4654841, 0.4659115, 0.4497107, 0.4666217, 0.4669855, 0.4506826, 0.4677879, 0.4682406, 0.4519413, 0.4689689, 0.4692431, 0.4529813, 0.4702598, 0.4705621, 0.4539826, 0.4543525, 0.4716072, 0.4549332, 0.4550847, 0.4725414, 0.4559626, 0.4561555, 0.4735686, 0.4565526, 0.4571887, 0.474662, 0.4576957, 0.45793, 0.4755207, 0.4758213, 0.4585142, 0.4589016, 0.4765039, 0.4597134, 0.4597104, 
1e-100, 1.569964e-06, 3.723561e-06, 6.165258e-06, 9.118671e-06, 1.261493e-05, 1.645e-05, 2.117361e-05, 2.643071e-05, 3.217975e-05, 3.922928e-05, 4.68556e-05, 5.525761e-05, 6.504675e-05, 7.556336e-05, 8.721637e-05, 0.0001006393, 0.0001147273, 0.0001304827, 0.0001482911, 0.0001669579, 0.0001876708, 0.0002108122, 0.0002348531, 0.0002613653, 0.0002908296, 0.0003220831, 0.0003551128, 0.0003934133, 0.0004329904, 0.000474375, 0.0005215246, 0.0005713106, 0.0006227467, 0.0006789715, 0.0007409635, 0.000804358, 0.0008720845, 0.0009477347, 0.001026033, 0.00110681, 0.00119647, 0.001290886, 0.001390283, 0.001494841, 0.001607846, 0.0017278, 0.001851206, 0.001978698, 0.002122102, 0.002269085, 0.002420628, 0.002582753, 0.00275831, 0.002935919, 0.003123457, 0.003317483, 0.003528264, 0.003739248, 0.003958273, 0.004191312, 0.004431241, 0.004685636, 0.004938987, 0.005209657, 0.005492569, 0.005782594, 0.006073574, 0.006387128, 0.006714081, 0.00704187, 0.007388289, 0.007741519, 0.008120284, 0.0084964, 0.008891276, 0.009300493, 0.009712124, 0.01014821, 0.01060115, 0.01106118, 0.01153448, 0.01202301, 0.01253485, 0.01304564, 0.01357218, 0.01412962, 0.01470319, 0.01527743, 0.01588978, 0.01648859, 0.01712063, 0.01776901, 0.01844997, 0.01910969, 0.01980475, 0.02048909, 0.02123307, 0.02196116, 0.02271484, 0.02349647, 0.02427474, 0.02510809, 0.02591535, 0.02676362, 0.02760533, 0.02848952, 0.0294025, 0.03030325, 0.03123056, 0.03216074, 0.0331197, 0.03412008, 0.03511937, 0.03614288, 0.0371712, 0.03822441, 0.03930142, 0.04038216, 0.04154274, 0.04262892, 0.04377927, 0.04490528, 0.04611605, 0.04733562, 0.04854275, 0.0497547, 0.05101916, 0.05227555, 0.05354762, 0.05486859, 0.05622648, 0.05754814, 0.05887007, 0.06021461, 0.06166414, 0.06309317, 0.06445235, 0.06587431, 0.06737077, 0.06882055, 0.0703258, 0.07185032, 0.07341811, 0.07494113, 0.0764875, 0.07805457, 0.0796368, 0.08131247, 0.08291497, 0.08450299, 0.08618007, 0.08789824, 0.08956334, 0.09127936, 0.09295494, 0.09470714, 0.09640553, 0.09819385, 0.09995477, 0.1017157, 0.1035234, 0.105403, 0.1071575, 0.1090335, 0.1109016, 0.1127682, 0.1146871, 0.1165844, 0.1185082, 0.1203586, 0.1223407, 0.1242434, 0.126229, 0.1282056, 0.1302266, 0.1321968, 0.1341454, 0.1361687, 0.1382538, 0.140239, 0.1423264, 0.1444182, 0.1464582, 0.1485471, 0.1506891, 0.1527366, 0.1548257, 0.1569683, 0.1590454, 0.1611833, 0.1633546, 0.1655327, 0.1676342, 0.1698609, 0.1720514, 0.1742095, 0.1763772, 0.1786756, 0.1806968, 0.1828438, 0.1850946, 0.1873909, 0.1895043, 0.1919448, 0.1940992, 0.1962785, 0.1984679, 0.2008611, 0.2030438, 0.2053441, 0.2076873, 0.2098962, 0.2121125, 0.2144368, 0.2165894, 0.2189703, 0.2213419, 0.2235929, 0.2259416, 0.2281853, 0.2303573, 0.2327023, 0.2350143, 0.2372948, 0.2396087, 0.2418518, 0.2441062, 0.2464718, 0.248808, 0.2510316, 0.2534703, 0.2556719, 0.2579026, 0.26017, 0.2625624, 0.2646872, 0.2671767, 0.2694169, 0.271629, 0.2739934, 0.2762836, 0.2784598, 0.2808088, 0.2829097, 0.2852937, 0.2876008, 0.2898813, 0.2920586, 0.2944466, 0.296616, 0.2989181, 0.3010548, 0.3032569, 0.305586, 0.3075921, 0.310019, 0.3121685, 0.314322, 0.316617, 0.318823, 0.3209262, 0.323317, 0.3253132, 0.3275672, 0.329762, 0.3317409, 0.3339534, 0.3358809, 0.3382944, 0.340318, 0.3424698, 0.3447327, 0.3467608, 0.3487572, 0.3510234, 0.3530334, 0.3550887, 0.357132, 0.3591227, 0.3611755, 0.3631979, 0.3655259, 0.3673497, 0.3693463, 0.371363, 0.3733594, 0.3753417, 0.3772798, 0.3793465, 0.3811148, 0.3832606, 0.3851621, 0.3871831, 0.3890653, 0.39103, 0.3928428, 0.3943448, 0.3964264, 0.3989055, 0.400431, 0.4022875, 0.4042581, 0.405966, 0.4077361, 0.4098164, 0.4114656, 0.4132096, 0.4150486, 0.4167447, 0.418629, 0.4204099, 0.4220988, 0.4236266, 0.4254902, 0.4271915, 0.4290045, 0.4306162, 0.4323471, 0.4340285, 0.4353216, 0.4374082, 0.4389268, 0.440667, 0.4422093, 0.4437442, 0.4453405, 0.4470195, 0.4487121, 0.4501071, 0.4516954, 0.453379, 0.4545973, 0.4564132, 0.4579455, 0.4594734, 0.4608672, 0.46198, 0.4638976, 0.4652011, 0.4668078, 0.4682747, 0.4695312, 0.4709127, 0.4723713, 0.4739655, 0.4752688, 0.4765754, 0.4779558, 0.4794542, 0.4807736, 0.482116, 0.483509, 0.4848543, 0.4857923, 0.4875603, 0.4886905, 0.4899946, 0.4914558, 0.4925034, 0.4936569, 0.4950172, 0.4961495, 0.4975544, 0.4985759, 0.499967, 0.5010148, 0.50225, 0.5034844, 0.5046044, 0.5058448, 0.506799, 0.5078515, 0.5091208, 0.5191501, 0.5113284, 0.5301834, 0.5134825, 0.5324056, 0.5157759, 0.534648, 0.5179121, 0.5368127, 0.5198879, 0.5392951, 0.522024, 0.5411681, 0.5240006, 0.5431614, 0.5258072, 0.5453129, 0.5280648, 0.5473897, 0.5296816, 0.5492597, 0.5317419, 0.5513214, 0.5342019, 0.5531518, 0.5348062, 0.5548853, 0.5368717, 0.5569405, 0.5389857, 0.558601, 0.5595592, 0.5411876, 0.5614212, 0.5432354, 0.5631746, 0.5447322, 0.5647973, 0.545961, 0.5472713, 0.567402, 0.5487511, 0.5689903, 0.5501823, 0.5705325, 0.5515708, 0.5525084, 0.5731681, 0.5540155, 0.5744053, 0.5552931, 0.5762321, 0.5769276, 0.55733, 0.5782096, 0.558708, 0.5794727, 0.5599958, 0.5607122, 0.5819389, 0.5621432, 0.5828538, 0.5840236, 0.5640248, 0.5852112, 0.5651326, 0.5659041, 0.5873087, 0.5666941, 0.5882745, 0.5681517, 0.5686899, 0.5901869, 0.5696961, 0.5704276, 0.592256, 0.5715841, 0.5933074, 0.5725547, 0.57309, 0.5947609, 0.5739407, 0.5747084, 0.5966174, 0.575447, 0.5758831, 0.5981869, 0.5771371, 0.5777825, 0.5998051, 0.578766, 0.6008551, 0.6012258, 0.5797465, 0.6022668, 0.602899, 0.5809762, 0.6034794, 0.6040004, 0.5822989, 0.605186, 0.6056118, 0.5837092, 0.6064897, 0.6068562, 0.5847652, 0.5851989, 0.6085117, 0.5858875, 0.5860353, 0.6093354, 0.5866718, 0.5873715, 0.6107395, 0.5880783, 0.588419, 0.6118021, 0.5888006, 0.5892921, 0.6128735, 0.6133815, 0.5901477, 0.5904089, 0.6137777, 0.5913215, 0.5915052, 
1e-100, 6.883165e-06, 1.460336e-05, 2.39298e-05, 3.410976e-05, 4.539251e-05, 5.845668e-05, 7.23953e-05, 8.824286e-05, 0.0001057902, 0.0001245517, 0.0001459588, 0.0001691307, 0.000193744, 0.0002218584, 0.0002517776, 0.0002835257, 0.0003195235, 0.0003578801, 0.0003983078, 0.0004439971, 0.0004923118, 0.0005429426, 0.0005999113, 0.0006602062, 0.0007226691, 0.0007929798, 0.0008665796, 0.0009431894, 0.00102866, 0.001117399, 0.001211928, 0.001312356, 0.00142121, 0.001534674, 0.001653679, 0.001782334, 0.00192069, 0.002060483, 0.002211149, 0.002373881, 0.002539047, 0.002716084, 0.002903175, 0.003102943, 0.003307701, 0.00351709, 0.003749988, 0.003992381, 0.004236429, 0.004495437, 0.004774835, 0.005063012, 0.005350635, 0.005672262, 0.005997973, 0.006334073, 0.006681886, 0.007053249, 0.007431666, 0.007818463, 0.008224054, 0.008635857, 0.009077897, 0.009522187, 0.009984658, 0.01045621, 0.01096341, 0.0114716, 0.01198794, 0.01252279, 0.01310491, 0.01367464, 0.01426412, 0.01487341, 0.01549881, 0.0161482, 0.01682091, 0.01749309, 0.01819624, 0.01892573, 0.01965849, 0.02042055, 0.02119164, 0.02201508, 0.02282726, 0.02368085, 0.02454047, 0.02539759, 0.02632308, 0.02724362, 0.02821173, 0.02915675, 0.03015865, 0.0311665, 0.0321974, 0.03324473, 0.03431589, 0.03541893, 0.03653391, 0.03770136, 0.03883338, 0.04001524, 0.04116409, 0.04245211, 0.04367797, 0.04495349, 0.04624582, 0.04751603, 0.04889097, 0.05022345, 0.05164596, 0.05301145, 0.05443888, 0.05590154, 0.05736846, 0.05887615, 0.06037098, 0.06189378, 0.06346843, 0.06502557, 0.06666087, 0.06828586, 0.06992785, 0.07155262, 0.07329097, 0.07503685, 0.07673909, 0.07845727, 0.08024623, 0.0820453, 0.08384922, 0.08567199, 0.08752791, 0.08942072, 0.09130765, 0.09320309, 0.09521637, 0.09714478, 0.09912048, 0.1011738, 0.1031477, 0.1052277, 0.1072943, 0.1093315, 0.1114739, 0.1136525, 0.115777, 0.1178745, 0.1200355, 0.1222846, 0.1245304, 0.126764, 0.1289739, 0.1311765, 0.1335155, 0.1358395, 0.1380893, 0.1404527, 0.1427963, 0.1451609, 0.147549, 0.1499235, 0.1524214, 0.154767, 0.1573208, 0.1598315, 0.1622462, 0.1646836, 0.1672248, 0.1697514, 0.1722706, 0.1748662, 0.1774141, 0.1798386, 0.1824659, 0.1850662, 0.1877011, 0.1903837, 0.1930004, 0.1957136, 0.1982472, 0.2010328, 0.2036724, 0.2063471, 0.2090374, 0.2118327, 0.2145098, 0.2172598, 0.2199858, 0.2228284, 0.2254916, 0.228343, 0.2311058, 0.2338592, 0.2366718, 0.239558, 0.242232, 0.2452557, 0.2479582, 0.2506381, 0.2534525, 0.2564117, 0.2592558, 0.2620699, 0.264837, 0.2677518, 0.2706474, 0.2735483, 0.2764922, 0.2793171, 0.2822786, 0.2851037, 0.2878576, 0.290771, 0.2937496, 0.2964872, 0.2996677, 0.3024334, 0.3052655, 0.3080394, 0.3111587, 0.3138745, 0.3169906, 0.3198752, 0.3226732, 0.3254783, 0.3286457, 0.331236, 0.3342562, 0.3372122, 0.3400696, 0.3428015, 0.3459564, 0.3486483, 0.3516872, 0.3544916, 0.3573453, 0.3602388, 0.36299, 0.3660112, 0.3688222, 0.3715198, 0.3743605, 0.37727, 0.3798979, 0.3830454, 0.3857585, 0.3885534, 0.3916403, 0.3943251, 0.3969161, 0.3999584, 0.4024812, 0.4053942, 0.4081561, 0.4108309, 0.4136289, 0.4163241, 0.4192989, 0.421869, 0.4244544, 0.4272751, 0.429994, 0.4324616, 0.4355457, 0.4378598, 0.4406696, 0.4432956, 0.4460389, 0.4485552, 0.4511255, 0.4539407, 0.4563823, 0.4590099, 0.4617173, 0.4641546, 0.4666545, 0.4694418, 0.4718189, 0.4743247, 0.4767871, 0.4794584, 0.4819504, 0.4842577, 0.4869905, 0.4891177, 0.4917467, 0.4939979, 0.496551, 0.4989185, 0.501397, 0.5037854, 0.5058433, 0.5086935, 0.5108344, 0.513196, 0.5154825, 0.5175256, 0.5206827, 0.5224042, 0.5249767, 0.5270518, 0.5291997, 0.5314404, 0.5337068, 0.5359955, 0.5381623, 0.5404936, 0.5422903, 0.5446699, 0.546757, 0.5489171, 0.5509353, 0.5532031, 0.5551532, 0.5569913, 0.559401, 0.5613448, 0.5634505, 0.5654348, 0.5673523, 0.5693566, 0.5711645, 0.5734615, 0.5752889, 0.5771704, 0.57928, 0.5810066, 0.5829516, 0.584982, 0.5869085, 0.5884115, 0.5902387, 0.5924336, 0.5941279, 0.5959808, 0.5977731, 0.5996342, 0.6012393, 0.6029322, 0.605104, 0.6065084, 0.6081927, 0.6099666, 0.611448, 0.6132896, 0.614988, 0.6164383, 0.6184566, 0.6195339, 0.6217156, 0.6231602, 0.6249016, 0.6263943, 0.6278128, 0.6294018, 0.6308842, 0.6327907, 0.6341927, 0.6354014, 0.6369194, 0.6382451, 0.639979, 0.6415014, 0.6430106, 0.6443927, 0.6456488, 0.6470122, 0.6488431, 0.6499402, 0.6513297, 0.6525975, 0.653804, 0.6671695, 0.6565727, 0.6821162, 0.6593702, 0.6844265, 0.6619665, 0.6874214, 0.6644809, 0.6900598, 0.6669939, 0.6927804, 0.6691514, 0.6951228, 0.6719616, 0.697712, 0.6739038, 0.7003282, 0.6761625, 0.7027127, 0.6786304, 0.7050399, 0.6811, 0.7072261, 0.6838219, 0.7099813, 0.6843904, 0.7118562, 0.6871364, 0.7138179, 0.6891878, 0.7161926, 0.7173581, 0.6924095, 0.7193313, 0.6942412, 0.7217469, 0.6963285, 0.7235331, 0.6982091, 0.6991483, 0.7265223, 0.7009942, 0.7286922, 0.7029234, 0.7305062, 0.7041425, 0.7055339, 0.7335106, 0.7074191, 0.7352654, 0.7089002, 0.7369567, 0.737709, 0.7111115, 0.739586, 0.7129013, 0.7409907, 0.7142259, 0.7148577, 0.7438912, 0.7169254, 0.7454811, 0.7463736, 0.7190623, 0.7473259, 0.7206258, 0.7211797, 0.7503593, 0.7223667, 0.7515702, 0.7235973, 0.7244374, 0.7539715, 0.7259556, 0.726811, 0.7559951, 0.7274776, 0.7575687, 0.7290776, 0.7298906, 0.7593461, 0.7307787, 0.7312412, 0.7612516, 0.7325391, 0.7332835, 0.7635206, 0.7344238, 0.7346182, 0.7653322, 0.7363781, 0.7663168, 0.7670353, 0.737511, 0.7678769, 0.7684166, 0.7391289, 0.7697686, 0.7705646, 0.7405663, 0.7710961, 0.7723359, 0.7423173, 0.7730498, 0.7736415, 0.7433544, 0.7437028, 0.7749472, 0.7447414, 0.7451324, 0.7765798, 0.7460145, 0.7460799, 0.7779386, 0.7473742, 0.7475574, 0.7793592, 0.7484353, 0.7486975, 0.7804789, 0.7809972, 0.7498138, 0.7501734, 0.7822341, 0.7513001, 0.7508597, 
1e-100, 1.796738e-05, 3.881252e-05, 6.129167e-05, 8.672122e-05, 0.000114418, 0.00014399, 0.0001777811, 0.0002136932, 0.0002527489, 0.0002960722, 0.0003421172, 0.0003925303, 0.0004470807, 0.0005048704, 0.0005688444, 0.0006367119, 0.0007085643, 0.0007877977, 0.0008723978, 0.0009604428, 0.001057195, 0.001159704, 0.001267556, 0.001384707, 0.001508741, 0.001638244, 0.00177658, 0.00192759, 0.002082554, 0.00224669, 0.002422533, 0.002607706, 0.002800478, 0.003007643, 0.003224128, 0.003452079, 0.003691576, 0.003944085, 0.004206989, 0.004482696, 0.004773532, 0.005085464, 0.005396955, 0.00572908, 0.006082398, 0.006448024, 0.00682206, 0.007219656, 0.007632273, 0.008062904, 0.008502641, 0.008980931, 0.009465222, 0.009968246, 0.0104873, 0.01103222, 0.01159803, 0.0121684, 0.01275897, 0.0133865, 0.01401753, 0.01467933, 0.0153464, 0.01605248, 0.01676822, 0.0175128, 0.01825823, 0.01904094, 0.01985532, 0.02067614, 0.02151261, 0.02239032, 0.02329037, 0.02420528, 0.02515186, 0.02611531, 0.02710996, 0.028115, 0.02915817, 0.03022676, 0.0313026, 0.03243306, 0.03358322, 0.0347434, 0.03593365, 0.03717035, 0.03843054, 0.03971082, 0.04101117, 0.04232996, 0.04370733, 0.0450852, 0.04650592, 0.04792234, 0.0493902, 0.05086078, 0.05241229, 0.05391698, 0.05550674, 0.05710285, 0.05873511, 0.06040025, 0.06206662, 0.06378247, 0.06547694, 0.06725339, 0.06907417, 0.07088023, 0.0727293, 0.07456844, 0.07648902, 0.07842239, 0.08039019, 0.08238392, 0.08432642, 0.0863997, 0.08842297, 0.09057499, 0.09268789, 0.09484183, 0.09702466, 0.09915137, 0.1014272, 0.1036719, 0.1059753, 0.1082635, 0.1105287, 0.1129611, 0.1153316, 0.1176745, 0.1200816, 0.1226169, 0.1250273, 0.127568, 0.130078, 0.1326921, 0.135247, 0.1378516, 0.1403931, 0.1431765, 0.1458554, 0.1485034, 0.1512233, 0.15399, 0.1567574, 0.1595743, 0.162322, 0.165122, 0.167985, 0.1709457, 0.1737346, 0.1766454, 0.1796264, 0.1826064, 0.1855725, 0.1885133, 0.1915792, 0.1945908, 0.1977032, 0.2007849, 0.2039076, 0.2069472, 0.2101253, 0.2133026, 0.216374, 0.2197134, 0.2227909, 0.2260356, 0.2291994, 0.2325452, 0.2357333, 0.2389749, 0.2422552, 0.2456221, 0.2489327, 0.2523114, 0.2556762, 0.259006, 0.2623128, 0.265838, 0.2691986, 0.2726274, 0.2761128, 0.2794832, 0.2828262, 0.2862926, 0.2898452, 0.2933681, 0.2968511, 0.3004505, 0.3038973, 0.3073369, 0.311043, 0.3144468, 0.3179007, 0.321493, 0.3252421, 0.3288132, 0.3323074, 0.3356257, 0.3392604, 0.3427664, 0.3465735, 0.3500675, 0.3537131, 0.3573146, 0.3610498, 0.3644466, 0.3682504, 0.3718223, 0.3753639, 0.3789649, 0.382727, 0.3862091, 0.3899309, 0.3935544, 0.3972219, 0.4006642, 0.4045326, 0.4080093, 0.4117527, 0.4153014, 0.4189719, 0.4226896, 0.4260057, 0.4297354, 0.4334583, 0.4370356, 0.4406267, 0.4443557, 0.4476432, 0.4512239, 0.4549702, 0.4585268, 0.462056, 0.4657627, 0.4692174, 0.4728283, 0.4763549, 0.4799521, 0.4832455, 0.4871232, 0.4903351, 0.4940727, 0.4975018, 0.501038, 0.5045454, 0.5079328, 0.5117564, 0.5150286, 0.5182406, 0.5217577, 0.5252868, 0.5284275, 0.5321786, 0.535433, 0.5387159, 0.5422413, 0.5456538, 0.5488198, 0.5523525, 0.5555714, 0.5589526, 0.5621703, 0.5656408, 0.568669, 0.5718459, 0.5753608, 0.5784948, 0.581702, 0.5849872, 0.5882068, 0.5914383, 0.5946072, 0.5978109, 0.6009122, 0.6040284, 0.6070538, 0.6100402, 0.6132213, 0.6165592, 0.6194706, 0.6223922, 0.6253618, 0.6285278, 0.6314331, 0.6344423, 0.6374897, 0.6399875, 0.6433014, 0.6462627, 0.6492334, 0.6521051, 0.6549378, 0.6575469, 0.6600859, 0.6634638, 0.6656389, 0.6695841, 0.6719349, 0.6746499, 0.6772601, 0.6799812, 0.6829298, 0.685304, 0.6879898, 0.6906062, 0.6933815, 0.6959332, 0.6987415, 0.7012359, 0.7033502, 0.7062055, 0.7088938, 0.7112799, 0.7139253, 0.7162571, 0.7187162, 0.7208124, 0.7238015, 0.7258857, 0.7285658, 0.7307058, 0.7331039, 0.7353232, 0.7377725, 0.7403475, 0.7426138, 0.7445216, 0.7470243, 0.7494528, 0.7513342, 0.7538231, 0.7559523, 0.7581705, 0.7600022, 0.7625155, 0.7645083, 0.7667498, 0.7688204, 0.7710534, 0.7727279, 0.7749167, 0.7771766, 0.7789863, 0.7808038, 0.7830833, 0.7847504, 0.7872095, 0.7890377, 0.7909989, 0.7928203, 0.7943859, 0.7968218, 0.7983508, 0.8003021, 0.8021097, 0.8040708, 0.8055725, 0.8075419, 0.8093447, 0.8110429, 0.8125479, 0.8145501, 0.816624, 0.818082, 0.8195421, 0.8215662, 0.8231509, 0.8244822, 0.8260359, 0.828142, 0.8294122, 0.8466512, 0.8326167, 0.8659052, 0.835796, 0.8692317, 0.8389683, 0.8726711, 0.8414302, 0.8769165, 0.8449111, 0.8795624, 0.8475412, 0.8826756, 0.8504131, 0.8852954, 0.8532398, 0.8878075, 0.8562005, 0.8901653, 0.8586707, 0.8935292, 0.8616012, 0.8961992, 0.8650224, 0.8990054, 0.8652492, 0.9013809, 0.8690895, 0.9040461, 0.8711064, 0.9066556, 0.9088599, 0.8749972, 0.9105753, 0.8774322, 0.9130275, 0.8794652, 0.9153219, 0.8818299, 0.8830334, 0.9199672, 0.8850193, 0.9212962, 0.8868936, 0.92371, 0.8893144, 0.8904631, 0.9271095, 0.8917515, 0.9295859, 0.8943258, 0.9311797, 0.933189, 0.8970747, 0.9342903, 0.898728, 0.9363187, 0.9009545, 0.9018279, 0.9392023, 0.903512, 0.9415791, 0.9432462, 0.9062413, 0.9441941, 0.9078986, 0.9085255, 0.9469636, 0.9103353, 0.9487029, 0.9117727, 0.9123678, 0.9513304, 0.9145273, 0.9155701, 0.9538864, 0.9165422, 0.9555388, 0.9178211, 0.9186457, 0.9579108, 0.9202243, 0.9206947, 0.9598827, 0.9223328, 0.9229728, 0.9628014, 0.9242901, 0.9249253, 0.9647124, 0.9262927, 0.9660921, 0.9669246, 0.9282416, 0.9680195, 0.9686135, 0.9294322, 0.970166, 0.9711771, 0.9315734, 0.972274, 0.9729366, 0.9331901, 0.9736098, 0.9746171, 0.9347622, 0.9351706, 0.9762433, 0.9358194, 0.9365574, 0.9785055, 0.9381543, 0.9381234, 0.9801331, 0.9391938, 0.9393624, 0.9813999, 0.9404911, 0.9411927, 0.9829357, 0.9833099, 0.9418941, 0.9424477, 0.9848088, 0.9440587, 0.9438875, 
1e-100, 3.917855e-05, 8.195039e-05, 0.0001304505, 0.0001818938, 0.0002382944, 0.0002998683, 0.0003651591, 0.0004376259, 0.0005149162, 0.000596854, 0.0006880643, 0.0007838169, 0.0008858299, 0.0009975146, 0.001115129, 0.001240392, 0.001375658, 0.001519465, 0.001670382, 0.001834326, 0.002004823, 0.002186284, 0.002380832, 0.002585265, 0.002797079, 0.003030369, 0.003270392, 0.003519949, 0.003791764, 0.004074135, 0.004364621, 0.00467814, 0.005006704, 0.005346542, 0.005703495, 0.00608419, 0.006470327, 0.006880469, 0.007312287, 0.007760688, 0.008226038, 0.008702825, 0.009219495, 0.009748609, 0.01029266, 0.01086685, 0.01146569, 0.01208829, 0.01272086, 0.01338404, 0.01407504, 0.01478803, 0.0155221, 0.01629484, 0.0170845, 0.01789997, 0.01873444, 0.0196122, 0.0205027, 0.02141777, 0.02236082, 0.02332618, 0.024347, 0.02536393, 0.02641017, 0.02748292, 0.02861178, 0.02974619, 0.030894, 0.03209315, 0.03333163, 0.03458536, 0.03587659, 0.03717361, 0.03853237, 0.03992462, 0.04133763, 0.04275259, 0.04421446, 0.04576488, 0.04728831, 0.04885323, 0.05045873, 0.05208534, 0.05379287, 0.05548616, 0.05720588, 0.05895991, 0.06080429, 0.06265505, 0.0645286, 0.06644447, 0.06837945, 0.07037349, 0.07241016, 0.07444229, 0.07647428, 0.07859007, 0.0807498, 0.08292762, 0.08513077, 0.08734533, 0.08963543, 0.09194065, 0.09428111, 0.09669344, 0.09903316, 0.1015181, 0.1039658, 0.1065326, 0.1090767, 0.1116182, 0.1142277, 0.1168618, 0.1195477, 0.1222863, 0.1250165, 0.1277825, 0.1305837, 0.133441, 0.1363609, 0.1392266, 0.1422044, 0.1451153, 0.1481147, 0.1511226, 0.1542357, 0.1573233, 0.1603946, 0.1635669, 0.1667571, 0.169984, 0.1732147, 0.1764718, 0.1798204, 0.1831149, 0.1864295, 0.1899091, 0.193361, 0.1967163, 0.2002281, 0.2036459, 0.20726, 0.2108296, 0.2143613, 0.2179108, 0.2215618, 0.2251595, 0.2287998, 0.2325876, 0.2363323, 0.2401247, 0.2438397, 0.2476683, 0.2514446, 0.2553236, 0.2592742, 0.2631328, 0.266937, 0.2709463, 0.2749953, 0.278955, 0.2829075, 0.2869715, 0.2909884, 0.2950562, 0.2991018, 0.3031122, 0.3071798, 0.3114239, 0.315508, 0.3196656, 0.3238892, 0.3279989, 0.3323119, 0.3363861, 0.3408504, 0.3450346, 0.3492468, 0.3534915, 0.3579279, 0.3622104, 0.3665217, 0.3708506, 0.3753463, 0.3795231, 0.384003, 0.3884335, 0.3927195, 0.3972582, 0.4016092, 0.4060252, 0.4103764, 0.4148278, 0.4193562, 0.4237659, 0.4283168, 0.4330375, 0.437177, 0.4413803, 0.445931, 0.4502497, 0.4548777, 0.4594663, 0.4640157, 0.4685622, 0.4730787, 0.4775092, 0.4819771, 0.4863893, 0.490811, 0.4954749, 0.4999908, 0.5045568, 0.5089626, 0.513609, 0.5179914, 0.5224107, 0.5270571, 0.5316316, 0.5358817, 0.5406154, 0.5449819, 0.5494939, 0.5539072, 0.5585454, 0.5628218, 0.5672572, 0.5715997, 0.5761537, 0.5806764, 0.5849425, 0.5894919, 0.5938966, 0.5980922, 0.6026114, 0.6070208, 0.6113653, 0.6158023, 0.620046, 0.6244641, 0.6287679, 0.6330909, 0.6372078, 0.6416281, 0.6458655, 0.6502105, 0.6545821, 0.6591536, 0.6628505, 0.6674164, 0.6713315, 0.6757435, 0.6796968, 0.6840288, 0.6879099, 0.6919787, 0.696522, 0.7004316, 0.7045769, 0.7086459, 0.7125858, 0.7165933, 0.7208507, 0.7248217, 0.7287137, 0.7325798, 0.7365798, 0.740519, 0.7442868, 0.7485558, 0.7521477, 0.7562251, 0.7599414, 0.7636527, 0.7675913, 0.7716459, 0.7751827, 0.7787592, 0.7825552, 0.7863236, 0.7900368, 0.7936389, 0.7976046, 0.8007295, 0.8045919, 0.8080387, 0.8118094, 0.8152153, 0.8188228, 0.8223331, 0.8253756, 0.8293746, 0.832878, 0.8361376, 0.8395188, 0.8430709, 0.8456728, 0.8503147, 0.8533636, 0.8564858, 0.8596513, 0.8629178, 0.8661277, 0.8695821, 0.8726692, 0.8757896, 0.878661, 0.8822673, 0.885249, 0.8885515, 0.8914622, 0.8945625, 0.8974175, 0.8998414, 0.90372, 0.9063937, 0.9094203, 0.9122061, 0.9151366, 0.9177974, 0.9208758, 0.9241818, 0.9266687, 0.9293524, 0.9322096, 0.9350028, 0.9377433, 0.9404456, 0.943255, 0.9460236, 0.9479497, 0.9514073, 0.9538587, 0.956231, 0.9590347, 0.9617073, 0.9636924, 0.9664175, 0.9695347, 0.9717275, 0.9737724, 0.976324, 0.9787543, 0.9811614, 0.9837028, 0.9861388, 0.987884, 0.990298, 0.993163, 0.9952777, 0.9976484, 0.9998569, 1.002022, 1.003855, 1.006214, 1.008713, 1.011002, 1.012539, 1.014963, 1.017067, 1.019199, 1.021274, 1.023492, 1.025832, 1.02712, 1.029613, 1.031582, 1.033471, 1.035147, 1.03719, 1.039011, 1.04073, 1.064208, 1.044812, 1.089229, 1.047822, 1.093093, 1.052054, 1.097064, 1.055375, 1.100804, 1.059046, 1.104272, 1.062371, 1.108317, 1.065671, 1.11113, 1.068859, 1.115129, 1.072524, 1.118716, 1.075476, 1.121811, 1.078137, 1.125391, 1.082756, 1.128421, 1.083183, 1.131349, 1.087014, 1.13422, 1.090216, 1.13763, 1.139284, 1.093966, 1.142445, 1.097526, 1.14492, 1.099885, 1.148008, 1.102422, 1.103554, 1.152076, 1.106593, 1.154926, 1.108584, 1.157599, 1.111422, 1.112604, 1.161706, 1.114946, 1.164241, 1.117085, 1.166565, 1.168032, 1.120806, 1.170148, 1.122306, 1.172595, 1.124713, 1.126141, 1.17618, 1.128243, 1.178522, 1.179616, 1.130952, 1.181801, 1.133497, 1.133977, 1.184923, 1.135647, 1.186912, 1.137921, 1.138773, 1.190669, 1.140788, 1.141489, 1.193606, 1.14343, 1.195445, 1.145018, 1.14561, 1.197972, 1.147034, 1.148066, 1.201151, 1.150078, 1.150694, 1.203385, 1.152055, 1.153149, 1.206266, 1.154696, 1.207566, 1.208715, 1.155979, 1.210081, 1.211206, 1.158722, 1.21294, 1.213517, 1.160465, 1.215236, 1.216005, 1.16256, 1.217116, 1.218016, 1.164129, 1.164481, 1.220085, 1.166308, 1.167011, 1.222205, 1.167173, 1.168678, 1.224423, 1.169415, 1.17011, 1.225865, 1.170719, 1.17113, 1.227874, 1.228607, 1.173088, 1.173673, 1.229782, 1.17489, 1.174992, 
1e-100, 7.466408e-05, 0.0001568645, 0.0002443591, 0.0003416547, 0.00044475, 0.0005553395, 0.0006762512, 0.0008034053, 0.0009415532, 0.001089914, 0.001245691, 0.001415393, 0.001594716, 0.001782271, 0.001987489, 0.002201996, 0.002427279, 0.002671705, 0.002926955, 0.003194009, 0.00348224, 0.003783484, 0.00409871, 0.004434483, 0.004784327, 0.005151921, 0.005537767, 0.005950146, 0.006372388, 0.006817539, 0.007284601, 0.007772447, 0.008281315, 0.00882023, 0.00937644, 0.009954898, 0.01055528, 0.01119207, 0.01184418, 0.01251635, 0.01322951, 0.01396373, 0.01472577, 0.01551729, 0.01633827, 0.01719322, 0.01805919, 0.01898106, 0.01993015, 0.02089776, 0.02190372, 0.02294867, 0.02401817, 0.02512404, 0.02627069, 0.0274562, 0.02866198, 0.02988447, 0.03116172, 0.03248268, 0.0338174, 0.03520749, 0.03661188, 0.03808361, 0.03955336, 0.04108307, 0.0426212, 0.04423104, 0.04587891, 0.04753144, 0.04921185, 0.05098716, 0.05278741, 0.0545898, 0.05644325, 0.05836878, 0.06030489, 0.0622954, 0.06428823, 0.06637513, 0.06847255, 0.0706101, 0.07283794, 0.07502695, 0.0773038, 0.07963088, 0.08202674, 0.08439379, 0.08684217, 0.08929107, 0.09184189, 0.09443338, 0.09699061, 0.09963014, 0.1023081, 0.1050367, 0.1078171, 0.1106542, 0.113494, 0.116397, 0.1193437, 0.1222831, 0.1253544, 0.1283563, 0.1314494, 0.1346071, 0.1377765, 0.1410048, 0.1442628, 0.1475398, 0.1509138, 0.1542356, 0.1577487, 0.1611084, 0.1646749, 0.1681304, 0.171804, 0.1754804, 0.1790785, 0.1827496, 0.1864768, 0.190278, 0.1940909, 0.1978629, 0.2018356, 0.2056891, 0.2096428, 0.2136941, 0.2177589, 0.2218416, 0.2259213, 0.2299674, 0.2341747, 0.2384414, 0.2426567, 0.246907, 0.251216, 0.2555996, 0.2599543, 0.2643686, 0.2688043, 0.2733247, 0.2777924, 0.2822323, 0.2868249, 0.2914684, 0.2961114, 0.3007294, 0.3053739, 0.3101911, 0.3149007, 0.3197097, 0.3244435, 0.3293704, 0.3342848, 0.3390737, 0.3438719, 0.3488416, 0.3539115, 0.3588649, 0.3639149, 0.3689023, 0.373812, 0.378801, 0.3839669, 0.3891356, 0.3942904, 0.3993969, 0.4047019, 0.409696, 0.4149921, 0.420368, 0.4255838, 0.430716, 0.4361775, 0.4413834, 0.446626, 0.4520881, 0.4573786, 0.4626441, 0.4681221, 0.4735007, 0.4788832, 0.4842993, 0.4898089, 0.4952327, 0.5006236, 0.5064196, 0.5116554, 0.5170969, 0.5224007, 0.5280669, 0.5335089, 0.5390668, 0.5446845, 0.5504258, 0.5556907, 0.5616772, 0.5667844, 0.5717093, 0.5774561, 0.5833886, 0.5887787, 0.594474, 0.6000698, 0.605629, 0.6109378, 0.6165776, 0.622045, 0.627683, 0.6334448, 0.6389901, 0.6444618, 0.6500292, 0.655369, 0.6610356, 0.6666209, 0.672244, 0.677774, 0.6833139, 0.68868, 0.6940133, 0.6995016, 0.7050959, 0.710594, 0.7161625, 0.7215918, 0.7267825, 0.7324146, 0.7375518, 0.7433644, 0.7484058, 0.7540024, 0.7593138, 0.7651024, 0.7699641, 0.7754664, 0.780935, 0.7859036, 0.7913372, 0.7966731, 0.8019389, 0.8073088, 0.8123716, 0.8175957, 0.8228549, 0.82794, 0.8332912, 0.8385231, 0.8435693, 0.8486357, 0.8537336, 0.8585932, 0.8640213, 0.8688088, 0.8737584, 0.8788769, 0.8837339, 0.8889373, 0.8940431, 0.8988342, 0.9036476, 0.9084439, 0.9132007, 0.9183285, 0.9225318, 0.9279179, 0.9322275, 0.937096, 0.9418739, 0.946675, 0.9514495, 0.9559446, 0.96062, 0.9650642, 0.9698392, 0.974284, 0.9786577, 0.9831502, 0.988124, 0.9921513, 0.9969301, 1.001011, 1.005573, 1.009597, 1.014013, 1.018582, 1.02214, 1.027088, 1.031019, 1.035437, 1.039558, 1.044054, 1.047787, 1.051583, 1.056075, 1.06002, 1.063247, 1.069186, 1.072238, 1.076358, 1.080167, 1.083943, 1.08743, 1.091673, 1.095327, 1.099328, 1.1031, 1.107065, 1.110479, 1.113655, 1.118111, 1.121359, 1.124973, 1.128514, 1.132247, 1.135407, 1.138974, 1.142911, 1.146078, 1.149482, 1.152843, 1.156285, 1.15972, 1.162843, 1.166697, 1.169777, 1.172571, 1.17626, 1.179305, 1.18264, 1.185183, 1.188852, 1.191815, 1.194747, 1.198211, 1.201079, 1.203861, 1.206996, 1.209886, 1.213132, 1.215649, 1.219135, 1.221691, 1.22377, 1.227529, 1.229946, 1.232219, 1.235619, 1.238525, 1.240728, 1.243437, 1.246639, 1.249025, 1.25127, 1.254226, 1.256379, 1.259493, 1.261708, 1.264414, 1.266814, 1.26851, 1.271919, 1.274067, 1.276897, 1.27885, 1.281236, 1.28353, 1.285772, 1.287875, 1.290753, 1.292494, 1.294113, 1.325262, 1.299294, 1.35838, 1.303745, 1.362823, 1.307794, 1.366726, 1.312213, 1.371751, 1.31586, 1.376052, 1.319826, 1.380258, 1.324078, 1.384265, 1.328161, 1.388924, 1.330936, 1.393178, 1.335332, 1.396675, 1.339019, 1.400648, 1.34357, 1.404317, 1.344445, 1.408245, 1.349157, 1.411003, 1.352733, 1.415608, 1.417654, 1.357542, 1.420859, 1.360955, 1.423631, 1.3637, 1.427695, 1.367193, 1.368059, 1.432569, 1.370975, 1.435905, 1.374409, 1.439237, 1.377478, 1.37877, 1.443588, 1.381803, 1.446853, 1.384132, 1.449423, 1.451284, 1.387868, 1.454018, 1.390737, 1.457003, 1.392957, 1.3939, 1.461559, 1.397291, 1.463776, 1.465424, 1.40037, 1.467157, 1.402525, 1.403829, 1.471674, 1.406141, 1.474102, 1.407557, 1.409741, 1.477916, 1.411655, 1.41273, 1.481739, 1.414525, 1.482959, 1.416271, 1.417575, 1.486996, 1.419212, 1.419539, 1.490216, 1.422419, 1.423549, 1.49361, 1.425067, 1.426046, 1.495965, 1.4275, 1.498397, 1.499617, 1.429847, 1.500721, 1.502364, 1.432931, 1.504472, 1.505181, 1.4348, 1.507065, 1.507757, 1.436725, 1.509614, 1.510891, 1.438715, 1.439057, 1.513193, 1.441014, 1.442035, 1.515672, 1.443235, 1.443819, 1.517971, 1.4445, 1.445492, 1.520253, 1.446571, 1.446935, 1.52167, 1.523068, 1.449218, 1.450203, 1.524594, 1.451214, 1.451489, 
1e-100, 0.0001295729, 0.0002702666, 0.0004237587, 0.0005855792, 0.0007620041, 0.0009495985, 0.001146908, 0.001363321, 0.001590349, 0.001831799, 0.002091517, 0.002363275, 0.002652033, 0.002961542, 0.00328352, 0.003625957, 0.003989211, 0.004367137, 0.004770707, 0.005194513, 0.005638972, 0.006104564, 0.006595398, 0.007106641, 0.00764013, 0.008211221, 0.008795626, 0.009409864, 0.01005703, 0.01072815, 0.01141927, 0.01215301, 0.01292191, 0.01369466, 0.01451447, 0.01537328, 0.01626534, 0.01717301, 0.01813084, 0.01912926, 0.02014418, 0.02121195, 0.02231294, 0.02344856, 0.02460999, 0.02583613, 0.02708681, 0.02837657, 0.0297087, 0.03108806, 0.03251429, 0.03397159, 0.03547183, 0.03701954, 0.03861928, 0.04026176, 0.04194694, 0.04366805, 0.04543598, 0.04721901, 0.04909254, 0.05099977, 0.05295468, 0.05491406, 0.05695481, 0.05903488, 0.06118568, 0.0633458, 0.06553127, 0.06782129, 0.07015835, 0.07252519, 0.07492729, 0.07738348, 0.07988051, 0.08247291, 0.08507877, 0.08771811, 0.0903932, 0.09322699, 0.09602464, 0.09885191, 0.1018185, 0.1047632, 0.1077694, 0.1108478, 0.1139374, 0.1171291, 0.1203841, 0.1236731, 0.1269619, 0.130365, 0.1337792, 0.1372953, 0.1408218, 0.1443843, 0.1480086, 0.1516852, 0.1554317, 0.1592211, 0.1630012, 0.1668629, 0.1708355, 0.1747911, 0.1788129, 0.1828561, 0.1870139, 0.1912377, 0.1954141, 0.1996916, 0.2040013, 0.2083301, 0.2127589, 0.2172234, 0.2217843, 0.226332, 0.2309427, 0.2356061, 0.240224, 0.2451221, 0.2498092, 0.2546972, 0.2594987, 0.2644327, 0.2694697, 0.2745128, 0.2795096, 0.2846674, 0.2897809, 0.2949469, 0.3001565, 0.3055501, 0.3108483, 0.3161029, 0.3213831, 0.3269689, 0.3325075, 0.3378408, 0.3433369, 0.3490849, 0.354587, 0.3602943, 0.3659929, 0.3718612, 0.3775833, 0.3833531, 0.3891846, 0.3949859, 0.4010236, 0.4070754, 0.4128104, 0.4189527, 0.4250529, 0.4310576, 0.4371795, 0.4432257, 0.4493285, 0.4555046, 0.4619092, 0.4679797, 0.4741993, 0.4803717, 0.4869753, 0.4931235, 0.4996178, 0.5061757, 0.5124195, 0.5188042, 0.5254302, 0.5318992, 0.5382036, 0.5449623, 0.5513037, 0.5577275, 0.5643544, 0.5711371, 0.5777075, 0.5842867, 0.5907965, 0.5977255, 0.6041175, 0.6109285, 0.6178167, 0.6245092, 0.631148, 0.6378372, 0.6446795, 0.6513127, 0.6578026, 0.6648544, 0.6714476, 0.6783276, 0.6852566, 0.6918862, 0.6987636, 0.7056877, 0.7125528, 0.7198352, 0.7259045, 0.7319895, 0.7393174, 0.746166, 0.7530497, 0.7597332, 0.7669172, 0.7734913, 0.780038, 0.7869085, 0.7940491, 0.8002493, 0.807663, 0.8142385, 0.8209628, 0.8277331, 0.8346882, 0.8412048, 0.8480398, 0.8548987, 0.8615853, 0.8679433, 0.8750053, 0.8813784, 0.8882668, 0.8953854, 0.9015243, 0.9082261, 0.9149058, 0.9211828, 0.9280393, 0.9347096, 0.9413559, 0.9476694, 0.9542194, 0.9607774, 0.9673867, 0.9738772, 0.980284, 0.9868991, 0.9927239, 0.9996895, 1.005985, 1.012201, 1.018493, 1.025172, 1.030918, 1.037398, 1.043816, 1.049832, 1.056588, 1.06257, 1.068509, 1.074582, 1.080524, 1.086771, 1.09282, 1.098669, 1.104945, 1.110694, 1.116898, 1.122894, 1.128499, 1.134168, 1.140679, 1.146153, 1.152259, 1.157471, 1.163387, 1.169249, 1.174828, 1.180316, 1.18577, 1.191826, 1.197082, 1.202667, 1.208316, 1.213826, 1.219391, 1.224381, 1.229792, 1.235361, 1.240263, 1.245767, 1.251546, 1.256079, 1.261808, 1.266448, 1.271941, 1.277088, 1.281942, 1.286845, 1.291329, 1.297042, 1.30189, 1.306784, 1.311871, 1.316624, 1.321238, 1.325032, 1.332265, 1.335567, 1.34043, 1.344911, 1.349692, 1.354196, 1.359214, 1.363485, 1.36719, 1.372435, 1.376921, 1.381318, 1.385731, 1.390063, 1.393828, 1.397666, 1.402715, 1.406828, 1.411051, 1.41518, 1.419243, 1.423193, 1.427254, 1.432089, 1.435385, 1.439062, 1.443283, 1.44743, 1.451406, 1.454909, 1.458989, 1.462499, 1.466325, 1.470082, 1.473666, 1.477326, 1.481196, 1.484789, 1.488095, 1.491408, 1.495847, 1.498696, 1.502278, 1.505419, 1.509177, 1.512568, 1.515612, 1.519133, 1.52239, 1.524953, 1.529705, 1.532159, 1.535334, 1.538574, 1.541708, 1.544347, 1.547543, 1.551073, 1.553597, 1.556702, 1.559959, 1.562974, 1.566038, 1.568728, 1.571442, 1.574231, 1.576593, 1.580105, 1.582599, 1.585125, 1.587862, 1.590569, 1.593115, 1.596164, 1.59833, 1.638003, 1.60363, 1.678594, 1.609248, 1.684535, 1.613822, 1.690247, 1.618354, 1.694991, 1.623228, 1.70073, 1.628008, 1.705354, 1.632672, 1.710126, 1.637317, 1.715667, 1.641789, 1.720177, 1.645671, 1.724466, 1.650604, 1.730037, 1.655788, 1.733491, 1.656525, 1.73758, 1.662114, 1.742575, 1.666511, 1.746693, 1.748036, 1.672368, 1.753132, 1.675878, 1.757005, 1.679358, 1.760846, 1.682891, 1.684773, 1.767127, 1.688575, 1.77042, 1.691531, 1.775202, 1.695177, 1.696569, 1.78081, 1.699753, 1.78315, 1.702692, 1.787095, 1.78917, 1.707417, 1.792084, 1.710361, 1.796021, 1.7135, 1.715418, 1.80081, 1.717877, 1.803036, 1.805684, 1.722335, 1.808118, 1.724203, 1.725358, 1.81249, 1.728491, 1.816264, 1.730983, 1.732526, 1.820281, 1.733885, 1.736964, 1.824952, 1.738505, 1.826571, 1.740398, 1.740934, 1.830942, 1.744142, 1.745103, 1.835225, 1.747018, 1.747732, 1.839056, 1.750209, 1.751899, 1.842303, 1.753494, 1.843985, 1.845707, 1.756223, 1.848157, 1.849449, 1.75904, 1.850674, 1.853415, 1.762103, 1.854474, 1.856133, 1.764276, 1.857332, 1.858345, 1.766497, 1.767719, 1.862501, 1.769281, 1.769432, 1.864757, 1.77169, 1.772175, 1.867828, 1.773443, 1.774332, 1.869552, 1.775254, 1.775924, 1.872772, 1.874145, 1.778274, 1.778234, 1.876114, 1.781158, 1.780845, 
1e-100, 0.0002128409, 0.0004409143, 0.0006833458, 0.0009465891, 0.001223488, 0.001519654, 0.00183589, 0.002166952, 0.002524066, 0.002900327, 0.003292821, 0.003717984, 0.004160067, 0.004623569, 0.005120097, 0.005636749, 0.00617769, 0.006755296, 0.007353899, 0.007976769, 0.008646584, 0.009330944, 0.01004869, 0.01080902, 0.01159489, 0.01241139, 0.01327174, 0.01416523, 0.01508848, 0.0160609, 0.01706733, 0.0181089, 0.01919546, 0.02032641, 0.02149007, 0.02268507, 0.02395533, 0.02525912, 0.02659243, 0.0279766, 0.02942114, 0.03090428, 0.03241339, 0.03401031, 0.03565165, 0.03732308, 0.03905714, 0.04084707, 0.04269103, 0.04457509, 0.04652335, 0.04852569, 0.05057997, 0.05268043, 0.05485671, 0.05708571, 0.05938182, 0.06167614, 0.06409228, 0.06652458, 0.06903198, 0.07156782, 0.07418219, 0.07688546, 0.07958579, 0.08238421, 0.08520659, 0.08810628, 0.0911017, 0.0940754, 0.0971348, 0.1003245, 0.1035325, 0.1067807, 0.1100875, 0.1134555, 0.1169319, 0.1204609, 0.1239896, 0.1276636, 0.1313533, 0.1350603, 0.1389074, 0.1427863, 0.1466988, 0.1507528, 0.15479, 0.1589919, 0.163093, 0.167409, 0.1717051, 0.1761476, 0.1805381, 0.1850913, 0.1896112, 0.1942987, 0.1989325, 0.2037195, 0.2084317, 0.2134913, 0.2184394, 0.2234084, 0.2283662, 0.2334026, 0.2386443, 0.243953, 0.2492095, 0.2545613, 0.2598197, 0.2654385, 0.2709177, 0.2765684, 0.2821232, 0.2878582, 0.293668, 0.2993844, 0.305343, 0.3112346, 0.3171746, 0.3232249, 0.3292278, 0.3354432, 0.3416879, 0.3479174, 0.3541293, 0.3605814, 0.3670786, 0.3734755, 0.3799395, 0.3864665, 0.3931328, 0.3997288, 0.4062557, 0.4130809, 0.4198974, 0.4267462, 0.4335709, 0.4406339, 0.4475061, 0.4545508, 0.4617677, 0.468748, 0.4759596, 0.4831548, 0.4902171, 0.497609, 0.5051316, 0.512417, 0.5195314, 0.5271371, 0.5345637, 0.5422323, 0.5498461, 0.5572126, 0.564603, 0.5723619, 0.5800429, 0.5875815, 0.5954882, 0.6031396, 0.6109477, 0.6187685, 0.6264461, 0.6345494, 0.6422282, 0.6504746, 0.658458, 0.6663402, 0.6741507, 0.6822609, 0.6901727, 0.6983202, 0.7064183, 0.7143133, 0.722336, 0.7305844, 0.7386848, 0.7467513, 0.7549734, 0.7633926, 0.7715321, 0.7795123, 0.7881462, 0.7961, 0.8044053, 0.812468, 0.8208266, 0.829087, 0.8373003, 0.8457085, 0.8540994, 0.8620116, 0.8706894, 0.8788689, 0.8872768, 0.8955911, 0.9040026, 0.9126079, 0.9202027, 0.9279311, 0.9368647, 0.9450061, 0.9535475, 0.9619493, 0.9698775, 0.9779862, 0.9867854, 0.9945232, 1.003352, 1.011369, 1.019623, 1.027888, 1.035798, 1.044022, 1.052287, 1.060621, 1.068841, 1.077189, 1.085317, 1.092906, 1.101099, 1.109296, 1.116848, 1.126067, 1.133616, 1.141129, 1.149299, 1.157515, 1.164867, 1.173623, 1.181211, 1.188679, 1.196905, 1.205145, 1.212587, 1.220417, 1.227966, 1.236013, 1.243421, 1.25124, 1.258966, 1.26643, 1.274051, 1.281636, 1.289362, 1.296754, 1.304579, 1.311281, 1.319774, 1.326752, 1.333734, 1.341241, 1.349171, 1.355676, 1.363177, 1.37014, 1.377211, 1.384505, 1.39177, 1.398967, 1.405592, 1.412795, 1.419702, 1.426513, 1.433816, 1.440708, 1.446916, 1.454023, 1.460416, 1.46749, 1.474364, 1.481147, 1.48761, 1.493839, 1.500738, 1.507271, 1.513568, 1.520004, 1.526663, 1.532609, 1.539216, 1.545618, 1.551815, 1.557789, 1.563715, 1.569987, 1.575398, 1.582322, 1.588042, 1.594046, 1.600146, 1.606108, 1.611519, 1.61742, 1.623533, 1.628838, 1.634463, 1.639282, 1.64739, 1.65206, 1.657659, 1.662889, 1.667718, 1.673942, 1.679282, 1.68494, 1.690096, 1.695415, 1.700284, 1.704744, 1.711157, 1.715849, 1.720522, 1.725884, 1.730947, 1.735575, 1.740425, 1.746316, 1.750352, 1.754847, 1.759983, 1.764775, 1.769518, 1.774039, 1.779151, 1.783281, 1.788351, 1.793115, 1.797215, 1.801552, 1.80654, 1.80964, 1.814484, 1.818825, 1.823729, 1.827265, 1.831352, 1.835523, 1.839913, 1.844342, 1.848387, 1.852529, 1.855896, 1.860251, 1.864533, 1.868142, 1.87162, 1.875312, 1.879965, 1.882646, 1.886724, 1.891382, 1.894434, 1.897547, 1.901249, 1.904974, 1.908634, 1.912752, 1.915654, 1.919083, 1.921311, 1.92599, 1.928945, 1.932303, 1.936513, 1.939091, 1.941825, 1.94535, 1.948955, 1.951766, 1.954513, 1.957733, 1.960607, 2.011457, 1.966901, 2.065243, 1.973252, 2.071068, 1.978628, 2.07786, 1.984561, 2.084098, 1.989855, 2.089675, 1.995145, 2.096106, 2.000893, 2.101897, 2.006049, 2.107501, 2.011925, 2.113525, 2.016231, 2.118608, 2.021156, 2.123777, 2.027795, 2.129439, 2.029105, 2.133908, 2.035185, 2.13977, 2.040016, 2.144473, 2.147357, 2.046709, 2.15161, 2.051077, 2.156856, 2.05561, 2.161269, 2.058593, 2.060848, 2.169064, 2.066157, 2.173135, 2.069718, 2.177819, 2.072894, 2.075399, 2.184265, 2.079277, 2.187966, 2.082366, 2.191356, 2.194158, 2.088057, 2.198445, 2.091763, 2.201993, 2.0938, 2.096947, 2.208054, 2.099925, 2.211361, 2.213856, 2.104161, 2.21672, 2.107937, 2.109581, 2.222824, 2.111932, 2.225722, 2.11551, 2.117033, 2.231103, 2.119635, 2.121303, 2.235528, 2.123564, 2.23905, 2.125965, 2.127797, 2.243373, 2.129516, 2.131995, 2.249139, 2.13371, 2.135412, 2.252913, 2.137032, 2.138412, 2.257028, 2.141476, 2.259271, 2.260964, 2.144178, 2.264029, 2.265651, 2.14769, 2.267802, 2.269318, 2.150321, 2.270774, 2.273049, 2.153812, 2.274733, 2.276301, 2.155746, 2.15665, 2.280699, 2.159669, 2.160252, 2.283739, 2.161745, 2.161802, 2.286227, 2.164198, 2.165079, 2.289688, 2.165825, 2.167611, 2.29296, 2.294734, 2.169753, 2.17052, 2.296612, 2.172325, 2.171857, 
1e-100, 0.0003239487, 0.000673572, 0.001046768, 0.001438265, 0.001862176, 0.002306509, 0.002775896, 0.003277656, 0.003801381, 0.004358571, 0.004948178, 0.005560918, 0.006213004, 0.006894598, 0.00760801, 0.008362572, 0.009151391, 0.009970419, 0.01084005, 0.01174612, 0.01267937, 0.01367407, 0.01470382, 0.01576341, 0.01687963, 0.01804709, 0.01924616, 0.02049686, 0.02181015, 0.02316079, 0.02455666, 0.02602002, 0.02752703, 0.02907998, 0.03070293, 0.03238056, 0.03409885, 0.0358759, 0.03773195, 0.03963333, 0.04158938, 0.04361708, 0.04571761, 0.0478629, 0.05005469, 0.05234365, 0.05469639, 0.05710237, 0.05955521, 0.06210825, 0.06474428, 0.06742281, 0.07018205, 0.07297049, 0.07589127, 0.07882886, 0.08188812, 0.08499083, 0.08814961, 0.09140967, 0.09470668, 0.09809317, 0.1015552, 0.1050315, 0.1086675, 0.1123357, 0.1160691, 0.119909, 0.1236909, 0.1277123, 0.1317461, 0.1358841, 0.1400354, 0.1443316, 0.1486338, 0.153063, 0.1575705, 0.1620887, 0.1666763, 0.1714648, 0.176245, 0.1810964, 0.1861277, 0.191109, 0.1961142, 0.2013562, 0.2065261, 0.2119488, 0.2173532, 0.222889, 0.2283813, 0.2340011, 0.239776, 0.2455542, 0.25134, 0.2573076, 0.2632364, 0.2694075, 0.2754386, 0.2817916, 0.2880115, 0.294302, 0.3008181, 0.3073214, 0.3139706, 0.3205309, 0.3272631, 0.3340608, 0.3409009, 0.3478867, 0.3547742, 0.3619204, 0.3690881, 0.3763223, 0.3836425, 0.3907119, 0.3983739, 0.4056755, 0.4134221, 0.4209618, 0.428688, 0.4364406, 0.4440295, 0.452003, 0.4599905, 0.4680028, 0.4760075, 0.4839742, 0.4922734, 0.5004663, 0.50859, 0.5168669, 0.5253007, 0.5337343, 0.5422115, 0.550767, 0.5594764, 0.5680863, 0.576813, 0.585365, 0.5944404, 0.6033303, 0.611995, 0.6210273, 0.6300691, 0.6390786, 0.6481924, 0.657166, 0.6660488, 0.6756104, 0.6849255, 0.693983, 0.7032589, 0.7126755, 0.7221114, 0.7315186, 0.7408946, 0.7504547, 0.7600254, 0.7694783, 0.7792101, 0.788829, 0.7982846, 0.808286, 0.8180209, 0.8277369, 0.8374192, 0.8471012, 0.856824, 0.8665945, 0.8765083, 0.8863101, 0.8961103, 0.9057537, 0.9162446, 0.9258866, 0.9358865, 0.9461347, 0.95596, 0.9658902, 0.9759204, 0.9859812, 0.9960165, 1.006082, 1.015993, 1.025806, 1.035989, 1.046277, 1.056407, 1.066514, 1.076718, 1.086667, 1.096716, 1.106795, 1.116773, 1.126759, 1.136968, 1.147798, 1.156505, 1.166052, 1.17657, 1.186592, 1.196509, 1.206808, 1.2167, 1.226719, 1.236898, 1.247064, 1.256489, 1.266853, 1.27653, 1.286502, 1.296303, 1.306771, 1.315978, 1.326279, 1.336045, 1.345304, 1.355296, 1.365454, 1.3748, 1.384854, 1.394354, 1.404353, 1.413345, 1.423145, 1.43271, 1.442592, 1.452028, 1.461623, 1.470873, 1.479579, 1.489795, 1.499325, 1.508575, 1.517306, 1.527074, 1.535929, 1.545739, 1.554803, 1.563625, 1.572852, 1.582507, 1.590768, 1.600022, 1.608942, 1.617958, 1.626806, 1.636418, 1.644646, 1.653033, 1.661995, 1.670723, 1.679509, 1.687692, 1.696505, 1.70457, 1.713591, 1.722007, 1.730473, 1.738671, 1.747377, 1.755706, 1.763363, 1.772257, 1.779709, 1.788082, 1.796011, 1.804071, 1.811144, 1.820471, 1.82795, 1.835832, 1.843911, 1.851392, 1.85917, 1.866447, 1.874452, 1.881673, 1.888862, 1.896528, 1.904619, 1.910984, 1.919383, 1.925802, 1.933441, 1.9405, 1.947509, 1.954239, 1.961276, 1.969257, 1.975671, 1.982693, 1.989584, 1.996031, 2.003, 2.009595, 2.015154, 2.024265, 2.029908, 2.036327, 2.043176, 2.049237, 2.0562, 2.061972, 2.067308, 2.075067, 2.080816, 2.086889, 2.09288, 2.099646, 2.104808, 2.110445, 2.117424, 2.122785, 2.128543, 2.134589, 2.140064, 2.14595, 2.151191, 2.157523, 2.162686, 2.168291, 2.173897, 2.17956, 2.184516, 2.190253, 2.195803, 2.200248, 2.206457, 2.211451, 2.215455, 2.220892, 2.226196, 2.230955, 2.236175, 2.241197, 2.246063, 2.250647, 2.25553, 2.260648, 2.265349, 2.26923, 2.274611, 2.279095, 2.283312, 2.287848, 2.292232, 2.297546, 2.301053, 2.305902, 2.309869, 2.314508, 2.318313, 2.323243, 2.327152, 2.329948, 2.336207, 2.339253, 2.343293, 2.347208, 2.351076, 2.354622, 2.358516, 2.363342, 2.366821, 2.369981, 2.374037, 2.377075, 2.381218, 2.385178, 2.388712, 2.391945, 2.454981, 2.398779, 2.523166, 2.406132, 2.529915, 2.412974, 2.537093, 2.418604, 2.544672, 2.426194, 2.551784, 2.431429, 2.558741, 2.438998, 2.565755, 2.444361, 2.572947, 2.450869, 2.578116, 2.456092, 2.585116, 2.462258, 2.591202, 2.4692, 2.596551, 2.471223, 2.60325, 2.47814, 2.609568, 2.483745, 2.615033, 2.618934, 2.491507, 2.624256, 2.496675, 2.629072, 2.500987, 2.635213, 2.505756, 2.508983, 2.643468, 2.512389, 2.649436, 2.517913, 2.653636, 2.522669, 2.524744, 2.661048, 2.528543, 2.666617, 2.533172, 2.670312, 2.673097, 2.539188, 2.678545, 2.544033, 2.682257, 2.547057, 2.549226, 2.688826, 2.552874, 2.693894, 2.696471, 2.558144, 2.698699, 2.56224, 2.56421, 2.707233, 2.566919, 2.710308, 2.570823, 2.571931, 2.715856, 2.575724, 2.578167, 2.721827, 2.579533, 2.724771, 2.582738, 2.585427, 2.731509, 2.588146, 2.589565, 2.736781, 2.592488, 2.593407, 2.742549, 2.59641, 2.596954, 2.746081, 2.599602, 2.749484, 2.752011, 2.604501, 2.754763, 2.755878, 2.607487, 2.759689, 2.760773, 2.611507, 2.762783, 2.765317, 2.613297, 2.767113, 2.769971, 2.618044, 2.618913, 2.773664, 2.620079, 2.622694, 2.777674, 2.624816, 2.625085, 2.781208, 2.625958, 2.626818, 2.784927, 2.630021, 2.63105, 2.788612, 2.789845, 2.633308, 2.634223, 2.792802, 2.636794, 2.636136, 
1e-100, 0.0004845086, 0.0009956963, 0.001539163, 0.002118289, 0.002724798, 0.003372553, 0.004052961, 0.00476527, 0.005528678, 0.006322712, 0.00715729, 0.008041536, 0.008960092, 0.009923574, 0.01094098, 0.01199559, 0.01310401, 0.01426615, 0.01546798, 0.01672868, 0.01804572, 0.01941507, 0.02083258, 0.02232267, 0.02385481, 0.02545005, 0.02712062, 0.02883093, 0.03061153, 0.03247069, 0.03439777, 0.03636283, 0.03842026, 0.04054366, 0.04272589, 0.04499338, 0.0473243, 0.04974068, 0.05220338, 0.05476708, 0.05741974, 0.0601187, 0.06289291, 0.06580884, 0.06874593, 0.07175433, 0.07486573, 0.07811502, 0.081361, 0.08473432, 0.08820842, 0.09173232, 0.09539554, 0.09908215, 0.1028877, 0.1068084, 0.1107925, 0.1148239, 0.1190143, 0.1232671, 0.1275615, 0.1319627, 0.1364643, 0.1411182, 0.1457734, 0.1505482, 0.1554072, 0.1603449, 0.1654322, 0.170513, 0.1756935, 0.1811418, 0.1865187, 0.191988, 0.1975425, 0.2032324, 0.2090842, 0.2149634, 0.2208584, 0.226965, 0.233121, 0.2393329, 0.2456905, 0.2520815, 0.258636, 0.2653175, 0.2719897, 0.278794, 0.2856665, 0.2926416, 0.2997284, 0.3069038, 0.314141, 0.3214656, 0.3290318, 0.336613, 0.34415, 0.3517725, 0.3594046, 0.3676678, 0.3755485, 0.3835751, 0.391765, 0.3998018, 0.4082061, 0.4164579, 0.4251231, 0.4336027, 0.4422332, 0.4510004, 0.4598452, 0.4688651, 0.4776098, 0.4865982, 0.4957801, 0.5048973, 0.5143648, 0.5236528, 0.5331621, 0.5426777, 0.5522261, 0.5619563, 0.5716191, 0.5814751, 0.5912956, 0.6010622, 0.6111943, 0.6213352, 0.6314808, 0.6415808, 0.6518473, 0.6622035, 0.6725295, 0.6831395, 0.6936167, 0.704409, 0.7149204, 0.7253686, 0.7364439, 0.7474388, 0.7580457, 0.7691437, 0.7797364, 0.7910483, 0.8021896, 0.8131658, 0.8242542, 0.8354852, 0.8465082, 0.8580069, 0.8693161, 0.8808209, 0.8924646, 0.9039002, 0.9153404, 0.9268595, 0.9382436, 0.9502599, 0.9620126, 0.9732823, 0.9852602, 0.9971947, 1.008864, 1.020566, 1.032607, 1.044649, 1.056176, 1.068175, 1.079951, 1.091855, 1.103988, 1.1161, 1.128031, 1.139887, 1.152293, 1.164244, 1.175931, 1.188756, 1.200313, 1.212714, 1.224447, 1.236895, 1.248937, 1.261041, 1.273253, 1.285819, 1.297667, 1.310035, 1.32223, 1.334391, 1.346367, 1.358699, 1.370645, 1.382785, 1.394966, 1.407582, 1.419383, 1.432392, 1.443292, 1.454174, 1.467129, 1.479419, 1.491024, 1.503464, 1.515621, 1.527797, 1.539496, 1.551496, 1.563703, 1.575224, 1.587216, 1.598851, 1.610988, 1.623133, 1.634735, 1.646577, 1.658431, 1.669402, 1.681926, 1.693141, 1.704998, 1.716292, 1.728135, 1.739774, 1.750933, 1.762701, 1.774364, 1.784865, 1.796995, 1.807999, 1.819237, 1.83033, 1.842321, 1.85305, 1.864287, 1.87497, 1.886231, 1.89736, 1.908185, 1.91951, 1.929443, 1.940995, 1.952121, 1.962508, 1.972449, 1.98379, 1.993691, 2.005099, 2.015937, 2.026029, 2.036216, 2.047443, 2.056906, 2.067122, 2.077081, 2.087081, 2.097352, 2.107948, 2.117759, 2.127207, 2.137557, 2.147178, 2.157126, 2.166536, 2.177642, 2.185547, 2.195591, 2.204811, 2.214835, 2.223689, 2.233911, 2.242813, 2.250981, 2.2616, 2.27013, 2.279461, 2.288547, 2.297764, 2.305828, 2.316001, 2.323844, 2.332979, 2.341442, 2.350008, 2.358186, 2.366102, 2.37571, 2.383907, 2.392085, 2.400297, 2.40865, 2.416752, 2.424543, 2.433401, 2.440507, 2.448368, 2.456234, 2.46297, 2.474223, 2.481073, 2.488356, 2.494725, 2.50339, 2.510639, 2.518205, 2.525807, 2.532881, 2.539893, 2.547195, 2.555071, 2.561733, 2.568389, 2.575251, 2.582262, 2.588746, 2.595934, 2.60346, 2.60931, 2.615906, 2.622441, 2.6294, 2.636252, 2.642422, 2.649206, 2.654584, 2.662236, 2.668396, 2.674072, 2.679906, 2.686361, 2.692157, 2.697292, 2.704514, 2.710516, 2.715129, 2.721361, 2.726214, 2.732721, 2.739028, 2.744576, 2.750643, 2.754327, 2.761528, 2.766634, 2.771061, 2.776625, 2.781759, 2.785842, 2.792853, 2.797386, 2.80293, 2.807707, 2.811927, 2.816993, 2.822459, 2.827489, 2.831786, 2.837111, 2.841131, 2.84535, 2.85095, 2.85485, 2.86002, 2.864028, 2.868884, 2.872594, 2.877446, 2.882073, 2.886141, 2.8894, 2.894504, 2.899024, 2.902608, 2.982802, 2.911002, 3.066742, 2.918363, 3.07551, 2.927349, 3.084463, 2.933404, 3.092268, 2.940877, 3.101073, 2.949121, 3.109673, 2.95679, 3.116163, 2.964275, 3.125464, 2.970257, 3.132376, 2.976817, 3.139338, 2.982712, 3.147326, 2.992512, 3.154285, 2.993118, 3.160834, 3.002907, 3.168421, 3.008696, 3.175388, 3.179099, 3.016293, 3.18484, 3.023584, 3.191105, 3.028489, 3.197306, 3.034159, 3.037594, 3.208039, 3.042947, 3.213705, 3.047586, 3.219525, 3.053409, 3.056471, 3.228251, 3.059194, 3.233777, 3.064071, 3.239433, 3.242998, 3.072812, 3.248598, 3.07609, 3.253182, 3.081558, 3.083761, 3.261217, 3.087948, 3.265196, 3.267821, 3.093936, 3.273313, 3.098728, 3.100279, 3.280245, 3.104576, 3.285962, 3.108224, 3.110189, 3.292036, 3.11288, 3.115053, 3.299621, 3.118502, 3.303249, 3.122419, 3.123047, 3.310771, 3.128314, 3.128997, 3.315534, 3.132689, 3.133919, 3.321447, 3.13711, 3.138557, 3.327569, 3.141973, 3.330171, 3.333972, 3.146478, 3.336907, 3.33889, 3.150268, 3.341745, 3.343111, 3.153971, 3.346829, 3.349409, 3.158165, 3.351199, 3.353824, 3.162118, 3.163875, 3.359277, 3.165205, 3.166444, 3.363151, 3.167828, 3.169779, 3.368268, 3.171399, 3.172479, 3.371359, 3.174609, 3.176395, 3.377048, 3.378208, 3.179083, 3.180075, 3.380594, 3.181596, 3.182593, 
1e-100, 0.0006826238, 0.001416437, 0.002189321, 0.003001559, 0.003868904, 0.00477393, 0.005730792, 0.00674112, 0.007791529, 0.00890739, 0.01007494, 0.01128477, 0.01257258, 0.01390719, 0.01529387, 0.01676585, 0.0182845, 0.01986129, 0.02153289, 0.02325273, 0.02503283, 0.02690759, 0.02884686, 0.03084399, 0.03294297, 0.03510586, 0.03733661, 0.03966596, 0.04206312, 0.0445578, 0.04711747, 0.04980478, 0.05253533, 0.05535109, 0.05829419, 0.06131464, 0.06441011, 0.06760441, 0.07090846, 0.07431506, 0.07776088, 0.08139026, 0.08507856, 0.08881935, 0.09274509, 0.09674327, 0.100838, 0.1050258, 0.109353, 0.1137843, 0.1182822, 0.1228816, 0.1276451, 0.1325504, 0.1374384, 0.1424922, 0.147713, 0.1529748, 0.1583794, 0.1638837, 0.1694841, 0.1751896, 0.1810669, 0.1869044, 0.1930092, 0.1992173, 0.2054139, 0.2118049, 0.2182479, 0.2248792, 0.2315485, 0.2384512, 0.2453931, 0.2524435, 0.2596192, 0.2668311, 0.2742171, 0.2817038, 0.2893531, 0.2971435, 0.3049741, 0.312851, 0.3210165, 0.3291437, 0.3373583, 0.3457109, 0.3542861, 0.3629386, 0.3717013, 0.3805873, 0.3895104, 0.3985736, 0.4077259, 0.4170356, 0.4264288, 0.4358169, 0.4453202, 0.4550565, 0.4647482, 0.4751565, 0.4847634, 0.4950102, 0.5052186, 0.5155007, 0.5259549, 0.5362744, 0.5469376, 0.5577136, 0.5685905, 0.5796047, 0.5903783, 0.6014493, 0.6124342, 0.6240415, 0.6352995, 0.6467392, 0.6582003, 0.6699791, 0.6818415, 0.6935872, 0.7054129, 0.7172307, 0.7294462, 0.741425, 0.7535431, 0.7662849, 0.7784239, 0.7908459, 0.8033691, 0.8161577, 0.8289431, 0.8416416, 0.8540263, 0.8672161, 0.8801994, 0.8932345, 0.9062557, 0.9193608, 0.9327284, 0.9459846, 0.9592487, 0.9726576, 0.9861496, 0.9996973, 1.012965, 1.026556, 1.040384, 1.054103, 1.067545, 1.082002, 1.095702, 1.109572, 1.123644, 1.137535, 1.151644, 1.165935, 1.180006, 1.193681, 1.208248, 1.222276, 1.236631, 1.251115, 1.265348, 1.279328, 1.293886, 1.3086, 1.322557, 1.337312, 1.351656, 1.366449, 1.380387, 1.395255, 1.410018, 1.424254, 1.4389, 1.453595, 1.468078, 1.482693, 1.497038, 1.512216, 1.526099, 1.541275, 1.555984, 1.570639, 1.585133, 1.599731, 1.614572, 1.629264, 1.644083, 1.658696, 1.672676, 1.687363, 1.701977, 1.716221, 1.731279, 1.74591, 1.760681, 1.775817, 1.789731, 1.802027, 1.816948, 1.831684, 1.846949, 1.860713, 1.875714, 1.889853, 1.903857, 1.918073, 1.932817, 1.946704, 1.960967, 1.975597, 1.989767, 2.003751, 2.017718, 2.031349, 2.045713, 2.059654, 2.074199, 2.087567, 2.101471, 2.114826, 2.128818, 2.142786, 2.156726, 2.169774, 2.183343, 2.196759, 2.210352, 2.223661, 2.236527, 2.250551, 2.263937, 2.276332, 2.289922, 2.303646, 2.315648, 2.329866, 2.341981, 2.355692, 2.368056, 2.381629, 2.393545, 2.406439, 2.418353, 2.431615, 2.444101, 2.456677, 2.470097, 2.480602, 2.493399, 2.506444, 2.517951, 2.529809, 2.542466, 2.552953, 2.566042, 2.577603, 2.589625, 2.601557, 2.612584, 2.624716, 2.635865, 2.647416, 2.659774, 2.669729, 2.681178, 2.692593, 2.70239, 2.714687, 2.725653, 2.736549, 2.747655, 2.758374, 2.768947, 2.77929, 2.79072, 2.80052, 2.810775, 2.821476, 2.831656, 2.841758, 2.852864, 2.862108, 2.872267, 2.881746, 2.892089, 2.90139, 2.911329, 2.921784, 2.930586, 2.940144, 2.949574, 2.959357, 2.968272, 2.977809, 2.987029, 2.993128, 3.007866, 3.015358, 3.024142, 3.032702, 3.041631, 3.049428, 3.057545, 3.067541, 3.075788, 3.084157, 3.092863, 3.100929, 3.108762, 3.11718, 3.125805, 3.133149, 3.140825, 3.148971, 3.157339, 3.165391, 3.172471, 3.180657, 3.187313, 3.195678, 3.203544, 3.210948, 3.218013, 3.225333, 3.233093, 3.238977, 3.247579, 3.254842, 3.260955, 3.266517, 3.274042, 3.280998, 3.28829, 3.29541, 3.301134, 3.307274, 3.314615, 3.321437, 3.328838, 3.333784, 3.339936, 3.346304, 3.351605, 3.359429, 3.365483, 3.369279, 3.376687, 3.382037, 3.387905, 3.394324, 3.400054, 3.406178, 3.410533, 3.416404, 3.423079, 3.428949, 3.432757, 3.439289, 3.444226, 3.447854, 3.45435, 3.4599, 3.465244, 3.469057, 3.474176, 3.47922, 3.484871, 3.489499, 3.494661, 3.499311, 3.502459, 3.602461, 3.513154, 3.704515, 3.521801, 3.714805, 3.530733, 3.724777, 3.539739, 3.734745, 3.548497, 3.742732, 3.557354, 3.753331, 3.565609, 3.762425, 3.573428, 3.771394, 3.580571, 3.779998, 3.589278, 3.788848, 3.595776, 3.796774, 3.60622, 3.807195, 3.608768, 3.813079, 3.618011, 3.82036, 3.625034, 3.829281, 3.834229, 3.634407, 3.840508, 3.640856, 3.847535, 3.647933, 3.855823, 3.654437, 3.657347, 3.866458, 3.664482, 3.874061, 3.67016, 3.880671, 3.675033, 3.677998, 3.890826, 3.683444, 3.897618, 3.690037, 3.902379, 3.907571, 3.697705, 3.91411, 3.702186, 3.919342, 3.707471, 3.709338, 3.928181, 3.715962, 3.933149, 3.936695, 3.722031, 3.942656, 3.727679, 3.7293, 3.951641, 3.733918, 3.956491, 3.737203, 3.740092, 3.965215, 3.743594, 3.746006, 3.971083, 3.750773, 3.97721, 3.755033, 3.757032, 3.984024, 3.760458, 3.761486, 3.991435, 3.766177, 3.767856, 3.997193, 3.769815, 3.771316, 4.00569, 3.777636, 4.008779, 4.011422, 3.781285, 4.014348, 4.018114, 3.785748, 4.021489, 4.022959, 3.789414, 4.025624, 4.029503, 3.795098, 4.033351, 4.035794, 3.79941, 3.800271, 4.041397, 3.802813, 3.803955, 4.046184, 3.805146, 3.80697, 4.049605, 3.809575, 3.811374, 4.056009, 3.81409, 3.814975, 4.060414, 4.063359, 3.817723, 3.81969, 4.066606, 3.821201, 3.821128, 
1e-100, 0.0009647657, 0.001976792, 0.003052725, 0.004185581, 0.005368608, 0.006631035, 0.007946256, 0.009324261, 0.01078226, 0.01229701, 0.01388884, 0.01555978, 0.01729106, 0.01911046, 0.02101165, 0.02297684, 0.02504521, 0.02718745, 0.02941237, 0.03173825, 0.0341533, 0.03664172, 0.03923732, 0.04194079, 0.04471652, 0.04760418, 0.05060639, 0.05368843, 0.05687566, 0.06019675, 0.0636053, 0.06711231, 0.07074711, 0.07451019, 0.07834524, 0.08231772, 0.08644549, 0.09064937, 0.09495207, 0.09942467, 0.1040007, 0.1086679, 0.1135024, 0.118472, 0.1235156, 0.1287012, 0.1340378, 0.1395536, 0.1451434, 0.150883, 0.1567552, 0.1627302, 0.1689065, 0.175121, 0.1815748, 0.1880941, 0.1947974, 0.2015353, 0.2085504, 0.2156468, 0.2227811, 0.230119, 0.2375649, 0.2452932, 0.2529623, 0.2608147, 0.2688586, 0.2770021, 0.2852845, 0.2936728, 0.3021822, 0.3109361, 0.3198709, 0.3287735, 0.3378292, 0.3470938, 0.3565095, 0.3660698, 0.3755792, 0.3853865, 0.3953877, 0.405414, 0.4156209, 0.4259064, 0.4363973, 0.447067, 0.4577477, 0.4683379, 0.4795012, 0.4905059, 0.501932, 0.5131741, 0.5246669, 0.536215, 0.5482159, 0.5598857, 0.5719713, 0.5841031, 0.5960124, 0.6090695, 0.6212343, 0.6339188, 0.6462982, 0.6593301, 0.6722438, 0.6853971, 0.698482, 0.711812, 0.7253062, 0.7387936, 0.7523319, 0.7660667, 0.7798415, 0.7938265, 0.8075962, 0.8221947, 0.8367021, 0.8509303, 0.8652096, 0.8794427, 0.8943724, 0.9089802, 0.923927, 0.9387647, 0.953864, 0.969023, 0.9843254, 0.9997028, 1.015039, 1.03027, 1.045902, 1.061523, 1.07741, 1.093326, 1.108733, 1.12457, 1.140795, 1.156946, 1.17267, 1.188853, 1.205557, 1.22156, 1.237992, 1.254365, 1.271203, 1.287696, 1.304293, 1.320835, 1.337047, 1.354811, 1.371579, 1.387989, 1.405235, 1.422466, 1.439155, 1.456303, 1.473211, 1.490459, 1.507859, 1.525274, 1.541907, 1.558855, 1.576064, 1.594074, 1.611055, 1.628757, 1.645957, 1.663292, 1.68046, 1.698631, 1.715873, 1.733461, 1.751002, 1.76802, 1.785577, 1.802539, 1.820652, 1.838137, 1.855182, 1.873188, 1.891093, 1.908052, 1.925863, 1.943641, 1.961378, 1.978073, 1.99623, 2.013591, 2.030406, 2.048029, 2.066074, 2.08276, 2.100905, 2.118464, 2.135447, 2.153133, 2.170513, 2.188882, 2.204429, 2.219756, 2.238254, 2.255555, 2.273136, 2.290038, 2.307668, 2.324101, 2.340898, 2.35756, 2.375432, 2.391789, 2.408337, 2.426344, 2.442054, 2.45876, 2.475546, 2.492645, 2.50854, 2.525895, 2.542252, 2.557989, 2.574439, 2.591524, 2.606843, 2.62387, 2.640112, 2.655681, 2.671559, 2.687793, 2.703339, 2.719068, 2.735032, 2.750544, 2.766217, 2.782859, 2.79765, 2.813437, 2.828948, 2.844338, 2.858554, 2.873848, 2.889656, 2.904607, 2.918834, 2.93398, 2.949404, 2.962941, 2.978914, 2.992938, 3.009171, 3.022772, 3.037173, 3.051117, 3.065327, 3.079423, 3.09357, 3.107077, 3.121823, 3.135352, 3.147806, 3.163026, 3.176661, 3.189498, 3.203595, 3.216693, 3.23052, 3.243457, 3.256251, 3.269646, 3.282701, 3.295495, 3.307557, 3.319652, 3.334409, 3.346181, 3.358471, 3.371468, 3.383468, 3.394919, 3.407749, 3.419676, 3.431861, 3.443652, 3.455253, 3.466332, 3.478402, 3.490858, 3.501747, 3.512357, 3.524249, 3.535697, 3.546917, 3.557152, 3.56909, 3.578965, 3.590556, 3.600706, 3.611792, 3.620243, 3.636321, 3.644501, 3.6531, 3.664677, 3.674273, 3.684499, 3.694709, 3.70487, 3.713889, 3.723423, 3.734457, 3.743763, 3.753039, 3.762738, 3.771588, 3.780983, 3.789635, 3.800141, 3.807852, 3.817514, 3.826925, 3.836122, 3.843823, 3.85305, 3.861637, 3.869641, 3.878948, 3.887082, 3.896787, 3.903806, 3.911967, 3.919196, 3.928683, 3.935306, 3.943895, 3.950101, 3.959477, 3.966681, 3.974956, 3.982235, 3.990375, 3.997737, 4.003977, 4.012618, 4.018798, 4.026687, 4.033107, 4.039818, 4.046001, 4.053002, 4.061237, 4.067366, 4.074094, 4.080535, 4.087194, 4.095162, 4.099948, 4.10714, 4.113047, 4.119357, 4.125126, 4.132097, 4.138162, 4.142928, 4.149058, 4.154488, 4.161402, 4.167162, 4.173276, 4.178731, 4.182852, 4.190429, 4.195129, 4.200624, 4.205719, 4.326015, 4.216116, 4.451963, 4.227222, 4.46366, 4.237004, 4.47365, 4.246208, 4.486773, 4.258036, 4.497555, 4.266472, 4.507825, 4.2752, 4.5191, 4.28489, 4.528393, 4.293178, 4.538644, 4.302021, 4.548237, 4.311437, 4.558557, 4.322129, 4.566363, 4.32628, 4.57611, 4.33472, 4.585861, 4.343425, 4.59406, 4.59936, 4.354984, 4.607786, 4.361876, 4.615507, 4.369583, 4.625863, 4.37764, 4.380042, 4.638188, 4.387826, 4.644448, 4.393587, 4.654001, 4.400353, 4.402261, 4.664094, 4.409472, 4.673423, 4.417231, 4.679742, 4.684016, 4.425832, 4.689623, 4.431943, 4.697703, 4.437051, 4.438894, 4.707651, 4.443548, 4.714176, 4.718513, 4.453455, 4.72401, 4.459258, 4.460162, 4.734345, 4.465707, 4.740317, 4.470331, 4.473305, 4.748049, 4.477563, 4.480738, 4.758487, 4.485585, 4.762869, 4.489434, 4.492965, 4.772839, 4.495385, 4.497617, 4.779739, 4.501407, 4.502922, 4.786946, 4.508243, 4.509842, 4.795099, 4.513597, 4.800235, 4.802613, 4.51901, 4.807825, 4.810407, 4.523481, 4.812958, 4.816872, 4.52894, 4.820604, 4.824409, 4.533023, 4.828285, 4.830937, 4.539755, 4.540436, 4.836641, 4.542585, 4.543071, 4.842119, 4.547345, 4.547681, 4.848774, 4.550341, 4.553096, 4.855744, 4.556648, 4.556172, 4.859638, 4.861847, 4.559271, 4.560018, 4.865874, 4.564583, 4.56426, 
1e-100, 0.001306942, 0.002707847, 0.004170931, 0.005711993, 0.007339061, 0.00903327, 0.01082652, 0.01270157, 0.01464954, 0.01671426, 0.0188606, 0.02108872, 0.02344068, 0.02587071, 0.02841206, 0.03106323, 0.03381043, 0.03667019, 0.03965032, 0.04273651, 0.04593222, 0.04926689, 0.05270765, 0.05625101, 0.05996986, 0.0637833, 0.06771081, 0.07181446, 0.07601649, 0.08035109, 0.08484785, 0.08947684, 0.09423623, 0.09912603, 0.1042003, 0.109403, 0.1146999, 0.1202325, 0.1259016, 0.1316562, 0.1376043, 0.14371, 0.149985, 0.1563236, 0.1630049, 0.1697092, 0.176587, 0.1836269, 0.1908976, 0.1983252, 0.2058545, 0.2135326, 0.2214397, 0.2295353, 0.2376573, 0.2460888, 0.2546319, 0.2633414, 0.2721608, 0.2812426, 0.2904214, 0.2997965, 0.3093164, 0.318975, 0.3288995, 0.3389724, 0.3490167, 0.3593713, 0.3699471, 0.3805355, 0.3913276, 0.402427, 0.4135552, 0.4249206, 0.4364601, 0.4480071, 0.4597817, 0.4718524, 0.4839939, 0.4963559, 0.5086673, 0.521297, 0.534107, 0.5470327, 0.5601241, 0.5732037, 0.5867458, 0.6003406, 0.6140584, 0.6278135, 0.6419804, 0.6559951, 0.670534, 0.6848491, 0.6994041, 0.7140084, 0.7290787, 0.7441239, 0.7590172, 0.7749544, 0.7902704, 0.8060463, 0.8216491, 0.8373257, 0.8531537, 0.8695651, 0.886146, 0.9024899, 0.9186443, 0.93504, 0.9521505, 0.969104, 0.9861972, 1.003225, 1.020481, 1.037887, 1.055408, 1.072806, 1.090721, 1.108316, 1.126369, 1.144244, 1.162513, 1.181074, 1.199508, 1.217536, 1.236448, 1.25525, 1.27382, 1.292575, 1.311516, 1.330696, 1.349209, 1.368311, 1.387482, 1.40715, 1.426452, 1.44579, 1.465767, 1.4852, 1.50491, 1.525199, 1.544438, 1.564776, 1.584771, 1.604036, 1.624702, 1.645225, 1.664956, 1.685525, 1.705713, 1.725985, 1.746752, 1.7673, 1.787365, 1.807326, 1.828302, 1.848691, 1.869465, 1.890055, 1.910773, 1.931766, 1.952401, 1.973152, 1.99367, 2.014609, 2.036024, 2.056812, 2.077528, 2.098069, 2.119304, 2.139815, 2.161023, 2.181098, 2.202498, 2.223051, 2.244218, 2.265227, 2.28604, 2.306981, 2.328035, 2.348913, 2.370124, 2.391042, 2.411795, 2.431955, 2.452523, 2.473818, 2.493563, 2.515607, 2.535693, 2.556694, 2.576362, 2.598335, 2.618598, 2.639225, 2.659883, 2.682232, 2.700046, 2.718044, 2.740286, 2.760646, 2.780829, 2.801787, 2.821768, 2.841349, 2.861592, 2.88195, 2.901835, 2.92135, 2.94193, 2.961805, 2.980715, 3.001454, 3.020082, 3.040341, 3.059986, 3.078941, 3.098871, 3.116674, 3.136843, 3.156173, 3.174831, 3.194621, 3.213408, 3.231368, 3.250606, 3.268879, 3.287498, 3.306224, 3.325432, 3.342425, 3.361, 3.380185, 3.398059, 3.415227, 3.434723, 3.451911, 3.469076, 3.487867, 3.504819, 3.522553, 3.540493, 3.557503, 3.574845, 3.592187, 3.609587, 3.626121, 3.643957, 3.660561, 3.676943, 3.692927, 3.709891, 3.727271, 3.742302, 3.759137, 3.775042, 3.791634, 3.807572, 3.823485, 3.8389, 3.854662, 3.871362, 3.886141, 3.902074, 3.915928, 3.931557, 3.945681, 3.962047, 3.975639, 3.991554, 4.006504, 4.021291, 4.035642, 4.049314, 4.064974, 4.078039, 4.093385, 4.106864, 4.120007, 4.133469, 4.149236, 4.162053, 4.17502, 4.188495, 4.201457, 4.21575, 4.227977, 4.24169, 4.254382, 4.266877, 4.278555, 4.292931, 4.305417, 4.317695, 4.330069, 4.340226, 4.352386, 4.369697, 4.379252, 4.391014, 4.403183, 4.413485, 4.425057, 4.437835, 4.448786, 4.459573, 4.471061, 4.481724, 4.493635, 4.504203, 4.515385, 4.525396, 4.535702, 4.54505, 4.557922, 4.568109, 4.578305, 4.588631, 4.597143, 4.610112, 4.618729, 4.627715, 4.638457, 4.647732, 4.65693, 4.666538, 4.676687, 4.685739, 4.693652, 4.704394, 4.710538, 4.722083, 4.73111, 4.740199, 4.748438, 4.755754, 4.765605, 4.774457, 4.783342, 4.79153, 4.799219, 4.806517, 4.815338, 4.823714, 4.831738, 4.837617, 4.845465, 4.853702, 4.862334, 4.870189, 4.877945, 4.886008, 4.891524, 4.900232, 4.907828, 4.91366, 4.922116, 4.928185, 4.934736, 4.941878, 4.948397, 4.95599, 4.962298, 4.96823, 4.974534, 4.982302, 4.989098, 4.995018, 5.000725, 5.006252, 5.012584, 5.01977, 5.16228, 5.029961, 5.313435, 5.042478, 5.325664, 5.054028, 5.340103, 5.065905, 5.352638, 5.075007, 5.366113, 5.086656, 5.37731, 5.098415, 5.38907, 5.106828, 5.400754, 5.117872, 5.411735, 5.127727, 5.423174, 5.137218, 5.434137, 5.149912, 5.443674, 5.152825, 5.452617, 5.164512, 5.464641, 5.174897, 5.472548, 5.479648, 5.185828, 5.489287, 5.195453, 5.499241, 5.203983, 5.50872, 5.209758, 5.216186, 5.524689, 5.223067, 5.531546, 5.230635, 5.539699, 5.23652, 5.240143, 5.556016, 5.248581, 5.561808, 5.255546, 5.569478, 5.576904, 5.265494, 5.582618, 5.27164, 5.588399, 5.277752, 5.281425, 5.601367, 5.287542, 5.607319, 5.615485, 5.297742, 5.620785, 5.302052, 5.3044, 5.633082, 5.309127, 5.636287, 5.316099, 5.318466, 5.647648, 5.322533, 5.326823, 5.661806, 5.331918, 5.663793, 5.336801, 5.33961, 5.673321, 5.343054, 5.345003, 5.682611, 5.348766, 5.350755, 5.690459, 5.3571, 5.360835, 5.703135, 5.363826, 5.705319, 5.710297, 5.369417, 5.714121, 5.719776, 5.373277, 5.719588, 5.725352, 5.379828, 5.729597, 5.73585, 5.38621, 5.736759, 5.741477, 5.389552, 5.392979, 5.748897, 5.395917, 5.396049, 5.752399, 5.397329, 5.401432, 5.761026, 5.405683, 5.406998, 5.766899, 5.409472, 5.411518, 5.772394, 5.77867, 5.41352, 5.413537, 5.778956, 5.416192, 5.419721, 
1e-100, 0.001788998, 0.003659265, 0.005648736, 0.00772599, 0.009901184, 0.01220112, 0.01459117, 0.01710272, 0.01973504, 0.02246498, 0.02534125, 0.0283279, 0.03141866, 0.03467873, 0.03805576, 0.04154009, 0.04521561, 0.04900219, 0.05289954, 0.05701123, 0.06123487, 0.06557443, 0.07014965, 0.07483669, 0.07966671, 0.08467183, 0.08986222, 0.09517423, 0.1006906, 0.106406, 0.112224, 0.1182115, 0.1244717, 0.1308757, 0.1373805, 0.1441435, 0.151082, 0.1581477, 0.1654313, 0.1729811, 0.180689, 0.1885027, 0.196611, 0.2049042, 0.2133136, 0.2219455, 0.2308365, 0.2398351, 0.2491342, 0.2586138, 0.2683026, 0.2781546, 0.288143, 0.2984695, 0.30896, 0.3195411, 0.3304012, 0.3414304, 0.352854, 0.3643019, 0.3758141, 0.3877563, 0.3997917, 0.4120974, 0.4245686, 0.4371778, 0.4500169, 0.463025, 0.4762718, 0.4896172, 0.503262, 0.5170311, 0.5311562, 0.5452848, 0.5595899, 0.5743242, 0.5892045, 0.6040135, 0.6191361, 0.6345724, 0.650179, 0.6659414, 0.6818031, 0.6978848, 0.7143076, 0.7306369, 0.7472207, 0.7639391, 0.7811289, 0.7983879, 0.8158201, 0.8331769, 0.8510055, 0.8689365, 0.8869566, 0.9049656, 0.9234011, 0.9415753, 0.9603852, 0.9804296, 0.998869, 1.017944, 1.036869, 1.056558, 1.076244, 1.096238, 1.116135, 1.135917, 1.156234, 1.17667, 1.197258, 1.217575, 1.238537, 1.259502, 1.280837, 1.301637, 1.322785, 1.34473, 1.366148, 1.387939, 1.409491, 1.431739, 1.453973, 1.475542, 1.498051, 1.520309, 1.543156, 1.565647, 1.587838, 1.611147, 1.633561, 1.656732, 1.679172, 1.703111, 1.726076, 1.749567, 1.772723, 1.796563, 1.820118, 1.844034, 1.867204, 1.891699, 1.9153, 1.938782, 1.963024, 1.987073, 2.01128, 2.035276, 2.05872, 2.083642, 2.108025, 2.132793, 2.157018, 2.181109, 2.205655, 2.230244, 2.254833, 2.27931, 2.304357, 2.328455, 2.353521, 2.378203, 2.403095, 2.427489, 2.452682, 2.478234, 2.501952, 2.526832, 2.551224, 2.576265, 2.601053, 2.626083, 2.65009, 2.674774, 2.699986, 2.724815, 2.74891, 2.774783, 2.79934, 2.824223, 2.848182, 2.873493, 2.897692, 2.922469, 2.946895, 2.971594, 2.995454, 3.020779, 3.045092, 3.069893, 3.093936, 3.118942, 3.142981, 3.167319, 3.19175, 3.21505, 3.239887, 3.265642, 3.287258, 3.308874, 3.334318, 3.357927, 3.382114, 3.40578, 3.429631, 3.452567, 3.476594, 3.500499, 3.523945, 3.546308, 3.570737, 3.592629, 3.616338, 3.639191, 3.662658, 3.683999, 3.70906, 3.730198, 3.753326, 3.77517, 3.798372, 3.820196, 3.842576, 3.865325, 3.885693, 3.908331, 3.930831, 3.951746, 3.973942, 3.995645, 4.016866, 4.036538, 4.060632, 4.080691, 4.10153, 4.122616, 4.142923, 4.163957, 4.183786, 4.206436, 4.225706, 4.245603, 4.266519, 4.286793, 4.305005, 4.326759, 4.345841, 4.365249, 4.386957, 4.405186, 4.423911, 4.443949, 4.462501, 4.481592, 4.500341, 4.51796, 4.537418, 4.554507, 4.57531, 4.593638, 4.610673, 4.629348, 4.647474, 4.664136, 4.685176, 4.699346, 4.717584, 4.735033, 4.752268, 4.768349, 4.785588, 4.803912, 4.820092, 4.836791, 4.852853, 4.869395, 4.886274, 4.902408, 4.918571, 4.934484, 4.949251, 4.965552, 4.980142, 4.997119, 5.012708, 5.026505, 5.041314, 5.056532, 5.072004, 5.086803, 5.101285, 5.116583, 5.129293, 5.144627, 5.158313, 5.173161, 5.187246, 5.197404, 5.21804, 5.227762, 5.242566, 5.25534, 5.268718, 5.281436, 5.294882, 5.307405, 5.319999, 5.334308, 5.345266, 5.357976, 5.370525, 5.383007, 5.394694, 5.406963, 5.419226, 5.430443, 5.443026, 5.454889, 5.466173, 5.477484, 5.488372, 5.500177, 5.510914, 5.522485, 5.533843, 5.544621, 5.555024, 5.565335, 5.575066, 5.586946, 5.597351, 5.60613, 5.614842, 5.627311, 5.637925, 5.64739, 5.65727, 5.665974, 5.676312, 5.684688, 5.695932, 5.7039, 5.713421, 5.723111, 5.730062, 5.73914, 5.748174, 5.756745, 5.76704, 5.774206, 5.783725, 5.792892, 5.802538, 5.808803, 5.816722, 5.825274, 5.831899, 5.842687, 5.849631, 5.855449, 5.862719, 5.870897, 5.877937, 5.886707, 5.895453, 5.901475, 5.909175, 5.914417, 5.923509, 5.929649, 5.936586, 5.943923, 5.949654, 6.119262, 5.962696, 6.299172, 5.975624, 6.313417, 5.988917, 6.328065, 6.00226, 6.343092, 6.014165, 6.357123, 6.024991, 6.370725, 6.038852, 6.38382, 6.048214, 6.396772, 6.05894, 6.40831, 6.07147, 6.422687, 6.083144, 6.435012, 6.09394, 6.44659, 6.100104, 6.457539, 6.112065, 6.468558, 6.121519, 6.479415, 6.485796, 6.135833, 6.498005, 6.144794, 6.508495, 6.15575, 6.520121, 6.164276, 6.167876, 6.534779, 6.175878, 6.544917, 6.184888, 6.553321, 6.190831, 6.197357, 6.57012, 6.206396, 6.579026, 6.212595, 6.587963, 6.592884, 6.222207, 6.601712, 6.231498, 6.609696, 6.236166, 6.239344, 6.623067, 6.249507, 6.632153, 6.636989, 6.257926, 6.64407, 6.264998, 6.267653, 6.657825, 6.272884, 6.662933, 6.277494, 6.28097, 6.675394, 6.28795, 6.291319, 6.687423, 6.294402, 6.69593, 6.302713, 6.303371, 6.703356, 6.309388, 6.311177, 6.712845, 6.316329, 6.319793, 6.725253, 6.323948, 6.325023, 6.734899, 6.330927, 6.739041, 6.744675, 6.336568, 6.747984, 6.752084, 6.342973, 6.756692, 6.761486, 6.350181, 6.765531, 6.77032, 6.35666, 6.775154, 6.777673, 6.361205, 6.362202, 6.783939, 6.365577, 6.36822, 6.792712, 6.371081, 6.372823, 6.800125, 6.378478, 6.379152, 6.807429, 6.380392, 6.382382, 6.813904, 6.816626, 6.386213, 6.387161, 6.824908, 6.391376, 6.392075, 
1e-100, 0.002387631, 0.004921496, 0.00756051, 0.0103457, 0.01326502, 0.01629655, 0.01950985, 0.02284079, 0.02631884, 0.02997422, 0.03375813, 0.03770896, 0.04183343, 0.04610006, 0.05055828, 0.05518812, 0.05996602, 0.06495161, 0.07012949, 0.07545231, 0.0810034, 0.08674294, 0.09264724, 0.09875789, 0.1051076, 0.1116086, 0.1183331, 0.1253082, 0.1324423, 0.1397608, 0.1473779, 0.1551767, 0.1631468, 0.1714557, 0.179905, 0.1885917, 0.1975019, 0.2066793, 0.2160597, 0.2256333, 0.2355818, 0.2456515, 0.2559655, 0.2664772, 0.2773579, 0.2884461, 0.2996784, 0.3111867, 0.3229669, 0.3351319, 0.3472389, 0.3597748, 0.3726837, 0.3856049, 0.3988295, 0.4122724, 0.4260013, 0.4399574, 0.4540821, 0.4686062, 0.4832701, 0.4980921, 0.5132902, 0.5286911, 0.5444135, 0.560298, 0.5762599, 0.5926354, 0.6092345, 0.6259365, 0.6428248, 0.6601347, 0.6776781, 0.6953616, 0.7133852, 0.7312467, 0.7498105, 0.768409, 0.7872874, 0.8062796, 0.825364, 0.8448667, 0.8647654, 0.8845124, 0.9046495, 0.9248846, 0.9455793, 0.9662464, 0.9872902, 1.008196, 1.029525, 1.051085, 1.072869, 1.094733, 1.116732, 1.139053, 1.161909, 1.184205, 1.206796, 1.230989, 1.253763, 1.277103, 1.300885, 1.324714, 1.348408, 1.372614, 1.396478, 1.421384, 1.445935, 1.470624, 1.49536, 1.520607, 1.546035, 1.57095, 1.596382, 1.62185, 1.647875, 1.674051, 1.699931, 1.726112, 1.752698, 1.778825, 1.805779, 1.832201, 1.859283, 1.886042, 1.912648, 1.940181, 1.967486, 1.995055, 2.022255, 2.049905, 2.076889, 2.105282, 2.133007, 2.161159, 2.189662, 2.217124, 2.24488, 2.273818, 2.302563, 2.330359, 2.359046, 2.38709, 2.41619, 2.444724, 2.47323, 2.50154, 2.52988, 2.560162, 2.58827, 2.61735, 2.646211, 2.675846, 2.704703, 2.733689, 2.762851, 2.79153, 2.822049, 2.850784, 2.879488, 2.909113, 2.937904, 2.96743, 2.996058, 3.027018, 3.05504, 3.083792, 3.113161, 3.142682, 3.170668, 3.200633, 3.230285, 3.258814, 3.288295, 3.317656, 3.346549, 3.376006, 3.404683, 3.433758, 3.462469, 3.491447, 3.520428, 3.54871, 3.577664, 3.606586, 3.636108, 3.663524, 3.693631, 3.722142, 3.749853, 3.777923, 3.806903, 3.835045, 3.863654, 3.892115, 3.920973, 3.950222, 3.97533, 3.999719, 4.02964, 4.057031, 4.086113, 4.11283, 4.141302, 4.1685, 4.195977, 4.222503, 4.249653, 4.276374, 4.303419, 4.331449, 4.357589, 4.384063, 4.411216, 4.437603, 4.463144, 4.489144, 4.515229, 4.541647, 4.567844, 4.593575, 4.619159, 4.644982, 4.668693, 4.695283, 4.720024, 4.744661, 4.770042, 4.79515, 4.819593, 4.842809, 4.869414, 4.892639, 4.917276, 4.942014, 4.96492, 4.989081, 5.012206, 5.03647, 5.058676, 5.083844, 5.106147, 5.12898, 5.152102, 5.174648, 5.196832, 5.219042, 5.24526, 5.264871, 5.286125, 5.309032, 5.330482, 5.35085, 5.37512, 5.394304, 5.415565, 5.436948, 5.459085, 5.479282, 5.499853, 5.521227, 5.541442, 5.5617, 5.582774, 5.600244, 5.620276, 5.641247, 5.659584, 5.679694, 5.699503, 5.718722, 5.737178, 5.757153, 5.776364, 5.792604, 5.812698, 5.830354, 5.847797, 5.866208, 5.885508, 5.902781, 5.919312, 5.937505, 5.954424, 5.971552, 5.98912, 6.005785, 6.022552, 6.039576, 6.056079, 6.07351, 6.088838, 6.104919, 6.120643, 6.134356, 6.152834, 6.164431, 6.18675, 6.199546, 6.215699, 6.228435, 6.242987, 6.259554, 6.273775, 6.287264, 6.302721, 6.316542, 6.331488, 6.344664, 6.359435, 6.372412, 6.386491, 6.399581, 6.413844, 6.426622, 6.439436, 6.453368, 6.464141, 6.480609, 6.49108, 6.503522, 6.516308, 6.528671, 6.540002, 6.552084, 6.56666, 6.577045, 6.587531, 6.600775, 6.611868, 6.622664, 6.635652, 6.645823, 6.657218, 6.666329, 6.679122, 6.689736, 6.700701, 6.710623, 6.720847, 6.729133, 6.74156, 6.752365, 6.760262, 6.769611, 6.781677, 6.789469, 6.803255, 6.810637, 6.820766, 6.829916, 6.837395, 6.849002, 6.856502, 6.866415, 6.874345, 6.884422, 6.890671, 6.900111, 6.9081, 6.916841, 6.926041, 6.934116, 6.943708, 6.9508, 6.958842, 6.966202, 6.974684, 6.979992, 6.989567, 6.997992, 7.003475, 7.201564, 7.018163, 7.409681, 7.033878, 7.428643, 7.048758, 7.445054, 7.06021, 7.459681, 7.076527, 7.475643, 7.088892, 7.490704, 7.100617, 7.50403, 7.113725, 7.519955, 7.124558, 7.533148, 7.137699, 7.548905, 7.151886, 7.562628, 7.163216, 7.574348, 7.168069, 7.585393, 7.182798, 7.599174, 7.190817, 7.610279, 7.618365, 7.210322, 7.633461, 7.219943, 7.64387, 7.229196, 7.654496, 7.239507, 7.243451, 7.673605, 7.252695, 7.683801, 7.260892, 7.693794, 7.271422, 7.275614, 7.711452, 7.282694, 7.722355, 7.29394, 7.731553, 7.737903, 7.302834, 7.745024, 7.310852, 7.756725, 7.318918, 7.323522, 7.770865, 7.328775, 7.781493, 7.787319, 7.341508, 7.793232, 7.348251, 7.351253, 7.805806, 7.357666, 7.816308, 7.364239, 7.368192, 7.828733, 7.374973, 7.378281, 7.841342, 7.383629, 7.848745, 7.389617, 7.391161, 7.858826, 7.397753, 7.398277, 7.870331, 7.405206, 7.407905, 7.883632, 7.414583, 7.417497, 7.894095, 7.420928, 7.897714, 7.903785, 7.428165, 7.907743, 7.912356, 7.432641, 7.918135, 7.925257, 7.443966, 7.928761, 7.933045, 7.447791, 7.93692, 7.94055, 7.453191, 7.456237, 7.950295, 7.457321, 7.459288, 7.95874, 7.465453, 7.46745, 7.967945, 7.47013, 7.472327, 7.974897, 7.47518, 7.476132, 7.981715, 7.983511, 7.478448, 7.481204, 7.99413, 7.488192, 7.488935, 
1e-100, 0.003198574, 0.006537884, 0.01008404, 0.01376646, 0.01763006, 0.02168579, 0.0258873, 0.03031703, 0.03491834, 0.03968247, 0.04471177, 0.04990652, 0.05528658, 0.06093476, 0.06675262, 0.07277501, 0.07908999, 0.08556875, 0.09227285, 0.09927525, 0.1064724, 0.1138924, 0.1216318, 0.1295636, 0.1377052, 0.1461789, 0.1549042, 0.16381, 0.1730639, 0.1826067, 0.1922873, 0.2023355, 0.2126654, 0.2231902, 0.2339901, 0.2452499, 0.256677, 0.268302, 0.2803269, 0.2926479, 0.3050766, 0.3179139, 0.331157, 0.3445707, 0.3582053, 0.3721935, 0.386565, 0.4012, 0.4159869, 0.431268, 0.4467559, 0.4624566, 0.4782883, 0.4948821, 0.5115134, 0.5282913, 0.5454708, 0.5630121, 0.5809817, 0.5989602, 0.6171959, 0.6358375, 0.6548273, 0.6739194, 0.6932973, 0.7132783, 0.7330177, 0.7534028, 0.773805, 0.7946679, 0.815841, 0.8371674, 0.858797, 0.8804812, 0.9025364, 0.9251686, 0.9478901, 0.970651, 0.993598, 1.017424, 1.041134, 1.064634, 1.089145, 1.113285, 1.138203, 1.162878, 1.187862, 1.213189, 1.238824, 1.264994, 1.290674, 1.316992, 1.343574, 1.370413, 1.397435, 1.424476, 1.451707, 1.479092, 1.506501, 1.535598, 1.563688, 1.591651, 1.62007, 1.649026, 1.677574, 1.707042, 1.735963, 1.765878, 1.794973, 1.825253, 1.855489, 1.885387, 1.91524, 1.945579, 1.977017, 2.007721, 2.038887, 2.06956, 2.100711, 2.132609, 2.163931, 2.195309, 2.227046, 2.259003, 2.290921, 2.322745, 2.355805, 2.388199, 2.420624, 2.452636, 2.486223, 2.519448, 2.55087, 2.584486, 2.617545, 2.651107, 2.684173, 2.717104, 2.750669, 2.784693, 2.818302, 2.851451, 2.884526, 2.919061, 2.952716, 2.986052, 3.020083, 3.053794, 3.087332, 3.122451, 3.156692, 3.190712, 3.225057, 3.259077, 3.293159, 3.326975, 3.36207, 3.395991, 3.429648, 3.464474, 3.497875, 3.532644, 3.5665, 3.601071, 3.634893, 3.669959, 3.703389, 3.737019, 3.771804, 3.806279, 3.84013, 3.87356, 3.90882, 3.941643, 3.976457, 4.009558, 4.044305, 4.078248, 4.111002, 4.145195, 4.178868, 4.211564, 4.246579, 4.279636, 4.313892, 4.345944, 4.380059, 4.413755, 4.445833, 4.479205, 4.512473, 4.544541, 4.578036, 4.610707, 4.64278, 4.676435, 4.709197, 4.743827, 4.772416, 4.80177, 4.834384, 4.867208, 4.898932, 4.93129, 4.961352, 4.994625, 5.025996, 5.056564, 5.087566, 5.119078, 5.149075, 5.181193, 5.211766, 5.241995, 5.271569, 5.303109, 5.333811, 5.362191, 5.393448, 5.422288, 5.452438, 5.481643, 5.509226, 5.540656, 5.569697, 5.597515, 5.627322, 5.654211, 5.68297, 5.711439, 5.740021, 5.767181, 5.796548, 5.823469, 5.851029, 5.877868, 5.904614, 5.932186, 5.958524, 5.986802, 6.0122, 6.039569, 6.065074, 6.090552, 6.117859, 6.142647, 6.167945, 6.193991, 6.22217, 6.244678, 6.271078, 6.294991, 6.318989, 6.343853, 6.368133, 6.391446, 6.415238, 6.440466, 6.464145, 6.486871, 6.510837, 6.533775, 6.556912, 6.58114, 6.603722, 6.623198, 6.646496, 6.669684, 6.690741, 6.711408, 6.735502, 6.756456, 6.777447, 6.79914, 6.81952, 6.840931, 6.861492, 6.882159, 6.902064, 6.922316, 6.941746, 6.963201, 6.982275, 7.002749, 7.020625, 7.039753, 7.059443, 7.078991, 7.097081, 7.116184, 7.135311, 7.150848, 7.171988, 7.188418, 7.206269, 7.224254, 7.240957, 7.254602, 7.280374, 7.294883, 7.31217, 7.326763, 7.344146, 7.359931, 7.377781, 7.39292, 7.408974, 7.423919, 7.440226, 7.456352, 7.47263, 7.487131, 7.50296, 7.518076, 7.529667, 7.548318, 7.560574, 7.575239, 7.590806, 7.604201, 7.617399, 7.631221, 7.645563, 7.659329, 7.671695, 7.6855, 7.698207, 7.712777, 7.727105, 7.739157, 7.752342, 7.759794, 7.778653, 7.789525, 7.800774, 7.813062, 7.825231, 7.835786, 7.847422, 7.861164, 7.872733, 7.8817, 7.894546, 7.903331, 7.915548, 7.928522, 7.938634, 7.947323, 7.959078, 7.970994, 7.981211, 7.991959, 8.001919, 8.011028, 8.019396, 8.031269, 8.041069, 8.050743, 8.056766, 8.067847, 8.07697, 8.088118, 8.097487, 8.105979, 8.116637, 8.122923, 8.132614, 8.142086, 8.148921, 8.157878, 8.166312, 8.173479, 8.180785, 8.411872, 8.197749, 8.652102, 8.213059, 8.668545, 8.230522, 8.688364, 8.244826, 8.705467, 8.258616, 8.720586, 8.273672, 8.739347, 8.287387, 8.7523, 8.301088, 8.769636, 8.315944, 8.786602, 8.329445, 8.802386, 8.339978, 8.817075, 8.355689, 8.830506, 8.36309, 8.844092, 8.375273, 8.856677, 8.386937, 8.87228, 8.879921, 8.40513, 8.893617, 8.418262, 8.907448, 8.427339, 8.918108, 8.436874, 8.440573, 8.938075, 8.453584, 8.949505, 8.462269, 8.959688, 8.472103, 8.478287, 8.981598, 8.486241, 8.991853, 8.495578, 9.000257, 9.008193, 8.509846, 9.018414, 8.514624, 9.027299, 8.52422, 8.529779, 9.046472, 8.53833, 9.055285, 9.061249, 8.548739, 9.070624, 8.555184, 8.559117, 9.083734, 8.56472, 9.091445, 8.573203, 8.576433, 9.108259, 8.584476, 8.586849, 9.123174, 8.595782, 9.129768, 8.599558, 8.600857, 9.142541, 8.60595, 8.610625, 9.156174, 8.618886, 8.620096, 9.166713, 8.624412, 8.629099, 9.178452, 8.633939, 9.18493, 9.190801, 8.639828, 9.195787, 9.200234, 8.648182, 9.208272, 9.213246, 8.655165, 9.218806, 9.223862, 8.662868, 9.227684, 9.229793, 8.666335, 8.669856, 9.242726, 8.67329, 8.676121, 9.25146, 8.680261, 8.683371, 9.261638, 8.686849, 8.687473, 9.268445, 8.690103, 8.691407, 9.276287, 9.279983, 8.699017, 8.699391, 9.290092, 8.701485, 8.705871, 
1e-100, 0.004237335, 0.008695223, 0.01332649, 0.01822614, 0.02331763, 0.02862295, 0.03420322, 0.03997557, 0.04601386, 0.05231792, 0.05881603, 0.06563407, 0.07268818, 0.07997623, 0.08758545, 0.0954773, 0.1035614, 0.1120561, 0.1207825, 0.1297413, 0.1390995, 0.1487083, 0.1585774, 0.1688441, 0.1793963, 0.1902161, 0.2014244, 0.2129005, 0.2246766, 0.2368462, 0.2493805, 0.2621856, 0.2752106, 0.2887732, 0.3026337, 0.3167038, 0.3312673, 0.3461164, 0.3612595, 0.3766954, 0.3926802, 0.4089142, 0.4255041, 0.4424002, 0.4596304, 0.4772233, 0.4950985, 0.5135058, 0.5321892, 0.5511795, 0.5702809, 0.5902892, 0.6104055, 0.6306603, 0.6512954, 0.6724707, 0.6938904, 0.7155582, 0.7375613, 0.7601224, 0.7826692, 0.8056772, 0.8290521, 0.8531189, 0.8770644, 0.9014704, 0.9257188, 0.9509128, 0.9763183, 1.001748, 1.027444, 1.053857, 1.080456, 1.107309, 1.134119, 1.161503, 1.18957, 1.217634, 1.245681, 1.27438, 1.302989, 1.33199, 1.361864, 1.391318, 1.421154, 1.451557, 1.481718, 1.512719, 1.543395, 1.57487, 1.60609, 1.638091, 1.669724, 1.702395, 1.734491, 1.767463, 1.800386, 1.833663, 1.865822, 1.901435, 1.93488, 1.96861, 2.002931, 2.037294, 2.07183, 2.107312, 2.142094, 2.177606, 2.212598, 2.248469, 2.284368, 2.320178, 2.356563, 2.392421, 2.429134, 2.464909, 2.502945, 2.540354, 2.577092, 2.613728, 2.65027, 2.688602, 2.726378, 2.763984, 2.802168, 2.839946, 2.878495, 2.916879, 2.955375, 2.993945, 3.032177, 3.071166, 3.108639, 3.148861, 3.188078, 3.226493, 3.265024, 3.304461, 3.343706, 3.382942, 3.422048, 3.462437, 3.500798, 3.541077, 3.580108, 3.619934, 3.659833, 3.699616, 3.738123, 3.779047, 3.819017, 3.858905, 3.897874, 3.938525, 3.978374, 4.018204, 4.058221, 4.097646, 4.137961, 4.177059, 4.217206, 4.25631, 4.295373, 4.334964, 4.375421, 4.416381, 4.454537, 4.494556, 4.53407, 4.573137, 4.612646, 4.652844, 4.691546, 4.730241, 4.769548, 4.808954, 4.848166, 4.886788, 4.924596, 4.964932, 5.002331, 5.041675, 5.079611, 5.119003, 5.157832, 5.195143, 5.233462, 5.272967, 5.309575, 5.34777, 5.385219, 5.424005, 5.460177, 5.498477, 5.536383, 5.573853, 5.611051, 5.649545, 5.681668, 5.714473, 5.754349, 5.792145, 5.827842, 5.864406, 5.90019, 5.935998, 5.971837, 6.006815, 6.042283, 6.077342, 6.113801, 6.148478, 6.183567, 6.217927, 6.250484, 6.286164, 6.320046, 6.353316, 6.390127, 6.421093, 6.454528, 6.488252, 6.521274, 6.55368, 6.587967, 6.619863, 6.652072, 6.683296, 6.716748, 6.747626, 6.779804, 6.809466, 6.842108, 6.874609, 6.906345, 6.936059, 6.965972, 6.997164, 7.026077, 7.057069, 7.087333, 7.117085, 7.146119, 7.174776, 7.203241, 7.232244, 7.261131, 7.289438, 7.321522, 7.347113, 7.374806, 7.402523, 7.430556, 7.459258, 7.484737, 7.511082, 7.538433, 7.565582, 7.590902, 7.61812, 7.644321, 7.668923, 7.696538, 7.72083, 7.74656, 7.769579, 7.795358, 7.82115, 7.843931, 7.868082, 7.892649, 7.91788, 7.94145, 7.964961, 7.985553, 8.011416, 8.032514, 8.056345, 8.078786, 8.102238, 8.124208, 8.145158, 8.167543, 8.189585, 8.210472, 8.230564, 8.253593, 8.273302, 8.295477, 8.315593, 8.337333, 8.356154, 8.375918, 8.395415, 8.414915, 8.435646, 8.452592, 8.468315, 8.496665, 8.512015, 8.530243, 8.547908, 8.567928, 8.583755, 8.602368, 8.618826, 8.638003, 8.656131, 8.673425, 8.690616, 8.706537, 8.724835, 8.740482, 8.758547, 8.773154, 8.789205, 8.804519, 8.821653, 8.837869, 8.852677, 8.866246, 8.882566, 8.896782, 8.911884, 8.92673, 8.942536, 8.957976, 8.969857, 8.985387, 8.999383, 9.013761, 9.026139, 9.040682, 9.053298, 9.067808, 9.080814, 9.092962, 9.106121, 9.118967, 9.131265, 9.143848, 9.156803, 9.168317, 9.181745, 9.190795, 9.205922, 9.21777, 9.227748, 9.241644, 9.251046, 9.26344, 9.273151, 9.285307, 9.297731, 9.306798, 9.31811, 9.326364, 9.337896, 9.347021, 9.359269, 9.368775, 9.376974, 9.388103, 9.399684, 9.410382, 9.417696, 9.427852, 9.436247, 9.444277, 9.454504, 9.465088, 9.471152, 9.478144, 9.743258, 9.497806, 10.01812, 9.516656, 10.03661, 9.532885, 10.05425, 9.549482, 10.0738, 9.565433, 10.09238, 9.580357, 10.10931, 9.595258, 10.12808, 9.610494, 10.14686, 9.623058, 10.16409, 9.639931, 10.18173, 9.652824, 10.19647, 9.665642, 10.20873, 9.676448, 10.22446, 9.688126, 10.23843, 9.703524, 10.2548, 10.26748, 9.721395, 10.27891, 9.733799, 10.29222, 9.743595, 10.30654, 9.756876, 9.761213, 10.32791, 9.769847, 10.3397, 9.783783, 10.35408, 9.793142, 9.798229, 10.37253, 9.809449, 10.38599, 9.818731, 10.39534, 10.40429, 9.830008, 10.41546, 9.839277, 10.42628, 9.849459, 9.853143, 10.44513, 9.864572, 10.45506, 10.46143, 9.874453, 10.47116, 9.879981, 9.885159, 10.48732, 9.893572, 10.49715, 9.900431, 9.903706, 10.51357, 9.912275, 9.915147, 10.52859, 9.92287, 10.53382, 9.927802, 9.930765, 10.54961, 9.93827, 9.940211, 10.56332, 9.949097, 9.95298, 10.57743, 9.956952, 9.959019, 10.58836, 9.962772, 10.5956, 10.60052, 9.97183, 10.60751, 10.61544, 9.982425, 10.62274, 10.62633, 9.988337, 10.6314, 10.63572, 9.99337, 10.64262, 10.64852, 10.00069, 10.00259, 10.65689, 10.00856, 10.01267, 10.66912, 10.01678, 10.0168, 10.67951, 10.01967, 10.02143, 10.68543, 10.02443, 10.02488, 10.69423, 10.69908, 10.03246, 10.03584, 10.71185, 10.03844, 10.04133, 
1e-100, 0.005579784, 0.01141804, 0.01756339, 0.02393259, 0.03062669, 0.0375999, 0.04481025, 0.05241462, 0.06026909, 0.06842879, 0.076964, 0.08574789, 0.09488955, 0.1043645, 0.1141438, 0.1243485, 0.1348995, 0.1457125, 0.1569356, 0.1685567, 0.1804522, 0.1927886, 0.2054828, 0.2185275, 0.2320017, 0.2458879, 0.2601297, 0.2746407, 0.2897573, 0.3051192, 0.320897, 0.3371993, 0.3537755, 0.3707961, 0.3882767, 0.4061801, 0.4243865, 0.4429722, 0.4621199, 0.4815231, 0.5015054, 0.5216908, 0.5424784, 0.5636257, 0.5850809, 0.6071968, 0.6295798, 0.6523872, 0.6754975, 0.6992856, 0.7232203, 0.7474472, 0.7721635, 0.7976079, 0.8231421, 0.849014, 0.8753353, 0.9023557, 0.9295249, 0.9568319, 0.984886, 1.013096, 1.042156, 1.070993, 1.100216, 1.130284, 1.160217, 1.190814, 1.221234, 1.252864, 1.284131, 1.316154, 1.348295, 1.380745, 1.413708, 1.446892, 1.480708, 1.514412, 1.548519, 1.583129, 1.618022, 1.652991, 1.689014, 1.724454, 1.760317, 1.796721, 1.833015, 1.870308, 1.907428, 1.945089, 1.982842, 2.02043, 2.058941, 2.098011, 2.136054, 2.175404, 2.21453, 2.254227, 2.29274, 2.335635, 2.374665, 2.414993, 2.456049, 2.496141, 2.537905, 2.578951, 2.620625, 2.662943, 2.704689, 2.746374, 2.788398, 2.8315, 2.874592, 2.916661, 2.960158, 3.002928, 3.046889, 3.090154, 3.133609, 3.177535, 3.220705, 3.265256, 3.309727, 3.35379, 3.398701, 3.443533, 3.487576, 3.533158, 3.577545, 3.622725, 3.667552, 3.712973, 3.757426, 3.802721, 3.847794, 3.893341, 3.939172, 3.984584, 4.030139, 4.076809, 4.121687, 4.167697, 4.213559, 4.258888, 4.305626, 4.350927, 4.396347, 4.442747, 4.487933, 4.535678, 4.58096, 4.626946, 4.673237, 4.719163, 4.764967, 4.80995, 4.855868, 4.901775, 4.948001, 4.992969, 5.038688, 5.084424, 5.130764, 5.175418, 5.221343, 5.268161, 5.311354, 5.357656, 5.402756, 5.448332, 5.491614, 5.539161, 5.582524, 5.6262, 5.672089, 5.716366, 5.760725, 5.805717, 5.849754, 5.8951, 5.937217, 5.981916, 6.026573, 6.069723, 6.113694, 6.156369, 6.19939, 6.243957, 6.286364, 6.328392, 6.372023, 6.415184, 6.457562, 6.499295, 6.543164, 6.584104, 6.626245, 6.671069, 6.708381, 6.744431, 6.789454, 6.830332, 6.871532, 6.912609, 6.953549, 6.993461, 7.032753, 7.072851, 7.113959, 7.15113, 7.192835, 7.231494, 7.270171, 7.308529, 7.348552, 7.386366, 7.425154, 7.46296, 7.501817, 7.537359, 7.576977, 7.612731, 7.65066, 7.688545, 7.723193, 7.759943, 7.79431, 7.831248, 7.867612, 7.902989, 7.939391, 7.972165, 8.007388, 8.043774, 8.076637, 8.111813, 8.145345, 8.180511, 8.211885, 8.246009, 8.278862, 8.311708, 8.344952, 8.377958, 8.408794, 8.441418, 8.474151, 8.50482, 8.538929, 8.569079, 8.600311, 8.630258, 8.659815, 8.689296, 8.719684, 8.748885, 8.778934, 8.80924, 8.837028, 8.866699, 8.893866, 8.924087, 8.953508, 8.980385, 9.00814, 9.034604, 9.061674, 9.089347, 9.116122, 9.143737, 9.168885, 9.195952, 9.221507, 9.248124, 9.273041, 9.299101, 9.324573, 9.348415, 9.373153, 9.39942, 9.423356, 9.44802, 9.471213, 9.494215, 9.518204, 9.542463, 9.564966, 9.587049, 9.609651, 9.632683, 9.65249, 9.676906, 9.698452, 9.720633, 9.74155, 9.762825, 9.781976, 9.799844, 9.831568, 9.847141, 9.865681, 9.885667, 9.905956, 9.926479, 9.945912, 9.964965, 9.982581, 10.00473, 10.02381, 10.04215, 10.05952, 10.07691, 10.09515, 10.11093, 10.13044, 10.14775, 10.1658, 10.18248, 10.19883, 10.21516, 10.23266, 10.24953, 10.26508, 10.28077, 10.29675, 10.31538, 10.33162, 10.34605, 10.36048, 10.37547, 10.3897, 10.40579, 10.42066, 10.43356, 10.44977, 10.46346, 10.47473, 10.48951, 10.50724, 10.51894, 10.5309, 10.54389, 10.55802, 10.57333, 10.58517, 10.59738, 10.61009, 10.61982, 10.63739, 10.64832, 10.66017, 10.67114, 10.68222, 10.69318, 10.70679, 10.71909, 10.7281, 10.73909, 10.7513, 10.76361, 10.77587, 10.78329, 10.79483, 10.80592, 10.81505, 10.82537, 10.83482, 10.84635, 10.85528, 10.86514, 10.87283, 10.88362, 10.89285, 11.18685, 10.91314, 11.49144, 10.93197, 11.51495, 10.95001, 11.53467, 10.96505, 11.55312, 10.98258, 11.57459, 10.99855, 11.59221, 11.01214, 11.61255, 11.03344, 11.634, 11.04753, 11.6492, 11.06169, 11.66533, 11.07789, 11.6857, 11.09168, 11.70002, 11.10252, 11.7144, 11.11482, 11.73374, 11.13041, 11.74956, 11.75781, 11.15167, 11.77569, 11.1634, 11.78931, 11.17301, 11.80131, 11.18535, 11.19297, 11.82601, 11.20306, 11.83915, 11.21452, 11.85651, 11.22798, 11.2303, 11.87452, 11.24124, 11.88765, 11.25037, 11.90123, 11.9094, 11.26584, 11.9193, 11.27674, 11.93437, 11.28728, 11.28971, 11.95302, 11.30042, 11.96265, 11.97009, 11.31063, 11.98026, 11.31818, 11.32193, 11.99495, 11.3336, 12.01079, 11.34252, 11.34536, 12.02652, 11.35107, 11.3563, 12.04095, 11.36298, 12.04797, 11.36829, 11.3694, 12.06255, 11.38251, 11.3846, 12.08033, 11.38985, 11.39228, 12.09491, 11.39794, 11.40326, 12.10895, 11.40667, 12.11378, 12.11906, 11.41592, 12.12927, 12.13542, 11.42398, 12.14138, 12.148, 11.43424, 12.15291, 12.15873, 11.43879, 12.16365, 12.16721, 11.44497, 11.44854, 12.18052, 11.4527, 11.45632, 12.19238, 11.46299, 11.46418, 12.20284, 11.46586, 11.4684, 12.21217, 11.46992, 11.47264, 12.22461, 12.22963, 11.47978, 11.47923, 12.23852, 11.48406, 11.48614, 
1e-100, 0.007327661, 0.01497102, 0.02293944, 0.03131698, 0.03998186, 0.04902603, 0.05847529, 0.0682203, 0.07841806, 0.08897456, 0.09988835, 0.1113104, 0.1230657, 0.1351871, 0.1478291, 0.1608471, 0.1742182, 0.1881573, 0.2024334, 0.2172027, 0.2324745, 0.248125, 0.2641414, 0.2808113, 0.2978745, 0.3152439, 0.3332821, 0.3517343, 0.3705541, 0.3900621, 0.4100466, 0.4302704, 0.4510248, 0.4724398, 0.494302, 0.51651, 0.5392654, 0.5625499, 0.5862008, 0.6105078, 0.6353964, 0.6604564, 0.6860442, 0.7122856, 0.7389684, 0.7661137, 0.7934691, 0.8217543, 0.8502254, 0.879363, 0.9086516, 0.938696, 0.9692358, 0.9999308, 1.031302, 1.063271, 1.09529, 1.127799, 1.160813, 1.19466, 1.228335, 1.262639, 1.297628, 1.332996, 1.368676, 1.404955, 1.440826, 1.477834, 1.515282, 1.552558, 1.59033, 1.629069, 1.668074, 1.706996, 1.74651, 1.786139, 1.826712, 1.867366, 1.908554, 1.949231, 1.991275, 2.033087, 2.07559, 2.117925, 2.160828, 2.20439, 2.24759, 2.291605, 2.335264, 2.379993, 2.424807, 2.470065, 2.515677, 2.561277, 2.606757, 2.653495, 2.699731, 2.745807, 2.791503, 2.842245, 2.888776, 2.935942, 2.983873, 3.03128, 3.08003, 3.128584, 3.177896, 3.226453, 3.274966, 3.32486, 3.374378, 3.424827, 3.474049, 3.523335, 3.57434, 3.624186, 3.674617, 3.724997, 3.776482, 3.82753, 3.878454, 3.928915, 3.981209, 4.032557, 4.084063, 4.135143, 4.187295, 4.239892, 4.291202, 4.342863, 4.395669, 4.447296, 4.497852, 4.551174, 4.60366, 4.655901, 4.708698, 4.760449, 4.813647, 4.866406, 4.919561, 4.971021, 5.023647, 5.076677, 5.128331, 5.18145, 5.233708, 5.286669, 5.33799, 5.39164, 5.443361, 5.495901, 5.549688, 5.601425, 5.652518, 5.704863, 5.758612, 5.809455, 5.862933, 5.913859, 5.965123, 6.018472, 6.069295, 6.121956, 6.171883, 6.226008, 6.276435, 6.32707, 6.379393, 6.429711, 6.479637, 6.530743, 6.581597, 6.630421, 6.68262, 6.732374, 6.782442, 6.831758, 6.883076, 6.933133, 6.982206, 7.030916, 7.081549, 7.130076, 7.178145, 7.227183, 7.27731, 7.325045, 7.372056, 7.421601, 7.47062, 7.516564, 7.564804, 7.612475, 7.659107, 7.705886, 7.755062, 7.803062, 7.845107, 7.886703, 7.934435, 7.981656, 8.027208, 8.072468, 8.118628, 8.161915, 8.207752, 8.251823, 8.29706, 8.339518, 8.38424, 8.427386, 8.469811, 8.513579, 8.556283, 8.599011, 8.641919, 8.684626, 8.725366, 8.766845, 8.808327, 8.850265, 8.888027, 8.93319, 8.971427, 9.01152, 9.051377, 9.090783, 9.130213, 9.169616, 9.207874, 9.245161, 9.286803, 9.323638, 9.362485, 9.398692, 9.436727, 9.473874, 9.509405, 9.548629, 9.584374, 9.619446, 9.656919, 9.689877, 9.725207, 9.760976, 9.796731, 9.829943, 9.867403, 9.899236, 9.932163, 9.967338, 10.00135, 10.03303, 10.06534, 10.09691, 10.13005, 10.16212, 10.19518, 10.22799, 10.25675, 10.28963, 10.31911, 10.34881, 10.38184, 10.40873, 10.43944, 10.46788, 10.49703, 10.52622, 10.55655, 10.58495, 10.61201, 10.63886, 10.66965, 10.69597, 10.72484, 10.75155, 10.77798, 10.80462, 10.82922, 10.85878, 10.8828, 10.90756, 10.93366, 10.95806, 10.98242, 11.00934, 11.03235, 11.05799, 11.0812, 11.10386, 11.12722, 11.15136, 11.17587, 11.19737, 11.21778, 11.2383, 11.26964, 11.28869, 11.30985, 11.3319, 11.35114, 11.37298, 11.39501, 11.41585, 11.43724, 11.4578, 11.47773, 11.49532, 11.51812, 11.53626, 11.55644, 11.57498, 11.59206, 11.61093, 11.63081, 11.65084, 11.66732, 11.68329, 11.7036, 11.71972, 11.73777, 11.75715, 11.77382, 11.79082, 11.80549, 11.82591, 11.83984, 11.85598, 11.87278, 11.88833, 11.90265, 11.91839, 11.93441, 11.95061, 11.96241, 11.97868, 11.99309, 12.00933, 12.02367, 12.03707, 12.05142, 12.06277, 12.08206, 12.09476, 12.10571, 12.11716, 12.1342, 12.14456, 12.15673, 12.17184, 12.18496, 12.19544, 12.20865, 12.21822, 12.23254, 12.24674, 12.2574, 12.26875, 12.27885, 12.29269, 12.30366, 12.31469, 12.32587, 12.33407, 12.34455, 12.3566, 12.3674, 12.37737, 12.38738, 12.3954, 12.40833, 12.73708, 12.42946, 13.07509, 12.4473, 13.09527, 12.46506, 13.11874, 12.48419, 13.14087, 12.50119, 13.16038, 12.51754, 13.18125, 12.53712, 13.2029, 12.55424, 13.22152, 12.57073, 13.24178, 12.58595, 13.26143, 12.60147, 13.27728, 12.61544, 13.29692, 12.62905, 13.31205, 12.64102, 13.33046, 12.65942, 13.34963, 13.35875, 12.67943, 13.37462, 12.69, 13.38976, 12.70717, 13.40575, 12.71532, 12.72132, 13.4295, 12.73698, 13.44673, 12.74816, 13.46187, 12.76039, 12.76567, 13.4842, 12.77635, 13.49688, 12.78409, 13.50798, 13.51649, 12.80149, 13.53235, 12.81115, 13.54517, 12.81922, 12.83022, 13.56731, 12.83616, 13.57624, 13.58305, 12.84885, 13.59391, 12.85706, 12.86358, 13.61374, 12.87165, 13.62376, 12.88076, 12.88335, 13.64089, 12.89238, 12.89455, 13.65761, 12.90065, 13.66794, 12.90925, 12.91142, 13.68405, 12.92143, 12.92467, 13.70412, 12.93293, 12.93424, 13.71672, 12.93921, 12.94093, 13.73134, 12.94766, 13.73633, 13.7445, 12.95742, 13.75261, 13.76018, 12.96914, 13.76675, 13.7717, 12.97607, 13.77754, 13.78329, 12.98292, 13.79248, 13.79688, 12.99089, 12.99458, 13.81232, 12.99906, 13.00033, 13.8222, 13.00645, 13.00779, 13.8328, 13.00941, 13.01234, 13.84021, 13.0155, 13.01865, 13.85428, 13.8623, 13.02711, 13.02751, 13.87402, 13.02774, 13.03234, 
1e-100, 0.009474781, 0.01939897, 0.02976054, 0.04047373, 0.05172954, 0.06336864, 0.07547818, 0.08810359, 0.1011036, 0.1146333, 0.1286815, 0.1431251, 0.1581233, 0.1736651, 0.1895709, 0.2062082, 0.2232419, 0.2407372, 0.2588991, 0.2775404, 0.2966667, 0.3163757, 0.3366225, 0.3573768, 0.3788039, 0.4007785, 0.4232224, 0.446224, 0.4698272, 0.493998, 0.518609, 0.5440453, 0.5697729, 0.5962409, 0.623285, 0.6509364, 0.6789328, 0.7076392, 0.737046, 0.7666992, 0.7969752, 0.8280131, 0.8596541, 0.8915366, 0.9240653, 0.9572688, 0.990918, 1.02512, 1.059908, 1.095387, 1.131223, 1.1673, 1.204268, 1.241614, 1.279572, 1.317983, 1.356867, 1.39625, 1.4362, 1.476304, 1.517296, 1.558854, 1.600867, 1.642974, 1.685301, 1.728996, 1.772015, 1.81645, 1.860512, 1.905807, 1.950929, 1.996631, 2.042953, 2.089636, 2.136465, 2.18397, 2.231621, 2.279832, 2.328091, 2.377082, 2.426474, 2.475768, 2.525857, 2.576339, 2.626721, 2.677811, 2.728751, 2.780549, 2.832554, 2.885234, 2.937091, 2.990178, 3.043189, 3.096451, 3.150561, 3.204208, 3.257394, 3.313, 3.366073, 3.424828, 3.478325, 3.533813, 3.589112, 3.645572, 3.702037, 3.758031, 3.814119, 3.870811, 3.92798, 3.985749, 4.0419, 4.100477, 4.157098, 4.215752, 4.273083, 4.330127, 4.389906, 4.44737, 4.506867, 4.564232, 4.623777, 4.683006, 4.741486, 4.799801, 4.85903, 4.919029, 4.977651, 5.037009, 5.096335, 5.15508, 5.214955, 5.275404, 5.332754, 5.393472, 5.452895, 5.51382, 5.57325, 5.632993, 5.691998, 5.751478, 5.812397, 5.872169, 5.930064, 5.990702, 6.049194, 6.109047, 6.169401, 6.227539, 6.284863, 6.347287, 6.406528, 6.464505, 6.524772, 6.583093, 6.64134, 6.701649, 6.759607, 6.819109, 6.876873, 6.935631, 6.994587, 7.052015, 7.109614, 7.169012, 7.225829, 7.285272, 7.342313, 7.39873, 7.455338, 7.512046, 7.569615, 7.626911, 7.682892, 7.738672, 7.796272, 7.85111, 7.907522, 7.964522, 8.017982, 8.075101, 8.129051, 8.183981, 8.239595, 8.293125, 8.346812, 8.401692, 8.455611, 8.509756, 8.562616, 8.617801, 8.670245, 8.723304, 8.774065, 8.827217, 8.877992, 8.932631, 8.984303, 9.039653, 9.083103, 9.130806, 9.182955, 9.23397, 9.285149, 9.336693, 9.385979, 9.433802, 9.484374, 9.532347, 9.580489, 9.629928, 9.677755, 9.726237, 9.773338, 9.822244, 9.867281, 9.914044, 9.961142, 10.00507, 10.05356, 10.09729, 10.14341, 10.18909, 10.23254, 10.27859, 10.32058, 10.36357, 10.41044, 10.45211, 10.49404, 10.53697, 10.57954, 10.61726, 10.66399, 10.70509, 10.74641, 10.78616, 10.82789, 10.86834, 10.90654, 10.94734, 10.98636, 11.02739, 11.06622, 11.10309, 11.14075, 11.18039, 11.21773, 11.25446, 11.29716, 11.33105, 11.36547, 11.40412, 11.43725, 11.47377, 11.50958, 11.54494, 11.57882, 11.61376, 11.64812, 11.68163, 11.71692, 11.74997, 11.78244, 11.81408, 11.85035, 11.87904, 11.91258, 11.94278, 11.97608, 12.00519, 12.03778, 12.06902, 12.0994, 12.13, 12.15938, 12.18966, 12.21662, 12.24758, 12.27553, 12.30484, 12.33277, 12.3627, 12.38708, 12.41946, 12.44433, 12.46894, 12.49614, 12.52331, 12.5486, 12.57318, 12.59986, 12.62677, 12.65172, 12.67718, 12.7001, 12.72446, 12.74817, 12.76929, 12.80251, 12.82139, 12.84402, 12.86915, 12.89179, 12.91423, 12.93721, 12.95692, 12.98297, 13.00266, 13.02435, 13.04469, 13.06665, 13.08621, 13.10497, 13.12622, 13.14808, 13.16766, 13.18778, 13.20628, 13.22503, 13.24561, 13.26433, 13.28474, 13.29939, 13.32012, 13.33951, 13.35745, 13.3753, 13.39304, 13.40804, 13.42821, 13.44618, 13.45855, 13.47554, 13.49252, 13.5075, 13.52568, 13.54082, 13.55758, 13.57238, 13.58676, 13.60324, 13.61916, 13.63392, 13.64887, 13.66325, 13.67753, 13.69237, 13.7046, 13.72044, 13.73275, 13.74604, 13.75925, 13.77346, 13.78638, 13.79832, 13.81301, 13.82256, 13.83804, 13.8519, 13.86387, 13.87574, 13.88647, 13.89865, 13.90873, 13.92199, 13.93549, 13.94466, 13.95473, 13.96513, 13.97499, 13.98773, 13.99928, 14.00935, 14.36822, 14.03027, 14.73947, 14.05142, 14.76209, 14.07025, 14.7833, 14.08677, 14.80655, 14.109, 14.82871, 14.12553, 14.85025, 14.1459, 14.87441, 14.16193, 14.89348, 14.18006, 14.91388, 14.19399, 14.93437, 14.21409, 14.95485, 14.22475, 14.96942, 14.24253, 14.99082, 14.25719, 15.00942, 14.27277, 15.02718, 15.03971, 14.29504, 15.0565, 14.30828, 15.07183, 14.32092, 15.08926, 14.33395, 14.34017, 15.11538, 14.35318, 15.13269, 14.36846, 15.14625, 14.37879, 14.38358, 15.17116, 14.3939, 15.18596, 14.40855, 15.198, 15.21042, 14.42345, 15.22489, 14.43332, 15.23684, 14.44428, 14.45069, 15.25945, 14.45945, 15.26886, 15.27657, 14.47214, 15.28877, 14.48118, 14.48784, 15.312, 14.49882, 15.32298, 14.50577, 14.51021, 15.34001, 14.5185, 14.52107, 15.35558, 14.52733, 15.36494, 14.53657, 14.54177, 15.38857, 14.55128, 14.55276, 15.40334, 14.55908, 14.56189, 15.41884, 14.56897, 14.57056, 15.4341, 14.57149, 15.44102, 15.45023, 14.58888, 15.46082, 15.46465, 14.59763, 15.47391, 15.47731, 14.6072, 15.4862, 15.48949, 14.6111, 15.4988, 15.50402, 14.62287, 14.62715, 15.51964, 14.63014, 14.6375, 15.53377, 14.6374, 14.63929, 15.54162, 14.64314, 14.64202, 15.55308, 14.65211, 14.65443, 15.56697, 15.5736, 14.65854, 14.65941, 15.58723, 14.6637, 14.66698, 
1e-100, 0.0122764, 0.02498988, 0.03826255, 0.05209633, 0.06637066, 0.08128825, 0.09672573, 0.1126355, 0.1292879, 0.1464007, 0.1641129, 0.1825186, 0.2013995, 0.2208707, 0.2410802, 0.2617901, 0.2831569, 0.3052031, 0.3277749, 0.3510439, 0.3750563, 0.3995117, 0.424695, 0.45063, 0.4770665, 0.5041016, 0.5320208, 0.5604446, 0.5894452, 0.6193619, 0.6497292, 0.6808161, 0.7123882, 0.7449966, 0.7778793, 0.811482, 0.8460202, 0.8809553, 0.9164731, 0.9527397, 0.9899197, 1.027321, 1.065317, 1.104226, 1.143677, 1.183596, 1.224427, 1.265717, 1.307394, 1.350063, 1.392988, 1.436536, 1.480833, 1.525605, 1.571006, 1.617, 1.663325, 1.710053, 1.757883, 1.806259, 1.85448, 1.903479, 1.953388, 2.00363, 2.054104, 2.105427, 2.156605, 2.208815, 2.261013, 2.314099, 2.367112, 2.421508, 2.475865, 2.530615, 2.585858, 2.641106, 2.697615, 2.754362, 2.81072, 2.867202, 2.925357, 2.983124, 3.041921, 3.100241, 3.159451, 3.218713, 3.278086, 3.337677, 3.398186, 3.458836, 3.519981, 3.580983, 3.642288, 3.704285, 3.76686, 3.829202, 3.891513, 3.954597, 4.015578, 4.083555, 4.144536, 4.209351, 4.273125, 4.336535, 4.401299, 4.465698, 4.531578, 4.595736, 4.661383, 4.726758, 4.791699, 4.85814, 4.923615, 4.988883, 5.054542, 5.121084, 5.18891, 5.25412, 5.320993, 5.386952, 5.453763, 5.521204, 5.587642, 5.654616, 5.721562, 5.787539, 5.855096, 5.923077, 5.990162, 6.057092, 6.12366, 6.191956, 6.25796, 6.325076, 6.392397, 6.460523, 6.526666, 6.593797, 6.660041, 6.728044, 6.795287, 6.86217, 6.927939, 6.994818, 7.061971, 7.128146, 7.194284, 7.260322, 7.324755, 7.394538, 7.459865, 7.526365, 7.592432, 7.65827, 7.723638, 7.788493, 7.85426, 7.920494, 7.985958, 8.048955, 8.115839, 8.178823, 8.242666, 8.306768, 8.371977, 8.437327, 8.498566, 8.562521, 8.626011, 8.688895, 8.751627, 8.814741, 8.876862, 8.938592, 9.001496, 9.06259, 9.123229, 9.187942, 9.246289, 9.308183, 9.368082, 9.428893, 9.489611, 9.548168, 9.609276, 9.669067, 9.727028, 9.786576, 9.844496, 9.904068, 9.96093, 10.01947, 10.07527, 10.13388, 10.19166, 10.24878, 10.30557, 10.36562, 10.41542, 10.4643, 10.52244, 10.57634, 10.6345, 10.68795, 10.74259, 10.79659, 10.84937, 10.90166, 10.95589, 11.00641, 11.06009, 11.10969, 11.16249, 11.21532, 11.2656, 11.31558, 11.36726, 11.41468, 11.46691, 11.51453, 11.56335, 11.61101, 11.66017, 11.71005, 11.75653, 11.80446, 11.85136, 11.89596, 11.94417, 11.98872, 12.03647, 12.07938, 12.12582, 12.17085, 12.21428, 12.25936, 12.30287, 12.34671, 12.38871, 12.43253, 12.47231, 12.51709, 12.5579, 12.599, 12.63924, 12.68118, 12.72113, 12.76019, 12.80326, 12.8416, 12.88134, 12.92, 12.95736, 12.99424, 13.03474, 13.06986, 13.10853, 13.14558, 13.18246, 13.21746, 13.25548, 13.28954, 13.32591, 13.36158, 13.39833, 13.43039, 13.46334, 13.49663, 13.53172, 13.56473, 13.5984, 13.63181, 13.66038, 13.69613, 13.72736, 13.75921, 13.79006, 13.82046, 13.84984, 13.88112, 13.91201, 13.94159, 13.97128, 14.00015, 14.02887, 14.05616, 14.08584, 14.11268, 14.14146, 14.16738, 14.19447, 14.22125, 14.25058, 14.27722, 14.3, 14.32516, 14.35141, 14.37339, 14.40962, 14.43232, 14.45608, 14.47853, 14.50282, 14.52739, 14.55154, 14.57608, 14.59746, 14.62017, 14.64326, 14.66608, 14.68781, 14.71055, 14.72952, 14.751, 14.77086, 14.79556, 14.81639, 14.83503, 14.85409, 14.87688, 14.89935, 14.9194, 14.93719, 14.95707, 14.97487, 14.99318, 15.01469, 15.03109, 15.05015, 15.06918, 15.08704, 15.10142, 15.11949, 15.13741, 15.15446, 15.17062, 15.18576, 15.20536, 15.22379, 15.23961, 15.25426, 15.26957, 15.28565, 15.30385, 15.31611, 15.33123, 15.34707, 15.35988, 15.37429, 15.38931, 15.40675, 15.41845, 15.43159, 15.446, 15.46175, 15.47622, 15.48901, 15.50185, 15.51353, 15.52559, 15.54083, 15.55178, 15.56509, 15.57513, 15.58716, 15.59953, 15.61396, 15.62346, 15.63437, 15.64566, 15.65804, 15.67218, 15.68194, 16.06944, 15.702, 16.47334, 15.72157, 16.49563, 15.74444, 16.52241, 15.7631, 16.54506, 15.78027, 16.56867, 15.80387, 16.59213, 15.82352, 16.61629, 15.84139, 16.64025, 15.85846, 16.65947, 15.87395, 16.67833, 15.89195, 16.70246, 15.90459, 16.71998, 15.92733, 16.7379, 15.94166, 16.76215, 15.95509, 16.77914, 16.7893, 15.98022, 16.80617, 15.99541, 16.82574, 16.0088, 16.84116, 16.0249, 16.02848, 16.87291, 16.04579, 16.88815, 16.05759, 16.9038, 16.07059, 16.07189, 16.93155, 16.08156, 16.94461, 16.09193, 16.95938, 16.97225, 16.11465, 16.98923, 16.12565, 17.00243, 16.13963, 16.14178, 17.02438, 16.15075, 17.03648, 17.04019, 16.16851, 17.0582, 16.17765, 16.18522, 17.07577, 16.1932, 17.08796, 16.19966, 16.20906, 17.10704, 16.21309, 16.21413, 17.1262, 16.22704, 17.13884, 16.23808, 16.23637, 17.16047, 16.25299, 16.2509, 17.17852, 16.25948, 16.25919, 17.19222, 16.2667, 16.27023, 17.21042, 16.27424, 17.21957, 17.22818, 16.29176, 17.24044, 17.24038, 16.29747, 17.2534, 17.25566, 16.30568, 17.26872, 17.27203, 16.31422, 17.2843, 17.28687, 16.32883, 16.33143, 17.30018, 16.33411, 16.33912, 17.3144, 16.33986, 16.33965, 17.32533, 16.3474, 16.34736, 17.33778, 16.35628, 16.35658, 17.35379, 17.36275, 16.36578, 16.36293, 17.37608, 16.3636, 16.36889, 
1e-100, 0.01551353, 0.0317588, 0.04858765, 0.06598909, 0.08417161, 0.1028862, 0.1223728, 0.1425462, 0.1632726, 0.1848185, 0.2070233, 0.2298249, 0.2535011, 0.2778028, 0.3027907, 0.3287072, 0.3551789, 0.3822738, 0.4103468, 0.4391513, 0.468424, 0.4988493, 0.5297021, 0.5613209, 0.5939714, 0.6271894, 0.6610334, 0.6958356, 0.7313203, 0.7674573, 0.8044013, 0.8423038, 0.8806528, 0.9196111, 0.9598902, 1.000596, 1.041849, 1.084052, 1.127007, 1.170461, 1.214856, 1.259842, 1.305709, 1.352064, 1.398953, 1.446904, 1.495427, 1.544477, 1.594157, 1.64509, 1.696103, 1.747773, 1.800145, 1.852901, 1.906841, 1.96082, 2.01583, 2.071306, 2.127189, 2.183641, 2.240894, 2.298743, 2.35692, 2.415751, 2.4751, 2.535415, 2.594366, 2.65546, 2.716822, 2.778845, 2.840844, 2.90346, 2.967025, 3.03034, 3.094716, 3.159098, 3.223336, 3.289284, 3.354825, 3.420865, 3.4877, 3.554318, 3.621746, 3.689175, 3.757305, 3.825176, 3.893641, 3.962859, 4.032438, 4.100858, 4.171308, 4.241431, 4.312329, 4.382517, 4.453955, 4.523827, 4.596337, 4.667496, 4.737633, 4.814335, 4.885018, 4.957703, 5.030554, 5.103441, 5.175696, 5.249272, 5.323183, 5.396413, 5.470766, 5.544153, 5.617927, 5.692028, 5.766554, 5.840691, 5.915523, 5.989987, 6.064583, 6.138886, 6.215042, 6.288968, 6.362604, 6.436636, 6.51344, 6.58874, 6.663153, 6.738072, 6.813605, 6.889246, 6.963356, 7.03895, 7.114896, 7.188996, 7.262138, 7.336639, 7.412581, 7.488282, 7.561737, 7.635559, 7.710793, 7.785377, 7.858501, 7.932814, 8.00767, 8.080518, 8.155841, 8.227579, 8.301999, 8.375264, 8.447009, 8.522847, 8.594643, 8.668824, 8.741474, 8.812218, 8.886465, 8.957874, 9.030505, 9.100597, 9.173638, 9.243498, 9.314674, 9.386331, 9.456213, 9.525887, 9.596067, 9.668253, 9.735809, 9.805791, 9.87513, 9.945767, 10.01171, 10.08195, 10.1495, 10.21587, 10.28693, 10.3523, 10.4194, 10.48669, 10.55155, 10.61743, 10.68451, 10.75011, 10.81544, 10.88029, 10.94701, 11.01033, 11.07472, 11.13787, 11.20193, 11.2657, 11.32671, 11.3899, 11.45246, 11.51279, 11.57723, 11.63863, 11.69976, 11.76281, 11.81855, 11.87058, 11.93395, 11.99428, 12.05445, 12.11272, 12.17119, 12.22837, 12.28463, 12.34266, 12.40238, 12.45545, 12.51215, 12.56856, 12.62337, 12.67883, 12.73309, 12.78626, 12.84097, 12.89487, 12.94899, 13.00078, 13.05352, 13.10396, 13.15609, 13.20761, 13.25796, 13.31088, 13.35926, 13.40841, 13.45822, 13.50814, 13.55419, 13.60464, 13.65324, 13.70074, 13.74774, 13.79499, 13.84051, 13.88815, 13.93201, 13.98006, 14.02531, 14.06873, 14.11381, 14.15397, 14.2002, 14.24457, 14.28795, 14.33102, 14.37623, 14.41391, 14.45686, 14.49678, 14.53847, 14.57884, 14.62007, 14.65869, 14.69678, 14.73613, 14.77766, 14.81644, 14.85297, 14.89106, 14.92617, 14.96632, 15.00525, 15.03864, 15.07493, 15.11119, 15.14453, 15.18031, 15.21613, 15.25216, 15.28492, 15.31914, 15.35221, 15.38425, 15.42015, 15.45277, 15.48273, 15.51572, 15.54678, 15.5773, 15.61154, 15.64069, 15.67018, 15.70011, 15.72969, 15.75674, 15.78845, 15.82094, 15.84671, 15.87348, 15.90336, 15.9301, 15.95904, 15.98366, 16.01145, 16.03483, 16.07146, 16.09395, 16.12195, 16.14599, 16.17116, 16.19555, 16.21787, 16.24723, 16.27058, 16.29286, 16.31749, 16.34069, 16.36198, 16.38607, 16.41139, 16.43257, 16.45269, 16.47642, 16.49747, 16.5203, 16.54178, 16.56286, 16.58439, 16.60371, 16.62596, 16.6489, 16.66646, 16.68726, 16.70618, 16.72373, 16.74393, 16.76468, 16.7824, 16.79739, 16.81795, 16.83487, 16.8549, 16.87119, 16.89073, 16.90772, 16.92509, 16.94426, 16.95943, 16.97399, 16.99318, 17.00958, 17.0235, 17.03918, 17.05447, 17.07047, 17.08459, 17.10127, 17.11449, 17.13035, 17.14628, 17.15985, 17.17437, 17.18769, 17.20623, 17.2188, 17.22935, 17.24246, 17.25718, 17.26846, 17.28062, 17.29562, 17.30944, 17.31899, 17.33235, 17.34466, 17.35856, 17.3716, 17.38384, 17.39659, 17.40577, 17.8252, 17.43015, 18.25321, 17.44964, 18.27873, 17.47135, 18.30461, 17.49292, 18.32903, 17.51459, 18.35506, 17.53463, 18.38214, 17.55685, 18.40448, 17.57359, 18.42871, 17.59071, 18.4499, 17.61317, 18.4723, 17.62696, 18.49447, 17.64068, 18.51672, 17.66835, 18.53569, 17.67931, 18.55976, 17.69604, 18.58289, 18.5892, 17.72458, 18.60562, 17.73776, 18.62441, 17.75458, 18.64474, 17.76966, 17.77204, 18.67222, 17.79271, 18.69378, 17.80551, 18.71043, 17.8189, 17.81542, 18.74209, 17.82995, 18.75722, 17.84511, 18.76899, 18.78475, 17.86349, 18.80231, 17.87264, 18.8173, 17.88655, 17.88872, 18.8384, 17.89975, 18.85267, 18.85622, 17.92353, 18.87152, 17.93006, 17.94179, 18.89502, 17.95083, 18.90584, 17.95177, 17.96348, 18.92632, 17.97067, 17.96349, 18.94436, 17.98447, 18.96189, 17.99943, 17.9935, 18.97935, 18.00917, 18.00457, 19.00272, 18.01849, 18.01627, 19.02289, 18.02209, 18.02441, 19.04139, 18.03418, 19.05354, 19.0614, 18.04846, 19.07108, 19.07047, 18.05598, 19.08631, 19.08354, 18.06413, 19.09829, 19.09616, 18.07495, 19.11856, 19.11759, 18.09466, 18.0965, 19.1342, 18.10296, 18.10319, 19.14749, 18.10683, 18.10296, 19.15669, 18.11274, 18.11447, 19.17504, 18.12417, 18.11992, 19.18797, 19.20273, 18.13301, 18.12827, 19.22176, 18.1343, 18.13495, 
1e-100, 0.01965965, 0.03993601, 0.06109678, 0.08298226, 0.1055111, 0.129024, 0.1532131, 0.1781463, 0.2040205, 0.2305551, 0.2580345, 0.2863783, 0.3153901, 0.3453398, 0.376162, 0.4076872, 0.4401233, 0.4735091, 0.5075045, 0.5426556, 0.5785636, 0.6151551, 0.6528312, 0.6913281, 0.7305161, 0.7706174, 0.8116905, 0.8535081, 0.8961849, 0.9398845, 0.9841868, 1.029221, 1.075401, 1.122378, 1.17015, 1.218629, 1.268035, 1.317998, 1.368983, 1.42092, 1.473421, 1.526697, 1.580796, 1.63582, 1.691298, 1.747429, 1.804797, 1.862615, 1.920926, 1.980462, 2.040441, 2.100776, 2.162329, 2.224123, 2.28679, 2.350404, 2.413996, 2.478446, 2.543652, 2.609888, 2.675789, 2.742591, 2.810304, 2.878579, 2.947023, 3.016098, 3.085466, 3.155731, 3.226412, 3.297438, 3.368633, 3.44102, 3.513918, 3.587008, 3.660268, 3.734101, 3.808787, 3.883587, 3.958258, 4.033349, 4.109959, 4.186068, 4.262829, 4.339976, 4.417093, 4.494415, 4.572226, 4.650276, 4.728807, 4.80767, 4.887449, 4.966404, 5.045994, 5.126073, 5.205789, 5.28586, 5.366737, 5.446689, 5.52581, 5.612868, 5.692258, 5.773399, 5.854129, 5.936452, 6.018415, 6.100972, 6.182992, 6.265382, 6.347706, 6.430013, 6.513489, 6.596379, 6.678799, 6.761801, 6.843699, 6.927353, 7.01004, 7.093574, 7.176455, 7.258844, 7.342371, 7.425001, 7.50825, 7.590719, 7.674661, 7.758213, 7.84035, 7.924034, 8.006188, 8.089007, 8.172365, 8.254119, 8.334658, 8.418647, 8.50038, 8.581772, 8.664641, 8.745446, 8.827758, 8.909544, 8.992775, 9.073228, 9.152618, 9.234652, 9.315502, 9.395943, 9.477126, 9.556587, 9.634724, 9.717436, 9.796968, 9.874809, 9.955452, 10.03402, 10.1121, 10.19094, 10.27046, 10.34606, 10.42493, 10.50183, 10.57887, 10.6559, 10.73112, 10.80957, 10.88414, 10.96236, 11.03674, 11.11213, 11.18649, 11.26213, 11.33566, 11.4095, 11.48465, 11.55414, 11.6291, 11.70005, 11.7731, 11.8461, 11.91654, 11.98925, 12.06003, 12.13013, 12.20036, 12.27042, 12.34062, 12.40886, 12.47769, 12.5471, 12.61347, 12.68136, 12.75063, 12.81639, 12.88328, 12.94965, 13.01637, 13.08275, 13.14834, 13.2146, 13.27351, 13.33075, 13.39741, 13.46339, 13.52702, 13.58906, 13.65124, 13.71334, 13.77278, 13.83438, 13.89508, 13.95395, 14.0134, 14.07464, 14.13238, 14.19083, 14.24841, 14.30733, 14.36285, 14.42037, 14.47755, 14.53198, 14.58881, 14.64673, 14.69834, 14.75472, 14.80819, 14.86083, 14.91342, 14.96702, 15.01879, 15.07158, 15.12282, 15.17221, 15.22459, 15.27465, 15.32808, 15.3761, 15.42372, 15.47389, 15.52007, 15.56678, 15.61748, 15.66474, 15.71198, 15.75728, 15.80364, 15.8487, 15.89442, 15.93813, 15.98812, 16.03103, 16.07243, 16.11739, 16.15651, 16.20243, 16.24562, 16.28597, 16.32836, 16.37033, 16.40894, 16.45132, 16.49084, 16.53162, 16.57103, 16.6089, 16.6499, 16.68485, 16.7242, 16.76336, 16.79875, 16.83551, 16.87287, 16.90552, 16.94646, 16.98135, 17.01456, 17.04946, 17.08546, 17.11862, 17.15354, 17.1863, 17.21802, 17.25498, 17.28628, 17.3167, 17.34798, 17.37953, 17.41253, 17.44121, 17.4723, 17.50284, 17.53481, 17.56394, 17.59326, 17.62197, 17.65013, 17.67952, 17.70836, 17.73433, 17.77102, 17.79561, 17.81874, 17.84867, 17.87371, 17.89924, 17.92668, 17.95283, 17.9774, 18.00246, 18.02834, 18.05357, 18.07575, 18.09951, 18.123, 18.14924, 18.17123, 18.19491, 18.21843, 18.24001, 18.2638, 18.28745, 18.3088, 18.33019, 18.35095, 18.37206, 18.39475, 18.41404, 18.434, 18.45584, 18.47609, 18.49341, 18.51488, 18.53329, 18.55187, 18.57182, 18.59054, 18.61002, 18.62884, 18.64883, 18.66645, 18.68519, 18.69945, 18.72164, 18.73777, 18.7518, 18.76855, 18.78561, 18.79954, 18.81689, 18.83522, 18.85138, 18.86573, 18.8814, 18.89793, 18.91439, 18.92841, 18.94383, 18.95821, 18.97284, 18.98717, 19.00117, 19.01538, 19.02739, 19.03993, 19.05417, 19.06824, 19.08138, 19.09574, 19.10886, 19.11947, 19.13556, 19.14864, 19.15953, 19.16948, 19.61452, 19.19515, 20.07286, 19.21644, 20.10252, 19.23936, 20.12671, 19.26112, 20.1551, 19.28749, 20.18273, 19.30755, 20.20778, 19.32826, 20.23422, 19.34661, 20.25876, 19.36595, 20.28009, 19.38484, 20.30501, 19.40499, 20.32837, 19.41696, 20.34868, 19.44944, 20.37282, 19.45611, 20.39747, 19.47513, 20.42551, 20.42569, 19.51199, 20.44948, 19.52752, 20.46893, 19.54528, 20.48957, 19.56095, 19.55631, 20.51866, 19.58318, 20.53748, 19.59701, 20.55494, 19.61103, 19.60292, 20.59329, 19.6173, 20.61603, 19.63648, 20.6231, 20.64454, 19.65656, 20.66066, 19.66677, 20.6756, 19.67997, 19.68186, 20.70044, 19.69195, 20.71312, 20.71307, 19.72708, 20.73145, 19.72309, 19.74076, 20.75286, 19.75058, 20.76717, 19.75164, 19.7687, 20.78791, 19.77689, 19.76513, 20.81061, 19.79389, 20.82454, 19.80296, 19.7948, 20.84531, 19.81922, 19.80707, 20.87639, 19.82474, 19.81715, 20.8974, 19.82547, 19.83053, 20.91832, 19.83499, 20.93214, 20.93701, 19.85033, 20.94709, 20.94408, 19.86319, 20.96406, 20.95882, 19.87321, 20.98265, 20.9756, 19.88304, 20.99981, 20.99367, 19.90739, 19.90995, 21.00952, 19.91563, 19.91758, 21.02305, 19.92386, 19.91278, 21.03941, 19.93322, 19.92231, 21.05661, 19.9452, 19.93243, 21.07088, 21.09019, 19.9507, 19.93622, 21.11072, 19.9516, 19.94714, 
1e-100, 0.02423932, 0.04957571, 0.07565572, 0.102627, 0.1305776, 0.1592974, 0.1891463, 0.2198371, 0.2513019, 0.2839658, 0.3174408, 0.3517927, 0.3872665, 0.4235363, 0.4608505, 0.4992555, 0.5383734, 0.5785223, 0.6198104, 0.6618752, 0.7049086, 0.749075, 0.7940092, 0.840051, 0.887109, 0.9350019, 0.9836361, 1.033602, 1.084423, 1.135967, 1.188513, 1.242088, 1.296272, 1.351937, 1.408259, 1.465508, 1.523411, 1.582505, 1.642299, 1.702924, 1.764184, 1.826931, 1.889993, 1.953896, 2.018843, 2.084482, 2.150746, 2.218028, 2.285909, 2.355325, 2.424027, 2.494201, 2.565153, 2.636709, 2.708697, 2.781867, 2.85556, 2.929602, 3.004494, 3.079974, 3.156188, 3.23271, 3.309981, 3.387932, 3.46597, 3.545502, 3.623916, 3.704387, 3.785003, 3.866213, 3.94701, 4.029199, 4.111225, 4.194633, 4.278564, 4.361457, 4.444972, 4.530056, 4.615272, 4.699786, 4.784952, 4.870514, 4.957668, 5.044098, 5.130969, 5.218077, 5.305524, 5.393617, 5.481062, 5.569464, 5.657803, 5.746244, 5.835202, 5.924211, 6.014041, 6.103469, 6.193512, 6.282948, 6.370572, 6.466718, 6.555147, 6.64528, 6.735954, 6.826571, 6.917027, 7.008394, 7.099317, 7.191063, 7.28144, 7.372712, 7.463439, 7.555915, 7.646417, 7.738512, 7.828675, 7.919639, 8.011704, 8.102421, 8.194064, 8.284915, 8.376432, 8.467243, 8.55782, 8.649137, 8.740564, 8.830602, 8.920699, 9.012333, 9.103185, 9.192626, 9.281972, 9.373994, 9.460245, 9.551257, 9.640132, 9.730222, 9.819851, 9.908427, 9.996527, 10.08583, 10.17494, 10.26361, 10.34885, 10.43792, 10.52553, 10.61114, 10.70018, 10.78554, 10.87063, 10.9614, 11.04514, 11.12989, 11.21496, 11.30197, 11.38535, 11.46914, 11.55313, 11.63862, 11.72166, 11.80416, 11.88809, 11.96885, 12.05179, 12.13373, 12.21542, 12.29693, 12.37817, 12.45848, 12.53898, 12.61811, 12.69765, 12.77769, 12.85506, 12.93203, 13.01192, 13.08818, 13.16517, 13.24292, 13.31757, 13.39688, 13.47126, 13.54668, 13.62112, 13.6968, 13.76949, 13.84216, 13.91554, 13.98973, 14.06052, 14.13392, 14.20579, 14.27622, 14.34473, 14.4179, 14.48721, 14.55804, 14.62748, 14.70056, 14.76153, 14.82313, 14.89255, 14.96119, 15.02879, 15.09442, 15.15861, 15.22497, 15.29102, 15.35474, 15.41952, 15.48137, 15.54426, 15.60829, 15.6685, 15.73127, 15.79301, 15.85343, 15.91414, 15.97368, 16.03429, 16.09126, 16.15091, 16.20831, 16.26704, 16.32576, 16.38243, 16.43547, 16.49339, 16.5488, 16.60344, 16.65718, 16.71321, 16.76445, 16.8205, 16.87436, 16.927, 16.97705, 17.03084, 17.0816, 17.1319, 17.18394, 17.23234, 17.2822, 17.33058, 17.38204, 17.42964, 17.47635, 17.52523, 17.57079, 17.61889, 17.66699, 17.71111, 17.75703, 17.8013, 17.84649, 17.89099, 17.93399, 17.9787, 18.02398, 18.06666, 18.10758, 18.14943, 18.18859, 18.23597, 18.27524, 18.31524, 18.35288, 18.39422, 18.43135, 18.47359, 18.50991, 18.54899, 18.58792, 18.62549, 18.6628, 18.70075, 18.73764, 18.77276, 18.80832, 18.8426, 18.87748, 18.91375, 18.95051, 18.98311, 19.01515, 19.04931, 19.08114, 19.11342, 19.1475, 19.18086, 19.21169, 19.24154, 19.27377, 19.30616, 19.33832, 19.367, 19.39638, 19.42439, 19.45338, 19.4916, 19.52028, 19.5463, 19.57274, 19.5997, 19.62729, 19.65703, 19.68112, 19.70886, 19.73457, 19.76013, 19.7864, 19.81312, 19.83962, 19.8632, 19.88678, 19.91159, 19.93924, 19.96327, 19.98656, 20.00949, 20.03395, 20.05728, 20.08171, 20.10283, 20.12461, 20.14697, 20.16721, 20.19076, 20.21256, 20.2331, 20.25148, 20.27287, 20.29232, 20.31233, 20.33611, 20.35665, 20.37633, 20.39356, 20.41459, 20.43284, 20.45168, 20.46953, 20.49013, 20.5074, 20.52513, 20.54128, 20.55937, 20.57494, 20.59229, 20.60831, 20.62823, 20.64686, 20.66215, 20.67842, 20.69401, 20.70969, 20.72814, 20.74172, 20.75548, 20.77282, 20.78749, 20.8006, 20.81588, 20.83165, 20.84472, 20.85776, 20.87512, 20.8893, 20.90468, 20.91558, 20.92848, 20.94267, 20.95483, 20.96783, 21.43397, 20.9946, 21.9169, 21.01688, 21.94514, 21.04086, 21.97679, 21.06474, 22.00501, 21.0871, 22.03576, 21.11301, 22.05861, 21.13361, 22.08682, 21.15366, 22.11164, 21.17516, 22.14002, 21.19661, 22.164, 21.22003, 22.19227, 21.23085, 22.2128, 21.26408, 22.23532, 21.27299, 22.2677, 21.29073, 22.30098, 22.29048, 21.33871, 22.3131, 21.35634, 22.33841, 21.37462, 22.36041, 21.39289, 21.38226, 22.39591, 21.41641, 22.4124, 21.43299, 22.43206, 21.44725, 21.43114, 22.48475, 21.44859, 22.50479, 21.46267, 22.50274, 22.53443, 21.48205, 22.55402, 21.49725, 22.56865, 21.5098, 21.51539, 22.59688, 21.528, 22.61319, 22.60161, 21.57266, 22.62314, 21.5633, 21.5923, 22.64335, 21.60094, 22.65694, 21.58695, 21.61606, 22.67975, 21.62729, 21.60946, 22.70894, 21.64953, 22.71967, 21.656, 21.63735, 22.74293, 21.67279, 21.65285, 22.79216, 21.68411, 21.66181, 22.80982, 21.67621, 21.68103, 22.83539, 21.68425, 22.84617, 22.85479, 21.70155, 22.86597, 22.84777, 21.71176, 22.88069, 22.86413, 21.72417, 22.90486, 22.88857, 21.73994, 22.92506, 22.90715, 21.77352, 21.77888, 22.92102, 21.78287, 21.78373, 22.93597, 21.7902, 21.76845, 22.95921, 21.806, 21.78213, 22.97308, 21.81363, 21.7886, 22.9876, 23.0218, 21.82349, 21.79743, 23.04347, 21.81819, 21.80728, 
1e-100, 0.02988188, 0.06064262, 0.09264608, 0.1255489, 0.159412, 0.194527, 0.2305184, 0.2676164, 0.3058319, 0.3449576, 0.3854015, 0.4268346, 0.469177, 0.5128216, 0.5575024, 0.6030533, 0.6499238, 0.6978475, 0.7466302, 0.7968802, 0.8479218, 0.9000252, 0.9533737, 1.007644, 1.062856, 1.119332, 1.176817, 1.235198, 1.294781, 1.355505, 1.416887, 1.47956, 1.543027, 1.607608, 1.672885, 1.739599, 1.806877, 1.875286, 1.944705, 2.01472, 2.086092, 2.157684, 2.230882, 2.304746, 2.379287, 2.454636, 2.530718, 2.607906, 2.685616, 2.76494, 2.844029, 2.923964, 3.005039, 3.08665, 3.169037, 3.251994, 3.335649, 3.419717, 3.504953, 3.590479, 3.676361, 3.7634, 3.850893, 3.939412, 4.02733, 4.116707, 4.205301, 4.295702, 4.386327, 4.477024, 4.568217, 4.660301, 4.753148, 4.845149, 4.939003, 5.032697, 5.126733, 5.220651, 5.315162, 5.410378, 5.505701, 5.601277, 5.697587, 5.79291, 5.889834, 5.985841, 6.082787, 6.180456, 6.277788, 6.376295, 6.473811, 6.571954, 6.670311, 6.768478, 6.867869, 6.966034, 7.064748, 7.162457, 7.259273, 7.365509, 7.462475, 7.56069, 7.660426, 7.758817, 7.859332, 7.959439, 8.058863, 8.158015, 8.257689, 8.358217, 8.457401, 8.556498, 8.655537, 8.754785, 8.854572, 8.952617, 9.052923, 9.151678, 9.251174, 9.350148, 9.448111, 9.547859, 9.646953, 9.744677, 9.843124, 9.940673, 10.03988, 10.13727, 10.23431, 10.3312, 10.42848, 10.52704, 10.62168, 10.7184, 10.81487, 10.912, 11.00786, 11.10374, 11.19875, 11.29465, 11.38863, 11.484, 11.57728, 11.67084, 11.76569, 11.85816, 11.95191, 12.04424, 12.13327, 12.22903, 12.31976, 12.41253, 12.50319, 12.5946, 12.68476, 12.77374, 12.86444, 12.95354, 13.04402, 13.13008, 13.21944, 13.30737, 13.3945, 13.48127, 13.56862, 13.65698, 13.74137, 13.82644, 13.91212, 13.99511, 14.07888, 14.16436, 14.24636, 14.32825, 14.41218, 14.49519, 14.57532, 14.65786, 14.73784, 14.81971, 14.89898, 14.97804, 15.05781, 15.13485, 15.21294, 15.29257, 15.36966, 15.44693, 15.52207, 15.5995, 15.67498, 15.74934, 15.82463, 15.8984, 15.97221, 16.04473, 16.11951, 16.19414, 16.25967, 16.32486, 16.40095, 16.46988, 16.54, 16.61056, 16.68013, 16.74657, 16.8176, 16.88313, 16.95106, 17.01631, 17.08384, 17.1497, 17.21381, 17.28021, 17.344, 17.40625, 17.47193, 17.53417, 17.5984, 17.65901, 17.72227, 17.78173, 17.84155, 17.90337, 17.96056, 18.02055, 18.07976, 18.13777, 18.19593, 18.25199, 18.30984, 18.36322, 18.42287, 18.4799, 18.53341, 18.58713, 18.64171, 18.69596, 18.74571, 18.80358, 18.85525, 18.90662, 18.95872, 19.00987, 19.05882, 19.11104, 19.16018, 19.21057, 19.26325, 19.30846, 19.35641, 19.40361, 19.4545, 19.50003, 19.54464, 19.59072, 19.63716, 19.67979, 19.72936, 19.77197, 19.81547, 19.85967, 19.90394, 19.94628, 19.9902, 20.03158, 20.07409, 20.11422, 20.15484, 20.19502, 20.23547, 20.27822, 20.31485, 20.35536, 20.39441, 20.43271, 20.46863, 20.50808, 20.54717, 20.58194, 20.6186, 20.65481, 20.69014, 20.72662, 20.76249, 20.79609, 20.83061, 20.86397, 20.89781, 20.93389, 20.9671, 21.00001, 21.03125, 21.06535, 21.09908, 21.12952, 21.15761, 21.1873, 21.2266, 21.25194, 21.28443, 21.3122, 21.34223, 21.37143, 21.40094, 21.42743, 21.4567, 21.48634, 21.51252, 21.53914, 21.56677, 21.5926, 21.62046, 21.64796, 21.67344, 21.69857, 21.72607, 21.75381, 21.77691, 21.80092, 21.82531, 21.85022, 21.87065, 21.89811, 21.92059, 21.94435, 21.9663, 21.98955, 22.0109, 22.03492, 22.05646, 22.07722, 22.0981, 22.12155, 22.14459, 22.16567, 22.18534, 22.20633, 22.22553, 22.24461, 22.26643, 22.28524, 22.30316, 22.32157, 22.34129, 22.35856, 22.37856, 22.39627, 22.41787, 22.43386, 22.45203, 22.47072, 22.48777, 22.504, 22.52217, 22.53939, 22.55435, 22.56931, 22.58558, 22.6025, 22.61594, 22.63237, 22.64633, 22.66484, 22.68087, 22.6967, 22.71135, 22.72277, 22.74234, 22.75396, 22.76735, 22.78325, 22.79682, 23.28355, 22.82403, 23.79226, 22.85026, 23.82279, 22.8771, 23.85528, 22.90489, 23.88573, 22.9288, 23.91422, 22.95266, 23.94346, 22.97487, 23.97388, 22.99742, 23.99723, 23.01891, 24.02409, 23.04246, 24.05522, 23.06534, 24.08236, 23.07326, 24.10859, 23.11897, 24.13028, 23.12659, 24.17204, 23.1434, 24.2144, 24.19315, 23.21122, 24.21923, 23.22976, 24.24196, 23.25059, 24.26588, 23.26829, 23.24171, 24.30186, 23.29406, 24.31907, 23.31117, 24.34031, 23.32858, 23.29785, 24.41362, 23.32047, 24.43667, 23.33596, 24.41962, 24.4688, 23.35682, 24.48841, 23.37176, 24.5036, 23.38362, 23.3903, 24.53331, 23.40728, 24.55674, 24.52907, 23.46831, 24.54904, 23.44245, 23.48736, 24.57319, 23.50174, 24.58981, 23.47356, 23.51555, 24.61395, 23.53199, 23.49443, 24.64251, 23.54839, 24.65827, 23.56321, 23.52453, 24.68173, 23.57885, 23.54219, 24.74714, 23.59205, 23.55575, 24.7702, 23.56967, 23.5725, 24.79534, 23.57909, 24.81117, 24.8171, 23.59911, 24.82942, 24.79762, 23.61006, 24.84876, 24.81733, 23.62467, 24.87202, 24.83808, 23.63757, 24.8913, 24.85579, 23.69306, 23.69757, 24.87343, 23.70251, 23.70592, 24.89416, 23.71691, 23.67782, 24.91386, 23.7306, 23.68977, 24.93267, 23.73726, 23.69629, 24.94826, 24.99695, 23.75025, 23.70767, 25.025, 23.74941, 23.72005, 
1e-100, 0.03599334, 0.07332492, 0.1116383, 0.151215, 0.1919642, 0.2337745, 0.2770361, 0.3213044, 0.3667235, 0.4135402, 0.4613535, 0.5103979, 0.5608066, 0.6121119, 0.6649126, 0.7188879, 0.7738446, 0.8301786, 0.8876882, 0.9462394, 1.005965, 1.066984, 1.12899, 1.192328, 1.256927, 1.322396, 1.388885, 1.456876, 1.525667, 1.595422, 1.66679, 1.739061, 1.811898, 1.886426, 1.961921, 2.038077, 2.11507, 2.193904, 2.273144, 2.352867, 2.434385, 2.516181, 2.599379, 2.683018, 2.768189, 2.853717, 2.940305, 3.027545, 3.115717, 3.20548, 3.294518, 3.384734, 3.476215, 3.567891, 3.660329, 3.753907, 3.847747, 3.942521, 4.03776, 4.133433, 4.230495, 4.326824, 4.424571, 4.523319, 4.621373, 4.721036, 4.81976, 4.920652, 5.021016, 5.122304, 5.223271, 5.325709, 5.427838, 5.530597, 5.633673, 5.736796, 5.840939, 5.945265, 6.048915, 6.153518, 6.258649, 6.363461, 6.470581, 6.575025, 6.682152, 6.787028, 6.894148, 7.00041, 7.106902, 7.213524, 7.320113, 7.42814, 7.534856, 7.642673, 7.749624, 7.857468, 7.964994, 8.072479, 8.177955, 8.293041, 8.397136, 8.504924, 8.612914, 8.721011, 8.828283, 8.935998, 9.043169, 9.151242, 9.259252, 9.366547, 9.474369, 9.581724, 9.689822, 9.795219, 9.902689, 10.00909, 10.11716, 10.22428, 10.32934, 10.43599, 10.54233, 10.64813, 10.75383, 10.85852, 10.96574, 11.06881, 11.1741, 11.28018, 11.3835, 11.48739, 11.59258, 11.69428, 11.79667, 11.90149, 12.00306, 12.10808, 12.20981, 12.3102, 12.41226, 12.51387, 12.61462, 12.71596, 12.81486, 12.91493, 13.0142, 13.11378, 13.21257, 13.30977, 13.40565, 13.50745, 13.60359, 13.70037, 13.79764, 13.89309, 13.98892, 14.08427, 14.18025, 14.27319, 14.36946, 14.46182, 14.55479, 14.64811, 14.73956, 14.83065, 14.9218, 15.01733, 15.10524, 15.19467, 15.28417, 15.37549, 15.46187, 15.55071, 15.63957, 15.72552, 15.8144, 15.89986, 15.98614, 16.07177, 16.15618, 16.24062, 16.32346, 16.40711, 16.49155, 16.57289, 16.65565, 16.73696, 16.81924, 16.8993, 16.97881, 17.05972, 17.13786, 17.21782, 17.2969, 17.37358, 17.45127, 17.52924, 17.60808, 17.68406, 17.75354, 17.82332, 17.89862, 17.97275, 18.0492, 18.12161, 18.19512, 18.26607, 18.33715, 18.40665, 18.47805, 18.54692, 18.61758, 18.68569, 18.75469, 18.8234, 18.88936, 18.95768, 19.0234, 19.08945, 19.15451, 19.21888, 19.28558, 19.34946, 19.41099, 19.47638, 19.53608, 19.60054, 19.66118, 19.72127, 19.78276, 19.8442, 19.90397, 19.96068, 20.02208, 20.08012, 20.1385, 20.19591, 20.25217, 20.30974, 20.36407, 20.42013, 20.47428, 20.53145, 20.58695, 20.63843, 20.69068, 20.74348, 20.7969, 20.8454, 20.90432, 20.95378, 21.00229, 21.05371, 21.10348, 21.15175, 21.20335, 21.25021, 21.29922, 21.34655, 21.39326, 21.43923, 21.48706, 21.53685, 21.57904, 21.62436, 21.6707, 21.71294, 21.75476, 21.80184, 21.84338, 21.8872, 21.93034, 21.9741, 22.01274, 22.05634, 22.09647, 22.13832, 22.17635, 22.21669, 22.25405, 22.29687, 22.3355, 22.3733, 22.4101, 22.45005, 22.48504, 22.52323, 22.56035, 22.59929, 22.63341, 22.66973, 22.70533, 22.74082, 22.77471, 22.80934, 22.84364, 22.87451, 22.90958, 22.94035, 22.97966, 23.01071, 23.04119, 23.07164, 23.10328, 23.13595, 23.16578, 23.19423, 23.22543, 23.25416, 23.28514, 23.31389, 23.34554, 23.3732, 23.40175, 23.43046, 23.45924, 23.4857, 23.51371, 23.54122, 23.56706, 23.59557, 23.62066, 23.64664, 23.67068, 23.69701, 23.71909, 23.74681, 23.77196, 23.79474, 23.81793, 23.84248, 23.86839, 23.89419, 23.91628, 23.93984, 23.96201, 23.98142, 24.00764, 24.02737, 24.0493, 24.07297, 24.09357, 24.11205, 24.13357, 24.15359, 24.17389, 24.19474, 24.21404, 24.23569, 24.25797, 24.27833, 24.29603, 24.31509, 24.33185, 24.35519, 24.37227, 24.38701, 24.40558, 24.42222, 24.44097, 24.45847, 24.47548, 24.49292, 24.51125, 24.52501, 24.54338, 24.56006, 24.57617, 24.59175, 24.60774, 24.62353, 24.63813, 24.65262, 24.66959, 25.17162, 24.69726, 25.69702, 24.72366, 25.73294, 24.75694, 25.76803, 24.78428, 25.79929, 24.81331, 25.83475, 24.83821, 25.86427, 24.8624, 25.89373, 24.88836, 25.92462, 24.91273, 25.95508, 24.94032, 25.98627, 24.96225, 26.01158, 24.97194, 26.03984, 25.0168, 26.06528, 25.02591, 26.12006, 25.0461, 26.17518, 26.13079, 25.14408, 26.16231, 25.1664, 26.18999, 25.18693, 26.21287, 25.20801, 25.15576, 26.25108, 25.23341, 26.27182, 25.25209, 26.29571, 25.27406, 25.22035, 26.39859, 25.23883, 26.41832, 25.25801, 26.38023, 26.45772, 25.28242, 26.47602, 25.29851, 26.49652, 25.31609, 25.32448, 26.53254, 25.33815, 26.55461, 26.50332, 25.43136, 26.52367, 25.37854, 25.45034, 26.54581, 25.46334, 26.56893, 25.41377, 25.48544, 26.59531, 25.50079, 25.43819, 26.62613, 25.52106, 26.64223, 25.53569, 25.47109, 26.66762, 25.55549, 25.49158, 26.76225, 25.57073, 25.50979, 26.79471, 25.52434, 25.52777, 26.81913, 25.53348, 26.83154, 26.8391, 25.55357, 26.85467, 26.79431, 25.56778, 26.87699, 26.81956, 25.58948, 26.9044, 26.84367, 25.60027, 26.92526, 26.86345, 25.68474, 25.68984, 26.87889, 25.6931, 25.69834, 26.90349, 25.71343, 25.64757, 26.92594, 25.72478, 25.65562, 26.94218, 25.73309, 25.66537, 26.96102, 27.04288, 25.74771, 25.67728, 27.07512, 25.75122, 25.69572, 
1e-100, 0.04299332, 0.08719416, 0.1329355, 0.1797817, 0.2279626, 0.277581, 0.3282987, 0.3805319, 0.4339973, 0.4886906, 0.544995, 0.6024242, 0.6609849, 0.7212287, 0.7825809, 0.8451398, 0.9091529, 0.9744011, 1.040862, 1.108833, 1.177828, 1.248018, 1.319719, 1.392363, 1.466238, 1.541528, 1.618029, 1.695561, 1.774372, 1.854185, 1.935202, 2.01732, 2.100651, 2.185076, 2.270248, 2.357023, 2.444696, 2.532941, 2.622639, 2.713141, 2.804647, 2.897163, 2.990725, 3.085148, 3.179875, 3.275912, 3.372952, 3.470293, 3.568767, 3.669161, 3.768979, 3.869571, 3.971116, 4.073541, 4.176392, 4.280227, 4.384399, 4.489041, 4.594748, 4.700605, 4.807433, 4.914959, 5.022787, 5.130967, 5.239302, 5.349727, 5.458421, 5.568727, 5.679538, 5.790146, 5.901269, 6.013105, 6.125684, 6.237378, 6.350373, 6.463433, 6.577022, 6.690182, 6.803828, 6.91797, 7.032714, 7.146356, 7.262274, 7.37662, 7.492066, 7.606589, 7.721786, 7.837208, 7.952862, 8.06825, 8.184418, 8.300668, 8.415361, 8.532452, 8.648137, 8.763663, 8.879529, 8.994957, 9.108495, 9.2314, 9.344747, 9.460891, 9.576402, 9.692047, 9.807261, 9.923412, 10.03855, 10.15479, 10.26932, 10.38419, 10.49929, 10.61445, 10.72882, 10.8429, 10.95695, 11.07067, 11.18567, 11.29929, 11.41188, 11.52492, 11.63842, 11.75045, 11.86294, 11.97506, 12.0872, 12.19898, 12.30993, 12.42204, 12.53194, 12.64139, 12.75331, 12.86407, 12.97028, 13.08106, 13.18851, 13.29865, 13.40753, 13.51332, 13.62053, 13.72846, 13.8347, 13.94084, 14.04466, 14.15142, 14.2557, 14.36077, 14.46484, 14.56708, 14.66952, 14.77571, 14.87861, 14.97912, 15.08192, 15.183, 15.28282, 15.38326, 15.4845, 15.58251, 15.68216, 15.77992, 15.87837, 15.97582, 16.07268, 16.16815, 16.26461, 16.36236, 16.45604, 16.54933, 16.6442, 16.7385, 16.83129, 16.92242, 17.01557, 17.10623, 17.20003, 17.28825, 17.37913, 17.47016, 17.55697, 17.64611, 17.73517, 17.82152, 17.91008, 17.99724, 18.08461, 18.16859, 18.25395, 18.33854, 18.42257, 18.50662, 18.58822, 18.67124, 18.7528, 18.83446, 18.91411, 18.9977, 19.0802, 19.16099, 19.23329, 19.30744, 19.38518, 19.46482, 19.54279, 19.6201, 19.69322, 19.77118, 19.84611, 19.91985, 19.99515, 20.06753, 20.14041, 20.21285, 20.28438, 20.35598, 20.42589, 20.49895, 20.56745, 20.63611, 20.70643, 20.77298, 20.84152, 20.90871, 20.97446, 21.04318, 21.10748, 21.17398, 21.23875, 21.30168, 21.36633, 21.42925, 21.49356, 21.55423, 21.6206, 21.67829, 21.74132, 21.80096, 21.85989, 21.91984, 21.97853, 22.038, 22.09791, 22.15591, 22.21327, 22.26737, 22.32631, 22.38078, 22.43795, 22.49176, 22.54947, 22.60159, 22.65571, 22.71111, 22.7609, 22.81368, 22.86579, 22.91735, 22.96705, 23.01914, 23.06923, 23.12023, 23.17063, 23.21876, 23.26707, 23.3161, 23.3642, 23.41092, 23.45661, 23.50377, 23.54813, 23.59393, 23.64204, 23.68638, 23.72898, 23.77608, 23.81751, 23.86291, 23.90672, 23.94975, 23.9891, 24.03275, 24.07184, 24.11481, 24.1558, 24.19656, 24.23874, 24.27521, 24.31519, 24.35591, 24.39413, 24.43231, 24.47074, 24.50642, 24.54777, 24.58318, 24.61943, 24.65407, 24.68859, 24.72503, 24.76305, 24.79807, 24.83247, 24.86645, 24.90018, 24.93271, 24.9667, 24.99958, 25.03308, 25.06375, 25.09523, 25.12748, 25.16136, 25.19211, 25.22167, 25.25467, 25.28226, 25.31593, 25.34565, 25.37292, 25.40048, 25.43091, 25.45901, 25.48603, 25.51447, 25.54448, 25.56942, 25.59763, 25.62203, 25.64861, 25.67689, 25.70439, 25.73021, 25.75527, 25.78198, 25.80665, 25.83306, 25.856, 25.87944, 25.90206, 25.92779, 25.95154, 25.97366, 25.99445, 26.01978, 26.03891, 26.06275, 26.08603, 26.1099, 26.1336, 26.15156, 26.17674, 26.19458, 26.21424, 26.23645, 26.25819, 26.27678, 26.29699, 26.31402, 26.33294, 26.35195, 26.37114, 26.38952, 26.41337, 26.43055, 26.44972, 26.46622, 26.4863, 26.50124, 26.5235, 26.53793, 26.5559, 26.57091, 26.58638, 26.6053, 27.1245, 26.63827, 27.66756, 26.66733, 27.70559, 26.70238, 27.74165, 26.73143, 27.7763, 26.7626, 27.8102, 26.78882, 27.84348, 26.81611, 27.87607, 26.84355, 27.90653, 26.87359, 27.94479, 26.90196, 27.97962, 26.92836, 28.00919, 26.94197, 28.03562, 26.98919, 28.06549, 26.99955, 28.13708, 27.02205, 28.21341, 28.14144, 27.15365, 28.16989, 27.17898, 28.19938, 27.19921, 28.22911, 27.22341, 27.14146, 28.26796, 27.25448, 28.29394, 27.27984, 28.32062, 27.30224, 27.21672, 28.45914, 27.23842, 28.48281, 27.26003, 28.41052, 28.52085, 27.28749, 28.54476, 27.2998, 28.5661, 27.32584, 27.33448, 28.61058, 27.35257, 28.62869, 28.5451, 27.47346, 28.56906, 27.3937, 27.50182, 28.59681, 27.51701, 28.61635, 27.43238, 27.54181, 28.65069, 27.55777, 27.46233, 28.68212, 27.58274, 28.70057, 27.59941, 27.50049, 28.72953, 27.61916, 27.52157, 28.86425, 27.63941, 27.54291, 28.8944, 27.5558, 27.56684, 28.92564, 27.56895, 28.93894, 28.95027, 27.59484, 28.96452, 28.86999, 27.61327, 28.99279, 28.8959, 27.62662, 29.01849, 28.91887, 27.64402, 29.04331, 28.94171, 27.76595, 27.77273, 28.96559, 27.78443, 27.78996, 28.99211, 27.80252, 27.69972, 29.01513, 27.81596, 27.71345, 29.03358, 27.82836, 27.7202, 29.05143, 29.17446, 27.84463, 27.73927, 29.20867, 27.84424, 27.75684, 
1e-100, 0.05051695, 0.102551, 0.1558262, 0.2107528, 0.2669868, 0.3246681, 0.3839521, 0.4444529, 0.5064781, 0.5700489, 0.6348408, 0.7011523, 0.7689397, 0.8379168, 0.9086746, 0.9805183, 1.053725, 1.128541, 1.204597, 1.281875, 1.36065, 1.440692, 1.521917, 1.604604, 1.688811, 1.773673, 1.860081, 1.947875, 2.036767, 2.126759, 2.218209, 2.310463, 2.403668, 2.498727, 2.594438, 2.691049, 2.78914, 2.888173, 2.988152, 3.088798, 3.191026, 3.293789, 3.397469, 3.502066, 3.608125, 3.714507, 3.821349, 3.929429, 4.038487, 4.149355, 4.258621, 4.369703, 4.481997, 4.59402, 4.707395, 4.820902, 4.935596, 5.050593, 5.166219, 5.282386, 5.399503, 5.515877, 5.633777, 5.75198, 5.870532, 5.990206, 6.108058, 6.228592, 6.349266, 6.46894, 6.589997, 6.710774, 6.832946, 6.954501, 7.076789, 7.197901, 7.320955, 7.44379, 7.566418, 7.688815, 7.813167, 7.93571, 8.060066, 8.182926, 8.306549, 8.430288, 8.553492, 8.677822, 8.801316, 8.925992, 9.050051, 9.174374, 9.297534, 9.422141, 9.545911, 9.669179, 9.792261, 9.916355, 10.03698, 10.16914, 10.29029, 10.41282, 10.53603, 10.65707, 10.78135, 10.90481, 11.02751, 11.1497, 11.27138, 11.3946, 11.51482, 11.63784, 11.7588, 11.88059, 12.00188, 12.12078, 12.24245, 12.36241, 12.48257, 12.60113, 12.71994, 12.84061, 12.95924, 13.07652, 13.19581, 13.31351, 13.43155, 13.54796, 13.66454, 13.78065, 13.89777, 14.01409, 14.12501, 14.24195, 14.35705, 14.47077, 14.58396, 14.69718, 14.8092, 14.92169, 15.03436, 15.14772, 15.25673, 15.36691, 15.47814, 15.58665, 15.6972, 15.80594, 15.91171, 16.02371, 16.13042, 16.23733, 16.34337, 16.45053, 16.55509, 16.66055, 16.76624, 16.86887, 16.97405, 17.0757, 17.17859, 17.28097, 17.38183, 17.48328, 17.58518, 17.68701, 17.78516, 17.88372, 17.98349, 18.08116, 18.17833, 18.27571, 18.37368, 18.46807, 18.56613, 18.65763, 18.75271, 18.84753, 18.93891, 19.03479, 19.12636, 19.21937, 19.30887, 19.40131, 19.49126, 19.58032, 19.66873, 19.76039, 19.84768, 19.93635, 20.024, 20.1111, 20.19588, 20.28152, 20.36745, 20.45354, 20.53919, 20.62317, 20.70021, 20.7793, 20.86295, 20.94542, 21.02733, 21.10937, 21.18939, 21.2669, 21.34569, 21.42145, 21.50454, 21.5803, 21.65815, 21.73574, 21.81241, 21.88584, 21.96047, 22.03464, 22.109, 22.18134, 22.25868, 22.32806, 22.39892, 22.47042, 22.54115, 22.61086, 22.68399, 22.75094, 22.82093, 22.88869, 22.95906, 23.02377, 23.09091, 23.15665, 23.22319, 23.2893, 23.35418, 23.41921, 23.48291, 23.54526, 23.60803, 23.671, 23.73403, 23.79744, 23.85715, 23.918, 23.97796, 24.03771, 24.09588, 24.1549, 24.21597, 24.27552, 24.33158, 24.38919, 24.44419, 24.50157, 24.55657, 24.61231, 24.66576, 24.71975, 24.77422, 24.82887, 24.88306, 24.93564, 24.98795, 25.04021, 25.09059, 25.13931, 25.19231, 25.24354, 25.29285, 25.3422, 25.39042, 25.43955, 25.48658, 25.53706, 25.58259, 25.63037, 25.67548, 25.72248, 25.7695, 25.8135, 25.85882, 25.90249, 25.95013, 25.9948, 26.03869, 26.08065, 26.12394, 26.16422, 26.2085, 26.25003, 26.29031, 26.33185, 26.37177, 26.41251, 26.45056, 26.4915, 26.5269, 26.56937, 26.6086, 26.64703, 26.68276, 26.72034, 26.76027, 26.79235, 26.82919, 26.86663, 26.90137, 26.93789, 26.97488, 27.01071, 27.04218, 27.07736, 27.11188, 27.14602, 27.17934, 27.21113, 27.24441, 27.27636, 27.30728, 27.3394, 27.37228, 27.40034, 27.43104, 27.4613, 27.49402, 27.52201, 27.55155, 27.57862, 27.61309, 27.64248, 27.67149, 27.69952, 27.72332, 27.75252, 27.77715, 27.80707, 27.83062, 27.8594, 27.88561, 27.91107, 27.93582, 27.96305, 27.98796, 28.01214, 28.03675, 28.06263, 28.08831, 28.11525, 28.13664, 28.16102, 28.18263, 28.20309, 28.23069, 28.25154, 28.2711, 28.29321, 28.31614, 28.33371, 28.35807, 28.38005, 28.40223, 28.42186, 28.44147, 28.46614, 28.48439, 28.50308, 28.52411, 28.54395, 28.56053, 28.57875, 28.59805, 28.61476, 28.63502, 29.16485, 28.67254, 29.725, 28.70738, 29.76734, 28.74532, 29.80509, 28.77942, 29.84768, 28.8084, 29.88077, 28.84218, 29.91794, 28.87055, 29.95484, 28.90553, 29.99235, 28.93334, 30.02811, 28.96752, 30.06428, 28.99385, 30.09685, 29.00876, 30.12534, 29.05956, 30.16049, 29.07183, 30.25566, 29.10068, 30.3597, 30.24662, 29.28398, 30.27783, 29.31129, 30.30851, 29.33733, 30.33971, 29.36298, 29.2321, 30.38209, 29.39848, 30.4097, 29.42672, 30.44442, 29.45294, 29.31746, 30.62948, 29.34532, 30.6584, 29.36677, 30.54283, 30.7066, 29.39822, 30.73147, 29.41762, 30.75651, 29.44581, 29.45263, 30.80137, 29.47193, 30.82271, 30.69116, 29.65149, 30.71906, 29.51676, 29.67832, 30.75133, 29.69656, 30.77527, 29.56396, 29.72652, 30.81173, 29.74789, 29.59408, 30.84233, 29.77387, 30.8618, 29.78874, 29.64285, 30.89981, 29.81749, 29.66787, 31.08827, 29.84152, 29.69036, 31.12013, 29.70621, 29.71197, 31.15279, 29.72005, 31.16643, 31.18061, 29.74467, 31.19569, 31.05255, 29.76693, 31.22752, 31.07993, 29.78955, 31.25702, 31.10612, 29.80281, 31.2827, 31.13233, 29.98381, 29.99032, 31.1507, 30.0063, 30.01207, 31.18591, 30.02219, 29.86515, 31.20783, 30.03923, 29.87836, 31.22933, 30.05879, 29.89291, 31.24652, 31.42916, 30.07753, 29.91412, 31.46928, 30.07883, 29.92724, 
1e-100, 0.05854365, 0.118682, 0.1805004, 0.2436605, 0.3085283, 0.3749374, 0.4427007, 0.5122875, 0.5832255, 0.6557542, 0.7299474, 0.805433, 0.8824862, 0.9611754, 1.041059, 1.122539, 1.205568, 1.289782, 1.375731, 1.463088, 1.551576, 1.641592, 1.733148, 1.825525, 1.919715, 2.014922, 2.111423, 2.209313, 2.308505, 2.408936, 2.510235, 2.613037, 2.716644, 2.821859, 2.927888, 3.03521, 3.143054, 3.252339, 3.362845, 3.473867, 3.585998, 3.699048, 3.813436, 3.928329, 4.043818, 4.160661, 4.278085, 4.396065, 4.5149, 4.635858, 4.755883, 4.87657, 4.998071, 5.120621, 5.243734, 5.3673, 5.491046, 5.615911, 5.741596, 5.867007, 5.993266, 6.120398, 6.247405, 6.374814, 6.50257, 6.632016, 6.759024, 6.888804, 7.017405, 7.147296, 7.276229, 7.407272, 7.537525, 7.668208, 7.7993, 7.929638, 8.060568, 8.191909, 8.32347, 8.453977, 8.586229, 8.717576, 8.849701, 8.981562, 9.112869, 9.244848, 9.375541, 9.508842, 9.640138, 9.772437, 9.904047, 10.03494, 10.16655, 10.29856, 10.42947, 10.56021, 10.69138, 10.82261, 10.94956, 11.08993, 11.21704, 11.34695, 11.47659, 11.60637, 11.737, 11.86616, 11.99582, 12.12481, 12.25414, 12.3835, 12.51126, 12.63983, 12.76774, 12.89464, 13.02265, 13.14787, 13.27704, 13.4019, 13.5293, 13.65409, 13.7793, 13.90441, 14.0294, 14.15341, 14.27867, 14.40152, 14.52592, 14.64817, 14.771, 14.89389, 15.01438, 15.13621, 15.25468, 15.37595, 15.49547, 15.61619, 15.73454, 15.85394, 15.97161, 16.09079, 16.20856, 16.3251, 16.44112, 16.55716, 16.67371, 16.78815, 16.90215, 17.0174, 17.12746, 17.24632, 17.35751, 17.46893, 17.58165, 17.69223, 17.80313, 17.91321, 18.02308, 18.13156, 18.24351, 18.34846, 18.45731, 18.56414, 18.67099, 18.77853, 18.88399, 18.9907, 19.09215, 19.19867, 19.30169, 19.4047, 19.50663, 19.61116, 19.71149, 19.81118, 19.91361, 20.01495, 20.11299, 20.21296, 20.31072, 20.41, 20.50628, 20.60534, 20.70105, 20.79776, 20.89267, 20.98615, 21.08069, 21.17478, 21.26714, 21.36066, 21.45451, 21.5448, 21.63649, 21.72716, 21.81942, 21.90921, 21.99837, 22.08778, 22.17241, 22.25519, 22.34249, 22.4311, 22.51724, 22.60291, 22.68873, 22.77263, 22.85708, 22.94093, 23.02565, 23.10673, 23.18859, 23.26979, 23.34896, 23.43227, 23.51198, 23.59103, 23.67096, 23.74774, 23.82769, 23.90211, 23.98063, 24.05741, 24.13303, 24.2099, 24.28449, 24.35749, 24.42938, 24.50339, 24.57892, 24.65164, 24.72167, 24.79164, 24.8629, 24.93649, 25.00544, 25.07348, 25.14253, 25.21484, 25.27969, 25.34775, 25.4158, 25.48122, 25.54743, 25.61232, 25.67905, 25.74443, 25.80647, 25.87024, 25.9368, 25.99895, 26.06217, 26.12175, 26.18378, 26.24529, 26.30454, 26.36367, 26.42568, 26.48508, 26.54435, 26.6046, 26.65947, 26.71816, 26.77412, 26.83296, 26.88713, 26.94068, 26.99647, 27.05035, 27.10705, 27.16268, 27.21277, 27.26683, 27.31944, 27.37161, 27.42547, 27.4771, 27.52711, 27.57695, 27.62902, 27.67711, 27.7285, 27.77859, 27.82609, 27.87575, 27.92209, 27.9711, 28.01816, 28.06252, 28.11081, 28.15655, 28.19901, 28.2495, 28.29072, 28.33338, 28.37934, 28.42115, 28.46208, 28.5053, 28.55373, 28.59503, 28.63389, 28.676, 28.71642, 28.75722, 28.79574, 28.83924, 28.8774, 28.91755, 28.95491, 28.9962, 29.03174, 29.07146, 29.10749, 29.14426, 29.18635, 29.21852, 29.25349, 29.29193, 29.32455, 29.3581, 29.39377, 29.4304, 29.46248, 29.49467, 29.53186, 29.5626, 29.59752, 29.6318, 29.66417, 29.69665, 29.72422, 29.76032, 29.78653, 29.82048, 29.8493, 29.87884, 29.90571, 29.93714, 29.96801, 29.99386, 30.02048, 30.05372, 30.08001, 30.11039, 30.13919, 30.16763, 30.19265, 30.21663, 30.2476, 30.27137, 30.29739, 30.32446, 30.3497, 30.37164, 30.39711, 30.42173, 30.44489, 30.46817, 30.49447, 30.51822, 30.54553, 30.56578, 30.5905, 30.61296, 30.63337, 30.65642, 30.68044, 30.69941, 30.71897, 30.74141, 30.76059, 30.78095, 30.80259, 31.34517, 30.84664, 31.90706, 30.88591, 31.95403, 30.92213, 31.99617, 30.96179, 32.04041, 30.99587, 32.07871, 31.03287, 32.1212, 31.06527, 32.16306, 31.10314, 32.20832, 31.13813, 32.24411, 31.17336, 32.28022, 31.20192, 32.31665, 31.21798, 32.34955, 31.27627, 32.38606, 31.28905, 32.51654, 31.32215, 32.65488, 32.47702, 31.56594, 32.51709, 31.59631, 32.55361, 31.6262, 32.58217, 31.65341, 31.47251, 32.63408, 31.69864, 32.66785, 31.72851, 32.70104, 31.7614, 31.5708, 32.95577, 31.59674, 32.9837, 31.62187, 32.80987, 33.03372, 31.65528, 33.06119, 31.68496, 33.09478, 31.71341, 31.71627, 33.13924, 31.74293, 33.16802, 32.97244, 31.98519, 33.00816, 31.79484, 32.01633, 33.03941, 32.0386, 33.07141, 31.8515, 32.07653, 33.10939, 32.09618, 31.88493, 33.14436, 32.12745, 33.16732, 32.1428, 31.93318, 33.19971, 32.17613, 31.96405, 33.46339, 32.20385, 31.98737, 33.50044, 32.01181, 32.01281, 33.53109, 32.02124, 33.55411, 33.56107, 32.0533, 33.58674, 33.3781, 32.07875, 33.6156, 33.40123, 32.10161, 33.64866, 33.43426, 32.11949, 33.68295, 33.46191, 32.36999, 32.37722, 33.48709, 32.39357, 32.40292, 33.51475, 32.41718, 32.18963, 33.54478, 32.43336, 32.20266, 33.5681, 32.45196, 32.22281, 33.59333, 33.8423, 32.48078, 32.2419, 33.88399, 32.47983, 32.2626, 
1e-100, 0.06709994, 0.1358274, 0.2060685, 0.2782222, 0.351826, 0.4272287, 0.5042918, 0.5827771, 0.6631215, 0.7449923, 0.8282841, 0.9134311, 0.9999861, 1.088049, 1.177989, 1.269121, 1.361699, 1.456155, 1.55177, 1.648732, 1.747457, 1.847367, 1.948603, 2.05147, 2.155741, 2.26075, 2.367716, 2.475643, 2.584621, 2.695144, 2.806769, 2.919346, 3.03312, 3.148619, 3.264668, 3.382051, 3.50017, 3.619591, 3.739833, 3.860544, 3.983065, 4.106076, 4.230002, 4.355048, 4.480476, 4.606902, 4.733836, 4.862011, 4.990619, 5.120902, 5.249767, 5.380542, 5.512104, 5.643702, 5.776298, 5.909013, 6.042329, 6.176198, 6.31096, 6.445978, 6.581194, 6.716471, 6.852317, 6.989089, 7.125814, 7.2637, 7.399746, 7.537515, 7.675692, 7.813355, 7.951225, 8.090366, 8.229048, 8.367976, 8.50697, 8.645688, 8.784339, 8.92441, 9.063557, 9.201827, 9.341783, 9.481234, 9.621421, 9.760407, 9.899894, 10.03904, 10.17793, 10.31779, 10.45613, 10.59551, 10.73429, 10.87343, 11.0121, 11.15098, 11.28954, 11.42849, 11.56627, 11.70327, 11.83775, 11.98435, 12.11847, 12.25537, 12.39182, 12.52846, 12.66497, 12.80167, 12.93814, 13.07335, 13.20827, 13.34375, 13.47928, 13.61442, 13.74836, 13.88102, 14.01455, 14.14756, 14.28229, 14.41457, 14.54638, 14.67822, 14.80984, 14.94037, 15.07232, 15.20268, 15.33319, 15.46133, 15.59159, 15.72163, 15.84939, 15.97649, 16.10493, 16.2321, 16.35725, 16.48444, 16.61009, 16.7369, 16.86158, 16.98554, 17.10996, 17.23534, 17.35837, 17.48068, 17.60182, 17.72492, 17.8447, 17.9671, 18.08832, 18.20724, 18.32484, 18.44783, 18.56424, 18.68323, 18.80255, 18.91852, 19.03499, 19.15031, 19.26584, 19.38307, 19.49766, 19.6111, 19.72434, 19.83802, 19.95159, 20.0633, 20.17494, 20.28825, 20.3971, 20.50806, 20.61784, 20.72649, 20.83422, 20.9428, 21.05104, 21.15606, 21.26559, 21.37055, 21.4764, 21.58005, 21.68597, 21.79163, 21.89399, 21.99565, 22.10052, 22.20223, 22.30356, 22.40449, 22.50639, 22.60676, 22.70568, 22.80656, 22.9024, 23.00181, 23.09824, 23.19634, 23.29296, 23.38877, 23.48644, 23.58087, 23.67152, 23.76257, 23.85661, 23.95051, 24.04302, 24.13319, 24.22472, 24.31548, 24.40629, 24.49902, 24.58776, 24.67597, 24.76417, 24.85145, 24.93824, 25.0269, 25.11184, 25.19845, 25.28273, 25.36741, 25.45216, 25.53549, 25.62019, 25.70186, 25.78158, 25.86811, 25.94705, 26.0282, 26.10928, 26.18807, 26.26706, 26.34529, 26.42011, 26.49906, 26.57955, 26.65663, 26.73254, 26.80749, 26.88204, 26.95719, 27.03239, 27.10674, 27.17954, 27.24966, 27.3251, 27.39545, 27.4655, 27.53749, 27.6103, 27.67647, 27.74934, 27.81939, 27.88491, 27.95422, 28.021, 28.08925, 28.15534, 28.21889, 28.28448, 28.35077, 28.4151, 28.48192, 28.54164, 28.60749, 28.6693, 28.73217, 28.79371, 28.85738, 28.91574, 28.97528, 29.03714, 29.09623, 29.15692, 29.21358, 29.27247, 29.32973, 29.38873, 29.44493, 29.50243, 29.561, 29.61491, 29.67029, 29.72213, 29.78219, 29.83382, 29.88708, 29.93904, 29.99209, 30.04275, 30.09865, 30.14744, 30.19626, 30.2486, 30.29902, 30.34798, 30.39638, 30.44592, 30.49306, 30.5398, 30.58758, 30.63806, 30.68609, 30.73072, 30.77778, 30.82132, 30.8678, 30.91269, 30.95853, 31.00259, 31.04439, 31.09124, 31.1288, 31.17644, 31.21561, 31.2557, 31.29883, 31.34041, 31.37993, 31.41984, 31.46016, 31.49909, 31.53864, 31.57757, 31.61402, 31.65227, 31.69387, 31.72968, 31.76791, 31.79763, 31.843, 31.87716, 31.91389, 31.94782, 31.98516, 32.01553, 32.04844, 32.08543, 32.11862, 32.1488, 32.18582, 32.21484, 32.2483, 32.27981, 32.31387, 32.34781, 32.37207, 32.40912, 32.43871, 32.4707, 32.50091, 32.52995, 32.55629, 32.58613, 32.61551, 32.64303, 32.6676, 32.69933, 32.72219, 32.75192, 32.77828, 32.80887, 32.83564, 32.85736, 32.88594, 32.91186, 32.93304, 32.95866, 32.98744, 33.0116, 33.0327, 33.05581, 33.08229, 33.10603, 33.12512, 33.15344, 33.69721, 33.20168, 34.27214, 33.2412, 34.32055, 33.28172, 34.3715, 33.32802, 34.41171, 33.36249, 34.45943, 33.40513, 34.50327, 33.44676, 34.55156, 33.485, 34.58804, 33.5243, 34.63432, 33.55782, 34.67801, 33.59395, 34.71301, 33.61235, 34.75372, 33.67653, 34.79779, 33.69268, 34.96818, 33.72953, 35.15026, 34.89585, 34.05839, 34.93633, 34.0897, 34.97505, 34.128, 35.00935, 34.15541, 33.89771, 35.0679, 34.21199, 35.10626, 34.24489, 35.14308, 34.2771, 34.00348, 35.48296, 34.03867, 35.51609, 34.06685, 35.26601, 35.56899, 34.10493, 35.60719, 34.13617, 35.63932, 34.16515, 34.17755, 35.69296, 34.20134, 35.72144, 35.44167, 34.5319, 35.47552, 34.26398, 34.57275, 35.51994, 34.59824, 35.54317, 34.32784, 34.63746, 35.58837, 34.65795, 34.3632, 35.63247, 34.69207, 35.65666, 34.71767, 34.42207, 35.698, 34.74938, 34.45539, 36.04703, 34.78079, 34.47935, 36.08934, 34.5041, 34.50287, 36.12156, 34.51821, 36.14901, 36.16125, 34.55093, 36.18446, 35.88588, 34.58241, 36.2183, 35.91507, 34.60356, 36.25287, 35.94527, 34.62598, 36.28772, 35.97639, 34.97214, 34.98856, 36.00817, 35.00146, 35.00428, 36.03866, 35.02309, 34.70412, 36.06453, 35.04445, 34.72535, 36.09377, 35.06443, 34.74028, 36.12252, 36.47088, 35.09695, 34.76513, 36.51428, 35.09806, 34.78438, 
1e-100, 0.07563593, 0.1532063, 0.2325183, 0.3133745, 0.3962508, 0.4806957, 0.5667837, 0.6547589, 0.7442214, 0.835534, 0.9285282, 1.022966, 1.119222, 1.216967, 1.316183, 1.417141, 1.519565, 1.623363, 1.72896, 1.836077, 1.944182, 2.054115, 2.165338, 2.277693, 2.391722, 2.506883, 2.623225, 2.741192, 2.86025, 2.980468, 3.101588, 3.224438, 3.348041, 3.472685, 3.598604, 3.725689, 3.853522, 3.982319, 4.112508, 4.243134, 4.374855, 4.507069, 4.640726, 4.775145, 4.909573, 5.045416, 5.181928, 5.318809, 5.456644, 5.596279, 5.734515, 5.873887, 6.01423, 6.155104, 6.296407, 6.437902, 6.579248, 6.72219, 6.86541, 7.008934, 7.15266, 7.296972, 7.441572, 7.586443, 7.731583, 7.877949, 8.021532, 8.167291, 8.313307, 8.459236, 8.605179, 8.751731, 8.89847, 9.045066, 9.192235, 9.338837, 9.485074, 9.632008, 9.779389, 9.926359, 10.07359, 10.21976, 10.36695, 10.51312, 10.66071, 10.80615, 10.9523, 11.10004, 11.24565, 11.39184, 11.5385, 11.68418, 11.82973, 11.97531, 12.1212, 12.26584, 12.4101, 12.55454, 12.69718, 12.85026, 12.99107, 13.13412, 13.27848, 13.42218, 13.56462, 13.70813, 13.85003, 13.99316, 14.13516, 14.27796, 14.41779, 14.55997, 14.69995, 14.84067, 14.98025, 15.11991, 15.26071, 15.4002, 15.54003, 15.67683, 15.81573, 15.95285, 16.09197, 16.22779, 16.36419, 16.50123, 16.637, 16.77318, 16.90635, 17.04234, 17.17732, 17.31233, 17.44079, 17.57635, 17.70932, 17.84099, 17.97517, 18.1059, 18.23623, 18.36753, 18.49662, 18.6271, 18.75611, 18.88506, 19.01291, 19.14036, 19.26849, 19.39492, 19.51821, 19.6474, 19.7747, 19.89902, 20.02321, 20.14741, 20.2706, 20.39499, 20.5188, 20.64084, 20.76232, 20.88319, 21.00469, 21.1259, 21.2454, 21.36431, 21.48302, 21.60327, 21.72066, 21.83807, 21.95515, 22.07187, 22.18799, 22.3031, 22.41895, 22.53161, 22.64877, 22.7612, 22.87521, 22.9877, 23.10009, 23.21233, 23.32288, 23.43252, 23.54282, 23.6538, 23.76144, 23.87166, 23.97928, 24.08798, 24.19324, 24.30494, 24.40856, 24.51441, 24.61854, 24.72412, 24.82878, 24.93273, 25.0389, 25.14027, 25.23835, 25.34026, 25.43889, 25.53944, 25.64031, 25.74266, 25.84223, 25.94171, 26.03942, 26.13697, 26.23397, 26.32884, 26.42418, 26.52283, 26.61955, 26.71256, 26.80724, 26.89971, 26.99168, 27.08405, 27.1773, 27.27031, 27.36115, 27.45066, 27.53995, 27.63095, 27.71987, 27.80501, 27.8982, 27.98427, 28.0695, 28.15526, 28.24322, 28.32421, 28.41691, 28.49717, 28.58226, 28.66456, 28.7514, 28.83041, 28.91357, 28.99364, 29.07449, 29.15649, 29.23507, 29.31456, 29.39164, 29.47033, 29.54913, 29.6269, 29.70578, 29.78107, 29.85372, 29.93223, 30.00464, 30.07958, 30.15268, 30.22973, 30.29913, 30.37083, 30.44327, 30.51454, 30.58492, 30.65686, 30.72599, 30.79349, 30.86331, 30.92898, 30.99718, 31.06665, 31.13507, 31.1991, 31.2658, 31.32969, 31.39844, 31.4642, 31.52739, 31.59061, 31.64927, 31.71693, 31.7764, 31.8381, 31.89991, 31.95994, 32.02137, 32.08058, 32.14181, 32.19851, 32.25473, 32.314, 32.36941, 32.42696, 32.48554, 32.53838, 32.59239, 32.64898, 32.70126, 32.75538, 32.80716, 32.86639, 32.91713, 32.96996, 33.02151, 33.07355, 33.12614, 33.17468, 33.22659, 33.27513, 33.32761, 33.37381, 33.42691, 33.46959, 33.51871, 33.56444, 33.6105, 33.66096, 33.70302, 33.74712, 33.79536, 33.83845, 33.87982, 33.92408, 33.96965, 34.0121, 34.0528, 34.10029, 34.14022, 34.18149, 34.22248, 34.2633, 34.30331, 34.34154, 34.38283, 34.42217, 34.45968, 34.49724, 34.53429, 34.57034, 34.60953, 34.64763, 34.68384, 34.71599, 34.75906, 34.79353, 34.82999, 34.86191, 34.89735, 34.9297, 34.96071, 34.99992, 35.02937, 35.06311, 35.09648, 35.12486, 35.15719, 35.18965, 35.2218, 35.25166, 35.28013, 35.3119, 35.34282, 35.37564, 35.40232, 35.43082, 35.45918, 35.48421, 35.51561, 35.5433, 35.56714, 35.59274, 35.62236, 35.64648, 35.67393, 35.70133, 35.73002, 36.27149, 35.77166, 36.8497, 35.82113, 36.90117, 35.87457, 36.95837, 35.918, 37.00496, 35.96379, 37.05714, 36.01075, 37.10253, 36.04955, 37.15754, 36.09764, 37.20498, 36.13766, 37.24963, 36.17904, 37.29898, 36.21982, 37.33097, 36.23634, 37.38087, 36.30286, 37.42872, 36.33118, 37.64957, 36.36584, 37.87195, 37.53924, 36.79462, 37.5796, 36.83181, 37.6241, 36.8751, 37.66511, 36.91429, 36.55997, 37.72817, 36.9678, 37.76955, 37.01236, 37.81073, 37.04411, 36.6784, 38.24713, 36.70909, 38.28791, 36.75022, 37.94265, 38.34878, 36.78699, 38.38845, 36.82567, 38.42629, 36.85867, 36.86846, 38.48178, 36.89361, 38.51412, 38.14183, 37.33461, 38.17489, 36.96255, 37.37932, 38.22269, 37.40918, 38.2573, 37.03401, 37.45019, 38.29636, 37.47837, 37.078, 38.34536, 37.51201, 38.36959, 37.53663, 37.14622, 38.41932, 37.58478, 37.18106, 38.8776, 37.61125, 37.21425, 38.92071, 37.23547, 37.2347, 38.95948, 37.24965, 38.98754, 39.00491, 37.2952, 39.02912, 38.6218, 37.31885, 39.06803, 38.65478, 37.35186, 39.10813, 38.69057, 37.37259, 39.14006, 38.7239, 37.83692, 37.84703, 38.75516, 37.86237, 37.87069, 38.79089, 37.89109, 37.46273, 38.81894, 37.91407, 37.47726, 38.84551, 37.9401, 37.50365, 38.87837, 39.33927, 37.96583, 37.53845, 39.39317, 37.97326, 37.5492, 
1e-100, 0.08464073, 0.1709566, 0.2590705, 0.34913, 0.4407847, 0.5344923, 0.6298697, 0.7268087, 0.8258221, 0.9263608, 1.028452, 1.132459, 1.237992, 1.345145, 1.45411, 1.564365, 1.676288, 1.789819, 1.904601, 2.021172, 2.138933, 2.258221, 2.37877, 2.500902, 2.624577, 2.748924, 2.875055, 3.002322, 3.130428, 3.260316, 3.390879, 3.522768, 3.655535, 3.790153, 3.925008, 4.061036, 4.198316, 4.336295, 4.475222, 4.614775, 4.755706, 4.896727, 5.039273, 5.182393, 5.326138, 5.470531, 5.615698, 5.761637, 5.907943, 6.05613, 6.202716, 6.351067, 6.499459, 6.648418, 6.798199, 6.948504, 7.098328, 7.248803, 7.400416, 7.552548, 7.703916, 7.855899, 8.008286, 8.161661, 8.314329, 8.467674, 8.619613, 8.774082, 8.927559, 9.081215, 9.234405, 9.389073, 9.543062, 9.697985, 9.851731, 10.00584, 10.16073, 10.3154, 10.46869, 10.6224, 10.77715, 10.93142, 11.08679, 11.24001, 11.39438, 11.54722, 11.70134, 11.85443, 12.0084, 12.16218, 12.31541, 12.46865, 12.62075, 12.77438, 12.92767, 13.07914, 13.23168, 13.38351, 13.53177, 13.69243, 13.8403, 13.99067, 14.14203, 14.29244, 14.44222, 14.59292, 14.74272, 14.89316, 15.04237, 15.1923, 15.34107, 15.48939, 15.63654, 15.78488, 15.93239, 16.07932, 16.22817, 16.37397, 16.52069, 16.66577, 16.81148, 16.95681, 17.10276, 17.246, 17.39063, 17.53541, 17.67963, 17.82305, 17.96542, 18.10856, 18.24877, 18.39236, 18.53232, 18.6745, 18.81534, 18.9556, 19.09637, 19.2354, 19.37373, 19.51265, 19.65177, 19.78919, 19.92514, 20.0615, 20.20053, 20.33604, 20.47295, 20.60772, 20.73946, 20.87836, 21.01241, 21.14688, 21.27826, 21.41373, 21.54487, 21.67702, 21.80909, 21.93955, 22.07149, 22.20136, 22.33159, 22.4594, 22.58699, 22.71737, 22.84523, 22.97529, 23.10181, 23.2276, 23.35184, 23.48036, 23.60472, 23.73017, 23.85265, 23.97814, 24.10244, 24.22338, 24.34661, 24.47103, 24.59065, 24.7129, 24.83244, 24.95384, 25.07195, 25.19318, 25.31364, 25.43027, 25.54954, 25.66841, 25.78272, 25.89879, 26.01338, 26.12899, 26.24505, 26.35986, 26.47497, 26.58957, 26.7034, 26.81309, 26.92294, 27.03239, 27.14238, 27.25242, 27.36436, 27.47477, 27.58404, 27.69109, 27.80251, 27.90558, 28.01454, 28.11955, 28.22819, 28.33024, 28.43796, 28.54015, 28.64336, 28.74637, 28.8518, 28.95197, 29.05539, 29.15558, 29.25561, 29.35589, 29.45716, 29.55506, 29.65218, 29.75006, 29.84919, 29.9443, 30.04163, 30.13654, 30.23261, 30.32577, 30.42234, 30.51595, 30.60599, 30.69942, 30.79539, 30.88379, 30.97286, 31.0651, 31.15168, 31.24647, 31.33072, 31.42099, 31.50901, 31.59942, 31.68137, 31.76808, 31.85466, 31.93874, 32.02489, 32.10832, 32.19108, 32.2731, 32.35454, 32.43753, 32.51768, 32.59645, 32.67961, 32.75576, 32.83686, 32.9139, 32.99353, 33.07075, 33.14697, 33.22363, 33.2963, 33.3754, 33.45127, 33.5244, 33.59558, 33.67145, 33.74042, 33.815, 33.8851, 33.95978, 34.0295, 34.09771, 34.16729, 34.23193, 34.30923, 34.37158, 34.43748, 34.50509, 34.57235, 34.63378, 34.70396, 34.76482, 34.83091, 34.89469, 34.95501, 35.01577, 35.07943, 35.14199, 35.20171, 35.26225, 35.31908, 35.3806, 35.44148, 35.50081, 35.5609, 35.61485, 35.67469, 35.7327, 35.78716, 35.83878, 35.89556, 35.95186, 36.0015, 36.06118, 36.1086, 36.16425, 36.21644, 36.26606, 36.31472, 36.36715, 36.41986, 36.46804, 36.51606, 36.56654, 36.61431, 36.66577, 36.71386, 36.75868, 36.80648, 36.84627, 36.90089, 36.94467, 36.98519, 37.03211, 37.07887, 37.11913, 37.15999, 37.21091, 37.24773, 37.28704, 37.33045, 37.36889, 37.41225, 37.45428, 37.49362, 37.53824, 37.56868, 37.61095, 37.65263, 37.68664, 37.72409, 37.76254, 37.79359, 37.8332, 37.8719, 37.90218, 37.93513, 37.97921, 38.00893, 38.04574, 38.07843, 38.11114, 38.14651, 38.17328, 38.21499, 38.24096, 38.26995, 38.30342, 38.33382, 38.36291, 38.39573, 38.42405, 38.45574, 38.48754, 38.50868, 38.54057, 38.57169, 39.11576, 38.62824, 39.69437, 38.67674, 39.75131, 38.73627, 39.80962, 38.77917, 39.86029, 38.83065, 39.91922, 38.8858, 39.97373, 38.92999, 40.02806, 38.97196, 40.08137, 39.023, 40.13044, 39.0694, 40.178, 39.1094, 40.22004, 39.13754, 40.27282, 39.20966, 40.32014, 39.23629, 40.6061, 39.28076, 40.89546, 40.44067, 39.82382, 40.4913, 39.8681, 40.53893, 39.91067, 40.57695, 39.95803, 39.49336, 40.65916, 40.0187, 40.69622, 40.06454, 40.73752, 40.09976, 39.624, 41.30442, 39.66138, 41.34699, 39.69466, 40.88458, 41.4139, 39.7498, 41.45625, 39.78531, 41.49406, 39.82501, 39.83186, 41.56216, 39.86135, 41.59506, 41.09991, 40.43004, 41.14437, 39.94235, 40.47838, 41.18944, 40.51055, 41.23369, 40.02282, 40.55608, 41.27715, 40.58964, 40.06535, 41.32245, 40.63035, 41.35887, 40.66208, 40.13742, 41.40265, 40.70832, 40.18276, 41.9923, 40.74478, 40.21208, 42.03585, 40.23463, 40.24349, 42.08369, 40.26019, 42.11564, 42.1266, 40.30542, 42.16323, 41.62704, 40.3354, 42.20094, 41.65731, 40.35848, 42.24562, 41.69891, 40.38711, 42.28312, 41.72664, 40.98747, 41.00457, 41.77045, 41.01883, 41.02755, 41.80364, 41.04597, 40.48067, 41.84002, 41.08179, 40.5137, 41.8689, 41.10708, 40.53962, 41.90249, 42.49746, 41.14162, 40.56641, 42.5525, 41.14329, 40.58111, 
1e-100, 0.09308848, 0.1883515, 0.2853451, 0.3840711, 0.4849154, 0.5873808, 0.6916672, 0.7978313, 0.9055696, 1.015317, 1.126672, 1.239541, 1.354326, 1.470635, 1.588412, 1.707993, 1.829012, 1.95151, 2.075748, 2.201285, 2.328033, 2.456595, 2.586339, 2.717116, 2.849835, 2.983317, 3.118076, 3.254679, 3.392067, 3.530347, 3.66986, 3.810687, 3.952326, 4.095227, 4.239122, 4.383936, 4.529501, 4.676262, 4.824097, 4.971741, 5.121317, 5.271178, 5.422247, 5.573463, 5.725343, 5.878167, 6.031564, 6.185225, 6.340116, 6.496592, 6.651541, 6.807196, 6.964547, 7.12138, 7.279118, 7.437209, 7.594829, 7.753738, 7.912629, 8.072144, 8.232145, 8.392033, 8.552409, 8.712965, 8.873794, 9.035358, 9.195097, 9.35599, 9.517938, 9.679493, 9.840425, 10.00238, 10.16396, 10.32604, 10.48859, 10.65009, 10.81116, 10.97379, 11.13646, 11.29742, 11.45945, 11.62097, 11.78363, 11.94513, 12.10683, 12.26823, 12.4298, 12.59166, 12.75285, 12.91367, 13.07486, 13.23522, 13.39682, 13.55811, 13.71863, 13.87759, 14.03875, 14.1985, 14.35505, 14.52267, 14.68023, 14.83977, 14.99895, 15.15625, 15.31495, 15.47418, 15.63325, 15.79064, 15.94861, 16.10584, 16.26366, 16.42093, 16.57838, 16.73321, 16.89123, 17.04643, 17.20282, 17.35779, 17.51428, 17.66883, 17.82303, 17.97632, 18.13224, 18.2851, 18.44073, 18.59213, 18.74675, 18.90015, 19.05065, 19.20372, 19.35612, 19.50853, 19.65743, 19.80801, 19.95797, 20.11074, 20.25909, 20.40796, 20.55687, 20.70531, 20.85441, 21.00259, 21.15096, 21.29777, 21.44513, 21.59089, 21.73621, 21.88363, 22.02828, 22.17656, 22.31987, 22.4651, 22.60935, 22.75207, 22.89642, 23.03767, 23.18206, 23.32311, 23.46504, 23.60436, 23.74708, 23.88761, 24.02675, 24.16781, 24.30553, 24.44674, 24.58389, 24.7226, 24.86107, 24.9984, 25.13262, 25.26944, 25.4053, 25.53972, 25.6778, 25.81058, 25.94425, 26.07766, 26.21081, 26.34272, 26.47617, 26.60846, 26.7399, 26.86797, 27.00361, 27.13064, 27.26018, 27.38887, 27.51588, 27.64586, 27.7735, 27.90129, 28.02786, 28.15226, 28.28054, 28.40565, 28.53208, 28.65268, 28.77726, 28.89562, 29.01799, 29.14093, 29.26431, 29.38423, 29.50506, 29.62633, 29.74314, 29.86177, 29.98084, 30.09729, 30.21954, 30.33472, 30.44979, 30.56585, 30.67682, 30.79528, 30.90743, 31.0231, 31.13742, 31.25094, 31.35891, 31.47065, 31.57956, 31.69138, 31.80059, 31.91068, 32.01631, 32.12276, 32.23237, 32.3412, 32.44218, 32.5506, 32.65454, 32.75906, 32.86335, 32.96663, 33.06823, 33.17084, 33.26887, 33.37144, 33.4725, 33.57221, 33.67166, 33.76574, 33.8646, 33.96538, 34.05986, 34.15321, 34.25462, 34.34397, 34.44177, 34.5326, 34.62733, 34.71973, 34.81544, 34.90312, 34.9935, 35.0809, 35.17282, 35.26046, 35.34911, 35.43764, 35.52108, 35.6126, 35.69463, 35.78004, 35.86488, 35.95015, 36.03116, 36.11465, 36.19885, 36.28175, 36.36007, 36.44073, 36.52033, 36.59339, 36.67899, 36.7564, 36.83094, 36.90801, 36.98508, 37.0628, 37.13774, 37.21054, 37.2824, 37.35678, 37.4311, 37.49715, 37.57022, 37.64343, 37.70982, 37.7761, 37.85162, 37.91629, 37.98552, 38.05118, 38.12281, 38.18582, 38.25683, 38.3206, 38.38494, 38.44835, 38.51062, 38.57716, 38.63494, 38.70211, 38.7582, 38.82207, 38.88067, 38.93909, 38.99465, 39.0552, 39.11687, 39.17021, 39.22347, 39.28572, 39.33706, 39.39477, 39.45157, 39.50657, 39.56381, 39.60789, 39.66927, 39.71893, 39.76918, 39.82136, 39.87379, 39.92397, 39.97079, 40.02227, 40.073, 40.1171, 40.16593, 40.21279, 40.26062, 40.30527, 40.35348, 40.39781, 40.43725, 40.49341, 40.53206, 40.57353, 40.61686, 40.66419, 40.70279, 40.74208, 40.78747, 40.83122, 40.86751, 40.91113, 40.94473, 40.98481, 41.02671, 41.06555, 41.10569, 41.13488, 41.18162, 41.21585, 41.25549, 41.28994, 41.32254, 41.35631, 41.39444, 41.42908, 41.46027, 41.48854, 41.52306, 41.55524, 41.59463, 41.62704, 41.65639, 41.68763, 41.72056, 42.25636, 41.78339, 42.83266, 41.84144, 42.89602, 41.89392, 42.95881, 41.95398, 43.0219, 42.00783, 43.08055, 42.0639, 43.14234, 42.12125, 43.19947, 42.16706, 43.25502, 42.21237, 43.30656, 42.26781, 43.36146, 42.30308, 43.40535, 42.34294, 43.46684, 42.41972, 43.5192, 42.45146, 43.86763, 42.48941, 44.23188, 43.64635, 43.18534, 43.69798, 43.23303, 43.75071, 43.27697, 43.8023, 43.32439, 42.73038, 43.87055, 43.40333, 43.92444, 43.45101, 43.96989, 43.49165, 42.87083, 44.67846, 42.92061, 44.72332, 42.95708, 44.1265, 44.80036, 43.02001, 44.85307, 43.05011, 44.88457, 43.09552, 43.10351, 44.95631, 43.14063, 45.00782, 44.36657, 43.84953, 44.41083, 43.22915, 43.91619, 44.45845, 43.94974, 44.503, 43.30627, 43.99651, 44.55177, 44.03369, 43.36067, 44.59849, 44.078, 44.63757, 44.12371, 43.44393, 44.68986, 44.1662, 43.48194, 45.43001, 44.20516, 43.52293, 45.48216, 43.54311, 43.55372, 45.53333, 43.58021, 45.56869, 45.58663, 43.62397, 45.61588, 44.92289, 43.65586, 45.66077, 44.96533, 43.68243, 45.70706, 44.99382, 43.71319, 45.75286, 45.04393, 44.47328, 44.49249, 45.07659, 44.50942, 44.52211, 45.12161, 44.53928, 43.81879, 45.14823, 44.56571, 43.85167, 45.18875, 44.60672, 43.87935, 45.21732, 45.9819, 44.65496, 43.90678, 46.03467, 44.64847, 43.9253, 
1e-100, 0.1018574, 0.205439, 0.3110168, 0.4184772, 0.5276906, 0.6390297, 0.7519897, 0.8666715, 0.9833654, 1.10166, 1.221648, 1.343391, 1.466669, 1.591791, 1.718453, 1.846532, 1.976362, 2.10758, 2.240243, 2.374501, 2.510034, 2.646862, 2.785332, 2.924979, 3.066102, 3.208114, 3.351612, 3.496207, 3.642006, 3.788933, 3.937027, 4.08592, 4.235994, 4.387436, 4.539442, 4.692467, 4.846508, 5.001363, 5.156842, 5.313267, 5.470952, 5.62861, 5.78778, 5.947319, 6.107471, 6.267884, 6.429765, 6.591641, 6.754243, 6.918448, 7.081451, 7.245565, 7.410102, 7.575033, 7.740811, 7.906844, 8.072862, 8.239197, 8.406189, 8.574009, 8.741684, 8.909792, 9.077557, 9.246458, 9.414896, 9.584355, 9.751621, 9.921914, 10.09108, 10.26038, 10.43007, 10.59998, 10.7702, 10.94027, 11.11015, 11.28063, 11.45072, 11.62084, 11.79079, 11.96126, 12.13258, 12.30307, 12.47365, 12.64308, 12.81414, 12.98369, 13.15389, 13.32469, 13.49427, 13.66488, 13.8353, 14.00459, 14.17528, 14.346, 14.51503, 14.68449, 14.854, 15.02201, 15.18984, 15.36638, 15.53302, 15.7019, 15.8692, 16.03902, 16.20731, 16.37747, 16.54464, 16.71256, 16.88045, 17.04838, 17.21649, 17.38401, 17.5518, 17.71852, 17.88549, 18.05194, 18.21824, 18.38615, 18.55074, 18.71659, 18.883, 19.04926, 19.21478, 19.37923, 19.54475, 19.70927, 19.87475, 20.0399, 20.20214, 20.36709, 20.53077, 20.69379, 20.85435, 21.01887, 21.18068, 21.3431, 21.50371, 21.6659, 21.82779, 21.98974, 22.15, 22.31198, 22.47064, 22.63182, 22.78999, 22.95058, 23.11054, 23.26745, 23.4241, 23.58469, 23.74521, 23.90188, 24.05784, 24.21452, 24.37171, 24.52704, 24.68218, 24.8359, 24.99103, 25.14763, 25.2993, 25.45367, 25.60873, 25.76254, 25.91534, 26.07072, 26.22102, 26.37055, 26.52309, 26.67165, 26.82309, 26.97359, 27.12229, 27.27053, 27.41958, 27.56691, 27.71686, 27.86294, 28.0101, 28.1581, 28.30505, 28.44739, 28.59395, 28.73862, 28.88294, 29.02769, 29.16797, 29.31115, 29.45215, 29.59618, 29.73789, 29.88035, 30.01922, 30.16022, 30.29821, 30.44038, 30.5761, 30.71154, 30.84592, 30.9825, 31.11152, 31.25309, 31.38929, 31.52258, 31.65494, 31.79034, 31.92079, 32.05578, 32.18832, 32.31833, 32.44526, 32.57908, 32.7034, 32.83536, 32.96328, 33.09283, 33.21908, 33.34493, 33.47003, 33.59101, 33.71882, 33.84377, 33.96453, 34.09085, 34.21141, 34.33102, 34.44985, 34.57335, 34.69049, 34.80889, 34.9282, 35.0416, 35.16209, 35.27414, 35.39553, 35.5093, 35.6183, 35.73292, 35.85022, 35.95582, 36.07387, 36.18346, 36.29449, 36.40315, 36.51531, 36.61855, 36.72639, 36.82914, 36.94118, 37.04585, 37.15234, 37.25531, 37.35564, 37.45971, 37.56297, 37.66224, 37.76324, 37.86612, 37.962, 38.06684, 38.16087, 38.26173, 38.35863, 38.45324, 38.54853, 38.63616, 38.7375, 38.82541, 38.91909, 39.01111, 39.10658, 39.19173, 39.28478, 39.37338, 39.46472, 39.55109, 39.63847, 39.72453, 39.8084, 39.90073, 39.97714, 40.06335, 40.14568, 40.22828, 40.30889, 40.39418, 40.47391, 40.55099, 40.63369, 40.70741, 40.78263, 40.86206, 40.93958, 41.01789, 41.08743, 41.16457, 41.24235, 41.31009, 41.38797, 41.46101, 41.5282, 41.60875, 41.67321, 41.74106, 41.80836, 41.87777, 41.9465, 42.01054, 42.08125, 42.14659, 42.20868, 42.27674, 42.33617, 42.39783, 42.46339, 42.52479, 42.59198, 42.64637, 42.70759, 42.77255, 42.83229, 42.88913, 42.94639, 43.00405, 43.06136, 43.12669, 43.1727, 43.22991, 43.28598, 43.33877, 43.38986, 43.44434, 43.5005, 43.54832, 43.59543, 43.65041, 43.70282, 43.75314, 43.80111, 43.8503, 43.90222, 43.94859, 44.00054, 44.04302, 44.08693, 44.13083, 44.17668, 44.22185, 44.27035, 44.31528, 44.35593, 44.39035, 44.44519, 44.48855, 44.52884, 44.56733, 44.6084, 44.64495, 44.68212, 44.72798, 44.76428, 44.80553, 44.84202, 44.87706, 44.911, 44.95113, 44.98839, 45.02455, 45.05771, 45.09581, 45.13503, 45.16837, 45.19913, 45.72982, 45.26567, 46.29456, 45.32815, 46.37177, 45.38573, 46.43113, 45.45151, 46.5005, 45.51275, 46.56915, 45.57106, 46.63436, 45.6234, 46.69288, 45.68502, 46.75294, 45.73646, 46.81509, 45.78776, 46.86717, 45.83516, 46.91907, 45.87937, 46.97671, 45.95254, 47.04337, 45.99334, 47.46008, 46.03619, 47.89437, 47.1637, 46.87579, 47.23344, 46.92656, 47.28508, 46.98265, 47.34528, 47.04417, 46.29461, 47.42121, 47.11484, 47.47249, 47.16468, 47.52519, 47.22086, 46.45153, 48.38218, 46.49735, 48.43646, 46.54667, 47.69776, 48.51782, 46.60995, 48.57136, 46.65017, 48.61867, 46.69522, 46.7038, 48.6843, 46.73754, 48.73055, 47.94813, 47.61782, 47.99444, 46.83632, 47.67097, 48.05973, 47.72157, 48.08909, 46.92045, 47.76772, 48.14474, 47.80368, 46.98154, 48.20238, 47.86428, 48.23998, 47.8909, 47.07082, 48.29361, 47.95102, 47.11825, 49.20155, 47.98896, 47.1491, 49.25525, 47.1801, 47.19267, 49.31435, 47.21038, 49.3551, 49.36963, 47.26438, 49.39951, 48.54285, 47.29442, 49.43662, 48.58458, 47.33095, 49.494, 48.62708, 47.35855, 49.54331, 48.68034, 48.3043, 48.30865, 48.70515, 48.33351, 48.34123, 48.74259, 48.37224, 47.48147, 48.78339, 48.40807, 47.50649, 48.82853, 48.44402, 47.53651, 48.85572, 49.79967, 48.48679, 47.56321, 49.8519, 48.48892, 47.58864, 
1e-100, 0.109807, 0.2218373, 0.3355839, 0.4512498, 0.5689194, 0.6882621, 0.8095891, 0.9326612, 1.057397, 1.184202, 1.312448, 1.442367, 1.574206, 1.707435, 1.842223, 1.978758, 2.116698, 2.256167, 2.397189, 2.539622, 2.683232, 2.828585, 2.975058, 3.122575, 3.271972, 3.422166, 3.573573, 3.726762, 3.880368, 4.035319, 4.191577, 4.34879, 4.506391, 4.665982, 4.826205, 4.987289, 5.149176, 5.312036, 5.475805, 5.640106, 5.805759, 5.971983, 6.138883, 6.306159, 6.474199, 6.643321, 6.812729, 6.982557, 7.15345, 7.325864, 7.497238, 7.669158, 7.842798, 8.015316, 8.189429, 8.363494, 8.537739, 8.713104, 8.888271, 9.064072, 9.241022, 9.41702, 9.593694, 9.770911, 9.948405, 10.12611, 10.3033, 10.48131, 10.65979, 10.83861, 11.01663, 11.19591, 11.375, 11.55502, 11.73487, 11.91341, 12.0928, 12.27313, 12.45356, 12.63292, 12.8129, 12.99349, 13.17467, 13.35509, 13.53519, 13.71532, 13.89657, 14.07735, 14.25722, 14.4374, 14.61836, 14.79983, 14.98038, 15.16117, 15.34161, 15.52264, 15.7042, 15.88373, 16.06225, 16.24962, 16.42822, 16.60879, 16.78969, 16.97109, 17.15043, 17.33188, 17.51123, 17.69296, 17.87292, 18.05392, 18.23168, 18.41395, 18.59534, 18.77337, 18.9535, 19.13223, 19.31196, 19.49392, 19.67207, 19.85128, 20.03148, 20.21042, 20.3893, 20.56911, 20.74896, 20.9276, 21.10392, 21.28299, 21.4615, 21.63907, 21.81624, 21.99422, 22.16857, 22.34739, 22.52454, 22.70282, 22.88002, 23.05481, 23.23138, 23.40834, 23.58585, 23.75975, 23.93379, 24.11105, 24.28648, 24.4582, 24.63485, 24.80635, 24.97936, 25.15793, 25.32835, 25.4989, 25.67351, 25.84623, 26.01777, 26.18869, 26.36082, 26.53119, 26.7047, 26.87336, 27.04323, 27.21147, 27.38425, 27.55271, 27.7196, 27.88957, 28.05728, 28.22381, 28.39053, 28.55647, 28.72202, 28.88664, 29.05415, 29.21736, 29.38251, 29.54866, 29.71106, 29.87299, 30.03233, 30.19974, 30.36146, 30.5207, 30.68109, 30.84562, 31.00285, 31.16171, 31.31983, 31.48261, 31.63586, 31.79613, 31.95381, 32.1091, 32.26391, 32.4182, 32.572, 32.72671, 32.88418, 33.03738, 33.18288, 33.33337, 33.48217, 33.6332, 33.78633, 33.93384, 34.08099, 34.23045, 34.38048, 34.52395, 34.67318, 34.81868, 34.96012, 35.10612, 35.24899, 35.39236, 35.53825, 35.67929, 35.81767, 35.95782, 36.09715, 36.23493, 36.37826, 36.51266, 36.64723, 36.78474, 36.91984, 37.0555, 37.18917, 37.32379, 37.455, 37.58506, 37.72005, 37.84445, 37.97767, 38.10562, 38.23488, 38.36198, 38.49002, 38.61364, 38.74004, 38.86434, 38.98866, 39.11256, 39.23117, 39.35575, 39.47597, 39.5929, 39.71249, 39.83463, 39.94894, 40.07075, 40.18437, 40.30245, 40.41683, 40.5285, 40.64252, 40.75434, 40.86663, 40.97508, 41.08683, 41.19571, 41.3021, 41.40622, 41.52039, 41.62578, 41.73149, 41.83463, 41.93891, 42.0388, 42.14962, 42.24567, 42.3475, 42.44663, 42.54449, 42.64454, 42.73589, 42.83983, 42.9352, 43.03099, 43.12711, 43.21961, 43.3129, 43.40614, 43.4957, 43.5861, 43.67592, 43.76449, 43.85134, 43.94299, 44.03035, 44.11426, 44.19509, 44.28681, 44.36901, 44.45515, 44.53646, 44.62212, 44.69538, 44.7793, 44.86616, 44.94356, 45.01898, 45.09758, 45.17317, 45.24602, 45.32805, 45.39994, 45.47649, 45.54567, 45.61725, 45.68716, 45.76453, 45.83653, 45.90226, 45.96761, 46.0439, 46.1088, 46.18154, 46.24769, 46.31111, 46.37868, 46.44109, 46.50965, 46.57016, 46.63325, 46.69887, 46.75981, 46.81713, 46.87788, 46.9392, 46.99599, 47.04956, 47.11014, 47.16736, 47.22753, 47.28338, 47.33887, 47.38971, 47.44979, 47.51111, 47.55632, 47.60604, 47.66083, 47.71454, 47.75769, 47.81065, 47.86754, 47.91095, 47.95683, 48.0081, 48.05322, 48.10433, 48.15177, 48.19718, 48.24889, 48.28104, 48.33124, 48.37789, 48.41954, 48.46028, 48.50149, 48.53744, 48.58596, 48.63096, 48.66702, 48.70178, 48.74617, 48.78915, 48.83131, 48.86578, 48.90332, 48.94387, 48.9748, 49.01195, 49.5374, 49.09048, 50.09893, 49.1553, 50.17088, 49.22013, 50.24645, 49.28771, 50.31664, 49.34932, 50.38957, 49.41941, 50.46044, 49.476, 50.51933, 49.53453, 50.57544, 49.59233, 50.64486, 49.64578, 50.69865, 49.70601, 50.76102, 49.75686, 50.81967, 49.82279, 50.88748, 49.8607, 51.39143, 49.91622, 51.92303, 51.02762, 50.93146, 51.08657, 50.98769, 51.1562, 51.04956, 51.21039, 51.1112, 50.19148, 51.29811, 51.19277, 51.35263, 51.24344, 51.4005, 51.29761, 50.36775, 52.44928, 50.41669, 52.50695, 50.45703, 51.59509, 52.59458, 50.52934, 52.64267, 50.57413, 52.69977, 50.61731, 50.63241, 52.77668, 50.67489, 52.82661, 51.85507, 51.74172, 51.90388, 50.77601, 51.80494, 51.96748, 51.84809, 52.00826, 50.86738, 51.90765, 52.06206, 51.94473, 50.93055, 52.12345, 52.0098, 52.17073, 52.04881, 51.02841, 52.22262, 52.09677, 51.0787, 53.32903, 52.14428, 51.11045, 53.37894, 51.14435, 51.1644, 53.45166, 51.18633, 53.48594, 53.49928, 51.22721, 53.53076, 52.48716, 51.27234, 53.58489, 52.52465, 51.29885, 53.64031, 52.58065, 51.34558, 53.69914, 52.62339, 52.48225, 52.49897, 52.65636, 52.52556, 52.52796, 52.69342, 52.5546, 51.45533, 52.74139, 52.60379, 51.50077, 52.7805, 52.64046, 51.52233, 52.82044, 53.95986, 52.68381, 51.55577, 54.01932, 52.68221, 51.57751, 
1e-100, 0.1178169, 0.2374178, 0.3591412, 0.482655, 0.6080581, 0.7355114, 0.8645828, 0.9955552, 1.128377, 1.262769, 1.39904, 1.536927, 1.676376, 1.817726, 1.960491, 2.10468, 2.250738, 2.398108, 2.54678, 2.69721, 2.848791, 3.001706, 3.156397, 3.312069, 3.46905, 3.62736, 3.786895, 3.94747, 4.10928, 4.272229, 4.436178, 4.601421, 4.767504, 4.934951, 5.103, 5.272408, 5.442558, 5.61372, 5.785368, 5.958044, 6.131733, 6.305859, 6.481402, 6.657419, 6.833859, 7.011004, 7.189581, 7.368431, 7.547167, 7.72804, 7.908488, 8.089296, 8.271121, 8.453434, 8.636416, 8.820227, 9.003381, 9.187417, 9.372327, 9.558017, 9.743792, 9.929618, 10.11628, 10.30319, 10.48988, 10.67892, 10.86503, 11.05437, 11.24296, 11.43114, 11.62072, 11.8108, 12.00063, 12.18956, 12.3813, 12.57266, 12.76239, 12.95377, 13.14475, 13.33706, 13.52938, 13.72039, 13.9133, 14.1049, 14.29776, 14.49023, 14.68303, 14.87672, 15.06931, 15.2637, 15.45681, 15.65133, 15.84467, 16.03998, 16.23423, 16.42867, 16.62194, 16.81504, 17.00695, 17.20905, 17.40334, 17.59578, 17.79149, 17.98573, 18.18066, 18.37692, 18.57063, 18.76669, 18.96144, 19.15845, 19.35529, 19.5508, 19.74443, 19.93775, 20.13598, 20.33035, 20.52745, 20.72039, 20.91753, 21.1138, 21.30728, 21.5027, 21.69834, 21.89403, 22.08907, 22.28312, 22.48029, 22.67514, 22.86987, 23.06251, 23.25932, 23.45561, 23.64785, 23.84011, 24.03697, 24.2315, 24.42343, 24.61931, 24.81482, 25.00632, 25.20052, 25.39052, 25.58375, 25.77795, 25.97111, 26.16262, 26.35506, 26.54585, 26.73469, 26.92848, 27.11902, 27.31495, 27.5029, 27.69253, 27.8805, 28.0708, 28.26277, 28.45137, 28.64063, 28.8261, 29.01632, 29.20621, 29.39096, 29.57929, 29.76303, 29.95058, 30.13528, 30.32269, 30.50391, 30.68913, 30.87211, 31.05831, 31.24128, 31.42181, 31.60796, 31.78898, 31.96781, 32.15004, 32.32842, 32.50963, 32.69074, 32.86769, 33.04486, 33.22529, 33.40415, 33.57776, 33.75345, 33.93283, 34.10624, 34.28076, 34.45771, 34.62773, 34.79791, 34.97128, 35.14261, 35.3129, 35.48688, 35.65356, 35.81905, 35.98319, 36.15127, 36.31807, 36.48248, 36.64956, 36.81737, 36.97615, 37.14578, 37.30515, 37.46893, 37.62623, 37.78948, 37.9476, 38.1085, 38.26876, 38.42595, 38.58146, 38.73905, 38.88692, 39.04789, 39.20277, 39.35215, 39.50713, 39.65529, 39.80862, 39.95675, 40.10644, 40.25676, 40.40469, 40.54869, 40.69232, 40.83649, 40.9793, 41.1257, 41.27291, 41.4097, 41.54687, 41.68921, 41.82908, 41.96095, 42.10889, 42.23984, 42.37609, 42.5123, 42.64319, 42.7764, 42.9117, 43.03956, 43.17177, 43.30521, 43.43381, 43.561, 43.68015, 43.81236, 43.9373, 44.05877, 44.1838, 44.30824, 44.42451, 44.55569, 44.66844, 44.78947, 44.90877, 45.02602, 45.14269, 45.25683, 45.37515, 45.48335, 45.59741, 45.7103, 45.82079, 45.92772, 46.04356, 46.15011, 46.26031, 46.36503, 46.47441, 46.57246, 46.68141, 46.78814, 46.88503, 46.98494, 47.08807, 47.18701, 47.28665, 47.38825, 47.48454, 47.57917, 47.68058, 47.76869, 47.8644, 47.9603, 48.05233, 48.1484, 48.23621, 48.32339, 48.42004, 48.50058, 48.58719, 48.68456, 48.7616, 48.86062, 48.93608, 49.01632, 49.10381, 49.18364, 49.26164, 49.34189, 49.42867, 49.50422, 49.57717, 49.66218, 49.73528, 49.81004, 49.88718, 49.9626, 50.04161, 50.10558, 50.18399, 50.25741, 50.32622, 50.39753, 50.46575, 50.52984, 50.60253, 50.67272, 50.73534, 50.79582, 50.8685, 50.92727, 50.99381, 51.05568, 51.12113, 51.18034, 51.2434, 51.30462, 51.36709, 51.42002, 51.48225, 51.5436, 51.59675, 51.65566, 51.71018, 51.7711, 51.81365, 51.87365, 51.92391, 51.98448, 52.03199, 52.08611, 52.13512, 52.18477, 52.24388, 52.29109, 52.33411, 52.38179, 52.43427, 52.47476, 52.52163, 52.5734, 52.61932, 52.65953, 52.70971, 52.74431, 52.79272, 52.83856, 52.88071, 52.9227, 52.96029, 53.00711, 53.04905, 53.08846, 53.12241, 53.16066, 53.66935, 53.23726, 54.2291, 53.31377, 54.30445, 53.37583, 54.38604, 53.46094, 54.46483, 53.52718, 54.53065, 53.59175, 54.59374, 53.66075, 54.67398, 53.72076, 54.73382, 53.78007, 54.7995, 53.83892, 54.87021, 53.89717, 54.92293, 53.94212, 54.99639, 54.02496, 55.06203, 54.06708, 55.65257, 54.1151, 56.26676, 55.21302, 55.31439, 55.27678, 55.37265, 55.34241, 55.45381, 55.40923, 55.50495, 54.41458, 55.4923, 55.59032, 55.54425, 55.65965, 55.60899, 55.71151, 54.5964, 56.83412, 54.65782, 56.89984, 54.7116, 55.80458, 56.98649, 54.77401, 57.04778, 54.82439, 57.10217, 54.86759, 54.88296, 57.16917, 54.93318, 57.2398, 56.08441, 56.18662, 56.13445, 55.03806, 56.26343, 56.19784, 56.2997, 56.23886, 55.13922, 56.35565, 56.30514, 56.4133, 55.21816, 56.368, 56.47044, 56.40684, 56.51839, 55.30713, 56.476, 56.57289, 55.36264, 57.76551, 56.62763, 55.4002, 57.83505, 55.44159, 55.4487, 57.89755, 55.47897, 57.94085, 57.96018, 55.52663, 57.9908, 56.73848, 55.56051, 58.04918, 56.7875, 55.60234, 58.10586, 56.8334, 55.65262, 58.16842, 56.88452, 56.98212, 57.00542, 56.92251, 57.02298, 57.03986, 56.97075, 57.07757, 55.77488, 57.01049, 57.11124, 55.81249, 57.05597, 57.15383, 55.84184, 57.09324, 58.43816, 57.1981, 55.87449, 58.50883, 57.20852, 55.89764, 
1e-100, 0.125154, 0.2524619, 0.3815019, 0.5126104, 0.6456076, 0.7803335, 0.9171446, 1.055573, 1.195857, 1.338118, 1.481843, 1.627343, 1.774646, 1.923344, 2.073826, 2.225853, 2.379257, 2.534521, 2.69118, 2.849138, 3.008646, 3.169792, 3.331791, 3.495327, 3.660509, 3.826475, 3.993982, 4.163128, 4.332906, 4.503932, 4.676534, 4.849917, 5.023871, 5.19991, 5.37664, 5.554387, 5.732969, 5.912953, 6.093592, 6.27469, 6.457397, 6.641033, 6.825487, 7.010109, 7.196119, 7.382893, 7.570115, 7.758712, 7.947707, 8.138344, 8.328166, 8.519027, 8.711597, 8.903763, 9.096908, 9.29074, 9.484958, 9.680181, 9.875719, 10.0722, 10.26945, 10.46666, 10.66385, 10.8627, 11.06249, 11.26156, 11.46044, 11.66093, 11.86241, 12.06326, 12.26454, 12.46669, 12.6701, 12.87336, 13.07682, 13.27996, 13.48498, 13.68988, 13.89446, 14.09954, 14.3056, 14.51156, 14.71973, 14.9255, 15.13233, 15.34005, 15.54792, 15.75651, 15.96375, 16.17313, 16.38124, 16.59227, 16.80154, 17.01116, 17.22044, 17.43144, 17.64096, 17.85265, 18.06169, 18.27884, 18.48975, 18.70028, 18.9132, 19.1235, 19.33582, 19.54813, 19.76189, 19.97429, 20.18827, 20.40213, 20.61333, 20.82756, 21.04122, 21.25318, 21.46807, 21.68012, 21.89761, 22.11262, 22.32568, 22.53961, 22.75281, 22.96776, 23.18183, 23.39504, 23.61136, 23.82421, 24.03965, 24.25321, 24.46778, 24.68309, 24.89768, 25.1104, 25.32317, 25.54085, 25.75272, 25.96784, 26.18182, 26.39631, 26.60977, 26.82049, 27.03425, 27.24971, 27.46133, 27.67131, 27.88188, 28.09837, 28.30829, 28.5222, 28.73021, 28.94519, 29.15756, 29.36727, 29.57862, 29.7871, 30.00051, 30.21128, 30.41715, 30.62661, 30.83496, 31.04463, 31.25068, 31.46221, 31.66517, 31.87007, 32.07834, 32.28642, 32.48683, 32.69376, 32.89971, 33.10276, 33.30662, 33.51037, 33.71374, 33.9127, 34.1199, 34.31803, 34.51768, 34.71615, 34.91549, 35.11708, 35.31257, 35.51053, 35.71238, 35.9043, 36.10175, 36.29792, 36.49468, 36.68498, 36.87913, 37.07405, 37.26262, 37.45286, 37.64596, 37.83195, 38.02675, 38.21774, 38.40252, 38.59249, 38.77207, 38.95131, 39.13783, 39.32268, 39.50936, 39.69099, 39.87183, 40.05307, 40.23558, 40.4121, 40.59043, 40.76453, 40.94364, 41.11935, 41.29329, 41.47705, 41.64104, 41.81213, 41.98518, 42.15938, 42.32394, 42.5009, 42.66326, 42.8297, 42.99416, 43.1662, 43.32492, 43.49158, 43.65314, 43.81243, 43.97511, 44.13736, 44.29243, 44.44867, 44.6107, 44.76758, 44.91785, 45.07879, 45.22618, 45.37721, 45.53013, 45.68023, 45.82777, 45.97793, 46.13095, 46.27257, 46.4125, 46.56115, 46.70831, 46.84051, 46.99837, 47.13081, 47.26944, 47.40812, 47.54779, 47.68149, 47.82002, 47.95261, 48.08652, 48.21889, 48.34944, 48.48205, 48.60509, 48.74129, 48.87174, 48.99402, 49.12595, 49.24355, 49.3631, 49.49588, 49.61347, 49.72911, 49.85514, 49.97134, 50.08821, 50.20479, 50.3267, 50.44072, 50.55166, 50.66955, 50.77559, 50.88966, 51.00133, 51.11052, 51.21366, 51.32537, 51.42804, 51.53486, 51.64237, 51.74361, 51.8507, 51.93981, 52.05315, 52.1514, 52.2522, 52.34659, 52.44885, 52.5359, 52.64299, 52.72847, 52.82886, 52.92288, 53.01161, 53.1012, 53.19202, 53.28211, 53.36993, 53.45652, 53.54269, 53.62616, 53.71417, 53.8006, 53.88401, 53.96667, 54.03807, 54.12961, 54.20958, 54.28993, 54.36372, 54.44615, 54.51674, 54.59868, 54.6755, 54.74393, 54.81933, 54.89414, 54.96116, 55.0361, 55.10717, 55.17567, 55.24201, 55.3107, 55.37454, 55.44861, 55.51488, 55.57865, 55.64828, 55.7017, 55.77905, 55.84056, 55.89607, 55.95698, 56.01722, 56.07786, 56.13519, 56.19715, 56.25681, 56.30797, 56.3635, 56.41852, 56.47961, 56.53406, 56.58738, 56.64159, 56.6955, 56.74931, 56.80266, 56.84696, 56.89765, 56.94236, 56.9902, 57.03951, 57.09454, 57.13843, 57.18651, 57.22255, 57.2799, 57.3305, 57.37338, 57.41202, 57.4527, 57.4989, 57.53302, 57.58482, 57.62568, 58.12015, 57.70011, 58.66085, 57.77649, 58.74631, 57.85927, 58.83021, 57.93576, 58.89848, 58.00772, 58.98163, 58.08116, 59.0568, 58.14329, 59.12105, 58.21224, 59.20118, 58.27454, 59.26838, 58.33053, 59.33503, 58.40561, 59.39873, 58.45669, 59.46641, 58.53022, 59.52854, 58.58163, 60.23857, 58.63688, 60.9435, 59.69088, 60.02762, 59.77306, 60.10217, 59.83804, 60.17264, 59.90062, 60.24199, 58.95533, 59.99951, 60.33039, 60.04822, 60.39245, 60.11297, 60.4575, 59.14656, 61.55011, 59.20572, 61.60613, 59.26753, 60.32439, 61.70748, 59.33429, 61.77182, 59.37626, 61.82948, 59.4365, 59.45088, 61.92195, 59.49662, 61.97195, 60.62212, 60.9672, 60.67013, 59.61626, 61.04076, 60.72747, 61.08231, 60.78132, 59.71378, 61.15518, 60.84432, 61.20409, 59.79695, 60.91967, 61.27294, 60.95861, 61.3154, 59.88851, 61.01548, 61.3825, 59.94471, 62.52874, 61.42618, 59.99561, 62.61193, 60.0446, 60.04011, 62.66912, 60.07173, 62.71721, 62.73379, 60.12445, 62.78771, 61.30859, 60.15975, 62.83048, 61.35808, 60.21322, 62.90282, 61.40786, 60.2492, 62.96067, 61.45181, 61.82385, 61.83898, 61.49045, 61.85838, 61.87192, 61.53822, 61.92159, 60.39737, 61.59381, 61.96606, 60.42763, 61.62409, 62.00972, 60.45291, 61.66753, 63.25089, 62.05019, 60.48529, 63.32545, 62.07658, 60.52771, 
1e-100, 0.1323886, 0.2666167, 0.4031179, 0.5413326, 0.6816238, 0.8238796, 0.9677932, 1.113775, 1.261511, 1.410921, 1.562323, 1.715317, 1.869986, 2.02663, 2.184712, 2.344264, 2.505888, 2.668709, 2.83302, 2.9993, 3.166646, 3.335559, 3.506182, 3.67802, 3.851137, 4.02598, 4.201969, 4.379085, 4.557577, 4.737734, 4.91866, 5.101121, 5.284625, 5.46985, 5.655534, 5.8429, 6.031048, 6.220053, 6.410552, 6.601881, 6.794933, 6.988014, 7.182794, 7.378591, 7.575019, 7.772348, 7.970949, 8.170079, 8.369779, 8.572275, 8.773646, 8.976285, 9.17995, 9.384412, 9.590212, 9.796452, 10.00203, 10.21042, 10.41925, 10.62831, 10.83804, 11.04893, 11.26072, 11.47254, 11.68415, 11.8994, 12.11195, 12.32818, 12.54197, 12.75871, 12.97491, 13.19337, 13.41051, 13.62796, 13.84855, 14.06744, 14.28774, 14.50908, 14.73055, 14.9519, 15.17606, 15.39733, 15.62116, 15.84529, 16.07022, 16.29456, 16.51938, 16.74605, 16.97323, 17.19964, 17.42578, 17.65465, 17.88145, 18.11086, 18.3388, 18.56957, 18.79738, 19.02743, 19.25534, 19.49403, 19.72199, 19.95431, 20.18694, 20.41775, 20.64926, 20.88178, 21.1152, 21.35046, 21.58439, 21.81684, 22.04687, 22.28573, 22.51789, 22.7535, 22.98764, 23.22266, 23.46002, 23.69448, 23.9294, 24.16388, 24.4004, 24.63538, 24.86979, 25.10673, 25.34399, 25.58227, 25.81464, 26.05468, 26.29, 26.52243, 26.76059, 26.99797, 27.23214, 27.46945, 27.70456, 27.9409, 28.17793, 28.41215, 28.64498, 28.88304, 29.11618, 29.35264, 29.58558, 29.82413, 30.05683, 30.29085, 30.52494, 30.75536, 30.98987, 31.23107, 31.46099, 31.69213, 31.92671, 32.15798, 32.38641, 32.62221, 32.84971, 33.08218, 33.31192, 33.54094, 33.76749, 33.99704, 34.22768, 34.453, 34.6825, 34.91165, 35.13753, 35.35953, 35.58847, 35.81511, 36.03762, 36.26142, 36.48502, 36.70478, 36.92786, 37.14873, 37.3717, 37.58969, 37.8096, 38.0281, 38.24464, 38.46158, 38.68087, 38.89775, 39.1127, 39.32917, 39.5441, 39.75417, 39.96436, 40.17584, 40.38334, 40.60056, 40.806, 41.01725, 41.21867, 41.43509, 41.63636, 41.84406, 42.03975, 42.24061, 42.44368, 42.64341, 42.85017, 43.04919, 43.24381, 43.44809, 43.64167, 43.83681, 44.03293, 44.23064, 44.41932, 44.6121, 44.80563, 44.99377, 45.18133, 45.37552, 45.55621, 45.74898, 45.93486, 46.11403, 46.29786, 46.47289, 46.66149, 46.84379, 47.01993, 47.19808, 47.37774, 47.54519, 47.72667, 47.8964, 48.07034, 48.23871, 48.41396, 48.58565, 48.75065, 48.92104, 49.08284, 49.25102, 49.41678, 49.57565, 49.73472, 49.90175, 50.06104, 50.2145, 50.37933, 50.53134, 50.68546, 50.84139, 50.99704, 51.14792, 51.29494, 51.45441, 51.59745, 51.74263, 51.89299, 52.03466, 52.17294, 52.32975, 52.46402, 52.60683, 52.74439, 52.88914, 53.02239, 53.15672, 53.30196, 53.42585, 53.56366, 53.69524, 53.8244, 53.95376, 54.0853, 54.2139, 54.34042, 54.46467, 54.58884, 54.71548, 54.83457, 54.96195, 55.07335, 55.19569, 55.31231, 55.43357, 55.54724, 55.66843, 55.77948, 55.8876, 56.00642, 56.11602, 56.22703, 56.33528, 56.44534, 56.552, 56.65982, 56.76499, 56.87139, 56.96698, 57.07186, 57.16352, 57.27834, 57.38148, 57.46899, 57.5652, 57.6602, 57.75606, 57.84739, 57.94401, 58.0433, 58.13064, 58.21613, 58.30817, 58.40108, 58.49098, 58.5749, 58.66204, 58.74796, 58.83168, 58.91879, 58.99875, 59.07738, 59.16194, 59.23893, 59.31551, 59.40207, 59.4803, 59.54718, 59.62207, 59.69911, 59.7719, 59.84616, 59.92492, 59.99627, 60.05741, 60.134, 60.20647, 60.27629, 60.34395, 60.40816, 60.47513, 60.53555, 60.60685, 60.67511, 60.73162, 60.78964, 60.8513, 60.91401, 60.97831, 61.03663, 61.09597, 61.15271, 61.21373, 61.28131, 61.32938, 61.37983, 61.44159, 61.494, 61.53844, 61.59725, 61.65862, 61.70168, 61.74668, 61.79909, 61.84971, 61.90754, 61.95425, 62.00102, 62.04913, 62.08821, 62.14478, 62.18797, 62.23566, 62.27304, 62.31455, 62.35452, 62.85005, 62.44797, 63.37367, 62.52973, 63.46087, 62.61591, 63.55227, 62.69097, 63.62858, 62.76814, 63.70688, 62.83755, 63.78737, 62.91736, 63.86313, 62.97789, 63.92765, 63.05075, 64.01378, 63.11957, 64.07962, 63.183, 64.1396, 63.24293, 64.21447, 63.32166, 64.28397, 63.36731, 65.07706, 63.43456, 65.89374, 64.45725, 65.02245, 64.52138, 65.10364, 64.61087, 65.16967, 64.669, 65.24197, 63.76425, 64.76355, 65.33925, 64.83699, 65.41046, 64.89408, 65.47918, 63.97674, 66.53925, 64.03319, 66.60066, 64.09199, 65.11291, 66.70123, 64.16413, 66.76767, 64.2144, 66.8221, 64.26602, 64.29485, 66.93105, 64.3461, 66.99062, 65.42993, 66.01306, 65.48541, 64.46164, 66.09849, 65.54018, 66.14104, 65.5893, 64.5655, 66.22848, 65.6758, 66.27765, 64.65225, 65.73085, 66.3488, 65.77928, 66.39297, 64.76128, 65.8457, 66.45914, 64.8099, 67.57614, 66.52277, 64.8719, 67.65636, 64.90834, 64.92387, 67.73217, 64.94552, 67.77449, 67.7845, 64.99892, 67.83369, 66.15063, 65.04615, 67.89787, 66.20266, 65.09507, 67.96726, 66.25694, 65.14708, 68.02568, 66.29922, 66.93392, 66.95114, 66.33628, 66.98206, 67.00867, 66.39663, 67.04204, 65.28033, 66.44353, 67.08898, 65.31696, 66.48847, 67.13547, 65.34488, 66.51841, 68.33996, 67.18923, 65.39983, 68.428, 67.20697, 65.43087, 
1e-100, 0.1394246, 0.2809429, 0.4242783, 0.5699113, 0.7173677, 0.8666989, 1.018229, 1.171376, 1.326665, 1.483863, 1.642683, 1.80351, 1.96615, 2.130281, 2.296434, 2.464229, 2.633505, 2.805035, 2.977849, 3.152186, 3.32842, 3.506213, 3.685127, 3.866006, 4.048456, 4.232006, 4.417409, 4.604651, 4.792509, 4.982121, 5.173447, 5.365707, 5.559165, 5.754794, 5.951429, 6.148897, 6.348067, 6.548736, 6.750404, 6.952749, 7.15755, 7.362765, 7.568844, 7.776578, 7.985832, 8.195352, 8.406392, 8.618678, 8.832591, 9.047779, 9.262273, 9.478703, 9.696839, 9.915155, 10.13397, 10.35568, 10.57666, 10.79922, 11.02314, 11.24794, 11.47376, 11.69996, 11.92634, 12.15477, 12.38478, 12.6153, 12.84488, 13.07662, 13.31007, 13.54308, 13.77652, 14.01105, 14.2493, 14.4855, 14.72278, 14.95957, 15.19893, 15.44017, 15.68075, 15.91949, 16.16312, 16.40583, 16.64875, 16.89308, 17.13537, 17.38147, 17.62704, 17.87381, 18.12026, 18.36718, 18.61707, 18.86648, 19.11675, 19.36718, 19.61534, 19.86784, 20.11909, 20.37046, 20.61704, 20.88022, 21.13203, 21.38549, 21.63874, 21.89108, 22.14745, 22.40402, 22.66082, 22.91644, 23.17236, 23.42759, 23.6846, 23.94402, 24.20028, 24.45974, 24.71937, 24.97642, 25.23506, 25.49221, 25.75274, 26.01097, 26.26913, 26.52795, 26.79015, 27.05166, 27.30707, 27.56948, 27.8313, 28.08927, 28.35036, 28.60968, 28.87207, 29.13022, 29.38784, 29.64686, 29.908, 30.16786, 30.42795, 30.68376, 30.9442, 31.20611, 31.46543, 31.72197, 31.98057, 32.23986, 32.49715, 32.75174, 33.01231, 33.27093, 33.52192, 33.78513, 34.03689, 34.29613, 34.55272, 34.80647, 35.05922, 35.3123, 35.5685, 35.81935, 36.0717, 36.32114, 36.57644, 36.82685, 37.07561, 37.32807, 37.57526, 37.82673, 38.07734, 38.32373, 38.56898, 38.81652, 39.05836, 39.30486, 39.54768, 39.79545, 40.03581, 40.27481, 40.51527, 40.75962, 40.99653, 41.23682, 41.47915, 41.7149, 41.9496, 42.1869, 42.42148, 42.65522, 42.88703, 43.12268, 43.34668, 43.58291, 43.81418, 44.04035, 44.26732, 44.50023, 44.72436, 44.9487, 45.17656, 45.39486, 45.61481, 45.82569, 46.04803, 46.26947, 46.49056, 46.70477, 46.9212, 47.13425, 47.34862, 47.55789, 47.77463, 47.98482, 48.19447, 48.39999, 48.61037, 48.81103, 49.01794, 49.2231, 49.4283, 49.62117, 49.83212, 50.02324, 50.22514, 50.42122, 50.61716, 50.81332, 51.00345, 51.19419, 51.38639, 51.57939, 51.76776, 51.95174, 52.13689, 52.32059, 52.50647, 52.68979, 52.87481, 53.05664, 53.23269, 53.40785, 53.58009, 53.75661, 53.93056, 54.10641, 54.28028, 54.44449, 54.61565, 54.78327, 54.9483, 55.12125, 55.28601, 55.442, 55.60857, 55.76995, 55.92286, 56.09062, 56.24192, 56.39349, 56.55439, 56.70238, 56.85709, 57.00686, 57.15723, 57.30682, 57.45283, 57.59824, 57.74575, 57.88782, 58.03381, 58.17252, 58.30685, 58.44505, 58.58431, 58.7184, 58.85257, 58.99934, 59.12273, 59.25297, 59.38654, 59.51417, 59.64526, 59.77297, 59.89545, 60.01981, 60.14801, 60.26601, 60.39376, 60.51249, 60.63092, 60.74666, 60.85793, 60.98799, 61.0995, 61.20655, 61.32342, 61.43673, 61.53846, 61.65775, 61.75551, 61.86102, 61.98579, 62.08281, 62.17939, 62.28505, 62.38863, 62.49072, 62.58626, 62.68591, 62.78398, 62.88405, 62.98145, 63.07545, 63.16732, 63.25717, 63.36208, 63.45063, 63.53228, 63.62761, 63.71557, 63.79517, 63.89116, 63.97206, 64.05878, 64.13851, 64.2229, 64.30228, 64.38519, 64.46491, 64.54141, 64.62235, 64.70208, 64.7747, 64.86027, 64.92791, 65.00409, 65.07815, 65.14197, 65.22845, 65.29205, 65.36137, 65.42894, 65.4944, 65.56247, 65.6296, 65.6998, 65.76107, 65.82438, 65.89336, 65.95695, 66.02043, 66.07577, 66.1425, 66.20221, 66.2571, 66.31972, 66.37492, 66.43801, 66.48175, 66.54113, 66.59223, 66.65568, 66.70456, 66.76062, 66.8151, 66.85526, 66.92306, 66.97456, 67.01762, 67.0658, 67.117, 67.15965, 67.20589, 67.25713, 67.30755, 67.34709, 67.81242, 67.42691, 68.33513, 67.52911, 68.4261, 67.61199, 68.51849, 67.69952, 68.60437, 67.77255, 68.68527, 67.84662, 68.75684, 67.9212, 68.84854, 67.99997, 68.91992, 68.07236, 68.99531, 68.15831, 69.08081, 68.20908, 69.13324, 68.27073, 69.20844, 68.34976, 69.29089, 68.41202, 70.18695, 68.46755, 71.11053, 69.46882, 70.28771, 69.54563, 70.36049, 69.62187, 70.43936, 69.69449, 70.5137, 68.82983, 69.78183, 70.61626, 69.85441, 70.69375, 69.9324, 70.77131, 69.05286, 71.78695, 69.12153, 71.87279, 69.17172, 70.15744, 71.96268, 69.24909, 72.03059, 69.30931, 72.10349, 69.37383, 69.39046, 72.19888, 69.44794, 72.27044, 70.48172, 71.34974, 70.54713, 69.56216, 71.42069, 70.61115, 71.48397, 70.6654, 69.68473, 71.5653, 70.74623, 71.63381, 69.77978, 70.81693, 71.69977, 70.85648, 71.74105, 69.88271, 70.92656, 71.81829, 69.94389, 72.90571, 71.89559, 70.00963, 72.98778, 70.04708, 70.05552, 73.06239, 70.08064, 73.10744, 73.13265, 70.14665, 73.1736, 71.25086, 70.19895, 73.25908, 71.31162, 70.25036, 73.32331, 71.36473, 70.29003, 73.38099, 71.40687, 72.33789, 72.35323, 71.44944, 72.39287, 72.41509, 71.52075, 72.46292, 70.45758, 71.56696, 72.51165, 70.49025, 71.60393, 72.56023, 70.51979, 71.63697, 73.71861, 72.61656, 70.58477, 73.82568, 72.64316, 70.60145, 
1e-100, 0.1464122, 0.2949169, 0.4457408, 0.5983979, 0.7534711, 0.910495, 1.069314, 1.230566, 1.393562, 1.558482, 1.725686, 1.894499, 2.065361, 2.23839, 2.412918, 2.589406, 2.768131, 2.94815, 3.130176, 3.314229, 3.499445, 3.686994, 3.876312, 4.066985, 4.259351, 4.453944, 4.649661, 4.847022, 5.046462, 5.247492, 5.449339, 5.653817, 5.859489, 6.06671, 6.275083, 6.485967, 6.697632, 6.910537, 7.12565, 7.341795, 7.558875, 7.778325, 7.999495, 8.221441, 8.444114, 8.668704, 8.89505, 9.122189, 9.349923, 9.581428, 9.812427, 10.04521, 10.27892, 10.5143, 10.7512, 10.98869, 11.22623, 11.46716, 11.70926, 11.95175, 12.19425, 12.43993, 12.68599, 12.93345, 13.18014, 13.43108, 13.68028, 13.93185, 14.18368, 14.43825, 14.69254, 14.94899, 15.20452, 15.46166, 15.72159, 15.98082, 16.23977, 16.50211, 16.76459, 17.0284, 17.29296, 17.55636, 17.82161, 18.08987, 18.35742, 18.62514, 18.89256, 19.16326, 19.4344, 19.70607, 19.97762, 20.24843, 20.52101, 20.79682, 21.07148, 21.34561, 21.61958, 21.89703, 22.17152, 22.45793, 22.73342, 23.00909, 23.28889, 23.57031, 23.85101, 24.12965, 24.40878, 24.68889, 24.97394, 25.25575, 25.53755, 25.82084, 26.10146, 26.38659, 26.67054, 26.95149, 27.23747, 27.51955, 27.80617, 28.0904, 28.37594, 28.66088, 28.94523, 29.22946, 29.51607, 29.80378, 30.08795, 30.37383, 30.65452, 30.94217, 31.22664, 31.51146, 31.79358, 32.07905, 32.36381, 32.64864, 32.93338, 33.2181, 33.50297, 33.78492, 34.07124, 34.35517, 34.63667, 34.91677, 35.20084, 35.48329, 35.76203, 36.04479, 36.32044, 36.60677, 36.88639, 37.16351, 37.44056, 37.71765, 37.99755, 38.27489, 38.55275, 38.82721, 39.10706, 39.37706, 39.65324, 39.92719, 40.19626, 40.47402, 40.74433, 41.01955, 41.27995, 41.55735, 41.81984, 42.09055, 42.35454, 42.62346, 42.88257, 43.147, 43.41505, 43.67714, 43.94114, 44.20292, 44.45924, 44.72204, 44.97934, 45.23435, 45.49116, 45.74555, 46.00289, 46.25215, 46.50424, 46.75844, 47.00971, 47.25232, 47.51028, 47.75555, 47.99602, 48.23855, 48.48889, 48.73129, 48.97483, 49.22129, 49.45121, 49.67486, 49.91724, 50.15475, 50.39097, 50.62044, 50.86245, 51.08712, 51.32043, 51.54809, 51.77653, 52.00084, 52.22216, 52.44727, 52.66985, 52.89575, 53.1136, 53.33238, 53.55044, 53.7617, 53.98026, 54.18832, 54.40347, 54.61969, 54.82679, 55.03501, 55.23869, 55.44309, 55.64183, 55.85364, 56.04934, 56.25116, 56.44966, 56.64654, 56.84404, 57.0397, 57.23193, 57.42201, 57.61476, 57.81039, 57.98917, 58.17947, 58.36652, 58.55204, 58.73128, 58.9154, 59.09515, 59.27382, 59.45012, 59.62378, 59.80476, 59.9733, 60.14672, 60.32301, 60.48192, 60.65321, 60.82032, 60.98571, 61.15477, 61.31123, 61.47569, 61.63578, 61.79418, 61.94666, 62.10536, 62.26738, 62.4186, 62.56375, 62.71684, 62.86323, 63.01241, 63.16133, 63.3102, 63.44975, 63.59392, 63.73473, 63.8824, 64.01695, 64.15548, 64.28819, 64.42458, 64.56228, 64.69807, 64.83275, 64.96466, 65.09248, 65.2108, 65.35041, 65.46963, 65.59571, 65.71966, 65.8451, 65.96042, 66.08927, 66.20546, 66.32391, 66.43832, 66.55155, 66.65861, 66.76945, 66.90305, 67.00232, 67.10437, 67.21921, 67.32185, 67.43042, 67.54063, 67.64926, 67.75089, 67.84323, 67.95431, 68.05097, 68.15187, 68.25043, 68.35018, 68.44044, 68.5432, 68.63596, 68.72781, 68.81527, 68.90456, 68.99117, 69.08578, 69.17899, 69.25926, 69.33997, 69.4276, 69.51297, 69.60885, 69.6869, 69.76842, 69.84997, 69.91789, 70.00893, 70.08191, 70.15681, 70.2399, 70.316, 70.38325, 70.45582, 70.52998, 70.60355, 70.66828, 70.73909, 70.80704, 70.88427, 70.95413, 71.02063, 71.08885, 71.14499, 71.22974, 71.28923, 71.34516, 71.4123, 71.46923, 71.52807, 71.58998, 71.65448, 71.71268, 71.76599, 71.82148, 71.88328, 71.94009, 72.00113, 72.05497, 72.10752, 72.15825, 72.20938, 72.26579, 72.31328, 72.36595, 72.41348, 72.45919, 72.50848, 72.56324, 73.01488, 72.65464, 73.51197, 72.74952, 73.6279, 72.84628, 73.70364, 72.92727, 73.7993, 73.00387, 73.88435, 73.10303, 73.96756, 73.16993, 74.05007, 73.25755, 74.14193, 73.33824, 74.22553, 73.41363, 74.30056, 73.48501, 74.37178, 73.54529, 74.44594, 73.62563, 74.51713, 73.6873, 75.53055, 73.75756, 76.56308, 74.71622, 75.7811, 74.80985, 75.86473, 74.88175, 75.9513, 74.95625, 76.02232, 74.13825, 75.06448, 76.14653, 75.14252, 76.22509, 75.20846, 76.30397, 74.37959, 77.29112, 74.44817, 77.37153, 74.50641, 75.44926, 77.47896, 74.58969, 77.55085, 74.65463, 77.62813, 74.73435, 74.74953, 77.74427, 74.80362, 77.80865, 75.80192, 76.93409, 75.87442, 74.93007, 77.0112, 75.93629, 77.0817, 76.01132, 75.07455, 77.17842, 76.08742, 77.23607, 75.15786, 76.15953, 77.31696, 76.21264, 77.37097, 75.28134, 76.26845, 77.45563, 75.35624, 78.49697, 77.53137, 75.41529, 78.58332, 75.46497, 75.46999, 78.66766, 75.50031, 78.71087, 78.73274, 75.56573, 78.79303, 76.63233, 75.63081, 78.87452, 76.68942, 75.69398, 78.95588, 76.74862, 75.73131, 79.01897, 76.80413, 78.02053, 78.04962, 76.86489, 78.09984, 78.11444, 76.91405, 78.15578, 75.91236, 76.96718, 78.21904, 75.9488, 77.01971, 78.27176, 75.99443, 77.06379, 79.40175, 78.35255, 76.05313, 79.49265, 78.37406, 76.09168, 
1e-100, 0.1539907, 0.3101024, 0.4682501, 0.6290866, 0.7918279, 0.9568332, 1.124121, 1.293313, 1.465169, 1.638971, 1.814583, 1.992872, 2.172946, 2.354814, 2.539389, 2.725672, 2.91388, 3.104721, 3.296938, 3.491315, 3.687922, 3.886441, 4.086666, 4.289203, 4.493407, 4.699335, 4.907806, 5.117926, 5.329484, 5.543522, 5.759134, 5.976932, 6.195718, 6.417366, 6.640246, 6.864356, 7.091827, 7.319827, 7.549671, 7.781233, 8.015473, 8.250352, 8.486705, 8.725175, 8.965862, 9.20732, 9.450143, 9.695972, 9.942118, 10.19137, 10.43995, 10.69162, 10.945, 11.19905, 11.45408, 11.71222, 11.97048, 12.23055, 12.49266, 12.75607, 13.02125, 13.28669, 13.55314, 13.82178, 14.09264, 14.3633, 14.63495, 14.91031, 15.18449, 15.4614, 15.73718, 16.01663, 16.29721, 16.57964, 16.86197, 17.14327, 17.42895, 17.71499, 18.0015, 18.28727, 18.5766, 18.8672, 19.15822, 19.44901, 19.74257, 20.03471, 20.32943, 20.62453, 20.92028, 21.21691, 21.51515, 21.81402, 22.11348, 22.41113, 22.71371, 23.0162, 23.31763, 23.61929, 23.91907, 24.23267, 24.53261, 24.84036, 25.14395, 25.44936, 25.75651, 26.06269, 26.37203, 26.6773, 26.98703, 27.29347, 27.60617, 27.91752, 28.22651, 28.53136, 28.83872, 29.15299, 29.4658, 29.77617, 30.08445, 30.397, 30.70876, 31.01985, 31.33124, 31.64221, 31.95558, 32.26344, 32.57416, 32.88795, 33.19718, 33.51009, 33.81803, 34.12956, 34.43892, 34.75175, 35.05883, 35.37412, 35.6826, 35.98665, 36.30057, 36.61066, 36.91735, 37.2274, 37.52972, 37.8361, 38.14433, 38.45179, 38.75584, 39.06348, 39.36403, 39.67008, 39.97167, 40.27274, 40.58148, 40.8782, 41.18034, 41.48048, 41.77698, 42.08235, 42.38019, 42.67784, 42.96819, 43.26653, 43.56099, 43.85545, 44.15193, 44.44741, 44.73272, 45.02435, 45.31367, 45.59939, 45.88771, 46.17482, 46.46242, 46.73952, 47.03177, 47.31243, 47.59542, 47.8741, 48.15269, 48.43383, 48.71037, 48.98306, 49.26555, 49.53341, 49.81097, 50.08046, 50.35507, 50.62014, 50.8886, 51.15753, 51.42339, 51.68816, 51.94999, 52.20901, 52.4747, 52.73397, 52.99397, 53.2628, 53.50348, 53.7434, 54.00383, 54.26048, 54.50507, 54.76026, 55.00412, 55.25094, 55.49019, 55.74118, 55.98659, 56.22228, 56.45986, 56.70414, 56.93404, 57.17322, 57.40902, 57.63873, 57.86945, 58.10224, 58.33093, 58.55528, 58.78703, 59.0097, 59.22646, 59.45453, 59.6683, 59.89268, 60.10921, 60.31861, 60.53647, 60.73754, 60.95648, 61.16199, 61.37574, 61.58451, 61.79068, 61.98707, 62.19589, 62.38948, 62.58365, 62.78756, 62.98636, 63.18153, 63.37237, 63.56297, 63.74937, 63.94541, 64.13654, 64.31585, 64.50732, 64.68825, 64.87085, 65.04857, 65.23378, 65.41434, 65.58613, 65.75587, 65.92628, 66.09629, 66.27375, 66.44655, 66.61144, 66.77067, 66.94175, 67.09889, 67.27403, 67.42683, 67.58412, 67.73403, 67.89995, 68.04887, 68.20527, 68.362, 68.51735, 68.66181, 68.80649, 68.95205, 69.10136, 69.24868, 69.3917, 69.53532, 69.6654, 69.81247, 69.95446, 70.09298, 70.22487, 70.35516, 70.48813, 70.61972, 70.76034, 70.87945, 71.00641, 71.13644, 71.26297, 71.38189, 71.5091, 71.63593, 71.7505, 71.86409, 71.99805, 72.10424, 72.22327, 72.33609, 72.45284, 72.55581, 72.67286, 72.78848, 72.8979, 73.00304, 73.11054, 73.21264, 73.32174, 73.43349, 73.52805, 73.62495, 73.72778, 73.83127, 73.92659, 74.02666, 74.12911, 74.21618, 74.30988, 74.39687, 74.49229, 74.58159, 74.68051, 74.76786, 74.84919, 74.94043, 75.03189, 75.12522, 75.20798, 75.28421, 75.36758, 75.45011, 75.53678, 75.61693, 75.68548, 75.76665, 75.84129, 75.91801, 75.99299, 76.07971, 76.14866, 76.21715, 76.29057, 76.365, 76.44097, 76.50733, 76.57733, 76.6477, 76.7144, 76.78, 76.85165, 76.908, 76.97086, 77.03065, 77.09668, 77.16544, 77.22542, 77.28672, 77.34624, 77.40898, 77.47676, 77.52642, 77.57574, 77.6419, 77.69639, 77.7416, 77.80051, 77.86513, 77.90996, 77.95936, 78.01819, 78.4551, 78.12564, 78.94592, 78.22456, 79.05138, 78.31979, 79.15132, 78.41242, 79.25057, 78.49931, 79.33375, 78.58644, 79.43028, 78.68056, 79.52564, 78.77018, 79.61148, 78.86151, 79.70963, 78.93221, 79.78376, 79.01452, 79.86515, 79.0719, 79.9483, 79.18464, 80.03106, 79.2389, 81.16517, 79.32552, 82.3187, 80.24354, 81.57609, 80.33626, 81.67032, 80.41807, 81.76076, 80.50158, 81.84614, 79.72687, 80.60386, 81.98981, 80.70613, 82.08058, 80.78738, 82.16178, 79.99713, 83.11606, 80.07407, 83.19432, 80.13971, 81.04279, 83.30637, 80.23614, 83.40194, 80.31538, 83.48792, 80.38399, 80.40211, 83.60179, 80.4684, 83.6816, 81.43192, 82.84872, 81.50159, 80.61521, 82.95754, 81.58985, 83.03727, 81.6549, 80.78179, 83.1426, 81.75273, 83.20113, 80.87444, 81.82743, 83.2813, 81.88425, 83.35489, 81.01065, 81.97361, 83.44456, 81.09182, 84.46022, 83.54502, 81.15943, 84.54151, 81.21334, 81.22059, 84.63235, 81.26845, 84.70317, 84.7247, 81.35111, 84.79036, 82.3766, 81.41368, 84.87503, 82.43392, 81.47527, 84.9629, 82.49561, 81.52076, 85.03473, 82.5492, 84.0965, 84.13224, 82.61973, 84.18614, 84.20895, 82.69481, 84.2572, 81.73144, 82.74539, 84.31471, 81.78232, 82.79317, 84.37491, 81.82088, 82.86373, 85.47688, 84.48359, 81.8953, 85.5686, 84.49942, 81.93586, 
1e-100, 0.1616786, 0.3258957, 0.4926105, 0.6614981, 0.833323, 1.007218, 1.183316, 1.36235, 1.543323, 1.726809, 1.912691, 2.100715, 2.291498, 2.484351, 2.67948, 2.877116, 3.077081, 3.278942, 3.483423, 3.690159, 3.898748, 4.110256, 4.323791, 4.539102, 4.757068, 4.977492, 5.199498, 5.423642, 5.650662, 5.879345, 6.11031, 6.3437, 6.578691, 6.815625, 7.055608, 7.297502, 7.540301, 7.785379, 8.033935, 8.283082, 8.534542, 8.78826, 9.044502, 9.301391, 9.560485, 9.822922, 10.08634, 10.35089, 10.61857, 10.88756, 11.15806, 11.43007, 11.70432, 11.98098, 12.25872, 12.53782, 12.81792, 13.10259, 13.38715, 13.67347, 13.96043, 14.25119, 14.54114, 14.83415, 15.12784, 15.42435, 15.72027, 16.0182, 16.31925, 16.6216, 16.9241, 17.22754, 17.5323, 17.83986, 18.14975, 18.45854, 18.76707, 19.08063, 19.39482, 19.70921, 20.02338, 20.33976, 20.65784, 20.97626, 21.29614, 21.61529, 21.93753, 22.25804, 22.58376, 22.90557, 23.2322, 23.55756, 23.887, 24.21231, 24.54312, 24.8698, 25.19974, 25.53022, 25.85806, 26.19922, 26.52698, 26.86035, 27.19221, 27.52746, 27.8606, 28.19517, 28.53211, 28.86747, 29.20371, 29.53936, 29.87715, 30.21194, 30.54714, 30.88526, 31.22226, 31.56022, 31.89738, 32.23631, 32.57796, 32.91315, 33.24818, 33.58676, 33.92663, 34.26264, 34.59997, 34.93859, 35.27481, 35.61468, 35.95098, 36.28615, 36.62587, 36.96185, 37.29091, 37.63202, 37.9678, 38.30215, 38.63812, 38.97109, 39.30737, 39.64253, 39.96949, 40.30382, 40.635, 40.96607, 41.29474, 41.62032, 41.95587, 42.28175, 42.60556, 42.94171, 43.26134, 43.59493, 43.91629, 44.24229, 44.56001, 44.88759, 45.21187, 45.52862, 45.85093, 46.16995, 46.48697, 46.805, 47.12398, 47.43267, 47.74584, 48.06262, 48.37308, 48.68353, 48.99463, 49.30398, 49.61007, 49.91467, 50.22279, 50.52574, 50.8321, 51.13227, 51.43161, 51.73232, 52.02729, 52.33009, 52.61729, 52.91674, 53.20938, 53.50263, 53.79069, 54.08378, 54.37312, 54.65701, 54.93764, 55.23286, 55.51094, 55.79087, 56.06813, 56.3541, 56.62694, 56.90834, 57.18532, 57.46721, 57.72076, 57.97195, 58.24682, 58.51633, 58.78841, 59.05577, 59.31389, 59.57511, 59.8354, 60.09478, 60.35178, 60.60279, 60.85855, 61.11291, 61.36227, 61.61664, 61.86152, 62.10621, 62.34914, 62.59363, 62.83213, 63.07184, 63.32258, 63.54965, 63.78336, 64.02034, 64.2506, 64.47129, 64.71714, 64.93605, 65.16265, 65.3865, 65.61228, 65.82652, 66.05163, 66.268, 66.48447, 66.70453, 66.91526, 67.12681, 67.33045, 67.54794, 67.75149, 67.95946, 68.16902, 68.36847, 68.56375, 68.76826, 68.9615, 69.15518, 69.36265, 69.55393, 69.74207, 69.93321, 70.1172, 70.30285, 70.49716, 70.68198, 70.85855, 71.0375, 71.2256, 71.40106, 71.5786, 71.75435, 71.93231, 72.09865, 72.28091, 72.43763, 72.60841, 72.77766, 72.94092, 73.10881, 73.26326, 73.42794, 73.59581, 73.75486, 73.90605, 74.06284, 74.21041, 74.37505, 74.52972, 74.68144, 74.82474, 74.97195, 75.11534, 75.25594, 75.40919, 75.54802, 75.68653, 75.82882, 75.96754, 76.09841, 76.24086, 76.37186, 76.51206, 76.63698, 76.76665, 76.89154, 77.02703, 77.14706, 77.28734, 77.39955, 77.52894, 77.64717, 77.77531, 77.89796, 78.01914, 78.12945, 78.247, 78.36684, 78.47826, 78.58678, 78.70566, 78.81503, 78.91983, 79.03509, 79.13752, 79.24862, 79.34737, 79.45379, 79.55328, 79.66703, 79.76485, 79.85916, 79.95317, 80.06862, 80.16282, 80.27006, 80.36, 80.44721, 80.53719, 80.62599, 80.72759, 80.80789, 80.90677, 80.99088, 81.07785, 81.16227, 81.24966, 81.33873, 81.41816, 81.49815, 81.58788, 81.67065, 81.7641, 81.83322, 81.91193, 81.98803, 82.05849, 82.15408, 82.22028, 82.29566, 82.36767, 82.43376, 82.50911, 82.58284, 82.65316, 82.72257, 82.79353, 82.86193, 82.93764, 83.00255, 83.0648, 83.13211, 83.1986, 83.26365, 83.32037, 83.38711, 83.45154, 83.50071, 83.56547, 83.61824, 83.68428, 83.75382, 83.80928, 84.23234, 83.91513, 84.72147, 84.03857, 84.82843, 84.14591, 84.93847, 84.24689, 85.04327, 84.35053, 85.15759, 84.45396, 85.2542, 84.55722, 85.37027, 84.65669, 85.46197, 84.74601, 85.56596, 84.83767, 85.65472, 84.92532, 85.74422, 84.99023, 85.82905, 85.10705, 85.92188, 85.18888, 87.18995, 85.2818, 88.46813, 86.17859, 87.80653, 86.26769, 87.89709, 86.36303, 87.9975, 86.45778, 88.10208, 85.74937, 86.60082, 88.25442, 86.68709, 88.3669, 86.78579, 88.45256, 86.05483, 89.37486, 86.13492, 89.46431, 86.21502, 87.08372, 89.60964, 86.33452, 89.70845, 86.42024, 89.81064, 86.50137, 86.53186, 89.93953, 86.59853, 90.02201, 87.52962, 89.25611, 87.59584, 86.77008, 89.37871, 87.72542, 89.47944, 87.79543, 86.95291, 89.5803, 87.88963, 89.65105, 87.07556, 87.9901, 89.76532, 88.04988, 89.8294, 87.23553, 88.16463, 89.9695, 87.34185, 90.92596, 90.06222, 87.42233, 91.03478, 87.48648, 87.49826, 91.12812, 87.53801, 91.19776, 91.24012, 87.64282, 91.31683, 88.62789, 87.71674, 91.42346, 88.70922, 87.78816, 91.51309, 88.76979, 87.85997, 91.59882, 88.84771, 90.73892, 90.77756, 88.92957, 90.82457, 90.84565, 89.0093, 90.90815, 88.11617, 89.07599, 90.9969, 88.17419, 89.13908, 91.06465, 88.23269, 89.2187, 92.12083, 91.17703, 88.3185, 92.24215, 91.21312, 88.36688, 
1e-100, 0.1708842, 0.3440949, 0.5199293, 0.6987872, 0.8798957, 1.063956, 1.250562, 1.439616, 1.631884, 1.826417, 2.023295, 2.223651, 2.426003, 2.630848, 2.838826, 3.049004, 3.262002, 3.477794, 3.695691, 3.916315, 4.139976, 4.365456, 4.593738, 4.825104, 5.058383, 5.294145, 5.532963, 5.773722, 6.017157, 6.263678, 6.512232, 6.762761, 7.015461, 7.272076, 7.53017, 7.790555, 8.053703, 8.318893, 8.586521, 8.856432, 9.128623, 9.40336, 9.679707, 9.958856, 10.23974, 10.52334, 10.8096, 11.09696, 11.38663, 11.67891, 11.97224, 12.26912, 12.56671, 12.86603, 13.16868, 13.47282, 13.77909, 14.08614, 14.3953, 14.70866, 15.02282, 15.33702, 15.6539, 15.97231, 16.29277, 16.61425, 16.93677, 17.26515, 17.59008, 17.91958, 18.2482, 18.58077, 18.91402, 19.24716, 19.58316, 19.92022, 20.26002, 20.59929, 20.93973, 21.28126, 21.62606, 21.97121, 22.31815, 22.66215, 23.01181, 23.35924, 23.70976, 24.0606, 24.41338, 24.76439, 25.12212, 25.47588, 25.82977, 26.18472, 26.54224, 26.89761, 27.25665, 27.6149, 27.96942, 28.34035, 28.69888, 29.05735, 29.4182, 29.78111, 30.14128, 30.50577, 30.86755, 31.23393, 31.59707, 31.96099, 32.32112, 32.68547, 33.05145, 33.4164, 33.78035, 34.14469, 34.51013, 34.87411, 35.24213, 35.60253, 35.96856, 36.33287, 36.69554, 37.05941, 37.42405, 37.7871, 38.15162, 38.51229, 38.87755, 39.24013, 39.59817, 39.9629, 40.31899, 40.6821, 41.04381, 41.40326, 41.76431, 42.12234, 42.47728, 42.83471, 43.19355, 43.54787, 43.89689, 44.25377, 44.61261, 44.9605, 45.31533, 45.66383, 46.00514, 46.36855, 46.71628, 47.06285, 47.40853, 47.75904, 48.10198, 48.44297, 48.79221, 49.12706, 49.47231, 49.80564, 50.14989, 50.48616, 50.81872, 51.15816, 51.48982, 51.82666, 52.155, 52.4858, 52.80887, 53.14788, 53.47076, 53.80044, 54.11604, 54.44253, 54.76513, 55.08227, 55.39804, 55.72276, 56.03356, 56.3529, 56.66636, 56.98342, 57.28764, 57.60043, 57.90635, 58.21572, 58.51961, 58.82219, 59.12252, 59.42475, 59.72435, 60.01628, 60.31567, 60.61286, 60.90548, 61.19905, 61.49915, 61.78892, 62.05491, 62.32502, 62.62086, 62.90096, 63.19197, 63.46828, 63.74233, 64.01771, 64.29685, 64.5668, 64.8393, 65.10822, 65.37755, 65.639, 65.90815, 66.17077, 66.43133, 66.69371, 66.95893, 67.20731, 67.46503, 67.71828, 67.97176, 68.2143, 68.4661, 68.71608, 68.96326, 69.20514, 69.4457, 69.6833, 69.91897, 70.1647, 70.39858, 70.62617, 70.86563, 71.09626, 71.31994, 71.55276, 71.77443, 72.00263, 72.21525, 72.44612, 72.66248, 72.88253, 73.09512, 73.31128, 73.52832, 73.74187, 73.94632, 74.15104, 74.37298, 74.57085, 74.77352, 74.98276, 75.18015, 75.37214, 75.57701, 75.7625, 75.96218, 76.15517, 76.35309, 76.54181, 76.72709, 76.909, 77.10019, 77.28804, 77.47528, 77.64256, 77.8204, 78.00934, 78.18079, 78.3595, 78.53737, 78.71436, 78.87639, 79.04802, 79.21756, 79.38895, 79.5505, 79.71494, 79.8798, 80.03152, 80.19888, 80.36339, 80.51497, 80.67252, 80.82385, 80.97241, 81.13482, 81.28543, 81.42881, 81.57851, 81.72637, 81.86376, 82.01234, 82.16107, 82.29973, 82.43505, 82.5768, 82.69771, 82.86605, 82.99398, 83.12868, 83.26415, 83.39219, 83.5204, 83.65368, 83.78169, 83.91325, 84.03712, 84.16028, 84.2968, 84.41373, 84.52938, 84.65291, 84.77126, 84.87924, 85.00361, 85.12104, 85.24045, 85.34635, 85.46212, 85.57105, 85.68901, 85.79845, 85.91311, 86.01556, 86.11519, 86.23343, 86.33604, 86.4413, 86.54267, 86.64332, 86.73667, 86.83801, 86.94577, 87.03683, 87.12659, 87.22855, 87.31803, 87.42567, 87.52727, 87.62047, 87.7046, 87.78616, 87.89134, 87.97053, 88.06145, 88.15073, 88.23604, 88.31913, 88.40337, 88.48565, 88.57241, 88.6446, 88.72905, 88.81468, 88.90109, 88.98427, 89.06257, 89.14047, 89.21167, 89.30089, 89.37821, 89.4397, 89.51746, 89.58694, 89.65426, 89.72946, 89.80104, 89.87586, 89.9412, 90.00622, 90.08471, 90.15212, 90.56327, 90.28787, 91.04153, 90.41411, 91.17649, 90.54026, 91.30966, 90.65466, 91.42293, 90.77353, 91.5448, 90.90033, 91.68081, 91.0194, 91.80312, 91.12002, 91.92966, 91.25292, 92.02882, 91.36042, 92.13991, 91.45889, 92.24951, 91.54733, 92.3559, 91.6922, 92.45742, 91.77297, 93.88438, 91.88582, 95.34125, 92.75049, 94.72728, 92.86547, 94.83892, 92.98504, 94.96424, 93.09196, 95.09108, 92.44615, 93.2713, 95.27916, 93.37578, 95.40038, 93.48026, 95.51459, 92.81088, 96.40826, 92.90161, 96.51117, 93.0037, 93.82733, 96.69573, 93.15995, 96.81181, 93.24999, 96.91045, 93.35886, 93.38414, 97.0771, 93.47616, 97.19017, 94.37166, 96.47779, 94.46995, 93.69108, 96.63878, 94.60525, 96.72872, 94.69662, 93.9073, 96.88038, 94.81424, 96.96384, 94.05212, 94.92656, 97.10473, 95.0146, 97.19591, 94.26672, 95.12761, 97.35348, 94.38715, 98.26139, 97.45519, 94.48516, 98.40051, 94.55769, 94.5843, 98.53305, 94.64827, 98.61608, 98.65035, 94.77756, 98.74344, 95.71348, 94.85845, 98.87242, 95.81211, 94.94339, 98.97544, 95.89664, 95.04381, 99.09957, 95.99417, 98.29382, 98.33285, 96.1027, 98.40741, 98.43815, 96.19087, 98.50487, 95.34995, 96.26634, 98.60798, 95.43209, 96.35186, 98.71082, 95.52484, 96.45114, 99.72776, 98.83624, 95.61381, 99.86199, 98.88847, 95.67623, 
1e-100, 0.180341, 0.3639509, 0.5503717, 0.7396347, 0.9323773, 1.127622, 1.325893, 1.527447, 1.731448, 1.938766, 2.149036, 2.361887, 2.578436, 2.797551, 3.019175, 3.244713, 3.472348, 3.702916, 3.936912, 4.173693, 4.412806, 4.656188, 4.901476, 5.149424, 5.400947, 5.654971, 5.911753, 6.171538, 6.43396, 6.699308, 6.967176, 7.239033, 7.512007, 7.788194, 8.067794, 8.349977, 8.633934, 8.921678, 9.21169, 9.503786, 9.79793, 10.09589, 10.39607, 10.69921, 11.0034, 11.31184, 11.62111, 11.93336, 12.24877, 12.56704, 12.88435, 13.20533, 13.53032, 13.8571, 14.18471, 14.51416, 14.84574, 15.18192, 15.51892, 15.8572, 16.19765, 16.54015, 16.88306, 17.23056, 17.57806, 17.9298, 18.27987, 18.63315, 18.9894, 19.3467, 19.70385, 20.06462, 20.42489, 20.78791, 21.15497, 21.51894, 21.88452, 22.2561, 22.62635, 22.99666, 23.36646, 23.73988, 24.11702, 24.49205, 24.86913, 25.24575, 25.62487, 26.00602, 26.38577, 26.76655, 27.14792, 27.53107, 27.91507, 28.30133, 28.68679, 29.07097, 29.45931, 29.8461, 30.22573, 30.62793, 31.01291, 31.3998, 31.79029, 32.17914, 32.56698, 32.9577, 33.34949, 33.7421, 34.13145, 34.5204, 34.91461, 35.30391, 35.69624, 36.08406, 36.47892, 36.86899, 37.26199, 37.65114, 38.03583, 38.43166, 38.82213, 39.20907, 39.5998, 39.98937, 40.38112, 40.76367, 41.15623, 41.54534, 41.92903, 42.31832, 42.70338, 43.09109, 43.47131, 43.85708, 44.23865, 44.62454, 45.00765, 45.3871, 45.76794, 46.14733, 46.52968, 46.90848, 47.2856, 47.66238, 48.03708, 48.41115, 48.78513, 49.15927, 49.52814, 49.9079, 50.27687, 50.64133, 51.01043, 51.3758, 51.74432, 52.10887, 52.4728, 52.83491, 53.19701, 53.55043, 53.91085, 54.26998, 54.62346, 54.9797, 55.33192, 55.69718, 56.03527, 56.38773, 56.73943, 57.08418, 57.42894, 57.77574, 58.11882, 58.45378, 58.79967, 59.1338, 59.47226, 59.80836, 60.14018, 60.47267, 60.80965, 61.13713, 61.46652, 61.788, 62.1218, 62.44114, 62.76134, 63.08511, 63.40584, 63.72268, 64.03792, 64.35399, 64.67171, 64.97674, 65.2876, 65.60243, 65.90991, 66.22207, 66.5111, 66.78859, 67.10107, 67.40576, 67.70365, 68.00029, 68.28869, 68.58395, 68.87833, 69.15933, 69.4556, 69.73788, 70.0236, 70.30527, 70.58441, 70.86759, 71.13846, 71.41773, 71.69154, 71.96531, 72.23844, 72.50867, 72.77079, 73.04049, 73.30309, 73.57045, 73.82113, 74.08713, 74.34779, 74.59805, 74.85621, 75.11332, 75.36073, 75.60333, 75.8573, 76.11163, 76.35579, 76.59495, 76.83791, 77.06713, 77.30961, 77.5552, 77.78896, 78.02724, 78.25504, 78.47836, 78.71433, 78.93766, 79.16743, 79.38629, 79.62694, 79.83855, 80.05635, 80.26959, 80.4906, 80.71367, 80.92237, 81.12702, 81.33578, 81.55509, 81.7614, 81.9714, 82.1737, 82.37677, 82.57806, 82.78683, 82.98819, 83.17588, 83.37932, 83.56563, 83.76746, 83.95188, 84.14326, 84.34367, 84.52454, 84.71186, 84.89961, 85.07143, 85.26781, 85.44433, 85.61986, 85.79843, 85.97808, 86.15384, 86.325, 86.50791, 86.66867, 86.84369, 87.0082, 87.17446, 87.34385, 87.50474, 87.67621, 87.83233, 87.99254, 88.15768, 88.31169, 88.47247, 88.61968, 88.77048, 88.94378, 89.09759, 89.24643, 89.4026, 89.54501, 89.69398, 89.83643, 89.98537, 90.13249, 90.26974, 90.41554, 90.55787, 90.69559, 90.83624, 90.97015, 91.11172, 91.23838, 91.37593, 91.50374, 91.63751, 91.77513, 91.90218, 92.03613, 92.14719, 92.29049, 92.41979, 92.54512, 92.66595, 92.78833, 92.90912, 93.02253, 93.14751, 93.2704, 93.38272, 93.50186, 93.61048, 93.72277, 93.84263, 93.95134, 94.07733, 94.17519, 94.29054, 94.40538, 94.51456, 94.62644, 94.72341, 94.83006, 94.93462, 95.04549, 95.13944, 95.23862, 95.34378, 95.43715, 95.53636, 95.6328, 95.74679, 95.83882, 95.92188, 96.02346, 96.11954, 96.21974, 96.30899, 96.40286, 96.49197, 96.58449, 96.66469, 96.75543, 96.83269, 96.92473, 97.00316, 97.10204, 97.18885, 97.27819, 97.3524, 97.43943, 97.83356, 97.61608, 98.33402, 97.75446, 98.49827, 97.90523, 98.65083, 98.06523, 98.79959, 98.2176, 98.95842, 98.36808, 99.12265, 98.51013, 99.25841, 98.64917, 99.40116, 98.78844, 99.54245, 98.92195, 99.67981, 99.04816, 99.81123, 99.15608, 99.94946, 99.33072, 100.0818, 99.4432, 101.7219, 99.56877, 103.3738, 100.4316, 102.8091, 100.5722, 102.9649, 100.7122, 103.1091, 100.8627, 103.2854, 100.2874, 101.0519, 103.5407, 101.201, 103.6743, 101.3427, 103.8262, 100.7259, 104.6582, 100.8676, 104.8115, 100.9868, 101.7871, 105.0612, 101.1731, 105.2105, 101.3061, 105.3085, 101.423, 101.4604, 105.5351, 101.5689, 105.669, 102.4531, 104.9899, 102.5917, 101.8646, 105.2392, 102.7434, 105.3161, 102.8644, 102.1318, 105.5278, 103.0211, 105.6269, 102.3311, 103.156, 105.8283, 103.2773, 105.9095, 102.6016, 103.4342, 106.1232, 102.7442, 106.9857, 106.2734, 102.8769, 107.1489, 102.9686, 102.9996, 107.3332, 103.0918, 107.4325, 107.5188, 103.2492, 107.5915, 104.1496, 103.3771, 107.7817, 104.2814, 103.4719, 107.9264, 104.3876, 103.5947, 108.0814, 104.5332, 107.3038, 107.39, 104.6413, 107.4332, 107.5202, 104.7704, 107.5774, 104.0053, 104.8804, 107.6901, 104.1074, 104.9958, 107.8243, 104.2191, 105.1126, 108.8128, 108.0358, 104.3414, 108.9836, 108.1146, 104.4354, 
1e-100, 0.1921179, 0.3870421, 0.5854711, 0.7872997, 0.9920075, 1.200595, 1.412106, 1.626938, 1.845581, 2.067053, 2.291914, 2.520431, 2.751771, 2.986759, 3.224949, 3.466087, 3.711164, 3.959717, 4.210721, 4.465569, 4.723607, 4.984852, 5.249098, 5.517007, 5.78769, 6.061641, 6.339181, 6.619748, 6.90286, 7.190265, 7.480138, 7.772979, 8.068706, 8.368062, 8.669819, 8.975489, 9.283923, 9.594725, 9.908082, 10.22474, 10.54517, 10.86701, 11.19142, 11.52033, 11.85081, 12.18396, 12.52101, 12.85907, 13.19977, 13.54403, 13.88973, 14.23957, 14.58914, 14.94197, 15.29856, 15.65779, 16.01635, 16.37863, 16.74313, 17.11006, 17.47955, 17.84899, 18.22417, 18.59854, 18.9744, 19.35416, 19.73253, 20.11609, 20.49997, 20.88567, 21.27153, 21.6624, 22.05309, 22.44272, 22.83738, 23.23394, 23.63022, 24.02484, 24.4246, 24.82471, 25.22468, 25.62962, 26.03377, 26.43556, 26.84249, 27.24815, 27.655, 28.06533, 28.47206, 28.88434, 29.29481, 29.70748, 30.12037, 30.53272, 30.94446, 31.35758, 31.77188, 32.18221, 32.59422, 33.02807, 33.43847, 33.85288, 34.26894, 34.68355, 35.10363, 35.52247, 35.93941, 36.35316, 36.76868, 37.18862, 37.60672, 38.02389, 38.44016, 38.85505, 39.27495, 39.68893, 40.10679, 40.52206, 40.93594, 41.35199, 41.76649, 42.18345, 42.60239, 43.00874, 43.42528, 43.83764, 44.25, 44.6655, 45.07445, 45.48072, 45.89346, 46.30322, 46.70361, 47.11936, 47.52152, 47.93069, 48.33466, 48.73813, 49.14384, 49.54909, 49.95074, 50.35162, 50.75012, 51.14864, 51.5453, 51.9412, 52.33998, 52.72857, 53.12291, 53.52054, 53.91285, 54.30609, 54.69371, 55.07828, 55.46316, 55.85096, 56.23366, 56.61812, 57.00221, 57.38141, 57.75831, 58.13582, 58.51386, 58.89026, 59.26386, 59.64654, 60.00982, 60.37839, 60.75206, 61.11516, 61.47798, 61.84361, 62.20773, 62.56557, 62.93141, 63.28569, 63.64303, 63.9993, 64.35674, 64.71375, 65.06009, 65.40953, 65.75596, 66.10633, 66.45391, 66.79273, 67.13248, 67.48033, 67.81399, 68.15244, 68.4894, 68.82858, 69.16016, 69.49394, 69.82118, 70.15432, 70.47797, 70.82005, 71.12093, 71.42543, 71.74905, 72.08179, 72.39666, 72.71257, 73.02801, 73.3443, 73.65679, 73.96096, 74.27995, 74.5782, 74.87664, 75.19352, 75.48701, 75.79413, 76.09223, 76.39357, 76.68838, 76.97429, 77.26992, 77.55777, 77.85329, 78.1377, 78.42765, 78.71996, 78.9881, 79.27147, 79.55386, 79.82965, 80.11819, 80.38456, 80.65791, 80.92168, 81.20447, 81.4629, 81.75067, 82.00742, 82.26995, 82.53209, 82.7934, 83.05628, 83.31539, 83.57202, 83.83005, 84.08944, 84.33396, 84.58969, 84.8277, 85.08344, 85.35102, 85.58864, 85.83312, 86.07374, 86.30658, 86.55525, 86.78629, 87.02288, 87.26077, 87.49716, 87.73548, 87.9646, 88.18604, 88.42247, 88.65073, 88.87399, 89.10931, 89.31146, 89.54855, 89.7683, 89.98968, 90.20164, 90.42455, 90.62955, 90.85448, 91.05643, 91.27529, 91.48548, 91.68536, 91.89992, 92.09834, 92.30967, 92.50948, 92.71737, 92.91392, 93.10794, 93.30962, 93.503, 93.7102, 93.88558, 94.08895, 94.27583, 94.46537, 94.6556, 94.84761, 95.03281, 95.2127, 95.40174, 95.57694, 95.75877, 95.97629, 96.13913, 96.31934, 96.49004, 96.66888, 96.84472, 97.01016, 97.19278, 97.35249, 97.52604, 97.70505, 97.86615, 98.03219, 98.19404, 98.3503, 98.51605, 98.67998, 98.8418, 99.00973, 99.16192, 99.33149, 99.48107, 99.63969, 99.7912, 99.94144, 100.094, 100.2315, 100.3983, 100.5386, 100.701, 100.8342, 100.9862, 101.1307, 101.2639, 101.412, 101.5437, 101.6838, 101.8352, 101.9788, 102.1207, 102.2575, 102.386, 102.5159, 102.6377, 102.7902, 102.9057, 103.0537, 103.1654, 103.2982, 103.4156, 103.5521, 103.6824, 103.8084, 103.9208, 104.0599, 104.1746, 104.31, 104.42, 104.5362, 104.6598, 104.7607, 104.8908, 105.0036, 105.1071, 105.2284, 105.3373, 105.4496, 105.5596, 105.6612, 105.7907, 105.8904, 105.9822, 106.1133, 106.2073, 106.3087, 106.7219, 106.5101, 107.2316, 106.7091, 107.4305, 106.9126, 107.6224, 107.0928, 107.8201, 107.3024, 108.0162, 107.4942, 108.213, 107.6846, 108.3958, 107.8714, 108.5789, 108.0272, 108.7665, 108.2055, 108.9359, 108.3734, 109.1083, 108.5392, 109.2811, 108.7172, 109.4593, 108.8676, 111.3555, 109.0443, 113.3017, 109.8867, 112.8447, 110.0771, 113.041, 110.2689, 113.2359, 110.4505, 113.4476, 109.9658, 110.7125, 113.7295, 110.8883, 113.9061, 111.0574, 114.1055, 110.5303, 114.9445, 110.7157, 115.1471, 110.8828, 111.6471, 115.4246, 111.1071, 115.6034, 111.2474, 115.7883, 111.4251, 111.4815, 116.0333, 111.6325, 116.1959, 112.4956, 115.6603, 112.6683, 112.0079, 115.8944, 112.8796, 116.0725, 113.0171, 112.3536, 116.2796, 113.2184, 116.4368, 112.6256, 113.4197, 116.6602, 113.5675, 116.7999, 112.9718, 113.7798, 117.0381, 113.1517, 117.9372, 117.242, 113.3284, 118.1405, 113.4499, 113.5111, 118.3552, 113.6024, 118.4943, 118.5701, 113.8127, 118.6832, 114.7076, 113.9655, 118.8757, 114.8511, 114.1217, 119.0751, 115.0124, 114.2724, 119.2699, 115.1916, 118.6149, 118.6789, 115.3395, 118.7835, 118.8386, 115.4725, 118.9546, 114.8162, 115.6375, 119.1229, 114.921, 115.7886, 119.2957, 115.0853, 115.9149, 120.255, 119.4972, 115.244, 120.4469, 119.6244, 115.3993, 
1e-100, 0.2041915, 0.4126272, 0.6243259, 0.8398766, 1.059386, 1.282031, 1.508791, 1.739172, 1.972839, 2.210911, 2.452223, 2.697014, 2.946403, 3.19873, 3.454751, 3.715086, 3.978408, 4.245307, 4.516559, 4.790952, 5.068863, 5.351056, 5.636434, 5.924743, 6.218082, 6.513908, 6.812863, 7.116322, 7.422417, 7.73203, 8.046028, 8.363127, 8.682775, 9.005766, 9.333588, 9.662878, 9.99529, 10.33241, 10.67197, 11.01278, 11.35975, 11.7082, 12.06026, 12.41427, 12.77123, 13.13263, 13.49604, 13.86034, 14.23014, 14.60242, 14.97432, 15.35099, 15.73062, 16.11213, 16.49638, 16.88196, 17.26983, 17.66176, 18.05493, 18.45019, 18.84782, 19.24749, 19.64858, 20.05232, 20.45876, 20.86702, 21.27559, 21.68551, 22.09979, 22.51487, 22.92962, 23.34771, 23.76743, 24.18815, 24.61111, 25.03381, 25.45884, 25.88651, 26.31412, 26.73941, 27.17007, 27.60108, 28.03611, 28.46766, 28.90374, 29.33565, 29.76963, 30.20799, 30.64365, 31.08154, 31.51935, 31.96088, 32.40064, 32.83946, 33.27861, 33.72052, 34.16269, 34.604, 35.04144, 35.50077, 35.93656, 36.37506, 36.82023, 37.26342, 37.70502, 38.14945, 38.59152, 39.03771, 39.47901, 39.92461, 40.36414, 40.80781, 41.25428, 41.6931, 42.13469, 42.57282, 43.01831, 43.46041, 43.90174, 44.33874, 44.78026, 45.21801, 45.66176, 46.09533, 46.53811, 46.9725, 47.40696, 47.84195, 48.28061, 48.71304, 49.1466, 49.58126, 50.00482, 50.44228, 50.87165, 51.30117, 51.73439, 52.16194, 52.58038, 53.01602, 53.44421, 53.86809, 54.28808, 54.70689, 55.13004, 55.54849, 55.97355, 56.38838, 56.79826, 57.22806, 57.63537, 58.04899, 58.45993, 58.87694, 59.28334, 59.69365, 60.10462, 60.50689, 60.92053, 61.31914, 61.72463, 62.11896, 62.52295, 62.92568, 63.32497, 63.72648, 64.11418, 64.51223, 64.90428, 65.29777, 65.68106, 66.07122, 66.46024, 66.84258, 67.2288, 67.62034, 67.9984, 68.38331, 68.7479, 69.14391, 69.51381, 69.89029, 70.26415, 70.63729, 71.00797, 71.37879, 71.75265, 72.12494, 72.48273, 72.85086, 73.21501, 73.57664, 73.93506, 74.29445, 74.64832, 75.01395, 75.37629, 75.74131, 76.06679, 76.39962, 76.75977, 77.10653, 77.45647, 77.79773, 78.15308, 78.49065, 78.83663, 79.17577, 79.52041, 79.84789, 80.18466, 80.52013, 80.85467, 81.17969, 81.52161, 81.85233, 82.1756, 82.49329, 82.83146, 83.14359, 83.47777, 83.79487, 84.10968, 84.42749, 84.74511, 85.06303, 85.38384, 85.69256, 86.00663, 86.30816, 86.62556, 86.91648, 87.24067, 87.54285, 87.84368, 88.14325, 88.44446, 88.75524, 89.05019, 89.34085, 89.63443, 89.92918, 90.20956, 90.51693, 90.80653, 91.09319, 91.38013, 91.66794, 91.96231, 92.23865, 92.50818, 92.80587, 93.08887, 93.36391, 93.63823, 93.89955, 94.19181, 94.46941, 94.73439, 94.99956, 95.28116, 95.53494, 95.81467, 96.08157, 96.35242, 96.59919, 96.86464, 97.12854, 97.38136, 97.64384, 97.91152, 98.14753, 98.41387, 98.67443, 98.9061, 99.17898, 99.41421, 99.66087, 99.90495, 100.157, 100.397, 100.6379, 100.8871, 101.135, 101.3654, 101.5962, 101.8319, 102.0745, 102.2962, 102.5402, 102.7697, 102.9991, 103.2294, 103.4587, 103.6958, 103.9002, 104.133, 104.3283, 104.6157, 104.8096, 105.0291, 105.2453, 105.4602, 105.6768, 105.9003, 106.1005, 106.3176, 106.5301, 106.7394, 106.9303, 107.1504, 107.3539, 107.5614, 107.7609, 107.959, 108.1627, 108.3633, 108.5853, 108.7585, 108.9646, 109.1283, 109.3687, 109.5477, 109.7369, 109.9219, 110.1199, 110.2942, 110.4723, 110.6789, 110.8478, 111.0354, 111.2121, 111.3971, 111.5689, 111.7617, 111.9386, 112.1281, 112.2809, 112.472, 112.6432, 112.8058, 112.9804, 113.1497, 113.3036, 113.4663, 113.6617, 113.8046, 113.9457, 114.1448, 114.2835, 114.4546, 114.6269, 114.7851, 114.9273, 115.0646, 115.2453, 115.376, 115.541, 115.7042, 115.8489, 115.9841, 116.1428, 116.2816, 116.4167, 116.5695, 116.7082, 116.8502, 117.017, 117.1399, 117.2843, 117.4171, 117.5699, 117.6751, 118.127, 117.9429, 118.671, 118.2181, 118.9243, 118.463, 119.1765, 118.735, 119.4573, 118.9601, 119.7083, 119.2354, 119.9376, 119.4628, 120.1984, 119.7061, 120.4048, 119.92, 120.6493, 120.1586, 120.8645, 120.3862, 121.1027, 120.6113, 121.316, 120.819, 121.5341, 121.0077, 123.8236, 121.2636, 126.1908, 122.1, 125.8193, 122.3534, 126.0863, 122.6015, 126.3452, 122.8461, 126.5891, 122.4031, 123.191, 126.9555, 123.3915, 127.2085, 123.6498, 127.4473, 123.2002, 128.3474, 123.434, 128.5895, 123.6379, 124.4053, 128.9475, 123.9586, 129.1564, 124.1531, 129.3802, 124.3498, 124.4558, 129.7233, 124.6577, 129.94, 125.4928, 129.5082, 125.733, 125.1261, 129.7937, 126.0055, 130.0125, 126.1761, 125.5634, 130.3098, 126.4654, 130.4811, 125.9544, 126.7347, 130.7914, 126.9161, 130.9927, 126.3847, 127.1692, 131.2776, 126.621, 132.1624, 131.5374, 126.8585, 132.4258, 127.0338, 127.1083, 132.7336, 127.2584, 132.8996, 132.9857, 127.4769, 133.1416, 128.3853, 127.7033, 133.3804, 128.5662, 127.9084, 133.6403, 128.8043, 128.1166, 133.8942, 129.0191, 133.2971, 133.3616, 129.2203, 133.5503, 133.6055, 129.3869, 133.7714, 128.7805, 129.582, 133.9741, 128.9977, 129.7949, 134.1932, 129.1635, 129.9704, 135.1675, 134.4588, 129.371, 135.3755, 134.6261, 129.5263, 
1e-100, 0.2189761, 0.4414849, 0.668591, 0.8995527, 1.134261, 1.373715, 1.616674, 1.864122, 2.115748, 2.370915, 2.630609, 2.894509, 3.161868, 3.43383, 3.70978, 3.98933, 4.273604, 4.561759, 4.853086, 5.149691, 5.449543, 5.752896, 6.060955, 6.372535, 6.687649, 7.007291, 7.330588, 7.657573, 7.988013, 8.322832, 8.660778, 9.001611, 9.347463, 9.6965, 10.0488, 10.40461, 10.76441, 11.12623, 11.49185, 11.86118, 12.23408, 12.60962, 12.9883, 13.37055, 13.75497, 14.14165, 14.53462, 14.92755, 15.32333, 15.72389, 16.12577, 16.53134, 16.93749, 17.34668, 17.75944, 18.17531, 18.59, 19.01059, 19.43167, 19.8565, 20.283, 20.71071, 21.14275, 21.57339, 22.0078, 22.44559, 22.88139, 23.3221, 23.76381, 24.20743, 24.64969, 25.09811, 25.54679, 25.99285, 26.44577, 26.89748, 27.35033, 27.80414, 28.25888, 28.71661, 29.17486, 29.63261, 30.09262, 30.5528, 31.01568, 31.47747, 31.94009, 32.4036, 32.86934, 33.3358, 33.8005, 34.26736, 34.73127, 35.20199, 35.66966, 36.13663, 36.60117, 37.07209, 37.53385, 38.02345, 38.48464, 38.95421, 39.42331, 39.89256, 40.36333, 40.8316, 41.30041, 41.77099, 42.24203, 42.71076, 43.17806, 43.65075, 44.11694, 44.58125, 45.057, 45.51993, 45.99144, 46.45548, 46.92585, 47.39221, 47.85792, 48.32265, 48.78726, 49.25256, 49.71569, 50.17892, 50.64156, 51.10296, 51.56565, 52.02556, 52.48588, 52.95236, 53.40164, 53.86147, 54.32174, 54.78125, 55.23374, 55.68894, 56.14584, 56.60482, 57.05812, 57.50681, 57.9512, 58.40552, 58.85643, 59.30272, 59.7464, 60.19763, 60.63484, 61.09442, 61.53598, 61.97844, 62.42499, 62.86265, 63.30302, 63.73482, 64.18406, 64.62159, 65.06421, 65.48873, 65.92476, 66.35889, 66.79439, 67.22757, 67.65565, 68.08623, 68.50586, 68.94119, 69.36393, 69.79553, 70.21623, 70.64533, 71.06092, 71.48037, 71.91154, 72.32228, 72.74489, 73.16796, 73.5779, 73.993, 74.40565, 74.81054, 75.22978, 75.6431, 76.05491, 76.45795, 76.87595, 77.28416, 77.6864, 78.09066, 78.49738, 78.89974, 79.29675, 79.70216, 80.10725, 80.50267, 80.90976, 81.31957, 81.68588, 82.0546, 82.46856, 82.85797, 83.25692, 83.65007, 84.03709, 84.42706, 84.82212, 85.2055, 85.58124, 85.96734, 86.3633, 86.73458, 87.118, 87.50491, 87.88149, 88.26969, 88.63863, 89.01013, 89.39454, 89.76139, 90.14683, 90.51136, 90.87405, 91.25169, 91.62191, 91.99241, 92.36376, 92.72889, 93.0766, 93.44922, 93.8122, 94.17165, 94.53402, 94.90337, 95.25598, 95.60529, 95.97055, 96.32931, 96.67154, 97.0355, 97.37778, 97.74767, 98.0937, 98.43446, 98.7805, 99.11855, 99.48057, 99.83119, 100.183, 100.5142, 100.8556, 101.181, 101.5466, 101.8614, 102.2058, 102.5336, 102.8702, 103.208, 103.5465, 103.8722, 104.2175, 104.5324, 104.8504, 105.1852, 105.4943, 105.8514, 106.1554, 106.4848, 106.8027, 107.1131, 107.435, 107.7722, 108.0559, 108.3753, 108.6833, 109.0048, 109.3196, 109.628, 109.938, 110.2384, 110.5461, 110.8579, 111.1449, 111.4702, 111.7534, 112.059, 112.3364, 112.6627, 112.9268, 113.2331, 113.528, 113.8144, 114.1192, 114.3931, 114.6981, 114.9712, 115.2579, 115.5306, 115.8584, 116.1353, 116.3993, 116.676, 116.9541, 117.215, 117.5072, 117.7605, 118.0472, 118.3186, 118.5839, 118.8343, 119.1188, 119.3637, 119.6294, 119.8965, 120.1507, 120.4238, 120.6764, 120.959, 121.1895, 121.435, 121.6932, 121.9481, 122.2008, 122.4231, 122.6859, 122.9207, 123.1757, 123.4038, 123.6568, 123.9033, 124.1308, 124.3719, 124.5922, 124.8119, 125.0758, 125.2953, 125.531, 125.7457, 125.9859, 126.2143, 126.4114, 126.6541, 126.8508, 127.0863, 127.297, 127.5161, 127.7096, 127.9223, 128.1687, 128.3925, 128.5651, 128.79, 128.9964, 129.2027, 129.4042, 129.5913, 129.7928, 129.9879, 130.2084, 130.3921, 130.5683, 130.7999, 130.9494, 131.1662, 131.3531, 131.5243, 131.7324, 131.8826, 132.1004, 132.2772, 132.4367, 132.6284, 132.8119, 133.2607, 133.1421, 133.8508, 133.4888, 134.2008, 133.8121, 134.55, 134.1546, 134.8919, 134.4749, 135.1997, 134.8136, 135.5265, 135.108, 135.8593, 135.4314, 136.1412, 135.7378, 136.4335, 136.0169, 136.7425, 136.3283, 137.0765, 136.5961, 137.3336, 136.8952, 137.588, 137.1563, 140.4075, 137.4456, 143.2352, 138.3524, 143.0159, 138.6623, 143.3477, 139.0048, 143.6861, 139.3176, 144.0225, 139.032, 139.7444, 144.4821, 140.0485, 144.8322, 140.3486, 145.1056, 139.9809, 146.0568, 140.3204, 146.3896, 140.5916, 141.3207, 146.8216, 140.971, 147.0661, 141.2701, 147.3884, 141.5111, 141.6524, 147.8033, 141.8755, 148.0977, 142.7948, 147.7513, 143.0474, 142.4875, 148.1648, 143.4035, 148.4547, 143.6422, 143.0979, 148.8015, 143.984, 149.0638, 143.5844, 144.3438, 149.4548, 144.5193, 149.7398, 144.1691, 144.889, 150.0457, 144.4717, 150.95, 150.3959, 144.761, 151.3105, 145.0203, 145.0784, 151.5916, 145.3089, 151.8761, 151.9728, 145.5844, 152.2179, 146.4431, 145.8364, 152.5359, 146.715, 146.1358, 152.8644, 146.9725, 146.4077, 153.1762, 147.2863, 152.6579, 152.7573, 147.5329, 152.9392, 153.0608, 147.7921, 153.2398, 147.2709, 147.9938, 153.5172, 147.5584, 148.2816, 153.7864, 147.7319, 148.5054, 154.7836, 154.1119, 148.0019, 155.0552, 154.394, 148.2507, 
1e-100, 0.2341118, 0.4732834, 0.7163695, 0.9644683, 1.216968, 1.47341, 1.735091, 2.000762, 2.270646, 2.54583, 2.824762, 3.108355, 3.396782, 3.689043, 3.985848, 4.287432, 4.592802, 4.902897, 5.217233, 5.535613, 5.858562, 6.186126, 6.517319, 6.852112, 7.192642, 7.535887, 7.883094, 8.23553, 8.591286, 8.950259, 9.314392, 9.68182, 10.0523, 10.42794, 10.80691, 11.18871, 11.57472, 11.96409, 12.35723, 12.75243, 13.15279, 13.55682, 13.9625, 14.37161, 14.78348, 15.20007, 15.61852, 16.03913, 16.4645, 16.89269, 17.32187, 17.75428, 18.1906, 18.62789, 19.06868, 19.51196, 19.95677, 20.40536, 20.85461, 21.30666, 21.76336, 22.21872, 22.67631, 23.13682, 23.60055, 24.06633, 24.52947, 24.9971, 25.46866, 25.93946, 26.4135, 26.88537, 27.36381, 27.84152, 28.31958, 28.79866, 29.28078, 29.76405, 30.24761, 30.72975, 31.21747, 31.70391, 32.19471, 32.6825, 33.1713, 33.66002, 34.15069, 34.64456, 35.13705, 35.62774, 36.1242, 36.62006, 37.11294, 37.61229, 38.10676, 38.5993, 39.09707, 39.59448, 40.08115, 40.60181, 41.0903, 41.59002, 42.08642, 42.58139, 43.08103, 43.58041, 44.08058, 44.57762, 45.0773, 45.57837, 46.07496, 46.57466, 47.07014, 47.56856, 48.06989, 48.56238, 49.06573, 49.56193, 50.06279, 50.55738, 51.05165, 51.5535, 52.0485, 52.54514, 53.03965, 53.53761, 54.03401, 54.53132, 55.02043, 55.52053, 56.01585, 56.50623, 56.99115, 57.49096, 57.98279, 58.47956, 58.96971, 59.45827, 59.95154, 60.44185, 60.93439, 61.4258, 61.90745, 62.39825, 62.8834, 63.37098, 63.86224, 64.3506, 64.82674, 65.31953, 65.8137, 66.29687, 66.78656, 67.26094, 67.74936, 68.23141, 68.71811, 69.19733, 69.67975, 70.15882, 70.63808, 71.12657, 71.59934, 72.07439, 72.54732, 73.04796, 73.51497, 73.99285, 74.47452, 74.94935, 75.41823, 75.90283, 76.37396, 76.84508, 77.33625, 77.79652, 78.26501, 78.74009, 79.20314, 79.68302, 80.16154, 80.63411, 81.10514, 81.57198, 82.04397, 82.51309, 82.98199, 83.46077, 83.91266, 84.38884, 84.85183, 85.31282, 85.78398, 86.24519, 86.71651, 87.18882, 87.66297, 88.13596, 88.56348, 88.98551, 89.46594, 89.93953, 90.41074, 90.88216, 91.3402, 91.79441, 92.24936, 92.69259, 93.16284, 93.62145, 94.08704, 94.5489, 95.00308, 95.45677, 95.90572, 96.35008, 96.81665, 97.27354, 97.72466, 98.19429, 98.62945, 99.07439, 99.51639, 99.9921, 100.4238, 100.8873, 101.3446, 101.7639, 102.2233, 102.6637, 103.1093, 103.5552, 103.9813, 104.4292, 104.89, 105.3138, 105.7853, 106.2132, 106.6281, 107.0716, 107.519, 107.9329, 108.4019, 108.8177, 109.2627, 109.702, 110.1292, 110.5472, 111.005, 111.4062, 111.8592, 112.2878, 112.6956, 113.1395, 113.5613, 113.9955, 114.3896, 114.8095, 115.2319, 115.6753, 116.0606, 116.4889, 116.9005, 117.3441, 117.7572, 118.1443, 118.5806, 118.9775, 119.3956, 119.7911, 120.1864, 120.6173, 120.9898, 121.4034, 121.8298, 122.1945, 122.6198, 122.9916, 123.3952, 123.8043, 124.1788, 124.595, 124.9647, 125.3801, 125.7472, 126.1363, 126.5247, 126.8895, 127.2874, 127.6616, 128.0536, 128.4196, 128.8037, 129.1977, 129.5237, 129.9165, 130.2669, 130.6373, 131.052, 131.4083, 131.7659, 132.122, 132.4725, 132.849, 133.2034, 133.5552, 133.9025, 134.245, 134.5884, 134.9567, 135.2892, 135.6489, 135.9681, 136.3268, 136.6506, 136.9993, 137.3645, 137.6732, 138.0242, 138.2967, 138.7022, 138.9923, 139.312, 139.6515, 139.9584, 140.2678, 140.573, 140.9259, 141.1984, 141.5321, 141.8455, 142.1525, 142.4484, 142.7561, 143.0844, 143.4022, 143.6627, 143.9809, 144.2679, 144.5485, 144.8639, 145.1431, 145.4165, 145.6797, 146.0181, 146.2572, 146.539, 146.8296, 147.1199, 147.3874, 147.6618, 147.9408, 148.1835, 148.4153, 148.7401, 148.9697, 149.2791, 149.5133, 149.7484, 150.0023, 150.276, 150.5243, 150.7735, 150.9959, 151.2904, 151.4913, 151.7817, 151.9808, 152.2472, 152.4767, 152.7106, 152.9085, 153.1803, 153.6313, 153.5993, 154.352, 154.0532, 154.7752, 154.5142, 155.2684, 154.9438, 155.62, 155.3933, 156.0825, 155.7632, 156.4981, 156.1972, 156.9085, 156.5823, 157.2915, 156.9834, 157.7019, 157.3383, 158.0993, 157.7413, 158.5128, 158.1426, 158.8216, 158.4231, 159.1586, 158.7773, 162.7268, 159.1895, 166.2266, 160.1836, 166.1741, 160.5899, 166.6169, 160.9728, 167.0727, 161.3797, 167.4389, 161.2347, 161.9254, 168.0727, 162.3532, 168.481, 162.7009, 168.9097, 162.507, 169.8752, 162.8988, 170.213, 163.262, 163.9586, 170.8416, 163.7202, 171.2141, 164.1036, 171.5477, 164.4537, 164.6117, 172.1226, 164.8984, 172.4893, 165.9037, 172.2817, 166.1652, 165.7216, 172.7887, 166.6476, 173.1015, 166.9519, 166.5094, 173.6093, 167.3612, 173.9471, 167.1356, 167.833, 174.4008, 168.1419, 174.7431, 167.8173, 168.5648, 175.1723, 168.2337, 176.1423, 175.6229, 168.6162, 176.5462, 168.9056, 168.9981, 177.0087, 169.2563, 177.2766, 177.3789, 169.631, 177.7372, 170.5127, 169.9259, 178.1028, 170.8449, 170.3301, 178.5185, 171.2291, 170.6627, 178.9053, 171.5068, 178.461, 178.6229, 171.8838, 178.863, 178.9637, 172.1742, 179.1961, 171.7877, 172.4759, 179.5866, 172.0711, 172.7798, 179.9069, 172.3571, 173.0781, 180.9394, 180.2898, 172.6585, 181.1748, 180.6735, 173.0056, 
1e-100, 0.2515362, 0.5074022, 0.769017, 1.034959, 1.305702, 1.581703, 1.861909, 2.147772, 2.438147, 2.732939, 3.033204, 3.337989, 3.646958, 3.961634, 4.280367, 4.603829, 4.932615, 5.265658, 5.602754, 5.945245, 6.291806, 6.642336, 6.998673, 7.358337, 7.722572, 8.091474, 8.464029, 8.841295, 9.222381, 9.608942, 9.997815, 10.3912, 10.78852, 11.18987, 11.59509, 12.00501, 12.41636, 12.83331, 13.25371, 13.67578, 14.10341, 14.53248, 14.96712, 15.40381, 15.84306, 16.2873, 16.73405, 17.18413, 17.63539, 18.09242, 18.54975, 19.01018, 19.47342, 19.93875, 20.40816, 20.88003, 21.35102, 21.82839, 22.30689, 22.78908, 23.27154, 23.75613, 24.24398, 24.7327, 25.22312, 25.71913, 26.21021, 26.70764, 27.20591, 27.70738, 28.20771, 28.71376, 29.21632, 29.72218, 30.23342, 30.74115, 31.25136, 31.76508, 32.27852, 32.78922, 33.30656, 33.82262, 34.34147, 34.86197, 35.38206, 35.90399, 36.42146, 36.94806, 37.47088, 37.99466, 38.51962, 39.04494, 39.57101, 40.09941, 40.62625, 41.15488, 41.68431, 42.21489, 42.73594, 43.2909, 43.81247, 44.34591, 44.87564, 45.4101, 45.94411, 46.47839, 47.01185, 47.54712, 48.08353, 48.6198, 49.15048, 49.69226, 50.22517, 50.76358, 51.30125, 51.83413, 52.38341, 52.91443, 53.45442, 53.98981, 54.53343, 55.07416, 55.61337, 56.14801, 56.69149, 57.2347, 57.76911, 58.31944, 58.85769, 59.3992, 59.94306, 60.48831, 61.02394, 61.57148, 62.1114, 62.65557, 63.19932, 63.74445, 64.28816, 64.83022, 65.3773, 65.92281, 66.46535, 67.01117, 67.56412, 68.10504, 68.65268, 69.18927, 69.73593, 70.30859, 70.84792, 71.39577, 71.94195, 72.48913, 73.04102, 73.59491, 74.14228, 74.68796, 75.25308, 75.78708, 76.33681, 76.88707, 77.43783, 77.99928, 78.55347, 79.11306, 79.65833, 80.20345, 80.76854, 81.33821, 81.88582, 82.44124, 82.99387, 83.54163, 84.10365, 84.66372, 85.21037, 85.7742, 86.32959, 86.89813, 87.4491, 88.01184, 88.58338, 89.13861, 89.69566, 90.27066, 90.82428, 91.37379, 91.92662, 92.50451, 93.06434, 93.63285, 94.20636, 94.76232, 95.32502, 95.90191, 96.46837, 97.04065, 97.57156, 98.12481, 98.69339, 99.27616, 99.8398, 100.4092, 100.9596, 101.5381, 102.0784, 102.663, 103.2595, 103.8138, 104.3762, 104.951, 105.5105, 106.0734, 106.6402, 107.2268, 107.7899, 108.3613, 108.9357, 109.4791, 110.0442, 110.606, 111.1757, 111.78, 112.334, 112.8727, 113.4452, 113.9935, 114.5851, 115.1519, 115.709, 116.262, 116.8346, 117.3914, 117.978, 118.5221, 119.0888, 119.6331, 120.2166, 120.7604, 121.3611, 121.9138, 122.444, 123.0109, 123.5492, 124.137, 124.6774, 125.2104, 125.8013, 126.3655, 126.8685, 127.4541, 128.0016, 128.5337, 129.1117, 129.6352, 130.1694, 130.7332, 131.2674, 131.8557, 132.3375, 132.8766, 133.4468, 133.9606, 134.5426, 135.0244, 135.578, 136.1314, 136.644, 137.169, 137.7206, 138.2212, 138.7438, 139.2695, 139.8031, 140.3242, 140.8204, 141.3915, 141.852, 142.3983, 142.8836, 143.4118, 143.9317, 144.4356, 144.9583, 145.4003, 145.9792, 146.4257, 146.949, 147.4461, 147.9157, 148.4261, 148.912, 149.4253, 149.883, 150.3698, 150.8631, 151.283, 151.8856, 152.318, 152.7819, 153.2779, 153.7067, 154.2033, 154.6271, 155.1201, 155.5738, 156.0248, 156.464, 156.9643, 157.383, 157.8489, 158.2771, 158.7226, 159.1871, 159.6084, 160.0983, 160.5068, 160.9052, 161.353, 161.789, 162.2185, 162.607, 163.0474, 163.443, 163.8687, 164.2925, 164.7125, 165.146, 165.5135, 165.9558, 166.2881, 166.7185, 167.134, 167.4926, 167.9047, 168.2926, 168.6764, 169.0544, 169.4283, 169.8425, 170.1701, 170.5579, 170.9197, 171.2777, 171.656, 171.9928, 172.3723, 172.7947, 173.1017, 173.4786, 173.841, 174.1463, 174.4992, 174.8531, 175.1482, 175.4877, 175.8894, 176.1688, 176.4678, 176.8838, 177.1345, 177.4932, 177.7959, 178.1462, 178.4577, 178.6867, 179.1075, 179.3618, 179.6166, 179.9524, 180.2734, 180.5339, 181.1077, 181.1003, 181.9222, 181.7543, 182.4074, 182.2618, 183.0123, 182.8388, 183.5886, 183.3492, 184.0974, 183.8871, 184.6847, 184.4671, 185.1086, 184.8978, 185.6821, 185.3972, 186.1556, 185.9221, 186.6871, 186.4052, 187.0615, 186.9933, 187.5816, 187.2185, 188.0078, 187.7479, 192.3807, 188.2147, 196.8686, 189.2604, 196.9181, 189.7371, 197.4845, 190.3287, 198.0097, 190.8219, 198.5486, 190.8391, 191.4689, 199.3222, 192.0537, 199.8349, 192.4612, 200.3285, 192.4814, 201.398, 192.8905, 201.8996, 193.3817, 194.1038, 202.5825, 194.0387, 203.1165, 194.4406, 203.5004, 194.8929, 195.0552, 204.2295, 195.4822, 204.6995, 196.4659, 204.6236, 196.8281, 196.4951, 205.2519, 197.4601, 205.6528, 197.8294, 197.4013, 206.2947, 198.4062, 206.7043, 198.2729, 198.8815, 207.299, 199.3643, 207.6897, 199.0926, 199.8037, 208.2476, 199.595, 209.2051, 208.8448, 200.1602, 209.7943, 200.3889, 200.5286, 210.3362, 200.9182, 210.6753, 210.8205, 201.3709, 211.1497, 202.2578, 201.8162, 211.7093, 202.7006, 202.2209, 212.1811, 203.1345, 202.6186, 212.7176, 203.5213, 212.367, 212.4964, 203.9218, 212.8692, 212.99, 204.343, 213.279, 203.9182, 204.796, 213.7579, 204.383, 205.0699, 214.1234, 204.7481, 205.402, 215.1675, 214.6715, 205.1999, 215.554, 215.0245, 205.4755, 
1e-100, 0.2693868, 0.5443023, 0.8239528, 1.109686, 1.400195, 1.695498, 1.996891, 2.302767, 2.613851, 2.930468, 3.251636, 3.578613, 3.910419, 4.246716, 4.588697, 4.935366, 5.286748, 5.6433, 6.004658, 6.370462, 6.741866, 7.117619, 7.497379, 7.882135, 8.272107, 8.665209, 9.063926, 9.466893, 9.873611, 10.28458, 10.70084, 11.12055, 11.54369, 11.97219, 12.40478, 12.83987, 13.2789, 13.72329, 14.16989, 14.61916, 15.07461, 15.53242, 15.99393, 16.45831, 16.92669, 17.39759, 17.87199, 18.34902, 18.8301, 19.31537, 19.79962, 20.28911, 20.78264, 21.27727, 21.77381, 22.27568, 22.77741, 23.28358, 23.79078, 24.30241, 24.81586, 25.33024, 25.84509, 26.36637, 26.88888, 27.41483, 27.93562, 28.46386, 28.99547, 29.52794, 30.06071, 30.59569, 31.13548, 31.67422, 32.21531, 32.75846, 33.30163, 33.84819, 34.39736, 34.94272, 35.49412, 36.04726, 36.6025, 37.15669, 37.71116, 38.27006, 38.82748, 39.3913, 39.95085, 40.51342, 41.08076, 41.64666, 42.21559, 42.78321, 43.35209, 43.92527, 44.49591, 45.07014, 45.63025, 46.23087, 46.79841, 47.37896, 47.9546, 48.53434, 49.1115, 49.70142, 50.28885, 50.87212, 51.45651, 52.04152, 52.63384, 53.22726, 53.81197, 54.40661, 54.9948, 55.59046, 56.19446, 56.78253, 57.3904, 57.98267, 58.58155, 59.18831, 59.79446, 60.4049, 61.00365, 61.61065, 62.21452, 62.8347, 63.4419, 64.05514, 64.66657, 65.28953, 65.88768, 66.50986, 67.13765, 67.75809, 68.3769, 68.9974, 69.62878, 70.26149, 70.88553, 71.52013, 72.14467, 72.78574, 73.40904, 74.04163, 74.68636, 75.33143, 75.95975, 76.61738, 77.23905, 77.89594, 78.55303, 79.194, 79.831, 80.47628, 81.13342, 81.78709, 82.45103, 83.10581, 83.76987, 84.43607, 85.09422, 85.75054, 86.41822, 87.09758, 87.76739, 88.44017, 89.11125, 89.79043, 90.4455, 91.11452, 91.80087, 92.46756, 93.15155, 93.84456, 94.5158, 95.21313, 95.87942, 96.59569, 97.28967, 97.98383, 98.68507, 99.36444, 100.0617, 100.7293, 101.4376, 102.1549, 102.8471, 103.5602, 104.2668, 104.9601, 105.6737, 106.3731, 107.0947, 107.8102, 108.5223, 109.2325, 109.9137, 110.612, 111.3443, 112.0609, 112.8171, 113.5083, 114.2064, 114.9282, 115.6679, 116.3792, 117.1427, 117.8279, 118.5579, 119.2942, 120.0116, 120.7325, 121.4557, 122.2267, 122.9346, 123.6264, 124.4027, 125.0983, 125.8451, 126.584, 127.3045, 128.0488, 128.772, 129.4524, 130.232, 130.9729, 131.6857, 132.4488, 133.157, 133.8171, 134.6048, 135.3533, 136.0736, 136.834, 137.5047, 138.2737, 138.9751, 139.7155, 140.4163, 141.1427, 141.8845, 142.623, 143.3415, 144.0656, 144.7613, 145.511, 146.2266, 146.9967, 147.6709, 148.3885, 149.1285, 149.789, 150.5951, 151.2084, 151.9578, 152.6591, 153.3938, 154.0631, 154.792, 155.4685, 156.2149, 156.9119, 157.6358, 158.3061, 158.961, 159.6776, 160.3516, 161.0354, 161.7862, 162.4062, 163.1034, 163.8256, 164.4734, 165.149, 165.817, 166.5069, 167.1799, 167.8075, 168.523, 169.1566, 169.863, 170.5296, 171.1693, 171.8272, 172.4673, 173.1304, 173.7374, 174.4282, 175.0493, 175.6949, 176.3678, 176.9451, 177.5917, 178.2306, 178.8484, 179.3415, 180.2079, 180.7298, 181.343, 181.9254, 182.5796, 183.1588, 183.7424, 184.3773, 184.9071, 185.524, 186.1645, 186.6824, 187.3706, 187.8544, 188.4499, 189.0249, 189.5481, 190.1616, 190.7139, 191.326, 191.7844, 192.4579, 192.9535, 193.4999, 194.0748, 194.5887, 195.1183, 195.6212, 196.2199, 196.6604, 197.2539, 197.7489, 198.2932, 198.761, 199.2926, 199.825, 200.333, 200.8405, 201.3011, 201.8038, 202.2399, 202.7876, 203.2722, 203.738, 204.1745, 204.7095, 205.1238, 205.6427, 206.0923, 206.4855, 207.004, 207.4043, 207.8904, 208.3462, 208.727, 209.2468, 209.6235, 210.06, 210.4952, 210.9447, 211.3048, 211.7616, 212.1312, 212.5706, 212.8897, 213.3949, 213.7725, 214.2317, 214.5784, 214.9703, 215.3355, 215.7534, 216.1648, 216.5243, 216.7776, 217.4042, 217.5717, 218.3698, 218.3312, 219.0908, 219.0324, 219.7729, 219.5833, 220.5186, 220.3197, 221.1569, 221.0683, 221.8064, 221.6975, 222.444, 222.2957, 223.1191, 222.998, 223.6319, 223.5572, 224.3292, 224.2055, 224.9347, 224.8407, 225.4581, 225.1853, 226.0693, 225.8679, 231.558, 226.3804, 237.2579, 227.6023, 237.5049, 228.2884, 238.1702, 228.9041, 238.7652, 229.5115, 239.5364, 229.6858, 230.3684, 240.5322, 230.987, 241.1198, 231.6245, 241.7635, 231.7353, 242.949, 232.2879, 243.5578, 232.8959, 233.6382, 244.4258, 233.6168, 244.9882, 234.2335, 245.6043, 234.7225, 234.9794, 246.3218, 235.5733, 247.0514, 236.4837, 247.0243, 236.9876, 236.6971, 247.7762, 237.6857, 248.3813, 238.2214, 237.8653, 248.9891, 238.9067, 249.5705, 238.8247, 239.4888, 250.3607, 240.0274, 250.7248, 239.8339, 240.6068, 251.4913, 240.5373, 252.4276, 252.1585, 241.1474, 253.3071, 241.5327, 241.7015, 253.8838, 242.0333, 254.2947, 254.5564, 242.6659, 254.903, 243.4486, 243.2152, 255.612, 244.1966, 243.6337, 256.0777, 244.6783, 244.1515, 256.7036, 245.1798, 256.6124, 256.7413, 245.5321, 257.1291, 257.2824, 246.0539, 257.6344, 245.8887, 246.4902, 258.1581, 246.2237, 246.9188, 258.6189, 246.6909, 247.3372, 259.6714, 259.1997, 247.224, 260.2531, 259.8042, 247.642, 
1e-100, 0.2884452, 0.5822376, 0.8822939, 1.187327, 1.498187, 1.814609, 2.136124, 2.463933, 2.796689, 3.134466, 3.478673, 3.827523, 4.181741, 4.54154, 4.906108, 5.276314, 5.651902, 6.032056, 6.417115, 6.808141, 7.203286, 7.60352, 8.009224, 8.41872, 8.833731, 9.253333, 9.677131, 10.10551, 10.53939, 10.97686, 11.41881, 11.86543, 12.31592, 12.7708, 13.22956, 13.69343, 14.16077, 14.63107, 15.10655, 15.58493, 16.06892, 16.55408, 17.0455, 17.53867, 18.03484, 18.5365, 19.04123, 19.54815, 20.05764, 20.57501, 21.09125, 21.61182, 22.13509, 22.66077, 23.19213, 23.72365, 24.25742, 24.79869, 25.33946, 25.88442, 26.43115, 26.98143, 27.5325, 28.08862, 28.64561, 29.2089, 29.76684, 30.33324, 30.90169, 31.47129, 32.04653, 32.62013, 33.19893, 33.77717, 34.36507, 34.94738, 35.53466, 36.12177, 36.71896, 37.31162, 37.90518, 38.50718, 39.10829, 39.71226, 40.32206, 40.9286, 41.5349, 42.15016, 42.76685, 43.38197, 44.00035, 44.61872, 45.24563, 45.87279, 46.50265, 47.12903, 47.7619, 48.39769, 49.03144, 49.70106, 50.33264, 50.9693, 51.61638, 52.2709, 52.9241, 53.57728, 54.23199, 54.88288, 55.5569, 56.22103, 56.88715, 57.55539, 58.22575, 58.90501, 59.58227, 60.26258, 60.94345, 61.62611, 62.32484, 63.00892, 63.7077, 64.40436, 65.10377, 65.80486, 66.50392, 67.22382, 67.93457, 68.63962, 69.34492, 70.0784, 70.80359, 71.53449, 72.2443, 72.9908, 73.71806, 74.46305, 75.20096, 75.9559, 76.70155, 77.44227, 78.19502, 78.96947, 79.71636, 80.48651, 81.23757, 82.00813, 82.78704, 83.56307, 84.32966, 85.11752, 85.92468, 86.69745, 87.47985, 88.26396, 89.08335, 89.88325, 90.69209, 91.50187, 92.31319, 93.11649, 93.94827, 94.76297, 95.57265, 96.41343, 97.24915, 98.09711, 98.91638, 99.75292, 100.602, 101.4333, 102.2683, 103.1386, 103.9735, 104.8317, 105.7142, 106.5681, 107.4325, 108.3382, 109.1847, 110.0603, 110.9541, 111.8033, 112.6873, 113.5664, 114.4492, 115.3407, 116.2262, 117.137, 118.0317, 118.937, 119.8689, 120.7612, 121.662, 122.5681, 123.4776, 124.3748, 125.3126, 126.2602, 127.1676, 128.0187, 128.9806, 129.9094, 130.8033, 131.7324, 132.7025, 133.6114, 134.5753, 135.5205, 136.4452, 137.3434, 138.2512, 139.2156, 140.1826, 141.1446, 142.0811, 143.0465, 143.9785, 144.8599, 145.8174, 146.7844, 147.7102, 148.7326, 149.6687, 150.5707, 151.5072, 152.5005, 153.3981, 154.4027, 155.3267, 156.2781, 157.2124, 158.1977, 159.1149, 160.102, 161.0674, 161.962, 162.9397, 163.8834, 164.8253, 165.7973, 166.6716, 167.6644, 168.6221, 169.4991, 170.4982, 171.3693, 172.3768, 173.2922, 174.3025, 175.141, 176.1095, 177.0026, 177.9854, 178.9003, 179.8103, 180.7395, 181.6595, 182.5843, 183.5375, 184.3703, 185.3095, 186.2325, 187.1068, 188.0682, 188.9515, 189.8427, 190.7389, 191.6552, 192.5527, 193.4516, 194.3336, 195.1929, 196.1071, 196.9695, 197.8244, 198.7078, 199.6365, 200.4629, 201.309, 202.1236, 203.034, 203.8891, 204.7376, 205.6302, 206.3364, 207.3222, 208.1273, 208.9395, 209.7724, 210.5994, 211.4346, 212.2397, 213.1252, 213.8465, 214.6626, 215.492, 216.2389, 216.998, 218.0032, 218.7014, 219.4779, 220.2109, 221.0129, 221.7132, 222.5529, 223.2609, 224.0576, 224.7138, 225.5951, 226.2711, 227.0384, 227.747, 228.4649, 229.2268, 229.8696, 230.7578, 231.3626, 232.0989, 232.8086, 233.493, 234.2046, 234.8596, 235.5548, 236.2119, 236.9216, 237.5351, 238.3215, 238.9742, 239.5553, 240.2741, 240.8621, 241.5948, 242.1978, 242.808, 243.4563, 244.0662, 244.7024, 245.2904, 245.9573, 246.6099, 247.1016, 247.7094, 248.2605, 248.9038, 249.4837, 250.0248, 250.6388, 251.1824, 251.8251, 252.3905, 252.9243, 253.4545, 254.0351, 254.5755, 255.0568, 255.6172, 256.2499, 256.7378, 257.1232, 257.8184, 258.2867, 258.8116, 259.3226, 259.8549, 260.3438, 260.6789, 261.3263, 261.7558, 262.2433, 262.6956, 263.154, 263.6099, 264.0553, 264.6763, 264.9653, 265.8319, 265.8278, 266.6821, 266.702, 267.5302, 267.5246, 268.3604, 268.2752, 269.167, 269.1972, 270.0046, 269.8743, 270.741, 270.7067, 271.5597, 271.5128, 272.2606, 272.2868, 272.9789, 272.9067, 273.788, 273.8818, 274.4839, 274.1679, 275.1459, 275.0255, 282.1507, 275.7667, 289.2521, 277.0541, 289.5021, 277.8614, 290.4939, 278.6684, 291.3008, 279.3399, 292.1096, 279.6584, 280.4809, 293.2848, 281.2177, 293.9279, 281.9979, 294.9058, 282.1697, 295.9991, 282.7849, 296.8699, 283.43, 284.1522, 297.9056, 284.4688, 298.5473, 285.0708, 299.2909, 285.7162, 286.0532, 300.3457, 286.6623, 301.0112, 287.7358, 301.2066, 288.3772, 288.148, 302.0991, 289.1975, 302.6793, 289.8763, 289.568, 303.7514, 290.6208, 304.3195, 290.6463, 291.51, 305.2774, 292.0183, 305.7319, 292.0492, 292.6566, 306.5847, 292.761, 307.7261, 307.4668, 293.3993, 308.373, 294.0058, 294.1167, 309.3638, 294.612, 309.8555, 309.9941, 295.2154, 310.5422, 296.2613, 295.874, 311.2907, 296.8207, 296.5344, 312.035, 297.5024, 297.0984, 312.597, 298.0058, 312.594, 312.8605, 298.7033, 313.3118, 313.5269, 299.126, 313.9128, 299.1629, 299.7675, 314.5031, 299.5216, 300.2056, 315.1346, 299.9921, 300.7043, 316.3619, 315.9381, 300.6709, 316.6982, 316.5754, 301.21, 
1e-100, 0.308125, 0.6221091, 0.9417185, 1.268068, 1.599621, 1.936974, 2.280623, 2.629481, 2.984473, 3.345037, 3.710801, 4.0833, 4.46085, 4.843252, 5.232649, 5.626468, 6.025461, 6.430802, 6.84045, 7.255663, 7.676853, 8.103003, 8.533108, 8.969727, 9.411434, 9.856317, 10.30745, 10.76333, 11.22348, 11.68953, 12.15997, 12.63455, 13.11201, 13.59791, 14.08583, 14.57778, 15.0756, 15.57707, 16.08265, 16.59111, 17.10733, 17.62554, 18.14736, 18.67413, 19.20619, 19.73861, 20.2777, 20.81939, 21.36633, 21.91873, 22.47099, 23.029, 23.59032, 24.15484, 24.72178, 25.29722, 25.87126, 26.45193, 27.0328, 27.62464, 28.2158, 28.80798, 29.40358, 30.00913, 30.61431, 31.22129, 31.82825, 32.44782, 33.06649, 33.69135, 34.3153, 34.94511, 35.58023, 36.21556, 36.85801, 37.5015, 38.15145, 38.8001, 39.45781, 40.11259, 40.77617, 41.44716, 42.11816, 42.78538, 43.4644, 44.144, 44.83152, 45.51878, 46.2113, 46.91108, 47.6132, 48.31916, 49.02723, 49.73232, 50.45055, 51.17247, 51.89941, 52.62409, 53.3395, 54.10554, 54.83542, 55.58308, 56.33157, 57.08838, 57.84014, 58.61093, 59.37647, 60.15056, 60.92733, 61.69559, 62.48743, 63.28269, 64.07924, 64.8751, 65.66396, 66.48555, 67.29763, 68.11803, 68.93993, 69.76595, 70.60761, 71.44246, 72.29023, 73.13182, 73.98722, 74.83731, 75.69279, 76.58933, 77.44679, 78.33153, 79.21057, 80.10298, 80.98763, 81.90034, 82.79123, 83.7175, 84.63748, 85.53847, 86.45726, 87.40799, 88.36324, 89.30738, 90.22251, 91.1572, 92.13813, 93.11131, 94.06375, 95.01786, 96.01842, 97.00353, 97.9934, 98.97852, 99.99052, 101.0141, 102.0065, 103.0085, 104.0167, 105.0637, 106.1452, 107.1502, 108.1855, 109.2271, 110.2588, 111.3208, 112.4039, 113.477, 114.5113, 115.5804, 116.6686, 117.7384, 118.8202, 119.9442, 121.0591, 122.1097, 123.246, 124.3884, 125.4749, 126.5943, 127.7292, 128.8817, 129.9819, 131.0957, 132.2536, 133.3574, 134.5295, 135.7151, 136.8819, 138.0414, 139.1974, 140.3579, 141.5312, 142.714, 143.8486, 145.0392, 146.2262, 147.442, 148.5938, 149.8305, 150.9613, 152.0976, 153.3002, 154.5309, 155.7309, 156.9638, 158.1791, 159.416, 160.5525, 161.8076, 163.032, 164.2033, 165.4402, 166.6843, 167.8734, 169.134, 170.3367, 171.5093, 172.7544, 174.0291, 175.2411, 176.464, 177.6931, 178.9168, 180.1396, 181.3879, 182.5879, 183.831, 185.0673, 186.3087, 187.4906, 188.6982, 189.9179, 191.1825, 192.4507, 193.5933, 194.8658, 196.0556, 197.3607, 198.5517, 199.6911, 200.9118, 202.1967, 203.3442, 204.6086, 205.8445, 206.9982, 208.243, 209.4164, 210.6078, 211.9078, 212.9686, 214.2394, 215.4287, 216.6026, 217.8622, 218.9414, 220.1686, 221.3367, 222.4689, 223.6097, 224.9056, 225.9624, 227.1729, 228.3056, 229.4662, 230.681, 231.752, 232.9052, 233.9912, 235.2001, 236.3246, 237.3898, 238.5888, 239.7105, 240.8177, 241.8756, 243.0254, 244.1448, 245.245, 246.2885, 247.4509, 248.4094, 249.6283, 250.6934, 251.7713, 252.8013, 253.8147, 254.8924, 255.8304, 257.0374, 258.0372, 258.9961, 260.1128, 261.081, 262.1055, 263.1452, 264.144, 264.9991, 266.0294, 267.2833, 268.1451, 269.082, 270.1088, 271.0479, 271.9351, 273.0136, 273.8586, 274.8987, 275.8547, 276.7561, 277.7998, 278.6553, 279.6018, 280.5002, 281.3674, 282.2787, 283.1713, 284.1029, 284.9416, 285.9151, 286.6228, 287.531, 288.4169, 289.2218, 290.06, 290.8971, 291.7554, 292.4787, 293.3922, 294.1703, 294.9735, 295.7616, 296.4713, 297.3168, 298.0324, 298.8749, 299.6293, 300.3327, 301.0963, 301.8484, 302.6446, 303.4121, 304.1332, 304.8227, 305.5563, 306.3693, 307.1264, 307.7541, 308.311, 309.1756, 309.8951, 310.3981, 311.2168, 311.8959, 312.5256, 313.1032, 313.8277, 314.4154, 315.143, 315.682, 316.3432, 317.0361, 317.4359, 318.141, 318.7467, 319.4479, 319.9348, 320.4758, 320.967, 321.6724, 322.2278, 322.7837, 323.2297, 323.809, 324.6317, 324.9974, 325.8505, 326.0288, 326.9559, 327.0324, 327.9419, 328.1358, 329.0237, 329.036, 330.0663, 330.0518, 330.9421, 330.9811, 331.969, 331.9618, 332.8272, 332.909, 333.7621, 333.8369, 334.58, 334.6737, 335.4299, 335.5399, 336.2876, 336.0695, 337.1017, 336.9891, 345.5227, 338.0143, 354.2145, 339.2434, 354.7914, 340.283, 355.8113, 341.1301, 356.8473, 342.2247, 357.8496, 342.4905, 343.3745, 359.1749, 344.2804, 360.2014, 345.1931, 361.024, 345.5448, 362.4598, 346.3643, 363.4789, 347.1371, 347.8933, 364.6021, 348.2629, 365.4444, 349.1339, 366.2089, 349.8025, 349.968, 367.5712, 350.9203, 368.2266, 352.027, 368.6367, 352.7784, 352.5483, 369.7227, 353.8884, 370.5407, 354.4064, 354.216, 371.6305, 355.5278, 372.2241, 355.7729, 356.3508, 373.2139, 356.9728, 374.043, 357.105, 357.8239, 374.964, 357.935, 376.1696, 376.0518, 358.8568, 377.0904, 359.4307, 359.5427, 377.9836, 360.3472, 378.6547, 378.8556, 360.871, 379.3507, 362.0654, 361.766, 380.376, 362.7753, 362.3882, 381.228, 363.2387, 363.0737, 382.0654, 364.0815, 381.9945, 382.2291, 364.7346, 382.8929, 383.0981, 365.3399, 383.5865, 365.3402, 365.9547, 384.2832, 365.932, 366.4833, 384.9832, 366.4697, 367.0548, 386.0964, 385.9104, 367.2684, 386.8347, 386.6015, 367.7041, 
1e-100, 0.328068, 0.6625148, 1.003524, 1.350154, 1.703598, 2.062854, 2.42806, 2.800233, 3.177661, 3.56121, 3.951318, 4.346594, 4.74808, 5.15587, 5.5687, 5.988284, 6.413396, 6.843594, 7.280302, 7.722843, 8.169409, 8.623226, 9.081849, 9.545515, 10.01604, 10.49089, 10.97155, 11.45709, 11.94899, 12.44521, 12.94592, 13.45395, 13.96521, 14.48285, 15.00591, 15.53431, 16.06564, 16.60254, 17.14624, 17.69384, 18.24584, 18.80354, 19.36684, 19.93328, 20.50421, 21.08165, 21.66669, 22.25193, 22.84049, 23.44492, 24.04432, 24.64981, 25.26222, 25.8795, 26.49999, 27.12462, 27.75734, 28.39489, 29.03471, 29.68249, 30.33377, 30.99028, 31.65008, 32.31839, 32.9884, 33.66975, 34.34722, 35.03324, 35.72844, 36.42597, 37.13562, 37.8412, 38.55192, 39.27067, 40.00883, 40.73313, 41.46772, 42.20706, 42.96649, 43.71578, 44.46714, 45.2305, 46.01244, 46.78853, 47.57776, 48.35838, 49.15218, 49.95644, 50.76495, 51.57846, 52.39801, 53.21964, 54.05742, 54.90418, 55.75297, 56.59808, 57.45854, 58.32413, 59.19285, 60.09874, 60.97692, 61.87647, 62.77301, 63.68371, 64.57565, 65.51995, 66.43539, 67.38327, 68.32125, 69.27164, 70.23955, 71.18161, 72.17438, 73.14786, 74.14684, 75.14859, 76.12862, 77.14514, 78.17124, 79.21209, 80.24531, 81.26018, 82.32957, 83.37075, 84.45306, 85.51367, 86.60706, 87.69688, 88.783, 89.87642, 91.02662, 92.15591, 93.25073, 94.38593, 95.529, 96.69299, 97.85942, 99.01784, 100.1981, 101.3761, 102.572, 103.7515, 104.9717, 106.1833, 107.4215, 108.6201, 109.8692, 111.1224, 112.3803, 113.6854, 114.9162, 116.2344, 117.5085, 118.7819, 120.0583, 121.3891, 122.7394, 124.0345, 125.3728, 126.6577, 127.9953, 129.3707, 130.7347, 132.0956, 133.4292, 134.8171, 136.2387, 137.6017, 139.0189, 140.4256, 141.858, 143.242, 144.6522, 146.0727, 147.5392, 149.0054, 150.429, 151.8458, 153.2456, 154.7335, 156.1699, 157.683, 159.1474, 160.6831, 162.1229, 163.6487, 165.1336, 166.5868, 168.0795, 169.6393, 171.1132, 172.644, 174.1977, 175.7153, 177.1939, 178.7749, 180.3209, 181.8713, 183.3725, 184.8155, 186.3248, 187.9331, 189.5146, 191.0683, 192.6219, 194.2148, 195.7379, 197.2527, 198.8556, 200.3384, 201.9805, 203.5676, 205.1453, 206.6924, 208.2987, 209.7792, 211.3952, 212.9141, 214.5333, 216.0497, 217.6791, 219.1909, 220.7126, 222.3645, 223.9554, 225.4199, 227.0668, 228.5208, 230.148, 231.7239, 233.3228, 234.796, 236.4266, 237.8862, 239.4911, 241.0799, 242.6068, 244.2092, 245.6739, 247.1945, 248.8089, 250.35, 251.8395, 253.333, 254.8083, 256.4745, 257.9547, 259.391, 261.001, 262.5595, 263.9024, 265.4909, 267.0333, 268.4612, 270.0316, 271.3819, 272.9637, 274.4076, 275.8462, 277.323, 278.7419, 280.1799, 281.8148, 283.1347, 284.7428, 286.0481, 287.5671, 289.0032, 290.4367, 291.826, 293.2114, 294.6785, 295.999, 297.3935, 298.7451, 300.1616, 301.495, 302.9146, 304.0946, 305.5698, 306.8712, 308.1873, 309.5497, 310.6798, 312.1485, 313.3851, 314.6745, 315.8906, 317.2426, 318.466, 319.7589, 321.1169, 322.3329, 323.6524, 324.8433, 326.0658, 327.3804, 328.3352, 330.0755, 331.0963, 332.2614, 333.4287, 334.5732, 335.8533, 336.9619, 338.1557, 339.0752, 340.4935, 341.4733, 342.6878, 343.6732, 344.9135, 345.9123, 347.0063, 348.1861, 349.1675, 350.2847, 351.3087, 352.3439, 353.349, 354.451, 355.4977, 356.521, 357.5755, 358.5093, 359.6444, 360.6372, 361.6318, 362.6254, 363.448, 364.6417, 365.445, 366.3199, 367.4148, 368.2568, 369.2261, 370.0213, 371.1453, 371.9515, 372.8166, 373.7613, 374.5588, 375.439, 376.2848, 377.1585, 378.0128, 378.7662, 379.6115, 380.5388, 381.3246, 382.0032, 382.9014, 383.6388, 384.5459, 385.2389, 386.0456, 386.8302, 387.4339, 388.4481, 389.0471, 389.7476, 390.4768, 391.3806, 391.886, 392.7808, 393.306, 394.1624, 394.6953, 395.5501, 396.0992, 396.7232, 397.3798, 398.1007, 398.8836, 399.3638, 400.2251, 400.7775, 401.5888, 401.7603, 402.7888, 403.0194, 403.846, 404.2552, 405.0897, 405.2602, 406.0968, 406.5248, 407.3216, 407.6328, 408.3774, 408.5997, 409.477, 409.4962, 410.4496, 410.6284, 411.5054, 411.8504, 412.5837, 412.179, 413.3202, 413.4938, 424.0788, 414.4241, 434.965, 416.0598, 435.6463, 417.2055, 436.9722, 418.1032, 438.0614, 419.1934, 439.2422, 419.9899, 420.676, 440.7668, 421.5294, 441.9936, 422.8579, 442.9793, 423.3674, 444.4373, 424.2965, 445.3946, 425.2526, 426.0937, 446.9731, 426.4589, 447.9497, 427.5232, 448.9064, 428.2756, 428.7141, 450.3486, 429.4565, 451.3502, 430.8269, 451.7631, 431.7076, 431.4621, 452.9458, 432.7558, 454.165, 433.6008, 433.4296, 455.2481, 434.5876, 456.1547, 435.1451, 435.779, 457.265, 436.3093, 458.0607, 436.6545, 437.5516, 459.2964, 437.7314, 460.4087, 460.2435, 438.7504, 461.4722, 439.3152, 439.5525, 462.5525, 440.2633, 463.1526, 463.6023, 441.1703, 464.2297, 442.0024, 442.0264, 465.2923, 442.9624, 442.6325, 466.4678, 443.7172, 443.3674, 467.0135, 444.5142, 467.4901, 467.6033, 445.2256, 468.2853, 468.411, 446.2087, 469.1551, 446.1194, 446.6802, 469.8448, 446.6895, 447.2688, 470.7585, 447.6326, 448.1518, 472.0661, 471.6331, 448.1199, 472.7137, 472.4605, 448.9248, 
1e-100, 0.3495867, 0.7053744, 1.067916, 1.437633, 1.813168, 2.195747, 2.584928, 2.979939, 3.382323, 3.790784, 4.205463, 4.62756, 5.055192, 5.48912, 5.930498, 6.377292, 6.830487, 7.290575, 7.756022, 8.228861, 8.707718, 9.192506, 9.683191, 10.1814, 10.68498, 11.19356, 11.71085, 12.23302, 12.76056, 13.29679, 13.83805, 14.38469, 14.93781, 15.49903, 16.06364, 16.63586, 17.21469, 17.80001, 18.38994, 18.98728, 19.59349, 20.2021, 20.81858, 21.44443, 22.07331, 22.70718, 23.3476, 24.00022, 24.65534, 25.3189, 25.98763, 26.66384, 27.34911, 28.03648, 28.73226, 29.43769, 30.14728, 30.86502, 31.59086, 32.32529, 33.06529, 33.80806, 34.56476, 35.33198, 36.10167, 36.87646, 37.65779, 38.45447, 39.2557, 40.06708, 40.88286, 41.70949, 42.54627, 43.38116, 44.23284, 45.09922, 45.96815, 46.83633, 47.7247, 48.61652, 49.52445, 50.44502, 51.35746, 52.27375, 53.22162, 54.17731, 55.13143, 56.09922, 57.07342, 58.06396, 59.07095, 60.07241, 61.08237, 62.10919, 63.14472, 64.18425, 65.24327, 66.30092, 67.37353, 68.48088, 69.5853, 70.68001, 71.79688, 72.91304, 74.05919, 75.22594, 76.37263, 77.55416, 78.72685, 79.91625, 81.13317, 82.33564, 83.57296, 84.807, 86.04035, 87.30324, 88.5728, 89.8825, 91.14539, 92.44253, 93.76059, 95.09113, 96.45662, 97.77566, 99.13687, 100.4915, 101.8916, 103.2826, 104.6901, 106.1342, 107.5381, 108.9543, 110.4031, 111.8882, 113.3452, 114.8205, 116.2916, 117.7931, 119.3199, 120.8589, 122.3627, 123.9397, 125.4864, 127.0484, 128.5972, 130.2273, 131.8454, 133.4243, 134.9817, 136.6861, 138.3476, 139.9689, 141.6209, 143.2812, 144.9576, 146.6567, 148.3702, 150.052, 151.7257, 153.5067, 155.1993, 156.9261, 158.6918, 160.4205, 162.2322, 164.0541, 165.7893, 167.563, 169.3545, 171.1787, 173.0409, 174.841, 176.6662, 178.4834, 180.3375, 182.1523, 183.9371, 185.8762, 187.6826, 189.6229, 191.5648, 193.4064, 195.2775, 197.1989, 199.0676, 201.0085, 202.8992, 204.808, 206.645, 208.6152, 210.5401, 212.4847, 214.466, 216.3961, 218.3568, 220.3049, 222.2591, 224.2188, 226.0356, 227.9453, 230.0141, 231.914, 233.9777, 235.9197, 237.8687, 239.8227, 241.7753, 243.7645, 245.7812, 247.7778, 249.7598, 251.7458, 253.696, 255.5714, 257.5653, 259.647, 261.6147, 263.5787, 265.6188, 267.516, 269.4969, 271.4456, 273.5278, 275.5322, 277.5228, 279.4164, 281.4008, 283.4181, 285.3332, 287.3444, 289.3151, 291.0929, 293.117, 295.1359, 296.9478, 299.0613, 300.7969, 302.7702, 304.6416, 306.7149, 308.3614, 310.3809, 312.2041, 314.1201, 315.936, 317.7666, 319.7163, 321.538, 323.4655, 325.4478, 327.307, 328.963, 331.0141, 332.7254, 334.7081, 336.3978, 338.2932, 340.0608, 342.0608, 343.7288, 345.4496, 347.2837, 349.1472, 350.8203, 352.6517, 354.293, 355.895, 357.8259, 359.4389, 361.0893, 362.9226, 364.6034, 366.2017, 367.7993, 369.6378, 371.2786, 372.8904, 374.6861, 376.2988, 377.8244, 379.506, 381.2302, 382.9154, 384.4437, 386.1533, 387.7032, 389.1848, 390.8264, 392.345, 393.719, 395.3091, 396.7893, 398.4458, 399.9485, 401.3912, 402.6006, 404.3873, 405.3711, 407.4966, 408.6659, 410.1335, 411.5494, 412.7279, 414.2976, 415.6096, 417.0806, 418.4876, 419.8965, 421.3201, 422.5767, 423.952, 425.3231, 426.4331, 427.8433, 429.0863, 430.61, 431.8008, 433.1621, 434.1552, 435.4986, 436.6511, 438.0251, 439.1467, 440.547, 441.5963, 442.786, 443.9637, 445.0971, 446.2054, 447.4407, 448.5564, 449.3941, 450.7518, 451.7183, 452.9855, 453.991, 454.9401, 456.1778, 457.1184, 458.2765, 459.41, 460.4414, 461.3344, 462.5649, 463.4414, 464.3945, 465.4027, 466.2304, 467.3445, 468.2961, 469.3676, 470.312, 471.1406, 471.9565, 472.8472, 473.7825, 474.7945, 475.7139, 476.3324, 477.4321, 478.1747, 478.9987, 479.9639, 480.6304, 481.6487, 482.3229, 482.8113, 483.9794, 484.6653, 485.3654, 486.0412, 487.0537, 487.7471, 488.6825, 489.2488, 490.1843, 490.6306, 491.3823, 492.0045, 493.0929, 493.4425, 494.2653, 494.8454, 495.7166, 495.9643, 497.2277, 497.4025, 498.1013, 498.6521, 499.62, 499.9017, 500.9221, 500.9167, 501.8821, 502.089, 502.9851, 503.8084, 504.1433, 503.9913, 505.1372, 505.4894, 518.0291, 506.6908, 530.7943, 508.1769, 531.8218, 509.2702, 533.2426, 510.7335, 534.633, 512.0206, 535.8328, 512.6289, 513.5494, 537.9096, 514.8827, 538.9433, 516.0959, 540.3787, 516.671, 541.7176, 518.007, 542.9834, 518.7877, 519.5085, 544.9071, 520.4466, 545.8567, 521.4442, 546.9238, 522.4469, 522.9132, 548.5683, 523.7667, 549.6345, 524.9878, 550.2598, 526.0186, 526.0725, 551.8052, 527.3613, 552.7914, 528.2111, 528.0993, 554.1773, 529.6475, 555.1368, 529.8715, 530.4968, 556.6312, 531.3913, 557.2331, 532.1278, 532.5691, 558.651, 532.923, 560.4158, 560.1175, 534.0034, 561.1204, 534.4284, 535.2206, 562.5384, 535.7089, 563.2423, 563.5226, 536.7281, 564.5512, 538.0137, 537.5095, 565.5194, 538.8016, 538.4181, 566.5812, 539.8054, 539.4021, 567.6773, 540.607, 567.8363, 568.4198, 541.4019, 568.9114, 569.3108, 542.0825, 569.8221, 542.3668, 542.7894, 570.8736, 543.0026, 543.7897, 571.6525, 543.5191, 544.4704, 573.4973, 573.0262, 544.6405, 573.6695, 573.528, 545.4877, 
1e-100, 0.3717751, 0.7513665, 1.138049, 1.53128, 1.932883, 2.340722, 2.756287, 3.179388, 3.608672, 4.045965, 4.490479, 4.941353, 5.400707, 5.866435, 6.339139, 6.821043, 7.309038, 7.803293, 8.306737, 8.816929, 9.333707, 9.859239, 10.39191, 10.93053, 11.4798, 12.03549, 12.59666, 13.16827, 13.74674, 14.33211, 14.92528, 15.52957, 16.13803, 16.75451, 17.38406, 18.01752, 18.66009, 19.30926, 19.9706, 20.63919, 21.31151, 22.00048, 22.6935, 23.39189, 24.10642, 24.83134, 25.5585, 26.29727, 27.04421, 27.80971, 28.5725, 29.34556, 30.134, 30.93441, 31.73456, 32.55158, 33.38296, 34.22227, 35.06695, 35.92704, 36.7989, 37.67588, 38.56331, 39.46691, 40.37751, 41.30793, 42.2354, 43.18053, 44.14961, 45.11845, 46.11016, 47.09595, 48.10046, 49.1204, 50.15934, 51.1993, 52.255, 53.32409, 54.4203, 55.5076, 56.60261, 57.72621, 58.8806, 60.03307, 61.20003, 62.35629, 63.53326, 64.74802, 65.95714, 67.19603, 68.42666, 69.70168, 70.96579, 72.2552, 73.54909, 74.86995, 76.21973, 77.54493, 78.90927, 80.29791, 81.67661, 83.06451, 84.49985, 85.91755, 87.37064, 88.82661, 90.29313, 91.81373, 93.30686, 94.83047, 96.33823, 97.91387, 99.5049, 101.0537, 102.6645, 104.2668, 105.8921, 107.537, 109.168, 110.8958, 112.547, 114.2739, 115.9457, 117.7211, 119.4896, 121.2084, 122.9888, 124.7697, 126.6033, 128.4331, 130.238, 132.0907, 133.9594, 135.7789, 137.6658, 139.6486, 141.5426, 143.4677, 145.4056, 147.3708, 149.3622, 151.3533, 153.3317, 155.3627, 157.4276, 159.4368, 161.448, 163.499, 165.6138, 167.764, 169.856, 171.8942, 174.0411, 176.1563, 178.3332, 180.433, 182.5896, 184.7692, 186.9647, 189.1845, 191.3812, 193.6296, 195.8261, 198.1083, 200.3895, 202.6596, 204.8603, 207.2013, 209.479, 211.8143, 214.1191, 216.3841, 218.6598, 220.9448, 223.315, 225.6415, 227.972, 230.3814, 232.7583, 235.0549, 237.4774, 239.895, 242.3036, 244.6646, 247.0974, 249.4585, 251.8514, 254.3424, 256.8267, 259.2267, 261.748, 264.2731, 266.6056, 269.1237, 271.5698, 273.9626, 276.4525, 278.9545, 281.2659, 283.5218, 286.0542, 288.4966, 290.9909, 293.3701, 295.7791, 298.1012, 300.6021, 303.0259, 305.38, 307.9846, 310.3244, 312.8402, 315.3868, 317.9246, 320.2932, 322.938, 325.2953, 327.7196, 330.131, 332.6173, 334.9799, 337.3872, 339.9212, 342.3208, 344.5074, 347.0894, 349.3606, 351.8573, 354.1754, 356.442, 358.8726, 361.2152, 363.5724, 366.0065, 368.3357, 370.696, 373.1559, 375.3475, 377.753, 380.3146, 382.5736, 384.7579, 387.2295, 389.4178, 391.8057, 393.8914, 396.3502, 398.5305, 401.0985, 403.1071, 405.2256, 407.5462, 409.8233, 411.8762, 414.2459, 416.3215, 418.469, 420.5696, 422.8489, 424.9237, 427.0035, 429.2583, 431.4774, 433.5914, 435.7135, 437.9094, 439.8121, 442.0548, 444.163, 445.9788, 448.2838, 450.1157, 452.2343, 454.2731, 456.1676, 458.1912, 460.3232, 462.1919, 464.1896, 466.1619, 467.8784, 469.8717, 471.6232, 473.7509, 475.4615, 477.1764, 479.1502, 481.0558, 482.663, 484.7813, 486.6553, 488.3047, 490.0785, 491.856, 493.568, 495.4873, 497.0612, 498.6866, 500.827, 502.2475, 504.1212, 505.5864, 507.5096, 508.9268, 510.6585, 512.1528, 514.0581, 515.4407, 516.937, 518.5407, 520.2337, 521.4709, 523.3309, 524.8539, 526.1548, 527.6758, 529.182, 530.5307, 531.9737, 533.5712, 534.8461, 536.4502, 537.7177, 539.4459, 540.7563, 541.9899, 543.5185, 544.6328, 546.0809, 547.4656, 548.8074, 549.7211, 551.3792, 552.3843, 553.7493, 555.0185, 556.4993, 557.6317, 558.7678, 560.0155, 561.3664, 562.301, 563.5793, 564.664, 565.5555, 566.9814, 568.0173, 569.0058, 570.4214, 571.1529, 572.4977, 573.4532, 574.601, 575.7025, 576.4937, 577.7177, 578.5065, 580.0035, 580.7831, 581.54, 582.6332, 583.6019, 584.3962, 585.5657, 586.5901, 587.4057, 588.2134, 589.2669, 589.9789, 590.9982, 591.9765, 592.7319, 593.2177, 594.6694, 595.1562, 596.2171, 597.1443, 597.7917, 598.3735, 599.3823, 600.1797, 601.0923, 601.2986, 602.5946, 602.776, 604.032, 604.3511, 605.4605, 605.9966, 606.7141, 607.2128, 607.994, 608.5592, 609.5436, 609.8376, 610.7202, 611.4617, 611.8143, 611.975, 613.2913, 613.4217, 628.4159, 614.9675, 644.3761, 616.3084, 645.8355, 618.0426, 646.9784, 619.1614, 648.4993, 621.0275, 650.185, 621.547, 622.5061, 652.2922, 624.2012, 653.7977, 625.0864, 655.0637, 626.1019, 656.5609, 627.353, 658.2628, 628.3995, 629.4215, 659.982, 630.0715, 661.1865, 631.2621, 662.7494, 632.3396, 632.8405, 664.4218, 633.9799, 665.4873, 635.3448, 666.6157, 636.2513, 636.2735, 668.0788, 638.2966, 669.408, 638.9335, 638.6348, 670.7852, 640.5105, 671.9399, 640.8949, 641.2321, 673.0929, 642.3648, 674.6651, 643.1549, 643.744, 675.6085, 644.3301, 677.1806, 677.1622, 645.3593, 678.5248, 646.2156, 646.4618, 679.7176, 647.4233, 680.888, 681.3052, 648.2448, 681.9819, 649.7179, 649.3744, 683.2609, 650.4853, 650.4129, 684.4152, 651.2772, 650.9062, 685.9793, 652.5519, 686.1806, 686.4751, 653.1772, 687.635, 687.8644, 654.1036, 688.207, 654.3784, 655.0374, 689.2918, 655.2213, 655.9737, 690.444, 655.9048, 656.6616, 691.9413, 691.5974, 656.9041, 692.5557, 692.7549, 657.9003, 
1e-100, 0.3991513, 0.805305, 1.220385, 1.643622, 2.073963, 2.513862, 2.96125, 3.416379, 3.881207, 4.353427, 4.834487, 5.324941, 5.822492, 6.329454, 6.845994, 7.369565, 7.903174, 8.44539, 8.99595, 9.557718, 10.12865, 10.70701, 11.296, 11.89565, 12.50305, 13.12017, 13.74941, 14.38748, 15.03445, 15.69554, 16.36303, 17.04153, 17.73277, 18.43687, 19.14686, 19.86737, 20.60691, 21.35208, 22.10695, 22.879, 23.66082, 24.45143, 25.25733, 26.07695, 26.90857, 27.74876, 28.60248, 29.47505, 30.35708, 31.25031, 32.16274, 33.09315, 34.0248, 34.96683, 35.93847, 36.9194, 37.91212, 38.91971, 39.94844, 40.98435, 42.04283, 43.10041, 44.20025, 45.30201, 46.41771, 47.54029, 48.69692, 49.87651, 51.04781, 52.24693, 53.46744, 54.69535, 55.95415, 57.20739, 58.50279, 59.82167, 61.13208, 62.46173, 63.82145, 65.20082, 66.59243, 68.0048, 69.43408, 70.85931, 72.34353, 73.84673, 75.32949, 76.85393, 78.38412, 79.98229, 81.54941, 83.16235, 84.7545, 86.39719, 88.04774, 89.74792, 91.4242, 93.13657, 94.84225, 96.67107, 98.44752, 100.2161, 102.0104, 103.8427, 105.716, 107.6023, 109.4825, 111.38, 113.2913, 115.2672, 117.267, 119.259, 121.2606, 123.2316, 125.3387, 127.3494, 129.4745, 131.5589, 133.6914, 135.8427, 137.9882, 140.1911, 142.4064, 144.6008, 146.841, 149.041, 151.3532, 153.6606, 155.9275, 158.1834, 160.5912, 162.9437, 165.3132, 167.7258, 170.1506, 172.5787, 175.029, 177.4766, 180.0095, 182.5127, 184.9797, 187.5292, 190.0692, 192.6089, 195.2267, 197.7619, 200.3644, 202.9947, 205.5849, 208.2904, 210.9008, 213.6591, 216.3064, 218.9702, 221.6177, 224.4616, 227.2598, 230.0191, 232.8927, 235.6166, 238.4749, 241.3145, 244.1172, 246.8616, 249.7415, 252.54, 255.3169, 258.1886, 261.1306, 263.9403, 266.6963, 269.5506, 272.4735, 275.1881, 278.1626, 281.0644, 284.0787, 286.9978, 289.9795, 292.9749, 296.08, 299.118, 302.014, 304.8626, 307.866, 310.8656, 313.7712, 316.8521, 319.7641, 322.8392, 325.686, 328.7762, 331.5382, 334.4752, 337.5009, 340.559, 343.4051, 346.5635, 349.5165, 352.3402, 355.4647, 358.4816, 361.5286, 364.5733, 367.6377, 370.6981, 373.6375, 376.5499, 379.6383, 382.4398, 385.5361, 388.5447, 391.3111, 394.4145, 397.3376, 400.3647, 403.112, 405.9316, 409.0592, 412.0214, 414.7986, 417.9686, 420.7479, 423.6388, 426.7025, 429.5434, 432.3822, 435.324, 438.3912, 441.1037, 443.8467, 446.7431, 449.7622, 452.5474, 455.4118, 458.053, 461.014, 463.6051, 466.5291, 469.2618, 471.9285, 474.7324, 477.5892, 479.9598, 482.8925, 485.6112, 488.0977, 490.8664, 493.9896, 496.6344, 498.9949, 501.6918, 504.3408, 507.0817, 509.5556, 512.1971, 514.6938, 517.4353, 520.0479, 522.2816, 525.0213, 527.4135, 529.9975, 532.642, 535.0839, 537.318, 539.7015, 542.339, 544.4915, 547.1469, 549.3009, 551.7349, 553.9327, 556.6367, 558.8007, 561.1534, 563.495, 565.8171, 567.8132, 570.2367, 572.5019, 574.6534, 576.9696, 579.3766, 581.3172, 583.6029, 585.611, 587.9773, 589.7564, 592.09, 593.9776, 596.2223, 598.3526, 600.3648, 602.141, 604.5603, 606.28, 607.869, 610.6386, 612.3968, 614.2183, 615.8833, 618.1076, 619.7772, 621.7569, 623.6806, 625.4122, 627.4979, 629.1802, 630.6622, 633.0062, 634.3827, 636.1928, 637.861, 639.734, 641.5875, 643.2495, 644.771, 646.5144, 647.9965, 649.8388, 651.4507, 653.1489, 654.5815, 656.174, 657.6716, 659.1872, 660.8119, 662.3443, 663.7111, 665.1269, 666.6064, 668.0205, 669.574, 671.1494, 672.1652, 673.7721, 675.4664, 676.8109, 678.0051, 679.5876, 680.6493, 682.1982, 683.2258, 684.512, 686.0505, 687.1813, 688.0325, 689.8179, 691.0359, 692.2846, 693.4785, 694.444, 695.6066, 696.7698, 698.1513, 699.1309, 700.1452, 701.5601, 702.4608, 703.1871, 705.0812, 705.673, 706.7002, 707.5825, 708.981, 709.9254, 710.6494, 711.7531, 712.606, 713.9288, 714.656, 715.5284, 716.8611, 717.5517, 718.6135, 719.4151, 720.55, 721.1235, 722.2658, 722.8574, 723.8363, 724.5663, 725.7807, 726.2859, 727.0577, 728.1301, 728.8052, 729.2007, 730.3999, 731.2279, 731.7271, 732.0107, 733.0077, 734.0861, 734.8286, 735.5108, 735.9177, 735.6218, 737.6075, 737.8359, 755.3334, 739.2414, 773.3266, 741.132, 774.9273, 742.8069, 776.6784, 744.2216, 778.5132, 745.6161, 780.176, 746.8474, 747.9302, 782.3365, 749.233, 784.2086, 750.6769, 785.3788, 751.5754, 787.4119, 752.9185, 788.5288, 754.5184, 755.4331, 791.0747, 756.0103, 792.33, 757.6378, 794.0681, 758.5566, 759.0775, 795.5173, 760.0537, 797.0749, 761.7285, 798.1372, 762.8501, 762.9616, 800.1279, 764.6016, 801.1724, 765.3285, 765.7857, 802.821, 767.162, 804.0083, 767.871, 768.2532, 805.6519, 769.3564, 807.1782, 770.3582, 770.7231, 808.4087, 771.5831, 809.7775, 809.9896, 772.5015, 811.5863, 773.8075, 773.6823, 812.7903, 774.7935, 813.6319, 814.7957, 776.2521, 815.1597, 776.9767, 776.9355, 816.8159, 778.2105, 777.9475, 817.9553, 778.9896, 778.9964, 819.5571, 780.5126, 819.7408, 819.9543, 781.2584, 820.9879, 821.713, 782.366, 822.0561, 782.4636, 782.6871, 823.1375, 783.5752, 784.1304, 824.4087, 784.5293, 785.0373, 825.7721, 825.8176, 785.172, 827.0248, 827.0955, 786.0826, 
1e-100, 0.4306009, 0.8723351, 1.322554, 1.782083, 2.252395, 2.731033, 3.221349, 3.72127, 4.229663, 4.751031, 5.281414, 5.821054, 6.374503, 6.936418, 7.509534, 8.096586, 8.692815, 9.299576, 9.922148, 10.55292, 11.19463, 11.85427, 12.52111, 13.20135, 13.8996, 14.60777, 15.32422, 16.06167, 16.81206, 17.57143, 18.34911, 19.13868, 19.94389, 20.76555, 21.60382, 22.45766, 23.31851, 24.20417, 25.10717, 26.01634, 26.94204, 27.89413, 28.86065, 29.83977, 30.84497, 31.86532, 32.90089, 33.9466, 35.02479, 36.12514, 37.22802, 38.35132, 39.5038, 40.67683, 41.85097, 43.05964, 44.29517, 45.54833, 46.81178, 48.10561, 49.42411, 50.76342, 52.09426, 53.47258, 54.88444, 56.31457, 57.73821, 59.19364, 60.69985, 62.2113, 63.7618, 65.30779, 66.88003, 68.47725, 70.10903, 71.75548, 73.45074, 75.14249, 76.86972, 78.60959, 80.36551, 82.18547, 84.02715, 85.85356, 87.72753, 89.6001, 91.50742, 93.50227, 95.44968, 97.4606, 99.43821, 101.5233, 103.5972, 105.682, 107.7778, 109.934, 112.0743, 114.316, 116.5051, 118.7829, 121.0084, 123.2929, 125.6031, 127.9382, 130.3215, 132.6697, 135.1029, 137.5342, 140.0275, 142.5108, 144.9824, 147.5631, 150.1002, 152.6751, 155.3217, 157.8771, 160.5711, 163.2588, 165.9761, 168.6542, 171.4262, 174.2013, 176.9472, 179.7858, 182.6166, 185.4408, 188.3036, 191.2132, 194.2546, 197.1793, 200.2274, 203.2492, 206.2727, 209.3032, 212.3721, 215.3616, 218.4565, 221.5832, 224.6022, 227.7179, 230.8182, 234.0381, 237.1311, 240.2426, 243.3648, 246.621, 250.0371, 253.2748, 256.4982, 259.8697, 263.2567, 266.5611, 269.996, 273.3413, 276.7516, 280.098, 283.4797, 286.7741, 290.2661, 293.6021, 297.0884, 300.5032, 303.9478, 307.3284, 310.7739, 314.3208, 317.7739, 321.3293, 324.9335, 328.4941, 332.111, 335.7239, 339.2163, 342.7121, 346.2634, 349.918, 353.5153, 356.9953, 360.4515, 364.0841, 367.5734, 371.1716, 374.783, 378.3793, 381.8314, 385.5962, 389.2049, 392.8194, 396.4569, 400.297, 403.859, 407.3365, 411.0499, 414.5681, 418.3558, 422.1036, 425.5442, 429.4708, 432.7117, 435.7252, 439.6027, 443.2857, 446.9818, 450.2846, 453.9358, 457.5737, 461.0917, 464.6098, 468.2259, 471.5566, 475.392, 478.8997, 482.5454, 486.136, 489.7663, 492.9989, 496.6645, 500.1366, 503.5896, 507.16, 510.6292, 514.0436, 517.2037, 520.9259, 524.144, 527.7239, 531.1285, 534.2448, 537.8376, 541.0745, 544.4901, 547.7513, 551.0321, 554.3662, 557.8556, 561.1741, 564.2905, 567.6404, 571.1229, 574.0377, 577.491, 580.6752, 583.9464, 587.2004, 590.1344, 593.4262, 596.5681, 599.6827, 602.9892, 605.8415, 608.9732, 611.9268, 614.9351, 618.0724, 621.1867, 623.8152, 626.6517, 630.079, 632.6037, 635.8322, 638.4879, 641.7104, 644.7074, 647.5044, 650.5228, 653.1412, 656.0844, 658.5797, 661.3867, 664.3796, 666.9446, 669.739, 672.516, 674.9355, 677.8562, 680.4201, 683.1742, 685.8013, 688.2539, 690.8386, 693.2856, 695.968, 698.4597, 701.0777, 703.1808, 705.6618, 708.2469, 710.5243, 713.1875, 715.3992, 717.9696, 720.1658, 722.4692, 724.9559, 727.3708, 729.3388, 732.2076, 733.5392, 736.4366, 738.5268, 740.7481, 742.9393, 744.8206, 746.9599, 749.2643, 751.3687, 753.4239, 755.4166, 757.3087, 759.4594, 761.582, 763.6116, 765.3893, 767.4497, 769.1862, 771.0546, 772.9647, 775.0124, 776.6465, 778.4813, 780.6443, 782.1613, 784.1358, 786.2099, 787.665, 789.4098, 790.9684, 793.097, 794.5511, 796.314, 797.6071, 799.299, 800.9243, 802.6132, 804.5503, 805.9528, 807.6321, 808.9959, 810.7167, 812.0216, 813.5407, 815.0379, 816.4639, 818.0952, 819.3707, 820.5741, 822.0648, 823.417, 824.6783, 825.907, 827.462, 829.3534, 830.2561, 831.5392, 832.7618, 834.3136, 835.5367, 836.6705, 837.9568, 838.7155, 840.307, 841.2496, 842.7266, 843.5315, 844.9718, 846.1018, 847.1418, 848.0923, 849.4459, 850.3888, 851.1369, 852.9133, 853.5992, 854.4228, 855.3757, 856.9707, 857.4326, 858.4864, 859.1558, 860.8856, 861.1994, 862.2067, 862.8684, 864.0414, 865.2544, 866.3714, 866.4906, 867.6061, 868.4037, 869.5306, 869.9683, 871.2076, 871.3323, 872.5728, 873.3679, 874.422, 875.3297, 875.6721, 875.6053, 876.7861, 877.9811, 898.4357, 878.8472, 919.4247, 881.3288, 921.3686, 882.6783, 922.9361, 884.7868, 925.1533, 886.0591, 926.9374, 887.4437, 888.9227, 929.6293, 889.8843, 930.9661, 891.3435, 932.8787, 892.5979, 934.801, 893.9948, 936.1071, 895.5627, 896.5347, 938.8895, 897.4227, 940.3016, 898.8786, 942.2196, 900.288, 900.6909, 943.924, 901.5982, 945.1878, 903.7286, 946.6673, 904.6754, 905.0431, 948.6328, 906.3076, 949.7715, 907.4158, 908.1336, 951.9122, 909.0372, 952.5899, 909.8582, 911.0182, 954.7809, 912.2661, 956.1853, 912.4186, 913.3634, 957.9714, 913.8002, 959.0835, 959.3674, 915.3243, 961.1392, 916.3586, 916.6014, 962.4747, 917.7054, 963.7122, 964.2341, 919.2889, 965.4802, 919.8483, 919.6939, 966.6203, 921.0252, 920.8739, 967.8909, 922.329, 922.3726, 969.6807, 923.5171, 970.1179, 970.5902, 924.919, 971.6446, 971.8047, 924.9252, 972.4459, 925.8677, 926.154, 974.1067, 927.1996, 927.3292, 975.4943, 927.9262, 928.1445, 976.4508, 976.7399, 928.6915, 977.8665, 978.2159, 929.622, 
1e-100, 0.4758805, 0.9613627, 1.461181, 1.97165, 2.492479, 3.029115, 3.575697, 4.135174, 4.70873, 5.293305, 5.894012, 6.508708, 7.1339, 7.776434, 8.433603, 9.102147, 9.790111, 10.49037, 11.20594, 11.94129, 12.6916, 13.45433, 14.23859, 15.04048, 15.85382, 16.68997, 17.5446, 18.41337, 19.3048, 20.21854, 21.14709, 22.09161, 23.06639, 24.05598, 25.0606, 26.09432, 27.14425, 28.21614, 29.30921, 30.43681, 31.58148, 32.73783, 33.92692, 35.1486, 36.37918, 37.62704, 38.91671, 40.22951, 41.55223, 42.91221, 44.30774, 45.7261, 47.1435, 48.60131, 50.11237, 51.62346, 53.16021, 54.73687, 56.34077, 57.96168, 59.61833, 61.29463, 63.02667, 64.76969, 66.53727, 68.32703, 70.17418, 72.04038, 73.91421, 75.82451, 77.79675, 79.79534, 81.77715, 83.81263, 85.8935, 88.00755, 90.12931, 92.28355, 94.47667, 96.69645, 98.94969, 101.221, 103.5414, 105.8816, 108.2888, 110.7116, 113.1504, 115.6223, 118.1169, 120.6733, 123.2252, 125.8092, 128.404, 131.0728, 133.7793, 136.5129, 139.2108, 141.9929, 144.7962, 147.7347, 150.5479, 153.5051, 156.4584, 159.4656, 162.4668, 165.504, 168.6093, 171.7272, 174.8288, 177.9996, 181.0527, 184.3097, 187.4657, 190.6668, 193.8759, 197.1639, 200.4775, 203.7125, 207.1071, 210.5764, 214.0391, 217.5954, 221.0167, 224.6434, 228.1087, 231.7201, 235.2689, 238.9621, 242.483, 246.2267, 249.774, 253.4242, 257.1663, 260.8307, 264.5308, 268.3314, 272.159, 275.9908, 279.9166, 283.881, 287.7277, 291.7458, 295.5609, 299.4451, 303.4834, 307.3204, 311.3947, 315.4067, 319.3418, 323.4762, 327.4419, 331.4477, 335.4722, 339.528, 343.6355, 347.7524, 351.9642, 356.2393, 360.4156, 364.5357, 368.6206, 372.7774, 377.054, 381.3503, 385.5011, 389.7816, 393.9374, 398.1532, 402.1891, 406.5255, 410.8008, 414.9815, 419.2329, 423.3592, 427.6528, 431.9464, 436.29, 440.4828, 444.8755, 449.3007, 453.6319, 457.7049, 462.0815, 466.4531, 470.7374, 474.9493, 479.4171, 483.4941, 487.729, 492.0679, 496.2716, 500.6221, 504.9364, 509.1424, 513.2306, 517.6435, 521.8731, 526.4005, 530.2274, 534.133, 538.7118, 542.9541, 547.428, 551.4461, 555.8148, 559.9647, 563.9251, 568.1353, 572.4502, 576.6533, 580.5225, 584.7997, 588.9815, 593.1149, 596.8778, 601.1146, 604.8857, 609.2765, 613.3916, 617.1765, 621.3862, 625.4203, 629.2194, 633.5709, 637.2351, 641.0282, 645.413, 648.9738, 653.055, 656.6291, 660.9275, 664.3477, 668.3252, 672.4177, 676.0391, 679.8706, 683.5181, 687.3667, 690.8817, 694.7356, 698.4889, 701.8618, 705.6622, 709.3273, 713.0767, 716.3599, 719.8502, 723.9139, 727.2151, 730.8597, 734.2694, 737.9687, 741.4203, 744.6657, 748.3424, 751.5262, 754.9071, 758.3437, 761.5931, 765.1652, 768.2222, 771.6624, 774.8502, 777.6124, 781.569, 784.1496, 787.6322, 790.7233, 793.7836, 796.6569, 799.8155, 802.9537, 806.1272, 808.9523, 812.3025, 815.1185, 817.9385, 821.0619, 823.7831, 826.8787, 829.5882, 832.5564, 835.2248, 838.3935, 840.9253, 843.3332, 846.7269, 848.9311, 851.5557, 854.3075, 857.4707, 859.7015, 862.3599, 864.7805, 867.541, 869.8926, 872.3392, 874.3923, 877.9209, 879.8725, 882.1743, 884.7615, 886.9337, 889.5819, 891.6942, 893.8376, 896.7127, 898.5497, 900.8061, 903.4885, 905.2811, 907.4819, 909.6161, 912.3476, 914.0256, 916.2785, 918.1471, 920.454, 922.3655, 924.4586, 926.4764, 928.315, 930.3664, 932.3591, 934.1661, 936.0141, 938.1968, 939.9637, 941.6091, 943.7212, 945.3793, 946.6064, 949.2875, 950.7232, 952.0562, 954.2127, 956.1662, 957.5751, 959.2308, 960.8062, 963.0081, 964.3442, 965.6751, 967.0853, 968.716, 970.517, 971.7706, 973.1264, 974.9179, 976.5415, 978.1749, 979.3749, 980.7447, 982.036, 983.6457, 984.6124, 986.4497, 988.0116, 988.5776, 989.8196, 991.7311, 992.7573, 993.8473, 995.0934, 996.9481, 997.6701, 998.612, 999.9196, 1001.371, 1002.848, 1003.847, 1004.437, 1005.542, 1007.225, 1008.009, 1009.161, 1010.322, 1011.22, 1012.177, 1013.617, 1014.358, 1015.509, 1016.438, 1017.428, 1018.214, 1019.966, 1020.22, 1021.03, 1021.679, 1022.998, 1023.515, 1024.773, 1025.412, 1026.588, 1026.925, 1027.67, 1029.047, 1029.61, 1030.029, 1031.258, 1031.391, 1054.958, 1033.334, 1079.384, 1035.243, 1081.198, 1036.966, 1083.31, 1039.359, 1085.472, 1040.668, 1086.977, 1042.36, 1043.184, 1090.321, 1045.253, 1091.662, 1046.18, 1093.382, 1048.102, 1095.806, 1049.297, 1097.525, 1051.214, 1051.688, 1099.891, 1052.61, 1102.252, 1054.741, 1103.153, 1055.354, 1055.99, 1105.651, 1057.572, 1107.205, 1059.793, 1108.336, 1060.95, 1061.001, 1110.451, 1062.161, 1111.997, 1063.634, 1064.411, 1114.381, 1065.45, 1115.113, 1066.194, 1066.849, 1117.558, 1068.432, 1118.59, 1069.365, 1069.953, 1120.393, 1070.217, 1121.946, 1122.204, 1072.058, 1123.487, 1072.779, 1073.181, 1125.612, 1074.559, 1127.227, 1127.426, 1075.538, 1128.376, 1076.546, 1076.512, 1129.691, 1077.891, 1078.151, 1131.717, 1079.859, 1079.414, 1132.985, 1079.97, 1133.585, 1134.093, 1081.601, 1135.332, 1135.649, 1082.206, 1136.193, 1082.597, 1083.688, 1137.988, 1083.72, 1084.965, 1139.589, 1084.975, 1085.34, 1140.32, 1140.671, 1086.291, 1141.887, 1141.587, 1086.758, 
1e-100, 0.5340527, 1.0862, 1.650431, 2.23084, 2.827406, 3.437409, 4.068464, 4.712968, 5.37212, 6.05467, 6.750809, 7.463925, 8.200444, 8.950952, 9.723884, 10.51747, 11.32695, 12.15844, 13.0142, 13.88798, 14.78349, 15.70494, 16.64397, 17.61165, 18.60299, 19.61657, 20.64498, 21.71592, 22.80344, 23.91222, 25.05374, 26.22109, 27.40866, 28.63379, 29.88599, 31.16077, 32.46136, 33.80014, 35.16608, 36.55931, 37.98686, 39.43962, 40.92273, 42.43239, 43.99101, 45.58196, 47.19079, 48.83522, 50.51831, 52.23972, 53.97367, 55.75486, 57.57956, 59.42606, 61.29717, 63.22368, 65.1819, 67.17526, 69.19256, 71.26242, 73.35856, 75.48692, 77.63298, 79.85964, 82.12326, 84.42955, 86.70605, 89.0666, 91.48035, 93.9256, 96.40781, 98.89435, 101.453, 104.0236, 106.6436, 109.3085, 111.9997, 114.7475, 117.5482, 120.3407, 123.2602, 126.1986, 129.1529, 132.1837, 135.1628, 138.2018, 141.2878, 144.3925, 147.5287, 150.6908, 153.8906, 157.1327, 160.3805, 163.695, 167.0214, 170.4814, 173.9831, 177.5423, 180.9948, 184.7693, 188.2564, 191.8735, 195.552, 199.1596, 202.9565, 206.6004, 210.3761, 214.1972, 217.916, 221.7637, 225.7251, 229.7238, 233.7453, 237.9106, 242.0103, 246.0633, 250.22, 254.2417, 258.4748, 262.6378, 266.783, 270.9643, 275.329, 279.506, 283.8614, 288.0707, 292.6118, 297.027, 301.5019, 305.9791, 310.5917, 315.0597, 319.4775, 324.0345, 328.6442, 333.3447, 337.8684, 342.4117, 347.135, 351.7509, 356.4528, 361.1665, 365.7972, 370.5298, 375.2169, 380.0631, 384.8663, 389.9395, 394.5494, 399.5197, 404.3777, 409.4262, 414.2258, 419.0427, 423.9154, 428.7533, 433.6853, 438.7144, 443.4953, 448.4201, 453.3419, 458.3632, 463.0389, 468.1479, 473.231, 478.3472, 483.4789, 488.4608, 493.3571, 498.4607, 503.476, 508.4522, 513.4052, 518.6151, 523.7272, 528.5123, 533.3147, 538.5088, 543.4568, 548.5304, 553.5959, 558.3032, 563.4736, 568.4593, 573.2604, 578.2229, 583.3888, 588.5133, 593.217, 598.4464, 603.3418, 608.473, 613.2131, 618.4668, 623.2245, 628.2275, 633.2982, 638.4692, 642.703, 647.3629, 652.3445, 656.9809, 662.0732, 666.7595, 671.9837, 676.2655, 681.1758, 686.0392, 690.8478, 695.6798, 700.2752, 705.0413, 709.7275, 714.4233, 719.2624, 724.0761, 728.7481, 732.9938, 738.0652, 742.3057, 747.0125, 751.4052, 756.1949, 760.6308, 765.0022, 769.3882, 773.8933, 778.4481, 782.4467, 787.2338, 791.3703, 795.4946, 800.1695, 804.3262, 808.8904, 813.1484, 817.1126, 821.7043, 825.4013, 829.8808, 834.0826, 838.2581, 842.4732, 846.1646, 850.4563, 854.453, 858.3683, 862.3599, 866.7468, 870.4712, 874.1468, 878.2931, 881.9822, 886.2991, 889.7184, 893.0062, 897.1885, 900.7569, 904.6496, 908.0442, 912.0066, 915.785, 919.3252, 922.92, 926.7572, 930.0291, 933.3866, 936.9617, 940.6216, 943.7592, 947.4714, 950.7125, 953.7293, 957.6143, 960.4426, 963.9825, 967.5972, 970.346, 973.4184, 976.7667, 980.1139, 983.2267, 986.1886, 989.1056, 992.2176, 995.2787, 998.3467, 1001.34, 1004.02, 1007.185, 1009.907, 1012.818, 1015.869, 1018.757, 1021.111, 1024.719, 1026.848, 1028.666, 1033.168, 1035.242, 1037.769, 1040.012, 1042.844, 1045.582, 1048.077, 1050.413, 1053.133, 1055.36, 1058.077, 1060.345, 1062.953, 1065.104, 1067.557, 1069.848, 1072.026, 1074.47, 1076.832, 1078.66, 1081.703, 1083.607, 1085.28, 1087.65, 1089.966, 1092.202, 1093.887, 1095.996, 1098.534, 1100.206, 1101.942, 1104.048, 1106.547, 1107.912, 1109.978, 1111.957, 1113.451, 1115.709, 1117.519, 1119.371, 1120.961, 1122.875, 1124.558, 1126.447, 1128.193, 1129.535, 1130.992, 1133.019, 1134.769, 1135.838, 1137.714, 1139.846, 1141.103, 1142.347, 1144.069, 1145.378, 1147.554, 1148.628, 1149.685, 1151.493, 1152.633, 1154.245, 1155.389, 1157.02, 1157.993, 1159.48, 1161.411, 1162.369, 1163.444, 1164.789, 1166.073, 1167.275, 1169.026, 1169.878, 1170.605, 1172.12, 1172.978, 1174.436, 1175.519, 1176.597, 1177.747, 1178.904, 1179.872, 1180.902, 1181.898, 1183.598, 1184.029, 1184.832, 1185.9, 1187.155, 1187.736, 1188.734, 1189.564, 1191.017, 1191.82, 1192.919, 1193.426, 1194.21, 1195.481, 1195.807, 1196.014, 1198.247, 1198.407, 1224.861, 1199.626, 1252.471, 1202.035, 1254.403, 1204.066, 1256.827, 1206.001, 1259.16, 1207.813, 1261.097, 1210.091, 1210.5, 1263.733, 1212.103, 1265.596, 1213.945, 1267.712, 1215.431, 1270.204, 1217.126, 1272.437, 1218.512, 1219.076, 1273.98, 1220.485, 1275.933, 1222.561, 1278.093, 1223.52, 1224.078, 1279.927, 1226.117, 1282.413, 1227.155, 1284.037, 1229.278, 1229.013, 1285.713, 1230.157, 1287.701, 1232.058, 1231.984, 1289.144, 1233.543, 1290.736, 1235.064, 1236.079, 1293.327, 1236.644, 1294.361, 1237.557, 1237.8, 1295.746, 1239.035, 1297.665, 1298.336, 1241.036, 1300.337, 1241.752, 1241.656, 1301.676, 1242.49, 1303.186, 1303.989, 1244.414, 1304.475, 1245.04, 1245.611, 1306.568, 1247.056, 1246.778, 1308.122, 1248.52, 1248.016, 1309.644, 1248.876, 1310.923, 1311.046, 1250.138, 1311.657, 1312.074, 1251.569, 1313.3, 1252.497, 1252.958, 1314.988, 1253.0, 1253.648, 1315.952, 1253.653, 1254.461, 1317.479, 1318.094, 1255.551, 1319.635, 1319.398, 1256.246, 
1e-100, 0.6204027, 1.257363, 1.918933, 2.596926, 3.294547, 4.016903, 4.75653, 5.520285, 6.306285, 7.112096, 7.949052, 8.805765, 9.681804, 10.59383, 11.52521, 12.47963, 13.47103, 14.48047, 15.51951, 16.59736, 17.69354, 18.81623, 19.98236, 21.17059, 22.3846, 23.64391, 24.93214, 26.24566, 27.60197, 28.99015, 30.4104, 31.86108, 33.36048, 34.88974, 36.44199, 38.05354, 39.70047, 41.36618, 43.08443, 44.85675, 46.63444, 48.45515, 50.34068, 52.24761, 54.19133, 56.17863, 58.24588, 60.32022, 62.42404, 64.60068, 66.82363, 69.07279, 71.33579, 73.6919, 76.06751, 78.4959, 80.93996, 83.47012, 86.04588, 88.65158, 91.30096, 94.03805, 96.83773, 99.66099, 102.5305, 105.4095, 108.345, 111.2652, 114.2725, 117.2681, 120.3901, 123.4801, 126.5979, 129.8458, 133.2268, 136.6247, 140.0175, 143.5222, 147.0009, 150.4818, 154.0348, 157.6461, 161.2707, 164.9569, 168.6406, 172.3747, 176.1256, 180.0372, 184.0413, 188.0136, 192.0665, 196.11, 200.2683, 204.3577, 208.5464, 212.7565, 216.9668, 221.2332, 225.4867, 229.9364, 234.2077, 238.6617, 243.1752, 247.7346, 252.4651, 257.0775, 261.7226, 266.3476, 271.2198, 275.7724, 280.6263, 285.4091, 290.1346, 295.1493, 299.9822, 304.8445, 309.7127, 314.7248, 319.7614, 324.8952, 330.0407, 335.4073, 340.5038, 345.65, 350.9756, 356.2153, 361.3992, 366.7138, 371.94, 377.1725, 382.5456, 387.8953, 393.1031, 398.6903, 403.9764, 409.4279, 415.1768, 420.6749, 426.2687, 431.8069, 437.4135, 443.2713, 448.6914, 454.2058, 459.8071, 465.566, 471.0966, 476.7843, 482.3116, 488.1275, 493.6803, 499.4393, 504.9239, 510.6419, 516.4475, 522.3184, 528.2643, 533.9005, 539.9508, 545.4517, 551.4895, 557.2386, 562.8809, 568.7862, 574.6437, 580.5052, 586.0335, 592.0719, 597.6302, 603.5133, 609.131, 614.8105, 620.2855, 626.1907, 632.2526, 637.754, 643.7453, 649.8075, 655.3998, 661.0313, 666.9341, 672.5668, 678.65, 684.004, 689.9457, 695.5701, 701.2394, 706.9282, 712.6596, 718.2339, 724.0503, 729.6945, 735.123, 740.5553, 746.4043, 751.8359, 757.5758, 763.3345, 768.7299, 773.4125, 779.4749, 784.9847, 790.5385, 795.9916, 801.5236, 806.871, 812.3027, 817.5429, 823.0804, 828.27, 833.7607, 839.2524, 844.2556, 849.4573, 854.8153, 859.9211, 865.2356, 870.1777, 875.3845, 880.5287, 885.7227, 890.9073, 896.0438, 900.8029, 905.9844, 910.9402, 915.5846, 921.1182, 925.6366, 930.5525, 935.2582, 940.3361, 945.0197, 949.9323, 954.6146, 959.2439, 963.9967, 968.6075, 973.1724, 977.9415, 982.314, 987.1584, 991.4575, 996.1568, 1000.44, 1004.626, 1009.514, 1013.595, 1018.475, 1022.376, 1026.854, 1031.477, 1035.334, 1039.557, 1043.844, 1047.847, 1052.001, 1056.053, 1060.344, 1064.062, 1068.49, 1072.038, 1076.057, 1080.395, 1084.001, 1087.669, 1091.721, 1095.521, 1098.83, 1102.811, 1106.614, 1110.311, 1113.76, 1117.7, 1121.261, 1124.745, 1128.338, 1131.814, 1135.176, 1138.668, 1142.199, 1145.508, 1149.039, 1152.393, 1155.064, 1159.269, 1161.928, 1165.272, 1168.366, 1171.846, 1174.592, 1177.869, 1180.859, 1184.318, 1187.25, 1189.718, 1192.543, 1195.19, 1199.526, 1201.598, 1204.411, 1207.041, 1210.089, 1212.934, 1215.7, 1218.491, 1220.939, 1223.207, 1226.611, 1228.662, 1231.508, 1233.946, 1236.763, 1238.93, 1241.397, 1243.951, 1246.794, 1249.125, 1251.069, 1253.329, 1255.981, 1258.252, 1260.295, 1262.36, 1264.636, 1266.96, 1269.592, 1271.371, 1273.507, 1275.473, 1277.263, 1280.039, 1282.001, 1283.466, 1285.804, 1288.04, 1289.576, 1291.348, 1293.164, 1295.966, 1297.622, 1298.931, 1300.581, 1302.466, 1304.398, 1305.983, 1307.889, 1309.313, 1311.846, 1313.354, 1314.831, 1316.072, 1317.53, 1319.4, 1321.162, 1323.187, 1324.402, 1325.754, 1326.818, 1328.492, 1329.914, 1331.375, 1332.851, 1334.362, 1335.769, 1336.948, 1338.254, 1339.812, 1341.651, 1342.539, 1343.332, 1344.76, 1345.964, 1347.455, 1348.506, 1349.595, 1351.245, 1352.906, 1353.582, 1354.625, 1355.173, 1356.974, 1357.301, 1359.227, 1360.406, 1361.192, 1361.761, 1362.936, 1364.023, 1365.001, 1365.678, 1367.064, 1367.699, 1368.693, 1369.856, 1370.777, 1372.656, 1372.885, 1372.592, 1373.514, 1374.319, 1404.251, 1376.507, 1434.942, 1378.937, 1437.505, 1381.41, 1439.632, 1382.855, 1441.261, 1384.799, 1443.659, 1386.731, 1387.605, 1446.886, 1389.172, 1448.734, 1391.041, 1451.077, 1392.51, 1453.127, 1394.532, 1455.473, 1396.169, 1396.532, 1458.447, 1398.479, 1459.693, 1399.583, 1461.23, 1401.251, 1401.898, 1464.685, 1403.827, 1466.093, 1405.046, 1467.632, 1406.358, 1406.466, 1469.307, 1408.257, 1471.703, 1410.233, 1410.426, 1473.886, 1411.585, 1474.998, 1412.538, 1413.578, 1477.595, 1414.624, 1478.871, 1415.478, 1415.696, 1480.787, 1417.578, 1482.987, 1483.069, 1418.289, 1484.147, 1419.76, 1419.999, 1487.011, 1421.57, 1487.669, 1487.923, 1422.091, 1489.456, 1423.921, 1423.875, 1492.169, 1425.583, 1425.772, 1493.442, 1426.472, 1425.982, 1494.315, 1427.525, 1495.028, 1496.137, 1429.005, 1497.408, 1497.914, 1429.766, 1498.63, 1430.768, 1431.048, 1500.603, 1431.889, 1432.068, 1501.456, 1432.336, 1432.918, 1503.67, 1503.875, 1433.665, 1504.644, 1504.608, 1434.705, 
1e-100, 0.7353807, 1.499148, 2.283771, 3.098986, 3.939091, 4.805194, 5.706005, 6.629942, 7.583488, 8.573387, 9.586187, 10.6351, 11.71836, 12.82764, 13.98157, 15.16653, 16.37816, 17.63717, 18.9299, 20.25627, 21.62237, 23.03181, 24.47159, 25.96191, 27.48972, 29.05124, 30.65251, 32.30729, 33.99888, 35.72872, 37.52988, 39.3606, 41.22617, 43.16278, 45.13473, 47.13881, 49.19642, 51.31095, 53.45565, 55.66571, 57.93637, 60.25092, 62.61104, 65.04575, 67.54697, 70.10576, 72.68516, 75.30898, 77.94964, 80.66372, 83.42781, 86.21319, 89.05897, 91.93716, 94.90233, 97.99041, 101.1364, 104.3572, 107.5636, 110.8311, 114.1589, 117.4997, 120.8969, 124.3148, 127.8269, 131.3394, 134.9386, 138.7315, 142.5242, 146.4081, 150.2854, 154.1771, 158.1789, 162.1813, 166.267, 170.3619, 174.4631, 178.7063, 182.9214, 187.1654, 191.5854, 196.1334, 200.6804, 205.2598, 209.9616, 214.568, 219.2186, 223.9455, 228.7241, 233.4724, 238.3625, 243.2274, 248.1682, 253.0775, 258.1426, 263.4202, 268.6994, 273.7798, 279.1181, 284.5375, 289.7561, 295.2013, 300.494, 306.0307, 311.3191, 316.9222, 322.4616, 327.9488, 333.6413, 339.1676, 345.0045, 350.9122, 356.6964, 362.4594, 368.3132, 374.3467, 380.3573, 386.205, 392.1276, 398.1214, 404.2277, 410.1285, 416.2419, 422.3082, 428.2736, 434.4517, 440.5756, 447.0057, 453.4292, 459.6935, 465.8972, 472.3826, 478.6323, 485.1204, 491.3532, 497.8531, 504.2946, 510.4694, 516.9592, 523.4976, 529.8726, 536.2403, 542.6489, 549.1149, 555.6977, 562.3795, 568.8159, 575.4388, 582.1204, 588.7963, 595.365, 601.7395, 608.6011, 614.9362, 621.7101, 628.3791, 634.7769, 641.5359, 648.1541, 654.7512, 660.9522, 667.6689, 674.2157, 680.7765, 687.4596, 694.3221, 700.7643, 707.3113, 714.0888, 720.7097, 727.1204, 733.8056, 740.4518, 746.6404, 753.4486, 759.7973, 766.5475, 772.847, 779.3797, 785.9135, 792.2731, 798.691, 805.199, 811.3454, 817.997, 824.4298, 831.0886, 837.2377, 843.4116, 850.1008, 856.592, 862.8229, 869.1526, 875.1843, 881.709, 887.9906, 894.1726, 900.8667, 906.5247, 911.7335, 918.3168, 924.6878, 930.3039, 936.7833, 942.4927, 948.6001, 954.3611, 960.5448, 966.6549, 972.131, 978.3093, 984.3716, 990.0425, 995.963, 1002.126, 1007.582, 1013.164, 1019.174, 1024.526, 1030.096, 1036.094, 1041.726, 1046.938, 1052.562, 1057.658, 1063.494, 1069.025, 1074.083, 1079.747, 1084.803, 1090.223, 1095.388, 1100.715, 1106.047, 1111.319, 1116.035, 1121.866, 1126.467, 1131.445, 1136.458, 1141.696, 1146.782, 1151.593, 1156.186, 1160.99, 1166.047, 1170.9, 1175.542, 1180.605, 1184.712, 1189.804, 1194.222, 1199.349, 1203.602, 1207.79, 1212.324, 1216.67, 1221.312, 1225.602, 1230.335, 1234.505, 1238.587, 1242.886, 1247.081, 1251.788, 1255.454, 1259.509, 1263.559, 1267.626, 1271.757, 1275.833, 1279.795, 1283.825, 1287.305, 1291.165, 1295.187, 1298.779, 1302.67, 1306.537, 1310.24, 1313.434, 1317.229, 1320.669, 1324.49, 1327.971, 1331.396, 1335.011, 1338.369, 1342.123, 1345.056, 1348.312, 1351.505, 1355.19, 1358.101, 1361.997, 1364.778, 1367.553, 1370.47, 1374.653, 1377.074, 1380.066, 1383.42, 1386.359, 1388.831, 1391.493, 1394.616, 1397.589, 1400.296, 1403.032, 1406.043, 1408.759, 1411.805, 1414.175, 1416.391, 1418.756, 1421.925, 1424.387, 1427.418, 1429.693, 1431.851, 1434.31, 1436.693, 1439.012, 1441.18, 1444.14, 1446.309, 1448.241, 1450.09, 1452.685, 1455.39, 1457.778, 1459.399, 1461.734, 1463.939, 1466.082, 1467.834, 1469.736, 1471.699, 1474.054, 1476.095, 1478.104, 1479.97, 1481.351, 1483.111, 1485.301, 1487.34, 1488.906, 1490.415, 1492.974, 1494.496, 1496.071, 1497.637, 1500.247, 1501.485, 1502.732, 1503.907, 1505.71, 1507.673, 1508.968, 1510.516, 1512.193, 1513.927, 1515.775, 1516.694, 1518.075, 1518.988, 1520.848, 1521.692, 1523.7, 1525.332, 1526.161, 1527.216, 1528.615, 1530.312, 1531.209, 1532.494, 1534.073, 1534.849, 1536.231, 1537.528, 1538.325, 1540.554, 1541.023, 1541.974, 1542.523, 1543.977, 1545.274, 1546.205, 1547.308, 1548.739, 1549.622, 1550.974, 1551.244, 1552.074, 1552.747, 1554.214, 1555.107, 1556.372, 1556.427, 1557.259, 1557.741, 1591.206, 1560.702, 1625.394, 1562.415, 1627.088, 1564.752, 1629.709, 1566.713, 1631.684, 1569.31, 1634.53, 1570.251, 1570.656, 1637.252, 1573.395, 1639.716, 1575.523, 1642.317, 1576.747, 1644.266, 1578.519, 1645.962, 1579.564, 1580.569, 1648.587, 1582.639, 1651.098, 1584.126, 1653.184, 1585.553, 1586.22, 1655.519, 1588.012, 1657.315, 1589.943, 1659.064, 1590.727, 1591.021, 1661.356, 1593.097, 1663.541, 1594.357, 1594.583, 1665.312, 1596.575, 1667.18, 1597.69, 1597.774, 1668.773, 1598.804, 1670.674, 1599.809, 1601.357, 1673.01, 1602.238, 1675.525, 1675.811, 1603.873, 1676.63, 1604.443, 1604.693, 1678.758, 1606.498, 1680.011, 1680.343, 1607.241, 1682.084, 1609.091, 1608.958, 1683.732, 1609.957, 1610.608, 1685.989, 1611.684, 1611.253, 1687.925, 1612.924, 1688.172, 1688.71, 1613.539, 1690.136, 1690.542, 1615.601, 1692.056, 1615.855, 1615.741, 1692.801, 1616.267, 1616.216, 1693.798, 1617.262, 1618.667, 1696.789, 1697.114, 1619.104, 1698.12, 1697.826, 1620.018, 
1e-100, 0.8936582, 1.818928, 2.784003, 3.776661, 4.808141, 5.876495, 6.974882, 8.119902, 9.298476, 10.51446, 11.77972, 13.07784, 14.41431, 15.80236, 17.22672, 18.69681, 20.22554, 21.78687, 23.40413, 25.0699, 26.77878, 28.52828, 30.34127, 32.19327, 34.10526, 36.08096, 38.10183, 40.19029, 42.35199, 44.57111, 46.83573, 49.12195, 51.4572, 53.85157, 56.30466, 58.78903, 61.3194, 63.92478, 66.69879, 69.4986, 72.38265, 75.2706, 78.21162, 81.21316, 84.27173, 87.37844, 90.53843, 93.72832, 97.07253, 100.4991, 104.021, 107.5576, 111.1526, 114.8327, 118.4992, 122.2702, 126.0542, 129.9124, 133.8708, 137.8383, 141.9134, 146.149, 150.3773, 154.7364, 159.1084, 163.5147, 167.9992, 172.4831, 177.1014, 181.7067, 186.3902, 191.0091, 195.8049, 200.8206, 205.8805, 210.929, 216.0231, 221.2316, 226.4247, 231.6932, 236.9302, 242.3058, 247.6738, 253.1456, 258.5687, 264.0647, 269.7342, 275.4622, 281.3162, 287.0494, 292.9508, 298.9569, 305.0035, 310.8882, 316.9566, 322.9445, 329.0807, 335.2895, 341.3396, 347.7669, 353.9499, 360.274, 366.759, 373.2415, 379.7116, 386.2865, 392.9068, 399.5684, 406.2776, 412.914, 419.5972, 426.3195, 432.8936, 439.6723, 446.4771, 453.2624, 460.173, 467.2653, 474.2401, 481.3175, 488.3351, 495.2412, 502.5363, 509.6944, 516.7448, 523.7928, 531.0032, 538.1849, 545.4132, 552.4725, 559.9418, 567.14, 573.9626, 581.4459, 588.7669, 596.0483, 603.4703, 610.985, 618.3972, 625.8528, 633.1413, 640.605, 647.9858, 655.4244, 662.9118, 669.9933, 677.627, 685.006, 692.2249, 700.0326, 707.0778, 714.6628, 722.1482, 729.6456, 736.7382, 744.5207, 751.9917, 759.5307, 766.9165, 774.4802, 781.8408, 789.382, 796.7682, 803.9884, 811.4309, 818.6986, 826.2394, 833.46, 840.9779, 848.0604, 855.4599, 862.5762, 870.1426, 877.3865, 884.5647, 891.9727, 899.2178, 906.4255, 913.5763, 921.1251, 927.8808, 935.2255, 942.3357, 949.6034, 956.5955, 963.8054, 970.7499, 977.7335, 984.6323, 991.7595, 998.6867, 1005.549, 1012.231, 1019.541, 1026.067, 1033.262, 1040.157, 1047.735, 1053.434, 1059.469, 1066.516, 1073.305, 1080.194, 1086.676, 1093.312, 1100.077, 1106.339, 1112.957, 1119.472, 1125.635, 1132.481, 1138.482, 1144.883, 1151.322, 1157.678, 1163.891, 1169.891, 1176.271, 1182.031, 1188.429, 1194.841, 1200.663, 1206.309, 1212.776, 1218.548, 1223.862, 1230.83, 1236.021, 1241.729, 1247.477, 1253.7, 1258.705, 1264.659, 1270.152, 1275.736, 1281.692, 1286.701, 1292.52, 1297.392, 1303.192, 1308.395, 1313.877, 1319.329, 1324.234, 1329.001, 1334.666, 1339.543, 1344.697, 1350.041, 1355.277, 1360.135, 1365.005, 1369.52, 1374.402, 1379.496, 1383.994, 1388.676, 1393.423, 1398.287, 1403.1, 1407.725, 1412.229, 1416.796, 1420.854, 1426.09, 1429.683, 1434.13, 1438.388, 1443.094, 1447.233, 1451.169, 1455.103, 1459.644, 1463.814, 1467.69, 1471.877, 1475.682, 1479.919, 1484.159, 1487.957, 1491.389, 1495.41, 1499.054, 1502.601, 1506.91, 1510.079, 1513.69, 1517.709, 1521.237, 1524.358, 1528.011, 1531.461, 1535.569, 1538.634, 1541.425, 1544.862, 1548.558, 1550.598, 1555.616, 1558.13, 1561.19, 1564.331, 1567.751, 1570.474, 1573.633, 1576.333, 1579.315, 1582.807, 1585.406, 1588.049, 1591.447, 1594.036, 1596.471, 1599.362, 1601.994, 1605.427, 1607.827, 1609.974, 1612.701, 1615.428, 1617.628, 1620.152, 1622.693, 1625.515, 1628.005, 1630.762, 1632.868, 1634.826, 1637.113, 1639.392, 1642.148, 1643.686, 1646.168, 1648.552, 1650.666, 1652.752, 1654.341, 1657.382, 1659.303, 1660.897, 1662.944, 1664.739, 1667.386, 1668.956, 1671.021, 1673.152, 1674.45, 1676.624, 1678.447, 1679.996, 1681.416, 1683.552, 1685.125, 1687.421, 1688.95, 1690.288, 1691.846, 1693.82, 1695.728, 1697.048, 1698.224, 1699.902, 1701.525, 1702.938, 1704.572, 1705.744, 1707.968, 1708.956, 1710.005, 1711.019, 1712.638, 1714.602, 1715.252, 1716.944, 1718.573, 1719.884, 1721.284, 1722.29, 1722.977, 1724.023, 1725.351, 1726.829, 1728.276, 1729.359, 1729.967, 1730.743, 1732.528, 1733.736, 1734.637, 1735.599, 1736.42, 1738.047, 1738.411, 1739.544, 1740.724, 1742.439, 1742.875, 1743.12, 1743.427, 1744.794, 1745.892, 1781.277, 1748.311, 1818.417, 1750.702, 1820.739, 1752.144, 1822.397, 1753.603, 1824.642, 1756.219, 1827.567, 1758.115, 1758.881, 1830.693, 1760.683, 1833.344, 1763.153, 1834.867, 1764.857, 1837.985, 1766.632, 1839.747, 1768.176, 1769.48, 1842.663, 1770.433, 1844.128, 1772.353, 1846.814, 1773.922, 1775.009, 1849.631, 1775.748, 1850.947, 1777.366, 1852.077, 1778.134, 1779.167, 1855.545, 1781.82, 1857.871, 1782.918, 1783.226, 1859.716, 1784.487, 1861.058, 1786.151, 1786.61, 1863.602, 1787.215, 1864.664, 1787.963, 1789.643, 1867.508, 1790.239, 1868.52, 1869.017, 1791.742, 1871.01, 1793.961, 1793.884, 1873.327, 1794.051, 1874.087, 1875.072, 1795.9, 1876.754, 1798.344, 1797.925, 1878.78, 1798.906, 1798.542, 1879.954, 1799.163, 1799.672, 1881.439, 1800.781, 1883.002, 1884.128, 1802.305, 1884.578, 1885.463, 1803.957, 1886.699, 1804.752, 1804.854, 1888.222, 1805.049, 1805.176, 1889.03, 1806.04, 1806.891, 1891.141, 1891.215, 1806.825, 1892.549, 1892.27, 1808.384, 
1e-100, 1.09971, 2.241027, 3.422254, 4.660003, 5.935435, 7.259862, 8.634342, 10.04772, 11.51413, 13.03011, 14.59353, 16.22275, 17.89842, 19.63638, 21.45737, 23.32311, 25.23336, 27.18961, 29.19069, 31.246, 33.3629, 35.52147, 37.75021, 40.14152, 42.56644, 45.03679, 47.56067, 50.13146, 52.78266, 55.47148, 58.2352, 61.04322, 63.96359, 67.02397, 70.16531, 73.31613, 76.52383, 79.81418, 83.13504, 86.55359, 89.9789, 93.5253, 97.1244, 100.9119, 104.7268, 108.6531, 112.6636, 116.6686, 120.7301, 124.8854, 129.0654, 133.3945, 137.6921, 142.0808, 146.632, 151.3279, 156.083, 160.8679, 165.6602, 170.5905, 175.5828, 180.6066, 185.641, 190.7619, 195.9514, 201.1931, 206.604, 212.1763, 217.7644, 223.4768, 229.174, 234.852, 240.6875, 246.4661, 252.4531, 258.4074, 264.4445, 270.4934, 276.6088, 282.8041, 289.2181, 295.6526, 302.1983, 308.6428, 315.2057, 321.8079, 328.5103, 335.204, 341.9497, 348.7146, 355.6023, 362.456, 369.398, 376.2766, 383.5138, 390.641, 398.0093, 405.1727, 412.1662, 419.9255, 427.1873, 434.5029, 441.9229, 449.526, 456.8085, 464.3495, 471.9073, 479.4366, 487.1289, 494.8437, 502.5959, 510.2688, 518.3534, 526.1555, 533.99, 541.8737, 549.7028, 557.7044, 565.7977, 573.7227, 581.7513, 589.6288, 597.4777, 605.686, 613.6019, 621.8382, 629.9865, 638.0345, 646.3746, 654.6918, 662.6987, 671.1708, 679.307, 687.5386, 695.7608, 703.9515, 712.2284, 720.7544, 728.8472, 737.212, 745.2706, 753.645, 761.5208, 770.0042, 778.485, 786.6141, 795.2221, 803.3866, 811.3194, 820.151, 828.2745, 836.7417, 844.8288, 853.5275, 861.531, 869.8156, 878.0978, 886.2247, 894.448, 902.3419, 910.9852, 919.1253, 926.9875, 935.1502, 943.421, 951.6355, 959.7009, 967.8886, 975.7389, 984.1076, 992.3006, 1000.43, 1007.951, 1016.321, 1024.26, 1032.124, 1039.858, 1047.759, 1055.608, 1063.487, 1071.349, 1079.096, 1086.821, 1094.734, 1102.269, 1110.178, 1117.829, 1125.305, 1133.284, 1141.012, 1148.221, 1155.923, 1163.164, 1170.998, 1178.024, 1185.836, 1193.455, 1201.001, 1207.288, 1213.948, 1221.601, 1228.719, 1236.317, 1243.18, 1250.349, 1257.053, 1264.44, 1271.615, 1278.339, 1285.316, 1292.184, 1299.053, 1305.866, 1312.761, 1319.519, 1326.021, 1333.157, 1339.013, 1345.837, 1352.645, 1359.124, 1365.067, 1371.573, 1377.969, 1384.502, 1390.723, 1396.803, 1403.177, 1409.219, 1415.419, 1421.615, 1427.115, 1433.492, 1439.458, 1444.988, 1451.444, 1456.921, 1462.941, 1468.15, 1474.729, 1479.893, 1485.581, 1490.851, 1496.463, 1502.354, 1507.384, 1513.066, 1518.172, 1524.214, 1529.033, 1534.468, 1539.761, 1544.671, 1549.43, 1554.883, 1559.341, 1564.437, 1569.331, 1575.185, 1579.635, 1584.313, 1588.592, 1593.774, 1598.636, 1603.289, 1607.617, 1612.121, 1616.781, 1621.267, 1625.956, 1630.188, 1634.659, 1638.686, 1643.161, 1647.696, 1651.706, 1655.644, 1660.155, 1663.947, 1667.863, 1671.96, 1676.222, 1680.061, 1684.358, 1688.082, 1691.469, 1695.536, 1699.321, 1702.784, 1706.311, 1710.187, 1713.494, 1717.602, 1721.444, 1724.342, 1727.46, 1731.282, 1733.432, 1738.7, 1741.325, 1745.138, 1747.977, 1750.997, 1753.89, 1757.276, 1760.3, 1763.818, 1767.029, 1769.61, 1772.768, 1775.787, 1778.505, 1780.92, 1784.125, 1786.395, 1789.921, 1792.846, 1795.219, 1797.519, 1800.352, 1802.934, 1805.742, 1808.219, 1811.289, 1813.243, 1815.345, 1818.197, 1820.468, 1823.107, 1825.818, 1827.816, 1829.431, 1832.127, 1834.314, 1836.314, 1838.349, 1840.661, 1842.793, 1845.734, 1847.49, 1849.91, 1851.515, 1853.054, 1855.93, 1857.563, 1859.116, 1861.453, 1863.244, 1865.185, 1866.756, 1868.277, 1871.05, 1872.389, 1873.711, 1875.438, 1877.163, 1879.134, 1881.178, 1882.827, 1884.353, 1886.123, 1887.565, 1888.822, 1889.874, 1891.377, 1892.909, 1894.484, 1896.39, 1897.625, 1898.664, 1899.915, 1901.81, 1903.051, 1904.499, 1905.584, 1907.046, 1908.56, 1909.69, 1910.873, 1912.84, 1913.704, 1914.697, 1915.152, 1916.617, 1917.803, 1919.149, 1920.046, 1922.001, 1922.509, 1923.994, 1924.709, 1925.608, 1926.055, 1927.445, 1928.35, 1929.835, 1930.568, 1931.038, 1931.606, 1932.486, 1933.409, 1971.89, 1935.681, 2011.469, 1938.191, 2014.381, 1940.508, 2016.382, 1942.798, 2019.198, 1944.488, 2021.005, 1946.226, 1947.649, 2024.974, 1949.757, 2027.631, 1951.343, 2029.615, 1952.797, 2031.387, 1953.965, 2033.208, 1955.919, 1956.627, 2037.259, 1959.178, 2038.988, 1960.228, 2040.904, 1962.551, 1962.984, 2043.918, 1964.637, 2045.897, 1965.797, 2047.009, 1967.129, 1967.892, 2050.224, 1969.342, 2051.511, 1970.483, 1970.94, 2053.878, 1972.65, 2055.384, 1973.303, 1974.453, 2058.059, 1975.214, 2059.981, 1977.09, 1977.73, 2063.22, 1979.482, 2064.325, 2064.282, 1979.853, 2066.053, 1981.472, 1982.209, 2068.317, 1982.394, 2069.284, 2069.962, 1984.575, 2071.705, 1985.623, 1985.439, 2074.267, 1987.297, 1987.054, 2075.537, 1988.366, 1988.638, 2077.143, 1989.3, 2078.21, 2079.024, 1990.764, 2080.039, 2080.891, 1991.695, 2081.217, 1992.014, 1992.148, 2082.215, 1992.67, 1994.02, 2084.613, 1994.521, 1995.268, 2087.129, 2086.832, 1995.6, 2088.372, 2087.747, 1996.86, 
1e-100, 1.365041, 2.804756, 4.303239, 5.846166, 7.434175, 9.068491, 10.75856, 12.50097, 14.29235, 16.19827, 18.19161, 20.23358, 22.32685, 24.47449, 26.68262, 28.95505, 31.27725, 33.6653, 36.22724, 38.8565, 41.52935, 44.27539, 47.07544, 49.94744, 52.88407, 55.87678, 58.9511, 62.17397, 65.49735, 68.91849, 72.36615, 75.91564, 79.48652, 83.16114, 86.87856, 90.7193, 94.59355, 98.6325, 102.8019, 107.0496, 111.3241, 115.7095, 120.1462, 124.6679, 129.2467, 133.9376, 138.6158, 143.4439, 148.4575, 153.545, 158.7207, 163.9558, 169.2422, 174.6099, 180.0498, 185.5659, 191.0983, 196.758, 202.5269, 208.3107, 214.2722, 220.3737, 226.4643, 232.7295, 238.9762, 245.2833, 251.6785, 258.1313, 264.7339, 271.1928, 277.8257, 284.5207, 291.3337, 298.2773, 305.3861, 312.4684, 319.597, 326.9167, 334.044, 341.4167, 348.5934, 356.1178, 363.5464, 371.0446, 378.6082, 386.1568, 393.8872, 401.8431, 409.7005, 417.5858, 425.5884, 433.6317, 441.5454, 449.7442, 457.8477, 465.9337, 474.2992, 482.4748, 490.4859, 499.1479, 507.538, 515.8483, 524.3038, 533.0052, 541.4765, 550.0704, 558.8265, 567.5237, 576.119, 584.8441, 593.7159, 602.3693, 611.2086, 619.8642, 628.6164, 637.6408, 646.5125, 655.4285, 664.2522, 673.4565, 682.5174, 691.4957, 700.5123, 709.5972, 718.6814, 727.5745, 737.0052, 745.9912, 755.0453, 764.1665, 773.0786, 782.2944, 791.3833, 800.5724, 809.5318, 818.7901, 827.9703, 837.1579, 846.3485, 855.3929, 864.8339, 873.9201, 883.1148, 892.141, 901.0707, 910.3598, 919.401, 928.6597, 937.6002, 946.8168, 955.9429, 964.9393, 973.7655, 982.8405, 992.1125, 1001.272, 1010.325, 1019.194, 1028.282, 1036.934, 1045.827, 1054.871, 1063.694, 1072.713, 1081.471, 1090.777, 1098.919, 1107.712, 1116.661, 1125.402, 1133.967, 1142.889, 1151.379, 1159.657, 1168.613, 1177.306, 1185.619, 1194.323, 1202.73, 1211.136, 1219.751, 1227.967, 1236.511, 1244.455, 1253.099, 1261.236, 1269.286, 1277.801, 1285.767, 1293.944, 1301.996, 1310.014, 1318.227, 1325.761, 1333.882, 1342.162, 1350.057, 1358.056, 1365.147, 1371.768, 1380.377, 1388.4, 1395.671, 1403.462, 1410.802, 1417.978, 1425.928, 1432.938, 1440.6, 1447.932, 1454.937, 1462.389, 1469.689, 1476.774, 1483.82, 1490.697, 1497.87, 1504.658, 1511.898, 1518.522, 1525.394, 1532.255, 1539.078, 1545.94, 1552.025, 1559.262, 1565.574, 1571.923, 1578.371, 1585.192, 1591.157, 1597.147, 1603.725, 1610.427, 1616.553, 1622.453, 1628.735, 1634.453, 1640.484, 1646.969, 1652.484, 1658.54, 1664.168, 1669.544, 1675.823, 1681.282, 1686.906, 1692.328, 1699.092, 1703.819, 1709.1, 1714.241, 1719.736, 1725.499, 1730.331, 1735.501, 1740.737, 1746.047, 1751.254, 1756.493, 1761.26, 1766.092, 1770.805, 1776.298, 1781.124, 1785.459, 1790.245, 1795.32, 1799.709, 1804.02, 1808.448, 1813.359, 1817.775, 1822.265, 1826.792, 1831.041, 1835.395, 1840.002, 1843.953, 1847.878, 1852.114, 1856.267, 1860.086, 1864.711, 1867.953, 1872.342, 1876.42, 1879.918, 1883.603, 1887.146, 1891.773, 1895.203, 1898.617, 1901.956, 1905.395, 1909.215, 1912.201, 1915.762, 1920.409, 1923.49, 1927.024, 1930.269, 1933.355, 1936.491, 1939.372, 1942.689, 1946.309, 1948.823, 1951.928, 1955.697, 1958.422, 1961.179, 1963.741, 1967.596, 1969.998, 1972.771, 1975.188, 1977.981, 1980.888, 1983.661, 1986.529, 1988.674, 1991.844, 1994.772, 1997.208, 1999.415, 2001.801, 2004.299, 2006.624, 2009.347, 2011.878, 2013.675, 2016.388, 2018.068, 2020.378, 2022.619, 2025.106, 2027.29, 2029.178, 2031.187, 2033.62, 2035.687, 2038.114, 2040.043, 2041.717, 2043.837, 2045.65, 2046.889, 2049.094, 2050.839, 2052.913, 2055.243, 2056.684, 2058.619, 2060.059, 2061.648, 2063.466, 2065.561, 2066.838, 2068.395, 2070.341, 2072.047, 2073.551, 2074.73, 2077.073, 2078.2, 2079.226, 2080.318, 2082.193, 2083.825, 2084.976, 2086.725, 2088.468, 2089.539, 2091.122, 2092.165, 2093.127, 2094.027, 2095.491, 2096.774, 2098.425, 2099.495, 2100.143, 2101.509, 2102.599, 2104.321, 2105.109, 2106.064, 2106.969, 2108.757, 2109.536, 2110.582, 2111.073, 2113.002, 2113.583, 2114.467, 2114.545, 2116.103, 2117.095, 2117.871, 2119.081, 2160.144, 2120.831, 2201.514, 2123.12, 2203.236, 2124.561, 2205.606, 2126.863, 2208.463, 2129.468, 2210.7, 2131.536, 2131.994, 2215.266, 2134.007, 2217.489, 2136.546, 2220.47, 2138.054, 2221.447, 2139.357, 2223.715, 2141.432, 2142.177, 2227.367, 2143.896, 2229.274, 2145.609, 2231.086, 2147.092, 2147.387, 2233.876, 2148.803, 2235.26, 2150.185, 2236.933, 2152.607, 2152.897, 2241.082, 2154.609, 2241.812, 2155.894, 2156.121, 2244.639, 2157.853, 2246.077, 2158.824, 2159.151, 2248.51, 2160.52, 2249.848, 2162.098, 2162.413, 2252.72, 2163.604, 2254.118, 2255.031, 2165.42, 2256.52, 2166.479, 2166.357, 2258.397, 2167.416, 2259.878, 2261.502, 2170.104, 2262.435, 2170.952, 2170.95, 2264.551, 2171.673, 2171.393, 2266.113, 2172.524, 2172.909, 2268.664, 2174.678, 2268.478, 2270.118, 2176.089, 2270.546, 2271.942, 2177.016, 2272.225, 2177.434, 2177.522, 2273.106, 2178.009, 2178.837, 2275.061, 2179.07, 2178.986, 2276.421, 2277.718, 2180.193, 2279.452, 2279.985, 2182.127, 
1e-100, 1.698023, 3.449211, 5.266015, 7.141338, 9.076107, 11.1519, 13.30438, 15.513, 17.78209, 20.11207, 22.52161, 24.99117, 27.5358, 30.21291, 32.9976, 35.86627, 38.79211, 41.79944, 44.87371, 48.02403, 51.25385, 54.55933, 58.0106, 61.60973, 65.28528, 69.01961, 72.82161, 76.71739, 80.69984, 84.75135, 88.91143, 93.12622, 97.52665, 102.0783, 106.682, 111.38, 116.1441, 121.0002, 125.9081, 130.9344, 136.0367, 141.195, 146.5043, 152.0095, 157.6033, 163.2603, 169.0307, 174.7958, 180.6656, 186.615, 192.6505, 198.8105, 204.9632, 211.2567, 217.8121, 224.479, 231.0888, 237.8368, 244.5655, 251.4688, 258.4355, 265.4637, 272.5324, 279.7143, 286.8927, 294.2216, 301.7371, 309.3537, 317.0555, 324.7604, 332.5627, 340.3643, 348.3135, 356.1616, 364.3045, 372.3686, 380.5747, 388.5973, 396.9718, 405.426, 413.8038, 422.467, 431.1982, 439.6767, 448.4341, 457.1651, 465.9865, 474.8327, 483.6493, 492.7374, 501.571, 510.5703, 519.7638, 528.8297, 538.0117, 547.3027, 556.557, 565.8769, 575.0757, 585.1303, 594.4095, 603.8807, 613.4048, 622.7594, 632.4697, 641.9747, 651.6029, 661.0351, 670.8012, 680.5838, 690.2411, 700.2539, 710.0538, 719.8465, 729.7216, 739.5675, 749.3625, 759.2591, 769.157, 779.0026, 789.0167, 799.008, 808.9728, 818.766, 828.6734, 838.4402, 848.6076, 858.8991, 868.9411, 878.766, 888.8578, 898.6357, 908.5912, 918.8564, 928.6576, 938.8276, 948.725, 958.5194, 968.518, 978.4686, 988.6981, 998.545, 1008.301, 1018.102, 1028.24, 1038.181, 1047.994, 1057.661, 1067.796, 1077.795, 1087.65, 1097.372, 1107.157, 1116.646, 1126.386, 1136.147, 1145.688, 1155.512, 1165.116, 1174.915, 1184.214, 1193.474, 1203.15, 1212.887, 1222.214, 1232.352, 1241.175, 1250.734, 1260.013, 1269.534, 1278.641, 1287.859, 1297.42, 1306.258, 1315.637, 1324.346, 1333.516, 1342.518, 1351.748, 1361.052, 1369.441, 1378.398, 1387.351, 1396.027, 1405.001, 1413.645, 1422.017, 1431.372, 1439.337, 1448.148, 1456.746, 1465.212, 1473.659, 1482.007, 1490.419, 1498.893, 1507.136, 1515.915, 1523.346, 1530.526, 1538.675, 1547.537, 1555.311, 1563.187, 1571.032, 1579.268, 1587.053, 1594.536, 1602.525, 1609.98, 1617.431, 1625.23, 1632.945, 1640.437, 1647.712, 1655.174, 1662.331, 1669.506, 1677.135, 1683.912, 1690.983, 1697.958, 1705.407, 1712.605, 1718.983, 1725.855, 1732.784, 1739.511, 1746.412, 1752.839, 1759.198, 1765.613, 1772.456, 1778.408, 1785.942, 1791.723, 1798.164, 1804.353, 1810.637, 1816.626, 1822.537, 1828.87, 1834.863, 1841.072, 1846.506, 1852.807, 1858.233, 1864.207, 1870.851, 1875.859, 1881.435, 1886.787, 1891.919, 1897.978, 1903.149, 1908.358, 1913.749, 1919.665, 1924.595, 1929.602, 1934.628, 1939.892, 1945.237, 1949.805, 1955.336, 1959.521, 1964.613, 1969.85, 1974.483, 1978.979, 1983.77, 1988.011, 1993.216, 1997.464, 2001.73, 2006.451, 2010.704, 2014.985, 2019.217, 2023.483, 2027.996, 2032.534, 2036.309, 2040.39, 2043.987, 2048.379, 2052.186, 2055.829, 2059.618, 2063.987, 2067.412, 2071.675, 2075.451, 2078.72, 2082.049, 2085.945, 2089.18, 2092.672, 2097.864, 2100.833, 2103.736, 2107.055, 2110.305, 2114.045, 2117.169, 2120.527, 2123.573, 2126.326, 2129.897, 2132.638, 2135.669, 2138.259, 2141.404, 2144.254, 2147.585, 2150.573, 2152.915, 2155.544, 2158.968, 2161.776, 2164.418, 2167.075, 2169.625, 2171.811, 2174.232, 2177.149, 2179.942, 2182.701, 2184.777, 2187.141, 2189.406, 2191.58, 2193.677, 2195.858, 2198.62, 2201.013, 2203.769, 2205.721, 2207.928, 2209.798, 2211.947, 2213.856, 2216.376, 2218.179, 2220.284, 2222.348, 2224.094, 2225.678, 2227.558, 2230.134, 2231.721, 2233.167, 2234.959, 2236.645, 2239.133, 2240.518, 2242.497, 2244.355, 2245.436, 2247.477, 2248.743, 2249.974, 2251.333, 2253.043, 2254.712, 2256.588, 2257.816, 2259.256, 2260.844, 2261.879, 2263.781, 2265.209, 2266.038, 2267.484, 2269.221, 2270.603, 2271.631, 2272.622, 2274.557, 2275.517, 2276.192, 2277.697, 2278.749, 2280.327, 2281.333, 2282.958, 2284.057, 2284.923, 2285.988, 2286.818, 2287.434, 2288.307, 2289.668, 2290.407, 2291.862, 2292.809, 2293.214, 2294.461, 2295.347, 2296.314, 2297.057, 2339.515, 2299.301, 2384.445, 2301.499, 2386.714, 2303.735, 2389.582, 2305.772, 2391.645, 2307.478, 2394.304, 2309.96, 2311.101, 2398.035, 2312.717, 2399.637, 2314.435, 2401.818, 2315.462, 2403.418, 2317.192, 2406.802, 2320.027, 2320.785, 2410.09, 2321.977, 2411.995, 2323.8, 2414.094, 2325.669, 2326.134, 2416.678, 2327.373, 2418.203, 2328.787, 2420.2, 2330.928, 2330.806, 2422.368, 2332.177, 2424.508, 2333.581, 2334.105, 2427.101, 2336.293, 2429.155, 2337.12, 2337.564, 2431.784, 2339.174, 2433.698, 2340.896, 2341.408, 2435.919, 2342.126, 2437.265, 2437.81, 2343.387, 2439.81, 2344.581, 2344.493, 2441.401, 2345.403, 2442.888, 2443.7, 2347.46, 2444.69, 2349.203, 2349.453, 2447.383, 2350.094, 2350.186, 2449.113, 2350.946, 2351.089, 2450.844, 2352.331, 2451.661, 2452.472, 2353.611, 2453.398, 2453.706, 2353.914, 2454.594, 2354.693, 2355.08, 2457.108, 2356.151, 2356.425, 2459.205, 2357.948, 2357.893, 2460.415, 2460.32, 2358.163, 2462.561, 2462.472, 2359.95, 
1e-100, 2.059245, 4.274909, 6.568597, 8.929901, 11.36124, 13.8635, 16.44962, 19.1126, 21.85213, 24.7799, 27.79095, 30.88009, 34.04716, 37.2934, 40.63355, 44.06389, 47.5674, 51.19084, 55.00826, 58.9176, 62.90418, 66.97422, 71.13236, 75.38602, 79.74655, 84.17243, 88.7017, 93.48238, 98.31637, 103.2587, 108.3033, 113.4457, 118.6704, 123.9777, 129.3793, 134.8592, 140.459, 146.302, 152.2385, 158.2049, 164.3504, 170.5533, 176.833, 183.1688, 189.6596, 196.2142, 202.8765, 209.5962, 216.5695, 223.6345, 230.7771, 238.0106, 245.2862, 252.6756, 260.191, 267.7518, 275.3331, 283.042, 290.8706, 298.7867, 306.859, 315.0797, 323.3774, 331.6178, 340.0339, 348.5155, 357.031, 365.5826, 374.3244, 382.9528, 391.7978, 400.6111, 409.7147, 418.858, 428.1021, 437.3746, 446.6744, 455.9852, 465.4694, 474.6511, 484.3055, 493.8388, 503.6455, 513.3059, 522.827, 532.6841, 542.5088, 552.4453, 562.5339, 572.6138, 582.7164, 592.818, 602.9313, 613.1473, 623.2876, 633.6519, 643.8822, 654.1226, 664.2275, 675.1666, 685.4502, 695.8929, 706.4928, 717.1758, 727.5834, 738.3033, 748.9638, 759.6278, 770.1181, 780.8862, 791.5764, 802.1778, 813.0124, 823.6868, 834.4024, 844.7923, 855.8581, 866.6331, 877.4211, 888.2186, 899.1624, 909.8963, 920.8979, 931.561, 942.5597, 953.318, 963.936, 974.7549, 985.5518, 996.2857, 1007.091, 1017.97, 1028.344, 1039.116, 1050.052, 1060.719, 1071.776, 1082.585, 1092.97, 1103.876, 1114.771, 1125.336, 1135.867, 1146.292, 1157.169, 1167.473, 1178.155, 1188.548, 1198.717, 1209.714, 1219.961, 1230.387, 1240.684, 1251.472, 1261.594, 1271.82, 1282.297, 1292.688, 1302.915, 1312.834, 1323.111, 1333.143, 1343.256, 1353.415, 1363.501, 1373.551, 1383.249, 1393.554, 1403.218, 1412.726, 1422.546, 1432.211, 1441.953, 1451.456, 1461.101, 1470.877, 1480.565, 1490.166, 1498.949, 1509.292, 1518.269, 1527.594, 1536.744, 1546.009, 1555.189, 1564.402, 1573.67, 1582.722, 1591.553, 1600.474, 1609.091, 1618.192, 1626.795, 1635.805, 1644.456, 1653.155, 1661.82, 1671.043, 1678.509, 1686.202, 1695.217, 1703.442, 1711.972, 1720.056, 1728.326, 1736.571, 1744.647, 1752.845, 1760.933, 1768.232, 1776.321, 1784.293, 1791.916, 1799.428, 1807.258, 1814.896, 1822.44, 1829.833, 1837.599, 1844.647, 1852.262, 1859.698, 1866.532, 1874.155, 1881.334, 1888.069, 1895.385, 1902.372, 1909.13, 1915.632, 1922.539, 1928.775, 1936.194, 1942.818, 1949.499, 1955.97, 1962.357, 1968.878, 1975.005, 1981.244, 1987.52, 1993.729, 1999.36, 2005.991, 2011.913, 2018.096, 2023.953, 2029.997, 2036.085, 2041.34, 2046.955, 2053.327, 2058.708, 2063.946, 2069.714, 2074.8, 2080.577, 2086.164, 2091.543, 2096.766, 2101.931, 2106.631, 2112.327, 2117.359, 2122.603, 2127.085, 2132.2, 2136.793, 2141.656, 2146.317, 2151.738, 2155.989, 2160.916, 2165.522, 2169.643, 2174.726, 2178.995, 2183.217, 2187.173, 2191.751, 2195.765, 2200.588, 2204.872, 2208.812, 2212.859, 2216.887, 2220.705, 2224.418, 2228.775, 2232.906, 2236.507, 2239.954, 2243.78, 2247.626, 2251.505, 2254.845, 2258.436, 2261.325, 2267.029, 2269.756, 2272.914, 2276.224, 2279.48, 2282.632, 2286.409, 2289.459, 2292.457, 2295.921, 2299.087, 2301.849, 2304.763, 2308.134, 2311.183, 2313.917, 2316.384, 2319.674, 2322.4, 2325.541, 2328.125, 2330.694, 2333.85, 2336.756, 2339.307, 2341.752, 2343.927, 2346.541, 2349.038, 2351.634, 2354.644, 2356.618, 2359.051, 2361.083, 2363.352, 2365.715, 2368.206, 2370.986, 2372.718, 2374.49, 2376.964, 2378.962, 2381.379, 2383.688, 2385.427, 2387.126, 2389.27, 2391.16, 2392.886, 2394.156, 2396.935, 2398.388, 2401.039, 2402.704, 2404.788, 2406.247, 2407.577, 2409.964, 2411.424, 2412.657, 2415.163, 2416.644, 2417.955, 2419.502, 2421.112, 2423.035, 2424.389, 2425.318, 2426.752, 2428.593, 2430.283, 2431.599, 2433.456, 2434.531, 2435.782, 2437.13, 2437.986, 2439.209, 2440.351, 2441.804, 2442.837, 2444.471, 2445.636, 2446.56, 2447.706, 2449.12, 2450.39, 2451.464, 2452.143, 2454.339, 2454.895, 2455.797, 2456.478, 2458.284, 2459.129, 2459.877, 2460.819, 2461.992, 2462.523, 2463.642, 2465.104, 2465.414, 2465.95, 2510.697, 2467.903, 2556.392, 2469.678, 2558.607, 2471.94, 2561.375, 2474.443, 2564.297, 2476.6, 2566.514, 2478.168, 2479.777, 2570.506, 2481.424, 2572.881, 2483.385, 2574.632, 2484.581, 2577.05, 2486.363, 2579.568, 2488.314, 2488.982, 2582.128, 2490.769, 2584.374, 2492.201, 2585.937, 2493.209, 2493.728, 2588.825, 2495.603, 2591.348, 2497.503, 2593.018, 2499.704, 2499.966, 2595.71, 2501.204, 2597.338, 2502.252, 2503.04, 2600.181, 2504.086, 2601.167, 2505.031, 2506.341, 2604.38, 2507.32, 2605.857, 2508.342, 2509.54, 2608.578, 2510.704, 2610.173, 2610.904, 2511.853, 2611.913, 2513.167, 2513.317, 2615.07, 2514.225, 2616.243, 2616.851, 2515.885, 2617.908, 2516.905, 2516.624, 2619.05, 2517.68, 2517.633, 2621.749, 2519.429, 2519.559, 2623.898, 2520.463, 2624.756, 2625.322, 2521.745, 2626.548, 2626.928, 2522.576, 2627.804, 2522.796, 2523.063, 2629.704, 2524.565, 2524.446, 2630.907, 2525.544, 2525.697, 2633.379, 2634.273, 2526.88, 2635.576, 2634.383, 2527.431, 
1e-100, 2.568388, 5.219036, 7.962298, 10.79198, 13.71929, 16.84217, 20.04994, 23.3484, 26.72829, 30.20297, 33.78205, 37.45119, 41.22406, 45.19118, 49.29894, 53.49253, 57.76446, 62.14763, 66.62299, 71.22299, 75.92048, 80.7024, 85.71658, 90.87339, 96.13333, 101.4739, 106.8903, 112.4626, 118.1243, 123.8807, 129.7529, 135.7293, 141.945, 148.2935, 154.7433, 161.2363, 167.8914, 174.5683, 181.3961, 188.2836, 195.347, 202.4614, 209.7408, 217.273, 224.8558, 232.4866, 240.2934, 248.136, 256.011, 264.0894, 272.1948, 280.4023, 288.6449, 297.0972, 305.7663, 314.5557, 323.2086, 332.1678, 341.0356, 350.1023, 359.1831, 368.413, 377.6459, 386.974, 396.3479, 405.9072, 415.5331, 425.3657, 435.1517, 445.1127, 455.0213, 465.0008, 475.143, 485.1669, 495.5156, 505.6849, 515.8587, 526.2156, 536.6225, 547.2936, 557.8612, 568.6941, 579.2566, 590.0932, 600.9979, 611.7938, 622.6246, 633.5796, 644.4993, 655.5403, 666.5806, 677.6767, 688.5935, 699.8836, 711.1785, 722.554, 733.6273, 745.1178, 756.2944, 768.1688, 779.3626, 790.7626, 802.1093, 813.4772, 825.0492, 836.4322, 847.7779, 859.4494, 871.0384, 882.5063, 894.2682, 905.8878, 917.257, 928.8176, 940.6639, 952.0487, 963.7687, 975.1715, 986.8809, 998.5833, 1009.966, 1021.534, 1033.204, 1044.649, 1056.144, 1067.764, 1079.403, 1090.976, 1102.358, 1114.072, 1125.506, 1137.305, 1148.23, 1159.94, 1171.36, 1182.903, 1194.008, 1205.454, 1216.9, 1228.34, 1239.593, 1250.847, 1261.809, 1273.168, 1284.447, 1295.438, 1306.634, 1317.915, 1328.595, 1340.244, 1350.941, 1362.172, 1372.934, 1383.739, 1394.725, 1405.359, 1416.182, 1426.982, 1437.952, 1448.192, 1458.711, 1469.439, 1480.099, 1490.792, 1501.21, 1511.59, 1521.877, 1532.395, 1542.455, 1552.966, 1562.942, 1573.428, 1583.09, 1592.877, 1603.359, 1613.235, 1623.041, 1633.182, 1642.654, 1652.423, 1662.114, 1671.77, 1681.245, 1691.137, 1700.459, 1709.626, 1719.566, 1728.734, 1738.069, 1747.161, 1756.577, 1765.833, 1774.491, 1783.646, 1792.759, 1801.711, 1810.81, 1820.364, 1827.789, 1835.279, 1845.019, 1853.443, 1861.968, 1870.505, 1879.265, 1887.234, 1895.563, 1903.557, 1911.911, 1920.059, 1928.354, 1936.195, 1944.165, 1951.882, 1959.726, 1967.561, 1975.348, 1983.077, 1990.724, 1997.948, 2006.11, 2013.41, 2020.496, 2027.873, 2035.045, 2041.93, 2049.608, 2056.383, 2063.288, 2070.2, 2077.046, 2083.71, 2090.625, 2097.862, 2104.338, 2110.818, 2117.54, 2124.202, 2130.189, 2136.852, 2142.953, 2149.704, 2155.806, 2161.953, 2168.308, 2174.407, 2180.507, 2186.378, 2192.607, 2198.188, 2203.763, 2209.146, 2215.464, 2220.803, 2226.422, 2231.717, 2237.492, 2242.586, 2248.079, 2253.423, 2258.97, 2263.924, 2269.344, 2274.548, 2279.546, 2284.398, 2289.462, 2294.32, 2299.033, 2303.727, 2308.512, 2313.809, 2318.173, 2322.439, 2327.354, 2332.055, 2336.277, 2340.736, 2344.874, 2349.951, 2353.798, 2358.058, 2362.039, 2366.363, 2370.342, 2374.09, 2377.949, 2382.253, 2386.148, 2390.458, 2394.131, 2398.032, 2401.739, 2405.226, 2409.113, 2412.709, 2416.251, 2419.082, 2424.313, 2427.084, 2430.458, 2433.826, 2437.277, 2440.917, 2443.824, 2446.909, 2450.18, 2453.311, 2456.189, 2459.008, 2462.477, 2465.29, 2468.633, 2471.396, 2474.298, 2477.177, 2479.785, 2483.162, 2485.553, 2488.219, 2490.817, 2493.485, 2496.318, 2498.958, 2501.413, 2504.67, 2506.898, 2509.204, 2511.522, 2513.975, 2516.388, 2518.188, 2521.234, 2523.097, 2525.816, 2528.412, 2530.309, 2532.118, 2534.521, 2536.699, 2538.619, 2540.756, 2542.692, 2544.819, 2546.818, 2548.365, 2550.209, 2551.87, 2554.477, 2556.18, 2557.685, 2559.494, 2561.558, 2563.231, 2565.263, 2567.014, 2568.419, 2570.064, 2571.566, 2572.918, 2574.268, 2575.785, 2577.437, 2579.37, 2581.132, 2582.232, 2583.975, 2585.213, 2586.741, 2588.076, 2589.342, 2590.817, 2592.606, 2593.727, 2595.115, 2596.163, 2597.732, 2598.839, 2599.615, 2601.011, 2601.967, 2603.939, 2604.516, 2606.156, 2607.186, 2607.958, 2609.014, 2610.03, 2610.71, 2611.501, 2612.772, 2613.594, 2614.843, 2616.16, 2616.79, 2618.441, 2618.039, 2619.588, 2620.81, 2621.069, 2622.678, 2668.668, 2624.325, 2715.556, 2626.184, 2718.172, 2628.436, 2720.415, 2630.619, 2723.128, 2632.527, 2725.489, 2634.767, 2635.162, 2728.483, 2636.725, 2730.678, 2638.206, 2732.314, 2640.108, 2735.427, 2642.798, 2738.33, 2644.391, 2644.773, 2740.9, 2646.37, 2742.832, 2648.197, 2744.927, 2649.486, 2649.866, 2747.175, 2650.903, 2749.657, 2653.046, 2751.001, 2654.599, 2654.565, 2753.757, 2656.399, 2756.543, 2658.622, 2658.797, 2758.403, 2659.612, 2760.251, 2660.932, 2661.538, 2763.261, 2663.081, 2764.935, 2664.377, 2664.589, 2766.566, 2665.387, 2768.107, 2768.374, 2666.821, 2770.354, 2667.751, 2668.067, 2772.428, 2669.041, 2773.741, 2774.335, 2671.52, 2776.463, 2672.315, 2672.15, 2777.681, 2673.032, 2673.375, 2779.934, 2674.162, 2674.381, 2781.306, 2675.253, 2782.364, 2782.936, 2675.781, 2783.491, 2784.401, 2677.34, 2785.578, 2678.321, 2678.675, 2788.314, 2680.347, 2680.578, 2790.271, 2680.75, 2680.741, 2791.874, 2792.06, 2681.915, 2794.03, 2792.529, 2682.747, 
1e-100, 3.076688, 6.354622, 9.7278, 13.19743, 16.76698, 20.44078, 24.22495, 28.11325, 32.1228, 36.36706, 40.71299, 45.15693, 49.70706, 54.37025, 59.15089, 64.05466, 69.04792, 74.25429, 79.62772, 85.13114, 90.71883, 96.44891, 102.2562, 108.1848, 114.2649, 120.4256, 126.7329, 133.3202, 139.9768, 146.7486, 153.6481, 160.6216, 167.6854, 174.926, 182.2826, 189.737, 197.2955, 205.101, 213.0253, 221.0837, 229.1861, 237.4319, 245.7546, 254.1955, 262.6479, 271.3212, 280.0416, 288.9277, 298.0186, 307.1674, 316.4099, 325.8257, 335.2271, 344.7047, 354.3914, 364.1033, 373.829, 383.7228, 393.6586, 403.7874, 414.0631, 424.4213, 434.7613, 445.2054, 455.7208, 466.4068, 476.9316, 487.6466, 498.4888, 509.3532, 520.37, 531.2707, 542.4512, 553.8033, 565.0166, 576.3725, 587.7576, 599.1323, 610.6113, 622.061, 633.7068, 645.2143, 657.007, 668.6674, 680.4405, 692.2423, 704.0161, 716.0375, 727.9885, 739.9307, 751.9905, 764.0856, 776.12, 788.2885, 800.4058, 812.3586, 824.5129, 836.5861, 848.5536, 861.4736, 873.3892, 885.7043, 898.0381, 910.2062, 922.5581, 934.9152, 947.2252, 959.4069, 971.923, 984.2747, 996.5553, 1008.887, 1021.131, 1033.366, 1045.584, 1057.951, 1070.352, 1082.744, 1095.125, 1107.447, 1119.629, 1131.91, 1144.165, 1156.526, 1168.585, 1180.87, 1193.102, 1205.321, 1217.308, 1229.546, 1241.652, 1253.512, 1265.267, 1277.638, 1289.418, 1301.704, 1313.602, 1325.513, 1337.339, 1349.359, 1361.116, 1372.975, 1384.421, 1396.407, 1407.911, 1419.614, 1431.199, 1442.925, 1453.964, 1465.768, 1477.396, 1488.779, 1500.376, 1511.334, 1522.619, 1533.911, 1545.23, 1556.328, 1567.603, 1578.385, 1589.388, 1600.427, 1611.217, 1622.094, 1632.842, 1644.145, 1654.561, 1664.968, 1675.696, 1686.415, 1696.504, 1707.346, 1717.376, 1727.786, 1738.379, 1748.63, 1758.597, 1768.787, 1778.956, 1788.949, 1798.786, 1808.625, 1818.77, 1828.53, 1838.284, 1847.929, 1857.559, 1867.29, 1876.308, 1886.049, 1895.381, 1904.5, 1914.086, 1923.131, 1932.13, 1941.473, 1950.977, 1960.484, 1968.332, 1975.941, 1985.198, 1994.147, 2003.068, 2011.676, 2020.381, 2028.771, 2037.276, 2045.426, 2053.76, 2061.926, 2070.178, 2078.282, 2086.117, 2094.284, 2102.472, 2110.013, 2118.073, 2125.561, 2133.452, 2141.153, 2148.754, 2156.443, 2163.568, 2171.498, 2178.099, 2185.579, 2193.047, 2199.854, 2206.864, 2213.81, 2221.102, 2228.009, 2234.965, 2241.674, 2248.397, 2255.233, 2261.858, 2268.447, 2274.259, 2280.938, 2287.232, 2293.576, 2300.285, 2306.191, 2312.208, 2318.441, 2324.575, 2330.254, 2336.844, 2342.228, 2348.245, 2354.308, 2359.651, 2365.465, 2371.2, 2376.786, 2382.008, 2387.33, 2392.861, 2398.109, 2403.068, 2408.872, 2413.923, 2418.735, 2424.031, 2428.859, 2433.852, 2438.891, 2443.526, 2448.363, 2453.638, 2457.981, 2462.52, 2467.296, 2471.83, 2476.199, 2480.468, 2484.658, 2489.249, 2493.793, 2498.071, 2502.692, 2506.473, 2510.563, 2514.572, 2518.561, 2522.461, 2526.795, 2530.579, 2534.434, 2537.994, 2541.814, 2545.955, 2549.563, 2553.211, 2557.149, 2560.473, 2563.757, 2568.357, 2571.131, 2574.236, 2577.749, 2580.835, 2584.429, 2587.996, 2590.897, 2593.985, 2597.014, 2600.09, 2602.99, 2606.458, 2609.188, 2612.023, 2614.715, 2617.45, 2621.079, 2623.572, 2626.78, 2629.472, 2632.059, 2634.896, 2637.141, 2639.469, 2641.936, 2644.578, 2647.014, 2649.905, 2652.351, 2654.52, 2656.756, 2659.18, 2661.287, 2664.2, 2666.244, 2668.51, 2670.377, 2672.514, 2674.895, 2677.171, 2679.377, 2681.164, 2683.002, 2684.887, 2686.899, 2688.597, 2690.252, 2692.637, 2694.654, 2697.099, 2698.689, 2700.528, 2702.195, 2703.854, 2705.487, 2707.583, 2709.144, 2711.107, 2712.507, 2713.943, 2715.448, 2716.723, 2718.683, 2719.939, 2721.399, 2722.724, 2724.347, 2726.118, 2727.401, 2729.167, 2730.241, 2731.292, 2732.631, 2733.794, 2734.962, 2735.837, 2737.674, 2738.31, 2740.011, 2741.276, 2743.032, 2743.748, 2744.738, 2746.115, 2747.215, 2748.02, 2749.982, 2750.893, 2751.695, 2752.529, 2753.624, 2754.553, 2755.421, 2756.289, 2757.536, 2758.581, 2759.734, 2759.43, 2761.302, 2762.52, 2762.037, 2763.034, 2809.993, 2764.534, 2858.07, 2766.748, 2861.337, 2769.615, 2863.984, 2771.831, 2866.253, 2773.275, 2868.848, 2775.576, 2776.191, 2872.172, 2777.871, 2874.174, 2779.47, 2876.217, 2781.084, 2878.797, 2782.972, 2880.357, 2784.479, 2785.158, 2883.444, 2786.64, 2885.391, 2787.773, 2887.665, 2790.041, 2790.165, 2890.509, 2792.227, 2892.721, 2794.028, 2894.76, 2795.642, 2795.406, 2897.022, 2796.919, 2898.872, 2799.01, 2799.155, 2901.12, 2799.769, 2903.133, 2801.321, 2801.9, 2905.053, 2802.881, 2907.394, 2804.61, 2804.9, 2909.468, 2805.91, 2911.457, 2911.56, 2807.313, 2913.08, 2808.57, 2808.685, 2915.666, 2809.37, 2916.607, 2916.959, 2810.66, 2918.182, 2811.451, 2811.586, 2920.149, 2813.125, 2813.521, 2922.77, 2814.404, 2814.601, 2924.35, 2815.295, 2925.098, 2925.901, 2816.196, 2926.664, 2926.986, 2817.205, 2928.675, 2818.23, 2818.16, 2930.504, 2819.072, 2819.329, 2932.332, 2819.983, 2820.491, 2934.173, 2934.071, 2821.039, 2936.096, 2935.201, 2822.742, 
1e-100, 3.713862, 7.543453, 11.49476, 15.56239, 19.79028, 24.22658, 28.77575, 33.43719, 38.21833, 43.12714, 48.15786, 53.31296, 58.59753, 64.15059, 69.82627, 75.61101, 81.52309, 87.56108, 93.71845, 100.0103, 106.4196, 112.9782, 119.7945, 126.7326, 133.8189, 140.9905, 148.2643, 155.7043, 163.259, 170.9535, 178.7323, 186.6888, 194.8869, 203.1859, 211.6048, 220.1518, 228.733, 237.4978, 246.3866, 255.3067, 264.4111, 273.5498, 283.0498, 292.5826, 302.2131, 311.9677, 321.8373, 331.7824, 341.8632, 352.0182, 362.1745, 372.5081, 382.8693, 393.5211, 404.2979, 415.1403, 425.8882, 436.9183, 447.9795, 459.1497, 470.2386, 481.5043, 492.875, 504.3675, 515.7718, 527.5264, 539.1509, 550.9755, 562.8049, 574.749, 586.6136, 598.7794, 610.8068, 622.9568, 635.132, 647.3335, 659.4895, 671.9408, 684.2746, 696.7202, 709.3283, 721.8481, 734.532, 747.1367, 759.8451, 772.6537, 785.2649, 797.9569, 810.7873, 823.4481, 836.284, 849.0862, 861.9974, 874.9658, 887.9416, 900.8283, 913.8867, 926.8094, 939.4688, 953.2863, 965.8215, 978.8663, 991.8869, 1004.88, 1017.819, 1030.747, 1043.792, 1056.8, 1069.919, 1083.028, 1095.726, 1109.124, 1121.908, 1134.75, 1147.835, 1160.582, 1173.82, 1186.692, 1199.684, 1212.292, 1225.217, 1237.85, 1250.883, 1263.329, 1276.294, 1289.072, 1301.623, 1314.537, 1327.273, 1339.731, 1352.464, 1365.124, 1377.284, 1389.988, 1402.288, 1414.931, 1427.351, 1439.646, 1451.862, 1464.09, 1476.429, 1488.497, 1500.715, 1512.645, 1524.923, 1536.92, 1549.006, 1560.552, 1572.237, 1584.746, 1596.37, 1608.158, 1619.97, 1631.386, 1643.132, 1654.689, 1666.222, 1677.587, 1689.105, 1700.032, 1711.598, 1722.68, 1733.844, 1745.042, 1756.174, 1767.388, 1778.037, 1788.98, 1799.705, 1810.762, 1821.284, 1831.971, 1842.591, 1853.265, 1863.865, 1874.094, 1884.316, 1894.628, 1904.909, 1915.125, 1924.851, 1935.095, 1945.233, 1955.087, 1964.843, 1975.242, 1984.685, 1994.534, 2003.767, 2013.936, 2023.299, 2032.672, 2041.952, 2051.254, 2060.486, 2069.923, 2079.175, 2088.948, 2096.827, 2105.022, 2114.033, 2122.857, 2131.82, 2140.815, 2149.196, 2157.617, 2165.802, 2174.545, 2183.241, 2191.216, 2199.465, 2207.735, 2215.544, 2223.935, 2231.792, 2239.958, 2247.844, 2255.329, 2263.25, 2271.103, 2278.577, 2286.004, 2293.124, 2301.103, 2307.998, 2315.342, 2322.554, 2329.466, 2336.8, 2343.894, 2350.615, 2357.461, 2364.983, 2371.58, 2378.466, 2385.175, 2391.633, 2398.093, 2404.273, 2410.95, 2417.494, 2423.969, 2430.065, 2436.118, 2442.422, 2448.618, 2454.548, 2460.051, 2466.514, 2472.041, 2477.651, 2483.835, 2489.263, 2494.64, 2500.436, 2505.971, 2511.085, 2516.783, 2522.113, 2527.662, 2533.057, 2537.95, 2543.106, 2548.349, 2553.59, 2557.776, 2562.732, 2567.566, 2572.552, 2577.165, 2582.341, 2587.006, 2591.542, 2596.153, 2600.545, 2605.052, 2609.226, 2614.155, 2618.304, 2622.286, 2626.297, 2630.757, 2634.822, 2638.796, 2642.765, 2646.89, 2651.007, 2654.814, 2659.073, 2662.935, 2666.417, 2670.084, 2673.582, 2677.677, 2681.031, 2684.315, 2687.9, 2690.988, 2695.665, 2698.739, 2701.983, 2705.59, 2708.631, 2711.49, 2714.648, 2718.009, 2720.921, 2723.717, 2727.121, 2730.31, 2733.491, 2736.715, 2739.274, 2742.175, 2744.795, 2747.509, 2750.647, 2753.281, 2755.555, 2758.549, 2761.225, 2763.876, 2766.309, 2768.9, 2771.68, 2773.785, 2776.102, 2778.476, 2780.807, 2783.532, 2785.672, 2787.768, 2790.2, 2792.685, 2794.948, 2796.966, 2798.843, 2800.697, 2802.857, 2805.114, 2807.207, 2808.96, 2810.97, 2812.782, 2814.56, 2816.499, 2818.134, 2820.541, 2822.375, 2823.822, 2825.659, 2827.604, 2829.186, 2831.017, 2832.411, 2833.795, 2835.534, 2836.95, 2838.237, 2839.795, 2841.534, 2842.784, 2845.044, 2846.333, 2847.918, 2849.063, 2850.272, 2852.202, 2853.845, 2854.585, 2856.404, 2857.654, 2858.65, 2859.845, 2861.042, 2862.513, 2863.269, 2864.584, 2865.703, 2867.194, 2868.168, 2869.644, 2870.562, 2871.487, 2872.281, 2873.421, 2874.144, 2875.012, 2876.099, 2877.14, 2877.78, 2879.736, 2880.667, 2881.428, 2882.102, 2883.466, 2883.31, 2885.096, 2886.626, 2886.38, 2887.129, 2934.596, 2888.738, 2983.642, 2890.803, 2986.16, 2893.067, 2988.613, 2894.814, 2991.087, 2896.834, 2993.265, 2898.345, 2898.833, 2996.126, 2900.424, 2998.297, 2902.825, 3001.15, 2904.534, 3003.437, 2906.382, 3005.482, 2908.022, 2908.506, 3008.024, 2910.078, 3010.335, 2911.493, 3012.048, 2913.013, 2913.197, 3015.082, 2914.87, 3016.878, 2916.475, 3018.459, 2918.231, 2918.296, 3021.686, 2920.245, 3023.22, 2921.133, 2921.327, 3025.167, 2922.582, 3026.827, 2924.312, 2924.975, 3029.802, 2925.687, 3030.997, 2926.503, 2926.75, 3033.34, 2927.66, 3034.618, 3035.122, 2929.305, 3036.874, 2930.444, 2930.508, 3039.476, 2931.374, 3040.9, 3041.325, 2933.561, 3042.594, 2934.384, 2934.297, 3044.227, 2935.462, 2935.421, 3046.076, 2936.105, 2936.28, 3047.773, 2937.185, 3048.446, 3049.311, 2938.274, 3050.789, 3051.114, 2940.201, 3053.032, 2940.995, 2941.38, 3054.724, 2941.66, 2942.012, 3056.13, 2942.507, 2942.684, 3057.791, 3057.736, 2942.747, 3059.791, 3058.147, 2944.438, 
1e-100, 4.370449, 8.952857, 13.65895, 18.49373, 23.45677, 28.55416, 33.78933, 39.15617, 44.7105, 50.49558, 56.40743, 62.44945, 68.6239, 74.91889, 81.37337, 87.96395, 94.68773, 101.6526, 108.8019, 116.0811, 123.4738, 131.0149, 138.6711, 146.4627, 154.4101, 162.4829, 170.7476, 179.2126, 187.8351, 196.5376, 205.375, 214.3492, 223.4452, 232.6728, 242.0165, 251.4599, 261.0784, 270.9304, 280.9142, 290.907, 301.097, 311.34, 321.7233, 332.2564, 342.8003, 353.4641, 364.2877, 375.266, 386.3999, 397.6464, 408.8392, 420.2815, 431.7214, 443.2837, 454.916, 466.654, 478.4016, 490.3712, 502.269, 514.4482, 526.6485, 538.9581, 551.2174, 563.7239, 576.224, 588.86, 601.2488, 613.9004, 626.6433, 639.404, 652.1241, 665.1156, 678.1502, 691.2088, 704.2912, 717.4983, 730.6332, 743.8601, 757.0973, 770.3022, 783.6253, 796.9197, 810.305, 823.6393, 837.0235, 850.4472, 864.006, 877.7046, 891.1104, 904.7136, 918.3333, 931.8375, 945.6, 959.0896, 972.6839, 986.4163, 999.8346, 1013.519, 1026.595, 1041.177, 1054.594, 1068.083, 1081.731, 1095.395, 1108.878, 1122.635, 1136.319, 1149.876, 1163.333, 1176.956, 1190.439, 1204.117, 1217.299, 1230.757, 1244.335, 1257.684, 1271.167, 1284.536, 1298.069, 1311.293, 1324.63, 1337.848, 1351.288, 1364.443, 1377.581, 1390.771, 1403.92, 1417.113, 1430.095, 1443.019, 1455.905, 1469.314, 1481.668, 1494.538, 1507.527, 1520.135, 1533.121, 1545.663, 1558.234, 1571.003, 1583.586, 1595.899, 1608.477, 1620.812, 1633.335, 1645.265, 1657.854, 1669.839, 1681.609, 1694.183, 1706.085, 1718.047, 1730.022, 1741.879, 1753.701, 1765.4, 1777.212, 1788.79, 1800.386, 1811.825, 1823.279, 1834.72, 1846.077, 1857.291, 1868.755, 1880.04, 1890.603, 1901.801, 1912.66, 1923.739, 1934.466, 1945.235, 1956.041, 1966.508, 1977.189, 1987.788, 1998.116, 2008.838, 2018.631, 2029.277, 2039.33, 2049.782, 2059.79, 2069.614, 2079.619, 2089.523, 2099.289, 2109.271, 2118.63, 2128.511, 2138.033, 2147.638, 2156.863, 2165.964, 2175.492, 2185.084, 2194.149, 2203.718, 2211.628, 2219.814, 2228.83, 2238.127, 2247.162, 2255.697, 2264.276, 2272.889, 2281.216, 2289.465, 2298.433, 2306.197, 2314.641, 2322.676, 2330.655, 2339.29, 2346.878, 2354.847, 2362.76, 2370.115, 2378.122, 2385.824, 2393.486, 2400.877, 2408.074, 2415.668, 2422.661, 2429.948, 2437.232, 2444.195, 2451.44, 2458.693, 2465.484, 2471.894, 2479.055, 2485.754, 2492.703, 2499.135, 2505.477, 2512.029, 2518.219, 2524.752, 2530.949, 2537.645, 2543.514, 2549.602, 2556.117, 2562.039, 2567.64, 2573.523, 2580.028, 2585.699, 2591.477, 2596.99, 2602.379, 2608.181, 2613.825, 2618.828, 2624.089, 2629.389, 2634.95, 2639.975, 2645.528, 2650.305, 2655.795, 2660.923, 2665.727, 2670.212, 2675.077, 2679.943, 2685.04, 2689.489, 2693.96, 2698.496, 2703.171, 2707.709, 2711.992, 2716.148, 2721.083, 2725.048, 2729.707, 2733.843, 2737.758, 2741.815, 2745.853, 2749.933, 2754.069, 2758.008, 2762.205, 2765.742, 2769.342, 2773.323, 2776.917, 2780.483, 2784.52, 2788.222, 2791.363, 2794.969, 2798.263, 2801.08, 2805.478, 2808.43, 2811.624, 2815.086, 2818.615, 2821.411, 2824.407, 2827.663, 2830.366, 2833.604, 2836.51, 2839.75, 2842.754, 2845.312, 2848.265, 2851.011, 2853.603, 2856.75, 2859.499, 2861.993, 2864.63, 2867.221, 2869.467, 2871.842, 2874.536, 2876.962, 2879.556, 2882.044, 2884.197, 2886.237, 2888.734, 2891.404, 2893.593, 2895.611, 2897.678, 2899.746, 2901.836, 2904.194, 2906.007, 2908.243, 2910.237, 2912.094, 2913.824, 2916.056, 2918.169, 2919.779, 2921.931, 2923.849, 2925.736, 2927.699, 2929.574, 2930.968, 2932.762, 2934.434, 2936.079, 2937.701, 2939.503, 2940.984, 2942.558, 2943.704, 2945.25, 2946.795, 2948.195, 2949.747, 2951.234, 2952.594, 2953.94, 2955.588, 2956.977, 2958.271, 2959.334, 2960.581, 2961.711, 2962.771, 2963.862, 2965.283, 2966.379, 2968.422, 2969.258, 2970.628, 2971.372, 2972.842, 2973.784, 2975.267, 2976.04, 2977.29, 2978.198, 2979.177, 2979.949, 2980.849, 2981.746, 2982.695, 2983.699, 2984.809, 2985.678, 2986.472, 2987.561, 2988.758, 2987.723, 2989.561, 2991.349, 2990.689, 2991.264, 3038.96, 2993.434, 3088.879, 2996.15, 3091.818, 2998.103, 3093.719, 2999.884, 3096.281, 3001.611, 3098.507, 3003.25, 3003.888, 3101.48, 3005.29, 3103.486, 3007.393, 3105.88, 3008.54, 3107.522, 3010.206, 3109.499, 3011.805, 3012.496, 3112.344, 3014.192, 3115.294, 3015.857, 3116.806, 3017.78, 3018.138, 3119.744, 3019.853, 3121.729, 3020.786, 3123.129, 3022.257, 3022.497, 3125.775, 3024.03, 3127.468, 3024.86, 3025.342, 3129.542, 3026.461, 3131.054, 3027.272, 3027.944, 3133.7, 3029.459, 3135.249, 3030.594, 3031.08, 3137.925, 3031.763, 3138.958, 3139.165, 3033.127, 3140.999, 3034.369, 3034.56, 3142.954, 3034.303, 3143.75, 3144.139, 3035.857, 3145.549, 3037.289, 3037.893, 3147.98, 3038.666, 3038.604, 3149.798, 3039.482, 3039.556, 3151.324, 3040.436, 3152.689, 3153.102, 3041.256, 3153.838, 3154.591, 3042.807, 3156.085, 3043.126, 3043.181, 3157.8, 3044.757, 3045.006, 3159.59, 3045.512, 3045.485, 3160.738, 3161.404, 3045.983, 3164.002, 3161.816, 3047.69, 
1e-100, 5.084097, 10.3166, 15.69716, 21.22225, 26.96891, 32.92296, 39.0141, 45.24332, 51.62168, 58.1421, 64.81755, 71.64167, 78.64727, 85.91211, 93.31074, 100.8484, 108.534, 116.3515, 124.3236, 132.4278, 140.6912, 149.1293, 157.8299, 166.6459, 175.6087, 184.664, 193.8774, 203.2268, 212.7179, 222.3191, 232.0743, 242.0638, 252.1586, 262.4532, 272.8108, 283.3039, 293.9161, 304.6476, 315.4772, 326.3781, 337.4996, 348.7028, 360.1237, 371.6163, 383.2397, 394.9582, 406.7555, 418.6085, 430.6345, 442.7609, 454.848, 467.1429, 479.4607, 491.9814, 504.606, 517.3187, 529.9854, 542.8687, 555.7941, 568.7623, 581.7867, 594.7854, 607.9615, 621.2486, 634.4512, 647.9075, 661.3123, 674.857, 688.5101, 702.0293, 715.6937, 729.376, 743.1863, 756.9058, 770.7675, 784.527, 798.4042, 812.2651, 826.3094, 840.2248, 854.2766, 868.4607, 882.5804, 896.5686, 910.7414, 924.8217, 938.8334, 953.0942, 967.2771, 981.3735, 995.4113, 1009.545, 1023.759, 1038.037, 1052.148, 1066.289, 1080.48, 1094.736, 1108.558, 1123.726, 1137.579, 1151.538, 1165.602, 1179.745, 1193.723, 1207.882, 1221.771, 1235.752, 1249.785, 1263.804, 1277.958, 1291.835, 1305.838, 1319.513, 1333.559, 1347.15, 1361.126, 1374.833, 1388.516, 1402.294, 1415.908, 1429.567, 1443.041, 1456.537, 1469.978, 1483.435, 1497.164, 1510.613, 1523.695, 1536.992, 1550.173, 1563.49, 1576.377, 1589.59, 1602.788, 1615.79, 1628.608, 1641.469, 1654.291, 1666.999, 1679.906, 1692.615, 1704.975, 1717.732, 1730.163, 1742.466, 1755.023, 1767.452, 1779.537, 1792.215, 1804.172, 1816.264, 1828.217, 1840.302, 1852.143, 1864.045, 1875.9, 1887.54, 1899.368, 1910.789, 1922.153, 1933.72, 1944.986, 1956.507, 1967.627, 1979.35, 1990.167, 2001.243, 2012.332, 2023.312, 2033.99, 2045.082, 2055.822, 2066.08, 2076.992, 2087.543, 2097.866, 2108.643, 2118.594, 2129.222, 2139.381, 2149.383, 2159.48, 2169.35, 2179.497, 2189.051, 2199.059, 2208.878, 2218.485, 2228.066, 2237.559, 2246.885, 2256.322, 2265.632, 2275.02, 2284.468, 2293.751, 2303.116, 2311.262, 2318.929, 2328.264, 2337.396, 2346.036, 2354.512, 2363.233, 2371.576, 2380.308, 2388.571, 2396.877, 2404.888, 2413.167, 2421.364, 2429.294, 2437.274, 2445.372, 2453.225, 2460.649, 2468.32, 2476.256, 2483.487, 2491.099, 2498.6, 2505.826, 2513.394, 2520.306, 2527.569, 2534.609, 2541.419, 2548.456, 2555.754, 2562.455, 2569.06, 2575.77, 2582.509, 2589.252, 2595.578, 2601.881, 2608.287, 2615.019, 2621.336, 2627.343, 2633.559, 2639.443, 2645.767, 2651.634, 2657.356, 2663.26, 2668.9, 2675.159, 2680.559, 2686.396, 2691.652, 2697.565, 2703.068, 2708.283, 2713.412, 2718.622, 2724.114, 2729.623, 2734.585, 2739.327, 2744.247, 2749.55, 2754.42, 2759.158, 2763.505, 2768.654, 2773.123, 2778.327, 2782.955, 2787.338, 2791.689, 2796.338, 2800.502, 2804.779, 2809.143, 2813.7, 2817.707, 2821.532, 2825.448, 2829.824, 2834.071, 2837.963, 2842.043, 2845.571, 2849.614, 2853.729, 2857.398, 2860.989, 2864.397, 2867.961, 2871.533, 2875.445, 2878.734, 2882.114, 2885.957, 2888.987, 2892.02, 2896.077, 2899.609, 2902.526, 2905.453, 2908.329, 2911.378, 2914.506, 2917.357, 2920.763, 2923.463, 2926.601, 2929.424, 2932.515, 2935.065, 2937.532, 2940.242, 2942.814, 2946.052, 2948.324, 2951.344, 2953.779, 2956.173, 2958.568, 2960.785, 2963.587, 2965.898, 2967.988, 2970.036, 2972.614, 2975.163, 2977.18, 2979.429, 2982.199, 2984.359, 2986.222, 2988.158, 2990.13, 2992.021, 2993.927, 2995.937, 2998.153, 3000.089, 3001.683, 3003.47, 3005.295, 3007.141, 3009.334, 3010.997, 3012.642, 3014.138, 3016.274, 3017.732, 3019.454, 3021.109, 3022.58, 3023.974, 3025.508, 3026.943, 3028.323, 3029.764, 3031.517, 3033.422, 3035.005, 3036.462, 3037.758, 3038.865, 3040.374, 3041.695, 3043.105, 3044.174, 3045.549, 3046.724, 3047.943, 3049.03, 3050.272, 3051.288, 3052.618, 3054.008, 3054.745, 3056.037, 3056.753, 3058.299, 3059.347, 3060.192, 3060.934, 3061.906, 3062.74, 3063.602, 3064.557, 3065.914, 3066.991, 3067.815, 3069.137, 3069.908, 3070.553, 3071.72, 3072.489, 3073.478, 3072.504, 3074.728, 3076.739, 3075.9, 3076.435, 3123.899, 3078.258, 3173.014, 3080.346, 3175.136, 3082.086, 3177.367, 3083.939, 3179.799, 3085.548, 3181.476, 3086.81, 3088.039, 3184.871, 3089.879, 3187.309, 3091.447, 3189.611, 3092.948, 3191.694, 3094.773, 3193.239, 3096.197, 3096.767, 3196.302, 3098.043, 3197.888, 3099.015, 3200.173, 3101.409, 3101.578, 3202.381, 3102.428, 3203.98, 3104.26, 3205.982, 3106.01, 3105.983, 3208.441, 3107.168, 3209.98, 3108.263, 3109.156, 3212.801, 3110.502, 3214.54, 3110.943, 3111.541, 3216.283, 3112.371, 3217.566, 3113.472, 3113.593, 3219.877, 3114.822, 3221.732, 3221.702, 3115.944, 3223.497, 3117.59, 3117.574, 3225.494, 3117.679, 3226.74, 3227.034, 3119.292, 3228.109, 3120.449, 3120.439, 3229.864, 3121.16, 3121.267, 3231.579, 3122.013, 3122.138, 3233.414, 3123.967, 3235.015, 3235.457, 3124.504, 3236.638, 3237.081, 3125.99, 3238.678, 3126.37, 3126.634, 3239.933, 3127.152, 3127.19, 3241.184, 3127.456, 3127.635, 3242.73, 3242.744, 3128.394, 3245.426, 3243.113, 3129.57, 
1e-100, 5.84816, 11.90383, 18.11083, 24.47145, 30.98599, 37.65995, 44.49398, 51.4863, 58.7194, 66.15739, 73.74563, 81.48481, 89.38162, 97.4233, 105.6363, 113.9887, 122.5078, 131.3015, 140.2294, 149.3146, 158.5339, 167.9355, 177.431, 187.0846, 196.9045, 206.8332, 217.0156, 227.3556, 237.8286, 248.4413, 259.1916, 270.0326, 281.0045, 292.1316, 303.3337, 314.6755, 326.2562, 337.9539, 349.7429, 361.6285, 373.6821, 385.7862, 398.0556, 410.3753, 422.8318, 435.2552, 447.942, 460.7185, 473.6306, 486.6739, 499.7423, 512.8244, 526.0527, 539.3143, 552.6801, 566.1417, 579.5358, 593.1825, 606.7542, 620.667, 634.4896, 648.408, 662.2833, 676.3742, 690.3901, 704.5617, 718.516, 732.7136, 746.9019, 761.1295, 775.3689, 789.7943, 804.2419, 818.6328, 833.1678, 847.6927, 862.1838, 876.6767, 891.3196, 905.8343, 920.4777, 935.0251, 949.6754, 964.2077, 978.9421, 993.4999, 1008.183, 1022.901, 1037.653, 1052.319, 1066.935, 1081.519, 1096.243, 1110.922, 1125.535, 1139.976, 1154.597, 1169.087, 1183.382, 1198.765, 1212.968, 1227.469, 1241.996, 1256.457, 1270.915, 1285.413, 1299.699, 1314.061, 1328.508, 1342.75, 1357.014, 1371.204, 1385.456, 1399.454, 1413.445, 1427.6, 1441.763, 1455.648, 1469.781, 1483.603, 1497.547, 1511.421, 1525.154, 1538.923, 1552.636, 1566.033, 1579.819, 1593.541, 1606.853, 1620.319, 1633.75, 1647.045, 1659.772, 1673.394, 1686.466, 1699.814, 1712.807, 1725.682, 1738.7, 1751.518, 1764.666, 1777.502, 1789.835, 1802.605, 1815.178, 1827.617, 1840.096, 1852.291, 1864.391, 1877.068, 1889.184, 1901.371, 1913.612, 1925.642, 1937.272, 1949.208, 1961.064, 1972.757, 1984.657, 1995.978, 2007.707, 2019.091, 2030.363, 2041.849, 2053.117, 2064.392, 2075.398, 2086.44, 2097.473, 2108.455, 2118.958, 2129.904, 2140.679, 2151.043, 2161.734, 2172.239, 2182.541, 2192.966, 2203.176, 2213.771, 2223.667, 2233.606, 2243.713, 2253.414, 2263.623, 2273.403, 2283.206, 2292.973, 2302.233, 2311.818, 2321.237, 2330.596, 2339.65, 2349.026, 2358.199, 2367.378, 2376.425, 2385.823, 2393.905, 2401.853, 2410.993, 2419.755, 2428.358, 2436.923, 2445.404, 2453.778, 2461.791, 2470.55, 2478.784, 2486.624, 2494.697, 2503.024, 2510.725, 2518.521, 2526.046, 2533.8, 2541.584, 2549.145, 2556.724, 2563.903, 2571.298, 2578.446, 2585.651, 2592.973, 2600.201, 2606.965, 2614.09, 2621.082, 2627.866, 2634.509, 2641.124, 2647.737, 2654.453, 2660.998, 2667.427, 2673.581, 2680.115, 2686.535, 2692.39, 2698.861, 2705.142, 2710.91, 2717.011, 2722.918, 2728.609, 2734.296, 2739.942, 2746.026, 2752.084, 2757.229, 2762.579, 2767.782, 2773.531, 2778.746, 2783.867, 2788.994, 2794.281, 2799.23, 2804.665, 2809.827, 2814.522, 2819.431, 2824.439, 2828.982, 2833.736, 2838.053, 2843.182, 2847.695, 2852.084, 2856.318, 2860.86, 2865.318, 2869.47, 2873.935, 2878.096, 2882.382, 2886.678, 2890.982, 2894.826, 2898.565, 2902.631, 2906.836, 2910.646, 2914.552, 2918.091, 2921.781, 2925.382, 2928.977, 2932.512, 2936.444, 2939.815, 2943.432, 2946.735, 2949.86, 2953.362, 2956.481, 2959.483, 2962.533, 2966.52, 2969.481, 2972.949, 2975.812, 2978.507, 2981.387, 2984.482, 2987.493, 2990.568, 2993.346, 2996.175, 2998.8, 3001.642, 3004.096, 3007.021, 3009.639, 3012.323, 3014.724, 3017.106, 3019.721, 3021.986, 3024.315, 3026.785, 3029.317, 3031.206, 3033.985, 3036.086, 3038.304, 3040.695, 3042.862, 3045.1, 3047.208, 3049.096, 3051.308, 3053.37, 3055.423, 3057.358, 3059.277, 3061.7, 3063.362, 3065.082, 3066.89, 3068.592, 3070.342, 3072.422, 3074.117, 3076.087, 3077.993, 3079.544, 3081.209, 3082.791, 3084.192, 3085.831, 3087.345, 3088.973, 3090.308, 3091.899, 3093.365, 3094.626, 3095.885, 3097.541, 3099.105, 3100.437, 3101.518, 3102.878, 3104.184, 3105.639, 3107.065, 3108.055, 3109.023, 3110.2, 3111.323, 3112.38, 3113.922, 3115.362, 3116.672, 3117.664, 3119.146, 3119.692, 3120.763, 3121.831, 3123.176, 3124.006, 3124.973, 3125.909, 3126.974, 3127.738, 3128.508, 3129.158, 3130.479, 3131.278, 3132.079, 3133.022, 3133.821, 3134.779, 3135.727, 3136.334, 3137.051, 3135.787, 3138.067, 3140.272, 3139.48, 3140.496, 3187.372, 3142.271, 3235.47, 3144.257, 3237.806, 3146.054, 3239.769, 3147.775, 3242.129, 3149.341, 3243.89, 3150.743, 3151.261, 3246.582, 3153.293, 3249.144, 3154.58, 3250.994, 3155.844, 3253.204, 3157.808, 3255.422, 3159.818, 3160.237, 3257.862, 3161.059, 3259.743, 3162.653, 3261.452, 3164.58, 3164.821, 3264.122, 3165.878, 3265.702, 3167.046, 3267.283, 3168.701, 3168.262, 3269.397, 3169.496, 3270.916, 3170.739, 3171.083, 3273.249, 3172.397, 3274.465, 3173.114, 3173.987, 3277.314, 3174.936, 3278.488, 3175.841, 3175.966, 3280.547, 3177.199, 3281.923, 3282.121, 3178.273, 3283.862, 3179.12, 3179.065, 3285.225, 3178.822, 3286.164, 3286.689, 3181.474, 3288.527, 3182.384, 3182.424, 3290.188, 3183.165, 3183.596, 3292.145, 3184.151, 3184.268, 3293.584, 3184.92, 3294.27, 3295.092, 3186.28, 3296.609, 3296.948, 3186.906, 3298.174, 3188.091, 3188.432, 3299.9, 3188.756, 3188.781, 3301.25, 3189.21, 3189.445, 3302.404, 3302.859, 3190.064, 3305.609, 3303.028, 3191.132
} ;


/**
 * Interpolates a value using a 4-point interpolation algorithm.
 *
 * @param ee The input value for E.
 * @param zz The input value for z.
 * @param print_flag A flag indicating whether to print debug information.
 * @return The interpolated value.
 */
double Interpolate4(double ee, double zz, int print_flag) {
      int ei = 1, zi = 1;
      int col, row;
      
      double x1, y1, x2, y2;
      double z11, z12, z21, z22;
      double a, b, c, d, t;

      if ((ee >= E4[0]) && (ee <= E4[EDIM_Finke-1]) &&
          (zz >= Z4[0]) && (zz <= Z4[ZDIM_Finke-1])) {

        if (print_flag) fprintf(stderr, "Input: E=%7.5f z=%7.5f\n\n", ee, zz);
      
        for (col=0; col <= ZDIM_Finke-2; col++) {
           if ((zz >= Z4[col]) && (zz <= Z4[col+1])) {
                 zi = col;
             break;
           }
        }
            
        for (row=0; row <= EDIM_Finke-2; row++) {
           if ((ee >= E4[row]) && (ee <= E4[row+1])) {
             ei = row;
             break;
           }
        }
        
        if (print_flag) fprintf(stderr, "Found: E_min=%7.5f E_max=%7.5f ei=%d\n",
E4[ei], E4[ei+1], ei);
        if (print_flag) fprintf(stderr, "Found: z_min=%7.5f z_max=%7.5f zi=%d\n\n",
Z4[zi], Z4[zi+1], zi);
      
        x1  = Z4[zi];
        x2  = Z4[zi+1];
        y1  = E4[ei];
        y2  = E4[ei+1];
      
        z11 = T4[ZDIM_Finke*ei+zi];           // [ei][zi];
        z12 = T4[ZDIM_Finke*(ei+1)+zi];       // [ei+1][zi];
        z21 = T4[ZDIM_Finke*ei+zi+1];         // [ei][zi+1];
        z22 = T4[ZDIM_Finke*(ei+1)+zi+1];     // [ei+1][zi+1];
      
        a = (z11 - z12 - z21 + z22)             / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        b = (y1*z12 - y2*z11 - y1*z22 + y2*z21) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        c = (x1*z21 - x2*z11 - x1*z22 + x2*z12) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        d = (x1*y1*z22 - x1*y2*z21 - x2*y1*z12 + x2*y2*z11) / (x1*y1 - x1*y2 - x2*y1
+ x2*y2);
      
        // the function
        t = a*ee*zz + b*zz + c*ee + d;
      } else {
        t = 0;

/*
        fprintf(stderr, "Incorrect range:\n");
        fprintf(stderr, "E_min=%e  E=%e  E_max=%e\n", E1[0], ee, E1[LDIM-1]);
        fprintf(stderr, "Z_min=%e  Z=%e  Z_max=%e\n", Z1[0], zz, Z1[LDIM-1]);*/
      }
      /*
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y1=%7.5f)    (x2=%7.5f, y1=%7.5f)\n", x1, y1, x2, y1);
      if (print_flag) fprintf(stderr, "(x=%7.5f, y=%7.5f)\n",
zz, ee);
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y2=%7.5f)    (x2=%7.5f, y2=%7.5f)\n\n", x1, y2, x2, y2);
     
    if (print_flag) fprintf(stderr, "(z11=%7.3f, z11_dim=%d)   (z21=%7.3f,z12_dim=%d)\n", z11, ZDIM_Finke*ei+zi,z21,ZDIM_Finke*ei+zi+1 );
      if (print_flag) fprintf(stderr, "             (z=%7.3f)\n", t);
      if (print_flag) fprintf(stderr, "(z12=%7.3f, z12_dim=%d)    (z22=%7.3f,z12_dim=%d)\n", z12, ZDIM_Finke*(ei+1)+zi,  z22,ZDIM_Finke*(ei+1)+zi+1);
         */
      return t;
}

/**
 * Calculates the optical depth, tau, using the IRA_Finke interpolation method.
 *
 * @param nu The frequency in Hz
 * @param zz The redshift value
 * @return The calculated optical depth
 */
double tau_IRA_Finke(double nu, double zz) {
      double E_keV, E_GeV, E_TeV, tau;
      
      E_keV = nu / keV;
      E_GeV = E_keV / 1.0e+6;
      E_TeV = E_keV / 1.0e+9;
      
        tau = Interpolate4(E_TeV, zz, 0);
        //fprintf(stderr, "E_TeV=%7.5f zz=%7.5f tau=%7.5f\n", E_TeV, zz, tau);
      
      return tau;
}



/**************************************************************************
*  IR absorbtion 
*
*  franceschini 2017
*
*
****************************************************************************/

const int ZDIM_Franceschini17 = 10;
const int EDIM_Franceschini17 = 56;
const int MDIM_Franceschini17 = 560;

double Z5[ZDIM_Franceschini17] = {0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};


double E5[EDIM_Franceschini17] =
{0.0052,0.00631,0.00767,0.00932,0.01132,0.01375,0.01671,0.0203,0.02466,0.02997,0.03641,0.04423,0.05374,0.06529,0.07932,0.09636,0.11708,0.14224,0.17281,0.20995,0.25507,0.30989,0.3765,
0.45742,0.55573,0.67516,0.82027,0.99657,1.21076,1.47098,1.78712,2.17122,2.63787,3.20481,3.8936,4.73042,5.7471,6.9823,8.48296,10.3061,12.5211,15.2122,18.4817,22.4539,27.2798,33.1429,
40.2661,48.9203,59.4344,72.2083,87.7276,106.582,129.489,157.319,191.131,232.21};


double T5[MDIM_Franceschini17] = 
{0,0,0,0,0,0,0,0.00001,0.00098,0.0025,
0,0,0,0,0,0,0,0.00041,0.00384,0.00706,
0,0,0,0,0,0,0.00002,0.00238,0.01031,0.01614,
0,0,0,0,0,0,0.00052,0.0077,0.02253,0.03325,
0,0,0,0,0,0,0.00269,0.01862,0.04384,0.06308,
0,0,0,0,0,0.0002,0.00834,0.03824,0.07863,0.10985,
0,0,0,0,0,0.00144,0.0201,0.07046,0.13069,0.17634,
0,0,0,0,0,0.00526,0.04148,0.11902,0.20242,0.26437,
0,0,0,0,0.00026,0.01411,0.0764,0.18647,0.29553,0.37522,
0,0,0,0.00006,0.00168,0.03132,0.12823,0.27492,0.41267,0.51353,
0,0,0,0.00075,0.00574,0.06037,0.19977,0.38853,0.56127,0.68797,
0,0.00001,0.00013,0.0032,0.01442,0.10498,0.29586,0.53788,0.75593,0.91357,
0.00004,0.00015,0.0009,0.00881,0.03012,0.17,0.42821,0.74227,1.01767,1.2124,
0.00017,0.00059,0.00279,0.01897,0.05553,0.26482,0.61799,1.02783,1.37303,1.60979,
0.00044,0.00146,0.00625,0.03557,0.09497,0.40704,0.89424,1.42494,1.85166,2.13489,
0.0009,0.00292,0.01196,0.06187,0.15657,0.62209,1.28834,1.96569,2.48455,2.81862,
0.00164,0.0053,0.02132,0.10376,0.25361,0.93967,1.83323,2.68098,3.30228,3.69168,
0.00288,0.00923,0.03649,0.17004,0.40341,1.38927,2.56014,3.59959,4.33198,4.78158,
0.00485,0.01543,0.0603,0.2716,0.62492,1.99735,3.49334,4.74466,5.59439,6.10519,
0.00785,0.02488,0.09631,0.42039,0.93622,2.78417,4.65017,6.12472,7.0919,7.66666,
0.01228,0.03875,0.14856,0.62743,1.35031,3.75962,6.02682,7.72167,8.80131,9.44462,
0.01849,0.05813,0.22045,0.89889,1.87116,4.9114,7.58655,9.48565,10.6756,11.3903,
0.02671,0.08371,0.31387,1.23409,2.49233,6.20194,9.26186,11.3434,12.6433,13.4434,
0.03695,0.11541,0.42716,1.62395,3.19434,7.56774,10.966,13.2189,14.6472,15.5431,
0.04889,0.15213,0.55633,2.05458,3.94206,8.93147,12.6258,15.0605,16.6643,17.7073,
0.06195,0.19214,0.69563,2.50232,4.69382,10.2273,14.2003,16.903,18.7701,20.0347,
0.07565,0.23399,0.8387,2.9418,5.40766,11.4154,15.7408,18.8822,21.1241,22.7128,
0.08899,0.27442,0.97492,3.34713,6.05113,12.5223,17.399,21.1526,23.9144,25.9438,
0.10138,0.31194,1.09952,3.70442,6.61381,13.64,19.3428,23.8785,27.3181,29.8841,
0.11221,0.3444,1.20615,4.00661,7.11168,14.9226,21.7192,27.2691,31.5636,34.8158,
0.1214,0.37194,1.29529,4.26781,7.60245,16.5345,24.7027,31.5446,37.0071,41.2633,
0.12901,0.3947,1.37052,4.52208,8.16693,18.5987,28.5221,37.193,44.2853,50.0201,
0.13601,0.41568,1.44294,4.8208,8.91296,21.2999,33.6172,44.8184,54.3879,62.3975,
0.14377,0.43949,1.53147,5.22872,9.95027,24.9165,40.6821,55.7377,69.2046,80.6985,
0.15419,0.47162,1.65382,5.81803,11.3558,29.901,50.8562,71.9941,91.3676,108.244,
0.16911,0.51756,1.83234,6.62927,13.2598,37.152,66.2602,96.9755,125.369,150.666,
0.19065,0.58408,2.08465,7.71468,15.8888,47.9635,90.2827,135.775,177.7,218.111,
0.21876,0.67023,2.41055,9.20728,19.6965,64.8492,128.281,195.456,257.948,338.714,
0.25607,0.78567,2.85808,11.3604,25.5316,92.0801,187.29,285.133,384.066,621.861,
0.30752,0.94942,3.50961,14.703,34.7408,135.263,275.826,415.094,615.195,1451.43,
0.38586,1.19616,4.50601,20.0272,49.645,200.724,399.269,608.09,1175.933,3855.73,
0.51166,1.59762,6.14335,28.6023,73.4936,293.144,563.797,952.984,2741.434,9943.42,
0.71277,2.23537,8.72773,42.356,109.697,414.81,786.935,1776.48,6.84E+03,1.00E+04,
1.03624,3.25844,12.901,63.5006,160.938,565.925,1156.41,3.96E+03,1.00E+04,1.00E+04,
1.55756,4.9016,19.4448,93.6672,227.793,749.334,2000.75,9.32E+03,1.00E+04,1.00E+04,
2.3373,7.33367,28.8458,133.104,308.522,1013.7,4.18E+03,1.00E+04,1.00E+04,1.00E+04,
3.40014,10.6261,41.1443,180.335,400.08,1548.18,9.40E+03,1.00E+04,1.00E+04,1.00E+04,
4.72388,14.6948,55.8835,233.223,505.738,2.91E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
6.2396,19.311,72.197,290.373,660.765,6.28E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
7.82189,24.1678,89.2541,364.018,1.01E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
9.50683,29.3591,109.08,522.73,2.00E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
11.8962,37.2086,147.522,996.47,4.51E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
18.3001,59.3783,268.848,2318.77,9.99E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
39.6796,132.463,642.223,5.41E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
99.778,333.153,1.58E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,
235.983,777.945,3.50E+03,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04,1.00E+04};





/**
 * This function performs interpolation to estimate the value of a function at a given point (ee, zz).
 * The function uses bilinear interpolation to estimate the value based on a 2D table of values.
 * The interpolation is performed within a specific range defined by the input parameters.
 *
 * @param ee      The value of the first input parameter.
 * @param zz      The value of the second input parameter.
 * @param print_flag  Determines if additional debug information should be printed.
 *                    Set to non-zero to enable printing, zero to disable printing.
 * @return        The estimated value of the function at the given point (ee, zz).
 *                If the inputs are outside the defined range, the function returns 0.
 *
 * Example usage:
 *
 * double result = Interpolate5(0.01, 0.5, 1);
 *
 * This function requires the following external variables to be defined:
 *
 * - E5          : An array representing the range of possible values for the first input parameter.
 *                 Defined as double E5[EDIM_Franceschini17].
 * - EDIM_Franceschini17    : The number of elements in the array E5. Defined as an integer constant.
 * - Z5          : An array representing the range of possible values for the second input parameter.
 *                 Defined as double Z5[ZDIM_Franceschini17].
 * - ZDIM_Franceschini17    : The number of elements in the array Z5. Defined as an integer constant.
 * - T5          : A 2D array representing the values of the function to be interpolated.
 *                 Defined as double T5[EDIM_Franceschini17 * ZDIM_Franceschini17].
 *                 The values in the array represent the function values at specific (ee, zz) coordinate pairs.
 *                 The indexing of the array is as follows:
 *                 - T5[ei][zi] represents the value at E5[ei] and Z5[zi].
 *
 * For correct usage, make sure to define the external variables before calling this function.
 */
double Interpolate5(double ee, double zz, int print_flag) {
      int ei = 1, zi = 1;
      int col, row;
      
      double x1, y1, x2, y2;
      double z11, z12, z21, z22;
      double a, b, c, d, t;

      if ((ee >= E5[0]) && (ee <= E5[EDIM_Franceschini17-1]) &&
          (zz >= Z5[0]) && (zz <= Z5[ZDIM_Franceschini17-1])) {

        if (print_flag) fprintf(stderr, "Input: E=%7.5f z=%7.5f\n\n", ee, zz);
      
        for (col=0; col <= ZDIM_Franceschini17-2; col++) {
           if ((zz >= Z5[col]) && (zz <= Z5[col+1])) {
                 zi = col;
             break;
           }
        }
            
        for (row=0; row <= EDIM_Franceschini17-2; row++) {
           if ((ee >= E5[row]) && (ee <= E5[row+1])) {
             ei = row;
             break;
           }
        }
            
        if (print_flag) fprintf(stderr, "Found: E_min=%7.5f E_max=%7.5f ei=%d\n",
E5[ei], E5[ei+1], ei);
        if (print_flag) fprintf(stderr, "Found: z_min=%7.5f z_max=%7.5f zi=%d\n\n",
Z5[zi], Z5[zi+1], zi);
      
        x1  = Z5[zi];
        x2  = Z5[zi+1];
        y1  = E5[ei];
        y2  = E5[ei+1];
      
        z11 = T5[ZDIM_Franceschini17*ei+zi];           // [ei][zi];
        z12 = T5[ZDIM_Franceschini17*(ei+1)+zi];       // [ei+1][zi];
        z21 = T5[ZDIM_Franceschini17*ei+zi+1];         // [ei][zi+1];
        z22 = T5[ZDIM_Franceschini17*(ei+1)+zi+1];     // [ei+1][zi+1];
      
        a = (z11 - z12 - z21 + z22)             / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        b = (y1*z12 - y2*z11 - y1*z22 + y2*z21) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        c = (x1*z21 - x2*z11 - x1*z22 + x2*z12) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        d = (x1*y1*z22 - x1*y2*z21 - x2*y1*z12 + x2*y2*z11) / (x1*y1 - x1*y2 - x2*y1
+ x2*y2);
      
        // the function
        t = a*ee*zz + b*zz + c*ee + d;
      } else {
        t = 0;


        //fprintf(stderr, "Incorrect range:\n");
        //fprintf(stderr, "E_min=%e  E=%e  E_max=%e\n", E1[0], ee, E1[LDIM-1]);
        //fprintf(stderr, "Z_min=%e  Z=%e  Z_max=%e\n", Z1[0], zz, Z1[LDIM-1]);
      }
      
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y1=%7.5f)    (x2=%7.5f, y1=%7.5f)\n", x1, y1, x2, y1);
      if (print_flag) fprintf(stderr, "(x=%7.5f, y=%7.5f)\n",
zz, ee);
      if (print_flag) fprintf(stderr, "(x1=%7.5f, y2=%7.5f)    (x2=%7.5f, y2=%7.5f)\n\n", x1, y2, x2, y2);
     
    if (print_flag) fprintf(stderr, "(z11=%7.3f, z11_dim=%d)   (z21=%7.3f,z12_dim=%d)\n", z11, ZDIM_Franceschini17*ei+zi,z21,ZDIM_Franceschini17*ei+zi+1 );
      if (print_flag) fprintf(stderr, "             (z=%7.3f)\n", t);
      if (print_flag) fprintf(stderr, "(z12=%7.3f, z12_dim=%d)    (z22=%7.3f,z12_dim=%d)\n", z12, ZDIM_Franceschini17*(ei+1)+zi,  z22,ZDIM_Franceschini17*(ei+1)+zi+1);
            
      return t;
}

/**
 * Calculates the absorption coefficient (tau) for a given frequency (nu) and redshift (zz) using the IRA_Franceschini17 model.
 * The absorption coefficient represents the amount of radiation that is absorbed as it travels through the universe.
 * The function performs interpolation to estimate the value of tau based on a 2D table of function values.
 *
 * @param nu The frequency in Hz.
 * @param zz The redshift.
 * @return The absorption coefficient tau.
 *
 * Example usage:
 *
 * double tau = tau_IRA_Franceschini17(4.8e+14, 0.5);
 */
double tau_IRA_Franceschini17(double nu, double zz) {
      double E_keV, E_GeV, E_TeV, tau;
      
      E_keV = nu / keV;
      E_GeV = E_keV / 1.0e+6;
      E_TeV = E_keV / 1.0e+9;
      
        tau = Interpolate5(E_TeV, zz, 0);
        //fprintf(stderr, "E_TeV=%7.5f zz=%7.5f tau=%7.5f\n", E_TeV, zz, tau);
      return tau;
}





/*
********************************************************************************
*
* PARAMETER NECESSARY FOR CALCULATIONS OF IIR ABSORPTION COEFFICIENT (Stecker & Jager 1998)
*
********************************************************************************
*/
/**
 * Calculates the absorption coefficient using the IIR model (Stecker & Jager 1998).
 *
 * @param i_i The index value to determine the row in the absorption coefficient matrix.
 * @param z_i The value for which the absorption coefficient is calculated.
 * @param level The radiation field level. Set to non-zero for a higher IR radiation field, and zero for a lower IR radiation field.
 * @return The absorption coefficient for the given parameters.
 */
double a_i(int i_i, double z_i, int level) {
      int j;
      double a_ij[4][3];
      double a_i;
      
      if (level) {
      
        //for a higher IR radiation field
	
        a_ij[0][0] = 1.46; a_ij[1][0] =  0.10; a_ij[2][0] = 0.42; a_ij[3][0] =  0.07;
        a_ij[0][1] = 1.46; a_ij[1][1] = -1.03; a_ij[2][1] = 1.66; a_ij[3][1] = -0.56;
        a_ij[0][2] = 0.15; a_ij[1][2] = -0.35; a_ij[2][2] = 0.58; a_ij[3][2] = -0.20;
	
      } else {	
      
        //for a lower IR radiation field
      
        a_ij[0][0] = 1.11; a_ij[1][0] = -0.26; a_ij[2][0] = 1.17; a_ij[3][0] = -0.24;
        a_ij[0][1] = 1.15; a_ij[1][1] = -1.24; a_ij[2][1] = 2.28; a_ij[3][1] = -0.88;
        a_ij[0][2] = 0.00; a_ij[1][2] = -0.41; a_ij[2][2] = 0.78; a_ij[3][2] = -0.31;
      }	
      
      a_i = 0.0;
      
      for (j=0; j<=2; j++) {
         a_i = a_i + (a_ij[i_i][j] * pow(log10(z_i), (double)j));
      }
      
      return a_i;
}

/*
********************************************************************************
*
* ABSORPTION COEFFICIENT WHICH DESCRIBES ABSORPTION OF GAMMA-RAYS BY
* INTERGALACTIC INFRARED MEDIUM (Stecker & Jager 1998)
*
* nu    - observed frequency
* z     - redshift
* level - absorption level (1 - high, 0 -low)
*
********************************************************************************
*/
/**
 * Calculates the absorption coefficient that describes the absorption of gamma-rays by the intergalactic infrared medium.
 *
 * @param nu The observed frequency of the gamma-ray.
 * @param z The redshift value.
 * @param level The absorption level. Use 1 for high absorption level and 0 for low absorption level.
 * @return The absorption coefficient for the given parameters.
 */
double tau_IIR(double nu, double z, int level) {

      char name[32] = "tau_IIR";

      if ((nu < MIN_FREQ) ||
          (nu > MAX_FREQ)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu'\n", name);
	return 0.0;
      }
      if ((z <= 0.0) || 
         (z >  0.3)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'z'\n", name);
	return 0.0;
      }	

      int i;
      double E_TeV;
      double log_tau;         
     
      E_TeV = nu / 2.418e+26;  
     
      log_tau = 0.0;
     
      for (i = 0; i <= 3; i++){
         log_tau = log_tau + (a_i(i, z, level) * pow(log10(E_TeV), (double)(i)));
      }
      
      return 1.0 * pow(10.0, log_tau);
}

/*
********************************************************************************
*
* ABSORPTION COEFFICIENT FOR PAIR PRODUCTION BETWEEN VHE GAMMA RAYS AND 
* SYNCHROTRON RADIATION (for SSC) (Coppi & Blandford 1990, Inoue & Takahara 1996)
*
* elec_spec   - function which defines electrons energy spectrum 
* gamma_min   - minimal energy of electrons (Lorentz factor)
* gamma_max   - maximal energy of electrons (Lorentz factor)
* nu_c        - frequency for which coefficient will be calculated
* B           - value of magnetic field
* r           - radius or length of emitting region
* sph_cyl     - 0 for cylindrical geometry 1 for spherical geometry
* prec1       - integration precision for trapezoid integration
* prec2       - integration precision for Gauss-Legendre
*
********************************************************************************
*/
/*
double gg_abs_ssc(double (*elec_spec)(double), 
	     //double NU[],
             double gamma_min, 
	     double gamma_max, 
	     double nu_c, 
	     double B, 
             double r,
             int sph_cyl,
	     int prec1,
	     int prec2
	    ) {

      char name[32] = "gg_abs_ssc";

      if ((gamma_min < 1.0) || 
          (gamma_min > 1.0e+10) ||
	  (gamma_min > gamma_max)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_min'\n", name);
	return 0.0;
      }	
      if ((gamma_max < 1.0) ||
          (gamma_max > 1.0e+10)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_max'\n", name);
	return 0.0;
      }	
      if ((nu_c < MIN_FREQ) ||
          (nu_c > MAX_FREQ)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_c'\n", name);
	return 0.0;
      }	
      if ((B < 1.0e-10) ||
          (B > 1.0e+5)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'B'\n", name);	
	return 0.0;
      }	
      if (r < 0.0) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'r'\n", name);	
	return 0.0;
      }	
      if ((prec1 < 3) ||
          (prec1 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec1'\n", name);	
	return 0.0;
      }	
      if ((prec2 < 3) ||
          (prec2 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec2'\n", name);	
	return 0.0;
      }

      double eps_c, eps_s, nu_s, n_s, jj, I_rad;
      
      eps_c  = nu_c * h / (m_e * c * c);
      eps_s  = 1.0 / eps_c;
      nu_s   = (eps_s * m_e * c * c) / h;
            
      if ((nu_s < MIN_FREQ) || (nu_s > MAX_FREQ)){
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_s'\n", name);	
	return 0.0;
      }

      
      jj  = j_syn(elec_spec, gamma_min, gamma_max, nu_s, B, prec1, prec2); 
      
      if (sph_cyl == 0)
        I_rad  = jj * r;
      else
        I_rad  = 0.75 * jj * r;  // cf. Kataoka et al., 1999
      
      n_s =  (4.0 * M_PI * I_rad) / (h * c * eps_s);
            
      return (0.2 * sig_T * n_s) / eps_c;
}

*/

/*
********************************************************************************
*
* COMPLETE PHOTON-PHOTON CROSS-SECTION (cf. Aharonian et al 2008, MNRAS 387 1206)
*
* nu_s       - target photon frequency
* nu_c       - High-Energy Photon Frequency
*
********************************************************************************
*/
/**
* Calculates the complete photon-photon cross-section according to the formula given by Aharonian et al. (2008)
*
* @param nu_s The target photon frequency
* @param nu_c The High-Energy Photon Frequency
* @return The cross-section value
*/
double sigma_gg(double nu_s, double nu_c){
  double eps_s, eps_c, s, sigma_gg;

  eps_s = h * nu_s;
  eps_c = h * nu_c;

  s = eps_s * eps_c / (m_e*m_e *c*c*c*c);

  if (s<= 1.0) {return 0.0;}
  else {
    sigma_gg = 3.*sig_T/(2.*s*s)*((s+0.5*log(s)-1./6.+1./(2.*s))*log(sqrt(s)+sqrt(s-1.))-(s+4./9.-1./(9.*s))*sqrt(1.-1./s));
    //printf("%e\n",sigma_gg);
    return sigma_gg;
  }
}

/*
********************************************************************************
*
* ABSORPTION COEFFICIENT FOR PAIR PRODUCTION BETWEEN VHE GAMMA RAYS AND
* SYNCHROTRON RADIATION
*
* nu_c        - frequency for which coefficient will be calculated
* I_c         - matrix which contains intensity of radiation field
* nu_dim      - dimension of matrix I_rad
* nu_min      - minimum frequency of radiation field photons
* nu_max      - maximum frequency of radiation field photons
* prec1       - integration precision for trapezoid integration
* prec2       - integration precision for Gauss-Legendre
*
********************************************************************************
*/

/**
 * Calculates the absorption coefficient for pair production between VHE gamma rays
 * and synchrotron radiation.
 *
 * @param nu_c     Frequency for which the coefficient will be calculated.
 * @param I_c      Array which contains intensity of radiation field.
 * @param nu_dim   Dimension of matrix I_rad.
 * @param nu_min   Minimum frequency of radiation field photons.
 * @param nu_max   Maximum frequency of radiation field photons.
 * @param prec1    Integration precision for trapezoid integration.
 * @param prec2    Integration precision for Gauss-Legendre.
 * @return         The calculated absorption coefficient.
 */
double gg_abs(double nu_c, 
		  double I_c[],
		  int nu_dim,
		  double nu_min,
		  double nu_max,
		  int prec1,
		  int prec2
		  ) {

      char name[32] = "gg_abs";

      if ((nu_c < nu_min ) ||
          (nu_c > nu_max )) {
        //fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_c': '%e', not in ['%e','%e']\n", name,nu_c,nu_min,nu_max);
	return 0.0;
      }		
      if ((prec1 < 3) ||
          (prec1 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec1'\n", name);	
	return 0.0;
      }	
      if ((prec2 < 3) ||
          (prec2 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec2'\n", name);	
	return 0.0;
      }

      double eps_c, n_s, I_rad_s;
            
      eps_c  = nu_c * h / (m_e * c * c);
                
                  
      double a, b;
      double stp, sum, val, prt;
      double x1, x2;
      const int x_dim = prec2*2;
      double x[x_dim+1], w[x_dim+1];
      int n, k, t; 



      if (eps_c < m_e*c*c/(h*nu_min)) {
	a = h * nu_min / (m_e * c * c);
      }
      
      else {
	a =  1. / eps_c;   
      }
      

      b = h * nu_max / (m_e * c * c);
      n = prec1;
      stp = (log10(b) - log10(a)) / (n - 1);
      val = log10(a);

      sum = 0.0;
	     
      for(k = 1; k <= n-1; k++) {
	       
	x1 = pow(10.,val);
	x2 = pow(10.,val + stp);
	gauleg(x1, x2, x, w, x_dim);
	    
	for(prt = 0.0, t = 1; t <= x_dim; t++) {

	  I_rad_s = linint(m_e*c*c*x[t]/h, I_c, nu_dim, nu_min, nu_max); //AZ changed nu_syn_dim+1 to nu_syn_dim 

	  n_s = (4.0 * M_PI * I_rad_s) / (h * c * x[t]);
	  
	  prt += w[t]*n_s*sigma_gg(m_e*c*c*x[t]/h ,nu_c);
	  
	}
	
	sum += prt;
	val += stp;
      }
      
      return sum;
      
}


/*
********************************************************************************
*
* ABSORPTION COEFFICIENT FOR PAIR PRODUCTION BETWEEN FOR 2nd ORDER SSC VHE GAMMA* RAYS AND 
* IC RADIATION (for SSC) (Coppi & Blandford 1990, Inoue & Takahara 1996)
*
* elec_spec   - function which defines electrons energy spectrum 
* gamma_min   - minimal energy of electrons (Lorentz factor)
* gamma_max   - maximal energy of electrons (Lorentz factor)
* nu_c        - frequency for which coefficient will be calculated
* I_rad1st    - matrix which contains intensity of 1st order IC radiation field
* nu_rad_min  - minimum frequency of radiation field photons
* nu_rad_max  - maximum frequency of radiation field photons
* nu_rad_dim  - dimension of matrix I_rad
* r           - radius or length of emitting region
* sph_cyl     - 0 for cylindrical geometry 1 for spherical geometry
* prec1       - integration precision for trapezoid integration
* prec2       - integration precision for Gauss-Legendre
*
********************************************************************************
*/
/*
double gg_abs_ssc2nd(double (*elec_spec)(double), 
		     double gamma_min, 
		     double gamma_max, 
		     double nu_c, 
		     double I_rad1st[],
		     double nu_rad_min,
		     double nu_rad_max,
		     int    nu_rad_dim,
		     double r,
		     int    sph_cyl,
		     int    prec1,
		     int    prec2
	    ) {

      char name[32] = "gg_abs_ssc2nd";

      if ((gamma_min < 1.0) || 
          (gamma_min > 1.0e+10) ||
	  (gamma_min > gamma_max)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_min'\n", name);
	return 0.0;
      }	
      if ((gamma_max < 1.0) ||
          (gamma_max > 1.0e+10)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_max'\n", name);
	return 0.0;
      }	
      if ((nu_c < MIN_FREQ) ||
          (nu_c > MAX_FREQ)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_c'\n", name);
	return 0.0;
      }	
      if (r < 0.0) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'r'\n", name);	
	return 0.0;
      }	
      if ((prec1 < 3) ||
          (prec1 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec1'\n", name);	
	return 0.0;
      }	
      if ((prec2 < 3) ||
          (prec2 > 1.0e+4)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec2'\n", name);	
	return 0.0;
      }

      double eps_c, eps_s, nu_s, n_s, jj, I_rad;
      
      eps_c  = nu_c * h / (m_e * c * c);
      eps_s  = 1.0 / eps_c;
      nu_s   = (eps_s * m_e * c * c) / h;
            
      if ((nu_s < MIN_FREQ) || (nu_s > MAX_FREQ)){
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_s'\n", name);	
	return 0.0;
      }
      
      jj = j_com(elec_spec, gamma_min, gamma_max, I_rad1st, nu_rad_min, nu_rad_max, nu_rad_dim, nu_s, prec1, prec2);

      if (sph_cyl == 0)
        I_rad  = jj * r;
      else
        I_rad  = 0.75 * jj * r;  // cf. Kataoka et al., 1999
      
      n_s =  (4.0 * M_PI * I_rad) / (h * c * eps_s);
            
      return (0.2 * sig_T * n_s) / eps_c;
      //return (sig_gg(nu_s,nu_c) * n_s);
}
*/


/*
********************************************************************************
*
* ABSORPTION COEFFICIENT FOR PAIR PRODUCTION BETWEEN VHE GAMMA RAYS AND 
* BLACKBODY RADIATION (EIC) (Inoue & Takahara 1996)
*
* elec_spec   - function which defines electrons energy spectrum 
* gamma_min   - minimal energy of electrons (Lorentz factor)
* gamma_max   - maximal energy of electrons (Lorentz factor)
* I_rad       - matrix which contains intensity of radiation field
* nu_rad_min  - minimum frequency of radiation field photons
* nu_rad_max  - maximum frequency of radiation field photons
* nu_rad_dim  - dimension of matrix I_rad
* nu_c        - frequency for which coefficient will be calculated
* r           - radius or length of emitting region
* sph_cyl     - 0 for cylindrical geometry 1 for spherical geometry
* prec1       - integration precision for trapezoid integration
* prec2       - integration precision for Gauss-Legendre
*
********************************************************************************
*/
/*
double gg_abs_eic(double (*elec_spec)(double), 
		  double gamma_min, 
		  double gamma_max, 
		  double I_rad[], 
		  double nu_rad_min, 
		  double nu_rad_max,
		  int    nu_rad_dim,
		  double nu_c,
		  double r,
		  int sph_cyl,
		  int prec1,
		  int prec2
		  ) {
  
  char name[32] = "gg_abs_eic";
  
  if ((gamma_min < 1.0) || 
      (gamma_min > 1.0e+10) ||
      (gamma_min > gamma_max)) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_min'\n", name);
    return 0.0;
  }	
  if ((gamma_max < 1.0) ||
      (gamma_max > 1.0e+10)) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'gamma_max'\n", name);
    return 0.0;
  }	
  if ((nu_c < MIN_FREQ) ||
      (nu_c > MAX_FREQ)) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_c'\n", name);
    return 0.0;
      }	
  if (r < 0.0) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'r'\n", name);	
    return 0.0;
  }	
  if ((prec1 < 3) ||
      (prec1 > 1.0e+4)) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec1'\n", name);	
    return 0.0;
  }	
  if ((prec2 < 3) ||
      (prec2 > 1.0e+4)) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'prec2'\n", name);	
    return 0.0;
  }
  
  double eps_c, n_c, jj, I_r;
  
  eps_c  = nu_c * h / (m_e * c * c);
  
  jj = j_com(elec_spec, gamma_min, gamma_max, I_rad, 
	     nu_rad_min, nu_rad_max, nu_rad_dim, nu_c,
	     prec1, prec2); 
  //if(DEBUG) cout << "DEBUG: SUBROUTINE gg_abs_eic: j_com=" << jj << endl;
  
  
  if (sph_cyl == 0)
    I_r  = jj * r;
  else
    I_r  = 0.75 * jj * r;  // cf. Kataoka et al., 1999
  
  n_c =  (4.0 * M_PI * I_r * eps_c) / (h * c);
  //if(DEBUG) cout << "DEBUG: SUBROUTINE gg_abs_eic: sig_gg=" << sig_gg(nu_c,nu_c) << endl;

  return (0.2 * sig_T * n_c) / eps_c; // OLD !!!
  //return (sig_gg(nu_c,nu_c) * n_c); // sometimes return NaN...
}
*/

/*
********************************************************************************
*
* COMPLETE PHOTON-PHOTON CROSS-SECTION, WITH FULL THOMSON + KLEIN-NISHINA
* REGIMES (cf. Aharonian et al 2008, MNRAS 387 1206)
*
* nu_s       - target photon frequency
* nu_c       - frequency of the photon tapping the target photon
*
********************************************************************************
*//*
double sig_gg(double nu_s, double nu_c,double NU[]){
  int i;
  double eps_s, eps_c, s, sig_gg, SIG;
  char   stmp[36];
  
  FILE*  stream_dat1;
  sprintf(stmp, "./data/s_matteo.dat");
  stream_dat1 = fopen(stmp, "a+");

  eps_s = h * nu_s;
  eps_c = h * nu_c;

  //s = eps_s * eps_c / (m_e*m_e *c*c*c*c);
  //if(DEBUG) cout << "DEBUG: SUBROUTINE sig_gg: s=" << s << endl;
  //debug:
  //if (s < 1.) s = 1.;
  //cout<<s<<endl;
  sig_gg = 0.;
  for (i = 1; i <= NU_DIM; i++) {
    s = h*nu_c * h*NU[i] / (m_e*m_e *c*c*c*c);
    SIG = 3.*sig_T/(2.*s*s)*(
			    (s+0.5*log(s)-1./6.+1./(2.*s))*log(sqrt(s)+sqrt(s-1))
			    -(s+4./9.-1./(9.*s))*sqrt(1.-1./s)
			    );
    if (isnan(SIG)) SIG = 0.;
    //on prend le max de SIG
    if (SIG > sig_gg) sig_gg = SIG;
    
    //simple integration
    //sig_gg += SIG/NU_DIM;
  }
  */
 /* 
  sig_gg = 3.*sig_T/(2.*s*s)*(
			      (s+0.5*log(s)-1./6.+1./(2.*s))*log(sqrt(s)+sqrt(s-1))
			      -(s+4./9.-1./(9.*s))*sqrt(1.-1./s)
			      );*//*
  //if (sig_gg < 0.0) sig_gg = 0.0;
  fprintf(stream_dat1,"%e %e\n",s,sig_gg);
  fclose(stream_dat1);

  
  sprintf(stmp, "%e", sig_gg);
  if(strcmp(stmp,"nan")  == 0 ||
     strcmp(stmp,"-nan") == 0 ||
     strcmp(stmp,"inf")  == 0 ||
     strcmp(stmp,"-inf") == 0) return 0.0;
  return sig_gg;
}
*/

/*
********************************************************************************
*
* PLANCK LAW, TO TAKE INTO ACCOUNT THE EXTERNAL INVERSE COMPTON
* RADIATION OF THE INFRARED PHOTONS REPROCESSED BY THE BROAD LINE REGIONS
*
* nu_BB      - frequency for which Planck's law will be calculated
* T_BB       - Temperature of the black body
*
********************************************************************************
*/
/**
 * Calculate the Planck function at a given frequency and temperature. Planck law, to take into account the external inverse compton. Radiation of the infrared photons reprocessed by the broad line regions
 *
 * @param nu_BB Frequency for which the Planck function will be calculated.
 * @param T_BB  Temperature of the black body.
 *
 * @return The intensity of the black body radiation in erg/s/sr/m²/Hz.
 *         Returns 0.0 if the input values are outside the valid range.
 */
double Planck(double nu_BB, double T_BB) {
  double I_BB = 0.0;
  
  char name[32] = "Planck";

  if ((nu_BB < MIN_FREQ) ||
      (nu_BB > MAX_FREQ)) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu_BB': %e\n", name,nu_BB);
    return 0.0;
  }	
  if ((T_BB < 1.0) ||
      (T_BB > 1.e+9)) {
    fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'T_BB'\n", name);
    return 0.0;
  }
  
  I_BB = 2.0 * h * pow(nu_BB , 3.0) / ( c * c ) *
    1.0 / ( exp( h * nu_BB / (k_B * T_BB) ) -1.0 );

  return I_BB; // in erg/s/sr/m2/Hz
}


/*
********************************************************************************
*
* ANALYTICAL SOLUTION OF THE TRANSFER EQUATION FOR CYLINDRICAL GEOMETRY
* CONSTANT EMISSION AND ABSORPTION ALONG EMITTING REGION
*
* I_inp       - input intensity
* jj          - emission coefficient
* kk          - absorption coefficient
* ll          - length of emitting region
*
********************************************************************************
*/
/**
 * \brief Analytical solution of the transfer equation for cylindrical geometry with constant emission and absorption coefficients along the emitting region.
 *
 * \param I_inp The input intensity.
 * \param jj The emission coefficient.
 * \param kk The absorption coefficient.
 * \param ll The length of the emitting region.
 *
 * \return The combined intensity after absorption and radiation.
 *
 * This function calculates the combined intensity after absorption and radiation for a given input intensity, emission coefficient, absorption coefficient, and length of the emitting region. It first calculates the optical depth using the length of the emitting region and the absorption coefficient. Based on the value of the optical depth, it calculates the absorbed intensity using the input intensity and the exponential decay factor. If the emission coefficient is greater than a very small threshold, it also calculates the radiated intensity using the emission and absorption coefficients.
 *
 * The implementation of this function accounts for three different cases depending on the value of the optical depth:
 * - If the optical depth is below a very small threshold, the absorbed intensity is equal to the input intensity. In this case, the radiated intensity is not considered.
 * - If the optical depth is above a large threshold, both the absorbed and radiated intensities are set to zero.
 * - If the optical depth is between the two thresholds, the absorbed intensity is calculated by multiplying the input intensity with the exponential decay factor, and the radiated intensity is calculated using the emission and absorption coefficients.
 */
double CylTransfEquat(double I_inp, double jj, double kk, double ll) {
      double tau = 0.0, I_abs = 0.0, I_rad = 0.0;
      
      tau  = ll * kk;
      
      if (tau  <  1.0e-10) {
        I_abs = I_inp;
      } else if (tau  >  7.0e+2) {                      
        I_abs = 0.0;
      } else if ((tau <= 7.0e+2)  && (tau >= 1.0e-10)) {
        I_abs = I_inp * exp(-tau);
      } 
      
      if (jj > 1.0e-300) {
        if  (tau < 1.0e-10) {                      
           I_rad = jj * ll;  
        } else if  (tau > 7.00e+2) {                       
           I_rad = jj / kk; 
        } else if ((tau <= 7.00e+2) && (tau >= 1.0e-10)) {
           I_rad = (jj / kk) * (1.0 - exp(-tau));	   
        }   
      }
      
      return I_abs + I_rad;
}

/*
********************************************************************************
*
* ANALYTICAL SOLUTION OF THE TRANSFER EQUATION FOR SPHERICAL GEOMETRY
* CONSTANT EMISSION AND ABSORPTION ALONG EMITTING REGION
*
* jj          - emission coefficient
* kk          - absorption coefficient
* ll          - length of emitting region
*
********************************************************************************
*/
/**
 * @brief Analytical solution of the transfer equation for spherical geometry.
 *
 * This function calculates the intensity of radiation given the emission and absorption coefficients
 * along with the length of the emitting region. The solution assumes constant emission and absorption
 * along the emitting region.
 *
 * @param jj The emission coefficient.
 * @param kk The absorption coefficient.
 * @param ll The length of the emitting region.
 * @return The intensity of radiation.
 *
 * @note The emission coefficient and absorption coefficient must be positive numbers.
 *
 * @warning The function does not handle cases where the emission coefficient is negative or zero,
 * or where the absorption coefficient is zero.
 *
 * @remark The intensity of radiation is computed using the following conditions:
 * - If the emission coefficient is less than 1.0e-300 or has a value of "-inf", "inf", or "nan",
 *   the intensity is set to 0.0.
 * - Otherwise, the intensity is calculated based on the length of the emitting region and the
 *   absorption coefficient. If the length is very large (tau > 7.00e+2), the intensity is
 *   set to jj divided by kk. If the length is very small (tau < 1.0e-4), the intensity is set
 *   to (4/3) times the emission coefficient times the length. Otherwise, the intensity is
 *   calculated using a formula that considers the relation between emission and absorption.
 *
 * @see SphTransfEquat2()
 */
double SphTransfEquat(double jj, double kk, double ll) {
      double tau, I_rad;
      char   stmp[36];
      
      sprintf(stmp, "%e", jj);
      if ((strcmp(stmp, "-inf") == 0) ||
          (strcmp(stmp, "inf")  == 0) ||
          (strcmp(stmp, "nan")  == 0) ||
	  (jj < 1.0e-300)) {
	   I_rad = 0.0;
      } else {
	I_rad = 0.0;
	tau  = 2.0 * ll * kk;	   
	
	if  (tau > 7.00e+2)                       I_rad = jj / kk; 
	if  (tau < 1.0e-4)                        I_rad = (4.0/3.0) * jj * ll;     
	if ((tau < 7.00e+2) && (tau > 1.0e-4))    I_rad = (jj / kk) * (1.0 - (2.0/(tau*tau)) * 
	                                                  (1.0 - exp(-tau)*(tau+1.0)));
      }
      
      return I_rad;
}




/********************************************************************************
*                                                                               *
* TRANSFORMATION OF SOURCE SURFACE INTENSITY TO OBSERVED FLUX FOR A DISK        *
* (CYLINDRICAL) SOURCE                                                          *
*                                                                               *
* Intens      - source surface intensity                                        *
* Radius      - radius of a source                                              *
* Doppler     - Doppler factor                                                  *
* z           - redshift                                                        *
* Hubble      - Hubble constant                                                 *
*                                                                               *
********************************************************************************/
/**
 * Calculates the flux from the intensity based on given parameters.
 *
 * @param Intens The intensity.
 * @param Radius The radius.
 * @param Doppler The Doppler shift.
 * @param z The redshift.
 * @param Hubble The Hubble constant.
 * @return The calculated flux.
 *
 * @pre The value of 'Radius' must be between 1.0 and 1.0e+100 (inclusive).
 * @pre The value of 'Doppler' must be between 1.0 and 1.0e+2 (inclusive).
 * @pre The value of 'z' must be greater than 0.0 and less than or equal to 5.0.
 * @pre The value of 'Hubble' must be between 50.0 and 100.0 (inclusive).
 *
 * @post If the input parameters are outside the valid ranges, the function
 *       prints an error message and returns 0.0.
 *
 * The following constants are referenced in the calculation:
 * - c: The speed of light in centimeters per second.
 *   Defined as: const double c = 2.997924 * 1.0e+10;
 *
 * The function calculates the Hubble constant H_0 by dividing 'Hubble' by 3.086
 * and multiplying by 1.0e-19. It then calculates the luminosity distance D_L
 * using the cosmological formula and the Hubble constant. Finally, it calculates
 * the flux using the given formula and returns the result.
 */
double Intens2Flux(double Intens, double Radius, 
                   double Doppler, double z, double Hubble) {

      char name[32] = "Intens2Flux";

      if ((Radius < 1.0) || 
          (Radius > 1.0e+100)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Radius'\n", name);
	return 0.0;
      }	      
      if ((Doppler < 1.0) || 
          (Doppler > 1.0e+2)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Doppler'\n", name);
	return 0.0;
      }	      
      if ((z <= 0.0) || 
          (z > 5.0)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'z'\n", name);
	return 0.0;
      }	      
      if ((Hubble < 50.0) || 
          (Hubble > 100.0)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Hubble'\n", name);
	return 0.0;
      }	

      double H_0 = (Hubble / 3.086) * 1.0e-19;
      double D_L = (2.0 * c * (z + 1.0 - sqrt(z + 1.0))) / H_0; // [cm]
      return M_PI * ((Radius*Radius) / (D_L*D_L)) * (1.0 + z) * 
             pow(Doppler,3.0) * Intens;
}




/********************************************************************************
*                                                                               *
* TRANSFORMATION OF SOURCE SURFACE INTENSITY TO OBSERVED FLUX FOR A CYLINDRICAL *
* SOURCE                                                                        *
*                                                                               *
* Intens1     - source surface intensity of the base                            *
* Intens2     - source surface intensity of the edge                            *
* Radius      - source radius                                                   *
* Length      - source length                                                   *
* Doppler     - Doppler factor                                                  *
* z           - redshift                                                        *
* Hubble      - Hubble constant                                                 *
* Theta       - angle with the line of sight                                    *
*                                                                               *
********************************************************************************/

/**
 * @brief Converts light intensity values to flux.
 *
 * This function takes the following parameters:
 *  - Intens1: Light intensity value 1.
 *  - Intens2: Light intensity value 2.
 *  - Radius: Radius of the cylinder.
 *  - Length: Length of the cylinder.
 *  - Doppler: Doppler value.
 *  - z: Redshift value.
 *  - Hubble: Hubble constant value.
 *  - Theta: Angle value in degrees.
 *
 * The function first validates the input parameters to ensure they fall within the expected ranges.
 * It validates the Radius should be greater than or equal to 1.0 and less than or equal to 1.0e+100.
 * It validates the Doppler should be greater than or equal to 1.0 and less than or equal to 1.0e+2.
 * It validates the z should be greater than 0.0 and less than or equal to 5.0.
 * It validates the Hubble should be greater than or equal to 50.0 and less than or equal to 100.0.
 * It validates the Theta should be greater than or equal to 0.0 and less than or equal to 90.0.
 *
 * If any of the parameters fail the validation, an error message will be printed to the standard error output,
 * and the function will return 0.0.
 *
 * If the parameters pass the validation, the function calculates the value of H_0 based on the Hubble constant,
 * and D_L based on the redshift value. Then it calculates and returns the flux value using the input parameters
 * and intermediate calculations.
 *
 * @param Intens1 Light intensity value 1.
 * @param Intens2 Light intensity value 2.
 * @param Radius Radius of the cylinder.
 * @param Length Length of the cylinder.
 * @param Doppler Doppler value.
 * @param z Redshift value.
 * @param Hubble Hubble constant value.
 * @param Theta Angle value in degrees.
 *
 * @return The calculated flux based on the input parameters.
 *
 * @note The value of c is a constant equal to 2.997924 * 1.0e+10.
 * @note The function assumes that the M_PI constant is defined and represents the value of pi.
 */
double CylIntens2Flux(double Intens1, double Intens2, double Radius, double Length,
                   double Doppler, double z, double Hubble, double Theta) {

      char name[32] = "Intens2Flux";

      if ((Radius < 1.0) || 
          (Radius > 1.0e+100)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Radius'\n", name);
	return 0.0;
      }	      
      if ((Doppler < 1.0) || 
          (Doppler > 1.0e+2)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Doppler'\n", name);
	return 0.0;
      }	      
      if ((z <= 0.0) || 
          (z > 5.0)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'z'\n", name);
	return 0.0;
      }	      
      if ((Hubble < 50.0) || 
          (Hubble > 100.0)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Hubble'\n", name);
	return 0.0;
      }	
      if ((Theta < 0.0) || 
          (Theta > 90.0)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Theta'\n", name);
	return 0.0;
      }	

      double H_0 = (Hubble / 3.086) * 1.0e-19;
      double D_L = (2.0 * c * (z + 1.0 - sqrt(z + 1.0))) / H_0; // [cm]
      return ((fabs((cos(Theta * M_PI / 180.0))) * M_PI * (Radius*Radius) * Intens1) +
	     (sin(Theta * M_PI / 180.0) * 2 *Radius * Length * Intens2)) * (1.0 + z) *
	     pow(Doppler,3.0);

}





/********************************************************************************
*                                                                               *
* TRANSFORMATION OF SOURCE SURFACE INTENSITY TO OBSERVED FLUX FOR RING LIKE     *
* SOURCE                                                                        *
*                                                                               *
* Intens      - source surface intensity                                        *
* InnRadius   - inner radius                                                    *
* OutRadius   - outer radius                                                    *
* Doppler     - Doppler factor                                                  *
* z           - redshift                                                        *
* Hubble      - Hubble constant                                                 *
* check       - (0 or 1) for check = 1 routine will carefully check input       *
*               parameters                                                      *
*                                                                               *
********************************************************************************/
/**
 * @brief Convert ring intensity to flux.
 *
 * This function calculates the flux from a ring intensity based on various parameters.
 *
 * @param Intens The intensity of the ring.
 * @param InnRadius The inner radius of the ring.
 * @param OutRadius The outer radius of the ring.
 * @param Doppler The Doppler value.
 * @param z The redshift value.
 * @param Hubble The Hubble value.
 * @param check A flag indicating if error checking should be performed.
 *
 * @return The flux calculated from the given ring intensity.
 *
 * @note The function checks if the given values for InnRadius, OutRadius, Doppler, z, and Hubble are within valid ranges and outputs an error message if any of them are invalid.
 *
 * @see c
 */
double RingIntens2Flux(double Intens, double InnRadius, double OutRadius, 
                       double Doppler, double z, double Hubble, int check) {

      if (check) {

        char name[32] = "RingIntens2Flux";

        if ((InnRadius < 1.0) || (InnRadius > 1.0e+100) ||
            (OutRadius < 1.0) || (OutRadius > 1.0e+100) ||
            (InnRadius >= OutRadius)
	   ) {
          fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'InnRadius' or 'OutRadius'\n", name);
	  return 0.0;
        }	      
        if ((Doppler < 1.0) || 
            (Doppler > 1.0e+2)
	   ) {
          fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Doppler'\n", name);
	  return 0.0;
        }	      
        if ((z <= 0.0) || 
            (z > 5.0)
	   ) {
          fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'z'\n", name);
	  return 0.0;
        }	      
        if ((Hubble < 50.0) || 
            (Hubble > 100.0)
	   ) {
          fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Hubble'\n", name);
	  return 0.0;
        }	

      }
     
      double H_0 = (Hubble / 3.086) * 1.0e-19;
      double D_L = (2.0 * c * (z + 1.0 - sqrt(z + 1.0))) / H_0; // [cm]
      return M_PI * ((OutRadius*OutRadius - InnRadius*InnRadius) / (D_L*D_L)) * 
             (1.0 + z) * pow(Doppler,3.0) * Intens;
}

/*
********************************************************************************
*
* FREQUENCY TRANSFORMATION FROM SOURCE FRAME TO OBSERVER FRAME
*
* nu          - source frequency
* Doppler     - Doppler factor
* z           - redshift
*
********************************************************************************
*/
/**
 * @brief Perform frequency transformation from the source frame to the observer frame.
 *
 * This function takes the source frequency, Doppler factor, and redshift as input, and calculates the transformed frequency.
 * The transformed frequency is obtained by multiplying the Doppler factor with the source frequency and dividing it by the sum of 1 and the redshift.
 * If the input values do not fall within the specified range, an error message is printed and 0.0 is returned.
 *
 * @param nu The source frequency.
 * @param Doppler The Doppler factor.
 * @param z The redshift.
 * @return The transformed frequency.
 *
 * @note The function assumes that the source frequency is in the range of 1.0e+5 to 1.0e+35.
 * @note The function assumes that the Doppler factor is in the range of 0.0 to 1.0e+2.
 * @note The function assumes that the redshift is in the range of 0.0 to 5.0.
 */
double FreqTransS2O(double nu, double Doppler, double z) {

      char name[32] = "FreqTransS2O";

      if ((nu < 1.0e+5) ||
          (nu > 1.0e+35)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu'\n", name);
	return 0.0;
      }	
      if ((Doppler < 0.0) || 
          (Doppler > 1.0e+2)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Doppler'\n", name);
	return 0.0;
      }	      
      if ((z <= 0.0) || 
          (z > 5.0)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'z'\n", name);
	return 0.0;
      }	

      return (Doppler * nu) / (1.0 + z);             
}

/*
********************************************************************************
*
* FREQUENCY TRANSFORMATION FROM OBSERVER FRAME TO SOURCE FRAME
*
* nu          - source frequency
* Doppler     - Doppler factor
* z           - redshift
*
********************************************************************************
*/
/**
 * @brief This function performs a frequency transformation from the observer frame to the source frame.
 *
 * The function takes in three parameters:
 * - `nu` represents the source frequency.
 * - `Doppler` is the Doppler factor.
 * - `z` is the redshift.
 *
 * @param nu The source frequency.
 * @param Doppler The Doppler factor.
 * @param z The redshift.
 *
 * @return The transformed frequency in the source frame.
 *
 * @note The function expects `nu` to be in the range between 1.0e+5 and 1.0e+35, `Doppler` to be in the range between 1.0 and 1.0e+2, and `z` to be in the range between 0.0 and 5.0. If any of the parameters fall outside these ranges, an error message is printed to stderr and 0.0 is returned.
 */
double FreqTransO2S(double nu, double Doppler, double z) {

      char name[32] = "FreqTransO2S";
      
      if ((nu < 1.0e+5) ||
          (nu > 1.0e+35)) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'nu'\n", name);
	return 0.0;
      }	
      if ((Doppler < 1.0) || 
          (Doppler > 1.0e+2)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'Doppler'\n", name);
	return 0.0;
      }	      
      if ((z < 0.0) || 
          (z > 5.0)
	 ) {
        fprintf(stderr, "SUBROUTINE: '%s' ERROR: Wrong value of 'z'\n", name);
	return 0.0;
      }	

      return (nu * (1.0 + z)) / Doppler;             
}


/**
 * Calculates the luminosity distance based on the redshift value.
 *
 * The calculation is based on the following formula:
 *
 * dl = c/H0 * ( z + (z*z*(1.-q0))/(1.+q0*z+sqrt(1.+2.*q0*z)) )
 *
 * where:
 * - dl is the luminosity distance
 * - c is the speed of light constant
 * - H0 is the Hubble constant multiplied by the relevant conversion factors for distance units
 * - z is the redshift value
 * - q0 is the deceleration parameter
 *
 * @param z The redshift value
 * @return The calculated luminosity distance
 */
double LuminDist(const double z){
  double dl = 0.;
  dl = c/H0 * ( z + (z*z*(1.-q0))/(1.+q0*z+sqrt(1.+2.*q0*z)) );
  return dl;
}


/*
********************************************************************************
*
* DISTANCE LUMINOSITY CALCULATOR
*
* adapted from a script by Ned Wright
* https://www.astro.ucla.edu/~wright/CosmoCalc.html
*
* z           - redshift
* H0          - Hubble constant
* WM          - Omega matter
*
* return distance luminosity in cm
********************************************************************************
*/
/**
 * @brief Calculates the distance luminosity in cm.
 *
 * This function calculates the distance luminosity in centimeters based on the redshift (z),
 * Hubble constant (H0), and Omega matter (WM). The function is adapted from a script by Ned Wright.
 * https://www.astro.ucla.edu/~wright/CosmoCalc.html
 *
 * @param z The redshift.
 * @param H0 The Hubble constant.
 * @param WM The Omega matter.
 * @return The distance luminosity in cm.
 */
double Distance_Luminosity(double z, double H0, double WM) {

    // initialize constants
    
      double WV = 1.-WM;    		// consider flat universe
      double h = H0/100.;
      double WR = 4.165E-5/(h*h);	// Omega(radiation), includes 3 massless neutrino species, T0 = 2.72528
      double WK = 1-WM-WR-WV;		// Omega curvaturve = 1-Omega(total)
      double c = 299792.458;		// velocity of light in km/sec
      double DCMR = 0.0;     		// comoving radial distance in units of c/H0
      double DA = 0.0;       		// angular size distance
      double DL = 0.0;      		// luminosity distance
      double DL_Mpc = 0.0;		// luminosity distance in Mpc
      double a = 1.0;			// 1/(1+z), the scale factor of the Universe
      double az = 1.0/(1+1.0*z);
      double ratio = 1.0;
      double adot = 0.;
      int n = 1000;        		// number of points in integrals
      int i = 0;
      double x,y,DCMT;


    // do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
      for(i = 1; i <= n; i++){
        a = az+(1-az)*(i+0.5)/n;
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
        DCMR = DCMR + 1./(a*adot);
      }
      DCMR = (1.-az)*DCMR/n;
    
    // tangential comoving distance
      x = sqrt(abs(WK))*DCMR;
      if(x > 0.1){
        if(WK > 0){
          ratio =  0.5*(exp(x)-exp(-x))/x;
        }
        else ratio = sin(x)/x;
      }
      else {
        y = x*x;
        if (WK < 0) y = -y;
        ratio = 1. + y/6. + y*y/120.;
       }
      DCMT = ratio*DCMR;
      DA = az*DCMT;
      DL = DA/(az*az);
      DL_Mpc = (c/H0)*DL;

      return DL_Mpc * 1.0e6 * pc;
}

/*
********************************************************************************
*
*    FUNCTIONS FOR HADRONIC INTERACTIONS (PION DECAY)
*
********************************************************************************
*/
/**
 * @brief This function sets the parameters for hadronic interactions, specifically for pion decay.
 *
 * The function calculates and sets various parameters needed for the hadronic interactions' calculations. It computes the energetics of electrons to normalize the proton spectrum, calculates the normalization of the proton spectrum from ENERGY50, and sets other related parameters.
 *
 * The function does not take any input parameters and does not return any values. It modifies the global variables ENERGY50, fEnergy_prot, fNorm_prot, and fNb_prot.
 */
void setHadronicParameters(){


  double t,t1,t2,step,per_decade,interval,spec_e,number_e,power_e,gamma;
  per_decade=1.0e+03;

  // Compute energetics of electrons, needed for normalization of proton spectrum
  t1=log10(GAMMA_MIN); t2=log10(GAMMA_MAX);
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  number_e=0.;
  power_e=0.;
  for(t=t1+step/2;t<t2;t+=step){
    gamma=pow(10.,t);
    spec_e= N_e(gamma)*gamma*interval;
    number_e+=spec_e;
    power_e+=gamma*spec_e;
  }
  double Volume=4./3.*M_PI*pow(R_src,3);
  // ratio = nb_p / nb_e
  // cf. energetics:
  // nb_p = (E50*1.e50*ERG/eV)/MPROTON
  // nb_e = m_e*c*c*power_e*ERG/eV*Volume/(511.e3)
  ENERGY50=(RATIO*MPROTON*m_e*c*c*power_e*Volume)/(511.e3*1.0e+50);
  



  
  // calculate normalization of proton spectrum from ENERGY50
  
  double e,number,power,spec;
  
  t1=log10(EPMIN); t2=log10(EPMAX);
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  fNorm_prot=1.;
  number = 0;
  power = 0;


  for(t=t1+step/2;t<t2;t+=step)
    {
      e=pow(10.,t);
      spec= N_p(e)*e*interval;
      number+=spec;
      power+=e*spec;
    }

  fEnergy_prot = power/number;
  fNorm_prot = (ENERGY50*1.e50*ERG/eV)/power;
  fNb_prot = (ENERGY50*1.e50*ERG/eV)/fEnergy_prot;
  
  if(DEBUG){
    cout << endl;
    cout << "DEBUG: Kp = " << Kp << endl;
    cout << "DEBUG: Energy50 = " << ENERGY50*1.e50*ERG/eV << " eV" << endl;
    cout << "DEBUG: fEnergy_prot = " << fEnergy_prot << " eV" << endl;
    cout << "DEBUG: number = " << number << endl;
    cout << "DEBUG: power = " << power << endl;
    cout << "DEBUG: fNorm_prot = " << fNorm_prot << endl;
    cout << "DEBUG: Number protons = " << fNb_prot << endl;
    cout << "DEBUG: Volume  = " << 4./3.*M_PI*pow(R_src,3) << " cm^3" << endl;
    cout << "DEBUG: Test Kp = " << (ENERGY50*1.e50*ERG/eV)/(power*4./3.*M_PI*pow(R_src,3)) << " cm^-3" << endl;
    cout << "DEBUG: Test Kp = " << fNb_prot/(4./3.*M_PI*pow(R_src,3)*number) << " cm^-3" << endl;
    cout << "DEBUG: Test Kp = " << fNorm_prot/(4./3.*M_PI*pow(R_src,3)) << " cm^-3" << endl;
    cout << "DEBUG: Test Norm = " << 4./3.*M_PI*pow(R_src,3)*Kp << endl;
  }
  
  return;
}





/**
 * Calculates the flux of gamma ray spectra from the decay of pions.
 * The calculations are based on equations from Kelner, Aharonian & Bugayov, PhRvD, 74, 034018, 2006.
 * The function takes the energy of the gamma ray as input and returns the calculated flux.
 *
 * @param energy The energy of the gamma ray in eV.
 * @return The flux of gamma ray spectra in erg cm<sup>-2</sup> s<sup>-1</sup>.
 */
double piondecay(double energy)
{
  // All equations are taken from Kelner, Aharonian & Bugayov, PhRvD, 74, 034018, 2006: 2006PhRvD..74c4018K
  
  DISTANCE = LuminDist(z); // redshift -> cm
  fFluxfact = 1./(4*M_PI)*pow(DISTANCE,-2.)*pow(DOP_B,4); // cm^-2

  // Energy unit: eV
  double per_decade=1.0e+04; // number of steps in integration
  double step,interval,t_proton,t_pion,energy_pion,energy_used,energy_proton,x,L,nde,t1,t2,flux,sigmapp,fgamma,bgamma,betagamma,kgamma,xbeta,eseuil,fac;
  
  t1=log10(energy); t2=log10(EPMAX);
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  flux = 0; // flux=dNgamma/dEgamma

  double Kpi=0.17;

  
  if(energy < 5.e9){ // Egamma < 5 GeV
    t1 = log10(energy + MPION*MPION/(4*energy));
    for(t_pion=t1+step/2.;t_pion<log10(100e9);t_pion+=step) // integral from E to 100 GeV
      {
	energy_pion=pow(10.,t_pion);
	energy_used = MPROTON + energy_pion/Kpi; // mp + Epion/Kpi, Kpi=0.17
	L = log(energy_used/1e12); // energy_used: proton energy in TeV
	eseuil = pow(1.22e9/energy_used,4); // threshold energy for pi0 production
	sigmapp = (34.3 + 1.88*L + 0.25*L*L)*pow(1-eseuil,2)*MILLIBARN; // cf. eq 79
	if(energy_pion/Kpi< 1e9) // approx not valid below 1 GeV -> sigmapp=0
	  sigmapp = 0;
	nde=N_p(energy_used)*interval*energy_pion;
	fac = -0.48*NP + 2.06; // ??? facteur eta ???
	//fac = 1.0;// !!! TEST, REMOVE THIS LINE !!! IN FACT, THIS SEEMS TO BE BETTER RACCORDING LOW AND HIGH ENERGY WITH fac=1 THAN WITH THE FROMULA ABOVE...
	flux += 2*nde*sigmapp/sqrt(energy_pion*energy_pion - MPION*MPION)*fac/Kpi; // cf. eq 77 & 78, Kpion=0.17
      }  
  }
  //if(energy >= 5.e9){ // Egamma > 5 GeV
  else{  // Egamma > 5 GeV
    for(t_proton=t1+step/2.;t_proton<t2;t_proton+=step) // integral from E to Emax_prot
      {
	energy_proton = pow(10,t_proton);
	L = log(energy_proton/1e12);
	nde=N_p(energy_proton)*interval;
	eseuil = pow(1.22e9/energy_proton,4);
	if(energy_proton < 1.22e9) continue; // 1.22 GeV = threshold energy for pi0 production
	sigmapp = (34.3 + 1.88*L + 0.25*L*L)*pow(1-eseuil,2)*MILLIBARN; // cf. eq. 79
	bgamma = 1.30 + 0.14*L + 0.011*L*L; // cf. eq 59
	betagamma = 1./(1.79+0.11*L+0.008*L*L); // cf. eq 60
	kgamma = 1./(0.801+0.049*L+0.014*L*L); // cf. eq 61
	x = energy/energy_proton; // Egamma/Eproton
	xbeta = pow(x,betagamma);
	fgamma = bgamma*log(x)/x*pow((1-xbeta)/(1+kgamma*xbeta*(1-xbeta)),4)*(1./log(x) - 4*betagamma*xbeta/(1-xbeta) - (4*kgamma*betagamma*xbeta*(1-2*xbeta))/(1+kgamma*xbeta*(1-xbeta))); // cf. eq 58
	
	if(energy_proton > 10e9) // Eproton > 10 GeV
	  flux += sigmapp*fgamma*nde; // cf. eq 71
      }
  }
  flux *= 1.45*fFluxfact*ATOMICDENSITY*LIGHTSPEED; // fFluxfact: emissivity -> observed flux, dN/dE
  //The factor 1.45 comes from Dermer et al. 1986 (take into account interaction between proton and nuclei)
  // cf. eq 71
  
  //cout << " flux = " << flux << endl;
  // OLD: return flux*energy*energy; // in eV m-2 s-1

  double result=flux*energy*energy; // E^2 dN/dE

  //if(DEBUG) cout << "flux=" << flux << "eV-1 cm-2 s-1" << endl;
  //if(DEBUG) cout << "sigmapp=" << sigmapp << endl;
  //if(DEBUG) cout << "fFluxfact=" << fFluxfact << endl;

  return result*eV/ERG; // eV cm-2 s-1 -> erg cm-2 s-1
}





/**
 * Calculates the hadronic parameters for testing purposes.
 *
 * This function calculates the normalization of the proton spectrum and sets
 * the values of the fNorm_prot and Kp variables. It performs an integration
 * from 1 TeV to infinity to calculate the normalization. The integration is
 * done using a step size of 1/(per_decade) and the result is stored in the
 * fNorm_prot variable. The Kp variable is also set to the same value. The
 * integration is done using the N_p function, which calculates the value of
 * the N_p function for a given energy.
 *
 * @return None.
 */
void setHadronicParametersTest(){

  
  // calculate normalization of proton spectrum
  
  double per_decade=1e06;
  double step,e,t1,t2,t,number,power,spec,interval;
  
  t1=log10(1.0e12); t2=log10(EPMAX*1.e3); // integration from 1 TeV to infinity
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  fNorm_prot=1.;
  Kp=1.;
  number = 0;
  power = 0;


  for(t=t1+step/2;t<t2;t+=step)
    {
      e=pow(10.,t);
      spec= N_p(e)*e*interval;
      number+=spec;
      power+=e*spec;
    }

  
  double A=(1.0*ERG/eV)/power;
  //double KpAndrea=107.; // cm^-3
  fNorm_prot=A;
  Kp=A;
  cout << "DEBUG: A=" << A << endl;

  
  return;
}



/**
 * Calculates the decay test for a given energy.
 *
 * This function calculates the decay test for a given energy using equations from the paper "Kelner, Aharonian & Bugayov, PhRvD, 74, 034018, 2006: 2006PhRvD..74c4018K".
 * The function takes the energy in electron volts (eV) as input and returns the result of the decay test.
 *
 * @param energy The energy in eV for which to calculate the decay test.
 * @return The result of the decay test for the given energy.
 */
double piondecaytest(double energy)
{
  // All equations are taken from Kelner, Aharonian & Bugayov, PhRvD, 74, 034018, 2006: 2006PhRvD..74c4018K
  

  // Energy unit: eV
  double per_decade=1.0e+04;
  double step,interval,t_proton,t_pion,energy_pion,energy_used,energy_proton,x,L,nde,t1,t2,flux,sigmapp,fgamma,bgamma,betagamma,kgamma,xbeta,eseuil,fac;
  
  t1=log10(energy); t2=log10(EPMAX);
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  flux = 0; // flux=dNgamma/dEgamma

  double Kpi=0.17;

  
  if(energy < 5.e9){ // Egamma < 5 GeV
    t1 = log10(energy + MPION*MPION/(4*energy));
    for(t_pion=t1+step/2.;t_pion<log10(100e9);t_pion+=step) // integral from E to 100 GeV
      {
	energy_pion=pow(10.,t_pion);
	energy_used = MPROTON + energy_pion/Kpi; // mp + Epion/Kpi, Kpi=0.17
	L = log(energy_used/1e12); // energy_used: proton energy in TeV
	eseuil = pow(1.22e9/energy_used,4); // threshold energy for pi0 production
	sigmapp = (34.3 + 1.88*L + 0.25*L*L)*pow(1-eseuil,2)*MILLIBARN; // cf. eq 79
	if(energy_pion/Kpi< 1e9) // approx not valid below 1 GeV -> sigmapp=0
	  sigmapp = 0;
	nde=N_p(energy_used)*interval*energy_pion;
	fac = -0.48*NP + 2.06; // ??? facteur eta ???
	//fac = 1.0;// !!! TEST, REMOVE THIS LINE !!! IN FACT, THIS SEEMS TO BE BETTER RACCORDING LOW AND HIGH ENERGY WITH fac=1 THAN WITH THE FROMULA ABOVE...
	flux += 2*nde*sigmapp/sqrt(energy_pion*energy_pion - MPION*MPION)*fac/Kpi; // cf. eq 77 & 78, Kpion=0.17
      }  
  }
  //if(energy >= 5.e9){ // Egamma > 5 GeV
  else{  // Egamma > 5 GeV
    for(t_proton=t1+step/2.;t_proton<t2;t_proton+=step) // integral from E to Emax_prot
      {
	energy_proton = pow(10,t_proton);
	L = log(energy_proton/1e12);
	nde=N_p(energy_proton)*interval;
	eseuil = pow(1.22e9/energy_proton,4);
	if(energy_proton < 1.22e9) continue; // 1.22 GeV = threshold energy for pi0 production
	sigmapp = (34.3 + 1.88*L + 0.25*L*L)*pow(1-eseuil,2)*MILLIBARN; // cf. eq. 79
	bgamma = 1.30 + 0.14*L + 0.011*L*L; // cf. eq 59
	betagamma = 1./(1.79+0.11*L+0.008*L*L); // cf. eq 60
	kgamma = 1./(0.801+0.049*L+0.014*L*L); // cf. eq 61
	x = energy/energy_proton; // Egamma/Eproton
	xbeta = pow(x,betagamma);
	fgamma = bgamma*log(x)/x*pow((1-xbeta)/(1+kgamma*xbeta*(1-xbeta)),4)*(1./log(x) - 4*betagamma*xbeta/(1-xbeta) - (4*kgamma*betagamma*xbeta*(1-2*xbeta))/(1+kgamma*xbeta*(1-xbeta))); // cf. eq 58
	
	if(energy_proton > 10e9) // Eproton > 10 GeV
	  flux += sigmapp*fgamma*nde; // cf. eq 71
      }
  }
  flux *= ATOMICDENSITY*LIGHTSPEED; // fFluxfact: emissivity -> observed flux, dN/dE
  //The factor 1.45 comes from Dermer et al. 1986 (take into account interaction between proton and nuclei)
  // cf. eq 71
  
  //cout << " flux = " << flux << endl;
  // OLD: return flux*energy*energy; // in eV m-2 s-1

  double result=flux*energy*energy; // E^2 dN/dE

  //if(DEBUG) cout << "flux=" << flux << "eV-1 cm-2 s-1" << endl;
  //if(DEBUG) cout << "sigmapp=" << sigmapp << endl;
  //if(DEBUG) cout << "fFluxfact=" << fFluxfact << endl;

  return result/1.e12; // eV s-1 -> TeV s-1
}






//-----------------------------------------------------------------------------

/**
 * Calculates the energy spectrum for electron emission from pion decay.
 *
 * This function calculates the energy spectrum for electron emission resulting
 * from pion decay. The input energy parameter represents the electron energy in
 * electron volts (eV). The output is the calculated energy spectrum for the
 * given electron energy.
 *
 * @param energy The energy of the electrons in eV for which to calculate the energy spectrum.
 * @return The calculated energy spectrum for the given electron energy.
 */
double piondecay_ahaprecise_electron(double energy)
{
  double per_decade=100;
  double step,interval,t_proton,t_pion,energy_pion,energy_used,energy_proton,x,L,nde,t1,t2,flux,sigmapp,felec,belec,betaelec,kelec,eseuil;

  t1=log10(energy); t2=log10(EPMAX);
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  flux = 0;
  if(energy < 4e10){
    t1 = log10(energy + MPION*MPION/(4*energy));
    for(t_pion=t1+step/2.;t_pion<t2;t_pion+=step)
      {
	energy_pion=pow(10.,t_pion);
	energy_used = MPROTON + energy_pion/0.17;
	L = log(energy_used/1e12);
	eseuil = pow(1.22e9/energy_used,4);
	sigmapp = (34.3 + 1.88*L + 0.25*L*L)*pow(1-eseuil,2)*MILLIBARN;
	if(energy_pion/(0.17)< 1e9)
	  sigmapp = 0;
	nde=N_p(energy_used)*interval*energy_pion;
	flux += 0.77*nde*sigmapp/sqrt(energy_pion*energy_pion - MPION*MPION)/0.17;
      }
  }
  else{
    for(t_proton=t1+step/2.;t_proton<t2;t_proton+=step)
      {
	energy_proton = pow(10,t_proton);
	L = log(energy_proton/1e12);
	nde=N_p(energy_proton)*interval;
	eseuil = pow(1.22e9/energy_proton,4);
	belec = 1./(69.5 + 2.65*L + 0.3*L*L);
	betaelec = 1./pow(0.201 + 0.062*L + 0.00042*L*L,0.25);
	kelec = (0.279+0.141*L+0.0172*L*L)/(0.3+pow(2.3+L,2));
	x = energy/energy_proton;
	felec = belec*pow(1+kelec*pow(log(x),2),3)*pow(-log(x),5)/(x*(1+0.3/pow(x,betaelec)));
	sigmapp = (34.3 + 1.88*L + 0.25*L*L)*pow(1-eseuil,2)*MILLIBARN;
	
	flux += sigmapp*felec*nde;
      }
 }
  flux *= 1.45*fFluxfact*ATOMICDENSITY*LIGHTSPEED;
  
  //if(DEBUG) cout << " flux = " << flux*energy*energy << " energy = " << energy << endl;
  return flux*energy*energy;
  
}


/**
 * Calculate the energetics of the system.
 *
 * This function calculates the energetics of the electrons and protons in the system. It computes various parameters
 * such as the number of particles, energy density, equipartition parameter, mean energy, and total energy. It also
 * prints out the derived parameters if the PRINT flag is set to 1.
 *
 * @return void
 */
void energetics(){
  double t,t1,t2,step,per_decade,interval,elecEnergyDensity,spec_e,number_e,power_e,gamma,magEnergyDensity,nbElectron,eB,meanElectronEnergy,totalElectronEnergy,spec_p,number_p,power_p,Volume,energy;
  per_decade=1.0e+03;

  // Compute energetics of electrons
  t1=log10(GAMMA_MIN); t2=log10(GAMMA_MAX);
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  number_e=0.;
  power_e=0.;
  //integration
  for(t=t1+step/2;t<t2;t+=step){
    gamma=pow(10.,t);
    spec_e= N_e(gamma)*gamma*interval;
    number_e+=spec_e;
    power_e+=gamma*spec_e;
  }
  Volume=4./3.*M_PI*pow(R_src,3);
  elecEnergyDensity=m_e*c*c*power_e; // erg cm^-3
  magEnergyDensity=B*B/(8.*M_PI); // erg cm^-3
  eB=magEnergyDensity/elecEnergyDensity;
  meanElectronEnergy=power_e/number_e*m_e*c*c*ERG/eV; // eV
  totalElectronEnergy=elecEnergyDensity*Volume; // erg
  //nbElectron=elecEnergyDensity/(m_e*c*c)*Volume;
  nbElectron=elecEnergyDensity*ERG/eV/(511.e3)*Volume;
  //densityElectron=nbElectron/Volume; // cm^-3

  //fEnergy_prot = power/number;
  //fNb_prot = (ENERGY50*1.e50*ERG/eV)/fEnergy_prot;


  // Compute energetics of protons
  t1=log10(EPMIN); t2=log10(EPMAX);
  step=1./per_decade;
  interval=pow(10.,step/2)-pow(10.,-step/2);
  number_p=0.;
  power_p=0.;
  for(t=t1+step/2;t<t2;t+=step){
    energy=pow(10.,t);
    spec_p= N_p(energy)*energy*interval;
    number_p+=spec_p;
    power_p+=energy*spec_p;
  }
  double protonEnergyDensity,nbProton,meanProtonEnergy,densityProton,totalEnergyProton;
  totalEnergyProton=ENERGY50*1.0e+50; // erg
  protonEnergyDensity=totalEnergyProton/Volume; // erg cm^-3
  meanProtonEnergy=power_p/number_p; // eV
  nbProton=protonEnergyDensity*ERG/eV/MPROTON*Volume;
  densityProton=nbProton/Volume;

  if(PRINT){
    cout << endl
	 << "DERIVED PARAMETERS FROM THE BLOB" << endl
	 << "----------" << endl
	 << endl
	 << "Volume of source        = " << Volume              << " cm^3"      << endl
	 << "Radius of source        = " << R_src/pc            << " pc"        << endl
	 << endl
	 << "Electrons:"                                                        << endl
	 << "----------"                                                        << endl
	 << "Number of electrons     = " << nbElectron                          << endl
      //<< " TEST nb electron = " << totalElectronEnergy*ERG/eV/meanElectronEnergy << endl
	 << "Electron energy density = " << elecEnergyDensity   << " erg cm^-3" << endl
	 << "Magnetic energy density = " << magEnergyDensity    << " erg cm^-3" << endl
	 << "Equipartition parameter = " << eB                                  << endl
	 << "Mean electron energy    = " << meanElectronEnergy  << " eV"        << endl
	 << "Total electron energy   = " << totalElectronEnergy << " erg"       << endl
	 //<< "Density of electron     = " << densityElectron     << " cm^-3"     << endl
	 << endl;
    if(CASE_PION){
    cout << "Protons:"                                                          << endl
	 << "--------"                                                          << endl
	 << "Number of protons       = " << nbProton                            << endl
      //<< " TEST nb proton = " << (ENERGY50*1.0e+50*ERG/eV)/meanProtonEnergy << endl
	 << "Proton energy density   = " << protonEnergyDensity << " erg cm^-3" << endl
	 << "Mean proton energy      = " << meanProtonEnergy    << " eV"        << endl
	 << "Total proton energy     = " << totalEnergyProton   << " erg"       << endl
	 << "Density of proton       = " << densityProton       << " cm^-3"     << endl
	 << "" << endl
      ;
    }
  }

  return;
}
