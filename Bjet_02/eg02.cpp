// PROGRAM CALIBRATE ELLIPTICAL GALAXY SPECTRUM FOR GIVEN
// DISTANCE, GALAXY MASS AND GALAXY AGE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <iostream.h> //standard
#include <iostream> //OH
#include <errno.h>

const int    NU_DIM_MAX = 300;
      int    NU_DIM     = 1;
     char    InpFile[256];
     char    NameSS[256];
     
   double    xs[NU_DIM_MAX+1];
   double    ys[NU_DIM_MAX+1];
   double    xg[1300];
   double    yg[1300];     

// CONSTANTS

const double c      = 2.997924 * 1.0e+10; // [cm/s]
const double pc     = 3.086    * 1.0e+18; // [cm]

double x[2000], y[2000];

// PHYSICAL PARAMETERS

double z;
double H_0;
double d_l;
double M_G;
int    AGE;

int load_params(char* name); 

double galflux(double nu) {
      int i;
      if ((nu > pow(10.0, 13.47)) && (nu < pow(10.0, 16.10))) { 
        for (i=1; i <= 1298-1; i++) {
          if ((xg[i] > nu) && (xg[i+1] < nu)) {
            return yg[i];
          }
        }
      } else {
        return 0.0;
      }
      return 0.0;
}

double CosmoCalc(double z, double H0){
     // Adapted from Wright (2006)

    // initialize constants
    double   WM = 0.27;     //Omega(matter)
    double   WV = 0.73;     //Omega(vacuum) or lambda
    double   WR = 0.;        // Omega(radiation)
    double   WK = 0.;        // Omega curvaturve = 1-Omega(total)
    double   c = 299792.458; // velocity of light in km/sec
    double   DCMR = 0.0;     // comoving radial distance in units of c/H0
    double   DA = 0.0;       // angular size distance
    double   DL = 0.0;       // luminosity distance
    double   DL_Mpc = 0.0;
    double   a = 1.0;        // 1/(1+z), the scale factor of the Universe
    double   az = 0.5;       // 1/(1+z(object))
    double   h = H0/100.;
    double   adot = 0.;
    int      n=1000;         // number of points in integrals
    int      i;
    WR = 4.165E-5/(h*h);   // includes 3 massless neutrino species, T0 = 2.72528
    WK = 1-WM-WR-WV;
    az = 1.0/(1+1.0*z);
    
    for (i=1; i <= n; i++){
        a = az*(i+0.5)/n;
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
    }
    DCMR = 0.0;

    // do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
    for (i=1; i <= n; i++){
        a = az+(1-az)*(i+0.5)/n;
        adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
        DCMR = DCMR + 1./(a*adot);
    }
    DCMR = (1.-az)*DCMR/n;
   
    // tangential comoving distance
    double ratio = 1.00;
    double x = sqrt(abs(WK))*DCMR;
    if (x > 0.1){
        if (WK > 0){
            ratio =  0.5*(exp(x)-exp(-x))/x;
        }else {
            ratio = sin(x)/x;
        }
    }else {
        double y = x*x;
        if (WK < 0){
            y = -y;
        }
        ratio = 1. + y/6. + y*y/120.;
    }
    double DCMT = ratio*DCMR;
    DA = az*DCMT;
    DL = DA/(az*az);
    DL_Mpc = (c/H0)*DL;

    return DL_Mpc * 1.0e6 *pc;
}

int main(int argc, char** argv) {
   int    i, j;
   FILE   *stream;
   double freq;
   char   name[256];
   double dtmp, gf;

   // READ FILES

   if (argc > 1) {
     strcpy(InpFile, argv[1]);
     //strcpy(OutFile, argv[2]);
   }
   else {
     printf("usage: parameters.gal\n");
     return 0;
   }     

   if (load_params(InpFile)); else return 0;
   
   printf("reading file with frequency matrix ... ");

   errno  = 0;
   stream = fopen("freq.dat", "r");

   if (errno == 0) {
     j = 1;
     for (i=1; i <= 1298; i++) {
        if (j == 5) {
          fscanf(stream,"%lf\n", &x[i]);
          j = 1;
        } 
        else {
          fscanf(stream,"%lf", &x[i]);
          j++;
        }
     } 
     fclose(stream);
   } else {
     printf("file: 'freq.dat' not found !!!\n");
     return 0;
   }
     
   printf("done\n");
   printf("reading normalized galaxy spectrum for age %d [Gy] ... ", AGE);

   errno  = 0;
   sprintf(name, "./spectra/EG_%dGy.dat", AGE);
   stream = fopen(name, "r");

   if (errno == 0) {
     j = 1;
     for (i=1; i <= 1298; i++) {
        if (j == 5) {
          fscanf(stream,"%lf\n", &y[i]);
          j = 1;
        } 
        else {
          fscanf(stream,"%lf", &y[i]);
          j++;
        }
     }
     fclose(stream);
   } else {
     printf("galaxy spectrum: '%s' not found !!!\n", name);     
     return 0;
   }

   printf("done\n");

   // CALCULATE PARAMETERS

   fprintf(stderr, "\nDERIVED PARAMETERS:\n");

   //d_l      = (2.0 * c * (z + 1.0 - sqrt(z + 1.0))) / 
   //           ((H_0 * 1.0e-19) / 3.086);   
   d_l = CosmoCalc(z, H_0);

   fprintf(stderr, "luminosity distance:       %e [M pc]\n", d_l / (1.0e+6 * pc));

   // SAVE GALAXY SPECTRUM

   sprintf(name, "data/HGS_%d.dat", AGE);
   printf("\nsaving galaxy spectrum: '%s' ... ", name);

   stream = fopen(name, "w+");
   
   for (i=1; i <= 1298; i++) {
      freq  =  c/(x[i]*1.0e-8);
      xg[i] = freq;
      yg[i] = M_G * x[i] * y[i] / (4.0*M_PI*d_l*d_l*freq);
      //fprintf(stream, "%e %e %e\n", xg[i], yg[i], freq * yg[i]); //standard
      fprintf(stream, "%e %e %e %e %e\n", xg[i], yg[i], freq * yg[i], log10(freq), log10(freq * yg[i])); //OH
   }
   //to avoid plmodel bugs
   fprintf(stream, "%e %e %e %e %e\n", 3.091281e+12, 6.987882e-27, 2.160151e-13, 12.0, -20.0);
   fprintf(stream, "%e %e %e %e %e\n", 3.091281e+11, 6.987882e-27, 2.160151e-13, 11.0, -20.0);
   fprintf(stream, "%e %e %e %e %e\n", 3.091281e+10, 6.987882e-27, 2.160151e-13, 10.0, -20.0);
   fprintf(stream, "%e %e %e %e %e\n", 3.091281e+09, 6.987882e-27, 2.160151e-13, 09.0, -20.0);
   fprintf(stream, "%e %e %e %e %e\n", 3.091281e+08, 6.987882e-27, 2.160151e-13, 08.0, -20.0);
   fprintf(stream, "%e %e %e %e %e", 3.091281e+07, 6.987882e-27, 2.160151e-13, 07.0, -20.0);
   
   fclose(stream);
   
   printf("done\n");        

   printf("reading jet spectrum ... ");   
   errno  = 0;
   stream = fopen(NameSS, "r");

   i = 0;

   if (errno == 0) {
     do {   
       i++;
       if (i > NU_DIM_MAX) {
         printf("I can not load more than %d spectral points !!!", NU_DIM_MAX);
         fclose(stream);
         return 0;
       }    
       fscanf(stream,"%lf %lf %lf %lf %lf %lf\n", 
              &xs[i], &ys[i], &dtmp, &dtmp, &dtmp, &dtmp);
     } while(!feof(stream));
     fclose(stream);
   } else {
     printf("file: '%s' not found !!!\n", NameSS);
     return 0;
   }

   printf("done (%d points read)\n", i);
  
   printf("saving total spectrum into: 'data/F_gal_syn.dat' ... ");
   stream = fopen("data/F_gal_syn.dat", "w+");
   
   for (j = 1; j <= i; j++) {
      gf =  galflux(xs[j]);
      fprintf(stream,"%e %e %e\n", xs[j], gf + ys[j], xs[j] * (gf + ys[j]));
   }
   
   fclose(stream);
     
   printf("done\n");
  
   return 0;
}

int load_params(char* name) {   
    int    DEBUG = 0;
    char   stmp[256];
    double dtmp;
    FILE   *stream;
    
    errno  = 0;
    stream = fopen(name, "r");
      
    if (errno == 0) {
      
      fscanf(stream, "%s\n", stmp);
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      fscanf(stream, "%s\n", stmp);
      if (DEBUG) fprintf(stderr, "%s\n", stmp);      
      fscanf(stream, "%s\n", stmp);
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
    
      fscanf(stream, "%lf  %s\n", &dtmp, stmp);
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 0.0) && (dtmp <= 6.0)) z = dtmp; 
      else {
        fprintf(stderr, "wrong value of redshift: %e (0..6) !!!\n", dtmp); 
        return 0; 
      } 
      
      fscanf(stream, "%lf  %s\n", &dtmp, stmp);
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 50.0) && (dtmp <= 100.0)) H_0 = dtmp; 
      else {
        fprintf(stderr, "wrong value of Hubble constant: %e (50...100 [km/s/Mpc]) !!!\n", dtmp); 
        return 0; 
      }
      
      fscanf(stream, "%lf  %s\n", &dtmp, stmp);
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+6) && (dtmp <= 1.0e+13)) M_G = dtmp; 
      else {
        fprintf(stderr, "wrong value of galaxy mass: %e (1e+6...1e+13 [of Sun mass]) !!!\n", dtmp); 
        return 0; 
      }
      
      fscanf(stream, "%lf  %s\n", &dtmp, stmp);
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 0.0) && (dtmp <= 19.0)) AGE = (int)(dtmp); 
      else {
        fprintf(stderr, "wrong value of galaxy age: %e (0...19 [Gy]) !!!\n", dtmp); 
        return 0; 
      }
      
      fscanf(stream, "%s  %s\n", NameSS, stmp);
      if (DEBUG) fprintf(stderr, "%s  %s\n", NameSS, stmp);
      
      fclose(stream);
      fprintf(stderr, "PARAMETERS: '%s' LOADED !\n", name);
      
      return 1;
      
    } else {
      fprintf(stderr, "ERROR (%d) %s\n", errno, strerror(errno));
      fprintf(stderr, "PARAMETERS: '%s' NOT LOADED !!!\n", name);
      return 0;
    }
    return 0;      
}