// SSC MODEL FOR SPHERICAL HOMOGENEOUS BLOB AND A STRATIFIED JET
// AND EXTERNAL INVERSE COMPTON FROM A BLACKBODY RADIATION FIELD (HOT CORONA / ACCRETION DISK / DUST TORUS)
// AND SSC OF SECOND ORDER
// AND A SELF CONSISTENT BLOB-JET INTERACTION
//
// Nov    2018: Added: Direct EIC from the blob on the disk radiation
// August 2018: Added: IIR absorbption by Franceschini 2017
// March  2016: Added: Consistent relativistic aberration for jet slices and blob emission in the jet
//              BUT still not consistent for the entire jet, advised to stay at THETA <~= PHI (half jet aperture) and THETA <~= 1/Lorentz jet
// August 2013: Added: IIR absorbption by Franceschini 2008, an extended jet emission (synchrotron, SSC), EIC on jet
// May    2010: Added: SSC 2nd order
// July   2008: Added: IIR ABSORPTION BY KNEISKE ET AL 2004

//Author: Olivier Hervet, SSC core of the code based on the work of Katarzynski et al. 2001

//If any comments or questions, please send a mail to ohervet@ucsc.edu

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <float.h>


#include "processes_supp_v02.h"
#include "bj02.h"

using namespace std;
using namespace bj02;

namespace bj02{

int  CASE_JET;
int  CASE_X;
int  CASE_EIC;
const int EBLFLAG = 3; // 0: Kneiske et al. (2002,2004), 1: Franceschini (2008), 2: Finke (2010), 3: Franceschini (2017)//

// TIME VARIABLES

time_t       T_START, T_END;
char         PARA_FN[256];

// GLOBAL VARIABLES

int    IIR_level=0;
int    NU_DIM=50;      // CURRENT NUMBER OF SPECTRAL POINTS
double NU_STR=1.0e+10; // START FREQUENCY
double NU_END=1.0e+25; // END FREQUENCY

//double GAMMA_MIN1      = 200.0; //Minimal gamma for the jet particles
double Utot_e          = 0.0;
double Utot_B          = 0.0;
double Utot_p          = 0.0;
double I_BASE	       = 0.0;
double I_JET	       = 0.0;



// PHYSICAL PARAMETERS

double       z;
double       H_0;
double       THETA, D_L;
double       DOP_B, LOR_B, V_B, V_B_APP;//DOP_BB, 
double       R_src,R_blr;
double       L_src, L_nuc;
double       B;
double       K1;
double       K2;
double       N1;
double       N2;
double       GAMMA_MIN;
double       GAMMA_BRK;
double       GAMMA_MAX;
double       T_BB;
double       tau;  // fraction of L_nuc scattered/rerocessed isotropically (EIC)
double       B_0;
double       N_0;
double       n_n;
double       n_N;
double       n_G;
double       D_b, D_b_src_j, D_BJ;
double       U_B,U_e;
double       jj1,kk1;
double       L_jet_eic;
double       L_tor, T_BB_tor;
double       Tcool_synch, Tcool_vhe, Tcool_radio;

// TRANSFORMATION PARAMETERS

double DOP_JET, LOR_JET, V_JET, V_JET_APP, DOP_B_J, LOR_B_J, V_B_J;

// GEOMETRY PARAMETERS


int          SL_DIM;
double       PHI, PHI_src;
double       THETA_src_j;
double       J_LEN, J_LEN_src;
double       A;
double       X_MIN;
double       X_MAX;
double       X_OUT;
double       Y_MIN;
double       Y_MAX;	
double       R_MOY;
double       Xcut;
double       antisym;
double       Sprev;
double       Sum_S;

//const int           SL_DIM_MAX = 200;//

// OTHER PARAMETERS

int          SL_CUR;   
double       GAMMA_MIN1;
double       GAMMA_MAX_0;
int          null0; 
double       V_exp, DEL_Tph, DEL_Tj;
char  prefix[256];


// VECTORS

double 	     NU[NU_MAX+1];
double       I_rad[NU_MAX+1];
double       I_rad1st[NU_MAX+1]; // for 2nd order SSC
double       I_CMB[NU_MAX+1];

double       X_VAL[SL_DIM_MAX+2];
double       Y_VAL[SL_DIM_MAX+2];
double       DEL_X[SL_DIM_MAX+2];
double       Sum_DEL[SL_DIM_MAX+2];
double       Stot[SL_DIM_MAX];
double       Sr_TOT[SL_DIM_MAX];

double       B_VAL[SL_DIM_MAX+1];
double       N_VAL[SL_DIM_MAX+1];
double       G_VAL[SL_DIM_MAX+1];
double       UJ_SLICE[SL_DIM_MAX+1];
double       PJ_SLICE[SL_DIM_MAX+1];

double       I_SYN_EXT[NU_DIM_MAX+1];
double       L_BB_nuc[NU_DIM_MAX+1];
double       I_eic_jet[NU_DIM_MAX+1];

double       I_rad_syn[NU_MAX+1];
double       I_rad_com[NU_MAX+1];
double       I_rad2nd[NU_MAX+1];
double       I_rad_ext[NU_MAX+1];
double       I_rad_ext1[NU_MAX+1];
double       I_com_ext[NU_MAX+1];
double       I_com_ext1[NU_MAX+1];
double       I_com_disc[NU_MAX+1];
double       I_rad_ext_s[NU_MAX+1];

double       C_e[NU_MAX+1];
double       C_a[NU_MAX+1];
double       C_e1[NU_MAX+1];
double       C_a1[NU_MAX+1];

double       F_IC_disk[NU_MAX+1];
double       F_IC_tot[NU_MAX+1];
double       NU_IC_disk[NU_MAX+1];




// MATRIXES

double       I_SYN_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       I_SYN_JET_EDGE[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       I_SYN_JET_BASE[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       J_SYN_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       K_ESA_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       F_SYN_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];

double       I_COM_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       I_COM_JET_EDGE[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       I_COM_JET_BASE[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       J_COM_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       K_ABS_SSC_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       F_COM_JET[SL_DIM_MAX+1][NU_DIM_MAX+1];

double       I_SYN_TOT[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       I_COM_TOT[SL_DIM_MAX+1][NU_DIM_MAX+1];
double       F_SYN_BLOB[SL_DIM_MAX+1][NU_DIM_MAX+1];

double       Sl[SL_DIM_MAX+1][SL_DIM_MAX+1];
double       Sr[SL_DIM_MAX+1][SL_DIM_MAX+1];
double       S[SL_DIM_MAX+1][SL_DIM_MAX+1];
double       Dcut[SL_DIM_MAX+1][SL_DIM_MAX+1];



// SUBROUTINES ******************************************************************

int ifexist(char *fname) {
  FILE* stream;
  errno = 0;
  stream = fopen(fname, "r");
  if (errno == 0) {
    fclose(stream);
    return 1;
  }
  else {
    return 0;
  }
}

// ELECTRON SPECTRUM
double N_e_BknPowLaw(double gamma) { // Broken power-law (possibly with exp cutoff at high energies)
  K2 = K1  * pow(GAMMA_BRK, N2-N1);
  
  if (gamma < GAMMA_BRK) {
    return K1 * pow(gamma, -N1);
  }
  if (gamma >= GAMMA_BRK) {
    return K2 * pow(gamma, -N2);// * exp(-gamma/GAMMA_MAX);
  } 
  
  return 1.0e-100;
}

double N_e_SmoothBknPowLaw(double gamma) { //cf. Tavecchio et al. 2001: 2001ApJ...554..725T, + exp cutoff at high energies
  return(K1 * pow(gamma, -N1) * pow((1.+gamma/GAMMA_BRK), (N1-N2))) * exp(-gamma/GAMMA_MAX);
}

double N_e_PileUp(double gamma){
  return(K1 * gamma*gamma * exp(-2.0 * gamma/GAMMA_BRK));
}

// Switch between different electron spectrum shapes
double N_e(double gamma){
  double foo=0.;
  //foo=N_e_PileUp(gamma);
  foo=N_e_SmoothBknPowLaw(gamma);
  //foo=N_e_BknPowLaw(gamma);
  return foo;
}


double N_e_Jet(double gg) {
      if ((gg >= GAMMA_MIN1) && (G_VAL[SL_CUR])) return N_VAL[SL_CUR] * pow(gg, -n_n);     
      return 0.0;
}

double ftr(double F) {
      return F * pow(DOP_JET, 3) * (1.0 + z);
}

//analytical integration over the electron spectrum

// simple power law (jet)
double int_spec_j(double N_0){
   if (n_n == 2){
      return N_0*(log(GAMMA_MAX_0)-log(GAMMA_MIN1));
   }
   else{           
      return N_0/(2-n_n)*(pow(GAMMA_MAX_0,2-n_n)-pow(GAMMA_MIN1,2-n_n));
   }
}

  // power law with exponential cut-off (blob)
double int_spec_b(){
   if (N1 == 2 and N2 == 2){
      return K1*((log(GAMMA_BRK)-log(GAMMA_MIN)) + pow(GAMMA_BRK, N2-N1)*(log(GAMMA_MAX)-log(GAMMA_BRK)));
   }
   if (N1 == 2 and N2 != 2){
      return K1*((log(GAMMA_BRK)-log(GAMMA_MIN)) + pow(GAMMA_BRK,(N2-N1))/(2-N2)*(pow(GAMMA_MAX,(2-N2))-pow(GAMMA_BRK,2-N2)));
   }
   if (N1 != 2 and N2 == 2){
      return K1*(1/(2-N1)*(pow(GAMMA_BRK,2-N1)-pow(GAMMA_MIN,2-N1)) + pow(GAMMA_BRK,N2-N1)*(log(GAMMA_MAX)-log(GAMMA_BRK)));
   }
   else{           
      return K1*(1/(2-N1)*(pow(GAMMA_BRK,2-N1)-pow(GAMMA_MIN,2-N1)) + pow(GAMMA_BRK,(N2-N1))/(2-N2)*(pow(GAMMA_MAX,2-N2)-pow(GAMMA_BRK,2-N2)));
   }
}

}



void description(){
  cout << "This code computes the radiative output from a homogeneous blob and a stratified jet given ther particle energy distribution" << endl
       << "The processes accounted for are:" << endl
       << endl
       << " - synchrotron radiation (blob and jet)" << endl
       << " - SSC radiation (blob and jet)" << endl
       << " - 2nd order SSC radiation (blob)" << endl
       << " - external inverse Compton on a disk/corona radiation field" << endl
       << " - external inverse Compton on the CMB" << endl
       << " - external inverse Compton on the jet synchrotoron radiation field" << endl
       << endl;
}




int load_params(char* name);



int main(int argc, char** argv) {
   int    i, j, k, l;
   char   stmp[256];
   char   comma[256];
   char   name[256];
   char   stmp1[256], stmp2[256];
   double tmp_min, tmp_max, tmp_stp, tmp_val, tmp_cur;
   double nu_tmp, fx_tmp, nu_tmp1, fx_tmp1;
   double jj, kk, tt, I_syn, I_SYN, I_com, I_com2, I_COM; 
   double dtmp, dtmp1, dtmp2, dtmp3;;
   double Ub_e, Uj_e, Ub_syn, Ub_ssc, Ub_ssc2, Ub_eicd, Ub_eicj, Ub_B;
   double Lb_B, Lj_B, Lb_e, Lj_e, Lb_r, Lj_r;
   double Pj_syn, Pj_ssc, elec_spec_tot_b, elec_spec_tot_j;
   FILE*  stream_dat;
   FILE*  stream_dat1;
   FILE*  stream_tau;

   T_START = time(NULL);


   cout << endl;
   cout << "SSC MODEL FOR SPHERICAL HOMOGENEOUS BLOB" << endl
	<< "EXTERNAL INVERSE COMPTON ON THE NUCLEUS RADIATION FIELD (HOT CORONA AND/OR ACCRETION DISK)" << endl
	<< "EXTERNAL INVERSE COMPTON ON CMB" << endl 
	<< "SSC MODEL FOR A STATIFIED JET" << endl 
	<< "EXTERNAL INVERSE COMPTON ON THE STATIFIED JET SYNCHROTRON RADIATION FIELD" << endl<< endl;

   sprintf(PARA_FN, "%s", "none");
   

   // LOAD PARAMETERS

   if (argc == 2) {
     if(strcmp(argv[1],"--help")==0){
       description();
       return 0;
     }
     strcpy(PARA_FN, argv[1]);
   }
   else {
     cout << "usage: bj parameters.par" << endl << endl;
     return 0;
   }    
      
   if (load_params(PARA_FN)); else return 0;


   double Gamma[G_DIM+1];
   double elec_spec[G_DIM+1];


  
   D_L       = (2.0 * c * (z + 1.0 - sqrt(z + 1.0))) / 
     ((H_0 * 1.0e-19) / 3.086);
   


   V_B       = c * (DOP_B*DOP_B * cos(THETA * M_PI / 180.0) - 
		    sqrt(1.0 - DOP_B*DOP_B * sin(THETA * M_PI / 180.0) * sin(THETA * M_PI / 180.0)))
     / (1.0 + DOP_B*DOP_B * cos(THETA * M_PI / 180.0) * cos(THETA * M_PI / 180.0));
   
   LOR_B   = 1.0 / sqrt(1.0 - pow(V_B / c, 2.0));	 	     

   //only for NGC1275
   /*
   LOR_B   = 9.86; //Doppler unbeamed case (too large angle)
   V_B = sqrt(1. - 1./(LOR_B*LOR_B)) *c;
   */
   
   
   
   if( DOP_B * sin((THETA)*M_PI/180.0)> 1.0){
     fprintf(stderr, "we must have DOP_B * sin(THETA) < 1, else the speed of the blob is NOT REAL...\n\nProgram aborted...");
     return 0;
   }
  
   V_B_APP = ((V_B) * sin(THETA * M_PI / 180)) / 
     (1.0 - (V_B/c)*cos(THETA * M_PI / 180));
   
   Ub_B = B*B/(8*M_PI);
   Lb_B = M_PI*R_src*R_src*LOR_B*LOR_B*c*Ub_B;
     
 
   if(CASE_JET){  
     
      V_JET    = c * (2.0 * DOP_JET*DOP_JET * cos(THETA * M_PI / 180.0) - 2.0 * 
              sqrt(1.0 - DOP_JET*DOP_JET + DOP_JET*DOP_JET*cos(THETA * M_PI / 180.0)*
              cos(THETA * M_PI / 180.0))) / (2.0*(1.0 + DOP_JET*DOP_JET * 
              cos(THETA * M_PI / 180.0) * cos(THETA * M_PI / 180.0)));              	 	                                 

      LOR_JET  = 1.0 / sqrt(1.0 - pow(V_JET/c, 2.0));

      V_JET_APP = ((V_JET) * sin(THETA * M_PI / 180.0)) / 
               (1.0 - (V_JET/c)*cos(THETA * M_PI / 180.0));
      
      J_LEN_src = J_LEN / LOR_JET;
      PHI_src =  atan(LOR_JET* tan(PHI*M_PI/180.0)) * 180.0/M_PI;
      D_b_src_j = D_b / LOR_JET;
	       
      //relativistic aberration [deg]
      THETA_src_j = acos((cos(THETA * M_PI / 180.0) - V_JET/c) / (1.0 - V_JET/c * cos(THETA * M_PI / 180.0))) * 180.0/M_PI;
      
      //jet lateral expansion speed (jet comoving frame)
      V_exp = V_JET * tan(PHI_src*M_PI/180.0);
      

      // transformation blob to jet frame
      V_B_J = (V_B - V_JET)/(1.0 - (V_B*V_JET/(c*c)));
      LOR_B_J = 1.0 / (sqrt(1.0 - pow(V_B_J/c,2.0)));
      DOP_B_J = 1.0 / (LOR_B_J * (1.0 - V_B_J/c * cos(THETA * M_PI / 180.0)));
      
      // distance blob - first slice (jet frame)
      D_BJ = (Y_MIN / tan(PHI_src * M_PI / 180.0)) - D_b_src_j;
   }
   //
   if (PRINT) {                
     fprintf(stderr, "Hubble constant:        %6.3f\n", H_0);
     fprintf(stderr, "redshift                %6.3f\n", z);
     fprintf(stderr, "\nBlob parameters:\n");
     fprintf(stderr, "---------------\n");
     fprintf(stderr, "Doppler factor:         %6.3f\n", DOP_B);
     fprintf(stderr, "theta:                  %6.3e\n", THETA);
     fprintf(stderr, "V blob:                 %6.3e\n", V_B/c);
     fprintf(stderr, "Lorentz factor:         %6.3e\n", LOR_B);            
     fprintf(stderr, "V apparent:             %6.3e\n", V_B_APP/c); 
     fprintf(stderr, "K_1:                    %6.3e\n", K1);
     fprintf(stderr, "n_1:                    %6.3f\n", N1);
     fprintf(stderr, "n_2:                    %6.3f\n", N2);
     fprintf(stderr, "gamma_min:              %6.3e\n", GAMMA_MIN);
     fprintf(stderr, "gamma_brk:              %6.3e\n", GAMMA_BRK);
     fprintf(stderr, "gamma_max:              %6.3e\n", GAMMA_MAX);  
     fprintf(stderr, "radius:              %6.3e\n", R_src); 
     if(CASE_EIC){
       fprintf(stderr, "\nEIC parameters:\n");
       fprintf(stderr, "---------------\n");
       fprintf(stderr, "disk blackbody temp.:   %6.3e\n",T_BB);
       fprintf(stderr, "tore blackbody temp.:   %6.3e\n",T_BB_tor);
       fprintf(stderr, "nuclear lumin.:         %6.3e\n",L_nuc);
       fprintf(stderr, "reprocessing fraction:  %6.3e\n",tau);
       fprintf(stderr, "R_blr:                  %6.3e\n",R_blr);
     }
     if(CASE_JET){
       fprintf(stderr, "\nJET parameters:\n");
       fprintf(stderr, "---------------\n");
       fprintf(stderr, "Doppler factor:          %6.3f\n", DOP_JET);
       fprintf(stderr, "V jet:                   %6.3e\n", V_JET/c);
       fprintf(stderr, "Lorentz factor:          %6.3e\n", LOR_JET);
       fprintf(stderr, "V apparent:              %6.3e\n", V_JET_APP/c);
       fprintf(stderr, "Doppler factor blob-jet: %6.3e\n", DOP_B_J);
       fprintf(stderr, "V blob-jet:              %6.3e\n", V_B_J/c);
       fprintf(stderr, "Lorentz factor blob-jet: %6.3e\n", LOR_B_J);
       fprintf(stderr, "K_0:                     %6.3e\n", N_0);
       fprintf(stderr, "n:                       %6.3e\n", n_n);
       fprintf(stderr, "n_K:                     %6.3e\n", n_N);
       fprintf(stderr, "gamma_max_0:             %6.3e\n", GAMMA_MAX_0);
       fprintf(stderr, "n_gamma_max:             %6.3e\n", n_G);
       fprintf(stderr, "B_0:                     %6.3e\n", B_0);
       fprintf(stderr, "n_B:                     %6.3e\n", n_B);
       fprintf(stderr, "R_0:                     %6.3e\n", Y_MIN);
       fprintf(stderr, "theta (comob frame):              %6.3e\n", THETA_src_j);
       fprintf(stderr, "jet_length (comov frame):         %6.3e\n", J_LEN_src);
       fprintf(stderr, "jet_length (host frame):           %6.3e\n", J_LEN);
       fprintf(stderr, "half-opening_angle (comov frame): %6.3e\n", PHI_src);
       fprintf(stderr, "half-opening_angle (host frame):   %6.3e\n", PHI);
       fprintf(stderr, "Distance SMBH-first slice (comov frame):   %e [cm] (%f [pc])\n", Y_MIN / tan(PHI_src * M_PI / 180.0), Y_MIN / tan(PHI_src * M_PI / 180.0)/pc); 
       fprintf(stderr, "Distance SMBH-first slice (host frame):    %e [cm] (%f [pc])\n", Y_MIN / tan(PHI * M_PI / 180.0), Y_MIN / tan(PHI * M_PI / 180.0)/pc);
       fprintf(stderr, "Distance SMBH-blob (jet comov frame):      %e [cm] (%f [pc])\n", D_b_src_j, D_b_src_j/pc);
       fprintf(stderr, "Distance SMBH-blob (host frame):           %e [cm] (%f [pc])\n", D_b, D_b/pc);
       fprintf(stderr, "nb_slices:                        %6.3d\n\n", SL_DIM);
       
       /*
       if( Y_MIN / tan(PHI_src * M_PI / 180.0) > D_b_src_j){
         fprintf(stderr, "we must have R_jet / tan(PHI_src) < D_blob_src_j, else the blob is NOT INSIDE THE JET...\n\nProgram aborted...\n\n");
         return 0;
       }*/
     }
     else cout << endl;
   }     
   



   

   // SAVING ELECTRONS ENERGY SPECTRUM
   
   sprintf(stmp, "./data/%s_es.dat", prefix);
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/%s_es.dat ./data/%s_prev_es.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
   stream_dat = fopen(stmp, "w+");   

   tmp_min   = log10(GAMMA_MIN);
   tmp_max   = log10(GAMMA_MAX);
   tmp_stp   = (tmp_max - tmp_min) / (G_DIM - 1);
   tmp_val   = tmp_min;
   
   for (i = 1; i <= G_DIM ; i++) {
      tmp_cur = pow(10.0, tmp_val);
      dtmp    = N_e(tmp_cur);
      Gamma[i] = tmp_cur;
      elec_spec[i] = dtmp*tmp_cur;
      if (dtmp > 1e-300) {
        fprintf(stream_dat, "%f %f %e %e\n", 
                log10(tmp_cur), log10(dtmp), tmp_cur, dtmp);  
      }                        
      tmp_val = tmp_val + tmp_stp;
   }   
   
   fclose(stream_dat);  
   
   elec_spec_tot_b = int_spec_b();
   Ub_e = m_e*c*c*elec_spec_tot_b;
   Lb_e = M_PI*R_src*R_src*LOR_B*LOR_B*c*Ub_e;




   // PREPARING FREQUENCY VECTOR
   
   tmp_min   = log10(FreqTransO2S(NU_STR, DOP_B, z));
   tmp_max   = log10(FreqTransO2S(NU_END, DOP_B, z));
   tmp_stp   = (tmp_max - tmp_min) / (NU_DIM);
   tmp_val   = tmp_min;
   
   
   for (i = 1; i <= NU_DIM; i++) {
      tmp_cur   = pow(10.0, tmp_val);
      NU[i] = tmp_cur;
      tmp_val   = tmp_val + tmp_stp; 
   }
   //Warning:: NU[1] < FreqTransO2S(NU_STR, DOP_B, z) !!
   NU[1] = FreqTransO2S(NU_STR, DOP_B, z);











   fprintf(stderr, "CALCULATING SYNCHROTRON SPECTRUM ... ");
  
  
   // CALCULATING SYNCHROTRON SPECTRUM
   if (CASE_JET == 0){
      sprintf(stmp, "./data/%s_ss.dat", prefix);
      if (ifexist(stmp)) {
        sprintf(comma, "mv ./data/%s_ss.dat ./data/%s_prev_ss.dat", prefix, prefix);
         if(system(comma)) return 1;
      }
      stream_dat = fopen(stmp, "w+");
   }


   for (i = 1; i <= NU_DIM; i++) {
      if (PRINT) fprintf(stderr, "%3i", i);
	 
      // emission & absorption coefficients
      jj    = j_syn(N_e, GAMMA_MIN, GAMMA_MAX, NU[i], B, SYN_PREC1, SYN_PREC2);
      kk    = k_esa(N_e, GAMMA_MIN, GAMMA_MAX, NU[i], B, ABS_PREC1, ABS_PREC2);            
      // Here the computation is made in the source frame
      


      // transfer equation//
      if (L_src == 0.0 ) {
        I_syn    = SphTransfEquat(jj, kk, R_src);
	// K. Katarzynski:
        // THIS IS TRICK !!! 
        // solution of the transfer for central point of blob is
        // the same like for cylindrical geometry      
        //I_rad[i] = 0.75 * CylTransfEquat(0.0, jj, kk, R_src);
        //OH, the coefficient 0.70 seems more accurate
        I_rad[i] = 0.70 * CylTransfEquat(0.0, jj, kk, R_src);
	I_rad_syn[i] = I_syn;
      } else {
        I_syn    = CylTransfEquat(0.0, jj, kk, L_src);
        I_rad[i] = I_syn;
      }
      
  
    
      // transformation to observer frame
      nu_tmp = FreqTransS2O(NU[i], DOP_B, z);
      fx_tmp = Intens2Flux(I_syn, R_src, DOP_B, z, H_0);
      
      if (CASE_JET == 0){
         if (fx_tmp > 1.0e-300) {
           fprintf(stream_dat, "%f %f %f\n", 
                log10(nu_tmp), 
                log10(fx_tmp), 
                log10(nu_tmp*fx_tmp)
               );
         } 
      }
         	   
      if (PRINT) fprintf(stderr, "\b\b\b");        
   }   
   
   if (CASE_JET == 0) fclose(stream_dat);
   
   fprintf(stderr, "DONE\n\n");
   
   //synch energy density
   Ub_syn = 4*M_PI/c*Simpson(I_rad, NU, NU_DIM, 1, NU_DIM);
   
   
   




       
   
   fprintf(stderr, "CALCULATING SSC SPECTRUM ... ");   
      
   // CALCULATING INVERSE-COMPTON SPECTRUM
   if (CASE_JET == 0){
      sprintf(stmp, "./data/%s_cs.dat", prefix);
      if (ifexist(stmp)) {
         sprintf(comma, "mv ./data/%s_cs.dat ./data/%s_prev_cs.dat", prefix, prefix);
         if(system(comma)) return 1;
      }
      stream_dat = fopen(stmp, "w+");

      sprintf(stmp, "./data/%s_tau.dat", prefix);
      if (ifexist(stmp)) {
         sprintf(comma, "mv ./data/%s_tau.dat ./data/%s_prev_tau.dat", prefix, prefix);
         if(system(comma)) return 1;
      }
      stream_tau = fopen(stmp, "w+");
   }
     
   for (i = 1; i <= NU_DIM; i++) {      
      if (PRINT) fprintf(stderr, "%3i", i);

      // emission & absorption coefficients
      jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad, 
                 FreqTransO2S(NU_STR, DOP_B, z), 
		 FreqTransO2S(NU_END, DOP_B, z), 
		 NU_DIM, NU[i], COM_PREC1, COM_PREC2); 

      kk = gg_abs(NU[i], I_rad, NU_DIM, 
		      FreqTransO2S(NU_STR, DOP_B, z), 
		      FreqTransO2S(NU_END, DOP_B, z), SYN_PREC1, SYN_PREC2);

      // absorption by pair production 
      if (L_src == 0.0) {
        I_com = SphTransfEquat(jj, kk, R_src);
	I_rad1st[i] = 0.70 * CylTransfEquat(0.0, jj, kk, R_src);
	I_rad_com[i] = I_com;
      } else {
        I_com = CylTransfEquat(0.0, jj, kk, L_src);
	I_rad1st[i] = I_com;
	I_rad_com[i] = I_com;
      }  
      
 
      if (CASE_JET == 0){
         // frequency transformation to observer frame
         nu_tmp = FreqTransS2O(NU[i], DOP_B, z);	 
      
         // optical depth for absorption of VHE gamma rays by IIR
         //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
        if(IIR_level==1){

          if(EBLFLAG==0){
             if (z <= 0.1){
	        tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	     }else{
	        tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	     }
           }else if(EBLFLAG==1){
             tt=tau_IRA_Franceschini(nu_tmp, z);
	   }else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
           }else if(EBLFLAG==3){
             tt=tau_IRA_Franceschini17(nu_tmp, z);
	   }

         }else tt=0.;

         fprintf(stream_tau,"%f %f\n",
	      nu_tmp*h*0.62415,  // Conversion Hz -> TeV
	      exp(-tt)
	      );


         // flux transformation to observer frame
         fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);
         dtmp   = fx_tmp;
      
         // absorption by IIR
         fx_tmp = fx_tmp * exp(-tt); 
         //fprintf(stderr, "log10 freq=%f Tau SSC EBL = %f \n",log10(nu_tmp),tt);
         
         //fprintf(stderr, "log10 freq=%f Tau SSC = %f \n",log10(nu_tmp),2.0 * R_src * kk);
         
      
         if (fx_tmp > 1.0e-300) {
           fprintf(stream_dat, "%f %f %f %f %f\n", 
                log10(nu_tmp), 
                log10(fx_tmp), 
                log10(nu_tmp*fx_tmp),
                log10(dtmp),
                log10(nu_tmp*dtmp) 
		);
         }
      }
      
      if (PRINT) fprintf(stderr, "\b\b\b");          
   }   
   
   if (CASE_JET == 0){
     fclose(stream_dat);
     fclose(stream_tau);
   }
  
   fprintf(stderr, "DONE\n\n");
   
   Ub_ssc = 4*M_PI/c*Simpson(I_rad_com, NU, NU_DIM, 1, NU_DIM);
   













   fprintf(stderr, "CALCULATING 2nd ORDER SSC SPECTRUM ... ");


   // CALCULATING 2nd ORDER INVERSE-COMPTON SPECTRUM
   
   if (CASE_JET == 0){
      sprintf(stmp, "./data/%s_cs2.dat", prefix);
      if (ifexist(stmp)) {
         sprintf(comma, "mv ./data/%s_cs2.dat ./data/%s_prev_cs2.dat", prefix, prefix);
         if(system(comma)) return 1;
      }
      stream_dat = fopen(stmp, "w+");
   
      sprintf(stmp, "./data/%s_tau.dat", prefix);
      if (ifexist(stmp)) {
         sprintf(comma, "mv ./data/%s_tau.dat ./data/%s_prev_tau.dat", prefix, prefix);
         if(system(comma)) return 1;
      }
      stream_tau = fopen(stmp, "w+");
   }
   
   sprintf(stmp1, "./data/%s_cs2_blob_frame.dat", prefix);
   stream_dat1 = fopen(stmp1, "w+");

   
   for (i = 1; i <= NU_DIM; i++) {      
      if (PRINT) fprintf(stderr, "%3i", i);

      // emission & absorption coefficients
      jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad1st, 
                 FreqTransO2S(NU_STR, DOP_B, z), 
		 FreqTransO2S(NU_END, DOP_B, z), 
		 NU_DIM, NU[i], COM_PREC1, COM_PREC2); 

      kk = gg_abs(NU[i], I_rad1st, NU_DIM, 
		      FreqTransO2S(NU_STR, DOP_B, z), 
		      FreqTransO2S(NU_END, DOP_B, z), COM_PREC1, COM_PREC2);
      //computed in the source frame


      // absorption by pair production 
      if (L_src == 0.0) {
        I_com2 = SphTransfEquat(jj, kk, R_src);
      } else {
        I_com2 = CylTransfEquat(0.0, jj, kk, L_src);
      }
      I_rad2nd[i] = I_com2;
      
      
      if (CASE_JET == 0){
         // frequency transformation to observer frame
         nu_tmp = FreqTransS2O(NU[i], DOP_B, z);

        if(IIR_level==1){

          if(EBLFLAG==0){
             if (z <= 0.1){
	        tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	     }else{
	        tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	     }
           }else if(EBLFLAG==1){
             tt=tau_IRA_Franceschini(nu_tmp, z);
	   }else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
           }else if(EBLFLAG==3){
             tt=tau_IRA_Franceschini17(nu_tmp, z);
	   }

         }else tt=0.;
      
         fprintf(stream_tau,"%f %f\n",
	      nu_tmp*h*0.62415,  // Conversion Hz -> TeV
	      exp(-tt)
	      );
      

         // flux transformation to observer frame
         fx_tmp = Intens2Flux(I_com2, R_src, DOP_B, z, H_0);
         dtmp   = fx_tmp;
      
         // absorption by IIR
         fx_tmp = fx_tmp * exp(-tt); 

         if (fx_tmp > 1.0e-300) {
           fprintf(stream_dat, "%f %f %f %f %f\n", 
                log10(nu_tmp), 
                log10(fx_tmp), 
                log10(nu_tmp*fx_tmp),
                log10(dtmp),
                log10(nu_tmp*dtmp) 
		);
         }
      }
      
      nu_tmp1 = FreqTransS2O(NU[i], 1, z);
      fx_tmp1 = Intens2Flux(I_com2, R_src, 1, z, H_0);
      if (fx_tmp1 > 1.0e-300) {
        fprintf(stream_dat1, "%f %f %f\n", 
             log10(nu_tmp1), 
             log10(fx_tmp1), 
             log10(nu_tmp1*fx_tmp1)
	      );
      }
      
      if (PRINT) fprintf(stderr, "\b\b\b");          
   }   
   
   if (CASE_JET == 0) fclose(stream_dat);
   fclose(stream_dat1);
  
   fprintf(stderr, "DONE\n\n");
   
   Ub_ssc2 = 4*M_PI/c*Simpson(I_rad2nd, NU, NU_DIM, 1, NU_DIM);





   // NUCLEAR LUMINOSITY

  sprintf(stmp, "./data/%s_nuc.dat", prefix);
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/%s_nuc.dat ./data/%s_prev_nuc.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
  
  sprintf(stmp1, "./data/%s_tor.dat", prefix);
   if (ifexist(stmp1)) {
      sprintf(comma, "mv ./data/%s_tor.dat ./data/%s_prev_tor.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
   stream_dat1 = fopen(stmp1, "w+");


   if(CASE_EIC){
   fprintf(stderr, "CALCULATING NUCLEAR LUMINOSITY SPECTRUM ... ");   
   stream_dat = fopen(stmp, "w+");



     double L_BB_disk;
     double T_BB_nuc=T_BB;
     double L_BB_tor;
     double L_x = 0;
     double L_BB_disk1;

     for (i = 1; i <= NU_DIM; i++) {
       if (PRINT) fprintf(stderr, "%3i", i);
	 
       L_BB_disk1 = L_BB_disk;
       
       L_BB_disk = L_nuc * Planck(NU[i],T_BB_nuc) / (sig/M_PI * pow(T_BB_nuc,4.0));
       
       // X corona (Ghisellini 2009)
       if (CASE_X && L_BB_disk <= L_BB_disk1) {
         //L_x = 0.5 * L_nuc * pow(NU[i],-1.) * exp(-NU[i]/3.0e+17);
         L_x = 0.2 * L_nuc * pow(NU[i],-1.0) * exp(-NU[i]/3.628e+19); //OJ 287
         //original
        //L_x = 0.3 * L_nuc * pow(NU[i],-1) * exp(-NU[i]/3.628e+19);
       }
       
       L_BB_tor = L_tor * Planck(NU[i],T_BB_tor) / (sig/M_PI * pow(T_BB_tor,4.0));
     
       // transformation to observer frame
       nu_tmp = FreqTransS2O(NU[i], 1.0, z);
       if (L_BB_disk >= L_x){
           L_BB_nuc[i] = L_BB_disk;
       }else{
            L_BB_nuc[i] = L_x;
       }
       //fprintf(stderr, "%6.3e\n", L_BB_nuc[i]);
       
       fx_tmp = (1.0+z) * L_BB_nuc[i] / (4.0 * M_PI * D_L * D_L);
       if (fx_tmp > 1.0e-300) {
	 fprintf(stream_dat, "%f %f %f\n", 
		 log10(nu_tmp), 
		 log10(fx_tmp), 
		 log10(nu_tmp*fx_tmp)
		 );
       }
       
       nu_tmp1 = FreqTransS2O(NU[i], 1.0, z);
       fx_tmp1 = (1.0+z) * L_BB_tor / (4.0 * M_PI * D_L * D_L);
      
       if (fx_tmp1 > 1.0e-300) {
	 fprintf(stream_dat1, "%f %f %f\n", 
		 log10(nu_tmp1), 
		 log10(fx_tmp1), 
		 log10(nu_tmp1*fx_tmp1)
		 ); 
       }
       if (PRINT) fprintf(stderr, "\b\b\b");        
     }   
   
   fclose(stream_dat1); 
   
   fprintf(stderr, "DONE\n\n");
   }
   




   sprintf(stmp, "./data/%s_ecdisc.dat", prefix);
   sprintf(stmp1, "./data/%s_prev_ecdisc.dat", prefix);
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/%s_ecdisc.dat ./data/%s_prev_ecdisc.dat", prefix, prefix);
      if(system(comma)) return 1;   
   } else if (ifexist(stmp1)){
      sprintf(comma, "rm ./data/%s_prev_ecdisc.dat", prefix);
      if(system(comma)) return 1;
   }
   

   if (CASE_EIC){

       fprintf(stderr, "CALCULATING EXT. INV. COMPTON SPECTRUM ON NUCLEAR RADIATION ... ");
       // Not finished yet, the code actually doesn't consider it
       // need to create a specific radiation file for it (merging with the IC BLR not a good idea)

       stream_dat = fopen(stmp, "w+");

       
       for (i = 1; i <= NU_DIM; i++) {
            if (PRINT) fprintf(stderr, "%3i", i);

            // radiation field from the disk in the blob frame
            I_rad_ext[i] = L_BB_nuc[i] / (4.0 * M_PI * D_b*D_b * 16.0*pow(LOR_B,4.0));  //Rybicky & lightman p.142 
            //fprintf(stderr, "%6.3e\n", I_rad_ext[i]);

            // emission & absorption coefficients
            
            //How the blob sees the direct disk radiation (redshift)
            jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad_ext, 
                        FreqTransO2S(NU_STR, DOP_B, z)/LOR_B,
                        FreqTransO2S(NU_END, DOP_B, z)/LOR_B,
                        NU_DIM, NU[i]/LOR_B, COM_PREC1, COM_PREC2); 

            // absorption by pair production on direct disc, not sure if it should be taken into account
            //pure anisotropic backward radiation, already abs on BLR
            /*
            kk = gg_abs(NU[i]/LOR_B, I_rad_ext, NU_DIM, 
                            FreqTransO2S(NU_STR, DOP_B, z)/LOR_B, 
                            FreqTransO2S(NU_END, DOP_B, z)/LOR_B, COM_PREC1, COM_PREC2);    
            */
            
            // absorption by pair production (on synch field)
            kk1 = gg_abs(NU[i], I_rad, NU_DIM, 
                        FreqTransO2S(NU_STR, DOP_B, z), 
                        FreqTransO2S(NU_END, DOP_B, z), SYN_PREC1, SYN_PREC2);
            
            if (L_src == 0.0) {
            I_com = SphTransfEquat(jj, kk1, R_src);
            } else {
            I_com = CylTransfEquat(0.0, jj, kk1, L_src);
            } 
            I_com_disc[i] = I_com;

            // frequency transformation to observer frame
            nu_tmp = FreqTransS2O(NU[i]/LOR_B, DOP_B, z);
            NU_IC_disk[i] = nu_tmp;


            // optical depth for absorption of VHE gamma rays by IIR
            if(IIR_level==1){

            if(EBLFLAG==0){
                if (z <= 0.1){
                    tt = tau_IRA_Kneiske(nu_tmp, z, 0);
                }else{
                    tt = tau_IRA_Kneiske(nu_tmp, z, 1);
                }
                }else if(EBLFLAG==1){
                    tt=tau_IRA_Franceschini(nu_tmp, z);
                }else if(EBLFLAG==2){
                    tt=tau_IRA_Finke(nu_tmp, z);
                }else if(EBLFLAG==3){
                    tt=tau_IRA_Franceschini17(nu_tmp, z);
                }
            }else tt=0.;
            
            
            // flux transformation to observer frame
            fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);
            //fx_tmp = 0.;
            //fprintf(stderr, "%6.3e\n", fx_tmp);
            dtmp   = fx_tmp;
            F_IC_disk[i] = dtmp;


            // absorption by IIR
            fx_tmp = fx_tmp * exp(-tt); 
            
            if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", 
                    log10(nu_tmp), 
                    log10(fx_tmp), 
                    log10(nu_tmp*fx_tmp),
                    log10(dtmp),
                    log10(nu_tmp*dtmp) 
                    );
            }

            if (PRINT) fprintf(stderr, "\b\b\b");  

       }
   
       fclose(stream_dat);

       fprintf(stderr, "DONE\n\n");

    if (CASE_JET == 0){
        sprintf(stmp, "./data/%s_ecs.dat", prefix);
        sprintf(stmp1, "./data/%s_prev_ecs.dat", prefix);
        if (ifexist(stmp)) {
            sprintf(comma, "mv ./data/%s_ecs.dat ./data/%s_prev_ecs.dat", prefix, prefix);
            if(system(comma)) return 1;   
        } else if (ifexist(stmp1)){
            sprintf(comma, "rm ./data/%s_prev_ecs.dat", prefix);
            if(system(comma)) return 1;
        }
    }
   
    if (D_b <= R_blr){

    fprintf(stderr, "CALCULATING EXT. INV. COMPTON SPECTRUM ON BLR ... ");   
    
  

    if (CASE_JET == 0) stream_dat = fopen(stmp, "w+");

    for (i = 1; i <= NU_DIM; i++) {

        if (PRINT) fprintf(stderr, "%3i", i);
        
        
        // radiation field from the disc reprocessed by the BLR

        I_rad_ext[i] = L_BB_nuc[i] * LOR_B * tau / (4.0 * M_PI * R_blr*R_blr);
        //fprintf(stderr, "%6.3e\n", L_BB_nuc[i]);

        // emission & absorption coefficients
        
        //How the blob sees the BLR radiation
        jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad_ext, 
                    FreqTransO2S(NU_STR, DOP_B, z)*LOR_B,
                    FreqTransO2S(NU_END, DOP_B, z)*LOR_B,
                    NU_DIM, NU[i]*LOR_B, COM_PREC1, COM_PREC2); 
/*
        kk = gg_abs(NU[i]*LOR_B, I_rad_ext, NU_DIM, 
                        FreqTransO2S(NU_STR, DOP_B, z)*LOR_B, 
                        FreqTransO2S(NU_END, DOP_B, z)*LOR_B, COM_PREC1, COM_PREC2);
        //absorption on BLR done further
        */
        
        // absorption by pair production (inside the blob, on synch field)
        kk1 = gg_abs(NU[i], I_rad, NU_DIM, 
		      FreqTransO2S(NU_STR, DOP_B, z), 
		      FreqTransO2S(NU_END, DOP_B, z), SYN_PREC1, SYN_PREC2);

        if (L_src == 0.0) {
        I_com = SphTransfEquat(jj, kk1, R_src);
        } else {
        I_com = CylTransfEquat(0.0, jj, kk1, L_src);
        } 
        
        //I_com_ext[i] = I_com_ext[i] + I_com;
        I_com_ext[i] = I_com;
        
        
        if (CASE_JET == 0){       
            // frequency transformation to observer frame
            nu_tmp = FreqTransS2O(NU[i]*LOR_B, DOP_B, z);
            //fprintf(stderr, "%6.3e\n", NU[i]);
            
            // optical depth for absorption of VHE gamma rays by IIR
            if(IIR_level==1){

            if(EBLFLAG==0){
                if (z <= 0.1){
                    tt = tau_IRA_Kneiske(nu_tmp, z, 0);
                }else{
                    tt = tau_IRA_Kneiske(nu_tmp, z, 1);
                }
                }else if(EBLFLAG==1){
                    tt=tau_IRA_Franceschini(nu_tmp, z);
                }else if(EBLFLAG==2){
                    tt=tau_IRA_Finke(nu_tmp, z);
                }else if(EBLFLAG==3){
                    tt=tau_IRA_Franceschini17(nu_tmp, z);
                }
            }else tt=0.;
        
            // flux transformation to observer frame
            dtmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);

            // absorption by EBL
            F_IC_tot[i] = dtmp * exp(-tt);
            
            if (F_IC_tot[i] > 1.0e-300) {
                //fprintf(stderr, "%6.3e %6.3e\n", nu_tmp, nu_tmp*F_IC_tot[i]);
            fprintf(stream_dat, "%f %f %f %f %f\n", 
                    log10(nu_tmp), 
                    log10(F_IC_tot[i]), 
                    log10(nu_tmp*F_IC_tot[i]),
                    log10(dtmp),
                    log10(nu_tmp*dtmp) 
                    );
            }
        }


        if (PRINT) fprintf(stderr, "\b\b\b");  

    }   
    if (CASE_JET == 0) fclose(stream_dat);
    fprintf(stderr, "DONE\n\n");
    
    //EIC from dust torus
    sprintf(stmp, "./data/%s_ecs1.dat", prefix);
    if (ifexist(stmp)) {
        sprintf(comma, "mv ./data/%s_ecs1.dat ./data/%s_prev_ecs1.dat", prefix, prefix);
        if(system(comma)) return 1;
    }
    stream_dat = fopen(stmp, "w+");


        for (i = 1; i <= NU_DIM; i++) {             
        if (PRINT == 0) fprintf(stderr, "%3i", i);


        // Blackbody radiation field from the dust torus
        I_rad[i] = LOR_B * tau * L_tor/(4.0 * M_PI * R_blr*R_blr) *
            Planck((NU[i]/LOR_B), T_BB_tor) / (sig/M_PI * pow(T_BB_tor,4.0));


        // emission & absorption coefficients
        jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad, 
                    FreqTransO2S(NU_STR, DOP_B, z)*LOR_B,
                    FreqTransO2S(NU_END, DOP_B, z)*LOR_B,
                    NU_DIM, NU[i]*LOR_B, COM_PREC1, COM_PREC2); 
        
        kk = gg_abs(NU[i]*LOR_B, I_rad, NU_DIM, 
                        FreqTransO2S(NU_STR, DOP_B, z)*LOR_B, 
                        FreqTransO2S(NU_END, DOP_B, z)*LOR_B, COM_PREC1, COM_PREC2);
        
        
        // absorption by pair production 
        if (L_src == 0.0) {
            I_com = SphTransfEquat(jj, kk, R_src);
        } else {
            I_com = CylTransfEquat(0.0, jj, kk, L_src);
        }    
            I_com_ext1[i] = I_com;
        
        // frequency transformation to observer frame
        nu_tmp = FreqTransS2O(NU[i]*LOR_B, DOP_B, z);
        
        // optical depth for absorption of VHE gamma rays by IIR
        // absorption by IIR by Kneiske et al
        if(IIR_level==1){

            if(EBLFLAG==0){
            if (z <= 0.1){
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
            }else{
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
            }
            }else if(EBLFLAG==1){
                tt=tau_IRA_Franceschini(nu_tmp, z);
            }else if(EBLFLAG==2){
                tt=tau_IRA_Finke(nu_tmp, z);
            }else if(EBLFLAG==3){
                tt=tau_IRA_Franceschini17(nu_tmp, z);
            }            
        }
        else tt=0.;
        

        // flux transformation to observer frame
        fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);
        dtmp   = fx_tmp;
        
        // absorption by IIR
        fx_tmp = fx_tmp * exp(-tt); 
        
        if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", 
                    log10(nu_tmp), 
                    log10(fx_tmp), 
                    log10(nu_tmp*fx_tmp),
                    log10(dtmp),
                    log10(nu_tmp*dtmp) 
                    );             
        }             	 	 
        if (PRINT == 0) fprintf(stderr, "\b\b\b");          
    }   
    
    fclose(stream_dat);
    
    
    Ub_eicd = 4*M_PI/c*Simpson(I_com_ext, NU, NU_DIM, 1, NU_DIM);
    
    
    

    }else if (CASE_EIC && D_b > R_blr){
            fprintf(stderr, "... BLOB OUTSIDE BLR, EXT. INV. COMPTON ON BLR NOT COMPUTED\n\n");
    }
   }


   
   
   
   
   
//JET PART

   
   sprintf(stmp, "./data/F_jet_syn.dat");
   sprintf(stmp1, "./data/F_jet_syn_prev.dat");
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/F_jet_syn.dat ./data/F_jet_syn_prev.dat");
      if(system(comma)) return 1;
   } else if (ifexist(stmp1)) {
      sprintf(comma, "rm ./data/F_jet_syn_prev.dat");
      if(system(comma)) return 1;
   }
 
   sprintf(stmp, "./data/F_jet_com.dat");
   sprintf(stmp1, "./data/F_jet_com_prev.dat");
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/F_jet_com.dat ./data/F_jet_com_prev.dat");
      if(system(comma)) return 1;
   } else if (ifexist(stmp1)) {
      sprintf(comma, "rm ./data/F_jet_com_prev.dat");
      if(system(comma)) return 1;
   }

   sprintf(stmp, "./data/%s_ecs_jet.dat", prefix);
   sprintf(stmp1, "./data/%s_prev_ecs_jet.dat", prefix);
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/%s_ecs_jet.dat ./data/%s_prev_ecs_jet.dat", prefix, prefix);
      if(system(comma)) return 1;
   } else if (ifexist(stmp1)) {
      sprintf(comma, "rm ./data/%s_prev_ecs_jet.dat", prefix);
      if(system(comma)) return 1;
   }  
   

   if(CASE_JET){

   fprintf(stderr, "CALCULATING JET SYNCHROTRON SPECTRUM... "); 
   

   // CALCULATE PARAMETERS   
   //ref comoving
     
   X_MIN = Y_MIN / tan(PHI_src * M_PI / 180.0);
   X_MAX = X_MIN + J_LEN_src * pc;   
   A     = Y_MIN / X_MIN;
   Y_MAX = A * X_MAX; 
   
   // PREPARING FREQUENCY VECTOR
   
   tmp_min   = log10(FreqTransO2S(NU_STR, DOP_JET, z));
   tmp_max   = log10(FreqTransO2S(NU_END, DOP_JET, z));
   tmp_stp   = (tmp_max - tmp_min) / (NU_DIM);
   tmp_val   = tmp_min;
   
   
   for (i = 1; i <= NU_DIM; i++) {
      tmp_cur   = pow(10.0, tmp_val);
      NU[i] = tmp_cur;
      tmp_val   = tmp_val + tmp_stp;
   }
   
   NU[1] = FreqTransO2S(NU_STR, DOP_JET, z);
   
   
   // SPLITTING JET IN SLICES
   
   tmp_min   = log10(X_MIN);
   tmp_max   = log10(X_MAX);
   tmp_stp   = (tmp_max - tmp_min) / (SL_DIM);
   tmp_val   = tmp_min; 

   Y_VAL[0] = 0.0;

   for (i = 1; i <= SL_DIM + 1; i++) {
      tmp_cur   = pow(10.0, tmp_val);
      X_VAL[i]  = tmp_cur;
      Y_VAL[i]  = A * X_VAL[i];
      tmp_val   = tmp_val + tmp_stp;
   }          
     
   stream_dat = fopen("./data/geometry.dat", "w+");
   
   for (i = 1; i <= SL_DIM; i++) {
      DEL_X[i] = X_VAL[i+1] - X_VAL[i];
      //print distance of each slice from the jet basis, and radius
      fprintf(stream_dat, "%e %e %f %f\n", X_VAL[i], Y_VAL[i], log10(X_VAL[i]), log10(Y_VAL[i]));
      B_VAL[i] = B_0         * pow(X_MIN/X_VAL[i], n_B);      
      N_VAL[i] = N_0* Y_MIN*Y_MIN / (Y_VAL[i]*Y_VAL[i]); //follow the jet gometry
      G_VAL[i] = GAMMA_MAX_0 * pow(X_MIN/X_VAL[i], n_G);
   }   
   
   fclose(stream_dat);
   
   
   // CALCULATING JET SYNCHROTRON SPECTRUM
   
   // CALCULATING ELECTRONS ENERGY SPECTRUM
   

   Lj_B = 0;
   Lj_e = 0;
   for (j = 1; j <= SL_DIM; j++) {
      sprintf(stmp, "./data/es_jet_slice_%d.dat",j);
      if (ifexist(stmp)) {
	  sprintf(comma, "mv ./data/es_jet_slice_%d.dat ./data/prev_es_jet_slice_%d.dat",j,j);
	  if(system(comma)) return 1;
      }
      stream_dat = fopen(stmp, "w+");  
      
      tmp_min   = log10(GAMMA_MIN1);
      tmp_max   = log10(GAMMA_MAX_0);
      tmp_stp   = (tmp_max - tmp_min) / (G_DIM - 1);
      tmp_val   = tmp_min;
      SL_CUR = j;
      for (i = 1; i <= G_DIM ; i++) {
         tmp_cur = pow(10.0, tmp_val);
         dtmp    = N_e_Jet(tmp_cur);
         //STORING VECTORS FOR THE COMPUTATION OF ENERGY DENSITIES
         Gamma[i] = tmp_cur;
         elec_spec[i] = dtmp*tmp_cur;
         tmp_val = tmp_val + tmp_stp;
	 
	 fprintf(stream_dat, "%e %e\n",tmp_cur , dtmp*tmp_cur);
      }
      fclose(stream_dat); 
      
      elec_spec_tot_j = Simpson(elec_spec, Gamma, G_DIM, 1, G_DIM);
      U_B = B_VAL[j]*B_VAL[j]/(8.*M_PI);
      Uj_e = m_e*c*c*elec_spec_tot_j;
      
      //We take powers through the first slice
      if (j == 1){
        Lj_e = M_PI*Y_VAL[j]*Y_VAL[j]*LOR_JET*LOR_JET*c*Uj_e;
	Lj_B = Lj_B + M_PI*Y_VAL[j]*Y_VAL[j]*LOR_JET*LOR_JET*c*U_B;
	Utot_e = Utot_e + Uj_e;
	Utot_B = Utot_B + U_B;
      }
      R_MOY = (Y_MAX - Y_MIN)/2;
    }
    

   // CALCULATE COEFFICIENTS OF JET
   
   sprintf(stmp, "./data/coeff_jet.dat");
   stream_dat1 = fopen(stmp, "w+");
   
   //fprintf(stderr, "Jet slices coefficients (comoving frame)\n");
   for (j = 1; j <= SL_DIM; j++) {
     if (PRINT) fprintf(stderr, "%3i", j);

      SL_CUR = j;
      //fprintf(stderr, "slice: %3d X: %e dX: %e N_0: %e N: %e  B: %e  gamma_max: %e\n", j, X_VAL[j], DEL_X[j], N_VAL[j], int_spec_j(N_VAL[j]), B_VAL[j], G_VAL[j]); 
      //fprintf(stream_dat1, "%d\t%e\t%e\t%e\t%e\t%e\t%e\n", j, X_VAL[j], DEL_X[j], N_VAL[j], int_spec_j(N_VAL[j]), B_VAL[j], G_VAL[j]);
      
      //Stot: Projected area of a slice, following the angle THETA_src_j
      Stot[j] = M_PI*Y_VAL[j]*Y_VAL[j]*cos(THETA_src_j*M_PI/180.0) + Y_VAL[j]*DEL_X[j]*sin(THETA_src_j*M_PI/180.0);
      
   
      sprintf(stmp, "./data/f_syn_slice_%d.dat", j);
      stream_dat = fopen(stmp, "w+");

      for (i = 1; i <= NU_DIM; i++) {      
         J_SYN_JET[j][i] = j_syn(N_e_Jet, GAMMA_MIN1, G_VAL[j], NU[i], B_VAL[j], SYN_PREC1, SYN_PREC2);
         K_ESA_JET[j][i] = k_esa(N_e_Jet, GAMMA_MIN1, G_VAL[j], NU[i], B_VAL[j], ABS_PREC1, ABS_PREC2);      

         I_SYN_JET_BASE[j][i] = CylTransfEquat(0.0, J_SYN_JET[j][i], K_ESA_JET[j][i], DEL_X[j]);
	 
	 //Edge surface intensity
	 if (THETA_src_j > PHI_src) {
	   I_SYN_JET_EDGE[j][i] = SphTransfEquat(J_SYN_JET[j][i], K_ESA_JET[j][i], Y_VAL[j]);
	 } else {
	   I_SYN_JET_EDGE[j][i] = 0.0;
	 }
	 
	 //Intensity of a slice following the angle of observation (comoving frame)
	 I_SYN_JET[j][i] = ((fabs((cos(THETA_src_j * M_PI / 180.0))) * M_PI * (Y_VAL[j]*Y_VAL[j]) * I_SYN_JET_BASE[j][i]) +
	     (sin(THETA_src_j * M_PI / 180.0) * 2 *Y_VAL[j] * DEL_X[j] * I_SYN_JET_EDGE[j][i])) / Stot[j];
	     
	 
	 // frequency transformation to observer frame
         nu_tmp = FreqTransS2O(NU[i], DOP_JET, z);

	 //misaligned jet
	 fx_tmp = CylIntens2Flux(I_SYN_JET_BASE[j][i], I_SYN_JET_EDGE[j][i], Y_VAL[j], DEL_X[j],
                   DOP_JET, z, H_0, THETA_src_j);
	 
      
         fprintf(stream_dat, "%e %e %e\n", nu_tmp, fx_tmp, nu_tmp*fx_tmp);      
      }        
      if (PRINT) fprintf(stderr, "\b\b\b");
      fclose(stream_dat);     
   } 
   fclose(stream_dat1);
   
   
   // CALCULATE TOTAL SYNCHROTRON SPECTRUM
   
  //new method 03_2016
  
  for (k = 1; k <= SL_DIM; k++) {
    //time gap before the central emission of a slice escape the jet by the edge
    DEL_Tph = Y_VAL[k] / (c * sin(THETA_src_j*M_PI/180.0) - V_exp);
    
    for (i = 1; i <= NU_DIM; i++) {
      I_SYN = 0.0;
      
      for (j = k; j <= SL_DIM; j++) {
	sprintf(stmp, "%e", J_SYN_JET[j][i]);
														
	if ((strcmp(stmp, "-inf") == 0) ||
	    (strcmp(stmp, "inf")  == 0) ||
	    (strcmp(stmp, "nan")  == 0) ||
	    (J_SYN_JET[j][i] < 1.0e-300)) {
	  I_SYN_JET[j][i] = 0.0;
	
	} else {
          //evolution time of a slice
          DEL_Tj = DEL_X[j]/V_JET;

	  //photons are crossing the slice before its evolution
	  if (DEL_Tj >= DEL_Tph){
	    dtmp = DEL_Tph*c * K_ESA_JET[j][i];
	  }
	  // The slice is evolving during the photons passage
	  else{
	    dtmp = DEL_Tj*c * K_ESA_JET[j][i];
	    DEL_Tph -= DEL_Tj;
	    if (DEL_Tph < 0.0) DEL_Tph = 0.0;
	  }	
	  
	  // thick
	  if (dtmp > 7.00e+2)                       I_SYN = I_SYN * 0.0        + I_SYN_JET[j][i];
	  // transparent
	  if (dtmp < 1.0e-10)                       I_SYN = I_SYN * 1.0        + I_SYN_JET[j][i];  
	  // thin
	  if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) I_SYN = I_SYN * exp(-dtmp) + I_SYN_JET[j][i];	  
	}
      }
      
      //intensity of a slice k at the end of its radiation transfert through the other slices downstream
      I_SYN_TOT[k][i] = I_SYN; 
      // transformation to observer frame (take the non-previously covered slice corona)
      F_SYN_JET[k][i] = M_PI * ((pow(Y_VAL[k], 2.0) - pow(Y_VAL[k-1], 2.0)) / (D_L * D_L)) * (1.0 + z) * I_SYN;
    }
  }
	
   sprintf(stmp, "./data/F_jet_syn.dat");
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/F_jet_syn.dat ./data/F_jet_syn_prev.dat");
      if(system(comma)) return 1;
   }
   errno  = 0;
   stream_dat = fopen(stmp, "w+");  
   
   sprintf(stmp1, "./data/F_jet_frame_syn.dat");
   stream_dat1 = fopen(stmp1, "w+");
   
   if (errno == 0) {              
        
     for (i = 1; i <= NU_DIM; i++) {
        dtmp = 0.0;
        for (j = 1; j <= SL_DIM; j++) {
           dtmp = dtmp + F_SYN_JET[j][i];
        }

        nu_tmp = FreqTransS2O(NU[i], DOP_JET, z);
        fx_tmp = ftr(dtmp);
	
	dtmp1 = log10(nu_tmp);
	dtmp2 = log10(fx_tmp);
	dtmp3 = log10(nu_tmp*fx_tmp);
	
	sprintf(stmp1, "%f", dtmp2);
	sprintf(stmp2, "%f", dtmp3);
	
        if ((strcmp(stmp1, "-inf") == 0) ||
            (strcmp(stmp1, "inf")  == 0) ||
            (strcmp(stmp1, "nan")  == 0) ||
	    (strcmp(stmp2, "-inf") == 0) ||
            (strcmp(stmp2, "inf")  == 0) ||
            (strcmp(stmp2, "nan")  == 0)) {  
        } else {
          fprintf(stream_dat, "%e %e %e %f %f %f\n",
                  nu_tmp, fx_tmp, nu_tmp * fx_tmp,
                  dtmp1,   dtmp2, dtmp3);
        } 
        
        nu_tmp1 = FreqTransS2O(NU[i], 1, z);	
        fx_tmp1 = dtmp * (1.0 + z);
        if (fx_tmp1 > 1.0e-300) {
           fprintf(stream_dat1, "%f %f %f\n", 
            log10(nu_tmp1), 
            log10(fx_tmp1), 
            log10(nu_tmp1*fx_tmp1)
	     );
        }
        
     }                 
     fclose(stream_dat);
     fclose(stream_dat1);
   } else {
     fprintf(stderr, "ERROR (%d) %s\n", errno, strerror(errno));
     fprintf(stderr, "TOTAL SYNCHROTRON SPECTRUM: '%s' NOT SAVED !!!\n", name);
     return 0;
   } 
   
   Pj_syn = 0.0;
   for (j = 1; j <= SL_DIM; j++) {
     UJ_SLICE[j] = 4*M_PI/c*Simpson(I_SYN_TOT[j], NU, NU_DIM, 1, NU_DIM);
     PJ_SLICE[j] = M_PI * (pow(Y_VAL[j], 2.0) - pow(Y_VAL[j-1], 2.0))  * LOR_JET*LOR_JET * c * UJ_SLICE[j];
     Pj_syn = Pj_syn + PJ_SLICE[j];
   }
   
   fprintf(stderr, " DONE\n\n");
   
   
  

     
     
   
   
   
   
   
   fprintf(stderr, "CALCULATING JET SSC SPECTRUM ... ");   
      
   // CALCULATING JET INVERSE-COMPTON SPECTRUM
      
   
   // CALCULATE COEFFICIENTS OF JET
   
   for (j = 1; j <= SL_DIM; j++) {
     if (PRINT) fprintf(stderr, "%3i", j);
      SL_CUR = j;
   
      sprintf(stmp, "./data/f_ssc_slice_%d.dat", j);
      stream_dat = fopen(stmp, "w+");

      for (i = 1; i <= NU_DIM; i++) { 
	
	 // emission & absorption coefficients
         J_COM_JET[j][i] = j_com(N_e_Jet, GAMMA_MIN1, G_VAL[j], I_SYN_JET[j], 
			FreqTransO2S(NU_STR, DOP_JET, z),
			FreqTransO2S(NU_END, DOP_JET, z),
			NU_DIM, NU[i], COM_PREC1, COM_PREC2);

	 K_ABS_SSC_JET[j][i] = gg_abs(NU[i], I_SYN_JET[j], NU_DIM, 
					  FreqTransO2S(NU_STR, DOP_JET, z), 
					  FreqTransO2S(NU_END, DOP_JET, z), 
					  SYN_PREC1, SYN_PREC2); 

	 
	 // absorption by pair production
         I_COM_JET_BASE[j][i] = CylTransfEquat(0.0, J_COM_JET[j][i], K_ABS_SSC_JET[j][i], DEL_X[j]);
	 
	 //Edge surface intensity
	 if (THETA > PHI) {
	   I_COM_JET_EDGE[j][i] = SphTransfEquat(J_COM_JET[j][i], K_ABS_SSC_JET[j][i], Y_VAL[j]);
	 } else {
	   I_COM_JET_EDGE[j][i] = 0.0;
	 }

	 // frequency transformation to observer frame
         nu_tmp = FreqTransS2O(NU[i], DOP_JET, z);
	 
	 // optical depth for absorption of VHE gamma rays by IIR
        if(IIR_level==1){

          if(EBLFLAG==0){
             if (z <= 0.1){
	        tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	     }else{
	        tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	     }
           }else if(EBLFLAG==1){
                tt=tau_IRA_Franceschini(nu_tmp, z);
	   }else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
           }else if(EBLFLAG==3){
             tt=tau_IRA_Franceschini17(nu_tmp, z);
	   }
         }else tt=0.;	 
	 
	 // flux transformation to observer frame
	 
	 I_COM_JET[j][i] = ((fabs((cos(THETA_src_j * M_PI / 180.0))) * M_PI * (Y_VAL[j]*Y_VAL[j]) * I_COM_JET_BASE[j][i]) +
	     (sin(THETA_src_j * M_PI / 180.0) * 2 *Y_VAL[j] * DEL_X[j] * I_COM_JET_EDGE[j][i])) / Stot[j];
	 
	 //misaligned jet
	 fx_tmp = CylIntens2Flux(I_COM_JET_BASE[j][i], I_COM_JET_EDGE[j][i], Y_VAL[j], DEL_X[j],
                   DOP_JET, z, H_0, THETA_src_j);
	 
         // absorption by IIR
         fx_tmp = fx_tmp * exp(-tt); 
      
         fprintf(stream_dat, "%e %e %e\n", nu_tmp, fx_tmp, nu_tmp*fx_tmp);      
      }          
      fclose(stream_dat);    
      
      if (PRINT) fprintf(stderr, "\b\b\b");
   }      
   
   // CALCULATE TOTAL COMPTON SPECTRUM
   
   
  for (k = 1; k <= SL_DIM; k++) {
    //time gap before the central emission of a slice escape the jet by the edge
    DEL_Tph = Y_VAL[k] / (c * sin(THETA_src_j*M_PI/180.0) - V_exp);
    
    for (i = 1; i <= NU_DIM; i++) {
      I_COM = 0.0;
      
      for (j = k; j <= SL_DIM; j++) {
	sprintf(stmp, "%e", J_COM_JET[j][i]);
														
	if ((strcmp(stmp, "-inf") == 0) ||
	    (strcmp(stmp, "inf")  == 0) ||
	    (strcmp(stmp, "nan")  == 0) ||
	    (J_COM_JET[j][i] < 1.0e-300)) {
	  J_COM_JET[j][i] = 0.0;
	
	} else {
          //evolution time of a slice
          DEL_Tj = DEL_X[j]/V_JET;

	  //photons are crossing the slice before its evolution
	  if (DEL_Tj >= DEL_Tph){
	    dtmp = DEL_Tph*c * K_ABS_SSC_JET[j][i];
	  }
	  // The slice is evolving during the photons passage
	  else{
	    dtmp = DEL_Tj*c * K_ABS_SSC_JET[j][i];
	    DEL_Tph -= DEL_Tj;
	    if (DEL_Tph < 0.0) DEL_Tph = 0.0;
	  }	
	  
	  if (dtmp > 7.00e+2)                       I_COM = I_COM * 0.0        + I_COM_JET[j][i];
	  if (dtmp < 1.0e-10)                       I_COM = I_COM * 1.0        + I_COM_JET[j][i];
	  if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) I_COM = I_COM * exp(-dtmp) + I_COM_JET[j][i];	  
	}
      }
      
      I_COM_TOT[k][i] = I_COM;
      // transformation to observer frame
      F_COM_JET[k][i] = M_PI * ((pow(Y_VAL[k], 2.0) - pow(Y_VAL[k-1], 2.0)) / (D_L * D_L)) * (1.0 + z) * I_COM;
    }
  }
   
   sprintf(stmp, "./data/F_jet_com.dat");
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/F_jet_com.dat ./data/F_jet_com_prev.dat");
      if(system(comma)) return 1;
   }
   errno  = 0;
   stream_dat = fopen(stmp, "w+");  
   
   sprintf(stmp1, "./data/F_jet_frame_com.dat");
   stream_dat1 = fopen(stmp1, "w+");
   
   if (errno == 0) {              
        
     for (i = 1; i <= NU_DIM; i++) {
        dtmp = 0.0;
        for (j = 1; j <= SL_DIM; j++) {
           dtmp = dtmp + F_COM_JET[j][i];
        }

        nu_tmp = FreqTransS2O(NU[i], DOP_JET, z);
        fx_tmp = ftr(dtmp);
	
	dtmp1 = log10(nu_tmp);
	dtmp2 = log10(fx_tmp);
	dtmp3 = log10(nu_tmp*fx_tmp);
	
	sprintf(stmp1, "%f", dtmp2);
	sprintf(stmp2, "%f", dtmp3);
	
        if ((strcmp(stmp1, "-inf") == 0) ||
            (strcmp(stmp1, "inf")  == 0) ||
            (strcmp(stmp1, "nan")  == 0) ||
	    (strcmp(stmp2, "-inf") == 0) ||
            (strcmp(stmp2, "inf")  == 0) ||
            (strcmp(stmp2, "nan")  == 0)) {
        } else {
          fprintf(stream_dat, "%e %e %e %f %f %f\n",
                  nu_tmp, fx_tmp, nu_tmp * fx_tmp,
                  dtmp1,   dtmp2, dtmp3);
        }  
        
        nu_tmp1 = FreqTransS2O(NU[i], 1, z);	
        fx_tmp1 = dtmp * (1.0 + z);
        if (fx_tmp1 > 1.0e-300) {
           fprintf(stream_dat1, "%f %f %f\n", 
            log10(nu_tmp1), 
            log10(fx_tmp1), 
            log10(nu_tmp1*fx_tmp1)
	     );
        }
        
     }                 
     fclose(stream_dat);
     fclose(stream_dat1);
     
   } else {
     fprintf(stderr, "ERROR (%d) %s\n", errno, strerror(errno));
     fprintf(stderr, "TOTAL SSC SPECTRUM: '%s' NOT SAVED !!!\n", name);
     return 0;
   } 
   
   Pj_ssc = 0.0;
   for (j = 1; j <= SL_DIM; j++) {
     UJ_SLICE[j] = 4*M_PI/c*Simpson(I_COM_TOT[j], NU, NU_DIM, 1, NU_DIM);
     PJ_SLICE[j] = M_PI * (pow(Y_VAL[j], 2.0) - pow(Y_VAL[j-1], 2.0))  *LOR_JET*LOR_JET * c * UJ_SLICE[j];
     Pj_ssc = Pj_ssc + PJ_SLICE[j];
   }
   
   fprintf(stderr, " DONE\n\n");
   
   
   

   
   
   
   
   
   
   
   fprintf(stderr, "CALCULATING EXT. INV. COMPTON SPECTRUM ON JET SYNCHROTRON... ");   
      
   // CALCULATING EXT. INV. COMPTON SPECTRUM ON JET SYNCHROTRON
   
    

   k = 1;
   while (X_VAL[k+1] <= D_b_src_j) {
     k += 1;
   }
   //if the blob is inside the jet
   if ( Y_MIN / tan(PHI_src * M_PI / 180.0) < D_b_src_j){
    //radiation transfert along the jet until the blob position, independent of the angle THETA
    //parallel plans
        for (l = 1; l <= k; l++) {
            //fprintf(stderr, "%e\n\n",l); 
            for (i = 1; i <= NU_DIM; i++) {
            I_SYN = 0.0;
            for (j = l; j <= k; j++) {
                sprintf(stmp, "%e", J_SYN_JET[j][i]);
                //fprintf(stderr, "%d\n",j); 								      
                if ((strcmp(stmp, "-inf") == 0) ||
                    (strcmp(stmp, "inf")  == 0) ||
                    (strcmp(stmp, "nan")  == 0) ||
                    (J_SYN_JET[j][i] < 1.0e-300)) {
                    I_SYN_JET[j][i] = 0.0;
                } else {
                    if (j < k){
                    dtmp = DEL_X[j] * K_ESA_JET[j][i];
                    // transparent
                    if (dtmp < 1.0e-10)                       I_SYN = I_SYN * 1.0      + J_SYN_JET[j][i] * DEL_X[j];
                    }else{
                    dtmp = (D_b_src_j - X_VAL[k]) * K_ESA_JET[j][i];
                    // transparent
                    if (dtmp < 1.0e-10)                       I_SYN = I_SYN * 1.0      + J_SYN_JET[j][i] * (D_b_src_j - X_VAL[k]);
                    }
                
                    // thick
                    if (dtmp > 7.00e+2)                       I_SYN = I_SYN * 0.0        + J_SYN_JET[j][i] / K_ESA_JET[j][i];
                    // thin
                    if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) I_SYN = I_SYN * exp(-dtmp) + I_SYN_JET[j][i];
                }
            }
            
            I_SYN_TOT[l][i] = I_SYN; 
            }
        }
   }
    
    
   for (i = 1; i <= NU_DIM; i++) {
     if (PRINT) fprintf(stderr, "%3i", i);
     
     I_rad_ext_s[i] = I_SYN_TOT[k][i];

     if (PRINT) fprintf(stderr, "\b\b\b");
   }
     
     
     
     
   //radiation transfert in the direction opposite of the jet propagation
   for (l = SL_DIM; l >= k; l--) {
     for (i = 1; i <= NU_DIM; i++) {
         I_SYN = 0.0;
         for (j = l; j >= k; j--) {
	    sprintf(stmp, "%e", J_SYN_JET[j][i]);                                                                                            	    
	    if ((strcmp(stmp, "-inf") == 0) ||
                (strcmp(stmp, "inf")  == 0) ||
                (strcmp(stmp, "nan")  == 0) ||
	        (J_SYN_JET[j][i] < 1.0e-300)) {
	      I_SYN_JET[j][i] = 0.0;
	    
	    } else {
	      if (j > k){ 
                dtmp = DEL_X[j] * K_ESA_JET[j][i];
	        // transparent
	        if (dtmp < 1.0e-10)                       I_SYN = I_SYN * 1.0        + J_SYN_JET[j][i] * DEL_X[j];  
		
	      } else {
		dtmp = (X_VAL[k+1] - D_b) * K_ESA_JET[j][i];
	        // transparent
	        if (dtmp < 1.0e-10)                       I_SYN = I_SYN * 1.0        + J_SYN_JET[j][i] * (X_VAL[k+1] - D_b);      
		
	      } 
		// thick
		if (dtmp > 7.00e+2)                       I_SYN = I_SYN * 0.0        + J_SYN_JET[j][i] / K_ESA_JET[j][i];
		// transparent
		//if (dtmp < 1.0e-10)                       I_SYN = I_SYN * 1.0        + J_SYN_JET[j][i] * DEL_X[j];  
		// thin
		if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) I_SYN = I_SYN * exp(-dtmp) + I_SYN_JET[j][i];	      
	    }
         }
	 
	 I_SYN_TOT[l][i] = I_SYN; 
     }
   }
   
   for (i = 1; i <= NU_DIM; i++) {
     if (PRINT) fprintf(stderr, "%3i", i);
     
     I_rad_ext1[i] = I_SYN_TOT[l][i];
     

     
     //if the blob is inside the jet
     if ( Y_MIN / tan(PHI_src * M_PI / 180.0) < D_b_src_j){
        // transformation jet to blob frame (blueshift)
        I_rad_ext1[i] = (I_rad_ext1[i] + I_rad_ext_s[i]) * LOR_B_J;
     } else {
         //if the blob is upstream the jet (i.e upstream radio core)
         I_rad_ext1[i] = (I_rad_ext1[i] + I_rad_ext_s[i]) * LOR_B_J / (4.0 * M_PI * (D_BJ + Y_MIN)*(D_BJ + Y_MIN) / (Y_MIN*Y_MIN));
     }
     
     // emission & absorption coefficients (how the blob sees the jet)
     C_e1[i] = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad_ext1,
		FreqTransO2S(NU_STR, DOP_JET, z)*LOR_B_J,
		FreqTransO2S(NU_END, DOP_JET, z)*LOR_B_J,
		NU_DIM, NU[i]*LOR_B_J, COM_PREC1, COM_PREC2);
     
     C_a1[i] = gg_abs(NU[i]*LOR_B_J, I_rad_ext1, NU_DIM, 
		      FreqTransO2S(NU_STR, DOP_JET, z)*LOR_B_J, 
		      FreqTransO2S(NU_END, DOP_JET, z)*LOR_B_J, 
		      COM_PREC1, COM_PREC2);
     
     if (PRINT) fprintf(stderr, "\b\b\b");
   }
   
   sprintf(stmp1, "./data/%s_ecs_jet_frame.dat", prefix);
   stream_dat1 = fopen(stmp1, "w+");

  
   for (i = 1; i <= NU_DIM; i++) {    
     jj = C_e1[i] ;
     kk = C_a1[i] ;
     
     
     // radiation transfert taking into account the pair absorption
     if (L_src == 0.0) {
       I_com = SphTransfEquat(jj, kk, R_src);
       I_eic_jet[i] = I_com;
     } else {
       I_com = CylTransfEquat(0.0, jj, kk, L_src);
       I_eic_jet[i] = I_com;
     } 
           
     //transformation in the observer frame
     nu_tmp1 = FreqTransS2O(NU[i]*LOR_B_J, DOP_B, z);
     fx_tmp1 = Intens2Flux(I_com, R_src, 1, z, H_0);
     if (fx_tmp1 > 1.0e-300) {
       fprintf(stream_dat1, "%f %f %f\n", 
            log10(nu_tmp1), 
            log10(fx_tmp1), 
            log10(nu_tmp1*fx_tmp1)
	     );
     } 
   }
   fclose(stream_dat1);

   fprintf(stderr, "DONE\n\n");
   
   Ub_eicj = 4*M_PI/c*Simpson(I_eic_jet, NU, NU_DIM, 1, NU_DIM);
   }

   
   
    if (D_b < R_blr){
    fprintf(stderr, "CALCULATING RADIATION ABSORPTION OF THE BLOB EMISSION THROUGH THE BLR ... ");

    // PREPARING FREQUENCY VECTOR
    tmp_min   = log10(FreqTransO2S(NU_STR, DOP_B, z)*LOR_B);
    tmp_max   = log10(FreqTransO2S(NU_END, DOP_B, z)*LOR_B);
    tmp_stp   = (tmp_max - tmp_min) / (NU_DIM);
    tmp_val   = tmp_min; 
    for (i = 1; i <= NU_DIM; i++) {
        tmp_cur   = pow(10.0, tmp_val);
        NU[i] = tmp_cur;
        tmp_val   = tmp_val + tmp_stp;
    }
    //Warning:: NU[1] < FreqTransO2S(NU_STR, DOP_B, z) !!
    NU[1] = FreqTransO2S(NU_STR, DOP_B, z)*LOR_B;
    
        
    //blob frame
    for (i = 1; i <= NU_DIM; i++) {
        //fprintf(stderr, "%6.3e\n", NU[i]);
        
        // absorption by pair production (on BLR soft photons)
        kk = gg_abs(NU[i], I_rad_ext, NU_DIM, 
                        FreqTransO2S(NU_STR, DOP_B, z), 
                        FreqTransO2S(NU_END, DOP_B, z), COM_PREC1, COM_PREC2);
        
        
        dtmp = (R_blr - D_b) *kk;
        //fprintf(stderr, "%6.3e\n", dtmp);
        // thick
        if (dtmp > 7.00e+2){
            I_rad_com[i] = 0.0;
            I_rad2nd[i] = 0.0;
        }
        // transparent
        //nothing change
        // thin
        if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)){
            I_rad_com[i] = I_rad_com[i] * exp(-dtmp);
            I_rad2nd[i] = I_rad2nd[i] * exp(-dtmp);
        }
    }
    
    // PREPARING FREQUENCY VECTOR
    tmp_min   = log10(FreqTransO2S(NU_STR, DOP_B, z));
    tmp_max   = log10(FreqTransO2S(NU_END, DOP_B, z));
    tmp_stp   = (tmp_max - tmp_min) / (NU_DIM);
    tmp_val   = tmp_min; 
    for (i = 1; i <= NU_DIM; i++) {
        tmp_cur   = pow(10.0, tmp_val);
        NU[i] = tmp_cur;
        tmp_val   = tmp_val + tmp_stp;
    }
    //Warning:: NU[1] < FreqTransO2S(NU_STR, DOP_B, z) !!
    NU[1] = FreqTransO2S(NU_STR, DOP_B, z);
    
        
    //blob frame
    for (i = 1; i <= NU_DIM; i++) {
        //fprintf(stderr, "%6.3e\n", NU[i]);
        
        // absorption by pair production (on BLR soft photons)
        kk = gg_abs(NU[i]*LOR_B, I_rad_ext, NU_DIM, 
                        FreqTransO2S(NU_STR, DOP_B, z)*LOR_B, 
                        FreqTransO2S(NU_END, DOP_B, z)*LOR_B, COM_PREC1, COM_PREC2);
        
        
        dtmp = (R_blr - D_b) *kk;
        //fprintf(stderr, "%6.3e\n", dtmp);
        // thick
        if (dtmp > 7.00e+2){
            I_com_ext[i] = 0.0;

        }
        // transparent
        //nothing change
        // thin
        if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)){
            I_com_ext[i] = I_com_ext[i] * exp(-dtmp);
        }
    }

    if (CASE_JET == 0){
        //SSC
        sprintf(stmp, "./data/%s_cs.dat", prefix);
        if (ifexist(stmp)) {
            sprintf(comma, "mv ./data/%s_cs.dat ./data/%s_prev_cs.dat", prefix, prefix);
            if(system(comma)) return 1;
        }
        stream_dat = fopen(stmp, "w+");
        
        
        for (i = 1; i <= NU_DIM; i++) {
        // frequency transformation to observer frame
        nu_tmp = FreqTransS2O(NU[i], DOP_B, z);	 

        // optical depth for absorption of VHE gamma rays by IIR
        //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
        if(IIR_level==1){
            if(EBLFLAG==0){
                if (z <= 0.1){
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
                }else{
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
                }
            }else if(EBLFLAG==1){
                tt=tau_IRA_Franceschini(nu_tmp, z);
            }else if(EBLFLAG==2){
                tt=tau_IRA_Finke(nu_tmp, z);
            }else if(EBLFLAG==3){
                tt=tau_IRA_Franceschini17(nu_tmp, z);
            }

        }else tt=0.;


        // flux transformation to observer frame
        fx_tmp = Intens2Flux(I_rad_com[i], R_src, DOP_B, z, H_0);
        dtmp   = fx_tmp;

        // absorption by IIR
        fx_tmp = fx_tmp * exp(-tt); 

        if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
            log10(nu_tmp), 
            log10(fx_tmp), 
            log10(nu_tmp*fx_tmp),
            log10(dtmp),
            log10(nu_tmp*dtmp) 
            );
        }
        }
    fclose(stream_dat);
    
    //SSC_2nd
    sprintf(stmp, "./data/%s_cs2.dat", prefix);
    if (ifexist(stmp)) {
        sprintf(comma, "mv ./data/%s_cs2.dat ./data/%s_prev_cs2.dat", prefix, prefix);
        if(system(comma)) return 1;
    }
    stream_dat = fopen(stmp, "w+");  
    
    for (i = 1; i <= NU_DIM; i++) {
        // frequency transformation to observer frame
        nu_tmp = FreqTransS2O(NU[i], DOP_B, z);	 

        // optical depth for absorption of VHE gamma rays by IIR
        //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
        if(IIR_level==1){
            if(EBLFLAG==0){
                if (z <= 0.1){
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
                }else{
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
                }
            }else if(EBLFLAG==1){
                tt=tau_IRA_Franceschini(nu_tmp, z);
            }else if(EBLFLAG==2){
                tt=tau_IRA_Finke(nu_tmp, z);
            }else if(EBLFLAG==3){
                tt=tau_IRA_Franceschini17(nu_tmp, z);
            }

        }else tt=0.;


        // flux transformation to observer frame
        fx_tmp = Intens2Flux(I_rad2nd[i], R_src, DOP_B, z, H_0);
        dtmp   = fx_tmp;

        // absorption by IIR
        fx_tmp = fx_tmp * exp(-tt); 

        if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
            log10(nu_tmp), 
            log10(fx_tmp), 
            log10(nu_tmp*fx_tmp),
            log10(dtmp),
            log10(nu_tmp*dtmp) 
            );
        }
        }
    fclose(stream_dat);  
    
    //EIC from disk-BLR
    sprintf(stmp, "./data/%s_ecs.dat", prefix);
    if (ifexist(stmp)) {
        sprintf(comma, "mv ./data/%s_ecs.dat ./data/%s_prev_ecs.dat", prefix, prefix);
        if(system(comma)) return 1;
    }
    stream_dat = fopen(stmp, "w+");
    
    for (i = 1; i <= NU_DIM; i++) {
        // frequency transformation to observer frame
        nu_tmp = FreqTransS2O(NU[i]*LOR_B, DOP_B, z); 

        // optical depth for absorption of VHE gamma rays by IIR
        //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
        if(IIR_level==1){
            if(EBLFLAG==0){
                if (z <= 0.1){
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
                }else{
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
                }
            }else if(EBLFLAG==1){
                tt=tau_IRA_Franceschini(nu_tmp, z);
            }else if(EBLFLAG==2){
                tt=tau_IRA_Finke(nu_tmp, z);
            }else if(EBLFLAG==3){
                tt=tau_IRA_Franceschini17(nu_tmp, z);
            }

        }else tt=0.;


        // flux transformation to observer frame
        fx_tmp = Intens2Flux(I_com_ext[i], R_src, DOP_B, z, H_0);
        dtmp   = fx_tmp;

        // absorption by IIR
        fx_tmp = fx_tmp * exp(-tt); 

        if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
            log10(nu_tmp), 
            log10(fx_tmp), 
            log10(nu_tmp*fx_tmp),
            log10(dtmp),
            log10(nu_tmp*dtmp) 
            );
        }
        }
    fclose(stream_dat);  
    
    }
    }
    
   
   
   
   if(CASE_JET){   
   
   fprintf(stderr, "CALCULATING RADIATION ABSORPTION OF THE BLOB EMISSION THROUGH THE JET ... ");
   
   // PREPARING FREQUENCY VECTOR
   
   tmp_min   = log10(FreqTransO2S(NU_STR, DOP_B, z));
   tmp_max   = log10(FreqTransO2S(NU_END, DOP_B, z));
   tmp_stp   = (tmp_max - tmp_min) / (NU_DIM);
   tmp_val   = tmp_min;
   
   for (i = 1; i <= NU_DIM; i++) {
      tmp_cur   = pow(10.0, tmp_val);
      NU[i] = tmp_cur;
      tmp_val   = tmp_val + tmp_stp;
   }

   j = 1;
   while (X_VAL[j] <= D_b_src_j) {
     j += 1;
     // j: first slice after the blob
   }
   
   //blob to jet frame
   for (i = 1; i <= NU_DIM; i++) {
      I_rad_syn[i] = I_rad_syn[i] * pow(DOP_B_J,3.0);
      I_rad_com[i] = I_rad_com[i] * pow(DOP_B_J,3.0);
      I_rad2nd[i] = I_rad2nd[i] * pow(DOP_B_J,3.0);
      I_com_ext[i] = I_com_ext[i] * pow(DOP_B_J,3.0);
      I_com_ext1[i] = I_com_ext1[i] * pow(DOP_B_J,3.0);
      I_eic_jet[i] = I_eic_jet[i] * pow(DOP_B_J,3.0);
   }
   
    //new method 03_2016
    //radiative transfert along the propagation of the jet
    
    //time gap before the blob emission escape the jet by the edge
    DEL_Tph = Y_VAL[j-1] / (c * sin(THETA_src_j*M_PI/180.0) - V_exp);
    k = j-1;
    
    while (k <= SL_DIM && DEL_Tph > 0.0 ){
      fprintf(stderr, "%3i", k);
      //evolution time of a slice
      if (k == j-1){
	DEL_Tj = (X_VAL[k+1] - D_b_src_j) / V_JET;
      }else{
	DEL_Tj = DEL_X[k]/V_JET;
      }
      for (i = 1; i <= NU_DIM; i++) {
      //photons cross the slice before its evolution
      if (DEL_Tj >= DEL_Tph){
	dtmp = DEL_Tph*c * K_ESA_JET[k][i];
	dtmp1 = DEL_Tph*c * K_ABS_SSC_JET[k][i];
      }
      // The slice is evolving during the photons passage
      else{
	dtmp = DEL_Tj*c * K_ESA_JET[k][i];
	dtmp1 = DEL_Tj*c * K_ABS_SSC_JET[k][i];
	DEL_Tph -= DEL_Tj;
      }
      // thick
      if (dtmp > 7.00e+2) I_rad_syn[i] = 0.0;
      
      if (dtmp1 > 7.00e+2){
	I_rad_com[i] = 0.0;
	I_rad2nd[i] = 0.0;
	I_eic_jet[i] = 0.0;
	I_com_ext[i] = 0.0;
        I_com_ext1[i] = 0.0;
      }
      // transparent
      //if (dtmp < 1.0e-10)                       I_rad[i] = I_rad[i] * 1.0;      
      // thin
      if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) I_rad_syn[i] = I_rad_syn[i] * exp(-dtmp);
      
      if ((dtmp1 < 7.00e+2) && (dtmp1 > 1.0e-10)){
	I_rad_com[i] = I_rad_com[i] * exp(-dtmp1);
	I_rad2nd[i] = I_rad2nd[i] * exp(-dtmp1);
	I_eic_jet[i] = I_eic_jet[i] * exp(-dtmp1);
	I_com_ext[i] = I_com_ext[i] * exp(-dtmp1);
        I_com_ext1[i] = I_com_ext1[i] * exp(-dtmp1);
      }
      }
      k += 1;
      fprintf(stderr, "\b\b\b");
    }
     
     
     
   // synchrotron  
   sprintf(stmp, "./data/%s_ss.dat", prefix);
   if (ifexist(stmp)) {
     sprintf(comma, "mv ./data/%s_ss.dat ./data/%s_prev_ss.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
   stream_dat = fopen(stmp, "w+");  
   
   // transformation to observer frame     
   for (i = 1; i <= NU_DIM; i++) {
      nu_tmp = FreqTransS2O(NU[i] * DOP_B_J, DOP_B/DOP_B_J , z);
      I_syn = I_rad_syn[i];
      fx_tmp = Intens2Flux(I_syn, R_src, DOP_B, z, H_0) / pow((DOP_B_J),3.0);

	
      if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f\n", 
                log10(nu_tmp), 
                log10(fx_tmp), 
                log10(nu_tmp*fx_tmp));
      }                   
   }     
   fclose(stream_dat); 
   
   
   //SSC
   sprintf(stmp, "./data/%s_cs.dat", prefix);
   if (ifexist(stmp)) {
     sprintf(comma, "mv ./data/%s_cs.dat ./data/%s_prev_cs.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
   stream_dat = fopen(stmp, "w+");  
   
   // transformation to observer frame      
   for (i = 1; i <= NU_DIM; i++) {
     nu_tmp = FreqTransS2O(NU[i]*DOP_B_J, DOP_B/DOP_B_J , z);
     I_com = I_rad_com[i];
     fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0) / pow((DOP_B_J),3.0);

	
     if(IIR_level==1){
        if(EBLFLAG==0){
           if (z <= 0.1){
	      tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	   }else{
	      tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	   }
        }else if(EBLFLAG==1){
           tt=tau_IRA_Franceschini(nu_tmp, z);
	}else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
        }else if(EBLFLAG==3){
            tt=tau_IRA_Franceschini17(nu_tmp, z);
        }
     }else tt=0.;
      
       
     dtmp   = fx_tmp;
     // absorption by IIR
     fx_tmp = fx_tmp * exp(-tt); 
      
     if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
             log10(nu_tmp), 
             log10(fx_tmp), 
             log10(nu_tmp*fx_tmp),
             log10(dtmp),
             log10(nu_tmp*dtmp));
      }                  
   }     
   fclose(stream_dat);
   
   
   //SSC_2nd
   sprintf(stmp, "./data/%s_cs2.dat", prefix);
   if (ifexist(stmp)) {
     sprintf(comma, "mv ./data/%s_cs2.dat ./data/%s_prev_cs2.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
   stream_dat = fopen(stmp, "w+");  
   
   // transformation to observer frame      
   for (i = 1; i <= NU_DIM; i++) {
     nu_tmp = FreqTransS2O(NU[i]*DOP_B_J, DOP_B/DOP_B_J , z);
     I_com2 = I_rad2nd[i];
     fx_tmp = Intens2Flux(I_com2, R_src, DOP_B, z, H_0) / pow((DOP_B_J),3.0);
	
     if(IIR_level==1){
        if(EBLFLAG==0){
           if (z <= 0.1){
	      tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	   }else{
	      tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	   }
        }else if(EBLFLAG==1){
           tt=tau_IRA_Franceschini(nu_tmp, z);
	}else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
        }else if(EBLFLAG==3){
            tt=tau_IRA_Franceschini17(nu_tmp, z);
        }       
     }else tt=0.;
      
       
     dtmp   = fx_tmp;
     // absorption by IIR
     fx_tmp = fx_tmp * exp(-tt); 
      
     if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
             log10(nu_tmp), 
             log10(fx_tmp), 
             log10(nu_tmp*fx_tmp),
             log10(dtmp),
             log10(nu_tmp*dtmp));
      }                  
   }     
   fclose(stream_dat);
   
   
   //eic from jet
   sprintf(stmp, "./data/%s_ecs_jet.dat", prefix);
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/%s_ecs_jet.dat ./data/%s_prev_ecs_jet.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
   stream_dat = fopen(stmp, "w+");
   
   // transformation to observer frame      
   for (i = 1; i <= NU_DIM; i++) {
     nu_tmp = FreqTransS2O(NU[i]*DOP_B_J*LOR_B_J, DOP_B/DOP_B_J , z);
     I_com = I_eic_jet[i];
     fx_tmp = M_PI * ((R_src*R_src) / (D_L*D_L)) * (1.0 + z) * (pow(DOP_B,3.0)/pow((DOP_B_J),3.0)) * I_com;
	
     if(IIR_level==1){
        if(EBLFLAG==0){
           if (z <= 0.1){
	      tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	   }else{
	      tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	   }
        }else if(EBLFLAG==1){
           tt=tau_IRA_Franceschini(nu_tmp, z);
	}else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
        }else if(EBLFLAG==3){
             tt=tau_IRA_Franceschini17(nu_tmp, z);
        }        
     }else tt=0.;
      
       
     dtmp   = fx_tmp;
     // absorption by IIR
     fx_tmp = fx_tmp * exp(-tt); 
      
     if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
             log10(nu_tmp), 
             log10(fx_tmp), 
             log10(nu_tmp*fx_tmp),
             log10(dtmp),
             log10(nu_tmp*dtmp));
      }
   }     
   fclose(stream_dat);
   
   
   //eic from disk-BLR
   sprintf(stmp, "./data/%s_ecs.dat", prefix);
   if (ifexist(stmp)) {
      sprintf(comma, "mv ./data/%s_ecs.dat ./data/%s_prev_ecs.dat", prefix, prefix);
      if(system(comma)) return 1;
   }
   stream_dat = fopen(stmp, "w+");
   
   // transformation to observer frame      
   for (i = 1; i <= NU_DIM; i++) {
     nu_tmp = FreqTransS2O(NU[i]*LOR_B*DOP_B_J, DOP_B/DOP_B_J, z);
     I_com = I_com_ext[i];
     fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0) / pow((DOP_B_J),3.0);

     if(IIR_level==1){
        if(EBLFLAG==0){
           if (z <= 0.1){
	      tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	   }else{
	      tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	   }
        }else if(EBLFLAG==1){
           tt=tau_IRA_Franceschini(nu_tmp, z);
	}else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
        }else if(EBLFLAG==3){
            tt=tau_IRA_Franceschini17(nu_tmp, z);
        }
     }else tt=0.;
      
       
     dtmp   = fx_tmp;
     // absorption by IIR
     fx_tmp = fx_tmp * exp(-tt); 
      
     if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
             log10(nu_tmp), 
             log10(fx_tmp), 
             log10(nu_tmp*fx_tmp),
             log10(dtmp),
             log10(nu_tmp*dtmp));
      }                  
   }     
   fclose(stream_dat);
   
   
   //eic from torus-BLR
    sprintf(stmp, "./data/%s_ecs1.dat", prefix);
    if (ifexist(stmp)) {
        sprintf(comma, "mv ./data/%s_ecs1.dat ./data/%s_prev_ecs1.dat", prefix, prefix);
        if(system(comma)) return 1;
    }
    stream_dat = fopen(stmp, "w+");
   
   // transformation to observer frame      
   for (i = 1; i <= NU_DIM; i++) {
     nu_tmp = FreqTransS2O(NU[i]*LOR_B*DOP_B_J, DOP_B/DOP_B_J, z);
     I_com = I_com_ext1[i];
     fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0) / pow((DOP_B_J),3.0);

     if(IIR_level==1){
        if(EBLFLAG==0){
           if (z <= 0.1){
	      tt = tau_IRA_Kneiske(nu_tmp, z, 0);
	   }else{
	      tt = tau_IRA_Kneiske(nu_tmp, z, 1);
	   }
        }else if(EBLFLAG==1){
           tt=tau_IRA_Franceschini(nu_tmp, z);
	}else if(EBLFLAG==2){
             tt=tau_IRA_Finke(nu_tmp, z);
        }else if(EBLFLAG==3){
             tt=tau_IRA_Franceschini17(nu_tmp, z);
        }        
     }else tt=0.;
      
       
     dtmp   = fx_tmp;
     // absorption by IIR
     fx_tmp = fx_tmp * exp(-tt); 
      
     if (fx_tmp > 1.0e-300) {
        fprintf(stream_dat, "%f %f %f %f %f\n", 
             log10(nu_tmp), 
             log10(fx_tmp), 
             log10(nu_tmp*fx_tmp),
             log10(dtmp),
             log10(nu_tmp*dtmp));
      }                  
   }     
   fclose(stream_dat);
   
   
   
   fprintf(stderr, " DONE\n\n");
   fprintf(stderr, "Jet slice where the blob takes place:  %6.3d\n",j-1);
   }
   
   
      
      
      


   
   
   
   
   fprintf(stderr, "\n\nDERIVED PARAMETERS FROM THE BLOB:\n");
   fprintf(stderr, "----------\n");    
   fprintf(stderr, "Velocity:                    %6.3e c \n", V_B/c);
   Tcool_synch = 3.0 * m_e * c /(4 * sig_T * GAMMA_MAX  * Ub_B);
   Tcool_vhe = 3.0 * m_e * c /(4 * sig_T * GAMMA_MAX  * Ub_syn);
   Tcool_radio = 3.0 * m_e * c /(4 * sig_T * GAMMA_MIN  * Ub_syn);
   if (Tcool_synch < Tcool_vhe){
       dtmp = 3.0 * m_e * c /(4 * sig_T * GAMMA_BRK  * Ub_B);
   } else{
       dtmp = 3.0 * m_e * c /(4 * sig_T * GAMMA_BRK  * Ub_syn);
   }
   fprintf(stderr, "Intrinsic synch Gamma_max cool. time:     %6.3e s \n", Tcool_synch);
   fprintf(stderr, "Intrinsic IC Gamma_max cool. time:        %6.3e s \n", Tcool_vhe);
   fprintf(stderr, "Intrinsic IC Gamma_min cool. time:        %6.3e s \n", Tcool_radio);
   fprintf(stderr, "Intrinsic Gamma_break cool. time:         %6.3e s \n", dtmp);
   //fprintf(stderr, "Intrinsic adiabatic cool. time:           %6.3e s \n", R_src/V_B); //not true, this is a lower limit
   fprintf(stderr, "Observed minimal variability:             %6.3e h \n", R_src*(1+z)/(c*DOP_B * 3600.));
   fprintf(stderr, "Electron energy density:     %6.3e cm^-3 \n", elec_spec_tot_b);
   fprintf(stderr, "Synchrotron energy density:  %6.3e erg cm^-3 \n", Ub_syn);
   fprintf(stderr, "Magnetic energy density:     %6.3e erg cm^-3 \n", Ub_B);
   fprintf(stderr, "U_B/U_e:                     %6.3e  \n\n\n", Ub_B/Ub_e);
   Lb_r = M_PI*R_src*R_src*LOR_B*LOR_B*c*(Ub_syn + Ub_ssc + Ub_ssc2 + Ub_eicd + Ub_eicj);
   Lj_r = Pj_syn + Pj_ssc;  
   fprintf(stderr, "Magnetic power:              %6.3e erg s^-1 \n", Lb_B); 
   fprintf(stderr, "Electrons power:             %6.3e erg s^-1 \n", Lb_e);
   fprintf(stderr, "Syn power:                   %6.3e erg s^-1 \n", M_PI*R_src*R_src*LOR_B*LOR_B*c*Ub_syn);
   fprintf(stderr, "Ssc power:                   %6.3e erg s^-1 \n", M_PI*R_src*R_src*LOR_B*LOR_B*c*Ub_ssc);
   fprintf(stderr, "Ssc2 power:                  %6.3e erg s^-1 \n", M_PI*R_src*R_src*LOR_B*LOR_B*c*Ub_ssc2);
   fprintf(stderr, "Eic disk power:              %6.3e erg s^-1 \n", M_PI*R_src*R_src*LOR_B*LOR_B*c*Ub_eicd);
   fprintf(stderr, "Eic jet power:               %6.3e erg s^-1 \n", M_PI*R_src*R_src*LOR_B*LOR_B*c*Ub_eicj);  
   fprintf(stderr, "Radiation power:             %6.3e erg s^-1 \n", Lb_r); 
   

   if(CASE_JET){
   fprintf(stderr, "\n\nDERIVED PARAMETERS FROM THE JET:\n");
   fprintf(stderr, "----------\n\n"); 
   fprintf(stderr, "velocity of jet:             %e [in c]\n", V_JET/c);
   fprintf(stderr, "Lorentz factor of jet:       %e\n", LOR_JET);
   fprintf(stderr, "apparent velocity of jet     %e [in c]\n", V_JET_APP / c);
   fprintf(stderr, "Minimal variability:         %6.3e h \n", Y_MIN*(1+z)/(c*DOP_JET * 3600.));
   fprintf(stderr, "outer radius of jet:         %e [cm] (%f [pc])\n", Y_MAX, Y_MAX/pc);
   fprintf(stderr, "theta in the jet frame:      %e [deg]\n", THETA_src_j);
   fprintf(stderr, "Magnetic energy density:     %6.3e erg cm^-3 \n", Utot_B);
   fprintf(stderr, "Electron energy density:     %6.3e erg cm^-3 \n", Utot_e);
   fprintf(stderr, "Proton energy density:       %6.3e erg cm^-3 \n", Utot_p);
   fprintf(stderr, "U_B/U_e:                     %6.3e  \n\n\n", Utot_B/Utot_e);   
   fprintf(stderr, "Magnetic power:              %6.3e erg s^-1 \n", Lj_B);   
   fprintf(stderr, "Electrons power:             %6.3e erg s^-1 \n", Lj_e);
   fprintf(stderr, "Syn power:                   %6.3e erg s^-1 \n", Pj_syn);
   fprintf(stderr, "Ssc power:                   %6.3e erg s^-1 \n", Pj_ssc);
   fprintf(stderr, "Radiation power:             %6.3e erg s^-1 \n", Lj_r);
   }
   
   fprintf(stderr, "\n\nGENERAL DERIVED PARAMETERS:\n");
   fprintf(stderr, "----------\n");
   fprintf(stderr, "luminosity distance:                                %e [M pc]\n", D_L / (1.0e+6 * pc));
   fprintf(stderr, "Total radiation power (nucleus included):           %6.3e erg s^-1 \n", Lb_r + Lj_r + L_nuc);
   fprintf(stderr, "Total power (without cold particles):               %6.3e erg s^-1 \n", Lb_r + Lj_r + L_nuc + Lb_B + Lj_B + Lb_e + Lj_e);
   fprintf(stderr, "Total power (without cold particles & nucleus):     %6.3e erg s^-1 \n", Lb_r + Lj_r + Lb_B + Lj_B + Lb_e + Lj_e);
   


   
   
   
   
   
   
   
   fprintf(stderr, "\nALL CALCULATIONS DONE !!! ");
   
   T_END = time(NULL);            
   
   if (difftime(T_END, T_START) > 60) 
     sprintf(stmp, "[time: %4.1f min]\n", difftime(T_END, T_START) / 60.0);
   else
     sprintf(stmp, "[time: %4.1f sec]\n", difftime(T_END, T_START));
   
   fprintf(stderr, "%s\n", stmp);
   
   fprintf(stderr, "SAVING FINAL SPECTRUM\n");
   
   return 1;
}





int load_params(char* name) {   
    int DEBUG = 0;
    char   stmp[256];
    double dtmp;
    FILE   *stream_par;
    errno  = 0;
    stream_par = fopen(name, "r");
   
      

    // General transformation parameters

    if (errno == 0) {
      
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);      
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);  

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 0.0) && (dtmp <= 6.0)) z = dtmp; 
      else {
        fprintf(stderr, "wrong value of redshift: %e (0..6) !!!\n", dtmp); 
        return 0; 
      } 
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 50.0) && (dtmp <= 100.0)) H_0 = dtmp; 
      else {
        fprintf(stderr, "wrong value of Hubble constant: %e (50...100 [km/s/Mpc]) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-03) && (dtmp <= 89.9)) THETA = dtmp;
      else {
        fprintf(stderr, "wrong value of angle: %e (1e-3...89.9 [degrees]) !!!\n", dtmp); 
        return 0; 
      }

      // Blob parameters
      
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);      
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0) && (dtmp <= 100)) DOP_B = dtmp;
      else {
        fprintf(stderr, "wrong value of Doppler factor: %e (1...100) !!!\n", dtmp);
        return 0; 
      }  

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-20) && (dtmp <= 1.0e+20)) K1 = dtmp; 
      else {
        fprintf(stderr, "wrong value of particle density: %e (1.0e-20...1.0e+20 [1/cm^3]) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0) && (dtmp <= 8.0)) N1 = dtmp; 
      else {
        fprintf(stderr, "wrong value of first slope\nfor particle spectrum: %e (1...8) !!!\n", dtmp); 
        return 0; 
      }
      //
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0) && (dtmp <= 8.0)) N2 = dtmp; 
      else {
        fprintf(stderr, "wrong value of second slope\nfor particle spectrum: %e (1...8) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0) && (dtmp <= 1.0e+5)) GAMMA_MIN = dtmp;
      else {
        fprintf(stderr, "wrong value of minimum\nelectrons energy: %e (1.0...1.0e+4) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+3) && (dtmp <= 1.0e+9)) GAMMA_MAX = dtmp;
      else {
        fprintf(stderr, "wrong value of maximum\nelectrons energy: %e (1.0e+3...1.0e+9) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= GAMMA_MIN) && (dtmp <= GAMMA_MAX)) GAMMA_BRK = dtmp;
      else {
        fprintf(stderr, "wrong value of break position in\nelectrons energy spectrum: %e (%e...%e) !!!\n", dtmp, GAMMA_MIN, GAMMA_MAX); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-5) && (dtmp <= 1.0e+5)) B = dtmp; 
      else {
        fprintf(stderr, "wrong value of magnetic field: %e (1.0e-5...1.0e+5 [G]) !!!\n", dtmp); 
        return 0; 
      }
            
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+05) && (dtmp <= 1.0e+25)) R_src = dtmp; 
      else {
        fprintf(stderr, "wrong value of radius of\nemitting region: %e (1.0e+05...1.0+25 [cm]) !!!\n", dtmp); 
        return 0; 
      }    
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 0.0) && (dtmp <= 1.0e+17)) L_src = dtmp; 
      else {
        fprintf(stderr, "wrong value of length of\nemitting region: %e (0.0...1.0+17 [cm]) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 0.0) && (dtmp <= 1.0)) IIR_level = (int)(dtmp); 
      else {
        fprintf(stderr, "wrong value of absorption level by\ninfrared intergalactic background: %e (0...1) !!!\n", dtmp); 
        return 0; 
      }

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+10) && (dtmp <= 1.0e+21)) D_b = dtmp; 
      else {
        fprintf(stderr, "wrong value of blob - SMBH distance : %e (1.0e+10...1.0e+21 [cm]) !!!\n", dtmp); 
        return 0; 
      }
      
      


      // Extern Inverse Compton parameters//

      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);      
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp == 0) || (dtmp == 1)) CASE_EIC = (int)dtmp; 
      else {
        fprintf(stderr, "wrong value of case EIC: %e (0 or 1) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp == 0) || (dtmp == 1)) CASE_X = (int)dtmp; 
      else {
        fprintf(stderr, "wrong value of case X: %e (0 or 1) !!!\n", dtmp); 
        return 0; 
      }

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 10.0) && (dtmp <= 1.0e+6)) T_BB = dtmp; 
      else {
        fprintf(stderr, "wrong value of disc black body\ntemperature: %e (10.0...1.0e+6 [K]) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 10.0) && (dtmp <= 1.0e+6)) T_BB_tor = dtmp; 
      else {
        fprintf(stderr, "wrong value of tore black body\ntemperature: %e (10.0...1.0e+6 [K]) !!!\n", dtmp); 
        return 0; 
      }

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+20) && (dtmp <= 1.0e+50)) L_nuc = dtmp; 
      else {
        fprintf(stderr, "wrong value of luminosity\nof the nucleus: %e (1.0e+20...1.0e+50 [erg/s]) !!!\n", dtmp); 
        return 0; 
      }
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-10) && (dtmp <= 1.0)) tau = dtmp; 
      else {
        fprintf(stderr, "wrong value of the fraction of L_nuc scattered/reprocessed isotropically: %e (1.0e-10...1) !!!\n", dtmp); 
        return 0; 
      }
      
            if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+20) && (dtmp <= 1.0e+50)) L_tor = dtmp; 
      else {
        fprintf(stderr, "wrong value of luminosity\nof the tore: %e (1.0e+20...1.0e+50 [erg/s]) !!!\n", dtmp); 
        return 0; 
      }
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-10) && (dtmp <= 1.0)) tau = dtmp; 
      else {
        fprintf(stderr, "wrong value of the fraction of L_tore scattered/reprocessed isotropically: %e (1.0e-10...1) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e10) && (dtmp <= 1.0e22)) R_blr = dtmp; 
      else {
        fprintf(stderr, "wrong value of the R_blr: %e (1.0e10...1.0e22) !!!\n", dtmp); 
        return 0; 
      }


      // Jet parameters

      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);      
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp == 0) || (dtmp == 1)) CASE_JET = (int)dtmp;
      else {
        fprintf(stderr, "wrong value of case JET: %e (0 or 1) !!!\n", dtmp); 
        return 0; 
      }

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0) && (dtmp <= 100)) DOP_JET = dtmp; 
      else {
        fprintf(stderr, "wrong value of jet Doppler factor: %e (1...100) !!!\n", dtmp); 
        return 0; 
      }    

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-5) && (dtmp <= 1.0e+5)) N_0 = dtmp; 
      else {
        fprintf(stderr, "wrong value of initial particle density: %e (1.0e-5...1.0e+5 [1/cm^3]) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0) && (dtmp <= 5.0)) n_n = dtmp; 
      else {
        fprintf(stderr, "wrong value of particle spectrum slope: %e (1...5) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0) && (dtmp <=  1.0e+5)) GAMMA_MIN1 = dtmp;
      else {
        fprintf(stderr, "wrong value of initial minimum electrons energy: %e (1.0...%e) !!!\n", dtmp, GAMMA_MIN1);
        return 0; 
      }      

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= GAMMA_MIN1) && (dtmp <= 1.0e+7)) GAMMA_MAX_0 = dtmp;
      else {
        fprintf(stderr, "wrong value of initial maximum electrons energy: %e (1.0e+3...1.0e+7) !!!\n", dtmp); 
        return 0; 
      }

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-5) && (dtmp <= 1.0e+5)) B_0 = dtmp; 
      else {
        fprintf(stderr, "wrong value of initial magnetic field: %e (1.0e-5...1.0e+5 [G]) !!!\n", dtmp); 
        return 0; 
      }

              
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+10) && (dtmp <= 1.0e+19)) Y_MIN = dtmp; 
      else {
        fprintf(stderr, "wrong value of inner jet radius: %e (1.0e+10...1.0+19 [cm]) !!!\n", dtmp); 
        return 0; 
      }   
      
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);      
      if ((dtmp >= 1.0) && (dtmp <= 1000.0)) J_LEN = dtmp; 
      else {
        fprintf(stderr, "wrong value of jet length: %e (1...1 000 [pc]) !!!\n", dtmp); 
        return 0; 
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e-3) && (dtmp <= 5.0)) PHI = dtmp;  // PHI =  atan(LOR_JET* tan(dtmp*M_PI/180.0)) * 180.0/M_PI;
      else {
        fprintf(stderr, "wrong value of jet opening angle: %e (1.0e-3...5.0 [deg]) !!!\n", dtmp); 
        return 0; 
      }
            
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 5.0) && (dtmp <= SL_DIM_MAX)) SL_DIM = (int)(dtmp);
      else {
        fprintf(stderr, "wrong number of slices: %d (5...%d) !!!\n", (int)(dtmp), SL_DIM_MAX); 
        return 0; 
      }

      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);      
      if(!fscanf(stream_par, "%s\n", stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s\n", stmp);
       
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 10.0) && (dtmp <= NU_DIM_MAX)) NU_DIM = (int)(dtmp);
      else {
        fprintf(stderr, "wrong number of spectral points: %d (10...%d [Hz]) !!!\n", (int)(dtmp), NU_DIM_MAX); 
        return 0; 
      }

      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+7) && (dtmp <= 1.0e+15)) NU_STR = dtmp;
      else {
        fprintf(stderr, "wrong value of minimum frequency: %e (1.0e+7...1.0e+15 [Hz]) !!!\n", dtmp);
        return 0;
      }
      
      if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
      if ((dtmp >= 1.0e+11) && (dtmp <= 1.0e+35)) NU_END = dtmp;
      else {
        fprintf(stderr, "wrong value of maximum frequency: %e (1.0e+11...1.0e+35 [Hz]) !!!\n", dtmp);
        return 0; 
      }

      if(!fscanf(stream_par, "%s  %s\n", prefix, stmp)) return 1;
      if (DEBUG) fprintf(stderr, "%s  %s\n", prefix, stmp);
      fprintf(stderr, "prefix name: %s \n", prefix);
      
      fclose(stream_par);
      fprintf(stderr, "PARAMETERS: '%s' LOADED !\n", name);
      
      return 1;
      
    } else {
      fprintf(stderr, "ERROR (%d) %s\n", errno, strerror(errno));
      fprintf(stderr, "PARAMETERS: '%s' NOT LOADED !!!\n", name);
      return 0;
    }
    return 0;      
}
