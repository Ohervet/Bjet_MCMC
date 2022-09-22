// Slightly modified version of bj02.cpp for use in MCMC code; modified by Sarah Youngquist, Feb 2022 (see below)

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


/*
Instructions:
Calling the executable:
- Usage 1: bj_mcmc --help
- Usage 2: bj_mcmc <parameter file>
- Usage 3: bj_mcmc 0 <parameter file>                             same execution as usage 2, _prev files made
- Usage 4: bj_mcmc 1 <parameter file>                             no _prev files made
- Usage 5: bj_mcmc 2 <data folder> <parameter file>               _prev files made
- Usage 6: bj_mcmc 3 <data folder> <parameter file>               no _prev files made
- Usage 7: bj_mcmc 0/1 <model type> <model params, at least 19 depending on model>
    ^params from command line, 0 = yes _prev files, 1 = no
- Usage 8: bj_mcmc 2/3 <data folder> <model type> <model params, at least 19 depending on model>
    ^params from command line, 2 = yes _prev files, 3 = no

Modes:
- 0: _prev files are made, save in default data folder
- 1: no _prev files are made, save in default data folder
- 2: _prev files are made, save in specified data folder
- 3: no _prev files are made, save in specified data folder

(default data directory is in the same directory as the executable, named "data")

Valid values for model type:
- 0: model with just blob
- 1: blob + EIC parameters
- Other models not yet implemented

Order of parameters in command line:
Model type 0:
 bj_mcmc <0, 1, 2, 3> <data folder if applicable> 0 (model type)
 z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src,
 L_src, IIR_level, D_b, NU_DIM, NU_STR, NU_END, prefix

 argc should be 22 or 23 depending on if data folder is listed

Model type 1:
 bj_mcmc <0, 1, 2, 3> <data folder if applicable> 1 (model type)
 z, H_0, THETA, DOP_B, K1, N1, N2, GAMMA_MIN, GAMMA_MAX, GAMMA_BRK, B, R_src,
 L_src, IIR_level, D_b, T_BB, TBB_tor, L_nuc, tau, L_tor, tau, NU_DIM, NU_STR, NU_END, prefix
 *Note that tau is present twice, this is a slight error in the bjet code. The second tau value is not used for
 anything, but it must be inputted.

 argc should be 28 or 29 depending on if data folder is listed

Example:
 ./bj_mcmc 3 /Users/sed_calculations 1 0.34 69.6 0.57 50.0 612.1 2.28 3.74 2816.9 1803000 44806 0.00236 5.94e+17 0 1 3.8e+15 2013 2.0e+4 1.7e+21 1.5e-10 5.5e+20 9.0e-5 99 50000000.0 1e+29 run

 ^ here, the 3 indicates that the data folder is specified and no prev file is made. 1 is the EIC model type. Then 0.34 is z (redshift) and then the rest of the parameters are enumerated.
*/

/*
Changes to the code by Sarah Youngquist in Feb 2022:
  - Allow for _prev files not to be created.
        Specified in the command line with the mode.
        Stored in INPUT_MODE
        Now, any time _prev_*.dat file would be made, it is conditional on INPUT_MODE being 0 (user entered 0 or 2)
  - Allow for specifying the data folder
        This is specified in the command line with the mode.
        Note: relative path to any data file
  - Use the model without a parameter file, inputting parameters from the
    command line.
  - Function load_params_from_list(argv, model_type, index, argc) added.
  - Most of main method moved into run_models
  - In main method, parse the input mode and set the input type, set data directory
     if applicable, parse model type if applicable, call either load_params or load_params_from_list
*/
//ToDo: integral form of smoothbkpowelaw+cutoff and PL+cutoff to allow more spectral shape with energetic calculation

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cerrno>
#include <float.h>


#include "processes_supp_mcmc.h"
#include "bj_mcmc.h"

using namespace std;
using namespace bj_mcmc02;

namespace bj_mcmc02 {
    int INPUT_MODE = 0; // Sarah Youngquist added;
    // 0 = normal
    // 1 = no prev files
    char DATA_DIRECTORY[512] = "data";
    int CASE_JET;
    int CASE_X;
    int CASE_EIC;
    const int EBLFLAG = 3; // 0: Kneiske et al. (2002,2004), 1: Franceschini (2008), 2: Finke (2010), 3: Franceschini (2017)//

// TIME VARIABLES

    time_t T_START, T_END;
    char PARA_FN[256];

// GLOBAL VARIABLES

    int IIR_level = 0;
    int NU_DIM = 50;      // CURRENT NUMBER OF SPECTRAL POINTS
    double NU_STR = 1.0e+10; // START FREQUENCY
    double NU_END = 1.0e+26; // END FREQUENCY

//double GAMMA_MIN1      = 200.0; //Minimal gamma for the jet particles
    double Utot_e = 0.0;
    double Utot_B = 0.0;
    double Utot_p = 0.0;
    double I_BASE = 0.0;
    double I_JET = 0.0;



// PHYSICAL PARAMETERS

    double z;
    double H_0;
    double THETA, D_L;
    double DOP_B, LOR_B, V_B, V_B_APP;//DOP_BB,
    double R_src, R_blr;
    double L_src, L_nuc;
    double B;
    double K1;
    double K2;
    double N1;
    double N2;
    double GAMMA_MIN;
    double GAMMA_BRK;
    double GAMMA_MAX;
    double T_BB;
    double tau;  // fraction of L_nuc scattered/rerocessed isotropically (EIC)
    double B_0;
    double N_0;
    double n_n;
    double n_N;
    double n_G;
    double D_b, D_b_src_j, D_BJ;
    double U_B, U_e;
    double jj1, kk1;
    double L_jet_eic;
    double L_tor, T_BB_tor;
    double Tcool_synch, Tcool_vhe, Tcool_radio, Tcool_synch_gbreak, Tcool_SSC_gbreak, Tcool_EIC_gbreak;

// TRANSFORMATION PARAMETERS

    double DOP_JET, LOR_JET, V_JET, V_JET_APP, DOP_B_J, LOR_B_J, V_B_J;

// GEOMETRY PARAMETERS


    int SL_DIM;
    double PHI, PHI_src;
    double THETA_src_j;
    double J_LEN, J_LEN_src;
    double A;
    double X_MIN;
    double X_MAX;
    double X_OUT;
    double Y_MIN;
    double Y_MAX;
    double R_MOY;
    double Xcut;
    double antisym;
    double Sprev;
    double Sum_S;

//const int           SL_DIM_MAX = 200;//

// OTHER PARAMETERS
    int SL_CUR;
    double GAMMA_MIN1;
    double GAMMA_MAX_0;
    int null0;
    double V_exp, DEL_Tph, DEL_Tj;
    char prefix[256];


// VECTORS

    double NU[NU_MAX + 1];
    double I_rad[NU_MAX + 1];
    double I_rad1st[NU_MAX + 1]; // for 2nd order SSC
    double I_CMB[NU_MAX + 1];

    double X_VAL[SL_DIM_MAX + 2];
    double Y_VAL[SL_DIM_MAX + 2];
    double DEL_X[SL_DIM_MAX + 2];
    double Sum_DEL[SL_DIM_MAX + 2];
    double Stot[SL_DIM_MAX];
    double Sr_TOT[SL_DIM_MAX];

    double B_VAL[SL_DIM_MAX + 1];
    double N_VAL[SL_DIM_MAX + 1];
    double G_VAL[SL_DIM_MAX + 1];
    double UJ_SLICE[SL_DIM_MAX + 1];
    double PJ_SLICE[SL_DIM_MAX + 1];

    double I_SYN_EXT[NU_DIM_MAX + 1];
    double L_BB_nuc[NU_DIM_MAX + 1];
    double I_eic_jet[NU_DIM_MAX + 1];

    double I_rad_syn[NU_MAX + 1];
    double I_rad_com[NU_MAX + 1];
    double I_rad2nd[NU_MAX + 1];
    double I_rad_ext[NU_MAX + 1];
    double I_rad_ext_Int[NU_MAX + 1];
    double I_rad_ext1[NU_MAX + 1];
    double I_com_ext[NU_MAX + 1];
    double I_com_ext1[NU_MAX + 1];
    double I_com_disc[NU_MAX + 1];
    double I_rad_ext_s[NU_MAX + 1];

    double I_rad_ext_D[200 + 1];
    double D[200 + 1];

    double C_e[NU_MAX + 1];
    double C_a[NU_MAX + 1];
    double C_e1[NU_MAX + 1];
    double C_a1[NU_MAX + 1];

    double F_IC_disk[NU_MAX + 1];
    double F_IC_tot[NU_MAX + 1];
    double NU_IC_disk[NU_MAX + 1];




// MATRIXES

    double I_SYN_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double I_SYN_JET_EDGE[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double I_SYN_JET_BASE[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double J_SYN_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double K_ESA_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double F_SYN_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];

    double I_COM_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double I_COM_JET_EDGE[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double I_COM_JET_BASE[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double J_COM_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double K_ABS_SSC_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double F_COM_JET[SL_DIM_MAX + 1][NU_DIM_MAX + 1];

    double I_SYN_TOT[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double I_COM_TOT[SL_DIM_MAX + 1][NU_DIM_MAX + 1];
    double F_SYN_BLOB[SL_DIM_MAX + 1][NU_DIM_MAX + 1];

    double Sl[SL_DIM_MAX + 1][SL_DIM_MAX + 1];
    double Sr[SL_DIM_MAX + 1][SL_DIM_MAX + 1];
    double S[SL_DIM_MAX + 1][SL_DIM_MAX + 1];
    double Dcut[SL_DIM_MAX + 1][SL_DIM_MAX + 1];



// SUBROUTINES ******************************************************************

    int ifexist(char *fname) {
      FILE *stream;
      errno = 0;
      stream = fopen(fname, "r");
      if (errno == 0) {
        fclose(stream);
        return 1;
      } else {
        return 0;
      }
    }

// ELECTRON SPECTRUM BLOB
    double N_e_BknPowLaw(double gamma) { // Broken power-law (possibly with exp cutoff at high energies)
      K2 = K1 * pow(GAMMA_BRK, N2 - N1);

      if (gamma < GAMMA_BRK) {
        return K1 * pow(gamma, -N1);
      }
      if (gamma >= GAMMA_BRK) {
        return K2 * pow(gamma, -N2);// * exp(-gamma/GAMMA_MAX);
      }

      return 1.0e-100;
    }

    double
    N_e_SmoothBknPowLaw(double gamma) { //cf. Tavecchio et al. 2001: 2001ApJ...554..725T, + exp cutoff at high energies
      return (K1 * pow(gamma, -N1) * pow((1. + gamma / GAMMA_BRK), (N1 - N2))) * exp(-gamma / GAMMA_MAX);
    }

    double N_e_PileUp(double gamma) {
      return (K1 * gamma * gamma * exp(-2.0 * gamma / GAMMA_BRK));
    }


// Switch between different electron spectrum shapes
    double N_e(double gamma) {
      double foo = 0.;
      //foo=N_e_PileUp(gamma);
      //foo=N_e_SmoothBknPowLaw(gamma);
      foo = N_e_BknPowLaw(gamma);
      return foo;
    }


// ELECTRON SPECTRUM JET
    double N_e_PowLaw(double gamma) { //simple power law
      return N_VAL[SL_CUR] * pow(gamma, -n_n);
    }

    double N_e_PowLawExpcut(double gamma) { //Power law with exponential cutoff
      return N_VAL[SL_CUR] * pow(gamma, -n_n) * exp(-gamma / GAMMA_MAX_0);
    }

/*
double N_e_Jet(double gg) {
      if ((gg >= GAMMA_MIN1) && (G_VAL[SL_CUR])) return N_VAL[SL_CUR] * pow(gg, -n_n);
      return 0.0;
}
*/
// Switch between different electron spectrum shapes for the jet
    double N_e_Jet(double gamma) {
      double foo = 0.;
      if ((gamma >= GAMMA_MIN1) && (G_VAL[SL_CUR])) {
        foo = N_e_PowLaw(gamma);
        //foo = N_e_PowLawExpcut(gamma);
      }
      return foo;
    }


    double ftr(double F) {
      return F * pow(DOP_JET, 3) * (1.0 + z);
    }

//ANALYTICAL INTEGRATION OVER THE ELECTRON SPECTRUM = integral(Ne(E) dE) or integral(E*Ne(E) dE)

// simple power law (jet)
    double int_spec_j(double N_0, int a) {
      //use a = 1 for particle density, a = 2 for energy density
      if (n_n == a) {
        return N_0 * (log(GAMMA_MAX_0) - log(GAMMA_MIN1));
      } else {
        return N_0 / (a - n_n) * (pow(GAMMA_MAX_0, a - n_n) - pow(GAMMA_MIN1, a - n_n));
      }
    }

// boken power-law (blob)
    double int_spec_b(int a) {
      //use a = 1 for particle density, a = 2 for energy density
      if (N1 == a and N2 == a) {
        return K1 * ((log(GAMMA_BRK) - log(GAMMA_MIN)) + pow(GAMMA_BRK, N2 - N1) * (log(GAMMA_MAX) - log(GAMMA_BRK)));
      }
      if (N1 == a and N2 != a) {
        return K1 * ((log(GAMMA_BRK) - log(GAMMA_MIN)) +
                     pow(GAMMA_BRK, (N2 - N1)) / (a - N2) * (pow(GAMMA_MAX, (a - N2)) - pow(GAMMA_BRK, a - N2)));
      }
      if (N1 != a and N2 == a) {
        return K1 * (1 / (a - N1) * (pow(GAMMA_BRK, a - N1) - pow(GAMMA_MIN, a - N1)) +
                     pow(GAMMA_BRK, N2 - N1) * (log(GAMMA_MAX) - log(GAMMA_BRK)));
      } else {
        return K1 * (1 / (a - N1) * (pow(GAMMA_BRK, a - N1) - pow(GAMMA_MIN, a - N1)) +
                     pow(GAMMA_BRK, (N2 - N1)) / (a - N2) * (pow(GAMMA_MAX, a - N2) - pow(GAMMA_BRK, a - N2)));
      }
    }


    // METHODS ========================================================================================================
    void description() {
      cout
          << "This code computes the radiative output from a homogeneous blob and a stratified jet given ther particle energy distribution"
          << endl << "The processes accounted for are:" << endl << endl << " - synchrotron radiation (blob and jet)"
          << endl << " - SSC radiation (blob and jet)" << endl << " - 2nd order SSC radiation (blob)" << endl
          << " - external inverse Compton on a disk/corona/dust torus radiation field" << endl
          << " - external inverse Compton on the CMB" << endl
          << " - external inverse Compton on the jet synchrotoron radiation field" << endl << endl;
    }

    int load_params(char *name) {
      // from parameter file
      int DEBUG = 0;
      char stmp[512];
      double dtmp;
      FILE *stream_par;
      errno = 0;
      stream_par = fopen(name, "r");


      // General transformation parameters

      if (errno == 0) {

        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 0.0) && (dtmp <= 6.0)) z = dtmp;
        else {
          fprintf(stderr, "wrong value of redshift: %e (0..6) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 50.0) && (dtmp <= 100.0)) H_0 = dtmp;
        else {
          fprintf(stderr, "wrong value of Hubble constant: %e (50...100 [km/s/Mpc]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-03) && (dtmp <= 89.9)) THETA = dtmp;
        else {
          fprintf(stderr, "wrong value of angle: %e (1e-3...89.9 [degrees]) !!!\n", dtmp);
          return 0;
        }

        // Blob parameters

        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 100)) DOP_B = dtmp;
        else {
          fprintf(stderr, "wrong value of Doppler factor: %e (1...100) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-20) && (dtmp <= 1.0e+20)) K1 = dtmp;
        else {
          fprintf(stderr, "wrong value of particle density: %e (1.0e-20...1.0e+20 [1/cm^3]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 8.0)) N1 = dtmp;
        else {
          fprintf(stderr, "wrong value of first slope\nfor particle spectrum: %e (1...8) !!!\n", dtmp);
          return 0;
        }
        //
        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 8.0)) N2 = dtmp;
        else {
          fprintf(stderr, "wrong value of second slope\nfor particle spectrum: %e (1...8) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 1.0e+5)) GAMMA_MIN = dtmp;
        else {
          fprintf(stderr, "wrong value of minimum\nelectrons energy: %e (1.0...1.0e+4) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+3) && (dtmp <= 1.0e+9)) GAMMA_MAX = dtmp;
        else {
          fprintf(stderr, "wrong value of maximum\nelectrons energy: %e (1.0e+3...1.0e+9) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= GAMMA_MIN) && (dtmp <= GAMMA_MAX)) GAMMA_BRK = dtmp;
        else {
          fprintf(stderr, "wrong value of break position in\nelectrons energy spectrum: %e (%e...%e) !!!\n", dtmp,
                  GAMMA_MIN, GAMMA_MAX);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-5) && (dtmp <= 1.0e+5)) B = dtmp;
        else {
          fprintf(stderr, "wrong value of magnetic field: %e (1.0e-5...1.0e+5 [G]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+05) && (dtmp <= 1.0e+25)) R_src = dtmp;
        else {
          fprintf(stderr, "wrong value of radius of\nemitting region: %e (1.0e+05...1.0+25 [cm]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 0.0) && (dtmp <= 1.0e+17)) L_src = dtmp;
        else {
          fprintf(stderr, "wrong value of length of\nemitting region: %e (0.0...1.0+17 [cm]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 0.0) && (dtmp <= 1.0)) IIR_level = (int) (dtmp);
        else {
          fprintf(stderr, "wrong value of absorption level by\ninfrared intergalactic background: %e (0...1) !!!\n",
                  dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+10) && (dtmp <= 1.0e+21)) D_b = dtmp;
        else {
          fprintf(stderr, "wrong value of blob - SMBH distance : %e (1.0e+10...1.0e+21 [cm]) !!!\n", dtmp);
          return 0;
        }




        // Extern Inverse Compton parameters//

        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp == 0) || (dtmp == 1)) CASE_EIC = (int) dtmp;
        else {
          fprintf(stderr, "wrong value of case EIC: %e (0 or 1) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp == 0) || (dtmp == 1)) CASE_X = (int) dtmp;
        else {
          fprintf(stderr, "wrong value of case X: %e (0 or 1) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 10.0) && (dtmp <= 1.0e+6)) T_BB = dtmp;
        else {
          fprintf(stderr, "wrong value of disc black body\ntemperature: %e (10.0...1.0e+6 [K]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 10.0) && (dtmp <= 1.0e+6)) T_BB_tor = dtmp;
        else {
          fprintf(stderr, "wrong value of tore black body\ntemperature: %e (10.0...1.0e+6 [K]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+20) && (dtmp <= 1.0e+50)) L_nuc = dtmp;
        else {
          fprintf(stderr, "wrong value of luminosity\nof the nucleus: %e (1.0e+20...1.0e+50 [erg/s]) !!!\n", dtmp);
          return 0;
        }
        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-10) && (dtmp <= 1.0)) tau = dtmp;
        else {
          fprintf(stderr,
                  "wrong value of the fraction of L_nuc scattered/reprocessed isotropically: %e (1.0e-10...1) !!!\n",
                  dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+20) && (dtmp <= 1.0e+50)) L_tor = dtmp;
        else {
          fprintf(stderr, "wrong value of luminosity\nof the tore: %e (1.0e+20...1.0e+50 [erg/s]) !!!\n", dtmp);
          return 0;
        }
        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-10) && (dtmp <= 1.0)){} // tau = dtmp;
        else {
          fprintf(stderr,
                  "wrong value of the fraction of L_tore scattered/reprocessed isotropically: %e (1.0e-10...1) !!!\n",
                  dtmp);
          return 0;
        }
        /*
        if(!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e10) && (dtmp <= 1.0e22)) R_blr = dtmp;
        else {
          fprintf(stderr, "wrong value of the R_blr: %e (1.0e10...1.0e22) !!!\n", dtmp);
          return 0;
        }
    */

        // Jet parameters

        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp == 0) || (dtmp == 1)) CASE_JET = (int) dtmp;
        else {
          fprintf(stderr, "wrong value of case JET: %e (0 or 1) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 100)) DOP_JET = dtmp;
        else {
          fprintf(stderr, "wrong value of jet Doppler factor: %e (1...100) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-5) && (dtmp <= 1.0e+7)) N_0 = dtmp;
        else {
          fprintf(stderr, "wrong value of initial particle density: %e (1.0e-5...1.0e+7 [1/cm^3]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 5.0)) n_n = dtmp;
        else {
          fprintf(stderr, "wrong value of particle spectrum slope: %e (1...5) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 1.0e+5)) GAMMA_MIN1 = dtmp;
        else {
          fprintf(stderr, "wrong value of initial minimum electrons energy: %e (1.0...%e) !!!\n", dtmp, GAMMA_MIN1);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= GAMMA_MIN1) && (dtmp <= 1.0e+7)) GAMMA_MAX_0 = dtmp;
        else {
          fprintf(stderr, "wrong value of initial maximum electrons energy: %e (1.0e+3...1.0e+7) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-5) && (dtmp <= 1.0e+5)) B_0 = dtmp;
        else {
          fprintf(stderr, "wrong value of initial magnetic field: %e (1.0e-5...1.0e+5 [G]) !!!\n", dtmp);
          return 0;
        }


        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+10) && (dtmp <= 1.0e+19)) Y_MIN = dtmp;
        else {
          fprintf(stderr, "wrong value of inner jet radius: %e (1.0e+10...1.0+19 [cm]) !!!\n", dtmp);
          return 0;
        }


        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0) && (dtmp <= 1000.0)) J_LEN = dtmp;
        else {
          fprintf(stderr, "wrong value of jet length: %e (1...1 000 [pc]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e-3) && (dtmp <= 5.0)) PHI = dtmp;  // PHI =  atan(LOR_JET* tan(dtmp*M_PI/180.0)) * 180.0/M_PI;
        else {
          fprintf(stderr, "wrong value of jet opening angle: %e (1.0e-3...5.0 [deg]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 5.0) && (dtmp <= SL_DIM_MAX)) SL_DIM = (int) (dtmp);
        else {
          fprintf(stderr, "wrong number of slices: %d (5...%d) !!!\n", (int) (dtmp), SL_DIM_MAX);
          return 0;
        }

        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);
        if (!fscanf(stream_par, "%s\n", stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s\n", stmp);


        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 10.0) && (dtmp <= NU_DIM_MAX)) NU_DIM = (int) (dtmp);
        else {
          fprintf(stderr, "wrong number of spectral points: %d (10...%d [Hz]) !!!\n", (int) (dtmp), NU_DIM_MAX);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+7) && (dtmp <= 1.0e+15)) NU_STR = dtmp;
        else {
          fprintf(stderr, "wrong value of minimum frequency: %e (1.0e+7...1.0e+15 [Hz]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%lf  %s\n", &dtmp, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%e  %s\n", dtmp, stmp);
        if ((dtmp >= 1.0e+11) && (dtmp <= 1.0e+35)) NU_END = dtmp;
        else {
          fprintf(stderr, "wrong value of maximum frequency: %e (1.0e+11...1.0e+35 [Hz]) !!!\n", dtmp);
          return 0;
        }

        if (!fscanf(stream_par, "%s  %s\n", prefix, stmp)) return 0;
        if (DEBUG) fprintf(stderr, "%s  %s\n", prefix, stmp);
        fprintf(stderr, "prefix name: %s \n", prefix);

        fclose(stream_par);
        fprintf(stderr, "PARAMETERS: '%s' LOADED !\n", name);

        return 1;  // everything went well

      } else {
        fprintf(stderr, "ERROR (%d) %s\n", errno, strerror(errno));
        fprintf(stderr, "PARAMETERS: '%s' NOT LOADED !!!\n", name);
        return 0;
      }
      return 0;
    }

    int load_params_from_list(char **list, int model_type, int starting_index, int list_length) {
      // start reading params
      double tmp[12];
      int list_index = starting_index;

      // -------------------------------------------------------------------------------------------------------
      // transformation params--redshift, hubble constant, angle
      for (int i = 0; i < 3; ++i) {
        tmp[i] = atof(list[list_index]);
        list_index++;
      }
      if (tmp[0] >= 0.0 && tmp[0] <= 6.0) {
        z = tmp[0];
      } else {
        fprintf(stderr, "wrong value of redshift: %e (0..6) !!!\n", tmp[0]);
        return 0;
      }

      if (tmp[1] >= 50.0 && tmp[1] <= 100.0) {
        H_0 = tmp[1];
      } else {
        fprintf(stderr, "wrong value of Hubble constant: %e (50...100 [km/s/Mpc]) !!!\n", tmp[1]);
        return 0;
      }

      if (tmp[2] >= 1.0e-03 && tmp[2] <= 89.8)
        THETA = tmp[2];
      else {
        fprintf(stderr, "wrong value of angle: %e (1e-3...89.9 [degrees]) !!!\n", tmp[2]);
        return 0;
      }

      // -------------------------------------------------------------------------------------------------------

      // blob params--doppler factor, particle density, slope1, slope2, gamma min, gamma max, gamma break, magnetic field strength,
      //    radius, emitting length, absorption level, blob distance
      for (int i = 0; i < 10; ++i) {
        tmp[i] = atof(list[list_index]);
        list_index++;
      }

      if (tmp[0] >= 1.0 && tmp[0] <= 100.0)
        DOP_B = tmp[0];
      else {
        fprintf(stderr, "wrong value of Doppler factor: %e (1...100) !!!\n", tmp[0]);
        return 0;
      }

      if (tmp[1] >= 1.0e-20 && tmp[1] <= 1.0e+20)
        K1 = tmp[1];
      else {
        fprintf(stderr, "wrong value of particle density: %e (1.0e-20...1.0e+20 [1/cm^3]) !!!\n", tmp[1]);
        return 0;
      }

      if (tmp[2] >= 1.0 && tmp[2] <= 8.0)
        N1 = tmp[2];
      else {
        fprintf(stderr, "wrong value of first slope\nfor particle spectrum: %e (1...8) !!!\n", tmp[2]);
        return 0;
      }

      if (tmp[3] >= 1.0 && tmp[3] <= 8.0)
        N2 = tmp[3];
      else {
        fprintf(stderr, "wrong value of second slope\nfor particle spectrum: %e (1...8) !!!\n", tmp[3]);
        return 0;
      }

      if (tmp[4] >= 1.0 && tmp[4] <= 1.0e+5)
        GAMMA_MIN = tmp[4];
      else {
        fprintf(stderr, "wrong value of minimum\nelectrons energy: %e (1.0...1.0e+4) !!!\n", tmp[4]);
        return 0;
      }

      if (tmp[5] >= 1.0e+3 && tmp[5] <= 1.0e+9)
        GAMMA_MAX = tmp[5];
      else {
        fprintf(stderr, "wrong value of maximum\nelectrons energy: %e (1.0e+3...1.0e+9) !!!\n", tmp[5]);
        return 0;
      }

      if (tmp[6] >= tmp[4] && tmp[6] <= tmp[5])
        GAMMA_BRK = tmp[6];
      else {
        fprintf(stderr, "wrong value of break position in\nelectrons energy spectrum: %e (%e...%e) !!!\n", tmp[6],
                GAMMA_MIN, GAMMA_MAX);
        return 0;
      }

      if (tmp[7] >= 1.0e-5 && tmp[7] <= 1.0e+5)
        B = tmp[7];
      else {
        fprintf(stderr, "wrong value of magnetic field: %e (1.0e-5...1.0e+5 [G]) !!!\n", tmp[7]);
        return 0;
      }

      if (tmp[8] >= 1.0e+05 && tmp[8] <= 1.0e+25)
        R_src = tmp[8];
      else {
        fprintf(stderr, "wrong value of radius of\nemitting region: %e (1.0e+05...1.0+25 [cm]) !!!\n", tmp[8]);
        return 0;
      }

      if (tmp[9] >= 0.0 && tmp[9] <= 1.0e+17)
        L_src = tmp[9];
      else {
        fprintf(stderr, "wrong value of length of\nemitting region: %e (0.0...1.0+17 [cm]) !!!\n", tmp[9]);
        return 0;
      }

      int tmp_iir = static_cast<int>(atof(list[list_index++]));
      double tmp_dist = atof(list[list_index++]);

      if (tmp_iir >= 0 && tmp_iir <= 1)
        IIR_level = tmp_iir;
      else {
        fprintf(stderr, "wrong value of absorption level by\ninfrared intergalactic background: %d (0...1) !!!\n",
                tmp_iir);
        return 0;
      }

      if (tmp_dist >= 1.0e+10 && tmp_dist <= 1.0e+21)
        D_b = tmp_dist;
      else {
        fprintf(stderr, "wrong value of blob - SMBH distance : %e (1.0e+10...1.0e+21 [cm]) !!!\n", tmp_dist);
        return 0;
      }

      // + NUCLEUS PARAMETERS
      if (model_type == 1) {
        // NUCLEUS PARAMETERS ----------------------------------------------------------------------------------------------
        CASE_EIC = 1;

        double temp = atof(list[list_index++]);
        if (temp >= 10.0 && temp <= 1.0e+6)
          T_BB = temp;
        else {
          fprintf(stderr, "wrong value of disc black body\ntemperature: %e (10.0...1.0e+6 [K]) !!!\n", temp);
          return 0;
        }

        temp = atof(list[list_index++]);
        if (temp >= 10.0 && temp <= 1.0e+6)
          T_BB_tor = temp;
        else {
          fprintf(stderr, "wrong value of tore black body\ntemperature: %e (10.0...1.0e+6 [K]) !!!\n", temp);
          return 0;
        }

        temp = atof(list[list_index++]);
        if (temp >= 1.0e+20 && temp <= 1.0e+50)
          L_nuc = temp;
        else {
          fprintf(stderr, "wrong value of luminosity\nof the nucleus: %e (1.0e+20...1.0e+50 [erg/s]) !!!\n", temp);
          return 0;
        }

        temp = atof(list[list_index++]);
        if (temp >= 1.0e-10 && temp <= 1.0)
          tau = temp;
        else {
          fprintf(stderr,
                  "wrong value of the fraction of L_nuc scattered/reprocessed isotropically: %e (1.0e-10...1) !!!\n",
                  temp);
          return 0;
        }

        temp = atof(list[list_index++]);
        if (temp >= 1.0e+20 && temp <= 1.0e+50)
          L_tor = temp;
        else {
          fprintf(stderr, "wrong value of luminosity\nof the tore: %e (1.0e+20...1.0e+50 [erg/s]) !!!\n", temp);
          return 0;
        }

        temp = atof(list[list_index++]);
        if (temp >= 1.0e-10 && temp <= 1.0) {}
          //tau = temp;
        else {
          fprintf(stderr,
                  "wrong value of the fraction of L_tore scattered/reprocessed isotropically: %e (1.0e-10...1) !!!\n",
                  temp);
          return 0;
        }
      }

      // -------------------------------------------------------------------------------------------------------

      // numerical parameters--number of pts, min freq, max freq, file prefix
      int tmp_pts = atoi(list[list_index++]);
      double tmp_min = atof(list[list_index++]);
      double tmp_max = atof(list[list_index++]);

      if (tmp_pts >= 10 && tmp_pts <= NU_DIM_MAX)
        NU_DIM = tmp_pts;
      else {
        fprintf(stderr, "wrong number of spectral points: %d (10...%d [Hz]) !!!\n", (int) (tmp_pts), NU_DIM_MAX);
        return 0;
      }

      if (tmp_min >= 1.0e+7 && tmp_min <= 1.0e+15)
        NU_STR = tmp_min;
      else {
        fprintf(stderr, "wrong value of minimum frequency: %e (1.0e+7...1.0e+15 [Hz]) !!!\n", tmp_min);
        return 0;
      }

      if (tmp_max >= 1.0e+11 && tmp_max <= 1.0e+35)
        NU_END = tmp_max;
      else {
        fprintf(stderr, "wrong value of maximum frequency: %e (1.0e+11...1.0e+35 [Hz]) !!!\n", tmp_max);
        return 0;
      }

      // get file_prefix
      sprintf(prefix, "%s", list[list_index]);
      // we're done and it worked!
      return 1;
    }

    int run_models() {
      // returns 0 if succeeds
      // returns 1 if system command fails
      // returns 2 if called incorrectly or file reading fails
      int i, j, k, l;
      char stmp[256];
      char comma[256];
      char name[256];
      char stmp1[256], stmp2[256];
      double tmp_min, tmp_max, tmp_stp, tmp_val, tmp_cur;
      double nu_tmp, fx_tmp, nu_tmp1, fx_tmp1;
      double jj, kk, tt, I_syn, I_SYN, I_com, I_com2, I_COM;
      double dtmp, dtmp1, dtmp2, dtmp3, dtmp4;
      double Ub_e, Uj_e, Ub_syn, Ub_ssc, Ub_ssc2, Ub_eicd, Ub_eicj, Ub_B, Ub_blr;
      double Lb_B, Lj_B, Lb_e, Lj_e, Lb_r, Lj_r, Lb_p, Lj_p;
      double Pj_syn, Pj_ssc, elec_spec_tot_b, elec_spec_tot_j;
      FILE *stream_dat;
      FILE *stream_dat1;
      FILE *stream_tau;

      T_START = time(NULL);

      double Gamma[G_DIM + 1];
      double elec_spec[G_DIM + 1];


      D_L = (2.0 * c * (z + 1.0 - sqrt(z + 1.0))) / ((H_0 * 1.0e-19) / 3.086);


      V_B = c * (DOP_B * DOP_B * cos(THETA * M_PI / 180.0) -
                 sqrt(1.0 - DOP_B * DOP_B * sin(THETA * M_PI / 180.0) * sin(THETA * M_PI / 180.0))) /
            (1.0 + DOP_B * DOP_B * cos(THETA * M_PI / 180.0) * cos(THETA * M_PI / 180.0));

      LOR_B = 1.0 / sqrt(1.0 - pow(V_B / c, 2.0));

      //only for NGC1275
      /*
      LOR_B   = 9.86; //Doppler unbeamed case (too large angle)
      V_B = sqrt(1. - 1./(LOR_B*LOR_B)) *c;
      */



      if (DOP_B * sin((THETA) * M_PI / 180.0) > 1.0) {
        fprintf(stderr,
                "we must have DOP_B * sin(THETA) < 1, else the speed of the blob is NOT REAL...\n\nProgram aborted...");
        return 2;
      }

      V_B_APP = ((V_B) * sin(THETA * M_PI / 180)) / (1.0 - (V_B / c) * cos(THETA * M_PI / 180));

      Ub_B = B * B / (8 * M_PI);
      Lb_B = M_PI * R_src * R_src * LOR_B * LOR_B * c * Ub_B;

      if (CASE_EIC) {
        R_blr = 0.1 * sqrt(L_nuc / 1.0e46) * 3.08568e18; //Sikora et al. 2009
      }


      if (CASE_JET) {

        V_JET = c * (2.0 * DOP_JET * DOP_JET * cos(THETA * M_PI / 180.0) - 2.0 * sqrt(
            1.0 - DOP_JET * DOP_JET + DOP_JET * DOP_JET * cos(THETA * M_PI / 180.0) * cos(THETA * M_PI / 180.0))) /
                (2.0 * (1.0 + DOP_JET * DOP_JET * cos(THETA * M_PI / 180.0) * cos(THETA * M_PI / 180.0)));

        LOR_JET = 1.0 / sqrt(1.0 - pow(V_JET / c, 2.0));

        V_JET_APP = ((V_JET) * sin(THETA * M_PI / 180.0)) / (1.0 - (V_JET / c) * cos(THETA * M_PI / 180.0));

        J_LEN_src = J_LEN / LOR_JET;
        PHI_src = atan(LOR_JET * tan(PHI * M_PI / 180.0)) * 180.0 / M_PI;
        D_b_src_j = D_b / LOR_JET;

        //relativistic aberration [deg]
        THETA_src_j =
            acos((cos(THETA * M_PI / 180.0) - V_JET / c) / (1.0 - V_JET / c * cos(THETA * M_PI / 180.0))) * 180.0 /
            M_PI;

        //jet lateral expansion speed (jet comoving frame)
        V_exp = V_JET * tan(PHI_src * M_PI / 180.0);


        // transformation blob to jet frame
        V_B_J = (V_B - V_JET) / (1.0 - (V_B * V_JET / (c * c)));
        LOR_B_J = 1.0 / (sqrt(1.0 - pow(V_B_J / c, 2.0)));
        DOP_B_J = 1.0 / (LOR_B_J * (1.0 - V_B_J / c * cos(THETA * M_PI / 180.0)));

        // distance blob - first slice (jet frame)
        D_BJ = X_MIN - D_b_src_j;
      }
      //
      if (PRINT) {
        fprintf(stderr, "Hubble constant:        %6.3f\n", H_0);
        fprintf(stderr, "redshift                %6.3f\n", z);
        fprintf(stderr, "\nBlob parameters:\n");
        fprintf(stderr, "---------------\n");
        fprintf(stderr, "Doppler factor:         %6.3f\n", DOP_B);
        fprintf(stderr, "theta:                  %6.3e\n", THETA);
        fprintf(stderr, "V blob:                 %6.3e\n", V_B / c);
        fprintf(stderr, "Lorentz factor:         %6.3e\n", LOR_B);
        fprintf(stderr, "V apparent:             %6.3e\n", V_B_APP / c);
        fprintf(stderr, "K_1:                    %6.3e\n", K1);
        fprintf(stderr, "n_1:                    %6.3f\n", N1);
        fprintf(stderr, "n_2:                    %6.3f\n", N2);
        fprintf(stderr, "gamma_min:              %6.3e\n", GAMMA_MIN);
        fprintf(stderr, "gamma_brk:              %6.3e\n", GAMMA_BRK);
        fprintf(stderr, "gamma_max:              %6.3e\n", GAMMA_MAX);
        fprintf(stderr, "radius:                 %6.3e\n", R_src);
        fprintf(stderr, "Distance SMBH-blob (host frame):           %e [cm] (%f [pc])\n", D_b, D_b / pc);
        if (CASE_EIC) {
          fprintf(stderr, "\nEIC parameters:\n");
          fprintf(stderr, "---------------\n");
          fprintf(stderr, "disk blackbody temp.:   %6.3e\n", T_BB);
          fprintf(stderr, "tore blackbody temp.:   %6.3e\n", T_BB_tor);
          fprintf(stderr, "nuclear lumin.:         %6.3e\n", L_nuc);
          fprintf(stderr, "reprocessing fraction:  %6.3e\n", tau);
          fprintf(stderr, "R_blr:                  %6.3e\n", R_blr);
        }
        if (CASE_JET) {
          fprintf(stderr, "\nJET parameters:\n");
          fprintf(stderr, "---------------\n");
          fprintf(stderr, "Doppler factor:          %6.3f\n", DOP_JET);
          fprintf(stderr, "V jet:                   %6.3e\n", V_JET / c);
          fprintf(stderr, "Lorentz factor:          %6.3e\n", LOR_JET);
          fprintf(stderr, "V apparent:              %6.3e\n", V_JET_APP / c);
          fprintf(stderr, "Doppler factor blob-jet: %6.3e\n", DOP_B_J);
          fprintf(stderr, "V blob-jet:              %6.3e\n", V_B_J / c);
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
          fprintf(stderr, "Distance SMBH-first slice (comov frame):   %e [cm] (%f [pc])\n", X_MIN, X_MIN / pc);
          fprintf(stderr, "Distance SMBH-first slice (host frame):    %e [cm] (%f [pc])\n",
                  Y_MIN / tan(PHI * M_PI / 180.0), Y_MIN / tan(PHI * M_PI / 180.0) / pc);
          fprintf(stderr, "Distance SMBH-blob (jet comov frame):      %e [cm] (%f [pc])\n", D_b_src_j, D_b_src_j / pc);
          fprintf(stderr, "nb_slices:                        %6.3d\n\n", SL_DIM);

          /*
          if( X_MIN > D_b_src_j){
            fprintf(stderr, "we must have R_jet / tan(PHI_src) < D_blob_src_j, else the blob is NOT INSIDE THE JET...\n\nProgram aborted...\n\n");
            return 2;
          }*/
        } else cout << endl;
      }

      // SAVING ELECTRONS ENERGY SPECTRUM

      sprintf(stmp, "%s/%s_es.dat", DATA_DIRECTORY, prefix);
      if (INPUT_MODE == 0) {

        if (ifexist(stmp)) {
          sprintf(comma, "mv %s/%s_es.dat %s/%s_prev_es.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
          if (system(comma)) return 1;
        }
      }
      stream_dat = fopen(stmp, "w+");
      tmp_min = log10(GAMMA_MIN);
      tmp_max = log10(GAMMA_MAX);
      tmp_stp = (tmp_max - tmp_min) / (G_DIM - 1);
      tmp_val = tmp_min;
      for (i = 1; i <= G_DIM; i++) {
        tmp_cur = pow(10.0, tmp_val);
        dtmp = N_e(tmp_cur);
        Gamma[i] = tmp_cur;
        elec_spec[i] = dtmp * tmp_cur;
        if (dtmp > 1e-300) {
          fprintf(stream_dat, "%f %f %e %e\n", log10(tmp_cur), log10(dtmp), tmp_cur, dtmp);
        }
        tmp_val = tmp_val + tmp_stp;
      }


      fclose(stream_dat);

      elec_spec_tot_b = int_spec_b(2);
      Ub_e = m_e * c * c * elec_spec_tot_b;
      Lb_e = M_PI * R_src * R_src * LOR_B * LOR_B * c * Ub_e;

      // 1 cold proton per electron
      Lb_p = M_PI * R_src * R_src * LOR_B * LOR_B * c * m_p * c * c * int_spec_b(1);



      // PREPARING FREQUENCY VECTOR

      tmp_min = log10(FreqTransO2S(NU_STR, DOP_B, z));
      tmp_max = log10(FreqTransO2S(NU_END, DOP_B, z));
      tmp_stp = (tmp_max - tmp_min) / (NU_DIM);
      tmp_val = tmp_min;

      for (i = 1; i <= NU_DIM; i++) {
        tmp_cur = pow(10.0, tmp_val);
        NU[i] = tmp_cur;
        tmp_val = tmp_val + tmp_stp;
      }
      //Warning:: NU[1] < FreqTransO2S(NU_STR, DOP_B, z) !!
      NU[1] = FreqTransO2S(NU_STR, DOP_B, z);


      fprintf(stderr, "CALCULATING SYNCHROTRON SPECTRUM ... ");


      // CALCULATING SYNCHROTRON SPECTRUM
      if (CASE_JET == 0) {
        sprintf(stmp, "%s/%s_ss.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ss.dat %s/%s_prev_ss.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");
      }


      for (i = 1; i <= NU_DIM; i++) {
        if (PRINT) fprintf(stderr, "%3i", i);

        // emission & absorption coefficients
        jj = j_syn(N_e, GAMMA_MIN, GAMMA_MAX, NU[i], B, SYN_PREC1, SYN_PREC2);
        kk = k_esa(N_e, GAMMA_MIN, GAMMA_MAX, NU[i], B, ABS_PREC1, ABS_PREC2);
        // Here the computation is made in the source frame



        // transfer equation//
        if (L_src == 0.0) {
          I_syn = SphTransfEquat(jj, kk, R_src);
          // K. Katarzynski:
          // THIS IS TRICK !!!
          // solution of the transfer for central point of blob is
          // the same like for cylindrical geometry
          //I_rad[i] = 0.75 * CylTransfEquat(0.0, jj, kk, R_src);
          //OH, the coefficient 0.70 seems more accurate
          I_rad[i] = 0.70 * CylTransfEquat(0.0, jj, kk, R_src);
          I_rad_syn[i] = I_syn;
        } else {
          I_syn = CylTransfEquat(0.0, jj, kk, L_src);
          I_rad[i] = I_syn;
        }



        // transformation to observer frame
        nu_tmp = FreqTransS2O(NU[i], DOP_B, z);
        fx_tmp = Intens2Flux(I_syn, R_src, DOP_B, z, H_0);

        if (CASE_JET == 0) {
          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp));
          }
        }

        if (PRINT) fprintf(stderr, "\b\b\b");
      }

      if (CASE_JET == 0) fclose(stream_dat);

      fprintf(stderr, "DONE\n\n");

      //synch energy density
      Ub_syn = 4 * M_PI / c * Simpson(I_rad, NU, NU_DIM, 1, NU_DIM);


      fprintf(stderr, "CALCULATING SSC SPECTRUM ... ");

      // CALCULATING INVERSE-COMPTON SPECTRUM
      if (CASE_JET == 0 && CASE_EIC == 0) {
        sprintf(stmp, "%s/%s_cs.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_cs.dat %s/%s_prev_cs.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        sprintf(stmp, "%s/%s_tau.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_tau.dat %s/%s_prev_tau.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_tau = fopen(stmp, "w+");
      }

      for (i = 1; i <= NU_DIM; i++) {
        if (PRINT) fprintf(stderr, "%3i", i);

        // emission & absorption coefficients
        jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z),
                   NU_DIM, NU[i], COM_PREC1, COM_PREC2);

        kk = gg_abs(NU[i], I_rad, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z), SYN_PREC1,
                    SYN_PREC2);

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


        if (CASE_JET == 0 && CASE_EIC == 0) {
          // frequency transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i], DOP_B, z);

          // optical depth for absorption of VHE gamma rays by IIR
          //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
          if (IIR_level == 1) {

            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }

          } else tt = 0.;

          fprintf(stream_tau, "%f %f\n", nu_tmp * h * 0.62415,  // Conversion Hz -> TeV
                  exp(-tt));


          // flux transformation to observer frame
          fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);
          dtmp = fx_tmp;

          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);
          //fprintf(stderr, "log10 freq=%f Tau SSC EBL = %f \n",log10(nu_tmp),tt);

          //fprintf(stderr, "log10 freq=%f Tau SSC = %f \n",log10(nu_tmp),2.0 * R_src * kk);


          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }

        if (PRINT) fprintf(stderr, "\b\b\b");
      }

      if (CASE_JET == 0 && CASE_EIC == 0) {
        fclose(stream_dat);
        fclose(stream_tau);
      }

      fprintf(stderr, "DONE\n\n");

      Ub_ssc = 4 * M_PI / c * Simpson(I_rad_com, NU, NU_DIM, 1, NU_DIM);


      fprintf(stderr, "CALCULATING 2nd ORDER SSC SPECTRUM ... ");


      // CALCULATING 2nd ORDER INVERSE-COMPTON SPECTRUM

      if (CASE_JET == 0) {
        sprintf(stmp, "%s/%s_cs2.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_cs2.dat %s/%s_prev_cs2.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        sprintf(stmp, "%s/%s_tau.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_tau.dat %s/%s_prev_tau.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_tau = fopen(stmp, "w+");
      }

      sprintf(stmp1, "%s/%s_cs2_blob_frame.dat", DATA_DIRECTORY, prefix);
      stream_dat1 = fopen(stmp1, "w+");


      for (i = 1; i <= NU_DIM; i++) {
        if (PRINT) fprintf(stderr, "%3i", i);

        // emission & absorption coefficients
        jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad1st, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z),
                   NU_DIM, NU[i], COM_PREC1, COM_PREC2);

        kk = gg_abs(NU[i], I_rad1st, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z), COM_PREC1,
                    COM_PREC2);
        //computed in the source frame


        // absorption by pair production
        if (L_src == 0.0) {
          I_com2 = SphTransfEquat(jj, kk, R_src);
        } else {
          I_com2 = CylTransfEquat(0.0, jj, kk, L_src);
        }
        I_rad2nd[i] = I_com2;


        if (CASE_JET == 0) {
          // frequency transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i], DOP_B, z);

          if (IIR_level == 1) {

            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }

          } else tt = 0.;

          fprintf(stream_tau, "%f %f\n", nu_tmp * h * 0.62415,  // Conversion Hz -> TeV
                  exp(-tt));


          // flux transformation to observer frame
          fx_tmp = Intens2Flux(I_com2, R_src, DOP_B, z, H_0);
          dtmp = fx_tmp;

          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }

        nu_tmp1 = FreqTransS2O(NU[i], 1, z);
        fx_tmp1 = Intens2Flux(I_com2, R_src, 1, z, H_0);
        if (fx_tmp1 > 1.0e-300) {
          fprintf(stream_dat1, "%f %f %f\n", log10(nu_tmp1), log10(fx_tmp1), log10(nu_tmp1 * fx_tmp1));
        }

        if (PRINT) fprintf(stderr, "\b\b\b");
      }

      if (CASE_JET == 0) fclose(stream_dat);
      fclose(stream_dat1);

      fprintf(stderr, "DONE\n\n");

      Ub_ssc2 = 4 * M_PI / c * Simpson(I_rad2nd, NU, NU_DIM, 1, NU_DIM);





      // NUCLEAR LUMINOSITY

      sprintf(stmp, "%s/%s_nuc.dat", DATA_DIRECTORY, prefix);
      if (INPUT_MODE == 0) {
        if (ifexist(stmp)) {
          sprintf(comma, "mv %s/%s_nuc.dat %s/%s_prev_nuc.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
          if (system(comma)) return 1;
        }
      }

      sprintf(stmp1, "%s/%s_tor.dat", DATA_DIRECTORY, prefix);
      if (INPUT_MODE == 0) {
        if (ifexist(stmp1)) {
          sprintf(comma, "mv %s/%s_tor.dat %s/%s_prev_tor.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
          if (system(comma)) return 1;
        }
      }
      stream_dat1 = fopen(stmp1, "w+");


      if (CASE_EIC) {
        fprintf(stderr, "CALCULATING NUCLEAR LUMINOSITY SPECTRUM ... ");
        stream_dat = fopen(stmp, "w+");


        double L_BB_disk;
        double T_BB_nuc = T_BB;
        double L_BB_tor;
        double L_x = 0;
        double L_BB_disk1;

        for (i = 1; i <= NU_DIM; i++) {
          if (PRINT) fprintf(stderr, "%3i", i);

          L_BB_disk1 = L_BB_disk;

          L_BB_disk = L_nuc * Planck(NU[i], T_BB_nuc) / (sig / M_PI * pow(T_BB_nuc, 4.0));

          // X corona (Ghisellini 2009)
          if (CASE_X && L_BB_disk <= L_BB_disk1) {
            //OJ 287 test
            //L_x = 0.2 * L_nuc * pow(NU[i],-1.0) * exp(-NU[i]/3.628e+19);
            //original
            L_x = 0.3 * L_nuc * pow(NU[i], -1) * exp(-NU[i] / 3.628e+19);
          }

          L_BB_tor = L_tor * Planck(NU[i], T_BB_tor) / (sig / M_PI * pow(T_BB_tor, 4.0));

          // transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i], 1.0, z);
          if (L_BB_disk >= L_x) {
            L_BB_nuc[i] = L_BB_disk;
          } else {
            L_BB_nuc[i] = L_x;
          }
          //fprintf(stderr, "%6.3e\n", L_BB_nuc[i]);

          fx_tmp = (1.0 + z) * L_BB_nuc[i] / (4.0 * M_PI * D_L * D_L);
          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp));
          }

          nu_tmp1 = FreqTransS2O(NU[i], 1.0, z);
          fx_tmp1 = (1.0 + z) * L_BB_tor / (4.0 * M_PI * D_L * D_L);

          if (fx_tmp1 > 1.0e-300) {
            fprintf(stream_dat1, "%f %f %f\n", log10(nu_tmp1), log10(fx_tmp1), log10(nu_tmp1 * fx_tmp1));
          }
          if (PRINT) fprintf(stderr, "\b\b\b");
        }

        fclose(stream_dat1);

        fprintf(stderr, "DONE\n\n");
        //}





        sprintf(stmp, "%s/%s_ecdisc.dat", DATA_DIRECTORY, prefix);
        sprintf(stmp1, "%s/%s_prev_ecdisc.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ecdisc.dat %s/%s_prev_ecdisc.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          } else if (ifexist(stmp1)) {
            sprintf(comma, "rm %s/%s_prev_ecdisc.dat", DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }

        //if (CASE_EIC){

        fprintf(stderr, "CALCULATING EXT. INV. COMPTON SPECTRUM ON NUCLEAR RADIATION ... ");

        stream_dat = fopen(stmp, "w+");


        for (i = 1; i <= NU_DIM; i++) {
          if (PRINT) fprintf(stderr, "%3i", i);

          // radiation field from the disk in the blob frame
          I_rad_ext[i] = L_BB_nuc[i] / (4.0 * M_PI * D_b * D_b * 16.0 * pow(LOR_B, 4.0));  //Rybicky & lightman p.142
          //fprintf(stderr, "%6.3e\n", I_rad_ext[i]);

          // emission & absorption coefficients

          //How the blob sees the direct disk radiation (redshift)
          jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad_ext, FreqTransO2S(NU_STR, DOP_B, z) / LOR_B,
                     FreqTransO2S(NU_END, DOP_B, z) / LOR_B, NU_DIM, NU[i] / LOR_B, COM_PREC1, COM_PREC2);

          // absorption by pair production on direct disc, not sure if it should be taken into account
          //pure anisotropic backward radiation, already abs on BLR
          /*
          kk = gg_abs(NU[i]/LOR_B, I_rad_ext, NU_DIM,
                          FreqTransO2S(NU_STR, DOP_B, z)/LOR_B,
                          FreqTransO2S(NU_END, DOP_B, z)/LOR_B, COM_PREC1, COM_PREC2);
          */

          // absorption by pair production (on synch field)
          kk1 = gg_abs(NU[i], I_rad, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z), SYN_PREC1,
                       SYN_PREC2);

          if (L_src == 0.0) {
            I_com = SphTransfEquat(jj, kk1, R_src);
          } else {
            I_com = CylTransfEquat(0.0, jj, kk1, L_src);
          }
          I_com_disc[i] = I_com;

          // frequency transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i] / LOR_B, DOP_B, z);
          NU_IC_disk[i] = nu_tmp;


          // optical depth for absorption of VHE gamma rays by IIR
          if (IIR_level == 1) {

            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }
          } else tt = 0.;


          // flux transformation to observer frame
          fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);
          //fx_tmp = 0.;
          //fprintf(stderr, "%6.3e\n", fx_tmp);
          dtmp = fx_tmp;
          F_IC_disk[i] = dtmp;


          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }

          if (PRINT) fprintf(stderr, "\b\b\b");

        }

        fclose(stream_dat);

        fprintf(stderr, "DONE\n\n");

        if (CASE_JET == 0) {
          sprintf(stmp, "%s/%s_ecs_disk.dat", DATA_DIRECTORY, prefix);
          sprintf(stmp1, "%s/%s_prev_ecs_disk.dat", DATA_DIRECTORY, prefix);
          if (INPUT_MODE == 0) {
            if (ifexist(stmp)) {
              sprintf(comma, "mv %s/%s_ecs_disk.dat %s/%s_prev_ecs_disk.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY,
                      prefix);
              if (system(comma)) return 1;
            } else if (ifexist(stmp1)) {
              sprintf(comma, "rm %s/%s_prev_ecs_disk.dat", DATA_DIRECTORY, prefix);
              if (system(comma)) return 1;
            }
          }
        }


        fprintf(stderr, "CALCULATING EXT. INV. COMPTON SPECTRUM ON BLR ... ");

        //double R_blr = 0.1*sqrt(L_nuc/1.0e46)*3.08568e18; //Sikora et al. 2009
        //fprintf(stderr, "%6.3e\n", R_blr);

        if (CASE_JET == 0) stream_dat = fopen(stmp, "w+");

        if (D_b >= 100. * R_blr) {
          fprintf(stderr, "Blob outside 100*R_BLR!\n");
        }

        for (i = 1; i <= NU_DIM; i++) {

          if (PRINT) fprintf(stderr, "%3i", i);

          //disk radiation and BLR density both decrease as 1/r2
          //I_rad_ext[i] = L_BB_nuc[i] * LOR_B / (4.0 * M_PI * D_b*D_b) * tau * R_blr*R_blr/(D_b*D_b);

          //BLR intensity in blob frame. Use the same as Nalewajko 2014:
          I_rad_ext[i] =
              L_BB_nuc[i] * LOR_B / (4.0 * M_PI * D_b * D_b) * tau * pow(D_b / R_blr, 2) / (1. + pow(D_b / R_blr, 4));



          //Case for a finite extension of the BLR (integration up to 10 RBLR)
          /*
          if (D_b >= 10.*R_blr){
              I_rad_ext[i] = 0.;
          }
          */
          //fprintf(stderr, "%6.3e\n", L_BB_nuc[i]);


          //need to integrate I_rad_ext[i] from the blob to inf calculate the BLR absorption
          //I_rad_ext_Int[i] = I_rad_ext[i] / 3. * D_b;


          //analytical integration for Rmax_BLR = 10*R_BLR
          //I_rad_ext_Int[i] = I_rad_ext[i] * D_b/3. * (1. - pow(D_b/(10.*R_blr),3) );

          //numerical integration (useful for more complex BLR profiles, up to 100*R_BLR)
          tmp_cur = D_b;
          tmp_stp = (100. * R_blr - D_b) / 200.;
          for (j = 1; j <= 200; j++) {
            D[j] = tmp_cur;
            //I_rad_ext_D[j] = L_BB_nuc[i] * LOR_B / (4.0 * M_PI * tmp_cur*tmp_cur) * tau *R_blr*R_blr/(tmp_cur*tmp_cur); //simple 1/r2 profile
            I_rad_ext_D[j] = L_BB_nuc[i] * LOR_B / (4.0 * M_PI * tmp_cur * tmp_cur) * tau * pow(tmp_cur / R_blr, 2) /
                             (1. + pow(tmp_cur / R_blr, 4)); //Nalewajko 14
            tmp_cur += tmp_stp;
          }
          I_rad_ext_Int[i] = Simpson(I_rad_ext_D, D, 200, 1, 200);


          // emission & absorption coefficients

          //How the blob sees the BLR radiation
          jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad_ext, FreqTransO2S(NU_STR, DOP_B, z) * LOR_B,
                     FreqTransO2S(NU_END, DOP_B, z) * LOR_B, NU_DIM, NU[i] * LOR_B, COM_PREC1, COM_PREC2);
/*
        kk = gg_abs(NU[i]*LOR_B, I_rad_ext, NU_DIM,
                        FreqTransO2S(NU_STR, DOP_B, z)*LOR_B,
                        FreqTransO2S(NU_END, DOP_B, z)*LOR_B, COM_PREC1, COM_PREC2);
        //absorption on BLR done further
        */

          // absorption by pair production (inside the blob, on synch field)
          kk1 = gg_abs(NU[i], I_rad, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z), SYN_PREC1,
                       SYN_PREC2);

          if (L_src == 0.0) {
            I_com = SphTransfEquat(jj, kk1, R_src);
          } else {
            I_com = CylTransfEquat(0.0, jj, kk1, L_src);
          }

          //I_com_ext[i] = I_com_ext[i] + I_com;
          I_com_ext[i] = I_com;


          if (CASE_JET == 0) {
            // frequency transformation to observer frame
            nu_tmp = FreqTransS2O(NU[i] * LOR_B, DOP_B, z);
            //fprintf(stderr, "%6.3e\n", NU[i]);

            // optical depth for absorption of VHE gamma rays by IIR
            if (IIR_level == 1) {

              if (EBLFLAG == 0) {
                if (z <= 0.1) {
                  tt = tau_IRA_Kneiske(nu_tmp, z, 0);
                } else {
                  tt = tau_IRA_Kneiske(nu_tmp, z, 1);
                }
              } else if (EBLFLAG == 1) {
                tt = tau_IRA_Franceschini(nu_tmp, z);
              } else if (EBLFLAG == 2) {
                tt = tau_IRA_Finke(nu_tmp, z);
              } else if (EBLFLAG == 3) {
                tt = tau_IRA_Franceschini17(nu_tmp, z);
              }
            } else tt = 0.;

            // flux transformation to observer frame
            dtmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);

            // absorption by EBL
            F_IC_tot[i] = dtmp * exp(-tt);

            if (F_IC_tot[i] > 1.0e-300) {
              //fprintf(stderr, "%6.3e %6.3e\n", nu_tmp, nu_tmp*F_IC_tot[i]);
              fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(F_IC_tot[i]), log10(nu_tmp * F_IC_tot[i]),
                      log10(dtmp), log10(nu_tmp * dtmp));
            }
          }


          if (PRINT) fprintf(stderr, "\b\b\b");

        }
        if (CASE_JET == 0) fclose(stream_dat);
        fprintf(stderr, "DONE\n\n");

        //EIC from dust torus
        sprintf(stmp, "%s/%s_ecs1.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ecs1.dat %s/%s_prev_ecs1.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");


        for (i = 1; i <= NU_DIM; i++) {
          if (PRINT == 0) fprintf(stderr, "%3i", i);


          // Blackbody radiation field from the dust torus
          I_rad[i] = LOR_B * tau * L_tor / (4.0 * M_PI * R_blr * R_blr) * Planck((NU[i] / LOR_B), T_BB_tor) /
                     (sig / M_PI * pow(T_BB_tor, 4.0));


          // emission & absorption coefficients
          jj = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad, FreqTransO2S(NU_STR, DOP_B, z) * LOR_B,
                     FreqTransO2S(NU_END, DOP_B, z) * LOR_B, NU_DIM, NU[i] * LOR_B, COM_PREC1, COM_PREC2);

          kk = gg_abs(NU[i] * LOR_B, I_rad, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z) * LOR_B,
                      FreqTransO2S(NU_END, DOP_B, z) * LOR_B, COM_PREC1, COM_PREC2);


          // absorption by pair production
          if (L_src == 0.0) {
            I_com = SphTransfEquat(jj, kk, R_src);
          } else {
            I_com = CylTransfEquat(0.0, jj, kk, L_src);
          }
          I_com_ext1[i] = I_com;

          // frequency transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i] * LOR_B, DOP_B, z);

          // optical depth for absorption of VHE gamma rays by IIR
          // absorption by IIR by Kneiske et al
          if (IIR_level == 1) {

            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }
          } else tt = 0.;


          // flux transformation to observer frame
          fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0);
          dtmp = fx_tmp;

          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
          if (PRINT == 0) fprintf(stderr, "\b\b\b");
        }

        fclose(stream_dat);


        //Ub_eicd = 4*M_PI/c*Simpson(I_com_ext, NU, NU_DIM, 1, NU_DIM);
        //this is calculated further considering all absorptions

        //BLR radiation density in the blob frame
        Ub_blr = 4 * M_PI / c * Simpson(I_rad_ext, NU, NU_DIM, 1, NU_DIM);
      }







//JET PART


      sprintf(stmp, "%s/F_jet_syn.dat", DATA_DIRECTORY);
      sprintf(stmp1, "%s/F_jet_syn_prev.dat", DATA_DIRECTORY);
      if (INPUT_MODE == 0) {
        if (ifexist(stmp)) {
          sprintf(comma, "mv %s/F_jet_syn.dat %s/F_jet_syn_prev.dat", DATA_DIRECTORY, DATA_DIRECTORY);
          if (system(comma)) return 1;
        } else if (ifexist(stmp1)) {
          sprintf(comma, "rm %s/F_jet_syn_prev.dat", DATA_DIRECTORY);
          if (system(comma)) return 1;
        }
      }

      sprintf(stmp, "%s/F_jet_com.dat", DATA_DIRECTORY);
      sprintf(stmp1, "%s/F_jet_com_prev.dat", DATA_DIRECTORY);
      if (INPUT_MODE == 0) {
        if (ifexist(stmp)) {
          sprintf(comma, "mv %s/F_jet_com.dat %s/F_jet_com_prev.dat", DATA_DIRECTORY, DATA_DIRECTORY);
          if (system(comma)) return 1;
        } else if (ifexist(stmp1)) {
          sprintf(comma, "rm %s/F_jet_com_prev.dat", DATA_DIRECTORY);
          if (system(comma)) return 1;
        }
      }

      sprintf(stmp, "%s/%s_ecs_jet.dat", DATA_DIRECTORY, prefix);
      sprintf(stmp1, "%s/%s_prev_ecs_jet.dat", DATA_DIRECTORY, prefix);

      if (INPUT_MODE == 0) {
        if (ifexist(stmp)) {
          sprintf(comma, "mv %s/%s_ecs_jet.dat %s/%s_prev_ecs_jet.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
          if (system(comma)) return 1;
        } else if (ifexist(stmp1)) {
          sprintf(comma, "rm %s/%s_prev_ecs_jet.dat", DATA_DIRECTORY, prefix);
          if (system(comma)) return 1;
        }
      }

      if (CASE_JET) {

        fprintf(stderr, "CALCULATING JET SYNCHROTRON SPECTRUM... ");


        // CALCULATE PARAMETERS
        //ref comoving

        X_MIN = Y_MIN / tan(PHI_src * M_PI / 180.0);
        fprintf(stream_dat, "%e XMIN: \n", X_MIN);
        X_MAX = X_MIN + J_LEN_src * pc;
        A = Y_MIN / X_MIN;
        Y_MAX = A * X_MAX;

        // PREPARING FREQUENCY VECTOR

        tmp_min = log10(FreqTransO2S(NU_STR, DOP_JET, z));
        tmp_max = log10(FreqTransO2S(NU_END, DOP_JET, z));
        tmp_stp = (tmp_max - tmp_min) / (NU_DIM);
        tmp_val = tmp_min;


        for (i = 1; i <= NU_DIM; i++) {
          tmp_cur = pow(10.0, tmp_val);
          NU[i] = tmp_cur;
          tmp_val = tmp_val + tmp_stp;
        }

        NU[1] = FreqTransO2S(NU_STR, DOP_JET, z);


        // SPLITTING JET IN SLICES

        tmp_min = log10(X_MIN);
        tmp_max = log10(X_MAX);
        tmp_stp = (tmp_max - tmp_min) / (SL_DIM);
        tmp_val = tmp_min;

        Y_VAL[0] = 0.0;

        for (i = 1; i <= SL_DIM + 1; i++) {
          tmp_cur = pow(10.0, tmp_val);
          X_VAL[i] = tmp_cur;
          Y_VAL[i] = A * X_VAL[i];
          tmp_val = tmp_val + tmp_stp;
        }
        sprintf(stmp, "%s/geometry.dat", DATA_DIRECTORY);
        stream_dat = fopen(stmp, "w+");

        for (i = 1; i <= SL_DIM; i++) {
          DEL_X[i] = X_VAL[i + 1] - X_VAL[i];
          //print distance of each slice from the jet basis, and radius
          fprintf(stream_dat, "%e %e %f %f\n", X_VAL[i], Y_VAL[i], log10(X_VAL[i]), log10(Y_VAL[i]));
          B_VAL[i] = B_0 * pow(X_MIN / X_VAL[i], n_B);
          N_VAL[i] = N_0 * Y_MIN * Y_MIN / (Y_VAL[i] * Y_VAL[i]); //follow the jet gometry
          G_VAL[i] = GAMMA_MAX_0 * pow(X_MIN / X_VAL[i], n_G);
        }

        fclose(stream_dat);


        // CALCULATING JET SYNCHROTRON SPECTRUM

        // CALCULATING ELECTRONS ENERGY SPECTRUM


        Lj_B = 0;
        Lj_e = 0;
        for (j = 1; j <= SL_DIM; j++) {
          sprintf(stmp, "%s/es_jet_slice_%d.dat", DATA_DIRECTORY, j);
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/es_jet_slice_%d.dat %s/prev_es_jet_slice_%d.dat", DATA_DIRECTORY, j, DATA_DIRECTORY,
                    j);
            if (system(comma)) return 1;
          }
          stream_dat = fopen(stmp, "w+");

          tmp_min = log10(GAMMA_MIN1);
          tmp_max = log10(GAMMA_MAX_0);
          tmp_stp = (tmp_max - tmp_min) / (G_DIM - 1);
          tmp_val = tmp_min;
          SL_CUR = j;
          for (i = 1; i <= G_DIM; i++) {
            tmp_cur = pow(10.0, tmp_val);
            dtmp = N_e_Jet(tmp_cur);
            //STORING VECTORS FOR THE COMPUTATION OF ENERGY DENSITIES
            Gamma[i] = tmp_cur;
            elec_spec[i] = dtmp * tmp_cur;
            tmp_val = tmp_val + tmp_stp;

            fprintf(stream_dat, "%e %e\n", tmp_cur, dtmp * tmp_cur);
          }
          fclose(stream_dat);

          elec_spec_tot_j = Simpson(elec_spec, Gamma, G_DIM, 1, G_DIM);
          U_B = B_VAL[j] * B_VAL[j] / (8. * M_PI);
          Uj_e = m_e * c * c * elec_spec_tot_j;

          //We take powers through the first slice
          if (j == 1) {
            Uj_e = m_e * c * c * elec_spec_tot_j; //exactly similar to analytical solution m_e*c*c*int_spec_j(N_0, 2);
            Lj_e = M_PI * Y_VAL[j] * Y_VAL[j] * LOR_JET * LOR_JET * c * Uj_e;
            Lj_B = Lj_B + M_PI * Y_VAL[j] * Y_VAL[j] * LOR_JET * LOR_JET * c * U_B;
            Utot_e = Utot_e + Uj_e;
            Utot_B = Utot_B + U_B;
            //Cold proton power (parity with electrons)
            Lj_p = M_PI * Y_VAL[j] * Y_VAL[j] * LOR_JET * LOR_JET * c * m_p * c * c * int_spec_j(N_0, 1);
          }
          R_MOY = (Y_MAX - Y_MIN) / 2;
        }


        // CALCULATE COEFFICIENTS OF JET

        sprintf(stmp, "%s/coeff_jet.dat", DATA_DIRECTORY);
        stream_dat1 = fopen(stmp, "w+");

        //fprintf(stderr, "Jet slices coefficients (comoving frame)\n");
        for (j = 1; j <= SL_DIM; j++) {
          if (PRINT) fprintf(stderr, "%3i", j);

          SL_CUR = j;
          //fprintf(stderr, "slice: %3d X: %e dX: %e N_0: %e N: %e  B: %e  gamma_max: %e\n", j, X_VAL[j], DEL_X[j], N_VAL[j], int_spec_j(N_VAL[j]), B_VAL[j], G_VAL[j]);
          //fprintf(stream_dat1, "%d\t%e\t%e\t%e\t%e\t%e\t%e\n", j, X_VAL[j], DEL_X[j], N_VAL[j], int_spec_j(N_VAL[j]), B_VAL[j], G_VAL[j]);

          //Stot: Projected area of a slice, following the angle THETA_src_j
          Stot[j] = M_PI * Y_VAL[j] * Y_VAL[j] * cos(THETA_src_j * M_PI / 180.0) +
                    Y_VAL[j] * DEL_X[j] * sin(THETA_src_j * M_PI / 180.0);


          sprintf(stmp, "%s/f_syn_slice_%d.dat", DATA_DIRECTORY, j);
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
            I_SYN_JET[j][i] =
                ((fabs((cos(THETA_src_j * M_PI / 180.0))) * M_PI * (Y_VAL[j] * Y_VAL[j]) * I_SYN_JET_BASE[j][i]) +
                 (sin(THETA_src_j * M_PI / 180.0) * 2 * Y_VAL[j] * DEL_X[j] * I_SYN_JET_EDGE[j][i])) / Stot[j];


            // frequency transformation to observer frame
            nu_tmp = FreqTransS2O(NU[i], DOP_JET, z);

            //misaligned jet
            fx_tmp = CylIntens2Flux(I_SYN_JET_BASE[j][i], I_SYN_JET_EDGE[j][i], Y_VAL[j], DEL_X[j], DOP_JET, z, H_0,
                                    THETA_src_j);


            fprintf(stream_dat, "%e %e %e\n", nu_tmp, fx_tmp, nu_tmp * fx_tmp);
          }
          if (PRINT) fprintf(stderr, "\b\b\b");
          fclose(stream_dat);
        }
        fclose(stream_dat1);


        // CALCULATE TOTAL SYNCHROTRON SPECTRUM

        //new method 03_2016

        for (k = 1; k <= SL_DIM; k++) {
          //time gap before the central emission of a slice escape the jet by the edge
          DEL_Tph = Y_VAL[k] / (c * sin(THETA_src_j * M_PI / 180.0) - V_exp);

          for (i = 1; i <= NU_DIM; i++) {
            I_SYN = 0.0;

            for (j = k; j <= SL_DIM; j++) {
              sprintf(stmp, "%e", J_SYN_JET[j][i]);

              if ((strcmp(stmp, "-inf") == 0) || (strcmp(stmp, "inf") == 0) || (strcmp(stmp, "nan") == 0) ||
                  (J_SYN_JET[j][i] < 1.0e-300)) {
                I_SYN_JET[j][i] = 0.0;

              } else {
                //evolution time of a slice
                DEL_Tj = DEL_X[j] / V_JET;

                //photons are crossing the slice before its evolution
                if (DEL_Tj >= DEL_Tph) {
                  dtmp = DEL_Tph * c * K_ESA_JET[j][i];
                }
                  // The slice is evolving during the photons passage
                else {
                  dtmp = DEL_Tj * c * K_ESA_JET[j][i];
                  DEL_Tph -= DEL_Tj;
                  if (DEL_Tph < 0.0) DEL_Tph = 0.0;
                }

                // thick
                if (dtmp > 7.00e+2) I_SYN = I_SYN * 0.0 + I_SYN_JET[j][i];
                // transparent
                if (dtmp < 1.0e-10) I_SYN = I_SYN * 1.0 + I_SYN_JET[j][i];
                // thin
                if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) I_SYN = I_SYN * exp(-dtmp) + I_SYN_JET[j][i];
              }
            }

            //intensity of a slice k at the end of its radiation transfert through the other slices downstream
            I_SYN_TOT[k][i] = I_SYN;
            // transformation to observer frame (take the non-previously covered slice corona)
            F_SYN_JET[k][i] = M_PI * ((pow(Y_VAL[k], 2.0) - pow(Y_VAL[k - 1], 2.0)) / (D_L * D_L)) * (1.0 + z) * I_SYN;
          }
        }

        sprintf(stmp, "%s/F_jet_syn.dat", DATA_DIRECTORY);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/F_jet_syn.dat %s/F_jet_syn_prev.dat", DATA_DIRECTORY, DATA_DIRECTORY);
            if (system(comma)) return 1;
          }
        }
        errno = 0;
        stream_dat = fopen(stmp, "w+");

        sprintf(stmp1, "%s/F_jet_frame_syn.dat", DATA_DIRECTORY);
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
            dtmp3 = log10(nu_tmp * fx_tmp);

            sprintf(stmp1, "%f", dtmp2);
            sprintf(stmp2, "%f", dtmp3);

            if ((strcmp(stmp1, "-inf") == 0) || (strcmp(stmp1, "inf") == 0) || (strcmp(stmp1, "nan") == 0) ||
                (strcmp(stmp2, "-inf") == 0) || (strcmp(stmp2, "inf") == 0) || (strcmp(stmp2, "nan") == 0)) {
            } else {
              fprintf(stream_dat, "%e %e %e %f %f %f\n", nu_tmp, fx_tmp, nu_tmp * fx_tmp, dtmp1, dtmp2, dtmp3);
            }

            nu_tmp1 = FreqTransS2O(NU[i], 1, z);
            fx_tmp1 = dtmp * (1.0 + z);
            if (fx_tmp1 > 1.0e-300) {
              fprintf(stream_dat1, "%f %f %f\n", log10(nu_tmp1), log10(fx_tmp1), log10(nu_tmp1 * fx_tmp1));
            }

          }
          fclose(stream_dat);
          fclose(stream_dat1);
        } else {
          fprintf(stderr, "ERROR (%d) %s\n", errno, strerror(errno));
          fprintf(stderr, "TOTAL SYNCHROTRON SPECTRUM: '%s' NOT SAVED !!!\n", name);
          return 2;
        }

        Pj_syn = 0.0;
        for (j = 1; j <= SL_DIM; j++) {
          UJ_SLICE[j] = 4 * M_PI / c * Simpson(I_SYN_TOT[j], NU, NU_DIM, 1, NU_DIM);
          PJ_SLICE[j] = M_PI * (pow(Y_VAL[j], 2.0) - pow(Y_VAL[j - 1], 2.0)) * LOR_JET * LOR_JET * c * UJ_SLICE[j];
          Pj_syn = Pj_syn + PJ_SLICE[j];
        }

        fprintf(stderr, " DONE\n\n");


        fprintf(stderr, "CALCULATING JET SSC SPECTRUM ... ");

        // CALCULATING JET INVERSE-COMPTON SPECTRUM


        // CALCULATE COEFFICIENTS OF JET

        for (j = 1; j <= SL_DIM; j++) {
          if (PRINT) fprintf(stderr, "%3i", j);
          SL_CUR = j;

          sprintf(stmp, "%s/f_ssc_slice_%d.dat", DATA_DIRECTORY, j);
          stream_dat = fopen(stmp, "w+");

          for (i = 1; i <= NU_DIM; i++) {

            // emission & absorption coefficients
            J_COM_JET[j][i] = j_com(N_e_Jet, GAMMA_MIN1, G_VAL[j], I_SYN_JET[j], FreqTransO2S(NU_STR, DOP_JET, z),
                                    FreqTransO2S(NU_END, DOP_JET, z), NU_DIM, NU[i], COM_PREC1, COM_PREC2);

            K_ABS_SSC_JET[j][i] = gg_abs(NU[i], I_SYN_JET[j], NU_DIM, FreqTransO2S(NU_STR, DOP_JET, z),
                                         FreqTransO2S(NU_END, DOP_JET, z), SYN_PREC1, SYN_PREC2);


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
            if (IIR_level == 1) {

              if (EBLFLAG == 0) {
                if (z <= 0.1) {
                  tt = tau_IRA_Kneiske(nu_tmp, z, 0);
                } else {
                  tt = tau_IRA_Kneiske(nu_tmp, z, 1);
                }
              } else if (EBLFLAG == 1) {
                tt = tau_IRA_Franceschini(nu_tmp, z);
              } else if (EBLFLAG == 2) {
                tt = tau_IRA_Finke(nu_tmp, z);
              } else if (EBLFLAG == 3) {
                tt = tau_IRA_Franceschini17(nu_tmp, z);
              }
            } else tt = 0.;

            // flux transformation to observer frame

            I_COM_JET[j][i] =
                ((fabs((cos(THETA_src_j * M_PI / 180.0))) * M_PI * (Y_VAL[j] * Y_VAL[j]) * I_COM_JET_BASE[j][i]) +
                 (sin(THETA_src_j * M_PI / 180.0) * 2 * Y_VAL[j] * DEL_X[j] * I_COM_JET_EDGE[j][i])) / Stot[j];

            //misaligned jet
            fx_tmp = CylIntens2Flux(I_COM_JET_BASE[j][i], I_COM_JET_EDGE[j][i], Y_VAL[j], DEL_X[j], DOP_JET, z, H_0,
                                    THETA_src_j);

            // absorption by IIR
            fx_tmp = fx_tmp * exp(-tt);

            fprintf(stream_dat, "%e %e %e\n", nu_tmp, fx_tmp, nu_tmp * fx_tmp);
          }
          fclose(stream_dat);

          if (PRINT) fprintf(stderr, "\b\b\b");
        }

        // CALCULATE TOTAL COMPTON SPECTRUM


        for (k = 1; k <= SL_DIM; k++) {
          //time gap before the central emission of a slice escape the jet by the edge
          DEL_Tph = Y_VAL[k] / (c * sin(THETA_src_j * M_PI / 180.0) - V_exp);

          for (i = 1; i <= NU_DIM; i++) {
            I_COM = 0.0;

            for (j = k; j <= SL_DIM; j++) {
              sprintf(stmp, "%e", J_COM_JET[j][i]);

              if ((strcmp(stmp, "-inf") == 0) || (strcmp(stmp, "inf") == 0) || (strcmp(stmp, "nan") == 0) ||
                  (J_COM_JET[j][i] < 1.0e-300)) {
                J_COM_JET[j][i] = 0.0;

              } else {
                //evolution time of a slice
                DEL_Tj = DEL_X[j] / V_JET;

                //photons are crossing the slice before its evolution
                if (DEL_Tj >= DEL_Tph) {
                  dtmp = DEL_Tph * c * K_ABS_SSC_JET[j][i];
                }
                  // The slice is evolving during the photons passage
                else {
                  dtmp = DEL_Tj * c * K_ABS_SSC_JET[j][i];
                  DEL_Tph -= DEL_Tj;
                  if (DEL_Tph < 0.0) DEL_Tph = 0.0;
                }

                if (dtmp > 7.00e+2) I_COM = I_COM * 0.0 + I_COM_JET[j][i];
                if (dtmp < 1.0e-10) I_COM = I_COM * 1.0 + I_COM_JET[j][i];
                if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) I_COM = I_COM * exp(-dtmp) + I_COM_JET[j][i];
              }
            }

            I_COM_TOT[k][i] = I_COM;
            // transformation to observer frame
            F_COM_JET[k][i] = M_PI * ((pow(Y_VAL[k], 2.0) - pow(Y_VAL[k - 1], 2.0)) / (D_L * D_L)) * (1.0 + z) * I_COM;
          }
        }

        sprintf(stmp, "%s/F_jet_com.dat", DATA_DIRECTORY);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/F_jet_com.dat %s/F_jet_com_prev.dat", DATA_DIRECTORY, DATA_DIRECTORY);
            if (system(comma)) return 1;
          }
        }
        errno = 0;
        stream_dat = fopen(stmp, "w+");

        sprintf(stmp1, "%s/F_jet_frame_com.dat", DATA_DIRECTORY);
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
            dtmp3 = log10(nu_tmp * fx_tmp);

            sprintf(stmp1, "%f", dtmp2);
            sprintf(stmp2, "%f", dtmp3);

            if ((strcmp(stmp1, "-inf") == 0) || (strcmp(stmp1, "inf") == 0) || (strcmp(stmp1, "nan") == 0) ||
                (strcmp(stmp2, "-inf") == 0) || (strcmp(stmp2, "inf") == 0) || (strcmp(stmp2, "nan") == 0)) {
            } else {
              fprintf(stream_dat, "%e %e %e %f %f %f\n", nu_tmp, fx_tmp, nu_tmp * fx_tmp, dtmp1, dtmp2, dtmp3);
            }

            nu_tmp1 = FreqTransS2O(NU[i], 1, z);
            fx_tmp1 = dtmp * (1.0 + z);
            if (fx_tmp1 > 1.0e-300) {
              fprintf(stream_dat1, "%f %f %f\n", log10(nu_tmp1), log10(fx_tmp1), log10(nu_tmp1 * fx_tmp1));
            }

          }
          fclose(stream_dat);
          fclose(stream_dat1);

        } else {
          fprintf(stderr, "ERROR (%d) %s\n", errno, strerror(errno));
          fprintf(stderr, "TOTAL SSC SPECTRUM: '%s' NOT SAVED !!!\n", name);
          return 2;
        }

        Pj_ssc = 0.0;
        for (j = 1; j <= SL_DIM; j++) {
          UJ_SLICE[j] = 4 * M_PI / c * Simpson(I_COM_TOT[j], NU, NU_DIM, 1, NU_DIM);
          PJ_SLICE[j] = M_PI * (pow(Y_VAL[j], 2.0) - pow(Y_VAL[j - 1], 2.0)) * LOR_JET * LOR_JET * c * UJ_SLICE[j];
          Pj_ssc = Pj_ssc + PJ_SLICE[j];
        }

        fprintf(stderr, " DONE\n\n");


        fprintf(stderr, "CALCULATING EXT. INV. COMPTON SPECTRUM ON JET SYNCHROTRON... ");

        // CALCULATING EXT. INV. COMPTON SPECTRUM ON JET SYNCHROTRON



        k = 1;
        while (X_VAL[k + 1] <= D_b_src_j) {
          k += 1;
        }
        //if the blob is inside the jet
        if (X_MIN < D_b_src_j) {
          //radiation transfert along the jet until the blob position, independent of the angle THETA
          //parallel plans
          for (l = 1; l <= k; l++) {
            //fprintf(stderr, "%e\n\n",l);
            for (i = 1; i <= NU_DIM; i++) {
              I_SYN = 0.0;
              for (j = l; j <= k; j++) {
                sprintf(stmp, "%e", J_SYN_JET[j][i]);
                //fprintf(stderr, "%d\n",j);
                if ((strcmp(stmp, "-inf") == 0) || (strcmp(stmp, "inf") == 0) || (strcmp(stmp, "nan") == 0) ||
                    (J_SYN_JET[j][i] < 1.0e-300)) {
                  I_SYN_JET[j][i] = 0.0;
                } else {
                  if (j < k) {
                    dtmp = DEL_X[j] * K_ESA_JET[j][i];
                    // transparent
                    if (dtmp < 1.0e-10) I_SYN = I_SYN * 1.0 + J_SYN_JET[j][i] * DEL_X[j];
                  } else {
                    dtmp = (D_b_src_j - X_VAL[k]) * K_ESA_JET[j][i];
                    // transparent
                    if (dtmp < 1.0e-10) I_SYN = I_SYN * 1.0 + J_SYN_JET[j][i] * (D_b_src_j - X_VAL[k]);
                  }

                  // thick
                  if (dtmp > 7.00e+2) I_SYN = I_SYN * 0.0 + J_SYN_JET[j][i] / K_ESA_JET[j][i];
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
              if ((strcmp(stmp, "-inf") == 0) || (strcmp(stmp, "inf") == 0) || (strcmp(stmp, "nan") == 0) ||
                  (J_SYN_JET[j][i] < 1.0e-300)) {
                I_SYN_JET[j][i] = 0.0;

              } else {
                if (j > k) {
                  dtmp = DEL_X[j] * K_ESA_JET[j][i];
                  // transparent
                  if (dtmp < 1.0e-10) I_SYN = I_SYN * 1.0 + J_SYN_JET[j][i] * DEL_X[j];

                } else {
                  dtmp = (X_VAL[k + 1] - D_b) * K_ESA_JET[j][i];
                  // transparent
                  if (dtmp < 1.0e-10) I_SYN = I_SYN * 1.0 + J_SYN_JET[j][i] * (X_VAL[k + 1] - D_b);

                }
                // thick
                if (dtmp > 7.00e+2) I_SYN = I_SYN * 0.0 + J_SYN_JET[j][i] / K_ESA_JET[j][i];
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



          //if the blob is inside the jet (outdated)
          /*
          if ( X_MIN < D_b_src_j){
             // transformation jet to blob frame (blueshift)
             I_rad_ext1[i] = (I_rad_ext1[i] + I_rad_ext_s[i]) * LOR_B_J;
          } else {
              //if the blob is upstream the jet (i.e upstream radio core)
              I_rad_ext1[i] = (I_rad_ext1[i] + I_rad_ext_s[i]) * LOR_B_J / (4.0 * M_PI * (D_BJ + Y_MIN)*(D_BJ + Y_MIN) / (Y_MIN*Y_MIN));
          }
        */
          //updated version: Hervet 2021
          //if the blob is at distance > Xmin+Ymin, no change
          if (D_b_src_j >= X_MIN + Y_MIN) {
            //solid angle full sphere
            dtmp = 4.0 * M_PI;
          } else if (D_b_src_j >= X_MIN and D_b_src_j < X_MIN + Y_MIN) {
            //height of the polar cap outside the jet
            dtmp1 = Y_MIN - (D_b_src_j - X_MIN);
            //solid angle of a sphere removed by the polar cap
            dtmp = 4.0 * M_PI - 2.0 * M_PI * dtmp1 / Y_MIN;
          } else if (D_b_src_j < X_MIN) {
            //radius of the sphere centered on the blob  passing through Xmin
            dtmp2 = X_MIN - D_b_src_j;
            dtmp3 = sqrt(dtmp2 * dtmp2 + Y_MIN * Y_MIN);
            //height of the polar cap inside the jet
            dtmp1 = dtmp3 - dtmp2;
            //solid angle of a sphere defined by the polar cap
            dtmp = 2.0 * M_PI * dtmp1 / Y_MIN;
          }
          //resulting intensity in the blob frame
          I_rad_ext1[i] = (I_rad_ext1[i] + I_rad_ext_s[i]) * dtmp / (4.0 * M_PI);



          // emission & absorption coefficients (how the blob sees the jet)
          C_e1[i] = j_com(N_e, GAMMA_MIN, GAMMA_MAX, I_rad_ext1, FreqTransO2S(NU_STR, DOP_JET, z) * LOR_B_J,
                          FreqTransO2S(NU_END, DOP_JET, z) * LOR_B_J, NU_DIM, NU[i] * LOR_B_J, COM_PREC1, COM_PREC2);

          C_a1[i] = gg_abs(NU[i] * LOR_B_J, I_rad_ext1, NU_DIM, FreqTransO2S(NU_STR, DOP_JET, z) * LOR_B_J,
                           FreqTransO2S(NU_END, DOP_JET, z) * LOR_B_J, COM_PREC1, COM_PREC2);

          if (PRINT) fprintf(stderr, "\b\b\b");
        }

        sprintf(stmp1, "%s/%s_ecs_jet_frame.dat", DATA_DIRECTORY, prefix);
        stream_dat1 = fopen(stmp1, "w+");


        for (i = 1; i <= NU_DIM; i++) {
          jj = C_e1[i];
          kk = C_a1[i];


          // radiation transfert taking into account the pair absorption
          if (L_src == 0.0) {
            I_com = SphTransfEquat(jj, kk, R_src);
            I_eic_jet[i] = I_com;
          } else {
            I_com = CylTransfEquat(0.0, jj, kk, L_src);
            I_eic_jet[i] = I_com;
          }

          //transformation in the observer frame
          nu_tmp1 = FreqTransS2O(NU[i] * LOR_B_J, DOP_B, z);
          fx_tmp1 = Intens2Flux(I_com, R_src, 1, z, H_0);
          if (fx_tmp1 > 1.0e-300) {
            fprintf(stream_dat1, "%f %f %f\n", log10(nu_tmp1), log10(fx_tmp1), log10(nu_tmp1 * fx_tmp1));
          }
        }
        fclose(stream_dat1);

        fprintf(stderr, "DONE\n\n");

        Ub_eicj = 4 * M_PI / c * Simpson(I_eic_jet, NU, NU_DIM, 1, NU_DIM);
      }


      fprintf(stderr, "CALCULATING RADIATION ABSORPTION OF THE BLOB EMISSION THROUGH THE BLR ... ");

      // PREPARING FREQUENCY VECTOR
      tmp_min = log10(FreqTransO2S(NU_STR, DOP_B, z) * LOR_B);
      tmp_max = log10(FreqTransO2S(NU_END, DOP_B, z) * LOR_B);
      tmp_stp = (tmp_max - tmp_min) / (NU_DIM);
      tmp_val = tmp_min;
      for (i = 1; i <= NU_DIM; i++) {
        tmp_cur = pow(10.0, tmp_val);
        NU[i] = tmp_cur;
        tmp_val = tmp_val + tmp_stp;
      }
      //Warning:: NU[1] < FreqTransO2S(NU_STR, DOP_B, z) !!
      NU[1] = FreqTransO2S(NU_STR, DOP_B, z) * LOR_B;
      fprintf(stderr, "%6.3e\n", NU_END);

      //absorption of SSC
      //blob frame
      for (i = 1; i <= NU_DIM; i++) {
        //fprintf(stderr, "%6.3e\n", NU[i]);

        // absorption by pair production on BLR soft photons, use of the already integrated rad field: I_rad_ext_Int
        kk = gg_abs(NU[i], I_rad_ext_Int, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z),
                    COM_PREC1, COM_PREC2);


        //dtmp = (R_blr - D_b) *kk;
        dtmp = kk; //(already integrated along the BLR extension)
        //remove BLR absorbption just for tests
        //dtmp = 0.;
        //small code to check the opacity at a given energy in the observer frame
        if (dtmp > 0.) {
          //fprintf(stderr, "BLR opacity: %6.3e  Eblob obs frame: %6.3e [TeV]\n", dtmp,  FreqTransS2O(NU[i], DOP_B, z)/LOR_B / (HZ_PER_EV*1e+12));
          dtmp1 = 0.266; //TeV
          //linear interpolation
          dtmp2 = FreqTransS2O(NU[i], DOP_B, z) / LOR_B / (HZ_PER_EV * 1e+12);
          dtmp3 = FreqTransS2O(NU[i + 1], DOP_B, z) / LOR_B / (HZ_PER_EV * 1e+12);
          kk = gg_abs(NU[i + 1], I_rad_ext_Int, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z), FreqTransO2S(NU_END, DOP_B, z),
                      COM_PREC1, COM_PREC2);
          if ((dtmp2 <= dtmp1) && (dtmp3 > dtmp1)) {
            dtmp4 = (kk - dtmp) * (dtmp1 - dtmp2) / (dtmp3 - dtmp2) + dtmp;
            fprintf(stderr, "BLR opacity: %6.3e  Eblob obs frame: %6.3e [TeV]\n", dtmp4, dtmp1);
          }
        }

        // thick
        if (dtmp > 7.00e+2) {
          I_rad_com[i] = 0.0;
          I_rad2nd[i] = 0.0;
        }
        // transparent
        //nothing change
        // thin
        if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) {
          I_rad_com[i] = I_rad_com[i] * exp(-dtmp);
          I_rad2nd[i] = I_rad2nd[i] * exp(-dtmp);
        }
      }

      // PREPARING FREQUENCY VECTOR
      tmp_min = log10(FreqTransO2S(NU_STR, DOP_B, z));
      tmp_max = log10(FreqTransO2S(NU_END, DOP_B, z));
      tmp_stp = (tmp_max - tmp_min) / (NU_DIM);
      tmp_val = tmp_min;
      for (i = 1; i <= NU_DIM; i++) {
        tmp_cur = pow(10.0, tmp_val);
        NU[i] = tmp_cur;
        tmp_val = tmp_val + tmp_stp;
      }
      //Warning:: NU[1] < FreqTransO2S(NU_STR, DOP_B, z) !!
      NU[1] = FreqTransO2S(NU_STR, DOP_B, z);


      //blob frame
      //absorption of EIC
      for (i = 1; i <= NU_DIM; i++) {
        //fprintf(stderr, "%6.3e\n", NU[i]);

        // absorption by pair production on BLR soft photons, use of the already integrated rad field: I_rad_ext_Int
        kk = gg_abs(NU[i] * LOR_B, I_rad_ext_Int, NU_DIM, FreqTransO2S(NU_STR, DOP_B, z) * LOR_B,
                    FreqTransO2S(NU_END, DOP_B, z) * LOR_B, COM_PREC1, COM_PREC2);

        dtmp = kk; //(already integrated along the BLR extension)
        //remove BLR absorbption just for tests
        //dtmp = 0.;

        /*
        // absorption by pair production (on BLR soft photons)
        kk = gg_abs(NU[i]*LOR_B, I_rad_ext, NU_DIM,
                        FreqTransO2S(NU_STR, DOP_B, z)*LOR_B,
                        FreqTransO2S(NU_END, DOP_B, z)*LOR_B, COM_PREC1, COM_PREC2);


        dtmp = (R_blr - D_b) *kk;
        */


        //fprintf(stderr, "%6.3e\n", dtmp);
        // thick
        if (dtmp > 7.00e+2) {
          I_com_ext[i] = 0.0;

        }
        // transparent
        //nothing change
        // thin
        if ((dtmp < 7.00e+2) && (dtmp > 1.0e-10)) {
          I_com_ext[i] = I_com_ext[i] * exp(-dtmp);
        }
      }

      if (CASE_JET == 0) {
        //SSC
        sprintf(stmp, "%s/%s_cs.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_cs.dat %s/%s_prev_cs.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");


        for (i = 1; i <= NU_DIM; i++) {
          // frequency transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i], DOP_B, z);

          // optical depth for absorption of VHE gamma rays by IIR
          //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }

          } else tt = 0.;


          // flux transformation to observer frame
          fx_tmp = Intens2Flux(I_rad_com[i], R_src, DOP_B, z, H_0);
          dtmp = fx_tmp;

          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);

        //SSC_2nd
        sprintf(stmp, "%s/%s_cs2.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_cs2.dat %s/%s_prev_cs2.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        for (i = 1; i <= NU_DIM; i++) {
          // frequency transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i], DOP_B, z);

          // optical depth for absorption of VHE gamma rays by IIR
          //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }

          } else tt = 0.;


          // flux transformation to observer frame
          fx_tmp = Intens2Flux(I_rad2nd[i], R_src, DOP_B, z, H_0);
          dtmp = fx_tmp;

          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);

        //EIC from disk-BLR
        sprintf(stmp, "%s/%s_ecs.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ecs.dat %s/%s_prev_ecs.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        for (i = 1; i <= NU_DIM; i++) {
          // frequency transformation to observer frame
          nu_tmp = FreqTransS2O(NU[i] * LOR_B, DOP_B, z);

          // optical depth for absorption of VHE gamma rays by IIR
          //tt     = tau_IIR(nu_tmp, z, IIR_level-1); // low level of IIR
          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }

          } else tt = 0.;


          // flux transformation to observer frame
          fx_tmp = Intens2Flux(I_com_ext[i], R_src, DOP_B, z, H_0);
          dtmp = fx_tmp;

          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);
        Ub_eicd = 4 * M_PI / c * Simpson(I_com_ext, NU, NU_DIM, 1, NU_DIM);

      }


      if (CASE_JET) {

        fprintf(stderr, "CALCULATING RADIATION ABSORPTION OF THE BLOB EMISSION THROUGH THE JET ... ");

        // PREPARING FREQUENCY VECTOR

        tmp_min = log10(FreqTransO2S(NU_STR, DOP_B, z));
        tmp_max = log10(FreqTransO2S(NU_END, DOP_B, z));
        tmp_stp = (tmp_max - tmp_min) / (NU_DIM);
        tmp_val = tmp_min;

        for (i = 1; i <= NU_DIM; i++) {
          tmp_cur = pow(10.0, tmp_val);
          NU[i] = tmp_cur;
          tmp_val = tmp_val + tmp_stp;
        }

        j = 1;
        while (X_VAL[j] <= D_b_src_j) {
          j += 1;
          // j: first slice after the blob
        }

        //blob to jet frame
        for (i = 1; i <= NU_DIM; i++) {
          I_rad_syn[i] = I_rad_syn[i] * pow(DOP_B_J, 3.0);
          I_rad_com[i] = I_rad_com[i] * pow(DOP_B_J, 3.0);
          I_rad2nd[i] = I_rad2nd[i] * pow(DOP_B_J, 3.0);
          I_com_ext[i] = I_com_ext[i] * pow(DOP_B_J, 3.0);
          I_com_ext1[i] = I_com_ext1[i] * pow(DOP_B_J, 3.0);
          I_eic_jet[i] = I_eic_jet[i] * pow(DOP_B_J, 3.0);
        }

        //new method 03_2016
        //radiative transfert along the propagation of the jet

        //time gap before the blob emission escape the jet by the edge
        DEL_Tph = Y_VAL[j - 1] / (c * sin(THETA_src_j * M_PI / 180.0) - V_exp);
        k = j - 1;
        dtmp2 = 0;
        //fprintf(stderr, "%6.3e DEL_Tph", DEL_Tph);
        while (k <= SL_DIM && DEL_Tph > 0.0) {
          fprintf(stderr, "%3i", k);
          //evolution time of a slice
          if (k == j - 1) {
            DEL_Tj = (X_VAL[k + 1] - D_b_src_j) / V_JET;
          } else {
            DEL_Tj = DEL_X[k] / V_JET;
          }
          for (i = 1; i <= NU_DIM; i++) {
            //photons cross the slice before its evolution
            if (DEL_Tj >= DEL_Tph) {
              dtmp = DEL_Tph * c * K_ESA_JET[k][i];
              dtmp1 = DEL_Tph * c * K_ABS_SSC_JET[k][i];
            }
              // The slice is evolving during the photons passage
            else {
              dtmp = DEL_Tj * c * K_ESA_JET[k][i];
              dtmp1 = DEL_Tj * c * K_ABS_SSC_JET[k][i];
              DEL_Tph -= DEL_Tj;
            }
            // thick
            if (dtmp > 7.00e+2) I_rad_syn[i] = 0.0;

            if (dtmp1 > 7.00e+2) {
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

            if ((dtmp1 < 7.00e+2) && (dtmp1 > 1.0e-10)) {
              I_rad_com[i] = I_rad_com[i] * exp(-dtmp1);
              I_rad2nd[i] = I_rad2nd[i] * exp(-dtmp1);
              I_eic_jet[i] = I_eic_jet[i] * exp(-dtmp1);
              I_com_ext[i] = I_com_ext[i] * exp(-dtmp1);
              I_com_ext1[i] = I_com_ext1[i] * exp(-dtmp1);
              //used to write the jet opacity at a given energy
              dtmp3 = (FreqTransS2O(NU[i] * DOP_B_J, DOP_B / DOP_B_J, z) / (HZ_PER_EV * 1e+12));
              //fprintf(stderr, "%6.3e Energy ", dtmp3);
              if ((0.4 > dtmp3) && (dtmp3 > 0.37)) {
                dtmp2 += dtmp1;
              }
            }
            //fprintf(stderr, "Jet opacity: %6.3e  Eblob obs frame: %6.3e [TeV]\n", dtmp1,  FreqTransS2O(NU[i] * DOP_B_J, DOP_B/DOP_B_J , z) / (HZ_PER_EV*1e+12));
          }
          k += 1;
          fprintf(stderr, "\b\b\b");
        }



        // synchrotron
        sprintf(stmp, "%s/%s_ss.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ss.dat %s/%s_prev_ss.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        // transformation to observer frame
        for (i = 1; i <= NU_DIM; i++) {
          nu_tmp = FreqTransS2O(NU[i] * DOP_B_J, DOP_B / DOP_B_J, z);
          I_syn = I_rad_syn[i];
          fx_tmp = Intens2Flux(I_syn, R_src, DOP_B, z, H_0) / pow((DOP_B_J), 3.0);


          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp));
          }
        }
        fclose(stream_dat);


        //SSC
        sprintf(stmp, "%s/%s_cs.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_cs.dat %s/%s_prev_cs.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        // transformation to observer frame
        for (i = 1; i <= NU_DIM; i++) {
          nu_tmp = FreqTransS2O(NU[i] * DOP_B_J, DOP_B / DOP_B_J, z);
          I_com = I_rad_com[i];
          fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0) / pow((DOP_B_J), 3.0);


          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }
          } else tt = 0.;


          dtmp = fx_tmp;
          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);


        //SSC_2nd
        sprintf(stmp, "%s/%s_cs2.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_cs2.dat %s/%s_prev_cs2.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        // transformation to observer frame
        for (i = 1; i <= NU_DIM; i++) {
          nu_tmp = FreqTransS2O(NU[i] * DOP_B_J, DOP_B / DOP_B_J, z);
          I_com2 = I_rad2nd[i];
          fx_tmp = Intens2Flux(I_com2, R_src, DOP_B, z, H_0) / pow((DOP_B_J), 3.0);

          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }
          } else tt = 0.;


          dtmp = fx_tmp;
          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);


        //eic from jet
        sprintf(stmp, "%s/%s_ecs_jet.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ecs_jet.dat %s/%s_prev_ecs_jet.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY,
                    prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        // transformation to observer frame
        for (i = 1; i <= NU_DIM; i++) {
          nu_tmp = FreqTransS2O(NU[i] * DOP_B_J * LOR_B_J, DOP_B / DOP_B_J, z);
          I_com = I_eic_jet[i];
          fx_tmp = M_PI * ((R_src * R_src) / (D_L * D_L)) * (1.0 + z) * (pow(DOP_B, 3.0) / pow((DOP_B_J), 3.0)) * I_com;

          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }
          } else tt = 0.;


          dtmp = fx_tmp;
          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);


        //eic from disk-BLR
        sprintf(stmp, "%s/%s_ecs.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ecs.dat %s/%s_prev_ecs.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        // transformation to observer frame
        for (i = 1; i <= NU_DIM; i++) {
          nu_tmp = FreqTransS2O(NU[i] * LOR_B * DOP_B_J, DOP_B / DOP_B_J, z);
          I_com = I_com_ext[i];
          fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0) / pow((DOP_B_J), 3.0);

          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }
          } else tt = 0.;


          dtmp = fx_tmp;
          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);


        //eic from torus-BLR
        sprintf(stmp, "%s/%s_ecs1.dat", DATA_DIRECTORY, prefix);
        if (INPUT_MODE == 0) {
          if (ifexist(stmp)) {
            sprintf(comma, "mv %s/%s_ecs1.dat %s/%s_prev_ecs1.dat", DATA_DIRECTORY, prefix, DATA_DIRECTORY, prefix);
            if (system(comma)) return 1;
          }
        }
        stream_dat = fopen(stmp, "w+");

        // transformation to observer frame
        for (i = 1; i <= NU_DIM; i++) {
          nu_tmp = FreqTransS2O(NU[i] * LOR_B * DOP_B_J, DOP_B / DOP_B_J, z);
          I_com = I_com_ext1[i];
          fx_tmp = Intens2Flux(I_com, R_src, DOP_B, z, H_0) / pow((DOP_B_J), 3.0);

          if (IIR_level == 1) {
            if (EBLFLAG == 0) {
              if (z <= 0.1) {
                tt = tau_IRA_Kneiske(nu_tmp, z, 0);
              } else {
                tt = tau_IRA_Kneiske(nu_tmp, z, 1);
              }
            } else if (EBLFLAG == 1) {
              tt = tau_IRA_Franceschini(nu_tmp, z);
            } else if (EBLFLAG == 2) {
              tt = tau_IRA_Finke(nu_tmp, z);
            } else if (EBLFLAG == 3) {
              tt = tau_IRA_Franceschini17(nu_tmp, z);
            }
          } else tt = 0.;


          dtmp = fx_tmp;
          // absorption by IIR
          fx_tmp = fx_tmp * exp(-tt);

          if (fx_tmp > 1.0e-300) {
            fprintf(stream_dat, "%f %f %f %f %f\n", log10(nu_tmp), log10(fx_tmp), log10(nu_tmp * fx_tmp), log10(dtmp),
                    log10(nu_tmp * dtmp));
          }
        }
        fclose(stream_dat);

        fprintf(stderr, "Jet opacity: %6.3e  Eblob observer frame: 3.703e-01 [TeV]\n",
                dtmp2); //,  FreqTransS2O(NU[i] * DOP_B_J, DOP_B/DOP_B_J , z) / (HZ_PER_EV*1e+12));



        fprintf(stderr, " DONE\n\n");
        fprintf(stderr, "Jet slice where the blob takes place:  %6.3d\n", j - 1);
      }


      fprintf(stderr, "\n\nDERIVED PARAMETERS FROM THE BLOB:\n");
      fprintf(stderr, "----------\n");
      fprintf(stderr, "Velocity:                    %6.3e c \n", V_B / c);
      Tcool_synch = 3.0 * m_e * c / (4 * sig_T * GAMMA_MAX * Ub_B);
      Tcool_vhe = 3.0 * m_e * c / (4 * sig_T * GAMMA_MAX * Ub_syn);
      Tcool_radio = 3.0 * m_e * c / (4 * sig_T * GAMMA_MIN * Ub_syn);
      Tcool_synch_gbreak = 3.0 * m_e * c / (4 * sig_T * GAMMA_BRK * Ub_B);
      Tcool_SSC_gbreak = 3.0 * m_e * c / (4 * sig_T * GAMMA_BRK * Ub_syn);
      Tcool_EIC_gbreak = 3.0 * m_e * c / (4 * sig_T * GAMMA_BRK * Ub_blr);


      fprintf(stderr, "Observed synch Gamma_max cool. time:     %6.3e s \n", Tcool_synch * (1 + z) / DOP_B);
      fprintf(stderr, "Observed IC Gamma_max cool. time:        %6.3e s \n", Tcool_vhe * (1 + z) / DOP_B);
      fprintf(stderr, "Observed IC Gamma_min cool. time:        %6.3e s \n", Tcool_radio * (1 + z) / DOP_B);
      fprintf(stderr, "Observed synch Gamma_break cool. time:    %6.3e s \n", Tcool_synch_gbreak * (1 + z) / DOP_B);
      fprintf(stderr, "Observed SSC Gamma_break cool. time:      %6.3e s \n", Tcool_SSC_gbreak * (1 + z) / DOP_B);
      fprintf(stderr, "Observed EIC BLR Gamma_break cool. time:  %6.3e s \n", Tcool_EIC_gbreak * (1 + z) / DOP_B);
      dtmp = 3.0 * m_e * c / (4 * sig_T * GAMMA_MAX * (Ub_blr + Ub_syn + Ub_B));
      fprintf(stderr, "Observed combined Gamma_max cool. time (synch+SSC+IC blr):  %6.3e h \n",
              dtmp * (1 + z) / (DOP_B * 3600.));
      dtmp = 3.0 * m_e * c / (4 * sig_T * GAMMA_BRK * (Ub_blr + Ub_syn + Ub_B));
      fprintf(stderr, "Observed combined Gamma_break cool. time (synch+SSC+IC blr):  %6.3e h \n",
              dtmp * (1 + z) / (DOP_B * 3600.));
      //fprintf(stderr, "Intrinsic adiabatic cool. time:           %6.3e s \n", R_src/V_B); //not true, this is a lower limit
      fprintf(stderr, "Observed minimal variability:             %6.3e h \n", R_src * (1 + z) / (c * DOP_B * 3600.));
      fprintf(stderr, "Electron particle density:   %6.3e cm^-3 \n", int_spec_b(1));
      fprintf(stderr, "Electron energy density:     %6.3e erg cm^-3 \n", Ub_e);
      fprintf(stderr, "Synchrotron energy density:  %6.3e erg cm^-3 \n", Ub_syn);
      fprintf(stderr, "Magnetic energy density:     %6.3e erg cm^-3 \n", Ub_B);
      fprintf(stderr, "U_B/U_e:                     %6.3e  \n\n\n", Ub_B / Ub_e);
      Lb_r = M_PI * R_src * R_src * LOR_B * LOR_B * c * (Ub_syn + Ub_ssc + Ub_ssc2 + Ub_eicd + Ub_eicj);
      Lj_r = Pj_syn + Pj_ssc;
      fprintf(stderr, "Magnetic power:              %6.3e erg s^-1 \n", Lb_B);
      fprintf(stderr, "Electrons power:             %6.3e erg s^-1 \n", Lb_e);
      fprintf(stderr, "Cold proton power (parity):  %6.3e erg s^-1 \n", Lb_p);
      fprintf(stderr, "Syn power:                   %6.3e erg s^-1 \n",
              M_PI * R_src * R_src * LOR_B * LOR_B * c * Ub_syn);
      fprintf(stderr, "Ssc power:                   %6.3e erg s^-1 \n",
              M_PI * R_src * R_src * LOR_B * LOR_B * c * Ub_ssc);
      fprintf(stderr, "Ssc2 power:                  %6.3e erg s^-1 \n",
              M_PI * R_src * R_src * LOR_B * LOR_B * c * Ub_ssc2);
      fprintf(stderr, "Eic disk power:              %6.3e erg s^-1 \n",
              M_PI * R_src * R_src * LOR_B * LOR_B * c * Ub_eicd);
      fprintf(stderr, "Eic jet power:               %6.3e erg s^-1 \n",
              M_PI * R_src * R_src * LOR_B * LOR_B * c * Ub_eicj);
      fprintf(stderr, "Radiation power (beware, not all abs. considered):             %6.3e erg s^-1 \n", Lb_r);

      if (CASE_EIC) {
        fprintf(stderr, "Blob distance:             %6.3e X R_BLR\n", D_b / R_blr);
      }


      if (CASE_JET) {
        fprintf(stderr, "\n\nDERIVED PARAMETERS FROM THE JET:\n");
        fprintf(stderr, "----------\n\n");
        fprintf(stderr, "velocity of jet:             %e [in c]\n", V_JET / c);
        fprintf(stderr, "Lorentz factor of jet:       %e\n", LOR_JET);
        fprintf(stderr, "apparent velocity of jet     %e [in c]\n", V_JET_APP / c);
        fprintf(stderr, "Minimal variability:         %6.3e h \n", Y_MIN * (1 + z) / (c * DOP_JET * 3600.));
        fprintf(stderr, "outer radius of jet:         %e [cm] (%f [pc])\n", Y_MAX, Y_MAX / pc);
        fprintf(stderr, "theta in the jet frame:      %e [deg]\n", THETA_src_j);
        fprintf(stderr, "Magnetic energy density:     %6.3e erg cm^-3 \n", Utot_B);
        fprintf(stderr, "Electron energy density:     %6.3e erg cm^-3 \n", Utot_e);
        //fprintf(stderr, "Proton energy density:       %6.3e erg cm^-3 \n", Utot_p);
        fprintf(stderr, "U_B/U_e:                     %6.3e  \n\n\n", Utot_B / Utot_e);
        fprintf(stderr, "Magnetic power:              %6.3e erg s^-1 \n", Lj_B);
        fprintf(stderr, "Electrons power:             %6.3e erg s^-1 \n", Lj_e);
        fprintf(stderr, "Cold proton power (parity):  %6.3e erg s^-1 \n", Lj_p);
        fprintf(stderr, "Syn power:                   %6.3e erg s^-1 \n", Pj_syn);
        fprintf(stderr, "Ssc power:                   %6.3e erg s^-1 \n", Pj_ssc);
        fprintf(stderr, "Radiation power:             %6.3e erg s^-1 \n", Lj_r);
      }

      fprintf(stderr, "\n\nGENERAL DERIVED PARAMETERS:\n");
      fprintf(stderr, "----------\n");
      fprintf(stderr, "luminosity distance:                                %e [M pc]\n", D_L / (1.0e+6 * pc));
      fprintf(stderr, "Total radiation power (nucleus included):           %6.3e erg s^-1 \n", Lb_r + Lj_r + L_nuc);
      fprintf(stderr, "Total power (without cold particles):               %6.3e erg s^-1 \n",
              Lb_r + Lj_r + L_nuc + Lb_B + Lj_B + Lb_e + Lj_e);
      fprintf(stderr, "Total power (without cold particles & nucleus):     %6.3e erg s^-1 \n",
              Lb_r + Lj_r + Lb_B + Lj_B + Lb_e + Lj_e);
      fprintf(stderr, "Total power (with cold particles & without nucleus):%6.3e erg s^-1 \n",
              Lb_r + Lj_r + Lb_B + Lj_B + Lb_e + Lj_e + Lb_p + Lj_p);


      fprintf(stderr, "\nALL CALCULATIONS DONE !!! ");

      T_END = time(NULL);

      if (difftime(T_END, T_START) > 60)
        sprintf(stmp, "[time: %4.1f min]\n", difftime(T_END, T_START) / 60.0);
      else
        sprintf(stmp, "[time: %4.1f sec]\n", difftime(T_END, T_START));

      fprintf(stderr, "%s\n", stmp);

      fprintf(stderr, "SAVING FINAL SPECTRUM\n");

      return 0;
    }

}


int main(int argc, char **argv) {
  int num_models = 2; // blob, blob + nuc
  sprintf(PARA_FN, "%s", "none");

  // LOAD PARAMETERS
  if (argc == 2) { // simple run with parameter file in default location
    if (strcmp(argv[1], "--help") == 0) {
      description();
      return 0;
    }
    strcpy(PARA_FN, argv[1]);
    if (load_params(PARA_FN)); else return 2;
  } else if (argc >= 3) { // input mode is specified

    int input_type = atoi(argv[1]);

    // set input mode
    if (input_type < 0 or input_type > 3) {
      std::cerr << "Valid input modes are 0-3\n";
      return 2;
    }
    bj_mcmc02::INPUT_MODE = input_type % 2;

    if (argc == 3 && input_type < 2) {  // called with input mode + param file
      // get and read param file
      strcpy(PARA_FN, argv[2]); // param file
      if (load_params(PARA_FN)); else return 2;
    } else if (argc == 4 && input_type >= 2) { // called with input mode + data file + param file
      strcpy(DATA_DIRECTORY, argv[2]);
      strcpy(PARA_FN, argv[3]);
      if (load_params(PARA_FN)); else return 2;
    } else if (argc >= 22) { // called with command line args
      int index = 2;
      // set data folder if necessary
      if (input_type >= 2) {
        strcpy(DATA_DIRECTORY, argv[index]);
        index++;
      }

      // model type
      int model_type = atoi(argv[index]);
      index++;
      // check if a model type number
      if (model_type < 0 || model_type >= num_models) {
        std::cerr << "Incorrect model specification number " << model_type << ". Valid values are < " << num_models
                  << endl;
        return 2;
      }
      // check if right number of args for model type
      if (model_type == 0 && ((input_type < 2 && argc != 22) || (input_type >= 2 & argc != 23))) {
        std::cerr << "Incorrect number of parameters for model.\n";
        return 2;
      }
      if (model_type == 1 && ((input_type < 2 && argc != 28) || (input_type >= 2 & argc != 29))) {
        std::cerr << "Incorrect number of parameters for model.\n";
        return 2;
      }

      // everything is correct, do the model
      if (!load_params_from_list(argv, model_type, index, argc))
        return 2;
    }
  } else {
    cout << "Valid usages:\n"
            "- Usage 1: bj_mcmc --help\n"
            "- Usage 2: bj_mcmc <parameter file>\n"
            "- Usage 3: bj_mcmc 0 <parameter file>                             same execution as usage 2, _prev files made\n"
            "- Usage 4: bj_mcmc 1 <parameter file>                             no _prev files made\n"
            "- Usage 5: bj_mcmc 2 <data folder> <parameter file>               _prev files made\n"
            "- Usage 6: bj_mcmc 3 <data folder> <parameter file>               no _prev files made\n"
            "- Usage 7: bj_mcmc 0/1 <model type> <model params, at least 19 depending on model>\n"
            "    ^params from command line, 0 = yes _prev files, 1 = no\n"
            "- Usage 8: bj_mcmc 2/3 <data folder> <model type> <model params, at least 19 depending on model>\n"
            "    ^params from command line, 2 = yes _prev files, 3 = no\n"
            "\n"
            "Valid values for model type:\n"
            "- 0: model with just blob\n"
            "- 1: blob + EIC parameters\n"
            "- Other models not yet implemented\n";
    return 2;
  }

  run_models();
}