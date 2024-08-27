/*
Slightly modified version of bj02.h for use in MCMC code; modified by Sarah Youngquist in Feb 2022
See bj_core.cpp for an explanation of the changes.
*/
/** @file */

#ifndef BLAZARS_MCMC_BJ02_COPY_H
#define BLAZARS_MCMC_BJ02_COPY_H


using namespace std;

#include <time.h>

namespace bj_core02 {



// FLAGS
    const int PRINT = 0; //!< 1=True, 0=False
    const int DEBUG = 0; //!< 1=True, 0=False
//#warning "DEUBG flag enabled"

// PHYSICAL CONSTANTS c.g.i.

    const double q0     = 0.5;
    const double pc     = 3.086    * 1.0e+18; //!< [cm]
    const double H0SI   = 70.;                //!< SI:  km/s/Mpc
    const double H0     = H0SI*1.e5/(1.e6*pc);//!< cgs: s-1
    const double sig_T  = 6.652453 * 1.0e-25; //!< [cm^2]
    const double m_e    = 9.109558 * 1.0e-28; //!< [g]
    const double c      = 2.997924 * 1.0e+10; //!< [cm / s]
    const double h      = 6.626296 * 1.0e-27; //!< [erg * s]
    const double sig    = 5.66961  * 1.0e-05; //!< [erg / (s * cm^2 * K^4)]
    const double m_p    = 1.672623 * 1.0e-24; //!< [g]
    const double e      = 4.803250 * 1.0e-10; //!< [esu]
    const double keV    = 2.417965 * 1.0e+17; //!< [Hz]
    const double G      = 6.6726   * 1.0e-08; //!< [dyn / (g^2 * cm^2)]
    const double k_B    = 1.380622 * 1.0e-16; //!< [erg / K]
    const double M_sol  = 1.989    * 1.0e+33; //!< [g]
    const double eV_K   = 11.6048  * 1.0e+3;  //!< [K]
    const double hour   = 3600.0;             //!< [s]
    const double day    = 24.0 * hour;        //!< [s]


// PHYSICAL CONSTANTS S.I.

    const double PLANCK = 6.6260693e-34; //!< J.s
    const double eV     = 1.602e-19;     //!< J
    const double PARSEC     = pc;	//!< cm
    const double ERG        = 1.e-7;	//!< J
    const double MELEC      = 0.511e6;   //!< eV
    const double HZ_PER_EV  = 2.145e14;
    const double LIGHTSPEED = c; //!< cm/s
    const double SIGMAT     = sig_T; //!<cm2
    const double HMASS      = 1.673e-27;	//!< kg
    const double MSOL = 1.989e30;	//!< kg
    const double ALPHA = 0.0073;	//!< fine structure const
    const double MILLIBARN = 1.e-27; //!< cm^2
    const double CHARGE = 1.602e-19; //!<Coulomb
    const double MELECKG = m_e; //!<g
    const double MPION = 135.e6;  //!<eV
    const double MPROTON = 938.272e6;  //!<eV

// CALCULATION PRECISION

    const int    SYN_PREC1       = 15;      //!< INTEGRATION PRECISION FOR SYNCH. SPEC.
    const int    SYN_PREC2       = 15;      //!< INTEGRATION PRECISION FOR SYNCH. SPEC.
    const int    ABS_PREC1       = 15;      //!< INTEGRATION PRECISION FOR ABSOR. SPEC.
    const int    ABS_PREC2       = 15;      //!< INTEGRATION PRECISION FOR ABSOR. SPEC.
    const int    COM_PREC1       = 15;      //!< INTEGRATION PRECISION FOR I.COM. SPEC.
    const int    COM_PREC2       = 15;      //!< INTEGRATION PRECISION FOR I.COM. SPEC.
    const int    G_MAX           = 1000;    //!< MAXIMUM NUMBER OF GAMMA POINTS
    const int    NU_MAX          = 5000;     //!< MAXIMUM NUMBER OF SPECTRAL POINTS


// TIME VARIABLES

    extern time_t       T_START, T_END;
    extern char         PARA_FN[256];

// GLOBAL VARIABLES


    const int    G_DIM=100;     //!< CURRENT NUMBER OF GAMMA POINTS
    const double N_MIN=1.0e-15; //!< MINIMAL ELECTRONS DENSITY
    const double N_MAX=50;  //!< MAXIMAL ELECTRONS DENSITY

//extern char  prefix[256];

    const int    NU_DIM_MAX      = 300;
    const int    SL_DIM_MAX      = 200;


// PHYSICAL PARAMETERS

    extern double       z;
    extern double       H_0;
    extern double       THETA, D_L;
    extern double       DOP_B, LOR_B, V_B, V_B_APP;//!<DOP_BB,
    extern double       R_src,R_blr;
    extern double       L_src, L_nuc;
    extern double       B;
    extern double       K1;
    extern double       K2;
    extern double       N1;
    extern double       N2;
    extern double       GAMMA_MIN;
    extern double       GAMMA_BRK;
    extern double       GAMMA_MAX;
    extern double       T_BB;
    extern double       tau;  //!< fraction of L_nuc scattered/rerocessed isotropically (EIC)
    extern double       B_0;
    const double        n_B = 1.0; //!< slope of the magnetic field along the jet
    extern double       N_0;
    extern double       n_n;
    extern double       n_N;
    extern double       n_G;
    extern double       D_b, D_b_src_j;
    extern double       U_B,U_e;
    extern double       jj1,kk1;
    extern double       L_jet_eic;
    extern double       L_tor, T_BB_tor;
    extern double       Tcool_vhe, Tcool_radio;

// TRANSFORMATION PARAMETERS

    extern double DOP_JET, LOR_JET, V_JET, V_JET_APP, DOP_B_J, LOR_B_J, V_B_J;

// GEOMETRY PARAMETERS

//const int           SL_DIM_MAX = 200;
    extern int          SL_DIM;
    extern double       PHI, PHI_src;
    extern double       THETA_src_j;
    extern double       J_LEN, J_LEN_src;
    extern double       A;
    extern double       X_MIN;
    extern double       X_MAX;
    extern double       X_OUT;
    extern double       Y_MIN;
    extern double       Y_MAX;
    extern double       R_MOY;
    extern double       Xcut;
    extern double       antisym;
    extern double       Sprev;
    extern double       Sum_S;

// OTHER PARAMETERS

    extern int          SL_CUR;
    extern double       GAMMA_MAX_0;
    extern int          null0;
    extern double       V_exp, DEL_Tph, DEL_Tj;



// SUBROUTINES ******************************************************************
    extern double N_e(double gamma);


// METHODS
    void description();
    // Sarah Youngquist added
    int run_models(); // moved prev code from main here
    int load_params(char* name);
    int load_params_from_list(char** list, int model_type, int starting_index, int list_length);
};

#endif //BLAZARS_MCMC_BJ02_COPY_H
