#ifndef SIMULATION_INCLUDED
#define SIMULATION_INCLUDED



#define INACTIVE  0x0001   // The ray is not being processed  
#define SHARP     0x0002   // The ray is sharp (close to normal) 
#define PENETR    0x0004   // Eps < 0; Critical surface search
#define WASREFL   0x0008   // First step after reflection
#define BADRAY    0x0010   // The ray is bad 
#define TRACED    0x0020   // The ray is traced (the ray traces are recorded)

#define PI 3.1415926535897931

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

#define BIT_IS_ON(PATTERN, iRay)     sim->Flags_I[iRay] &   (PATTERN)
#define BIT_IS_OFF(PATTERN, iRay)   ~sim->Flags_I[iRay] &   (PATTERN)
#define SET_BIT(PATTERN, iRay)       sim->Flags_I[iRay] |=  (PATTERN)
#define CLEAR_BIT(PATTERN, iRay)     sim->Flags_I[iRay] &= ~(PATTERN)



#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
printf("Error at %s:%d\n",__FILE__,__LINE__); \
exit(EXIT_FAILURE);}} while(0)


#define NPARAM 41
struct param {
  double DeltaS;          /* Initial ray path increment in solar radii */ 
  double Freq;            /* Frequency, Hz */
  double Omega;           /* Frequency, rad/s */
  double Freq2;           /* Frequency squared, Hz^2 */
  double Omega2;          /* Frequency squared, (rad/s)^2 */
  double k0;              /* k0 = Omega/c_light_cms, wave number in vacuum */
  double k0_Rsun;         /* k0 = Omega/c_light_Rsun, wave # in rad/Rsun */
  double RhoCr;           /* g/cm^3, critical density at given Freq */
  double RhoCrInv;        /* 1/RhoCr, RhoCr inverted */
  double Tol;             /* Tolerance - maximum radians per one step  */
  double Tol2;            /* Tol^2, tolerance squared */
  double AbsMinStep;      /* Limit to adaptive step decrease */
  double TolEps;          /* Tolerance for dielectric permittivity */
  double CntMax;          /* Upper limit to Newton iterations number */
  double h_chromo_km;     /* km, Chromosphere thickness */
  double r_chromo;        /* (Rsun_km + h_chromo_km - 1000.)/Rsun_km */
  double r_chro_cor;      /* (Rsun_km + h_chromo_km)/Rsun_km */
  double r_corona;        /* (Rsun_km + h_chromo_km + 1000.)/Rsun_km */
  double c_light_cms;     /* cm/s, Speed of light */
  double ProtonChargeCGSe;/* StatCoulombs, CGSe */
  double ProtonMass_g;    /* g */
  double ProtonMassInv;   /* g^(-1) */
  double ElectronMass_g;  /* g */
  double Boltzmannk_ergK; /* erg/K, k = 1.38065042424e-16 erg/K */
  double Rsun_cm;         /* cm, Solar radius */
  double Rsun_km;         /* km, Solar radius */
  double Te_corona_K;     /* 1.0e6 K */
  double Te_chromo_K;     /* 3.0e4 K */
  double AU_m;            /* 149597870700 m, Sun-earth distance in meters */
  double AU_Rsun;         /* 215.097 Rsun, Sun-earth distance in solar radii */
  double Msun_G;          /* 1.04 G/Rsun^3, solar dipole field at equator */
  double dirMx;           /* 0.0, Solar dipole x direction, CGI */
  double dirMy;           /* 0.0, Solar dipole y direction, CGI */
  double dirMz;           /* 1.0, Solar dipole z direction, CGI */
  double e_ovr_mcw2;      /* e/(m_e*c*omega), to calculate u_plasma */
  double e2_4pi_ovr_m;    /* 4*pi*e^/m_e, to calculate plasma freq^2 */
  double e2_4pi_ovr_mw2;  /* 4*pi*e^/(m_e*omega^2), to calculate v_plasma */
  double twokf2c2;        /* I = 2kf^2/c^2 * Tb, Rayleigh-Jeans factor */
  double lnLambda_13p7;   /* 13.7*Coulomb logarithm, ~13.7*20 for corona */
  double Cnu;             /* Coef. at Ginzburg's nu_eff */
  double callcount;       /* Counter of calls to advance_beam() */
};

#ifdef __cplusplus

#include "streamer.h"
#include <cstdlib>

void check_call(const char *func);

class Simulation {

    public:
        
        int nIter;
        int nRay;
	int nTracedPts;

	/*
	 * DEVICE pointers to data
	 */
        struct param *prm;
        double *Pos_ID; 
        double *Dir_ID; 
        double *DS_I; 
        short *Flags_I;
        double *Tbr_I; 
        double *TbrIQUV_IP; 
        double *OpDepth_I;
        double *Rho_I;
        double *GradRho_ID;
        double *Bfield_ID;
        double *PosPr_ID;
        double *DirPr_ID;
        double *DS_New_I;
        double *DistToCrSurf_I;
        int    *TracedPts_I;
	double *Trajectories_I;
	int    *LastStep_I;              

        int rtmode;
        int scattering;
        double rsph;
        Streamers_T oStreamers;

	/*
	 * HOST Pointers to data 
	 */
        double *Pos_ID_h, *Dir_ID_h, *DS_I_h, *Tbr_I_h;
        double *OpDepth_I_h;
        short  *Flags_I_h;
        double *Rho_I_h, *GradRho_ID_h;
        double *PosPr_ID_h, *DirPr_ID_h, *DS_New_I_h;
        double *TbrIQUV_IP_h, *Bfield_ID_h;
        double *DistToCrSurf_I_h;
        int    *TracedPts_I_h;
	double *Trajectories_I_h;
	int    *LastStep_I_h;              

        struct param *prm_h;

        Simulation(		
                int nIter,
                int nRay,
		int nTracedPts,
                struct param *prm,
                double *Pos_ID, 
                double *Dir_ID, 
                double *DS_I, 
                short *Flags_I,
                double *Tbr_I, 
                double *TbrIQUV_IP, 
                double *OpDepth_I,
                double *Rho_I,
                double *GradRho_ID,
                double *Bfield_ID,
                double *PosPr_ID,
                double *DirPr_ID,
                double *DS_New_I,
                double *DistToCrSurf_I,
		int    *TracedPts_I,
		double *Trajectories_I,
		int    *LastStep_I,
                int rtmode,
                int scattering,        
                double rsph,
                double *theta,
                double *phi,
                double *orientation,
                double *density,
                double *baseStrength,
                double *stalkStrength,
                double *scale);

        ~Simulation();

        void trace();

};
#endif

#ifdef __cplusplus
extern "C"{ 
#endif
void make_simulation(
        int nIter,
        int nRay,
	int nTracedPts,
        struct param *prm,
        double *Pos_ID, 
        double *Dir_ID, 
        double *DS_I, 
        short *Flags_I,
        double *Tbr_I, 
        double *TbrIQUV_IP, 
        double *OpDepth_I,
        double *Rho_I,
        double *GradRho_ID,
        double *Bfield_ID,
        double *PosPr_ID,
        double *DirPr_ID,
        double *DS_New_I,
        double *DistToCrSurf_I,
	int    *TracedPts_I,
	double *Trajectories_I,
	int    *LastStep_I,
        int rtmode,
        int scattering,
        double rsph,
        double *theta,
        double *phi,
        double *orientation,
        double *density,
        double *baseStrength,
        double *stalkStrength,
        double *scale);

#ifdef __cplusplus
};
#endif


#endif
