//
// raytrace.h
//    Definitions for beam_path(), beam_intensity(), and write_rwimage
//
#include <stdint.h>
#include <pthread.h>

// Bits for Flags_I
#define INACTIVE  0x0001   // The ray is not being processed  
#define SHARP     0x0002   // The ray is sharp (close to normal) 
#define PENETR    0x0004   // Eps < 0; Critical surface search
#define WASREFL   0x0008   // First step after reflection
#define BADRAY    0x0010   // The ray is bad 

#define pi 3.1415926535897931

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

#define BIT_IS_ON(PATTERN, iRay)     Flags_I[iRay] & (PATTERN)
#define BIT_IS_OFF(PATTERN, iRay)   ~Flags_I[iRay] & (PATTERN)
#define SET_BIT(PATTERN, iRay)    Flags_I[iRay] |= (PATTERN)
#define CLEAR_BIT(PATTERN, iRay)  Flags_I[iRay] &= ~(PATTERN)

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

struct TraceBits { /* Bit flags of the traced parameters */
  unsigned int Pos       : 1;
  unsigned int Dir       : 1;
  unsigned int Arclen    : 1;
  unsigned int AngBD     : 1;
  unsigned int Dist      : 1;
  unsigned int Rho       : 1;
  unsigned int gRho      : 1;
  unsigned int gRhoMagn  : 1;
  unsigned int Bfield    : 1;
  unsigned int BMagn     : 1;
  unsigned int Tbr       : 1;
  unsigned int Stokes    : 1;
};

struct BeamData {
  /* Array parameters */
  double *Param_I;
  double *Pos_ID;
  double *Dir_ID;
  double *DS_I; 
  double *Arclen_I;
  double *CosBD_I;
  double *Dist_I;
  double *Tbr_I; 
  double *OpDepth_I; 
  short  *Flags_I;
  double *Rho_I; 
  double *GradRho_ID;
  double *PosPr_ID; 
  double *DirPr_ID; 
  double *DS_New_I;
  double *TbrIQUV_IP; 
  double *Bfield_ID;
  double *DistToCrSurf_I;
  double *TrcPos_ID; 
  double *TrcDir_ID; 
  double *TrcArclen_I;
  double *TrcAngBD_I;
  double *TrcRho_I; 
  double *TrcGradRho_ID; 
  double *TrcGRhoMagn_I; 
  double *TrcDist_I;
  double *TrcBfield_ID;
  double *TrcBMagn_I;
  double *TrcTbrIQUV_IP; 
  double *TrcTbr_I;
  int32_t *TrcIrays_I; 
  int32_t *TrcLast_I;
  void (* plasma_parameters)();
  int *thiter;    /* Iteration counters for individual threads */
  int *thnactive; /* Numbers of active rays for individual threads */
  pthread_mutex_t *mutex;
  struct TraceBits *toTrace; /* Trace bit indicators */
  /* Elemental input variables */
  double RIntgSph;
  double pwr10; /* Negative power of 10 to compare with the nactive/nRay */
  double exp10; /* Negative exponential of 10  */
  double dexp; /* Delta of the 10's power */
  clock_t tstart; /* Start time in machine clocks */ 
  int diter; /* Delta of the iteration number between prints */
  int niterprinted;  /* Number of iterations already printed */
  int dtime; /* Time interval between diagnostic prints */
  int ntintervprinted; /* Number of time intervals printed */
  int nacthr;  /* Number of active threads */
  int printed;  /* 0|1: diagnostics has been printed */
  int nRay;
  int nIter;
  int rtmode;
  int scattering;
  int ntraj; /* Number of trajectories traced */
  int npmax; /* Max number of points, i.e. max trajectory length */
  /* int allThreadsHaveWorked; */
};


inline double dot_product(double a[3], double b[3]);
inline void cross_product(double a[3], double b[3],double c[3]);
double sum(double a[], int n);
double sum_squares(double vec[], int n);
void print1d(int n, double a[n]);
void print2d(int m, int n, double a[m][n]);
double minval(double a[], int n);
double maxval(double a[], int n);
inline double vmagn(double a[], int n); 
inline double v2magn(double a[]);
inline double v3magn(double a[]);
void linsolve(const int N, double a[N][N], double b[N], double x[N]);
void linsolve2(const int N, const int M, 
	       double a[N][N], double b[N][M], double x[N][M]);
void mvmul(int N, double a[N][N], double x[N], double b[N]);
void mmul(int const N, int const M, int const L, 
	  double a[N][M], double b[M][L], double r[N][L]);
void minv(const int N, double a[N][N]);
void minv2(const int N, double a[N][N], double ainv[N][N]);
void calc_rdir(double obs[3], double targ[2], double dir[3]);
int  calc_rpos(double obs[3], double dir[3], double rsph, double pos[3]);
void mcalc_rdir(double obs[3], double xruler[], double yruler[], 
		int nx, int ny,	double dir[][3]);
int mcalc_rsph(double obs[3], double dir[][3], double rsph, int nRay, 
	       short isec[], double pos[][3]);
void log10space(double start, double stop, long np, double sp[]);
void logspace(int base, double start, double stop, long np, double sp[]);


void advance_beam(
	int start,
	int stop,
	int thread_num,
	int nRay,
        double Param_I[],
        double Pos_ID[nRay][3], 
	double Dir_ID[nRay][3], 
	double Arclen_I[nRay], 
	double CosBD_I[nRay], 
	double DS_I[nRay], 
        void (*plasma_params)(),
        short Flags_I[nRay],
        double Tbr_I[nRay], 
        double TbrIQUV_IP[nRay][4], 
        double OpDepth_I[nRay],
        double Rho_I[nRay],
        double GradRho_ID[nRay][3],
        double Bfield_ID[nRay][3],
        double PosPr_ID[nRay][3],
        double DirPr_ID[nRay][3],
        double DS_New_I[nRay],
        double DistToCrSurf_I[nRay],
	int rtmode,
	int scattering);


void plasma_density(int nRay, 
		    double Param_I[], 
		    double Pos_ID[nRay][3], 
		    double Rho_I[nRay], 
                    double GradRho_ID[nRay][3], 
		    double DS_I[nRay], 
		    short Flags_I[nRay]);

void plasma_params(int start,
		   int stop,
		   int nRay, 
		   double Param_I[], 
		   double Pos_ID[nRay][3], 
		   double Rho_I[nRay], 
		   double GradRho_ID[nRay][3], 
		   double Bfield_ID[nRay][3],
		   double DS_I[nRay], 
		   short Flags_I[nRay]);

inline double calc_tbr(
                double Param_I[],
		double Rho,
		double Te,
		double ds,
		double Tbr,
		double *dtau);

inline void calc_tbriquv(
                double Param_I[],
		double Dir_D[3], 
		double Rho,
		double Bfield_D[3],
		double Te,
		double ds,
		double TbrIQUV_P[4],
		double *dtau,
		double *cosbd,
		int thread_num);


void beam_intens_traj(double XyzEarth_D[3], double RadioFrequency, 
		      double ImageRange_I[4], double rIntegration, 
		      int nXPixel, int nYPixel, 
		      double Intensity_II[nYPixel][nXPixel], 
		      int const nTrajMax, int const lenTrajMax, 
		      int nTraj, int SelTraj_II[2][nTrajMax], 
		      int LenTraj_I[nTrajMax], 
		      double Traj_DII[3][nTrajMax][lenTrajMax],
		      double Tb_II[nTrajMax][lenTrajMax]);

inline double saito(double r, double cosTh);

inline double density_only(int iRay, double Param_I[], 
			   double x, double y, double z);

