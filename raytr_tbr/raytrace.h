//
// raytrace.h
//    Definitions for beam_path(), beam_intensity(), and write_rwimage
//

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

double inline dot_product(double a[3], double b[3]);
void inline cross_product(double a[3], double b[3],double c[3]);
double sum(double a[], int n);
double sum_squares(double vec[], int n);
void print1d(int n, double a[n]);
void print2d(int m, int n, double a[m][n]);
double minval(double a[], int n);
double maxval(double a[], int n);

void linsolve(const int N, double a[N][N], double b[N], double x[N]);

void linsolve2(const int N, const int M, 
	       double a[N][N], double b[N][M], double x[N][M]);

void mvmul(int N, double a[N][N], double x[N], double b[N]);

void minv(const int N, double a[N][N]);

void beam_path(void (* Get_Plasma_Density)(), int nRay, 
	       short ExcludeRay_I[nRay], double Position_DI[3][nRay],
	       double Slope_DI[3][nRay], double DeltaS_I[nRay], 
	       double ToleranceInit, double DensityCr, 
	       double Intensity_I[nRay], short RayFlag_I[nRay], 
	       short *NewEntry, double RadioFrequency, 
	       double OpticalDepth_I[nRay], FILE *fh);

void plasma_density(int nRay, double Pos_DI[3][nRay], 
		    double Rho_I[nRay], double GradRho_DI[3][nRay],
		    double DS_I[nRay], short Flags_I[nRay]);

void beam_intensity(double XyzEarth_D[3], double RadioFrequency, 
		    double ImageRange_I[4], double rIntegration, 
		    int nXPixel, int nYPixel, 
		    double Intensity_II[nYPixel][nXPixel]);

void beam_intens_traj(double XyzEarth_D[3], double RadioFrequency, 
		      double ImageRange_I[4], double rIntegration, 
		      int nXPixel, int nYPixel, 
		      double Intensity_II[nYPixel][nXPixel], 
		      int const nTrajMax, int const lenTrajMax, 
		      int nTraj, int SelTraj_II[2][nTrajMax], 
		      int LenTraj_I[nTrajMax], 
		      double Traj_DII[3][nTrajMax][lenTrajMax],
		      double Tb_II[nTrajMax][lenTrajMax]);

double inline saito(double r, double cosTh);

double inline density_only(double x, double y,double z);
