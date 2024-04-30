#ifndef PLASMA_INCLUDED
#define PLASMA_INCLUDED

struct plasma_constants{

    double r_corm16;
    double r_corm6; 
    double r_corm2d5; 
    double r_corm17; 
    double r_corm7; 
    double r_corm3d5; 

    double A[16];
    double rhsa[4];
    double a[4];

};
__device__ void mvmul(int N, double *a, double x[], double b[]);
__global__ void plasma_parameters(Simulation *simulation);
__device__ double density_only(struct param *prm, Streamers_T oStreamers, 
        double x, double y, double z);
#endif
