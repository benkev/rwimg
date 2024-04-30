#ifndef CALC_INC
#define CALC_INC


double calc_tbr(
        struct param *prm, 
		double Rho,
		double Te,
		double ds,
		double Tbr,
		double *dtau) ;

void calc_tbriquv(
        struct param *prm,
		double Dir_D[3], 
		double Rho,
		double Bfield_D[3],
		double Te,
		double ds,
		double TbrIQUV_P[4],
		double *dtau,
		int thread_num) ;
#endif
