//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
#include <math.h>
#include "raytrace.h"
 
inline double calc_tbr(
                double Param_I[],
		double Rho,
		double Te,
		double ds,
		double Tbr,
		double *dtau) {
  /* 
   * Calculate brightness temperature
   */

  double nu_eff, Ne, Te32_inv, omega_p2, znsq;
  double alpha; 


  /* Access parameters by structure field names */
  struct param *prm = (void *) Param_I; 

  /*
   * Dimensionless plasma parameters
   */
  Ne =  Rho*(prm->ProtonMassInv);  /* Electron number density */
  Te32_inv = pow(Te,-1.5);            /* Te^(-3/2) to calc nu_eff */
  nu_eff = prm->Cnu*Ne*Te32_inv*log(220.0*Te*pow(Ne,-1./3.)); /*Ginzburg*/
  omega_p2 = (prm->e2_4pi_ovr_m)*Ne;       /* Plasma frequency squared */
  znsq = sqrt(prm->Omega2 - omega_p2);     /* sqrt(w^2 - wp^2) */

  /* 
   * Absorption coefficient 
   */
  alpha = omega_p2*nu_eff/(znsq*(prm->Omega)*(prm->c_light_cms)); 
 
  /*
   * Optical depth inctement
   */
  *dtau = alpha*ds*(prm->Rsun_cm);

  Tbr = Te + exp(-(*dtau))*(Tbr - Te);
  /* printf("alpha = %e, dtau = %e \n", alpha, *dtau); */
 
  return Tbr;
}
