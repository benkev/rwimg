//#include <iostream>
#include "simulation.h"
#include "advance_beam.h"
#include "streamer.h"
#include <cmath>
#include <stdio.h>
#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace std;

/*
* A quick function to check if there was an error with the last Kernel Call.
* Fill in the char * with a descriptive error message. This function terminates
* the program.
*/
void check_call(const char *func){
    cudaError_t error = cudaGetLastError();
    if(error!=cudaSuccess){
        printf("%s: %s\n",func,cudaGetErrorString(error));
        exit(1);
    }   
}

/*
* This destructor copies all of the data back from the GPU to the CPU. Then all
* of the memory on the GPU is cleaned up.
*/
Simulation::~Simulation(){

    Streamers_free(this->oStreamers);

    int nRay3 = this->nRay*3;

    CUDA_CALL(cudaMemcpy(this->Pos_ID_h, this->Pos_ID, sizeof(double)*nRay3,
                cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(this->Dir_ID_h, this->Dir_ID, sizeof(double)*nRay3,
                cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(this->DS_I_h, this->DS_I, sizeof(double)*this->nRay,
                cudaMemcpyDeviceToHost));


    if(this->rtmode>1)
        CUDA_CALL(cudaMemcpy(this->OpDepth_I_h, this->OpDepth_I,
                    sizeof(double)*this->nRay, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(this->Flags_I_h, this->Flags_I,
                sizeof(short)*this->nRay, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(this->Rho_I_h, this->Rho_I, sizeof(double)*this->nRay,
                cudaMemcpyDeviceToHost));

    CUDA_CALL(cudaMemcpy(this->GradRho_ID_h, this->GradRho_ID,
                sizeof(double)*nRay3, cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(this->PosPr_ID_h, this->PosPr_ID, sizeof(double)*nRay3,
                cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaMemcpy(this->DirPr_ID_h, this->DirPr_ID,
                sizeof(double)*this->nRay, cudaMemcpyDeviceToHost));

    if(this->rtmode==2)
        CUDA_CALL(cudaMemcpy(this->Tbr_I_h, this->Tbr_I, \
			     sizeof(double)*this->nRay,
			     cudaMemcpyDeviceToHost));

    CUDA_CALL(cudaMemcpy(this->DS_New_I_h, this->DS_New_I,
                sizeof(double)*this->nRay, cudaMemcpyDeviceToHost));

    if(this->rtmode==3){
        CUDA_CALL(cudaMemcpy(this->TbrIQUV_IP_h, this->TbrIQUV_IP, \
			     sizeof(double)
                    *this->nRay*4, cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaMemcpy(this->Bfield_ID_h, this->Bfield_ID,
                    sizeof(double)*nRay3, cudaMemcpyDeviceToHost));
    }

    CUDA_CALL(cudaMemcpy(this->DistToCrSurf_I_h, this->DistToCrSurf_I,
                sizeof(double)*this->nRay, cudaMemcpyDeviceToHost));

    if(nTracedPts > 0) {
      /* Traced Rays */

      CUDA_CALL(cudaMemcpy(TracedPts_I_h, TracedPts_I, 
	  sizeof(int)*this->nTracedPts, cudaMemcpyDeviceToHost));
      CUDA_CALL(cudaMemcpy(Trajectories_I_h, Trajectories_I, 
		 sizeof(double)*this->nTracedPts*3*this->nIter, 
			 cudaMemcpyDeviceToHost));
      CUDA_CALL(cudaMemcpy(LastStep_I_h, LastStep_I, 
		 sizeof(int)*this->nTracedPts, cudaMemcpyDeviceToHost));
      /* End Traced Rays */
    }

    CUDA_CALL(cudaFree(this->Pos_ID));
    CUDA_CALL(cudaFree(this->Dir_ID));
    CUDA_CALL(cudaFree(this->DS_I));
    if(this->rtmode==2)
        CUDA_CALL(cudaFree(this->Tbr_I));

    CUDA_CALL(cudaFree(this->OpDepth_I));
    CUDA_CALL(cudaFree(this->Flags_I));
    CUDA_CALL(cudaFree(this->Rho_I));

    CUDA_CALL(cudaFree(this->GradRho_ID));
    CUDA_CALL(cudaFree(this->PosPr_ID));
    CUDA_CALL(cudaFree(this->DirPr_ID));
    CUDA_CALL(cudaFree(this->DS_New_I));

    if(this->rtmode==3){
        CUDA_CALL(cudaFree(this->TbrIQUV_IP));
        CUDA_CALL(cudaFree(this->Bfield_ID));
    }
    CUDA_CALL(cudaFree(this->DistToCrSurf_I));

    CUDA_CALL(cudaFree(this->TracedPts_I));
    CUDA_CALL(cudaFree(this->Trajectories_I));
    CUDA_CALL(cudaFree(this->LastStep_I));
}

/*
* In the constructor, a copy of all the parameter fields is made on the device.
* Therefore there are two copies of each variable stored. The normal fields
* accessed by sim->field are the device copies. Fields accessed by
* sim_h->field_h are the fields on the CPU
*/
Simulation::Simulation(
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
        int rsph,
        double *theta,
        double *phi,
        double *orientation,
        double *density,
        double *baseStrength,
        double *stalkStrength,
        double *scale) {

  this->oStreamers = Streamers_new();        

  int i = 0;

  while (!isnan(theta[i])) {
    Streamers_makeStreamer(this->oStreamers,
			   1.0,
			   theta[i],
			   phi[i],
			   orientation[i],
			   density[i],
			   baseStrength[i],
			   stalkStrength[i],
			   scale[i],
			   scale[i],
			   scale[i]);
    i++;
  } 

  /*
   * Fill in the double ptraj[nTracedPts] table with the starting 
   * positions in the trajectory store array Trajectories_I[nTracedPts,nRay,3].
   */
  

  this->nIter=nIter;
  this->nRay=nRay;
  this->nTracedPts=nTracedPts;
  this->prm_h=prm;
  this->Pos_ID_h=Pos_ID; 
  this->Dir_ID_h=Dir_ID; 
  this->DS_I_h=DS_I; 
  this->Flags_I_h=Flags_I;
  this->Tbr_I_h=Tbr_I; 
  this->TbrIQUV_IP_h=TbrIQUV_IP; 
  this->OpDepth_I_h=OpDepth_I;
  this->Rho_I_h=Rho_I;
  this->GradRho_ID_h=GradRho_ID;
  this->Bfield_ID_h=Bfield_ID;
  this->PosPr_ID_h=PosPr_ID;
  this->DirPr_ID_h=DirPr_ID;
  this->DS_New_I_h=DS_New_I;
  this->DistToCrSurf_I_h=DistToCrSurf_I;
  this->TracedPts_I_h=TracedPts_I;
  this->Trajectories_I=Trajectories_I_h;
  this->LastStep_I=LastStep_I_h;
  this->rtmode=rtmode;
  this->scattering=scattering;
  this->rsph = rsph+1;

  int nRay3 = nRay*3;


  CUDA_CALL(cudaMalloc((void **) &this->Pos_ID, sizeof(double)*nRay3));
  CUDA_CALL(cudaMalloc((void **) &this->Dir_ID, sizeof(double)*nRay3));
  CUDA_CALL(cudaMalloc((void **) &this->DS_I, sizeof(double)*nRay));
  CUDA_CALL(cudaMalloc((void **) &this->Tbr_I, sizeof(double)*nRay));

  CUDA_CALL(cudaMalloc((void **) &this->OpDepth_I,
		       sizeof(double)*nRay));
  CUDA_CALL(cudaMalloc((void **) &this->Flags_I, sizeof(short)*nRay));
  CUDA_CALL(cudaMalloc((void **) &this->Rho_I, sizeof(double)*nRay));

  CUDA_CALL(cudaMalloc((void **) &this->GradRho_ID,
		       sizeof(double)*nRay3));
  CUDA_CALL(cudaMalloc((void **) &this->PosPr_ID,
		       sizeof(double)*nRay3));
  CUDA_CALL(cudaMalloc((void **) &this->DirPr_ID,
		       sizeof(double)*nRay3));
  CUDA_CALL(cudaMalloc((void **) &this->DS_New_I,
		       sizeof(double)*nRay));

  if(rtmode==3){
    CUDA_CALL(cudaMalloc((void **) &this->TbrIQUV_IP,
			 sizeof(double)*nRay*4));
    CUDA_CALL(cudaMalloc((void **) &this->Bfield_ID,
			 sizeof(double)*nRay3));
  }

  CUDA_CALL(cudaMalloc((void **) &this->DistToCrSurf_I,
		       sizeof(double)*nRay));

  if(nTracedPts > 0) {
    /* Traced Rays */
    CUDA_CALL(cudaMalloc((void **) &this->TracedPts_I,
			 sizeof(int)*nTracedPts));
    CUDA_CALL(cudaMalloc((void **) &this->Trajectories_I,
			 sizeof(double)*nTracedPts*nIter*3));
    CUDA_CALL(cudaMalloc((void **) &this->LastStep_I,
			 sizeof(int)*nTracedPts));
    /* End Traced Rays */
  }

  CUDA_CALL(cudaMalloc((void **) &this->prm, sizeof(struct param)));


  CUDA_CALL(cudaMemcpy(this->Pos_ID, this->Pos_ID_h,
		       sizeof(double)*nRay3, cudaMemcpyHostToDevice));


  CUDA_CALL(cudaMemcpy(this->Dir_ID, this->Dir_ID_h,
		       sizeof(double)*nRay3, cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(this->DS_I, this->DS_I_h, sizeof(double)*nRay,
		       cudaMemcpyHostToDevice));

  if(this->rtmode>1)
    CUDA_CALL(cudaMemcpy(this->OpDepth_I, this->OpDepth_I_h,
			 sizeof(double)*nRay, cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(this->Flags_I, this->Flags_I_h,
		       sizeof(short)*nRay, cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(this->Rho_I, this->Rho_I_h,
		       sizeof(double)*nRay, cudaMemcpyHostToDevice));

  CUDA_CALL(cudaMemcpy(this->GradRho_ID, this->GradRho_ID_h,
		       sizeof(double)*nRay3, cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(this->PosPr_ID, this->PosPr_ID_h,
		       sizeof(double)*nRay3, cudaMemcpyHostToDevice));
  CUDA_CALL(cudaMemcpy(this->DirPr_ID, this->DirPr_ID_h,
		       sizeof(double)*nRay, cudaMemcpyHostToDevice));

  if(this->rtmode==2)
    CUDA_CALL(cudaMemcpy(this->Tbr_I, this->Tbr_I_h,
			 sizeof(double)*nRay, cudaMemcpyHostToDevice));

  CUDA_CALL(cudaMemcpy(this->DS_New_I, this->DS_New_I_h,
		       sizeof(double)*nRay, cudaMemcpyHostToDevice));

  if(this->rtmode==3){
    CUDA_CALL(cudaMemcpy(this->TbrIQUV_IP, this->TbrIQUV_IP_h, 
			 sizeof(double*)*nRay*4, cudaMemcpyHostToDevice));

    CUDA_CALL(cudaMemcpy(this->Bfield_ID, this->Bfield_ID_h,
			 sizeof(double)*nRay3, cudaMemcpyHostToDevice));
  }

  CUDA_CALL(cudaMemcpy(this->DistToCrSurf_I, this->DistToCrSurf_I_h,
		       sizeof(double)*nRay, cudaMemcpyHostToDevice));

  if(nTracedPts > 0) {
    /* Traced Rays */
    CUDA_CALL(cudaMemcpy(this->TracedPts_I, this->TracedPts_I_h,
			 sizeof(int)*nTracedPts, cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(this->Trajectories_I, this->Trajectories_I_h,
			 sizeof(double)*nTracedPts*nIter*3, 
			 cudaMemcpyHostToDevice));
    CUDA_CALL(cudaMemcpy(this->LastStep_I, this->LastStep_I_h,
			 sizeof(int)*this->nTracedPts, cudaMemcpyHostToDevice));
    /* End Traced Rays */
  }
  CUDA_CALL(cudaMemcpy(this->prm, this->prm_h, sizeof(struct param),
		       cudaMemcpyHostToDevice));
}



/*
* This is the interface between C and C++ functions. The C function calls this
* which then creates the necessary simulation object and calls the methods to
* run the simulation.
*/
extern "C" void make_simulation(
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
        int rsph,
        double *theta,
        double *phi,
        double *orientation,
        double *density,
        double *baseStrength,
        double *stalkStrength,
        double *scale){




            Simulation sim (nIter,
                    nRay,
		    nTracedPts,
                    prm,
                    Pos_ID,
                    Dir_ID,
                    DS_I,
                    Flags_I,
                    Tbr_I,
                    TbrIQUV_IP,
                    OpDepth_I,
                    Rho_I,
                    GradRho_ID,
                    Bfield_ID,
                    PosPr_ID,
                    DirPr_ID,
                    DS_New_I,
                    DistToCrSurf_I,
		    TracedPts_I,
	            Trajectories_I,
	            LastStep_I,
                    rtmode,
                    scattering,
                    rsph,
                    theta,
                    phi,
                    orientation,
                    density,
                    baseStrength,
                    stalkStrength,
                    scale);





            sim.trace();

        }
/*
* This is the function that calls the advance_beam function. It makes a copy of
* the simulation object onto the GPU and then calls advance_beam with a
* reference to that GPU copy of simulation.
*/
void Simulation::trace(){
    Simulation *device_sim;

    CUDA_CALL(cudaMalloc((void **) &device_sim, sizeof(Simulation)));
    CUDA_CALL(cudaMemcpy(device_sim, this, sizeof(Simulation),
                cudaMemcpyHostToDevice));

    advance_beam(this, device_sim);

    CUDA_CALL(cudaMemcpy(this, device_sim, sizeof(Simulation),
                cudaMemcpyDeviceToHost));
    CUDA_CALL(cudaFree(device_sim));
}


