#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "raytrace.h"
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>

#define nTrajMax   120   // Maximum number of selected trajectories
#define lenTrajMax 5000  // Maximum lengths of selected trajectories
#define nXPixel    4   // Image plane X dimension in pixels
#define nYPixel    4   // Image plane Y dimension in pixels


double const cPi = 3.1415926535897931;
double const ProtonChargeSGSe = 4.8e-10;      // StatCoulombs, SGSe
double const cProtonMass = 1.672621636E-24;   //g
double const cElectronMass = 9.10938215E-28;  //g
double const cSolarRadius = 6.95E10;         //cm
double const h_chromo = 10000;         //cm


double const DeltaS = 0.1;         // in R_sun units

int main(int argc, char *argv[]) 
{
  //double REarth_D[3] = {150.0, 145.7, 50.0}; // SGI earth position in R_sol
  double REarth_D[3] = {215.0, 0.0, 0.0};   // SGI earth position in R_sol
  double Freq = 100.0E6;            // Hz
  /* ImageRange is (x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper) */
  double ImageRange_I[4] = {-2.0, -2.0, 2.0, 2.0}; 
  double rIntegration = 10.0;  // Radius of "integration sphere" in solar radii
  static double Intensity_II[nYPixel][nXPixel];  // Sun radioimage at Freq, Hz
  int i, j, nTraj;
  static int SelTraj_II[2][nTrajMax]; // (y,x) indices of selected trajectories
  static int LenTraj_I[nTrajMax];    // Lengths of calculated trajectories 
  static double Traj_DII[3][nTrajMax][lenTrajMax]; // The returned trajectories
  static double Tb_II[nTrajMax][lenTrajMax]; // Brightness temperatures 
  int iTraj, lenTraj;
  double dy;
  char fname[100]; 
  /* char fnpart[100]; */
  char cmd[100];
  FILE *fh;

  clock_t cstart, cend;
  double elapsed;
  time_t time1; 
  /*  time_t time2; */


  cstart = clock();
  time(&time1);

  printf("cstart = %g\n", (double)cstart);

  //SelTraj_II[0][0] = 3; SelTraj_II[1][0] = 4;   // i1, j1
  //SelTraj_II[0][1] = 10; SelTraj_II[1][1] = 10; // i2, j2
  printf("\nSelected pixels for trajectories:\n");
  nTraj = 0;
  //for (i = 0; i < nYPixel; i += 10) {
  //  SelTraj_II[0][nTraj] = i;         // Select i each 10th
  //  SelTraj_II[1][nTraj] = nXPixel/2; // Select j at the middle vertical line
  //  printf("%3i %3i\n", SelTraj_II[0][nTraj], SelTraj_II[1][nTraj]);
  //  if (++nTraj == nTrajMax) break;
  //}
  for (i = 0; i < 120; i++) { // Denser near the center
    SelTraj_II[0][nTraj] = i + 40;  //nYPixel/3;        
    SelTraj_II[1][nTraj] = nXPixel/2; // Select j at the middle vertical line
    printf("%3i %3i\n", SelTraj_II[0][nTraj], SelTraj_II[1][nTraj]);
    if (++nTraj == nTrajMax) break;
  }

  printf("nTraj = %i\n", nTraj);
  
  if (argc > 1) {
    Freq = 1e6*atof(argv[1]);
  }
  printf("\nFreq = %g\n", Freq);

  
  beam_intens_traj(REarth_D, Freq, ImageRange_I, 
		   rIntegration, nXPixel, nYPixel, Intensity_II, 
		   nTrajMax, lenTrajMax, nTraj, SelTraj_II, 
		   LenTraj_I, Traj_DII, Tb_II);

  fh = fopen("image.txt", "w");
  fprintf(fh, "#\n# Image generated via raytracing\n");
  fprintf(fh, "#\n# Radio wave frequency, Hz:\n");
  fprintf(fh, "%12.4e\n", Freq);
  fprintf(fh, "#\n# Image lower-left and upper-right points,\n");
  fprintf(fh, "# in the SGI coordinates, x0, y0, x1, y1,\n");
  fprintf(fh, "# solar radius units:\n");
  for (j = 0; j < 4; j++) fprintf(fh, "%g ", ImageRange_I[j]);
  fprintf(fh, "\n#\n");
  fprintf(fh, "# Image Y and X sizes in pixels:\n");
  fprintf(fh, "%i %i\n", nYPixel, nXPixel);
  fprintf(fh, "#\n# Data:\n");

  
  for (i = 0; i < nYPixel; i++) {
    for (j = 0; j < nXPixel; j++) {
      fprintf(fh, "%g ", Intensity_II[i][j]);
    }
    fprintf(fh, "\n");
  }

  fclose(fh);

  //return;
  
  //
  // Save the selected ray trajectories as a table
  // nTraj*SUM(LenTraj_I) lines by 3 columns of X Y Z coordinates
  // 
  //  lenTraj = maxval(LenTraj_I, nTraj);
  //		      int nTraj, int lenTrajMax, int SelTraj_II[2][nTraj], 
  //		      int LenTraj_I[nTraj], 
  //                  double Traj_DII[3][nTraj][lenTrajMax]) 
  printf("#\n# Ray trajectories generated via raytracing\n");
  printf("#\n# Radio wave frequency, Hz:\n");
  printf("%12.4e\n", Freq);
  printf("#\n# Lengths of rays stored:\n");
  for (j = 0; j < nTraj; j++) printf("%i ", LenTraj_I[j]);
  printf("\n#\n");
  printf("# Data:\n");

  fh = fopen("traj.txt", "w");

  fprintf(fh, "#\n# Ray trajectories generated via raytracing\n");
  fprintf(fh, "#\n# Radio wave frequency, Hz:\n");
  fprintf(fh, "%12.4e\n", Freq);
  fprintf(fh, "#\n# Lengths of rays stored:\n");
  for (j = 0; j < nTraj; j++) fprintf(fh, "%i ", LenTraj_I[j]);
  fprintf(fh, "\n#\n");
  fprintf(fh, "# Data:\n");

  for (iTraj = 0; iTraj < nTraj; iTraj++) {
    lenTraj = LenTraj_I[iTraj];
    for (i = 0; i < lenTraj; i++) { 
      for (j = 0; j < 3; j++) fprintf(fh, "%12.4e ", Traj_DII[j][iTraj][i]);
      fprintf(fh, "\n");
    }
  }

  fclose(fh);

  fh = fopen("tbr.txt", "w");

  fprintf(fh, "#\n# Brightness temperatures generated via raytracing\n");
  fprintf(fh, "#\n# Radio wave frequency, Hz:\n");
  fprintf(fh, "%12.4e\n", Freq);
  fprintf(fh, "#\n# Lengths of data stored:\n");
  for (j = 0; j < nTraj; j++) fprintf(fh, "%i ", LenTraj_I[j]);
  fprintf(fh, "\n#\n");
  fprintf(fh, "# Data:\n");

  for (iTraj = 0; iTraj < nTraj; iTraj++) {
    lenTraj = LenTraj_I[iTraj];
    for (i = 0; i < lenTraj; i++) { 
      fprintf(fh, "%12.4e ", Tb_II[iTraj][i]);
      fprintf(fh, "\n");
    }
  }

  fclose(fh);

  //=============================================
  //char fname[100], fnpart[100];
  //double dy;

  dy = (ImageRange_I[3] - ImageRange_I[1])/nYPixel;
  
  
  //strcpy(fname, "tb_");                   // fname = "tb_";
  //sprintf(fnpart, "%i", (int)(Freq/1e6)); // fnpart = (str)Freq/1e6; (in MHz)
  sprintf(fname, "tb_%i.txt", (int)(Freq/1e6)); // fname = Freq/1e6; (in MHz)
  //sprintf(fnpart, "_%i.txt", (int)(Freq/1e6)); // fnpart = Freq/1e6; (in MHz)
  //strcat(fname, fnpart);                  // fname = "tb_" + "200";
  //strcpy(fnpart, ".txt");                 // fnpart = ".txt";
  //strcat(fname, fnpart);                  // fname = "tb_200" + ".txt";

  fh = fopen(fname, "w");
  
  for (i = nYPixel/2; i < nYPixel; i++) {
    fprintf(fh, "%g %g\n", dy/2.0 + (float)(i-nYPixel/2)*dy, 
	    Intensity_II[i][nXPixel/2]);
  }
  //  fprintf(fh, "\n");
  
  fclose(fh);
  
 
  cend = clock();
  elapsed = ((double)(cend - cstart))/CLOCKS_PER_SEC;
  printf("\n\nCLOCKS_PER_SEC = %ld;\n", CLOCKS_PER_SEC);
  printf("Elapsed time: %3g sec = %02g:%02g\n", 
  	 elapsed, floor(elapsed/60.), fmod(elapsed,60.));

  //time(&time2);
  //elapsed = difftime (time2, time1);
  //printf("Elapsed time: %3g sec = %02g:%02g\n\n\n\n", 
  //	 elapsed, floor(elapsed/60.), fmod(elapsed,60.));


  sprintf(cmd, "cp image.txt image_%i.txt", (int)(Freq/1e6));
  //sprintf(cmd, "cp image.txt image_%ia.txt", (int)(Freq/1e6));
  system(cmd);
  sprintf(cmd, "cp traj.txt traj_%i.txt", (int)(Freq/1e6));
  system(cmd);
  
  return 0;
}
