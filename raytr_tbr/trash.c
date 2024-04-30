    if (Rho_I[iRay] >= RhoCr) {
      if (Flags_I[iRay] == FALSE) {
	// Tresspass happend first time for this ray:
	if (iRay == TrcRay) printf("++TRESPASS happend first time\n");
	//Pos_DI[0][iRay] = Dir_DI[0][iRay]*HalfDS; 
	//Pos_DI[1][iRay] = Dir_DI[1][iRay]*HalfDS;
	//Pos_DI[2][iRay] = Dir_DI[2][iRay]*HalfDS;

	Pos_DI[0][iRay] = PosPr_DI[iRay];
	Pos_DI[1][iRay] = PosPr_DI[nRay+iRay];
	Pos_DI[2][iRay] = PosPr_DI[n2Ry+iRay];
	
	Dir_DI[0][iRay] = DirPr_DI[iRay];
	Dir_DI[1][iRay] = DirPr_DI[nRay+iRay];
	Dir_DI[2][iRay] = DirPr_DI[n2Ry+iRay];

	Flags_I[iRay] = TRUE;   // Mark "bad ray"
      }
