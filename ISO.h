//**************************************************
// This function initializes the Properties for the Isotopologues
// The values are taken from HITRAN molparam.txt
// The ids must be adapted in order to be consistent with the line list files in HITRAN
// The parameters are:
// m.id: Id of molecule
// m.NL: Number of Lines in the Molecule file
// m.nISO: Number of Isotopologues in the Molecule file
// The last line specifies the name of the Molecule file.

//Author: Simon Grimm
//November 2014
// ****************************************************
__host__ void Init(Molecule &m, Param param){
	if(m.id == 1){//H2O
		m.nISO = 6;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){11,  161,	.997317E+00,    1.7464E+02,    1,     18.010565};
		m.ISO[1] = (Isotopologue){12,  181,	1.99983E-03,    1.7511E+02,    1,     20.014811};
		m.ISO[2] = (Isotopologue){13,  171,	3.71884E-04,    1.0479E+03,    6,     19.014780};
		m.ISO[3] = (Isotopologue){14,  162,	3.10693E-04,    8.5901E+02,    6,     19.016740};
		m.ISO[4] = (Isotopologue){15,  182,	6.23003E-07,    8.7519E+02,    6,     21.020985};
		m.ISO[5] = (Isotopologue){16,  172,	1.15853E-07,    5.2204E+03,   36,     20.020956};
		m.nFiles = 1;		//number of data files
		//HITRAN2012
		m.NL[0] = 224515;	//number of lines
		m.NLmax = 224515;	//same as the number of lines
		sprintf(m.dataFilename[0], "%s%s", param.path, "01_hit12.par");
		//HITRAN2008
		//m.NL[0] = 69201;	//number of lines
		//m.NLmax = 69201;	//same as the number of lines
		//sprintf(m.dataFilename[0], "%s%s", param.path, "01_hit08.par");
		//Test
		//m.NL[0] = 64123;	//number of lines
		//m.NLmax = 64123;	//same as the number of lines
		//sprintf(m.dataFilename[0], "%s%s", param.path, "01_H2O.bin");
		
		
		if(param.useHITEMP == 1){
			m.nFiles = 34;			//number of data files
			m.NL[ 0] = 2048387;		//number of lines per data file
			m.NL[ 1] = 3470074;
			m.NL[ 2] = 3625575;
			m.NL[ 3] = 3500494;
			m.NL[ 4] = 4972481;
			m.NL[ 5] = 3561630;
			m.NL[ 6] = 3795745;
			m.NL[ 7] = 3907835;
			m.NL[ 8] = 3969887;
			m.NL[ 9] = 3820310;
			m.NL[10] = 4960367;
			m.NL[11] = 4163192;
			m.NL[12] = 4363563;
			m.NL[13] = 4336535;
			m.NL[14] = 3758715;
			m.NL[15] = 3283806;
			m.NL[16] = 2971465;
			m.NL[17] = 2992331;
			m.NL[18] = 3076556;
			m.NL[19] = 2708657;
			m.NL[20] = 2299459;
			m.NL[21] = 4583911;
			m.NL[22] = 2421691;
			m.NL[23] = 2834046;
			m.NL[24] = 4079576;
			m.NL[25] = 3784869;
			m.NL[26] = 3452571;
			m.NL[27] = 2680328;
			m.NL[28] = 3154700;
			m.NL[29] = 2671165;
			m.NL[30] = 2220333;
			m.NL[31] = 1735773;
			m.NL[32] = 3573701;
			m.NL[33] = 1461436;

			m.NLmax = 4972481;		//The naximum of number of lines per file

			sprintf(m.dataFilename[ 0], "%s%s", param.path, "01_00000-00050_HITEMP2010.par");
			sprintf(m.dataFilename[ 1], "%s%s", param.path, "01_00050-00150_HITEMP2010.par");
			sprintf(m.dataFilename[ 2], "%s%s", param.path, "01_00150-00250_HITEMP2010.par");
			sprintf(m.dataFilename[ 3], "%s%s", param.path, "01_00250-00350_HITEMP2010.par");
			sprintf(m.dataFilename[ 4], "%s%s", param.path, "01_00350-00500_HITEMP2010.par");
			sprintf(m.dataFilename[ 5], "%s%s", param.path, "01_00500-00600_HITEMP2010.par");
			sprintf(m.dataFilename[ 6], "%s%s", param.path, "01_00600-00700_HITEMP2010.par");
			sprintf(m.dataFilename[ 7], "%s%s", param.path, "01_00700-00800_HITEMP2010.par");
			sprintf(m.dataFilename[ 8], "%s%s", param.path, "01_00800-00900_HITEMP2010.par");
			sprintf(m.dataFilename[ 9], "%s%s", param.path, "01_00900-01000_HITEMP2010.par");
			sprintf(m.dataFilename[10], "%s%s", param.path, "01_01000-01150_HITEMP2010.par");
			sprintf(m.dataFilename[11], "%s%s", param.path, "01_01150-01300_HITEMP2010.par");
			sprintf(m.dataFilename[12], "%s%s", param.path, "01_01300-01500_HITEMP2010.par");
			sprintf(m.dataFilename[13], "%s%s", param.path, "01_01500-01750_HITEMP2010.par");
			sprintf(m.dataFilename[14], "%s%s", param.path, "01_01750-02000_HITEMP2010.par");
			sprintf(m.dataFilename[15], "%s%s", param.path, "01_02000-02250_HITEMP2010.par");
			sprintf(m.dataFilename[16], "%s%s", param.path, "01_02250-02500_HITEMP2010.par");
			sprintf(m.dataFilename[17], "%s%s", param.path, "01_02500-02750_HITEMP2010.par");
			sprintf(m.dataFilename[18], "%s%s", param.path, "01_02750-03000_HITEMP2010.par");
			sprintf(m.dataFilename[19], "%s%s", param.path, "01_03000-03250_HITEMP2010.par");
			sprintf(m.dataFilename[20], "%s%s", param.path, "01_03250-03500_HITEMP2010.par");
			sprintf(m.dataFilename[21], "%s%s", param.path, "01_03500-04150_HITEMP2010.par");
			sprintf(m.dataFilename[22], "%s%s", param.path, "01_04150-04500_HITEMP2010.par");
			sprintf(m.dataFilename[23], "%s%s", param.path, "01_04500-05000_HITEMP2010.par");
			sprintf(m.dataFilename[24], "%s%s", param.path, "01_05000-05500_HITEMP2010.par");
			sprintf(m.dataFilename[25], "%s%s", param.path, "01_05500-06000_HITEMP2010.par");
			sprintf(m.dataFilename[26], "%s%s", param.path, "01_06000-06500_HITEMP2010.par");
			sprintf(m.dataFilename[27], "%s%s", param.path, "01_06500-07000_HITEMP2010.par");
			sprintf(m.dataFilename[28], "%s%s", param.path, "01_07000-07500_HITEMP2010.par");
			sprintf(m.dataFilename[29], "%s%s", param.path, "01_07500-08000_HITEMP2010.par");
			sprintf(m.dataFilename[30], "%s%s", param.path, "01_08000-08500_HITEMP2010.par");
			sprintf(m.dataFilename[31], "%s%s", param.path, "01_08500-09000_HITEMP2010.par");
			sprintf(m.dataFilename[32], "%s%s", param.path, "01_09000-11000_HITEMP2010.par");
			sprintf(m.dataFilename[33], "%s%s", param.path, "01_11000-30000_HITEMP2010.par");
		}
	}
	if(m.id == 2){//CO2
		m.nFiles = 1;		//number of data files
		m.NL[0] = 471847;
		m.NLmax = 471847;
		m.nISO = 10;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){21,  626,	.984204E+00,    2.8694E+02,    1,     43.989830};
		m.ISO[1] = (Isotopologue){22,  636,	1.10574E-02,    5.7841E+02,    2,     44.993185};
		m.ISO[2] = (Isotopologue){23,  628,	3.94707E-03,    6.0948E+02,    1,     45.994076};
		m.ISO[3] = (Isotopologue){24,  627,	7.33989E-04,    3.5527E+03,    6,     44.994045};
		m.ISO[4] = (Isotopologue){25,  638,	4.43446E-05,    1.2291E+03,    2,     46.997431};
		m.ISO[5] = (Isotopologue){26,  637,	8.24623E-06,    7.1629E+03,   12,     45.997400};
		m.ISO[6] = (Isotopologue){27,  828,	3.95734E-06,    3.2421E+02,    1,     47.998322};
		m.ISO[7] = (Isotopologue){28,  728,	1.47180E-06,    3.7764E+03,    6,     46.998291};
		m.ISO[8] = (Isotopologue){29,  727,	1.36847E-07,    1.1002E+04,    1,     45.998262};
		m.ISO[9] = (Isotopologue){20,  838,	4.44600E-08,    6.5350E+02,    2,     49.001675};
		//m.ISO[0] = (Isotopologue){31,  837	1.65354E-08,    7.6152E+03,   12,     48.001646};
		sprintf(m.dataFilename[0], "%s%s", param.path, "02_hit12.par");

		if(param.useHITEMP == 1){
			m.nFiles = 20;			//number of data files
			m.NL[ 0] = 259836;		//number of lines per data file
			m.NL[ 1] = 730397;
			m.NL[ 2] = 1140808;
			m.NL[ 3] = 943281;
			m.NL[ 4] = 378616;
			m.NL[ 5] = 456531;
			m.NL[ 6] = 576241;
			m.NL[ 7] = 1103020;
			m.NL[ 8] = 507509;
			m.NL[ 9] = 212373;
			m.NL[10] = 214133;
			m.NL[11] = 939432;
			m.NL[12] = 1046856;
			m.NL[13] = 213769;
			m.NL[14] = 138258;
			m.NL[15] = 866153;
			m.NL[16] = 450170;
			m.NL[17] = 177444;
			m.NL[18] = 344265;
			m.NL[19] = 494516;

			m.NLmax = 1140808;              //The naximum of number of lines per file

			sprintf(m.dataFilename[ 0], "%s%s", param.path, "02_00000-00500_HITEMP2010.par");
			sprintf(m.dataFilename[ 1], "%s%s", param.path, "02_00500-00625_HITEMP2010.par");
			sprintf(m.dataFilename[ 2], "%s%s", param.path, "02_00625-00750_HITEMP2010.par");
			sprintf(m.dataFilename[ 3], "%s%s", param.path, "02_00750-01000_HITEMP2010.par");
			sprintf(m.dataFilename[ 4], "%s%s", param.path, "02_01000-01500_HITEMP2010.par");
			sprintf(m.dataFilename[ 5], "%s%s", param.path, "02_01500-02000_HITEMP2010.par");
			sprintf(m.dataFilename[ 6], "%s%s", param.path, "02_02000-02125_HITEMP2010.par");
			sprintf(m.dataFilename[ 7], "%s%s", param.path, "02_02125-02250_HITEMP2010.par");
			sprintf(m.dataFilename[ 8], "%s%s", param.path, "02_02250-02500_HITEMP2010.par");
			sprintf(m.dataFilename[ 9], "%s%s", param.path, "02_02500-03000_HITEMP2010.par");
			sprintf(m.dataFilename[10], "%s%s", param.path, "02_03000-03250_HITEMP2010.par");
			sprintf(m.dataFilename[11], "%s%s", param.path, "02_03250-03500_HITEMP2010.par");
			sprintf(m.dataFilename[12], "%s%s", param.path, "02_03500-03750_HITEMP2010.par");
			sprintf(m.dataFilename[13], "%s%s", param.path, "02_03750-04000_HITEMP2010.par");
			sprintf(m.dataFilename[14], "%s%s", param.path, "02_04000-04500_HITEMP2010.par");
			sprintf(m.dataFilename[15], "%s%s", param.path, "02_04500-05000_HITEMP2010.par");
			sprintf(m.dataFilename[16], "%s%s", param.path, "02_05000-05500_HITEMP2010.par");
			sprintf(m.dataFilename[17], "%s%s", param.path, "02_05500-06000_HITEMP2010.par");
			sprintf(m.dataFilename[18], "%s%s", param.path, "02_06000-06500_HITEMP2010.par");
			sprintf(m.dataFilename[19], "%s%s", param.path, "02_06500-12785_HITEMP2010.par");
		}

	}
	if(m.id == 5){//CO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 4606;
		m.NLmax = 4606;
		m.nISO = 6;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){51,  26,	.986544E+00,    1.0712E+02,    1,     27.994915};
		m.ISO[1] = (Isotopologue){52,  36,	1.10836E-02,    2.2408E+02,    2,     28.998270};
		m.ISO[2] = (Isotopologue){53,  28,	1.97822E-03,    1.1247E+02,    1,     29.999161};
		m.ISO[3] = (Isotopologue){54,  27,	3.67867E-04,    6.5934E+02,    6,     28.999130};
		m.ISO[4] = (Isotopologue){55,  38,	2.22250E-05,    2.3582E+02,    2,     31.002516};
		m.ISO[5] = (Isotopologue){56,  37,	4.13292E-06,    1.3809E+03,   12,     30.002485};
		sprintf(m.dataFilename[0], "%s%s", param.path, "05_hit12.par");

		if(param.useHITEMP == 1){
			m.nFiles = 1;			//number of data files
			m.NL[ 0] = 113631;		//number of lines per data file

			m.NLmax = 113631;              //The naximum of number of lines per file

			sprintf(m.dataFilename[0], "%s%s", param.path, "05_HITEMP2010.par");
		}
	}
	if(m.id == 6){//CH4
		m.nFiles = 1;		//number of data files
		m.NL[0] = 468013;
		m.NLmax = 468013;
		m.nISO = 4;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){61,  211,	.988274E+00,    5.9052E+02,    1,     16.031300};
		m.ISO[1] = (Isotopologue){62,  311,	1.11031E-02,    1.1808E+03,    2,     17.034655};
		m.ISO[2] = (Isotopologue){63,  212,	6.15751E-04,    4.7954E+03,    3,     17.037475};
		m.ISO[3] = (Isotopologue){64,  312,	6.91785E-06,    9.5990E+03,    6,     18.040830};
		sprintf(m.dataFilename[0], "%s%s", param.path, "06_hit12.par");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
	}
}


__host__ int InitCia(Molecule &m, ciaSystem &cia, Param param){
	m.NLmax = 0;
	cia.Nsets = 0;
	cia.mass1 = 1.0;
	if(strcmp(param.ciaSystem, "H2-H2") == 0){
		cia.Nsets = 113;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H2_2011.cia");
		cia.mass1 = 2.0 * 1.00794; //mass of H2 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-H2_eq") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H2_eq_2011.cia");
		cia.mass1 = 2.0 * 1.00794; //mass of H2 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-H2_norm") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H2_norm_2011.cia");
		cia.mass1 = 2.0 * 1.00794; //mass of H2 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-He") == 0){
		cia.Nsets = 339;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-He_2011.cia");
		cia.mass1 = 2.0 * 4.002602; //mass of He in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-He_eq") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-He_eq_2011.cia");
		cia.mass1 = 2.0 * 4.002602; //mass of He in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-He_norm") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-He_norm_2011.cia");
		cia.mass1 = 2.0 * 4.002602; //mass of He in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-CH4_eq") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-CH4_eq_2011.cia");
		cia.mass1 = 16.04246; //mass of CH4 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-CH4_norm") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-CH4_norm_2011.cia");
		cia.mass1 = 16.04246; //mass of CH4 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-H") == 0){
		cia.Nsets = 4;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H_2011.cia");
		cia.mass1 = 1.00794; //mass of H in g / mol
	}
	else{
		printf("Error: cia System not found %s\n", param.ciaSystem);
		return 0;
	}
	return 1;
}

