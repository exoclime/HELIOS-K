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
__host__ void Init(Molecule &m){
	if(m.id == 1){//H2O
		m.NL = 224515;
		m.nISO = 6;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){11,  161,	.997317E+00,    1.7464E+02,    1,     18.010565};
		m.ISO[1] = (Isotopologue){12,  181,	1.99983E-03,    1.7511E+02,    1,     20.014811};
		m.ISO[2] = (Isotopologue){13,  171,	3.71884E-04,    1.0479E+03,    6,     19.014780};
		m.ISO[3] = (Isotopologue){14,  162,	3.10693E-04,    8.5901E+02,    6,     19.016740};
		m.ISO[4] = (Isotopologue){15,  182,	6.23003E-07,    8.7519E+02,    6,     21.020985};
		m.ISO[5] = (Isotopologue){16,  172,	1.15853E-07,    5.2204E+03,   36,     20.020956};
		sprintf(m.dataFilename, "%s", "01_hit12.par");
	}
	if(m.id == 2){//CO2
		m.NL = 471847;
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
		sprintf(m.dataFilename, "%s", "02_hit12.par");
	}
	if(m.id == 5){//CO
		m.NL = 4606;
		m.nISO = 6;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){51,  26,	.986544E+00,    1.0712E+02,    1,     27.994915};
		m.ISO[1] = (Isotopologue){52,  36,	1.10836E-02,    2.2408E+02,    2,     28.998270};
		m.ISO[2] = (Isotopologue){53,  28,	1.97822E-03,    1.1247E+02,    1,     29.999161};
		m.ISO[3] = (Isotopologue){54,  27,	3.67867E-04,    6.5934E+02,    6,     28.999130};
		m.ISO[4] = (Isotopologue){55,  38,	2.22250E-05,    2.3582E+02,    2,     31.002516};
		m.ISO[5] = (Isotopologue){56,  37,	4.13292E-06,    1.3809E+03,   12,     30.002485};
		sprintf(m.dataFilename, "%s", "05_hit12.par");
	}
	if(m.id == 6){//CH4
		m.NL = 468013;
		m.nISO = 4;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){61,  211,	.988274E+00,    5.9052E+02,    1,     16.031300};
		m.ISO[1] = (Isotopologue){62,  311,	1.11031E-02,    1.1808E+03,    2,     17.034655};
		m.ISO[2] = (Isotopologue){63,  212,	6.15751E-04,    4.7954E+03,    3,     17.037475};
		m.ISO[3] = (Isotopologue){64,  312,	6.91785E-06,    9.5990E+03,    6,     18.040830};
		sprintf(m.dataFilename, "%s", "06_hit12.par");
	}
}