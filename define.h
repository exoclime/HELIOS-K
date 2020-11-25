
#ifndef M_PI
#define _USE_MATH_DEFINES  //for Windows
#endif

#ifndef DEFINE_H
#define DEFINE_H


//Build Data

#ifndef GIT_DESCRIBE
#define GIT_DESCRIBE "Undefined"
#endif

#ifndef BUILD_DATE
#define BUILD_DATE "Undefined"
#endif

#ifndef BUILD_SYSTEM
#define BUILD_SYSTEM "Undefined"
#endif

#ifndef BUILD_PATH
#define BUILD_PATH "Undefined"
#endif

#ifndef BUILD_SM
#define BUILD_SM "Undefined"
#endif


#define VERSION 2.06


#define def_T0 296.0 		//Reference Temperature in K
#define def_PObar 0.986923	//Referecne Pressure 1 bar in atm for ExoMol
#define def_POatm 1.0		//Referecne Pressure 1 atm in atm for Hitran
#define def_kB 1.3806489e-16 	//Boltzmann constant in erg/K = cm^2 g / (s^2 K)
#define def_h 6.62606957e-27	//Planck costant in erg s
#define def_c 2.99792458e10 	//Speed of light cm/s
#define def_NA 6.0221412927e23	//Avogadro Constant  1/mol
#define def_amagat 2.6867774e19 // molecules cm^-3

#define M_PIf 3.14159265358979323846f

#define def_TOL 1.43e-17		//Tolerance in the Voigt function 3.58e-9 2.48e-12 1.43e-17 3.25e-27 1.69e-33
#define def_TOLF 2.48e-12f		//Tolerance in the Voigt function
#define def_nthmax 1048576   		//Maximum number of threads for line and sort kernels
#define def_nlmax 32768			//Maximum number of lines per kernel launch, to prevent from time out on Desktop machines
#define def_maxlines 1048576ll		//maximum number of lines stored on the GPU, Should not be less than maximum in HITEMP lines, must be long long int type
#define def_maxfiles 500		//maximum number of files per molecule

#define def_NmaxSample 100		//Maximum Number of resample coefficients for K(y)


#define def_nlA 1024			//split lines into blocks of this size, A Line kernels
#define def_nlB 1024			//split lines into blocks of this size, A Line kernels
#define def_nlC 1024			//split lines into blocks of this size, A Line kernels

#define def_rBs 256			//number of streams in readBuffer copy

#define def_KSn 32


//default values of parameters
#define def_qALPHA_L 0.5	//q value in the Lorentz half width q = Pself / P
#define def_gammaF 1.0		//scaling factor for Lorentzian halfwidth
#define def_PROFILE 1		//1 = Voigt, 2 = Lorentz, 3 = Gauss, 4 = cross section
#define def_doTuning 1		//use the self-tuning routines
#define def_removePlinth 0	//use to remove the plinth from cutted line wings


struct Isotopologue{
	int id;			//id in HITRAN notation
	double Ab;		//Abundance
	double Q;		//partition function reference value
	int g;			//degeneracy
	double m;		//mass
	char cid[4];		//id in HITRAN notation as string
};

struct Molecule{
	int id;			//Molecule number in HITRAN notation
	int nISO;		//Number of Isotopologues
	long long int NL[def_maxfiles];		//Number of Lines per file
	long long int NLmax;
	Isotopologue *ISO;
	char mName[160];	//name of states and trans files
	char dataFilename[def_maxfiles][160];
	int fileLimit[def_maxfiles + 1];
	int nFiles;
	int nStates;		//Number of states in EXOMOL linelist
	int ntcol;		//Number of columns in transition files
	int npfcol;		//number of columns in partition function file
	double defaultL;	//default value for Lorentz half width for EXOMOL
	double defaultn;	//default value for temperature exponent for EXOMOL
	int version;
};

struct Partition{
	int n;
	int *id;
	double *Q;
};

struct ciaSystem{
	int Nsets;
	char dataFilename[160];
	double mass1;
	double mass2;
};


struct Param{
	char name[160];
	double T;
	double P;
	char PFilename[160];
	int nP;
	int usePFile; 
	char ciaSystem[160];
	int useCia;
	char path[300];
	char pathK[300];
	char mParamFilename[400];
	char SpeciesFilename[160];
	int nSpecies;
	int useSpeciesFile; 
	int dataBase;
	double numin;
	double numax;
	double dnu;
	int Nxb;
	int cutMode;
	double cut;
	int dev;
	int nbins;
	int nedges;
	char bins[160];
	char edges[160];
	int doResampling;
	int nC;
	int doTransmission;
	int nTr;
	double dTr;
	int doStoreFullK;
	int doStoreK;
	double kmin;
	double qalphaL;
	int doMean;
	int useIndividualBins;
	int useOutputEdges;
	int units;
	int useIndividualX;
	int replaceFiles;
	int profile;
	double gammaF;
	int doTuning;
	int removePlinth;
	char subLorentzianFilename[160];
	int useSubLorentzian;
};

//Sub-Lorentzian coefficients
float sLchi_h[17];

struct Line{
	double *nu_h, *nu_d;		//Wavenumber
	double *S_h, *S_d;		//Intensity
	float *Sf_d;
	double *S1_d;			//modified Intensity
	float *S1f_d;
	double *A_h, *A_d;		//Einstein A coefficient
	double *delta_h, *delta_d;	//pressure induced line shift
	double *EL_h, *EL_d;		//Energy of lower state
	double *vy_h, *vy_d;		//Lorentz Halfwidth / Doppler Halfwidth
	float *vyf_d;
	float *va_d;			//(numin - nu) * ialphaD
	float *vb_d;			//dnu * ialphaD
	float *vcut2_d;			//(cut * ialphaD)^2
	double *ialphaD_h, *ialphaD_d;	//Doppler Halfwidth
	double *plinth_d;
	double *n_h, *n_d;		//temperature dependent exponent
	double *Sort_d;			//helper array used to sort the other arrays
	int *ID_d;			//line id used for sorting
	double *nuLimitsA0_d;	//limits for Line blocks, min
	double *nuLimitsA1_d;	//limits for Line blocks, max
	double *nuLimitsAL0_d;	//limits for Line blocks, min
	double *nuLimitsAL1_d;	//limits for Line blocks, max
	double *nuLimitsAR0_d;	//limits for Line blocks, min
	double *nuLimitsAR1_d;	//limits for Line blocks, max
	double *nuLimitsB0_d;	//limits for Line blocks, min
	double *nuLimitsB1_d;	//limits for Line blocks, max
	double *nuLimitsC0_d;	//limits for Line blocks, min
	double *nuLimitsC1_d;	//limits for Line blocks, max

	long long int *iiLimitsA0_h, *iiLimitsA0_d;	//limits for Line blocks, min
	long long int *iiLimitsA1_h, *iiLimitsA1_d;	//limits for Line blocks, max
	long long int *iiLimitsAL0_h, *iiLimitsAL0_d;	//limits for Line blocks, min
	long long int *iiLimitsAL1_h, *iiLimitsAL1_d;	//limits for Line blocks, max
	long long int *iiLimitsAR0_h, *iiLimitsAR0_d;	//limits for Line blocks, min
	long long int *iiLimitsAR1_h, *iiLimitsAR1_d;	//limits for Line blocks, max
	long long int *iiLimitsB0_h, *iiLimitsB0_d;	//limits for Line blocks, min
	long long int *iiLimitsB1_h, *iiLimitsB1_d;	//limits for Line blocks, max
	long long int *iiLimitsC0_h, *iiLimitsC0_d;	//limits for Line blocks, min
	long long int *iiLimitsC1_h, *iiLimitsC1_d;	//limits for Line blocks, max

	long long int *iiLimitsAT_m, *iiLimitsAT_d;	//limits for all blocks, mapped memory
	long long int *iiLimitsALT_m, *iiLimitsALT_d;	//limits for all blocks, mapped memory
	long long int *iiLimitsART_m, *iiLimitsART_d;	//limits for all blocks, mapped memory
	long long int *iiLimitsBT_m, *iiLimitsBT_d;	//limits for all blocks, mapped memory
	long long int *iiLimitsCT_m, *iiLimitsCT_d;	//limits for all blocks, mapped memory

};
#endif
