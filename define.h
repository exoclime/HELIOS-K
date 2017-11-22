#define VERSION 1.60

#define def_T0 296.0 		//Temperature in K
#define def_kB 1.3806489e-16 	//Boltzmann constant in erg/K
#define def_h 6.62606957e-27	//Planck costant in erg s
#define def_c 2.99792458e10 	//Speed of light cm/s
#define def_NA 6.0221412927e23	//Avogadro Constant  1/mol
#define def_amagat 2.6867774e19 // molecules cm^-3

#define M_PIf 3.14159265358979323846f

#define TOL 1.43e-17		//Tolerance in the Voigt function 3.58e-9 2.48e-12 1.43e-17 3.25e-27 1.69e-33
#define TOLF 2.48e-12f		//Tolerance in the Voigt function
#define NCheb 12		//Number of Chebychev coefficients in the q.dat file
#define nthmax 32768   		//Maximum number of threads in 2.0 Cards
#define nlmax 32768		//Maximum number of lines per kernel launch, to prevent from time out on Desktop machines
#define maxbins 1000		//Maximum number of bins
#define qALPHA_L 0.5		//q value in the Lorentz half width q = Pself / P
//#define maxlines 50000000	//maximum number of lines stored on the GPU, Should not be less than maximum in HITEMP lines
#define maxlines 8000000	//maximum number of lines stored on the GPU, Should not be less than maximum in HITEMP lines

#define PROFILE	1		//1 = Voigt, 2 = Lorentz, 3 = Gauss
#define NmaxSample 100		//Maximum Number of resample coefficients for K(y)
#define NXLOW 100000		//Linef versus Line2f

struct Isotopologue{
	int id;			//id in HITRAN notation
	int AFGL;		//id in AFGL notation
	double Ab;		//Abundance
	double Q;
	int g;
	double m;		//mass
};

struct Molecule{
	int id;			//Molecule number in HITRAN notation
	int nISO;		//Number of Isotopologues
	int NL[122];		//Number of Lines per file
	int NLmax;
	Isotopologue *ISO;
	char dataFilename[122][160];
	int fileLimit[123];
	int nFiles;
	int nStates;		//Number of states in EXOMOL linelist
	int ntcol;
	double meanMass;
	double defaultL;	//default value for Lorentz half width for EXOMOL
	double defaultn;		//default value for temperature exponent for EXOMOL
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
	int nMolecule;
	char ciaSystem[160];
	int useCia;
	char path[300];
	char pathK[300];
	int useHITEMP;
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
	int RLOW;
};

struct Line{
	double *nu_h, *nu_d;            //Wavenumber
	double *S_h, *S_d;        //Intensity
	float *Sf_d;
	double *S1_h, *S1_d;     //modified Intensity
	float *S1f_d;
	double *A_h, *A_d;              //Einstein A coefficient
	double *delta_h, *delta_d;      //pressure induced line shift
	double *EL_h, *EL_d;  		//Energy of lower state
	double *vy_h, *vy_d; 	//Lorentz Halfwidth / Doppler Halfwidth
	float *vyf_d;
	float *va_h, *va_d;	 	//(numin - nu) * ialphaD
	float *vb_h, *vb_d;	 	//dnu * ialphaD
	float *vcut2_h, *vcut2_d; 	//(cut * ialphaD)^2
	double *ialphaD_h, *ialphaD_d;  //Doppler Halfwidth
	double *n_h, *n_d;              //temperature dependent exponent
	double *Q_h, *Q_d;              //partition function
	int *ID_h, *ID_d;		//line id used for sorting

	
};
