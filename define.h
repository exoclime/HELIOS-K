#define VERSION 1.01

#define T0 296.0 		//Temperature in K
#define kB 1.3806489e-16 	//Boltzmann constant in erg/K
#define h 6.62606957e-27	//Planck costant in erg s
#define c 2.99792458e10 	//Speed of light cm/s
#define NA 6.0221412927e23	//Avogadro Constant  1/mol

#define TOL 1.43e-17		//Tolerance in the Voigt function 3.58e-9 2.48e-12 1.43e-17 3.25e-27 1.69e-33
#define NCheb 12		//Number of Chebychev coefficients in the q.dat file
#define nthmax 32768   		//Maximum number of threads in 2.0 Cards
#define nlmax 32768		//Maximum number of lines per kernel launch, to prevent from time out on Desktop machines
#define qALPHA_L 0.5		//q value in the Lorentz half width q = Pself / P

#define PROFILE	1	//1 = Voigt, 2 = Lorentz, 3 = Gauss

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
	int NL;			//Number of Lines
	Isotopologue *ISO;
	char dataFilename[160];
};

struct Partition{
	int n;
	int *id;
	double *Q;
};

struct Param{
	double T;
	double P;
	int nMolecule;
	double numin;
	double numax;
	double dnu;
	double dev;
};
