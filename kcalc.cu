#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "define.h"
#include "ISO.h"

// *********************************************
//This function calculates the Series Sigma1. Sigma2 and Sigma3 (Equations 27, 28, and 29) from Alg 916
//The parameter TOL sets a tolerance where to truncate the series
//It returns the values for sigma1 sigma2 and sigma3
//
// This implementation has still problems with loss of accuracy
//
//
//Author Simon Grimm, Adapted from Zaghloul & Ali, Algorithm 916
//November 2014
// **********************************************
__device__ void Sigma(double x, double y, double &s1, double &s2, double &s3, double a, double ex2, int id){

	s1 = 0.0;
	double sold1 = s1;
	s2 = 0.0;
	double sold2 = s2;
	s3 = 0.0;
	double sold3 = s3;

	double f, f3p, f3n;
	double an, an3p, an3n;

	double yy = y * y;

	int n0 = (int)(ceil(x / a)); //starting point for sigma3 series
	int n3p, n3n;

	int stop1 = 0;
	int stop2 = 0;
	int stop3 = 0;

	double ean2;

	double e2axn = exp(-2.0 * a * x);
	double e2axns = 1.0;

	double an0 = a * n0;
	double e3p = exp(-2.0 * a * (an0 - x - a));
	double e3n = exp(2.0 * a * (an0 - x));

	double e3ps = exp(-(an0 * an0 - 2.0 * a * an0 - 2.0 * an0 * x + x * x + a * a + 2.0 * a * x));
	double e3ns = exp(-(an0 * an0 - 2.0 * an0 * x + x * x));

	double st;

	for(int n = 1; n < 100; ++n){
		n3p = n0 + n - 1;
		n3n = n0 - n;
		an = a * n;
		an3p = a * n3p;
		an3n = a * n3n;

		ean2 = exp(-an * an);
		e2axns *= e2axn;	

		e3ps *= e3p;
		e3ns *= e3n;

		f = 1.0 / (an * an + yy);
		f3p = 1.0 / (an3p * an3p + yy);
		f3n = 1.0 / (an3n * an3n + yy);

		st = f * ean2 * ex2;

		s1 += st;
		s2 += st * e2axns;
		s3 += f3p * ean2 * e3ps;
		if(n3n >= 1) s3 += f3n * ean2 * e3ns;

		
		if(fabs(s1 - sold1) < TOL) stop1 = 1;		
		if(fabs(s2 - sold2) < TOL) stop2 = 1;		
		if(fabs(s3 - sold3) < TOL) stop3 = 1;		
		if(stop1 == 1 && stop2 ==1 && stop3 == 1) break;

		sold1 = s1;
		sold2 = s2;
		sold3 = s3;
//if(n >= 100-1) printf("Sigma Series did not converge\n");
	}
}
// *********************************************
//This function calculates the Series Sigma1. Sigma2 and Sigma3 (Equations 15, 16, and 17) from Alg 916
//The parameter TOL sets a tolerance where to truncate the series
//It returns the values for sigma1 sigma2 and sigma3

//Author Simon Grimm, Adapted from Zaghloul & Ali, Algorithm 916
//November 2014
// **********************************************
__device__ void Sigmab(double x, double y, double &s1, double &s2, double &s3, double a, double ex2, int id){

	s1 = 0.0;
	double sold1 = s1;
	s2 = 0.0;
	double sold2 = s2;
	s3 = 0.0;
	double sold3 = s3;

	double f, f3p, f3n;
	double an, an3p, an3n;

	double yy = y * y;

	if(x < 0.0) x = -x;

		int n0 = (int)(ceil(x / a)); //starting point for sigma3 series
		int n3p, n3n;

		int stop1 = 0;
		int stop2 = 0;
		int stop3 = 0;

		double e2axn = exp(-2.0 * a * x);

		for(int n = 1; n < 100; ++n){
		n3p = n0 + n - 1;
		n3n = n0 - n;
		an = a * n;
		an3p = a * n3p;
		an3n = a * n3n;

		f = 1.0 / (an * an + yy);
		f3p = 1.0 / (an3p * an3p + yy);
		f3n = 1.0 / (an3n * an3n + yy);

		s1 += f * exp(-(an * an + x * x));
		s2 += f * exp(-(an + x) * (an + x));
		s3 += f3p * exp(-(an3p - x) * (an3p - x));
		if(n3n >= 1) s3 += f3n * exp(-(an3n - x) * (an3n - x));

		if(fabs(s1 - sold1) < TOL) stop1 = 1;
		if(fabs(s2 - sold2) < TOL) stop2 = 1;
		if(fabs(s3 - sold3) < TOL) stop3 = 1;
		if(stop1 == 1 && stop2 ==1 && stop3 == 1) break;

		sold1 = s1;
		sold2 = s2;
		sold3 = s3;
//		if(n >= 100-1) printf("Sigma Series did not converge %d\n", id);
	}
}
// **************************************************
//This function calculates the Voigt profile V(x,y) as equation 13 from Zaghloul & Ali, Algorithm 916
//it calls the Sigma function
//The parameter TOL sets a tolerance where to truncate the series

//Author Simon Grimm, Adapted from Zaghloul & Ali, Algorithm 916
//November 2014
// **********************************************
__device__ double voigt_916(double x, double y, double a, int id){

	double s1, s2, s3;
	double ex2 = exp(-x * x);

	//Compute Sigma Series
	if(x != 0.0 && y != 0.0) Sigmab(x, y, s1, s2, s3, a, ex2, id);

	double xy = x * y;
	double a2ipi = 2.0 * a / M_PI;
	double cos2xy = cos(2.0 * xy);
	double sinxy = sin(xy);

	double t1;
	t1 = ex2 * erfcx(y) * cos2xy;
	t1 += a2ipi * x * sinxy * ex2 * sinxy / xy;
	t1 += a2ipi * y * (-cos2xy * s1 + 0.5 * (s2 + s3));
	
	if(x == 0) t1 = erfcx(y);
	if(y == 0) t1 = exp(-x * x);
	//if(x*x + y*y > 1.0e18) t1 = y / (sqrt(M_PI) * (x * x + y * y));
	
	return t1;

}
// **************************************************
//This kernel calculates the integrate line strength, the Lorentz and the Doppler halfwidths
//
//Author Simon Grimm
//November 2014
// **********************************************
__global__ void S_kernel(double *nu_d, double *S_d, double *A_d, double *EL_d, double *alphaL_d, double *alphaD_d, double *n_d, double *mass_d, double *gamma_d, double *delta_d, double *Q_d, int NL, double T, double P, int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){
		double c2 = c * c;
		double m = mass_d[id] / NA;			// mass in g

		double nu = nu_d[id] + delta_d[id] * P;		//read nu from alphaD
		nu_d[id] = nu;
		double S = S_d[id] / m;				//cm / g
		double EL = EL_d[id];  				//1/cm
		double Q = Q_d[id];				//Q0 / Q(T)
		double alphaL = alphaL_d[id];
		
		S_d[id] = S * Q * exp(-EL * h * c / (kB * T) + EL * h * c / (kB * T0)) * (1.0 - exp(-h * nu * c / (kB * T))) / (1.0 - exp(-h * nu * c / (kB * T0))); 
		alphaD_d[id] = 1.0 / (nu * sqrt(2.0 * kB * T / (m * c2)));	//inverse Doppler halfwith
		alphaL *= P * pow(T0 / T, gamma_d[id]);
		alphaL += A_d[id] / (4.0 * M_PI * c);				//1/cm
		alphaL_d[id] = alphaL;
	}
}

// **************************************************
//This kernel directly calls the Voigt function
//
//Author Simon Grimm
//November 2014
// **********************************************
__global__ void Voigt_line_kernel(double a, double dnu, double *K_d, double Nx, double xmax){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	if(id < Nx){
		double aTOL = M_PI * sqrt(-1.0 / log(TOL * 0.5));
		double x = fabs(-xmax + id * 2.0 * xmax / ((double)(Nx)));
		K_d[id] = voigt_916(x, a, aTOL, id);
	}
}

// **************************************************
//This kernel computes the line shape
// It uses patterns of shared memory the reduce global memory access
//
//Author Simon Grimm
//November 2014
// **********************************************
template <int NB, int nl>
__global__ void Line_kernel(double *nu_d, double *S_d, double *alphaL_d, double *alphaD_d, double *n_d, double *K_d, double dnu, double numin, int Nx, int NL, int ii, int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	__shared__ double nu_s[NB];
	__shared__ double S_s[NB];
	__shared__ double alphaL_s[NB];
	__shared__ double ialphaD_s[NB];
	__shared__ int xyFlag_s[2];

	double K = 0.0;

	double a = M_PI * sqrt(-1.0 / log(TOL * 0.5));
	double sqln2 = sqrt(log(2.0));

	double isqrtpi = 1.0 / sqrt(M_PI);
	double nu = numin + id * dnu;

	for(int i = 0; i < nl; i += NB){
		xyFlag_s[0] = 0;
		xyFlag_s[1] = 0;
		if(i + idx + ii < NL){
			nu_s[idx] = nu_d[i + idx + ii];
			S_s[idx] = S_d[i + idx + ii];
			alphaL_s[idx] = alphaL_d[i + idx + ii];
			ialphaD_s[idx] = alphaD_d[i + idx + ii];
		}
		else{
			nu_s[idx] = 0.0;
			S_s[idx] = 0.0;
			alphaL_s[idx] = 0.0;
			ialphaD_s[idx] = 0.0;
		}
		__syncthreads();

# if PROFILE == 1
		//Check smallest values for x and y
		for(int j = 0; j < NB; ++j){
			if(i + j < NL){
				double x = sqln2 * fabs((nu - nu_s[j]) * ialphaD_s[j]);
				double y = sqln2 * alphaL_s[j] * ialphaD_s[j];
				if(x * x + y * y < 1e6) xyFlag_s[1] = 1; 
				if(x * x + y * y < 100) xyFlag_s[0] = 1; 
			}
		}
		__syncthreads();
		if(xyFlag_s[0] == 1){
			for(int j = 0; j < NB; ++j){
				if(i + j < NL){
					double x = sqln2 * fabs((nu - nu_s[j]) * ialphaD_s[j]);
					double y = sqln2 * alphaL_s[j] * ialphaD_s[j];
					K += S_s[j] * voigt_916(x, y, a, id) * sqln2 * ialphaD_s[j] * isqrtpi;
//if(id < 32) printf("%d %.20g %.20g %.20g %.20g %.20g %.20g\n", id, nu, x, y, nu_s[j], 1.0 / ialphaD_s[j], alphaL_s[j]);
				}

			}
		}
		else if(xyFlag_s[1] == 1){
			for(int j = 0; j < NB; ++j){
				if(i + j < NL){
					//2nd order Gauss Hermite Quadrature
					double x = sqln2 * fabs((nu - nu_s[j]) * ialphaD_s[j]);
					double y = sqln2 *  alphaL_s[j] * ialphaD_s[j];
					double xxyy = x * x + y * y;
					double t = y / 3.0;
					double t1 = 2.0 * t / (M_PI * xxyy);
					double t2 = t * (xxyy + 1.5) / (M_PI * (xxyy + 1.5) * (xxyy + 1.5) - 4.0 * x * x * 1.5);
					K += S_s[j] * sqln2 * ialphaD_s[j] * (t1 + t2);
				}
			}
		}
		else{
			for(int j = 0; j < NB; ++j){
				if(i + j < NL){
					//1 order Gauss Hermite Quadrature
					double x = sqln2 * fabs((nu - nu_s[j]) * ialphaD_s[j]);
					double y = sqln2 *  alphaL_s[j] * ialphaD_s[j];
					K += S_s[j] * sqln2 * y * ialphaD_s[j] / (M_PI * (x * x + y * y));
				}
			}
		}
#endif
# if PROFILE == 2
		for(int j = 0; j < NB; ++j){
			if(i + j < NL){
				K += S_s[j] * alphaL_s[j] / (M_PI * ((nu - nu_s[j]) * (nu - nu_s[j]) + alphaL_s[j] * alphaL_s[j]));
			}	
		}
#endif
# if PROFILE == 3
		for(int j = 0; j < NB; ++j){
			if(i + j < NL){
				K += S_s[j] * ialphaD_s[j] * isqrtpi * exp(-x * x);
			}	
		}
#endif
		__syncthreads();
	}
	if(id < Nx) K_d[id] += K;
}

// ****************************************
// This function computes the Chebyshev polynomials as a function of T
// n must by greater than 1.
// *****************************************
__host__ void Chebyshev(double T, double *Cheb, int n){
	Cheb[0] = 1.0;
	Cheb[1] = T;
	for(int i = 2; i < n; ++i){
		Cheb[i] = 2.0 * T * Cheb[i - 1] - Cheb[i - 2];
	}
	
}
// ****************************************
// This function reads the file q.dat and computes for each Isotopologue
// the corresponding Partition function Q(T)

//Author: Simon Grimm
//November 2014
// *****************************************
__host__ int ChebCoeff(char *qFilename, Partition &part, double T){

	//Calculate Chebychev polynomial
	double Cheb[NCheb];
	Chebyshev(T, Cheb, NCheb);

	//Read Chebychev Coefficients from q file	
	FILE *qFile;
	//Check size of q.dat file
	qFile = fopen(qFilename, "r");
	if(qFile == NULL){
		printf("Error: q.dat file not found\n");
		return 0;
	}
	int j;
	for(j = 0; j < 100; ++j){
		int id;
		double coeff;
		int er = fscanf (qFile, "%d", &id);
		if (er <= 0) break;
		for(int i = 0; i < NCheb; ++i){
			fscanf (qFile, "%lf", &coeff);
		}
	}
	fclose(qFile);
	part.n = j;
	
	part.id = (int*)malloc(j * sizeof(int));
	part.Q = (double*)malloc(j * sizeof(double));
	
	qFile = fopen(qFilename, "r");
	for(j = 0; j < 100; ++j){
		int id;
		double coeff;
		double Q = 0.0;
		int er = fscanf (qFile, "%d", &id);
		if (er <= 0) break;
		for(int i = 0; i < NCheb; ++i){
			fscanf (qFile, "%lf", &coeff);
			Q += coeff * Cheb[i];
		}
		part.id[j] = id;
		part.Q[j] = Q;
	}
	fclose(qFile);
	return 1;
}

__host__ int read_parameters(Param &param, char *paramFilename, int argc, char*argv[]){
	//Read parameters from param.dat file
	FILE *paramFile;
	param.dev = 0;
	paramFile = fopen(paramFilename, "r");
		char skip[160];
		char skip2[160];
		//read T
		fgets(skip, 4, paramFile);
		fscanf (paramFile, "%lf", &param.T);
		fgets(skip2, 3, paramFile);
		//read P
		fgets(skip, 4, paramFile);
		fscanf (paramFile, "%lf", &param.P);
		fgets(skip2, 3, paramFile);
		//read Molecule
		fgets(skip, 11, paramFile);
		fscanf (paramFile, "%d", &param.nMolecule);
		fgets(skip2, 3, paramFile);
		//read numin
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%lf", &param.numin);
		fgets(skip2, 3, paramFile);
		//read numax
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%lf", &param.numax);
		fgets(skip2, 3, paramFile);
		//read dnu
		fgets(skip, 6, paramFile);
		fscanf (paramFile, "%lf", &param.dnu);
		fgets(skip2, 3, paramFile);

	fclose(paramFile);

	//Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-T") == 0){
			param.T = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-P") == 0){
			param.P = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-m") == 0){
			param.nMolecule = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-numin") == 0){
			param.numin = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-numax") == 0){
			param.numax = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dnu") == 0){
			param.dnu = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dev") == 0){
			param.dev = atoi(argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}

	}
	return 1;
}

int main(int argc, char*argv[]){

	cudaError_t error;
	int er;

	int devCount = 0;
	cudaGetDeviceCount(&devCount);

	if(devCount == 0){
		printf("Error: No valid cuda device!\n");
		return 0;
	}
	if(devCount == 1) printf("There is %d CUDA Device\n", devCount); 
	else printf("There are %d CUDA Devices\n", devCount); 


	FILE *dataFile;
	FILE *OutFile;
	FILE *Out2File;
	char qFilename[160];
	char paramFilename[160];
	char OutFilename[160];
	char Out2Filename[160];
	sprintf(qFilename, "%s", "q.dat");
	sprintf(paramFilename, "%s", "param.dat");
	sprintf(OutFilename, "%s", "Out.dat");
	sprintf(Out2Filename, "%s", "Out2.dat");

	//Read prameters
	Param param;
	er = read_parameters(param, paramFilename, argc, argv);
	if(er == 0){
		return 0;
	}
	if(param.dev > devCount || param.dev < 0){
		printf("Error: Devive Number is not allowed\n");
		return 0;
	}

	cudaSetDevice(param.dev);
	printf("Version: %g\n", VERSION);
	printf("T = %g\nP = %g\nMolecule = %d\nnumin = %g\nnumax = %g\ndnu = %g\n", param.T, param.P, param.nMolecule, param.numin, param.numax, param.dnu);
	printf("Profile = %d\n", PROFILE);
	printf("Using device %d\n", param.dev);

	int Nx = (int)((param.numax - param.numin) / param.dnu);

	Partition part;
	er = ChebCoeff(qFilename, part, param.T);
	if(er == 0){
		return 0;
	}

	Molecule m;
	m.NL = 0;
	m.id = param.nMolecule;	//1 = H2O, 2 = CO, 5 = CO, 6 = CH4
	m.nISO = 0;

	//Initialize the Isotopologue properties for ISO.h
	Init(m);

	if(m.NL == 0){
		printf("Molecule Id is not allowed\n");
		return 0;
	}
	
	
	timeval tt1;			//start time
	timeval tt2;			//end time
	long long times, timems;	//elapsed time in seconds and microseconds

	cudaDeviceSynchronize();
	gettimeofday(&tt1, NULL);
	times = 0.0;
	timems = 0.0;
	
	double *nu_h, *nu_d;		//Wavenumber
	double *S_h, *S_d;		//Intensity
	double *A_h, *A_d;		//Einstein A coefficient
	double *delta_h, *delta_d;	//line shift
	double *EL_h, *EL_d;		//lower state energy
	double *alphaL_h, *alphaL_d;	//Lorentz Halfwidth
	double *alphaD_h, *alphaD_d;	//Doppler Halfwidth
	double *n_h, *n_d;		//temperature dependent exponent
	double *gamma_h, *gamma_d;	//pressure shift coefficient
	double *mass_h, *mass_d;
	double *Q_h, *Q_d;		//partition function

	double *K_h, *K_d;
	
	nu_h = (double*)malloc(m.NL * sizeof(double));
	S_h = (double*)malloc(m.NL * sizeof(double));
	A_h = (double*)malloc(m.NL * sizeof(double));
	delta_h = (double*)malloc(m.NL * sizeof(double));
	EL_h = (double*)malloc(m.NL * sizeof(double));
	alphaL_h = (double*)malloc(m.NL * sizeof(double));
	alphaD_h = (double*)malloc(m.NL * sizeof(double));
	n_h = (double*)malloc(m.NL * sizeof(double));
	gamma_h = (double*)malloc(m.NL * sizeof(double));
	mass_h = (double*)malloc(m.NL * sizeof(double));
	Q_h = (double*)malloc(m.NL * sizeof(double));

	cudaMalloc((void **) &nu_d, m.NL * sizeof(double));
	cudaMalloc((void **) &S_d, m.NL * sizeof(double));
	cudaMalloc((void **) &A_d, m.NL * sizeof(double));
	cudaMalloc((void **) &delta_d, m.NL * sizeof(double));
	cudaMalloc((void **) &EL_d, m.NL * sizeof(double));
	cudaMalloc((void **) &alphaL_d, m.NL * sizeof(double));
	cudaMalloc((void **) &alphaD_d, m.NL * sizeof(double));
	cudaMalloc((void **) &n_d, m.NL * sizeof(double));
	cudaMalloc((void **) &gamma_d, m.NL * sizeof(double));
	cudaMalloc((void **) &mass_d, m.NL * sizeof(double));
	cudaMalloc((void **) &Q_d, m.NL * sizeof(double));


	K_h = (double*)malloc(Nx * sizeof(double));
	cudaMalloc((void **) &K_d, Nx * sizeof(double));


	for(int i = 0; i < Nx; ++i){
		K_h[i] = 0.0;
	}

	
	{
		dataFile  = fopen(m.dataFilename, "r");
		if(dataFile == NULL){
			printf("Error: line list file not found\n");
			return 0;
		}
		//read line list file		

		char c1[3];
		//char c2[2];
		char c3[13];
		char c4[11];
		char c5[11];
		char c6[6];
		char c7[6];
		char c8[11];
		char c9[5];
		char c10[9];
		char c11[16];
		char c12[16];
		char c13[16];
		char c14[16];
		char c15[7];
		char c16[13];
		char c17[2];
		char c18[8];
		char c19[8];
		
		char skip[5];

//int count[40];
//for(int cc = 0; cc < 40; ++cc){
//count[cc] = 0;	
//}
		for(int i = 0; i < m.NL; ++i){
		
			fgets(skip, 1, dataFile);
			//fgets(c1, 3, dataFile);
			//fgets(c2, 2, dataFile);
			fgets(c1, 4, dataFile);		//Use combined notation for Id (AFGL and molecule + abundance number
			fgets(c3, 13, dataFile);
			fgets(c4, 11, dataFile);
			fgets(c5, 11, dataFile);
			fgets(c6, 6, dataFile);
			fgets(c7, 6, dataFile);
			fgets(c8, 11, dataFile);
			fgets(c9, 5, dataFile);
			fgets(c10, 9, dataFile);
			fgets(c11, 16, dataFile);
			fgets(c12, 16, dataFile);
			fgets(c13, 16, dataFile);
			fgets(c14, 16, dataFile);
			fgets(c15, 7, dataFile);
			fgets(c16, 13, dataFile);
			fgets(c17, 2, dataFile);
			fgets(c18, 8, dataFile);
			fgets(c19, 8, dataFile);
			fgets(skip, 6, dataFile);
			
			nu_h[i] = strtod(c3, NULL);		
			S_h[i] = strtod(c4, NULL);		
			A_h[i] = strtod(c5, NULL);		
			delta_h[i] = strtod(c9, NULL);
			EL_h[i] = strtod(c8, NULL);		
			
			double gammaAir = strtod(c6, NULL);
			double gammaSelf = strtod(c7, NULL);
			alphaL_h[i] = (1.0 - qALPHA_L) * gammaAir + qALPHA_L * gammaSelf;
			n_h[i] = strtod(c9, NULL);
			alphaD_h[i] = n_h[i];
			
			gamma_h[i] = strtod(c10, NULL);
			int id= std::atoi(c1);
			int idAFGL;
//count[0] += 1;
			//Assign the Isotopologue properties
			for(int j = 0; j < m.nISO; ++j){
				if(id == m.ISO[j].id){
					mass_h[i] = m.ISO[j].m;
					Q_h[i] = m.ISO[j].Q;
					idAFGL = m.ISO[j].AFGL;
				}
			}
			double Q;
			//Assign the Partition function
			for(int j = 0; j < part.n; ++j){
				if(idAFGL == part.id[j]){
					Q = part.Q[j];
				}
			}
			Q_h[i] /= exp(Q);
			
//if(i < 10) printf("%d %d %d %g %g\n", i, id, idAFGL, exp(Q), Q_h[i]);
		}
//for(int cc = 0; cc < 40; ++cc){
//printf("%d %d\n", cc, count[cc]);	
//}
		fclose(dataFile);
	}
	
	
//for(int i = 0; i < m.NL; ++i){
//	printf("%d %g %g\n", i, nu_h[i], S_h[i]);

//}

	cudaMemcpy(nu_d, nu_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(S_d, S_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(A_d, A_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(delta_d, delta_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(EL_d, EL_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(alphaL_d, alphaL_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(alphaD_d, alphaD_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(n_d, n_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(gamma_d, gamma_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(mass_d, mass_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Q_d, Q_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);

	cudaMemcpy(K_d, K_h, Nx * sizeof(double), cudaMemcpyHostToDevice);


	for(int k = 0; k < m.NL; k += nthmax){
		int Nk = min(nthmax, m.NL);
		S_kernel <<< (Nk + 127) / 128, 128 >>> (nu_d, S_d, A_d, EL_d, alphaL_d, alphaD_d, n_d, mass_d, gamma_d, delta_d, Q_d, m.NL, param.T, param.P, k);
	}	

	cudaDeviceSynchronize();
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	printf("Time before Line_kernel: %g seconds\n", times + timems/1000000.0);
/*
	cudaMemcpy(nu_h, nu_d, m.NL * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(S_h, S_d, m.NL * sizeof(double), cudaMemcpyDeviceToHost);
	for(int i = 0; i < m.NL; ++i){
		printf("%d %.20g %.20g\n", i, nu_h[i], S_h[i]);

	}
*/
	cudaDeviceSynchronize();
	gettimeofday(&tt1, NULL);
	times = 0.0;
	timems = 0.0;

	for(int k = 0; k < Nx; k += nthmax){
		int Nk = min(nthmax, Nx);
		printf("Reached k =  %d, Total = %d\n", k, Nx);
		for(int i = 0; i < m.NL; i += nlmax){
			//This loop reduces the running time of the kernel to a few seconds
			//A longer running time of a single kernel can cause a time out
			Line_kernel <32, nlmax> <<< (Nk + 31) / 32, 32 >>> (nu_d, S_d, alphaL_d, alphaD_d, n_d, K_d, param.dnu, param.numin, Nx, m.NL, i, k);
		}
	}
	cudaDeviceSynchronize();
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	printf("Time for Line_kernel:    %g seconds\n", times + timems/1000000.0);

	gettimeofday(&tt1, NULL);
	times = 0.0;
	timems = 0.0;

	cudaMemcpy(K_h, K_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);

	OutFile = fopen(OutFilename, "w");
	for(int j = 0; j < Nx; ++j){
		double x = param.numin + j * param.dnu;
		fprintf(OutFile, "%.20g %.20g\n", x, K_h[j]);
	}
	fclose(OutFile);


	thrust::device_ptr<double> K_dt = thrust::device_pointer_cast(K_d);
	thrust::sort(K_dt, K_dt + Nx);

	cudaMemcpy(K_h, K_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);

	Out2File = fopen(Out2Filename, "w");
	for(int j = 0; j < Nx; ++j){
		double x = param.numin + j * param.dnu;
		fprintf(Out2File, "%.20g %.20g\n", x, K_h[j]);
	}
	fclose(Out2File);

	cudaDeviceSynchronize();
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	printf("Time after Line_kernel:  %g seconds\n", times + timems/1000000.0);

	error = cudaGetLastError();
	printf("error = %d = %s\n",error, cudaGetErrorString(error));

/*
for(int i = 0; i < 1; ++i){

 	double a = 0.0;
/*
	if(i == 0) a = 0.0;
	if(i == 1) a = 0.1;
	if(i == 2) a = 0.3;
	if(i == 3) a = 0.5;
	if(i == 4) a = 0.7;
	if(i == 5) a = 0.9;
	if(i == 6) a = 2.0;
	if(i == 7) a = 5.0;
	if(i == 8) a = 15.0;
* /
/*
	if(i == 0) a = 20.0;
	if(i == 1) a = 40.0;
	if(i == 2) a = 60.0;
	if(i == 3) a = 80.0;
	if(i == 4) a = 100.0;
	if(i == 5) a = 150.0;
	if(i == 6) a = 200.0;
* /	
	a = 200.0;

	double xmax = 100.0;
	Voigt_line_kernel <<< (Nx + 127) / 128, 128 >>> (a, dnu, K_d, Nx, xmax);

	cudaMemcpy(K_h, K_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);

for(int j = 0; j < Nx; ++j){
	double x = -xmax + j * 2.0 * xmax / ((double)(Nx));
	printf("%.20g %.20g %g\n", x, K_h[j], a);
}
}
*/

	return 0;
}
