
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
__device__ void Sigma(const double x, const double y, double &s1, double &s2, double &s3, const double a, const double ex2, const int id){

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
__device__ void Sigmab(double x, const double y, double &s1, double &s2, double &s3, const double a, const double ex2, const int id){

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
__device__ void Sigmabf(float x, const float y, float &s1, float &s2, float &s3, const float a, const float ex2, const int id){

	s1 = 0.0f;
	float sold1 = s1;
	s2 = 0.0f;
	float sold2 = s2;
	s3 = 0.0f;
	float sold3 = s3;

	float f, f3p, f3n;
	float an, an3p, an3n;

	float yy = y * y;

	x = fabsf(x);

	int n0 = (int)(ceil(x / a)); //starting point for sigma3 series
	int n3p, n3n;

	int stop1 = 0;
	int stop2 = 0;
	int stop3 = 0;

	float e2axn = expf(-2.0f * a * x);

	for(int n = 1; n < 100; ++n){
		n3p = n0 + n - 1;
		n3n = n0 - n;
		an = a * n;
		an3p = a * n3p;
		an3n = a * n3n;

		f = 1.0f / (an * an + yy);
		f3p = 1.0f / (an3p * an3p + yy);
		f3n = 1.0f / (an3n * an3n + yy);

		s1 += f * expf(-(an * an + x * x));
		s2 += f * expf(-(an + x) * (an + x));
		s3 += f3p * expf(-(an3p - x) * (an3p - x));
		if(n3n >= 1) s3 += f3n * expf(-(an3n - x) * (an3n - x));

		if(fabs(s1 - sold1) < TOLF) stop1 = 1;
		if(fabs(s2 - sold2) < TOLF) stop2 = 1;
		if(fabs(s3 - sold3) < TOLF) stop3 = 1;
		if(stop1 == 1 && stop2 ==1 && stop3 == 1) break;

		sold1 = s1;
		sold2 = s2;
		sold3 = s3;
//		if(n >= 100-1) printf("Sigma Series did not converge %d\n", id);
	}
}

// *************************************************
//This function calculates the Voigt profile V(x,y) as equation 13 from Zaghloul & Ali, Algorithm 916
//it calls the Sigma function
//The parameter TOL sets a tolerance where to truncate the series

//Author Simon Grimm, Adapted from Zaghloul & Ali, Algorithm 916
//November 2014
// **************************************************
__device__ double voigt_916(const double x, const double y, const double a, const int id){

	double s1, s2, s3;
	double ex2 = exp(-x * x);

	//Compute Sigma Series
	if(x != 0.0 && y != 0.0) Sigmab(x, y, s1, s2, s3, a, ex2, id);

	double xy = x * y;
	double a2ipi = 2.0 * a / M_PI;
	double cos2xy = cos(2.0 * xy);
	double sinxy = sin(xy);

	double t1 = ex2 * erfcx(y) * cos2xy;
	t1 += a2ipi * x * sinxy * ex2 * sinxy / xy;
	t1 += a2ipi * y * (-cos2xy * s1 + 0.5 * (s2 + s3));
	
	if(x == 0) t1 = erfcx(y);
	if(y == 0) t1 = ex2;
	//if(x*x + y*y > 1.0e18) t1 = y / (sqrt(M_PI) * (x * x + y * y));
	
	return t1;
}
__device__ float voigt_916f(const float x, const float y, const float a, const int id){

	float s1, s2, s3;
	float ex2 = expf(-x * x);

	//Compute Sigma Series
	      if(x != 0.0f && y != 0.0f) Sigmabf(x, y, s1, s2, s3, a, ex2, id);

	float xy = x * y;
	float a2ipi = 2.0f * a / M_PIf;
	float cos2xy = cosf(2.0f * xy);
	float sinxy = sinf(xy);

	float t1 = ex2 * erfcxf(y) * cos2xy;
	t1 += a2ipi * x * sinxy * ex2 * sinxy / xy;
	t1 += a2ipi * y * (-cos2xy * s1 + 0.5f * (s2 + s3));

	if(x == 0) t1 = erfcxf(y);
	if(y == 0) t1 = ex2;
	//if(x*x + y*y > 1.0e18) t1 = y / (sqrt(M_PI) * (x * x + y * y));

	return t1;
}


// *************************************************
//This kernel calculates the integrated line strength, the Lorentz and the Doppler halfwidths
//
//Author Simon Grimm
//November 2014
// *************************************************
__global__ void S_kernel(double *nu_d, double *S_d, double *A_d, double *EL_d, double *vy_d, double *ialphaD_d, double *n_d, double *delta_d, int *ID_d, int NL, double T, double P, int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){

		double nu = nu_d[id] + delta_d[id] * P;
		if(nu == 0.0) nu = 0.0000001;
		nu_d[id] = nu;
		double S = S_d[id];				//cm / g
//printf("%d %g %g %g %g %g\n", id, nu_d[id], S_d[id], m, mass_d[id], Q_d[id]);
		double EL = EL_d[id];  				//1/cm
		double alphaL = vy_d[id];
		
		S_d[id] = (float)(S * exp(-EL * def_h * def_c / (def_kB * T) + EL * def_h * def_c / (def_kB * def_T0)) * (1.0 - exp(-def_h * nu * def_c / (def_kB * T))) / (1.0 - exp(-def_h * nu * def_c / (def_kB * def_T0)))); 
		ialphaD_d[id] /= nu;	//inverse Doppler halfwidth
		alphaL *= P * pow(def_T0 / T, n_d[id]);
		alphaL += A_d[id] / (4.0 * M_PI * def_c);				//1/cm
		vy_d[id] = alphaL * ialphaD_d[id];
		ID_d[id] = id;
//if(id < 1000) printf("%d %g %g %g %g\n", id, nu_d[id], S_d[id], ialphaD_d[id], vy_d[id]);
	}
}
// *************************************************
//This kernel calculates the integrated line strength, the Lorentz and the Doppler halfwidths
//
//Author Simon Grimm
//October 2016
// *************************************************
__global__ void S2_kernel(double *nu_d, double *S_d, float *Sf_d, double *A_d, double *vy_d, float *vyf_d, double *ialphaD_d, double *n_d, double *delta_d, double *EL_d, int *ID_d, float *va_d, float *vb_d, float *vcut2_d, double *S1_d, float *S1f_d, const int NL, const double numin, const double dnu, const double cut, const int cutMode, int useIndividualX, const double T, const double P, const int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){

		double nu = nu_d[id] + delta_d[id] * P;
		nu_d[id] = nu;
		double alphaL = vy_d[id];
		double ialphaD = ialphaD_d[id] / nu;
		double EL = EL_d[id];
	
		double c = def_h * def_c / def_kB;
		S_d[id] *= exp(-EL * c / T + EL * c / def_T0) * (1.0 - exp(-c * nu / T)) / (1.0 - exp(-c * nu / def_T0)) * ialphaD;
	
		ialphaD_d[id] = ialphaD;	//inverse Doppler halfwidth
		alphaL *= P * pow(def_T0 / T, n_d[id]);
		alphaL += (A_d[id] / (4.0 * M_PI * def_c));				//1/cm
		vy_d[id] = alphaL * ialphaD;
		ID_d[id] = id;
              
		if(useIndividualX == 0){
			va_d[id] = (float)((numin - nu) * ialphaD);
			vb_d[id] = (float)(dnu * ialphaD);
		}
		else{
			va_d[id] = (float)(-nu * ialphaD);
			vb_d[id] = (float)(ialphaD);
		}

		vcut2_d[id] = (float)(cut * cut * ialphaD * ialphaD); //square of modified cut lenght
		if(cutMode == 1){
			vcut2_d[id] = (float)(cut * cut * vy_d[id] * vy_d[id]);
		}
		if(cutMode == 2){
			vcut2_d[id] = (float)(cut * cut);
		}

		S1_d[id] = S_d[id] * vy_d[id] / M_PI;
		if(nu == 0.0){
			S_d[id] = 0.0;
			S1_d[id] = 0.0;
			ialphaD_d[id] = 0.0;
			vy_d[id] = 0.0;
			va_d[id] = 0.0f;
			vb_d[id] = 0.0f;
			vcut2_d[id] = 0.0f;
		}
		vyf_d[id] = (float)(vy_d[id]); 
		Sf_d[id] = (float)(S_d[id]); 
		S1f_d[id] = (float)(S1_d[id]); 
//if(id < 1000) printf("%d %g %g %g %g %g\n", id, nu_d[id], S_d[id], ialphaD_d[id], EL, exp(-c * nu / T));



	}
}
__global__ void Sf_kernel(double *S_d, float *Sf_d, double *vy_d, float *vyf_d, double *S1_d, float *S1f_d, const int NL, const int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){
		vyf_d[id] = (float)(vy_d[id]); 
		Sf_d[id] = (float)(S_d[id]); 
		S1f_d[id] = (float)(S1_d[id]); 
	}
}

// *************************************************
//This kernel calls directly the Voigt function
// It is usefull to test the profile 
//
//Author Simon Grimm
//November 2014
// *************************************************
__global__ void Voigt_line_kernel(double a, double dnu, double *K_d, double Nx, double xmax){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	if(id < Nx){
		double aTOL = M_PI * sqrt(-1.0 / log(TOL * 0.5));
		double x = fabs(-xmax + id * 2.0 * xmax / ((double)(Nx)));
		K_d[id] = voigt_916(x, a, aTOL, id);
	}
}
// *************************************************
//This kernel initializes K_d with kmin
//
//Author Simon Grimm
//September 2016
// *************************************************
__global__ void InitialK_kernel(double *K_d, const double Nx, const double kmin, const int k){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + k;

	if(id < Nx){
		K_d[id] = kmin;
	}
}

// *************************************************
// This kernel initializes the location of nu 
//
// Author Simon Grimm
// September 2016
// *************************************************
__global__ void setX_kernel(double *x_d, const double Nx, const double numin, const double dnu, const int Nxb, int useIndividualX, double *binBoundaries_d, const int k){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + k;

	if(id < Nx){
		if(useIndividualX == 0){
			x_d[id] = numin + id * dnu;
		}
		else{
			int bin = id / Nxb;
			double dnu = (binBoundaries_d[bin + 1] - binBoundaries_d[bin]) / ((double)(Nxb));
			int start = bin * Nxb;
			x_d[id] = binBoundaries_d[bin] + (id - start) * dnu;
		}
	}
}

// *************************************************
//This kernel initializes the binkeys with the bin number
//
//Nx is the total number of points,
//Nxb is the number of points ber bin
//
//Author Simon Grimm
//Septermber 2016
// *************************************************
__global__ void binKey_kernel(int *binKey_d, const int Nx, const int Nxb, double *binBoundaries_d, const int nbins, const double numax, double *x_d, int useIndividualX, const int k){

	int id = blockIdx.x * blockDim.x + threadIdx.x + k;

	if(id < Nx){
		//int bin = id / Nxb;
		//binKey_d[id] = bin;

		//if(useIndividualX == 1){
			double nu = x_d[id];
			binKey_d[id] = 0;

			for(int i = 0; i < nbins; ++i){
				double nl = binBoundaries_d[i];
				double nr = binBoundaries_d[i + 1];

				if(nl <= nu && nu < nr){
					binKey_d[id] = i;
					break;
				}
			}
			if(nu >= numax){
				binKey_d[id] = nbins; 
			}
		//}
	}
}
// *************************************************
//This kernel determines the starting index of the bins
//
//Nx is the total number of points,
//nbins is the number of bins
//
//Author Simon Grimm
//September 2016
// *************************************************
__global__ void binIndex_kernel(int *binKey_d, int *binIndex_d, const int Nx, const int nbins, const int k){

	int id = blockIdx.x * blockDim.x + threadIdx.x + k;

	if(id < Nx - 1){
		int bin = binKey_d[id];
		int bin1 = binKey_d[id + 1];

		if(bin1 > bin){
			binIndex_d[bin1] = id + 1;
		}
	}
	if(id == 0){
		binIndex_d[0] = 0;
		binIndex_d[nbins] = Nx;
		binIndex_d[nbins + 1] = Nx;
	}
}



//*************************************************
//This kernel copies an array a to b
//NL is the size of the array
//
//Author Simon Grimm
//January 2015
//*************************************************
__global__ void Copy_kernel(double *a_d, double *b_d, int NL, int k){

	int id = blockIdx.x * blockDim.x + threadIdx.x + k;

	if(id < NL){
		b_d[id] = a_d[id];
	}
}
__global__ void Copyf_kernel(float *a_d, double *b_d, int NL, int k){

	int id = blockIdx.x * blockDim.x + threadIdx.x + k;

	if(id < NL){
		b_d[id] = a_d[id];
	}
}

//*************************************************
//This kernel sorts the array a_d acording to the key in ID_d 
//and writes the result into b_d
//NL is the size of the array
//
//Author Simon Grimm
//January 2015
//*************************************************
__global__ void Sort_kernel(double *a_d, double *b_d, int *ID_d, int NL, int k){

        int id = blockIdx.x * blockDim.x + threadIdx.x + k;

        if(id < NL){
		b_d[id] = a_d[ID_d[id]];
	}
}
__global__ void Sortf_kernel(double *a_d, float *b_d, int *ID_d, int NL, int k){

        int id = blockIdx.x * blockDim.x + threadIdx.x + k;

        if(id < NL){
		b_d[id] = a_d[ID_d[id]];
	}
}

//*************************************************
//This kernel initializes the Limits arrays
//
//n is the size of the Limits_d array
//NL is the number of Lines
//Author Simon Grimm
//January 2015
//*************************************************
__global__ void setLimits_kernel(int2 *Limits_d, int n, int NL, double cut){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if(id < n){
		if(cut != 0.0){
			Limits_d[id].x = NL;
			Limits_d[id].y = 0;
		}
		else{
			Limits_d[id].x = 0;
			Limits_d[id].y = NL;
		}
	}
}

//*************************************************
//This kernel computes the limits of line indexes for each block in the Line kernel
// bl is the number of treads per block in the Line kernel
// n is the size of the Limits_d array
//
// cutMode == 0: cut at absolute values
// cutMode == 1: cut at Lorentz factors
// cutMode == 2: cut at Doppler factors
//
//Author Simon Grimm
//January 2015
//*************************************************
__global__ void Cutoff_kernel(double *nu_d, int *ID_d, int2 *Limits_d, double *vy_d, double *ialphaD_d, int bl, double numin, double dnu, int NL, int n, double cut, int cutMode, int Nx, double *x_d, int useIndividualX){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if(id < NL){

		double nu = nu_d[id];

		if(cutMode == 1){
			cut *= vy_d[id] / ialphaD_d[id];
		}
		else if(cutMode == 2){
			cut /= ialphaD_d[id];
		}
	
		//Determine the block index in x of the limits	
		int il = ((nu - numin + cut) / dnu) / bl;
		int ir = ((nu - numin - cut) / dnu) / bl;
		//if an irregular spacing in x is used, a binary search for the index is needed
		if(useIndividualX == 1){
			il = 0;
			ir = 0;
			int l = 0;
			int r = Nx - 1;
			int i;
			int m;
			for(i = 0; i < 10000; ++i){
				m = l + (r - l) / 2;
				if(m == 0){
					il = 0;
					break;
				}
				if(m >= Nx - 2){
					il = m / bl;
					break;
				}
//if(id == 224514) printf("%d %d %d %d %d %g %g\n", i, l, r, m, Nx, nu + cut, x_d[m]);
				if(x_d[m] <= (nu + cut) && (nu + cut) < x_d[m + 1]){
					il = m / bl;
					break;		
				}
				else{
					if(x_d[m] > (nu + cut)) r = m;
					else l = m;
				}
			}
			if(i >= 10000 -1) printf("Error: binary search for limits did not converge for il %d %g %g\n", id, nu + cut, x_d[m]);
			l = 0;
			r = Nx - 1;
			for(i = 0; i < 10000; ++i){
				m = l + (r - l) / 2;
				if(m == 0){
					ir = 0;
					break;
				}
				if(m >= Nx - 2){
					ir = m / bl;
					break;
				}
//if(id == 223381) printf("%d %d %d %d %d %g %g\n", i, l, r, m, Nx, nu - cut, x_d[m]);
				if(x_d[m] <= (nu - cut) && (nu - cut) < x_d[m + 1]){
					ir = m / bl;
					break;		
				}
				else{
					if(x_d[m] > (nu - cut)) r = m;
					else l = m;
				}
			}
			if(i >= 10000 -1) printf("Error: binary search for limits did not converge for ir %d %g %g\n", id, nu - cut, x_d[m]);
		}

		il += 1;
		
		il = max(il, 0);
		ir = max(ir, 0);
		il = min(il, n);
		ir = min(ir, n);

		if(il != 0 || ir != 0){	
			for(int j = ir; j <= il; ++j){ 
				atomicMin(&Limits_d[j].x, id);
			}
			for(int j = ir; j <= il; ++j){ 
				atomicMax(&Limits_d[j].y, id);
			}
		}
//if(id > NL - 50 || id < 50) printf("%d %g %d %d %d\n", id, nu_d[id], il, ir, ID_d[id]);
	}
}

//*************************************************
//This kernel computes the maximum number of lines for all thread blocks in the Line kernel
//
//n is the size of the Limits_d array
//NL is the number of Lines
//Author Simon Grimm
//January 2015
//*************************************************
__global__ void MaxLimits_kernel(int2 *Limits_d, int *MaxLimits_d, int n, int NL){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if(id < n){
		int il = Limits_d[id].x;
		int ir = Limits_d[id].y;
		atomicMax(&MaxLimits_d[0], ir - il);
	}
}

// *************************************************
//This kernel computes the line shape
// It uses patterns of shared memory the reduce global memory access
//
// NB is the number of threads per block
// nl is the number of lines per kernel launch nlmax
//
//Author Simon Grimm
//November 2014
// *************************************************
template <int NB>
__global__ void Line_kernel(float *S_d, float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, double *x_d, const int Nx, const int NL, int2 *Limits_d, const int nl, const int ii, const int kk, const int useIndividualX, const int Nxb, double *binBoundaries_d, const float a, const float b, const float c){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	__shared__ float S_s[NB];
	__shared__ float S1_s[NB];
	__shared__ float vy_s[NB];
	__shared__ float va_s[NB];
	__shared__ float vb_s[NB];
	__shared__ float vcut2_s[NB];

	double K = 0.0;
	int2 Limits;
	if(blockIdx.x + kk / blockDim.x < (Nx + blockDim.x - 1) / blockDim.x){
		Limits= Limits_d[blockIdx.x + kk / blockDim.x];
	}
	else{
		Limits.x = 0;
		Limits.y = 0;
	}

	for(int i = 0; i < nl; i += NB){
		if(i + idx + ii + Limits.x < NL){
			vy_s[idx] = vy_d[i + idx + ii + Limits.x];
			va_s[idx] = va_d[i + idx + ii + Limits.x];
			vb_s[idx] = vb_d[i + idx + ii + Limits.x];
			S_s[idx] = S_d[i + idx + ii + Limits.x];
			S1_s[idx] = S1_d[i + idx + ii + Limits.x];
			vcut2_s[idx] = vcut2_d[i + idx + ii + Limits.x];
		}
		else{
			S_s[idx] = 0.0f;
			S1_s[idx] = 0.0f;
			vy_s[idx] = 0.0f;
			va_s[idx] = 0.0f;
			vb_s[idx] = 0.0f;
			vcut2_s[idx] = 0.0f;
		}
		__syncthreads();
		float x;
# if PROFILE == 1

		for(int k = 0; k < NB; ++k){
			//Check smallest values for x and y
			if(useIndividualX == 0){
				x = va_s[k] + id * vb_s[k];
			}
			else{
				int bin = id / Nxb;
				double dnu = (binBoundaries_d[bin + 1] - binBoundaries_d[bin]) / ((double)(Nxb));
				int bstart = bin * Nxb;
				x = (float)(binBoundaries_d[bin] * vb_s[k] + va_s[k] + (id - bstart) * dnu * vb_s[k]);
			}
			float t1 = x * x;
			if(i + k + ii + Limits.x < NL && t1 < vcut2_s[k]){
				float y = vy_s[k];
				float xxyy = t1 + y * y;
				if(xxyy < 100.0f){
//printf("%d %g\n", id, y);
					float s1, s2, s3;
					float ex2 = expf(-x * x);

					//Compute Sigma Series
					if(x != 0.0f && y != 0.0f) Sigmabf(x, y, s1, s2, s3, a, ex2, ii);

					float xy = x * y;
					float cos2xy = cosf(2.0f * xy);
					float sinxy = sinf(xy);

					float t1 = ex2 * erfcxf(y) * cos2xy;
					float t2 = sinxy * ex2 * sinxy / y;
					float t3 = y * (-cos2xy * s1 + 0.5f * (s2 + s3));
					t1 += c * (t2 + t3);
					
					if(x == 0.0f) t1 = erfcxf(y);
					if(y == 0.0f) t1 = ex2;

					K += S_s[k] * t1 * b;
				}
				else if(xxyy < 1.0e6f){
					t1 *= 6.0f;
					float t2 = xxyy + 1.5f;
					float t3 = M_PIf * t2;

					float t4 = (t3 * (2.0f * t2 + xxyy) - 2.0f * t1) / (3.0f * xxyy * (t3 * t2 - t1));
					K += S1_s[k] * t4;
				}
				else{
					//1 order Gauss Hermite Quadrature
					K += S1_s[k] / xxyy;
//if(i + k + ii < 100) printf("%g %g %g\n", S1_s[k], x, y);
				}
			}
		}

#endif
# if PROFILE == 2
		for(int k = 0; k < NB; ++k){
			if(useIndividualX == 0){
				x = va_s[k] + id * vb_s[k];
			}
			else{
				int bin = id / Nxb;
				double dnu = (binBoundaries_d[bin + 1] - binBoundaries_d[bin]) / ((double)(Nxb));
				int bstart = bin * Nxb;
				x = binBoundaries_d[bin] * vb_s[k] + va_s[k] + (id - bstart) * dnu * vb_s[k];
			}
			float t1 = x * x;
			if(i + k + ii + Limits.x < NL && t1 < vcut2_s[k]){
				float xxyy = t1 + y_s[k] * y_s[k];
				K += S1_s[k] / xxyy;
			}	
		}
#endif
# if PROFILE == 3
		for(int k = 0; k < NB; ++k){
			if(useIndividualX == 0){
				x = va_s[k] + id * vb_s[k];
			}
			else{
				int bin = id / Nxb;
				double dnu = (binBoundaries_d[bin + 1] - binBoundaries_d[bin]) / ((double)(Nxb));
				int bstart = bin * Nxb;
				x = binBoundaries_d[bin] * vb_s[k] + va_s[k] + (id - bstart) * dnu * vb_s[k];
			}
			float t1 = x * x;
			if(i + k + ii + Limits.x < NL && t1 < vcut2_s[k]){
				K += S_s[k] * b * expf(-x * x);
			}	
		}
#endif
		__syncthreads();
		if(i + ii + idx + Limits.x > Limits.y) break;

		__syncthreads();
	}
	if(id < Nx) K_d[id] += K;
}


// *************************************************
//This kernel computes the line shape
//
// E = 0 first order
// E = 1 third order
// E = 2 higher oder
// E = -1 first order with reduced resolution in x
// E = 10, 11 correct boundary to lower resolution
// E = 12, 13 correct boundary to cut off 
//
//Author Simon Grimm
//October 2016
// *************************************************
template <int NB, const int E>
__global__ void Line2f_kernel(float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int nl, const int useIndividualX, const int Nxb, double *binBoundaries_d, const float a, const float b, const float c){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;
	__shared__ float S1_s[NB];
	__shared__ float vy_s[NB];
	__shared__ float va_s[NB];
	__shared__ float vb_s[NB];
	__shared__ float vcut2_s[NB];

	if(idx < nl){ 
		S1_s[idx] = S1_d[il + idx];
		vy_s[idx] = vy_d[il + idx];
		va_s[idx] = va_d[il + idx];
		vb_s[idx] = vb_d[il + idx];
		vcut2_s[idx] = vcut2_d[il + idx];

	}
	__syncthreads();
	float x;
	double dnu;
	double bb;
	if(id < Nk){
		int ii = nstart + id;
		double K = 0.0;
		if(useIndividualX != 0){
			int bin = ii / Nxb;
			bb = binBoundaries_d[bin];
			dnu = (binBoundaries_d[bin + 1] - bb) / ((double)(Nxb));
			int bstart = bin * Nxb;
			if(E >= 0){
				dnu *= (ii - bstart);
			}
			else{
				dnu *= (ii * 10 - bstart);

			}
		}
		for(int ill = 0; ill < nl; ++ill){
			float y = vy_s[ill];
			if(useIndividualX == 0){
				if(E >= 0){
					x = va_s[ill] + ii * vb_s[ill];
				}
				else{
					x = va_s[ill] + ii * 10 * vb_s[ill];
				}
			}
			else{
				x = bb * vb_s[ill] + va_s[ill] + dnu * vb_s[ill];
			}
			float t1 = x * x;
			float xxyy = t1 + y * y;

			if(t1 < vcut2_s[ill]){	

//printf("%g %g %g\n", S1_s[ill], x, y);

				if(E <= 0 && xxyy >= 1.0e6f){
				//1 order Gauss Hermite Quadrature
					K += S1_s[ill] / xxyy;
				}
				if(E == 1 && xxyy >= 100.0f && xxyy < 1.0e6f){
				//2nd order Gauss Hermite Quadrature
					t1 *= 6.0f;
					float t2 = xxyy + 1.5f;
					float t3 = M_PIf * t2;
				
					float t4 = (t3 * (2.0f * t2 + xxyy) - 2.0f * t1) / (3.0f * xxyy * (t3 * t2 - t1));
					K += S1_s[ill] * t4;
				}
				if(E == 2 && xxyy < 100.0f){
					float s1, s2, s3;
					float ex2 = expf(-x * x);

					//Compute Sigma Series
					if(x != 0.0 && y != 0.0) Sigmabf(x, y, s1, s2, s3, a, ex2, ii);

					float xy = x * y;
					float cos2xy = cosf(2.0f * xy);
					float sinxy = sinf(xy);

					float t1 = ex2 * erfcxf(y) * cos2xy;
					float t2 = sinxy * ex2 * sinxy / y;
					float t3 = y * (-cos2xy * s1 + 0.5f * (s2 + s3));
					t1 += c * (t2 + t3);
					
					if(x == 0.0f) t1 = erfcxf(y);
					if(y == 0.0f) t1 = ex2;

					K += S1_s[ill] * t1 * b;
				}
				if(E == 10){ //correct interpolation

					int i0 = (int)((-sqrtf((1.0e6f - y * y)) - va_s[ill]) / vb_s[ill]);
					i0 = (i0 / 10) * 10;
//printf("%d %d %d %d\n", id, ii, nstart, i0);
					float x0 = va_s[ill] + i0 * vb_s[ill];
					float xxyy0 = x0 * x0 + y * y;
					float Kc0 = S1_s[ill] / xxyy0;
					
					float Kc = Kc0 - Kc0 / 10.0 * (ii - i0);


					if(Kc >= 0.0f && Kc <= Kc0){
						if(xxyy >= 1.0e6f) K += S1_s[ill] / xxyy;
						K -= Kc;
					}
				}
				if(E == 11){ //correct interpolation

					int i0 = (int)((sqrtf((1.0e6f - y * y)) - va_s[ill]) / vb_s[ill]);
					i0 = (i0 / 10) * 10;
//printf("%d %d %d %d\n", id, ii, nstart, i0);
					float x0 = va_s[ill] + (i0 + 10) * vb_s[ill];
					float xxyy0 = x0 * x0 + y * y;
					float Kc0 = S1_s[ill] / xxyy0;
					
					float Kc = Kc0 / 10.0 * (ii - i0);


					if(Kc >= 0.0f && Kc <= Kc0){
						if(xxyy >= 1.0e6f) K += S1_s[ill] / xxyy;
						K -= Kc;
					}
					
				}
			}
			if(E == 12){ //correct interpolation
				int i0 = (int)((sqrtf((vcut2_s[ill])) - va_s[ill]) / vb_s[ill]);
				i0 = (i0 / 10) * 10;
//printf("%d %d %d %d\n", id, ii, nstart, i0);
				float x0 = va_s[ill] + i0 * vb_s[ill];
				float xxyy0 = x0 * x0 + y * y;
				float Kc0 = S1_s[ill] / xxyy0;
				
				float Kc = Kc0 - Kc0 / 10.0 * (ii - i0);


				if(Kc >= 0.0f && Kc <= Kc0){
					if(t1 < vcut2_s[ill]) K += S1_s[ill] / xxyy;
					K -= Kc;
				}
			}
			if(E == 13){ //correct interpolation
				int i0 = (int)((-sqrtf((vcut2_s[ill])) - va_s[ill]) / vb_s[ill]);
				i0 = (i0 / 10) * 10;
//printf("%d %d %d %d\n", id, ii, nstart, i0);
				float x0 = va_s[ill] + (i0 + 10) * vb_s[ill];
				float xxyy0 = x0 * x0 + y * y;
				float Kc0 = S1_s[ill] / xxyy0;
				
				float Kc = Kc0 / 10.0 * (ii - i0);


				if(Kc >= 0.0f && Kc <= Kc0){
					if(t1 < vcut2_s[ill]) K += S1_s[ill] / xxyy;
					K -= Kc;
				}
			}
		}
		K_d[ii] += K;
	}
}
// *************************************************
// This kernel initializes the location of nu 
//
// Author Simon Grimm
// September 2016
// *************************************************
__global__ void InterpolateX1_kernel(double *K_d, double *K1_d, const double Nx, const int Nxb, int useIndividualX, double *binBoundaries_d, const int k){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + k;

	if(id < Nx){
		if(useIndividualX == 0){
			int i1 = id / 10;

			double K0 = K1_d[i1];
			double K1 = K1_d[i1 + 1];

			K_d[id] += (K0 + (K1 - K0) / 10.0 * (id - i1 * 10));

		}
		else{
//			int bin = id / Nxb;
//			double dnu = (binBoundaries_d[bin + 1] - binBoundaries_d[bin]) / ((double)(Nxb));
//			int start = bin * Nxb;
//			x_d[id] = binBoundaries_d[bin] + (id - start) * dnu;
		}
	}

}
__global__ void InterpolateX2_kernel(double *K_d, double *K1_d, const double Nx, const int Nxb, int useIndividualX, double *binBoundaries_d, const int k){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + k;


	if(id < Nx){
		K_d[id] += K1_d[id];
	}
}

