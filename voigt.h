
// *********************************************
//This function calculates the Series Sigma1. Sigma2 and Sigma3 (Equations 27, 28, and 29) from Alg 916
//The parameter def_TOL sets a tolerance where to truncate the series
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

		
		if(fabs(s1 - sold1) < def_TOL) stop1 = 1;
		if(fabs(s2 - sold2) < def_TOL) stop2 = 1;
		if(fabs(s3 - sold3) < def_TOL) stop3 = 1;
		if(stop1 == 1 && stop2 ==1 && stop3 == 1) break;

		sold1 = s1;
		sold2 = s2;
		sold3 = s3;
//if(n >= 100-1) printf("Sigma Series did not converge\n");
	}
}
// *********************************************
//This function calculates the Series Sigma1. Sigma2 and Sigma3 (Equations 15, 16, and 17) from Alg 916
//The parameter def_TOL sets a tolerance where to truncate the series
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

		if(fabs(s1 - sold1) < def_TOL) stop1 = 1;
		if(fabs(s2 - sold2) < def_TOL) stop2 = 1;
		if(fabs(s3 - sold3) < def_TOL) stop3 = 1;
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

		if(fabs(s1 - sold1) < def_TOLF) stop1 = 1;
		if(fabs(s2 - sold2) < def_TOLF) stop2 = 1;
		if(fabs(s3 - sold3) < def_TOLF) stop3 = 1;
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
//The parameter def_TOL sets a tolerance where to truncate the series

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
//October 2016
// *************************************************
__global__ void S2_kernel(double *nu_d, double *S_d, double *A_d, double *vy_d, double *ialphaD_d, double *n_d, double *delta_d, double *EL_d, int *ID_d, const int NL, const double T, const double P, const int kk){

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
		alphaL *= (P / def_POatm) * pow(def_T0 / T, n_d[id]);
		alphaL += (A_d[id] / (4.0 * M_PI * def_c));				//1/cm
		vy_d[id] = alphaL * ialphaD;
		ID_d[id] = id;
	
		if(nu == 0.0){
			S_d[id] = 0.0;
			ialphaD_d[id] = 0.0;
			vy_d[id] = 0.0;
		}
//if(id < 100) printf("%d %g %g %g %g %g\n", id, nu_d[id], S_d[id], ialphaD_d[id], EL, exp(-c * nu / T));



	}
}

__global__ void Sf_kernel(double *nu_d, double *S_d, double *A_d, double *vy_d, double *ialphaD_d, double *n_d,  double *EL_d, int *ID_d, const int NL, const double c, const double T1, const double P, const int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){

		double nu = nu_d[id];
		double ialphaD = ialphaD_d[id] / nu;
		double EL = EL_d[id];
	
		S_d[id] *= exp(-c * EL) * (1.0 - exp(-c * nu)) * ialphaD;	
		ialphaD_d[id] = ialphaD;	//inverse Doppler halfwidth

	
		double alphaL = vy_d[id];
		alphaL *= (P / def_PObar) * pow(T1, n_d[id]);
		alphaL += A_d[id];				//1/cm
		vy_d[id] = alphaL * ialphaD;
		ID_d[id] = id;

		if(nu == 0){
			S_d[id] = 0.0;
			ialphaD_d[id] = 0.0;
			vy_d[id] = 0.0;
		}
//if(id < 10) printf("ialphaD %d %g %g %g\n", id, nu_d[id], ialphaD_d[id], vy_d[id]);


//if( id < 100) printf("S %d %g\n", id, S_d[id]);
//if(id > 156000 && id < 158000) printf("S %d %g\n", id, S_d[id]);
	}
}

// *************************************************
//This kernel creates float arrays from double arrays and arrays which depends only on the pre-sorted other Line arrays
//It is called after the sorting routine to avoid sorting twice the same data
//Author: Simon Grimm
//Date: Oktober 2019
// *************************************************
__global__ void S3_kernel(double *nu_d, double *S_d, double *S1_d, double *vy_d, double *ialphaD_d, float* Sf_d, float *S1f_d, float *vyf_d, float *vcut2_d, float *va_d, float *vb_d, const double cut, const int cutMode, const int profile, const double numin, const double dnu, int useIndividualX, const int NL, const int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){
		double ialphaD = ialphaD_d[id];
		vcut2_d[id] = (float)(cut * cut * ialphaD * ialphaD); //square of modified cut lenght

		if(cutMode == 1){
			vcut2_d[id] = (float)(cut * cut * vy_d[id] * vy_d[id]);
		}
		if(cutMode == 2){
			vcut2_d[id] = (float)(cut * cut);
		}

		if(profile < 4){
			S1_d[id] = S_d[id] * vy_d[id] / M_PI;
		}
		else{
			S1_d[id] = S_d[id];
		}

		double nu = nu_d[id];
		if(useIndividualX == 0){
			va_d[id] = (float)((numin - nu) * ialphaD);
			vb_d[id] = (float)(dnu * ialphaD);
//if(id == 304040) printf("vab %d %.20g %.20g %.20g %.20g\n", id, numin, nu, ialphaD, dnu);
//if(id == 304040) printf("vab %d %.20g %.20g %.20g %.20g\n", id, va_d[id], (numin - nu) * ialphaD, vb_d[id], dnu * ialphaD);
		}
		else{
			va_d[id] = (float)(-nu * ialphaD);
			vb_d[id] = (float)(ialphaD);
		}

		vyf_d[id] = (float)(vy_d[id]); 
		Sf_d[id] = (float)(S_d[id]); 
		S1f_d[id] = (float)(S1_d[id]);
	}

}

__global__ void L_kernelExomol(double *readBuffer_d, double *nu_d, double *S_d, double *EL_d, double *ialphaD_d, double* A_d, double *vy_d, double *n_d, const double defaultL, const double defaultn, const double gammaF, const double mass, const double T, const double Q, const double Abundance, const double Sscale, const int NL, const int kk){

	int id =  blockIdx.x * blockDim.x + threadIdx.x + kk;

	if(id < NL){

		nu_d[id] = readBuffer_d[id * 4 + 0];
		double S = readBuffer_d[id * 4 + 1];
		EL_d[id] = readBuffer_d[id * 4 + 2];
		double A = readBuffer_d[id * 4 + 3];

//if(id < 10) printf("%d %g %g %g %g\n", id, nu_d[id], S, EL_d[id], A);

		ialphaD_d[id] = def_c * sqrt( mass / (2.0 * def_kB * T));
		A_d[id] = A / (4.0 * M_PI * def_c);
		vy_d[id] = defaultL * gammaF;
		n_d[id] = defaultn;
		S_d[id] = S * Abundance * Sscale / Q;
// if(id < 100) printf("%d %g %g %g %g %g %g %g\n", id, nu_d[id], S_d[id], ialphaD_d[id], EL_d[id], Q, 0.0, vy_d[id]);
	}
}
__global__ void L_kernelKurucz(double *readBuffer_d, double *nu_d, double *S_d, double *EL_d, double *ialphaD_d, double* A_d, double *vy_d, double *n_d, const double defaultL, const double defaultn, const double gammaF, const double mass, const double T, const double Q, const double Abundance, const double Sscale, const int NL, const int kk){
	int id =  blockIdx.x * blockDim.x + threadIdx.x + kk;

	if(id < NL){

		nu_d[id] = readBuffer_d[id * 5 + 0];
		double S = readBuffer_d[id * 5 + 1];
		EL_d[id] = readBuffer_d[id * 5 + 2];
		double A = readBuffer_d[id * 5 + 3];
		double GammaN = readBuffer_d[id * 5 + 4];

//if(id < 10) printf("%d %g %g %g %g %g\n", id, nu_d[id], S, EL_d[id], A, GammaN);

		ialphaD_d[id] = def_c * sqrt( mass / (2.0 * def_kB * T));
		A_d[id] = (A + GammaN) / (4.0 * M_PI * def_c);
		vy_d[id] = defaultL * gammaF;
		n_d[id] = defaultn;
		S_d[id] = S * Abundance * Sscale / Q;

// if(id < 100) printf("%d %g %g %g %g %g %g %g\n", id, nu_d[id], S_d[id], ialphaD_d[id], EL_d[id], Q, 0.0, vy_d[id]);
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
		double aTOL = M_PI * sqrt(-1.0 / log(def_TOL * 0.5));
		double x = fabs(-xmax + id * 2.0 * xmax / ((double)(Nx)));
		K_d[id] = voigt_916(x, a, aTOL, id);
	}
}
// *************************************************
//This kernel calls directly the Voigt function
// It is usefull to test the profile 
//
//Author Simon Grimm
//November 2014
// *************************************************
__global__ void Voigt_2d_kernel(const double a, const double b, const double c, double *K_d, int Nx, int Ny, size_t pitch, double xMax, double yMax){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if(idx < Nx && idy < Ny){
		double x = idx * xMax / double(Nx - 1);
		double y = idy * yMax / double(Ny - 1);
		double s1, s2, s3;
		double ex2 = expf(-x * x);

		//Compute Sigma Series
		if(x != 0.0 && y != 0.0) Sigmab(x, y, s1, s2, s3, a, ex2, idx);

		double xy = x * y;
		double cos2xy = cosf(2.0 * xy);
		double sinxy = sinf(xy);

		double t1 = ex2 * erfcx(y) * cos2xy;
		double t2 = sinxy * ex2 * sinxy / y;
		double t3 = y * (-cos2xy * s1 + 0.5 * (s2 + s3));
		t1 += c * (t2 + t3);
		
		if(x == 0.0) t1 = erfcx(y);
		if(y == 0.0) t1 = ex2;

		//K_d[idy * Nx + idx] = t1 * b;
		double *row = (double *)(((char *)K_d)+(idy*pitch));
		row[idx] = t1 * b;
//if(idy == 0) printf("a %d %d %g %g %g\n", idx, idy, x, y, float(idy * Nx + idx));
//printf("%g %g %g %g %g %g\n", x, y, s1, s2, s3, K_d[idy * Nx + idx]);
	}
}
__global__ void Voigt_2df_kernel(const float a, const float b, const float c, float *K_d, int Nx, int Ny, size_t pitch, float xMax, float yMax){

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if(idx < Nx && idy < Ny){
		float x = idx * xMax / float(Nx - 1);
		float y = idy * yMax / float(Ny - 1);
		float s1, s2, s3;
		float ex2 = expf(-x * x);

		//Compute Sigma Series
		if(x != 0.0f && y != 0.0f) Sigmabf(x, y, s1, s2, s3, a, ex2, idx);

		float xy = x * y;
		float cos2xy = cosf(2.0f * xy);
		float sinxy = sinf(xy);

		float t1 = ex2 * erfcxf(y) * cos2xy;
		float t2 = sinxy * ex2 * sinxy / y;
		float t3 = y * (-cos2xy * s1 + 0.5f * (s2 + s3));
		t1 += c * (t2 + t3);
		
		if(x == 0.0f) t1 = erfcxf(y);
		if(y == 0.0f) t1 = ex2;

		//K_d[idy * Nx + idx] = t1 * b;
		float *row = (float *)(((char *)K_d)+(idy*pitch));
		row[idx] = t1 * b;
		//row[idx] = float(idy * Nx + idx);
//if(idy == 0) printf("a %d %d %g %g %g\n", idx, idy, x, y, float(idy * Nx + idx));
//printf("%g %g %g %g %g %g\n", x, y, s1, s2, s3, K_d[idy * Nx + idx]);
	}
}



// *************************************************
//This kernel initializes K_d with kmin
//
//Author Simon Grimm
//September 2016
// *************************************************
__global__ void InitialK_kernel(double *K_d, const int Nx, const double kmin, const int k){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + k;

	if(id < Nx){
		K_d[id] = kmin;
	}
}

// *******************************************************
// This kernel adds the different streams together to K_d
// Author: Simon Grimm
// May 2019
// ********************************************************
__global__ void AddKStreams_kernel(double *K_d, double *KS_d, const int nStreams, const int NX){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	if(id < NX){
		double K = 0.0;
		for(int i = 0; i < nStreams; ++i){
			K += KS_d[i * NX + id];
//if(id == 0) printf("KS %d %d %g %g\n", id, i, K, KS_d[i * NX + id]);
			KS_d[i * NX + id] = 0.0;
		}
		K_d[id] += K;
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
//printf("x_d %d %.20g\n", id, x_d[id]);
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

__global__ void print_kernel(double *nu_d, double *ialphaD_d, double *vy_d, int *ID_d, int N, int EE){

	
	for(int id = 0; id < N; ++id){
		printf("ID %d %g %g %g %d\n", EE, nu_d[id], ialphaD_d[id], vy_d[id], ID_d[id]);
	}
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
__global__ void Line2f_kernel(float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int nl, const int useIndividualX, const int Nxb, double *binBoundaries_d, const float a, const float b, const float c, const int profile){

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
	int bstart;
	if(id < Nk){
		int ii = nstart + id;
		double K = 0.0;
		if(useIndividualX != 0){
			int bin = ii / Nxb;
			bb = binBoundaries_d[bin];
			dnu = (binBoundaries_d[bin + 1] - bb) / ((double)(Nxb));
			bstart = bin * Nxb;
//printf("a %d %d %d %g %g %d %g\n", id, ii, bin, bb, dnu, bstart, dnu);
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
				if(E >= 0){
					x = bb * vb_s[ill] + va_s[ill] + (dnu * (ii - bstart) * vb_s[ill]);
				}
				else{
					x = bb * vb_s[ill] + va_s[ill] + (dnu * (ii * 10 - bstart) * vb_s[ill]);
				}
			}
//printf("x %d %d %g %g %g\n", ill, id, x, vb_s[ill], va_s[ill]);
			float t1 = x * x;
			float xxyy = t1 + y * y;
			if(profile == 1){			
				if(t1 < vcut2_s[ill]){	
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
			else if(profile == 2){
				//Lorentz profile
				if(E <= 0 && t1 < vcut2_s[ill]){
					K += S1_s[ill] / xxyy;
				}
			}
			else if(profile == 3){
				//Doppler profile
				if(E <= 0 && t1 < vcut2_s[ill]){
					K += S1_s[ill] * b * expf(-x * x);
				}
			}
			else if(profile == 4){
				//crossection as in Hill et all 2012
				float xp, xm, dd;
				if(useIndividualX == 0){
					xp = x + 0.5f * vb_s[ill];
					xm = x - 0.5f * vb_s[ill];
					dd = 1.0f / vb_s[ill];
				}
				else{
					xp = x + 0.5f * dnu * vb_s[ill];
					xm = x - 0.5f * dnu * vb_s[ill];
					dd = 1.0f / (vb_s[ill] * dnu);
				}
				if(E <= 0 && t1 < vcut2_s[ill]){
					K += S1_s[ill] * dd * 0.5f * (erff(xp) - erff(xm));
				}
			}
			
		}
		K_d[ii] += K;
	}
}


// Case A
// Individual X = 0
template <const int profile>
__global__ void Line6fA_kernel(float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il1, long long int *iiLimitsA0_d, const int Nk, const int NL, const int NBx, const int NBy, const int fK, const int LR){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;
	int by = blockIdx.y;

	int il = il1 + by * def_nlA;

	extern __shared__ double K_s[];

	int nstart = Nk;
	if(il < NL){
		nstart = iiLimitsA0_d[il / def_nlB];
	}

	K_s[idx] = 0.0;
	__syncthreads();

	int ii, iii;
	float x, t1, xxyy;
	float S1, y, va, vb, vcut2;

	for(int iil = 0; iil < def_nlA; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			va = va_d[iL];
			vb = vb_d[iL];
			vcut2 = vcut2_d[iL];
		}
		for(int i = 0; i < NBx; ++i){
			if(iL < NL){ 
				ii = (i + idx) % NBx;
				iii = nstart + idy + ii;
				x = va + iii * vb;
				t1 = x * x;
				xxyy = t1 + y * y;
				if(profile == 1){  //Voigt
					if(t1 < vcut2 && xxyy >= 1.0e6f){
						//1 order Gauss Hermite Quadrature
						if(LR == 10 || (LR == 11 && x <= 0.0 ) || (LR == 12 && x > 0.0  )){
//if(iii == 24982) printf("xA %d   %d %d %.20g %.20g %.20g %.20g %.20g\n", LR, iL, iii, x, S1, y, vb, va);
							K_s[(i + idx) % (NBy * NBx)] += S1 / xxyy;
//if(iii <= 24985 && iii >= 24970) printf("xA %d   %d %d %.20g %.20g %.20g %.20g %.20g %.30g %.30g\n", LR, iL, iii, x, S1, y, vb, va, 1.0 / xxyy, K_s[(i + idx) % (NBy * NBx)]);
						}
					}
				}
				if(profile == 2){  //Lorentz
					if(t1 < vcut2){	
						K_s[(i + idx) % (NBy * NBx)] += S1 / xxyy;
					}
				}
				if(profile == 3){  //Doppler
					if(t1 < vcut2){	
						K_s[(i + idx) % (NBy * NBx)] += S1 * expf(-t1) / sqrtf(M_PI);
					}
				}
				if(profile == 4){  //crossection as in Hill et all 2012
					if(t1 < vcut2){
						float xp = x + 0.5f * vb;
						float xm = x - 0.5f * vb;
						float dd = 1.0f / vb;
						K_s[(i + idx) % (NBy * NBx)] += S1 * dd * 0.5f * (erff(xp) - erff(xm));
					}
				}
			}
			__syncthreads();
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx && fK == 1){
		K_d[nstart + idy + idx  + by * Nk] += K_s[idx];
	}

}
// Case A
// Individual X = 1
template <const int profile>
__global__ void Line6fAX_kernel(float *S1_d, float *vy_d, double *nu_d, double *ialphaD_d, float *vcut2_d, double *K_d, const int il1, long long int *iiLimitsA0_d, const int Nk, const int NL, const int NBx, const int NBy, double *x_d, const int fK, const int LR){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;
	int by = blockIdx.y;

	int il = il1 + by * def_nlA;

	extern __shared__ double s_s[];
	double *K_s = s_s; 				//size NBx * NBy
	double *x_s = (double*)&K_s[NBx * NBy];		//size NBx

	int nstart = Nk;
	if(il < NL){
		nstart = iiLimitsA0_d[il / def_nlB];
	}

	K_s[idx] = 0.0;
	if(nstart + idy + idx < Nk && idx < NBx){
		x_s[idx] = x_d[nstart + idy + idx];
	}
	else if(idx < NBx){
		x_s[idx] = 0.0;
	}

	__syncthreads();

	int ii;
	float x, t1, xxyy;
	float S1, y, vcut2;
	double nu, ialphaD;

	for(int iil = 0; iil < def_nlA; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			nu = nu_d[iL];
			ialphaD = ialphaD_d[iL];
			vcut2 = vcut2_d[iL];
		}
		for(int i = 0; i < NBx; ++i){
			if(iL < NL){ 
				ii = (i + idx) % NBx;
				x = float((x_s[ii] - nu) * ialphaD);
				t1 = x * x;
				xxyy = t1 + y * y;
				if(profile == 1){  //Voigt
					if(t1 < vcut2 && xxyy >= 1.0e6f){
						//1 order Gauss Hermite Quadrature
						if(LR == 10 || (LR == 11 && x <= 0.0 ) || (LR == 12 && x > 0.0  )){
//printf("x %d   %d %g %g %g %g %g %g\n", LR, iL, x_s[ii], x, S1, y, nu, ialphaD);
							K_s[(i + idx) % (NBy * NBx)] += S1 / xxyy;
						}
					}
				}
				if(profile == 2){  //Lorentz
					if(t1 < vcut2){	
						K_s[(i + idx) % (NBy * NBx)] += S1 / xxyy;
					}
				}
				if(profile == 3){  //Doppler
					if(t1 < vcut2){	
						K_s[(i + idx) % (NBy * NBx)] += S1 * expf(-t1) / sqrtf(M_PI);
					}
				}
				if(profile == 4){  //crossection as in Hill et all 2012
					if(t1 < vcut2){
						float dnu;
						if(ii < NBx-1) dnu = x_s[ii + 1] - x_s[ii];
						else dnu = dnu = x_s[ii] - x_s[ii - 1];

						float xp = x + 0.5f * dnu * ialphaD;
						float xm = x - 0.5f * dnu * ialphaD;
						float dd = 1.0f / (ialphaD * dnu);
						K_s[(i + idx) % (NBy * NBx)] += S1 * dd * 0.5f * (erff(xp) - erff(xm));
					}
				}
			}
			__syncthreads();
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx && fK == 1){
		K_d[nstart + idy + idx + by * Nk] += K_s[idx];
	}

}

// Case B
// Individual X = 0
__global__ void Line6fB_kernel(float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il1, long long int *iiLimitsB0_d, const int Nk, const int NL, const int NBx, const int NBy, const int fK){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;
	int by = blockIdx.y;

	int il = il1 + by * def_nlB;

	extern __shared__ double K_s[];

	int nstart = Nk;
	if(il < NL){
		nstart = iiLimitsB0_d[il / def_nlB];
	}

	K_s[idx] = 0.0;
	__syncthreads();

	int ii, iii;
	float x, t1, xxyy;
	float S1, y, va, vb, vcut2;

	for(int iil = 0; iil < def_nlB; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			va = va_d[iL];
			vb = vb_d[iL];
			vcut2 = vcut2_d[iL];
		}
		for(int i = 0; i < NBx; ++i){
			if(iL < NL){ 
				ii = (i + idx) % NBx;
				iii = nstart + idy + ii;
				x = va + iii * vb;
				t1 = x * x;
				xxyy = t1 + y * y;
				if(t1 < vcut2 && xxyy < 1.0e6f && xxyy >= 100.0f){	
//if(iii == 24982) printf("xB %d %d %.20g %.20g %.20g %.20g %.20g\n", iL, iii, x, S1, y, vb, va);
					//2nd order Gauss Hermite Quadrature
					t1 *= 6.0f;
					float t2 = xxyy + 1.5f;
					float t3 = M_PIf * t2;

					float t4 = (t3 * (2.0f * t2 + xxyy) - 2.0f * t1) / (3.0f * xxyy * (t3 * t2 - t1));
					
					K_s[(i + idx) % (NBy * NBx)] += S1 * t4;
//if(iii <= 24985 && iii >= 24970) printf("xB %d %d %.20g %.20g %.20g %.20g %.20g %.30g %.30g\n", iL, iii, x, S1, y, vb, va, t4, K_s[(i + idx) % (NBy * NBx)]);
				}
			}
			__syncthreads();
		}
	}


	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx && fK == 1){
		K_d[nstart + idy + idx + by * Nk] += K_s[idx];
	}

}
// Case B
// Individual X = 1
__global__ void Line6fBX_kernel(float *S1_d, float *vy_d, double *nu_d, double *ialphaD_d, float *vcut2_d, double *K_d, const int il1, long long int *iiLimitsB0_d, const int Nk, const int NL, const int NBx, const int NBy, double *x_d, const int fK){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;
	int by = blockIdx.y;

	int il = il1 + by * def_nlB;

	extern __shared__ double s_s[];
	double *K_s = s_s;                              //size NBx * NBy
	double *x_s = (double*)&K_s[NBx * NBy];         //size NBx

	int nstart = Nk;
	if(il < NL){
		nstart = iiLimitsB0_d[il / def_nlB];
	}
	K_s[idx] = 0.0;
	if(nstart + idy + idx < Nk && idx < NBx){
		x_s[idx] = x_d[nstart + idy + idx];
	}
	else if(idx < NBx){
		x_s[idx] = 0.0;
	}

	__syncthreads();

	int ii;
	float x, t1, xxyy;
	float S1, y, vcut2;
	double nu, ialphaD;

	for(int iil = 0; iil < def_nlB; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			nu = nu_d[iL];
			ialphaD = ialphaD_d[iL];
			vcut2 = vcut2_d[iL];
		}
		for(int i = 0; i < NBx; ++i){
			if(iL < NL){ 
				ii = (i + idx) % NBx;
				x = float((x_s[ii] - nu) * ialphaD);
				t1 = x * x;
				xxyy = t1 + y * y;
				if(t1 < vcut2 && xxyy < 1.0e6f && xxyy >= 100.0f){	
//if(x_s[ii] > 25119.0 && x_s[ii] < 25121.0)  printf("x %d %d %.20g %.20g %.20g %.20g\n", iL, nstart + idy + ii, x_s[ii], x, y, S1);
					//2nd order Gauss Hermite Quadrature
					t1 *= 6.0f;
					float t2 = xxyy + 1.5f;
					float t3 = M_PIf * t2;

					float t4 = (t3 * (2.0f * t2 + xxyy) - 2.0f * t1) / (3.0f * xxyy * (t3 * t2 - t1));
					
					K_s[(i + idx) % (NBy * NBx)] += S1 * t4;
				}
			}
			__syncthreads();
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx && fK == 1){
		K_d[nstart + idy + idx + by * Nk] += K_s[idx];
	}
}


// Case C
// Individual X = 0
__global__ void Line6fC_kernel(float *S_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il1, long long int *iiLimitsC0_d, const int Nk, const int NL, const int NBx, const int NBy, const float a, const float b, const float c, const int fK){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;
	int by = blockIdx.y;

	int il = il1 + by * def_nlC;

	extern __shared__ double K_s[];

	int nstart = Nk;
	if(il < NL){
		nstart = iiLimitsC0_d[il / def_nlC];
	}

	K_s[idx] = 0.0;

	int ii, iii;
	float x, t1, xxyy;
	float S, y, va, vb, vcut2;

	__syncthreads();

	for(int iil = 0; iil < def_nlC; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S = S_d[iL];
			y = vy_d[iL];
			va = va_d[iL];
			vb = vb_d[iL];
			vcut2 = vcut2_d[iL];
		}

		for(int i = 0; i < NBx; ++i){
			if(iL < NL){ 
				ii = (i + idx) % NBx;
				iii = nstart + idy + ii;
				x = va + iii * vb;
				t1 = x * x;
				xxyy = t1 + y * y;
				if(t1 < vcut2 && xxyy < 100.0f){	

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

					K_s[(i + idx) % (NBy * NBx)] += S * t1 * b;
//if(iii <= 26253 && iii >= 26253) printf("xC %d %d %.20g %.20g %.20g %.20g %.20g %.30g %.30g\n", iL, iii, x, S, y, vb, va, t1, K_s[(i + idx) % (NBy * NBx)]);
				}
			}
			__syncthreads();
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx && fK == 1){
		K_d[nstart + idy + idx + by * Nk] += K_s[idx];
	}

}

// Case C
// Individual X = 1
__global__ void Line6fCX_kernel(float *S_d, float *vy_d, double *nu_d, double *ialphaD_d, float *vcut2_d, double *K_d, const int il1, long long int *iiLimitsC0_d, const int Nk, const int NL, const int NBx, const int NBy, const float a, const float b, const float c, double *x_d, const int fK){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;
	int by = blockIdx.y;

	int il = il1 + by * def_nlC;

	extern __shared__ double s_s[];
	double *K_s = s_s;                              //size NBx * NBy
	double *x_s = (double*)&K_s[NBx * NBy];         //size NBx

	int nstart = Nk;
	if(il < NL){
 		nstart = iiLimitsC0_d[il / def_nlC];
	}

	K_s[idx] = 0.0;
	if(nstart + idy + idx < Nk && idx < NBx){
		x_s[idx] = x_d[nstart + idy + idx];
	}
	else if(idx < NBx){
		x_s[idx] = 0.0;
	}

	int ii;
	float x, t1, xxyy;
	float S, y, vcut2;
	double nu, ialphaD;

	__syncthreads();

	for(int iil = 0; iil < def_nlC; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S = S_d[iL];
			y = vy_d[iL];
			nu = nu_d[iL];
			ialphaD = ialphaD_d[iL];
			vcut2 = vcut2_d[iL];
		}

		for(int i = 0; i < NBx; ++i){
			if(iL < NL){ 
				ii = (i + idx) % NBx;
				x = float((x_s[ii] - nu) * ialphaD);
//printf("x %d %d %g %g %g\n", iL, iii, x, vb, va);
				t1 = x * x;
				xxyy = t1 + y * y;
				if(t1 < vcut2 && xxyy < 100.0f){	

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

					K_s[(i + idx) % (NBy * NBx)] += S * t1 * b;
				}
			}
			__syncthreads();
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx && fK == 1){
		K_d[nstart + idy + idx + by * Nk] += K_s[idx];
	}
}

//fk = 0 does not update K_d, is used only for the tuning performance test
//no streams, but def_KSn y blocks
__host__ int Line6A_Call(Line L, Param param, double* KS_d, double * x_d, int il1, int NL, int nntt, int nk, int Nx, cudaStream_t Stream, int LR, int fK){


	int nt1 = 0;
	for(int i = 0; i < def_KSn; ++i){
		int il = il1 + i * def_nlA;
		if(il < NL){
			long long int ii00, ii11;
			if(LR == 10){
				ii00 = L.iiLimitsA0_h[il / def_nlA];
				ii11 = L.iiLimitsA1_h[il / def_nlA];
			}
			if(LR == 11){
				ii00 = L.iiLimitsAL0_h[il / def_nlA];
				ii11 = L.iiLimitsAL1_h[il / def_nlA];
			}
			if(LR == 12){
				ii00 = L.iiLimitsAR0_h[il / def_nlA];
				ii11 = L.iiLimitsAR1_h[il / def_nlA];
			}

//printf("%d %lld %lld\n", LR, L.iiLimitsAL0_h[0], L.iiLimitsAR0_h[0]);
			int nt = ii11 - ii00;
			if(ii00 > ii11) nt = 0; //this prevents from int overflow
			nt1 = max(nt1, nt);
int nnkk = (nt + nk - 1) / nk;
if(LR == 10 && (il % 10000 == 0 || il < 5 * def_nlA)) printf("AA   %d %lld %lld %d | blocks %d threads %d stream %d \n",il, ii00, ii11, nt1, nnkk, nntt, i);
if(LR == 11 && (il % 10000 == 0 || il < 5 * def_nlA)) printf("AAL  %d %lld %lld %d | blocks %d threads %d stream %d \n",il, ii00, ii11, nt1, nnkk, nntt, i);
if(LR == 12 && (il % 10000 == 0 || il < 5 * def_nlA)) printf("AAR  %d %lld %lld %d | blocks %d threads %d stream %d \n",il, ii00, ii11, nt1, nnkk, nntt, i);

		}
	}

	int nnkk = (nt1 + nk - 1) / nk;
	if(nnkk < 0) nnkk = 0;
	if(nnkk > 0){


		if(param.profile == 1){
			if(LR == 10){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 1 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 1 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 11){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 1 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 1 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 12){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 1 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 1 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
		}
		if(param.profile == 2){
			if(LR == 10){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 2 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 2 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 11){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 2 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 2 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 12){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 2 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 2 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
		}
		if(param.profile == 3){
			if(LR == 10){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 3 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 3 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 11){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 3 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 3 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 12){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 3 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 3 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
		}
		if(param.profile == 4){
			if(LR == 10){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 4 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 4 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsA0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 11){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 4 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 4 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAL0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
			if(LR == 12){
				if(param.useIndividualX == 0){
					Line6fA_kernel < 4 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, fK, LR);
				}
				else{
					Line6fAX_kernel < 4 > <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsAR0_d, Nx, NL, nk, nntt / nk, x_d, fK, LR);
				}
			}
		}
	}
	return nt1; 
}

//fk = 0 does not update K_d, is used only for the tuning performance test
//no streams, but def_KSn y blocks
__host__ int Line6B_Call(Line L, Param param, double* KS_d, double * x_d, int il1, int NL, int nntt, int nk, int Nx, cudaStream_t Stream, int fK){


	int nt1 = 0;
	for(int i = 0; i < def_KSn; ++i){
		int il = il1 + i * def_nlB;
		if(il < NL){
			long long int ii00 = L.iiLimitsB0_h[il / def_nlB];
			long long int ii11 = L.iiLimitsB1_h[il / def_nlB];

			int nt = ii11 - ii00;
			if(ii00 > ii11) nt = 0; //this prevents from int overflow
			nt1 = max(nt1, nt);
int nnkk = (nt + nk - 1) / nk;
if(il % 10000 == 0 || il < 5 * def_nlB) printf("BB  %d %lld %lld %d | blocks %d threads %d stream %d \n",il, ii00, ii11, nt1, nnkk, nntt, i);

		}
	}

	int nnkk = (nt1 + nk - 1) / nk;
	if(nnkk < 0) nnkk = 0;
	if(nnkk > 0){


		if(param.useIndividualX == 0){
			Line6fB_kernel <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsB0_d, Nx, NL, nk, nntt / nk, fK);
		}
		else{
			Line6fBX_kernel <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.S1f_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsB0_d, Nx, NL, nk, nntt / nk, x_d, fK);
		}
	}
	return nt1; 
}

//fk = 0 does not update K_d, is used only for the tuning performance test
//no streams, but def_KSn y blocks
__host__ int Line6C_Call(Line L, Param param, double* KS_d, double * x_d, int il1, int NL, int nntt, int nk, int Nx, double a, double b, double c, cudaStream_t Stream, int fK){


	int nt1 = 0;
	for(int i = 0; i < def_KSn; ++i){
		int il = il1 + i * def_nlC;
		if(il < NL){
			long long int ii00 = L.iiLimitsC0_h[il / def_nlC];
			long long int ii11 = L.iiLimitsC1_h[il / def_nlC];

			int nt = ii11 - ii00;
			if(ii00 > ii11) nt = 0; //this prevents from int overflow
			nt1 = max(nt1, nt);
int nnkk = (nt + nk - 1) / nk;
if(il % 10000 == 0 || il < 5 * def_nlC) printf("CC  %d %lld %lld %d | blocks %d threads %d stream %d \n",il, ii00, ii11, nt1, nnkk, nntt, i);

		}
	}

	int nnkk = (nt1 + nk - 1) / nk;
	if(nnkk < 0) nnkk = 0;
	if(nnkk > 0){


		if(param.useIndividualX == 0){
			Line6fC_kernel <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , nntt * sizeof(double), Stream >>> (L.Sf_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, KS_d, il1, L.iiLimitsC0_d, Nx, NL, nk, nntt / nk, a, b, c, fK);
		}
		else{
			Line6fCX_kernel <<< dim3(nnkk, def_KSn, 1), dim3(nntt, 1, 1) , (nntt + nk) * sizeof(double), Stream >>> (L.Sf_d, L.vyf_d, L.nu_d, L.ialphaD_d, L.vcut2_d, KS_d, il1, L.iiLimitsC0_d, Nx, NL, nk, nntt / nk, a, b, c, x_d, fK);
		}
	}
	return nt1; 
}


// ***************************************
// This kernel removes overlapping regions A and Al, AR

//Author: Simon Grimm
//Date: February 2020
// **************************************
__global__ void iiLimitsCheck(long long int *iiLimitsA0_d,  long long int *iiLimitsA1_d, long long int *iiLimitsAL0_d,  long long int *iiLimitsAL1_d, long long int *iiLimitsAR0_d,  long long int *iiLimitsAR1_d, const int NLimits){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	if(id < NLimits){
		if((iiLimitsAL1_d[id] > iiLimitsAL0_d[id]  || iiLimitsAR1_d[id] > iiLimitsAR0_d[id] )  && iiLimitsA1_d[id] > iiLimitsA0_d[id]){
			iiLimitsA0_d[id] = min(iiLimitsAL0_d[id], iiLimitsA0_d[id]);
			iiLimitsA1_d[id] = max(iiLimitsAR1_d[id], iiLimitsA1_d[id]);
			iiLimitsAL0_d[id] = 1.0e8;
			iiLimitsAL1_d[id] = -1.0e8;
			iiLimitsAR0_d[id] = 1.0e8;
			iiLimitsAR1_d[id] = -1.0e8;
		}
	}
}


__global__ void iiLimits_kernel(double *nuLimits0_d, double* nuLimits1_d, long long int *iiLimits0_d, long long int *iiLimits1_d, double *binBoundaries_d, const int NLimits, const double numin, const double dnu, const int Nx, const int useIndividualX, const int nbins, const int Nxb, const int EE){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;
	
	if(id < NLimits){

		double nu00 = nuLimits0_d[id];
		double nu11 = nuLimits1_d[id];
		
		long long int ii11 = 0LL;
		long long int ii00 = (long long int)(Nx);


		if(nu00 <= nu11){
			if(useIndividualX == 0){
				ii11 = (long long int)((nu11 - numin) / dnu) + 2;
				ii00 = (long long int)((nu00 - numin) / dnu) - 1;
			}
			else{
				for(int bin = nbins -1; bin >= 0; --bin){
					if(binBoundaries_d[bin] < nu11){
						double dnu = (binBoundaries_d[bin + 1] - binBoundaries_d[bin]) / ((double)(Nxb));
						int bstart = bin * Nxb;
						ii11 = (nu11 - binBoundaries_d[bin]) / dnu + bstart + 2;
						break;
					}
				}
				for(int bin = 0; bin < nbins; ++bin){
					if(binBoundaries_d[bin + 1] > nu00){
						double dnu = (binBoundaries_d[bin + 1] - binBoundaries_d[bin]) / ((double)(Nxb));
						int bstart = bin * Nxb;
						ii00 = (nu00 - binBoundaries_d[bin]) / dnu + bstart - 1;
						break;
					}
				}
			}
		
			ii11 = min((long long int)(Nx), ii11);
			ii00 = max(0LL, ii00);
		}	

		iiLimits0_d[id] = ii00;
		iiLimits1_d[id] = ii11;

//if(id == 683) printf("iilimitsK %d %d %g %g %lld %lld\n", EE, id, nu00, nu11, ii00, ii11);
	}

}

// ********************************************************
// This kernel finds the minimum and maxmun of the iiLimits
// and stores them in iiLimitsT
// it uses a parallel reduction sum with only 1 thread block

//Author: Simon Grimm
//Date: May 2020
// ********************************************************
template <int nb>
__global__ void iiLimitsMax_kernel(long long int *iiLimits0_d, long long int *iiLimits1_d, long long int *iiLimitsT_d, const int Nx, const int nl){


	int idy = threadIdx.x;

	__shared__ long long int ii00_s[nb];
	__shared__ long long int ii11_s[nb];

	ii11_s[idy] = 0LL;
	ii00_s[idy] = (long long int)(Nx);

	__syncthreads();

	for(int k = 0; k < nl; k += nb){
		if(idy + k < nl){
			ii00_s[idy] = min(ii00_s[idy], iiLimits0_d[idy + k]);
			ii11_s[idy] = max(ii11_s[idy], iiLimits1_d[idy + k]);
		}
	}
	__syncthreads();

	if(nb >= 512){
		if(idy < 256){
			ii00_s[idy] = min(ii00_s[idy], ii00_s[idy + 256]);
			ii11_s[idy] = max(ii11_s[idy], ii11_s[idy + 256]);
		}
	}
	__syncthreads();

	if(nb >= 256){
		if(idy < 128){
			ii00_s[idy] = min(ii00_s[idy], ii00_s[idy + 128]);
			ii11_s[idy] = max(ii11_s[idy], ii11_s[idy + 128]);
		}
	}
	__syncthreads();

	if(nb >= 128){
		if(idy < 64){
			ii00_s[idy] = min(ii00_s[idy], ii00_s[idy + 64]);
			ii11_s[idy] = max(ii11_s[idy], ii11_s[idy + 64]);
		}
	}
	__syncthreads();
	if(idy < 32){
		volatile long long int *ii00 = ii00_s;
		volatile long long int *ii11 = ii11_s;

		ii00[idy] = min(ii00[idy], ii00[idy + 32]);
		ii11[idy] = max(ii11[idy], ii11[idy + 32]);

		ii00[idy] = min(ii00[idy], ii00[idy + 16]);
		ii11[idy] = max(ii11[idy], ii11[idy + 16]);

		ii00[idy] = min(ii00[idy], ii00[idy + 8]);
		ii11[idy] = max(ii11[idy], ii11[idy + 8]);

		ii00[idy] = min(ii00[idy], ii00[idy + 4]);
		ii11[idy] = max(ii11[idy], ii11[idy + 4]);

		ii00[idy] = min(ii00[idy], ii00[idy + 2]);
		ii11[idy] = max(ii11[idy], ii11[idy + 2]);

		ii00[idy] = min(ii00[idy], ii00[idy + 1]);
		ii11[idy] = max(ii11[idy], ii11[idy + 1]);
	}
	__syncthreads();

	if(idy == 0){
		iiLimitsT_d[0] = ii00_s[0];
		iiLimitsT_d[1] = ii11_s[0];
	}

}



//This kernel finds the minimal and maximal index in wavenumbers for each block of lines
//EE 10: A
//EE 11: AL
//EE 12: AR
//EE 20: B
//EE 30: C
__global__ void nuLimits_kernel(double *nu_d, double *ialphaD_d, double *vy_d, float *vcut2_d, double* nuLimits0_d, double* nuLimits1_d, const double numin, const double numax, const int nl, const int NL, const int profile, const int EE){

	int idx = threadIdx.x;
	int id = blockIdx.x * nl + idx;

	__shared__ double nu0_s[32];
	__shared__ double nu1_s[32];

	int lane = threadIdx.x % warpSize;
	int warp = threadIdx.x / warpSize;

	if(warp == 0){
		nu0_s[threadIdx.x] = 1.0e8;
		nu1_s[threadIdx.x] = -1.0e8;
	}

	double nu0 = 1.0e8;
	double nu1 = -1.0e8;

	for(int i = 0; i < nl; i+= blockDim.x){
		if(id + i < NL){	
			double ialphaD = ialphaD_d[id + i];
			double vy = vy_d[id + i];
			double aD2 = 1.0 / (ialphaD * ialphaD);
			double aL2 = vy * vy * aD2;
			double vcut2 = vcut2_d[id + i];
			double cut = sqrt(vcut2 / (ialphaD * ialphaD));
			double nu = nu_d[id + i];

			double DnuB2 = 1.0e6 * aD2 - aL2;
			double DnuB = 0.0;
			if(DnuB2 > 0.0){
				DnuB = sqrt(DnuB2);
				DnuB = fmin(DnuB, cut);
			}

			double DnuC2 = 1.0e2 * aD2 - aL2;
			double DnuC = 0.0;
			if(EE >= 20 && DnuC2 > 0.0){
				DnuC = sqrt(DnuC2);
				DnuC = fmin(DnuC, cut);
			}

			if(EE == 10){
			//A
				if(DnuB == 0.0 || profile != 1){
					nu0 = fmin(nu - cut, nu0);
					nu1 = fmax(nu + cut, nu1);
				}
//if(blockIdx.x == 12 && nu >= numin - 100 * cut && nu <= numax + 100 * cut) printf("A__ %d %d %g %g %g %g %g %g\n", id, idx, nu, cut, DnuB, DnuC, nu0, nu1);			
			}

			if(EE == 11){
			//AL
				if(DnuB != 0.0 && DnuB != cut){
					nu0 = fmin(nu - cut, nu0);
					nu1 = fmax(nu - DnuB, nu1);
				}
			}

			if(EE == 12){
			//AR
				if(DnuB != 0.0 && DnuB != cut){
					nu0 = fmin(nu + DnuB, nu0);
					nu1 = fmax(nu + cut, nu1);
				}
			}

			if(EE == 20){
			//B
				if(DnuB != 0.0 && DnuB != DnuC){
					nu0 = fmin(nu - DnuB, nu0);
					nu1 = fmax(nu + DnuB, nu1);
//if(blockIdx.x == 12) printf("B__ %d %d %g %g %g %g %g %g\n", id, idx, nu, cut, DnuB, DnuC, nu0, nu1);			
				}
			}
			if(EE == 30){
			//C
				if(DnuC != 0){
					nu0 = fmin(nu - DnuC, nu0);
					nu1 = fmax(nu + DnuC, nu1);
				}
			}
//if(idx < 10  && blockIdx.x == 0) printf("%d__ %d %g %g %g %g %g %g\n", EE, idx, nu, cut, DnuB, DnuC, nu0, nu1);			

		}
	}
	__syncthreads();
//if(EE == 10 && id >= 12458 && id < 12558 ) printf("AAa %d %d %d %d %g %g\n", blockIdx.x, idx, warp, lane, nu0, nu1); 

	for(int i = 16; i >= 1; i/=2){
		nu0 = fmin(__shfl_xor_sync(0xffffffff, nu0, i, 32), nu0);
		nu1 = fmax(__shfl_xor_sync(0xffffffff, nu1, i, 32), nu1);
	}

	__syncthreads();
//if(EE == 10 && id >= 12458 && id < 12558 ) printf("AAb %d %d %d %d %g %g\n", blockIdx.x, idx, warp, lane, nu0, nu1); 


	if(blockDim.x > warpSize){
		//reduce across warps
		if(lane == 0){
			nu0_s[warp] = nu0;
			nu1_s[warp] = nu1;
		}
		__syncthreads();
//if(EE == 10 && idx < 32 && blockIdx.x == 29) printf("AAc %d %g %g\n", idx, nu0_s[idx], nu1_s[idx]); 
//if(EE == 11 && idx < 32 && blockIdx.x == 29) printf("AALc %d %g %g\n", idx, nu0_s[idx], nu1_s[idx]); 
//if(EE == 12 && idx < 32 && blockIdx.x == 29) printf("ARRc %d %g %g\n", idx, nu0_s[idx], nu1_s[idx]); 

		//reduce previous warp results in the first warp
	
		if(warp == 0){
			nu0 = nu0_s[threadIdx.x];
			nu1 = nu1_s[threadIdx.x];
//if(EE == 10) printf("AAc %d %g %g\n", threadIdx.x, nu0, nu1); 

			for(int i = 16; i >= 1; i/=2){
				nu0 = fmin(__shfl_xor_sync(0xffffffff, nu0, i, 32), nu0);
				nu1 = fmax(__shfl_xor_sync(0xffffffff, nu1, i, 32), nu1);
			}
//if(EE == 10 && idx == 0 && blockIdx.x == 29) printf("AAd %d %d %g %g\n", blockIdx.x, idx, nu0, nu1); 
//if(EE == 11 && idx == 0 && blockIdx.x == 29) printf("AALd %d %d %g %g\n", blockIdx.x, idx, nu0, nu1); 
//if(EE == 12 && idx == 0 && blockIdx.x == 29) printf("AARd %d %d %g %g\n", blockIdx.x, idx, nu0, nu1); 
//if(EE == 20 && idx == 0 && blockIdx.x == 29) printf("BBd %d %d %g %g\n", blockIdx.x, idx, nu0, nu1); 
		}
	}

	__syncthreads();
	if(threadIdx.x == 0){
		nuLimits0_d[blockIdx.x] = fmax(nu0, numin);
		nuLimits1_d[blockIdx.x] = fmin(nu1, numax);
//if(EE == 10) printf("nulimits A  %d  %d %g %g\n", EE, blockIdx.x, nuLimits0_d[blockIdx.x], nuLimits1_d[blockIdx.x]);
//if(EE == 11) printf("nulimits AL %d  %d %g %g\n", EE, blockIdx.x, nuLimits0_d[blockIdx.x], nuLimits1_d[blockIdx.x]);
//if(EE == 12) printf("nulimits AR %d  %d %g %g\n", EE, blockIdx.x, nuLimits0_d[blockIdx.x], nuLimits1_d[blockIdx.x]);
	}
}

