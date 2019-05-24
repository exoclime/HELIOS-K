
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
__global__ void S2_kernel(double *nu_d, double *S_d, float *Sf_d, double *A_d, double *vy_d, float *vyf_d, double *ialphaD_d, double *n_d, double *delta_d, double *EL_d, int *ID_d, float *va_d, float *vb_d, float *vcut2_d, double *S1_d, float *S1f_d, const int NL, const double numin, const double dnu, const double cut, const int cutMode, int profile, int useIndividualX, const double T, const double P, const int kk){

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
		if(profile < 4){
			S1_d[id] = S_d[id] * vy_d[id] / M_PI;
		}
		else{
			S1_d[id] = S_d[id];
		}
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
//if(id < 100) printf("%d %g %g %g %g %g %g %g\n", id, nu_d[id], S_d[id], Sf_d[id], S1f_d[id], ialphaD_d[id], EL, exp(-c * nu / T));



	}
}
__global__ void Sf_kernel(double *nu_d, double *S_d, float *Sf_d, double *A_d, double *vy_d, float *vyf_d, double *ialphaD_d, double *n_d,  double *EL_d, double *S1_d, float *S1f_d, float *va_d, float *vb_d, float *vcut2_d, const int NL, const double numin, const double dnu, const double cut, const int cutMode, int profile, int useIndividualX, const double T, const double P, const int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){

		double nu = nu_d[id];
		double ialphaD = ialphaD_d[id] / nu;
		double EL = EL_d[id];
	
		double c = def_h * def_c / (def_kB * T);
		S_d[id] *= exp(-c * EL) * (1.0 - exp(-c * nu)) * ialphaD;	
		ialphaD_d[id] = ialphaD;	//inverse Doppler halfwidth

	
		double alphaL = vy_d[id];
		alphaL *= (P / def_PObar) * pow(def_T0 / T, n_d[id]);
		alphaL += A_d[id];				//1/cm
		vy_d[id] = alphaL * ialphaD;

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
		if(profile < 4){
			S1_d[id] = S_d[id] * vy_d[id] / M_PI;
		}
		else{
			S1_d[id] = S_d[id];
		}

		if(nu == 0){
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

//if( id < 100) printf("S %d %g %f %f %g\n", id, S_d[id], Sf_d[id], vcut2_d[id], cut);
//if(id > 156000 && id < 158000) printf("S %d %g %f\n", id, S_d[id], Sf_d[id]);
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
// profile = 1, voigt
// Individual X = 0
template <const int NBx, const int NBy, const int profile>
__global__ void Line4fA_kernel(float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int NL, const int nl){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;

	__shared__ double K_s[NBx * NBy];

	K_s[idx] = 0.0;
	__syncthreads();

	int ii, iii;
	float x, t1, xxyy;
	float S1, y, va, vb, vcut2;

	for(int iil = 0; iil < nl; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			va = va_d[iL];
			vb = vb_d[iL];
			vcut2 = vcut2_d[iL];

			for(int i = 0; i < NBx; ++i){
				ii = (i + idx) % NBx;
				iii = nstart + idy + ii;
				x = va + iii * vb;
//printf("x %d %d %g %g %g\n", ill, id, x, vb_s[ill], va_s[ill]);
				t1 = x * x;
				xxyy = t1 + y * y;
				if(profile == 1){  //voigt
					if(t1 < vcut2 && xxyy >= 1.0e6f){	
						//1 order Gauss Hermite Quadrature
						K_s[(i + idx) % (NBy * NBx)] += S1 / xxyy;
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
				__syncthreads();
			}
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx){
		K_d[nstart + idy + idx] += K_s[idx];
	}
}
// Case A
// profile = 1, voigt
// Individual X = 1
template <const int NBx, const int NBy, const int profile>
__global__ void Line4fAx_kernel(float *S1_d, float *vy_d, double *nu_d, double *ialphaD_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int NL, const int nl, double *x_d){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;

	__shared__ double K_s[NBx * NBy];
	__shared__ double x_s[NBx];

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

	for(int iil = 0; iil < nl; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			nu = nu_d[iL];
			ialphaD = ialphaD_d[iL];
			vcut2 = vcut2_d[iL];

			for(int i = 0; i < NBx; ++i){
				ii = (i + idx) % NBx;

				x = float((x_s[ii] - nu) * ialphaD);
//printf("x %d %d %g %g %g\n", ill, id, x, vb_s[ill], va_s[ill]);
				t1 = x * x;
				xxyy = t1 + y * y;
				if(profile == 1){  //voigt
					if(t1 < vcut2 && xxyy >= 1.0e6f){	
						//1 order Gauss Hermite Quadrature
						K_s[(i + idx) % (NBy * NBx)] += S1 / xxyy;
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
				__syncthreads();
			}
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx){
		K_d[nstart + idy + idx] += K_s[idx];
	}
}

// Case B
// profile = 1, voigt
// Individual X = 0
template <int NBx, int NBy>
__global__ void Line4fB_kernel(float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int NL, const int nl){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;

	__shared__ double K_s[NBx * NBy];

	K_s[idx] = 0.0;
	__syncthreads();

	int ii, iii;
	float x, t1, xxyy;
	float S1, y, va, vb, vcut2;

	for(int iil = 0; iil < nl; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			va = va_d[iL];
			vb = vb_d[iL];
			vcut2 = vcut2_d[iL];

			for(int i = 0; i < NBx; ++i){
				ii = (i + idx) % NBx;
				iii = nstart + idy + ii;
				x = va + iii * vb;
//printf("x %d %d %g %g %g\n", ill, id, x, vb_s[ill], va_s[ill]);
				t1 = x * x;
				xxyy = t1 + y * y;
				if(t1 < vcut2 && xxyy < 1.0e6f && xxyy >= 100.0f){	
					//2nd order Gauss Hermite Quadrature
					t1 *= 6.0f;
					float t2 = xxyy + 1.5f;
					float t3 = M_PIf * t2;

					float t4 = (t3 * (2.0f * t2 + xxyy) - 2.0f * t1) / (3.0f * xxyy * (t3 * t2 - t1));
					
					K_s[(i + idx) % (NBy * NBx)] += S1 * t4;
				}
				__syncthreads();
			}
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx){
		K_d[nstart + idy + idx] += K_s[idx];
	}
}

// Case B
// profile = 1, voigt
// Individual X = 1
template <int NBx, int NBy>
__global__ void Line4fBx_kernel(float *S1_d, float *vy_d, double *nu_d, double *ialphaD_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int NL, const int nl, double *x_d){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;

	__shared__ double K_s[NBx * NBy];
	__shared__ double x_s[NBx];

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

	for(int iil = 0; iil < nl; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			nu = nu_d[iL];
			ialphaD = ialphaD_d[iL];
			vcut2 = vcut2_d[iL];

			for(int i = 0; i < NBx; ++i){
				ii = (i + idx) % NBx;
	
				x = float((x_s[ii] - nu) * ialphaD);
//if(iL < 100) printf("%d %d %g %g %g\n", iL, nstart + idy + ii, nu, ialphaD, x);

//printf("x %d %d %g %g %g\n", ill, id, x, vb_s[ill], va_s[ill]);
				t1 = x * x;
				xxyy = t1 + y * y;
				if(t1 < vcut2 && xxyy < 1.0e6f && xxyy >= 100.0f){	
					//2nd order Gauss Hermite Quadrature
					t1 *= 6.0f;
					float t2 = xxyy + 1.5f;
					float t3 = M_PIf * t2;

					float t4 = (t3 * (2.0f * t2 + xxyy) - 2.0f * t1) / (3.0f * xxyy * (t3 * t2 - t1));
					
					K_s[(i + idx) % (NBy * NBx)] += S1 * t4;
				}
				__syncthreads();
			}
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx){
		K_d[nstart + idy + idx] += K_s[idx];
	}
}
// Case C
// profile = 1, voigt
// Individual X = 0
template <int NBx, int NBy>
__global__ void Line4fC_kernel(float *S1_d, float *vy_d, float *va_d, float *vb_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int NL, const int nl, const float a, const float b, const float c){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;

	__shared__ double K_s[NBx * NBy];

	K_s[idx] = 0.0;
	__syncthreads();

	int ii, iii;
	float x, t1, xxyy;
	float S1, y, va, vb, vcut2;

	for(int iil = 0; iil < nl; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			va = va_d[iL];
			vb = vb_d[iL];
			vcut2 = vcut2_d[iL];

			for(int i = 0; i < NBx; ++i){
				ii = (i + idx) % NBx;
				iii = nstart + idy + ii;
				x = va + iii * vb;
//printf("x %d %d %g %g %g\n", ill, id, x, vb_s[ill], va_s[ill]);
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

					K_s[(i + idx) % (NBy * NBx)] += S1 * t1 * b;
				}
				__syncthreads();
			}
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx){
		K_d[nstart + idy + idx] += K_s[idx];
	}
}

// Case C
// profile = 1, voigt
// Individual X = 1
template <int NBx, int NBy>
__global__ void Line4fCx_kernel(float *S1_d, float *vy_d, double *nu_d, double *ialphaD_d, float *vcut2_d, double *K_d, const int il, const int nstart, const int Nk, const int NL, const int nl, const float a, const float b, const float c, double *x_d){

	int idx = threadIdx.x;
	int idy = blockIdx.x * NBx;

	__shared__ double K_s[NBx * NBy];
	__shared__ double x_s[NBx];

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

	for(int iil = 0; iil < nl; iil += blockDim.x){
			
		int iL = iil + il + idx;
		if(iL < NL){ 
			S1 = S1_d[iL];
			y = vy_d[iL];
			nu = nu_d[iL];
			ialphaD = ialphaD_d[iL];
			vcut2 = vcut2_d[iL];

			for(int i = 0; i < NBx; ++i){
				ii = (i + idx) % NBx;

				x = float((x_s[ii] - nu) * ialphaD);

//printf("x %d %d %g %g %g\n", ill, id, x, vb_s[ill], va_s[ill]);
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

					K_s[(i + idx) % (NBy * NBx)] += S1 * t1 * b;
				}
				__syncthreads();
			}
		}
	}

	__syncthreads();

	if(idx < NBx){
		for(int j = 1; j < NBy; ++j){
			K_s[idx] += K_s[idx + j * NBx];
			__syncthreads();
		}
	}

	if(nstart + idy + idx < Nk && idx < NBx){
		K_d[nstart + idy + idx] += K_s[idx];
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

