
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

// *************************************************
//This function calculates the Voigt profile V(x,y) as equation 13 from Zaghloul & Ali, Algorithm 916
//it calls the Sigma function
//The parameter TOL sets a tolerance where to truncate the series

//Author Simon Grimm, Adapted from Zaghloul & Ali, Algorithm 916
//November 2014
// **************************************************
__device__ double voigt_916(double x, double y, double a, int id){

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
	if(y == 0) t1 = exp(-x * x);
	//if(x*x + y*y > 1.0e18) t1 = y / (sqrt(M_PI) * (x * x + y * y));
	
	return t1;
}

// *************************************************
//This kernel calculates the integrate line strength, the Lorentz and the Doppler halfwidths
//
//Author Simon Grimm
//November 2014
// *************************************************
__global__ void S_kernel(double *nu_d, double *S_d, double *A_d, double *EL_d, double *alphaL_d, double *alphaD_d, double *n_d, double *mass_d, double *delta_d, double *Q_d, int *ID_d, int NL, double T, double P, int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	if(id < NL){
		double m = mass_d[id] / NA;			// mass in g

		double nu = nu_d[id] + delta_d[id] * P;		//read nu from alphaD
		if(nu == 0.0) nu = 0.0000001;
		nu_d[id] = nu;
		double S = S_d[id] / m;				//cm / g
		double EL = EL_d[id];  				//1/cm
		double Q = Q_d[id];				//Q0 / Q(T)
		double alphaL = alphaL_d[id];
		
		S_d[id] = S * Q * exp(-EL * def_h * def_c / (def_kB * T) + EL * def_h * def_c / (def_kB * T0)) * (1.0 - exp(-def_h * nu * def_c / (def_kB * T))) / (1.0 - exp(-def_h * nu * def_c / (def_kB * T0))); 
		alphaD_d[id] = def_c / nu * sqrt( m / (2.0 * def_kB * T));	//inverse Doppler halfwith
		alphaL *= P * pow(T0 / T, n_d[id]);
		alphaL += A_d[id] / (4.0 * M_PI * def_c);				//1/cm
		alphaL_d[id] = alphaL;
		ID_d[id] = id;
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
//January 2015
// *************************************************
__global__ void InitialK_kernel(double *K_d, double Nx, double kmin){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx;

	if(id < Nx){
		K_d[id] = kmin;
	}
}

// *************************************************
//This kernel initializes the binkeys with the bin number
//
//Nx is the total number of points,
//Nxb is the number of points ber bin
//
//Author Simon Grimm
//January 2015
// *************************************************
__global__ void binKey_kernel(int *binKey_d, int Nx, int Nxb){

	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id < Nx){
		int bin = id / Nxb;
		binKey_d[id] = bin;
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
//This kernel comnputes the limits on line indexes for each block in the Line kernel
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
__global__ void Cutoff_kernel(double *nu_d, int *ID_d, int2 *Limits_d, double *alphaL_d, double *ialphaD_d, int bl, double numin, double dnu, int NL, int n, double cut, int cutMode){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if(id < NL){


		// cutMode == 0: cut absolute values 
		if(cutMode == 1){			//Cut factors of Lorentz halfwidth
			cut *= alphaL_d[id];
		}
		else if(cutMode == 2){			//cut fctors of Lorentz / Gauss halfwidth -> y parameter in Voigt Profile
			cut /= ialphaD_d[id];
		}
		
		int il = ((nu_d[id] - numin + cut) / dnu) / bl;
		int ir = ((nu_d[id] - numin - cut) / dnu) / bl;

		il += 1;
		
		il = max(il, 0);
		ir = max(ir, 0);
		il = min(il, n);
		ir = min(ir, n);
	
		for(int j = ir; j <= il; ++j){ 
			atomicMin(&Limits_d[j].x, id);
		}
		for(int j = ir; j <= il; ++j){ 
			atomicMax(&Limits_d[j].y, id);
		}
//if(id  < 200) printf("%d %g %d %d %d %d\n", id, nu_d[id], il, i, ir, ID_d[id]);
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
__global__ void Line_kernel(double *nu_d, double *S_d, double *alphaL_d, double *alphaD_d, double *K_d, double dnu, double numin, int Nx, int NL, int2 *Limits_d, double cut, int cutMode, int nl, int ii, int kk){

	int idx = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idx + kk;

	__shared__ double nu_s[NB];
	__shared__ double S_s[NB];
	__shared__ double ialphaD_s[NB];
	__shared__ double y_s[NB];
	__shared__ double cut_s[NB];

	double K = 0.0;
	int2 Limits;
	if(blockIdx.x + kk / blockDim.x < (Nx + blockDim.x - 1) / blockDim.x){
		Limits= Limits_d[blockIdx.x + kk / blockDim.x];
	}
	else{
		Limits.x = 0;
		Limits.y = 0;
	}

	double a = M_PI * sqrt(-1.0 / log(TOL * 0.5));
	double sqln2 = 1.0;//sqrt(log(2.0)); //This factor will cancel out, becuase alphaD = sqln2 * vo sqrt(2kt/mc2)
	double isqrtpi = 1.0 / sqrt(M_PI);
	
	double nu = numin + id * dnu;

	for(int i = 0; i < nl; i += NB){
		if(i + idx + ii + Limits.x < NL){
			nu_s[idx] = nu_d[i + idx + ii + Limits.x];
			ialphaD_s[idx] = sqln2 * alphaD_d[i + idx + ii + Limits.x];
			y_s[idx] = alphaL_d[i + idx + ii + Limits.x] * ialphaD_s[idx];
			S_s[idx] = S_d[i + idx + ii + Limits.x] * ialphaD_s[idx];
			if(cutMode == 0) cut_s[idx] = cut;
			else if (cutMode == 1) cut_s[idx] = cut * alphaL_d[i + idx + ii + Limits.x];
			else if (cutMode == 2) cut_s[idx] = cut * y_s[idx];
		}
		else{
			nu_s[idx] = 0.0;
			S_s[idx] = 0.0;
			ialphaD_s[idx] = 0.0;
			y_s[idx] = 0.0;
			cut_s[idx] = 0.0;
		}
		__syncthreads();

# if PROFILE == 1
		for(int k = 0; k < NB; ++k){
			//Check smallest values for x and y
			if(i + k + ii + Limits.x < NL && fabs(nu - nu_s[k]) < cut_s[k]){
				double x = fabs((nu - nu_s[k]) * ialphaD_s[k]);
				double xxyy = x * x + y_s[k] * y_s[k];
//				if(__any(xxyy < 100)){
				if(xxyy < 100.0){
					K += S_s[k] * voigt_916(x, y_s[k], a, id) * isqrtpi;
				}
				else if(xxyy < 1.0e6){
//				else if(__any(xxyy < 1e6)){
					//2nd order Gauss Hermite Quadrature
					double t = y_s[k] / 3.0;
					double t1 = 2.0 * t / (M_PI * xxyy);
					double t2 = t * (xxyy + 1.5) / (M_PI * (xxyy + 1.5) * (xxyy + 1.5) - 4.0 * x * x * 1.5);
					K += S_s[k] * (t1 + t2);
				}
				else{
					//1 order Gauss Hermite Quadrature
					K += S_s[k] * y_s[k] / (M_PI * xxyy);
				}
			}
		}
#endif
# if PROFILE == 2
		for(int k = 0; k < NB; ++k){
			if(i + k + ii + Limits.x < NL && fabs(nu - nu_s[k]) < cut_s[k]){
				double x = fabs((nu - nu_s[k]) * ialphaD_s[k]);
				double xxyy = x * x + y_s[k] * y_s[k];
				K += S_s[k] * y_s[k] / (M_PI * xxyy);
			}	
		}
#endif
# if PROFILE == 3
		for(int k = 0; k < NB; ++k){
			if(i + k + ii + Limits.x < NL && fabs(nu - nu_s[k]) < cut_s[k]){
				double x = fabs((nu - nu_s[k]) * ialphaD_s[k]);
				K += S_s[k] * isqrtpi * exp(-x * x);
			}	
		}
#endif
		__syncthreads();
		if(i + ii + idx + Limits.x > Limits.y) break;

		__syncthreads();
	}
	if(id < Nx) K_d[id] += K;
}

