
__constant__ double wS_c[3];


// ****************************************
// This kernel computes the Vandermonde like matrix for the least squares
// using Chebyshev polynomials
//
// nl is the bigger dimension, NC the smaller one
//
// Author: Simon Grimm
// December 2014
// *****************************************
__global__ void Vandermonde_kernel(double *A_d, double nl, int NC){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	int NL = (int)(nl);

	if(id < NL){

		double x = id * 2.0 / (nl - 1.0) - 1.0;
//printf("%d %g\n", id, x);
		A_d[id + NL * 0] = 1.0;
		A_d[id + NL * 1] = x;

		double t = x;
		double t0 = 1.0;

		for(int i = 2; i < NC; ++i){
			double tt = t;
			t = 2.0 * x * t - t0;
			t0 = tt;
			A_d[id + NL * i] = t;
		}
	}
}
// ****************************************
// This is just a helper function, to print a matrix to screen
//
// Author: Simon Grimm
// December 2014
// *****************************************
__device__ void printA(double *A_d, int NL, int NC){

	printf("A \n \n");

	for(int j = 0; j < NL; ++j){
		for(int i = 0; i < NC; ++i){
			printf("%g  ", A_d[j + NL * i]);

		}
		printf("\n");
	}
}



// ****************************************
// This kernel performs a QR decomposition of the Matrix V_d
// it overwrites the upper right part of it with R, except of
// the diagonal part, which is written in D_d. C_d contains 
// the Householder scalar c. The rest of V_d is filled with the 
// Householder vectors.
//
// NL is the bigger dimension of V_d, NC the smaller one
// nb is the number of threads per block
// The kernel must be launched with only one block
//
// Author: Simon Grimm
// December 2014
// *****************************************
template <int nb>
__global__ void QR_kernel(double *V_d, double *C_d, double *D_d, int NL, int NC){

	int idy = threadIdx.x;

	__shared__ double a_s[nb];
	double scale;

	for(int i = 0; i < NC; ++i){ // if V is square then i < NC - 1

		a_s[idy] = 0.0;

		for(int k = 0; k < NL; k += nb){
			if(idy + k < NL && idy + k >= i){
				a_s[idy] = fmax(a_s[idy], fabs(V_d[(idy + k) + NL * i]));
			}
		}
		__syncthreads();

		if(nb >= 512){
			if(idy < 256){
				a_s[idy] = fmax(a_s[idy], a_s[idy + 256]);
			}
		}
		__syncthreads();

		if(nb >= 256){
			if(idy < 128){
				a_s[idy] = fmax(a_s[idy], a_s[idy + 128]);
			}
		}
		__syncthreads();

		if(nb >= 128){
			if(idy < 64){
				a_s[idy] = fmax(a_s[idy], a_s[idy + 64]);
			}
		}
		__syncthreads();
		if(idy < 32){
			volatile double *a = a_s;
			a[idy] = fmax(a[idy], a[idy + 32]);
			a[idy] = fmax(a[idy], a[idy + 16]);
			a[idy] = fmax(a[idy], a[idy + 8]);
			a[idy] = fmax(a[idy], a[idy + 4]);
			a[idy] = fmax(a[idy], a[idy + 2]);
			a[idy] = fmax(a[idy], a[idy + 1]);
		}
		__syncthreads();
		scale = a_s[0];
		__syncthreads();
		a_s[idy] = 0.0;
		__syncthreads();
		for(int k = 0; k < NL; k += nb){
			if(idy + k < NL && idy + k >= i){
				double V = V_d[(idy + k) + NL * i] / scale;
				V_d[(idy + k) + NL * i] = V;
//if(i == 1) printf("s %d %d %d %g %g %d %d\n", (idy + k) + NL * i, k, i, V_d[(idy + k) + NL * i], scale, NL, NC);
				a_s[idy] += V * V;
			}
		}
		__syncthreads();

		if(nb >= 512){
			if(idy < 256){
				a_s[idy] += a_s[idy + 256];
			}
		}
		__syncthreads();

		if(nb >= 256){
			if(idy < 128){
				a_s[idy] += a_s[idy + 128];
			}
		}
		__syncthreads();

		if(nb >= 128){
			if(idy < 64){
				a_s[idy] += a_s[idy + 64];
			}
		}
		__syncthreads();

		if(idy < 32){
			volatile double *a = a_s;
			a[idy] += a[idy + 32];
			a[idy] += a[idy + 16];
			a[idy] += a[idy + 8];
			a[idy] += a[idy + 4];
			a[idy] += a[idy + 2];
			a[idy] += a[idy + 1];
		}
		if(idy == 0){
			double s = sqrt(a_s[0]);
			double V = V_d[i + NL * i];
			double sigma = s;
			if(V < 0.0) sigma = -s;
			V_d[i + NL * i] += sigma;
			C_d[i] = sigma * (V + sigma);
			D_d[i] = -scale * sigma;
		}
		__syncthreads();
		for(int j = i + 1; j < NC; ++j){
			a_s[idy] = 0.0;
			__syncthreads();
			for(int k = 0; k < NL; k += nb){
				if(idy + k < NL && idy + k >= i){
					a_s[idy] += V_d[(idy + k) + NL * i] * V_d[(idy + k) + NL * j];
				}
			}
			__syncthreads();

			if(nb >= 512){
				if(idy < 256){
					a_s[idy] += a_s[idy + 256];
				}
			}
			__syncthreads();

			if(nb >= 256){
				if(idy < 128){
					a_s[idy] += a_s[idy + 128];
				}
			}
			__syncthreads();

			if(nb >= 128){
				if(idy < 64){
					a_s[idy] += a_s[idy + 64];
				}
			}
			__syncthreads();

			if(idy < 32){
				volatile double *a = a_s;
				a[idy] += a[idy + 32];
				a[idy] += a[idy + 16];
				a[idy] += a[idy + 8];
				a[idy] += a[idy + 4];
				a[idy] += a[idy + 2];
				a[idy] += a[idy + 1];
			}
			__syncthreads();
			double t = a_s[0] / C_d[i];
			for(int k = 0; k < NL; k += nb){
				if(idy + k < NL && idy + k >= i){
					V_d[(idy + k) + NL * j] -= t * V_d[(idy + k) + NL * i];
				}
			}
			__syncthreads();
		}//end if j loop
	}// end of i loop
	__syncthreads();
}

// ****************************************
// This kernel performes a least square fit by using the previously computed 
// QR decomposition of the matrix V_d. It uses V_d, C_d and D_d from the QR_kernel
// the array b_d contains initially the right hand solution vector of V * x = b
// The kernel first computes Qt * b and solves then R * x = Qt * b
// The solution x is written in the beginning of b_d
//
// NL is the bigger dimension of V_d, NC the smaller
// nb is the number of threads per block
// The kernel must be launched with only one block
//
// Author: Simon Grimm
// December 2014
// *****************************************
template <int nb>
__global__ void leastSquare_kernel(double *V_d, double *C_d, double *D_d, double *b_d, int NL, int NC){
	int idy = threadIdx.x;
	int idx = blockIdx.x;

	__shared__ double a_s[nb];

	//Compute Qt b
	for(int j = 0; j < NC - 1; ++j){
		a_s[idy] = 0.0;
		__syncthreads();	
		for(int k = 0; k < NL; k += nb){
			if(idy + k < NL && idy + k >= j){
				double b = b_d[idy + k + idx * NL];
//printf("%d %d %d %g\n", (idy + k) + NL * j, k, j, V_d[(idy + k) + NL * j]);
				a_s[idy] += V_d[(idy + k) + NL * j] * b;
			}
		}
		__syncthreads();

		if(nb >= 512){
			if(idy < 256){
				a_s[idy] += a_s[idy + 256];
			}
		}
		__syncthreads();

		if(nb >= 256){
			if(idy < 128){
				a_s[idy] += a_s[idy + 128];
			}
		}
		__syncthreads();

		if(nb >= 128){
			if(idy < 64){
				a_s[idy] += a_s[idy + 64];
			}
		}
		__syncthreads();

		if(idy < 32){
			volatile double *a = a_s;
			a[idy] += a[idy + 32];
			a[idy] += a[idy + 16];
			a[idy] += a[idy + 8];
			a[idy] += a[idy + 4];
			a[idy] += a[idy + 2];
			a[idy] += a[idy + 1];
		}
		__syncthreads();
		double t = a_s[0] / C_d[j];
		for(int k = 0; k < NL; k += nb){
			if(idy + k < NL && idy + k >= j){
				b_d[idy + k + idx * NL] -= t * V_d[(idy + k) + NL * j];
			}
		}
		__syncthreads();
	}//end if j loop in Qt b

	//solve R x = Qt b

	if(idy == 0) b_d[NC - 1 + idx * NL] /= D_d[NC - 1];
	__syncthreads();

	for(int i = NC - 2; i >= 0; --i){
		a_s[idy] = 0.0;
		__syncthreads();
		//Assume NC < nb
		if(idy < NC && idy >= i + 1){
			a_s[idy] += V_d[i + NL * idy] * b_d[idy + idx * NL];
		}
		__syncthreads();

		if(NC >= 512){
			if(idy < 256){
				a_s[idy] += a_s[idy + 256];
			}
		}
		__syncthreads();

		if(NC >= 256){
			if(idy < 128){
				a_s[idy] += a_s[idy + 128];
			}
		}
		__syncthreads();

		if(NC >= 128){
			if(idy < 64){
				a_s[idy] += a_s[idy + 64];
			}
		}
		__syncthreads();

		if(idy < 32){
			volatile double *a = a_s;
			a[idy] += a[idy + 32];
			a[idy] += a[idy + 16];
			a[idy] += a[idy + 8];
			a[idy] += a[idy + 4];
			a[idy] += a[idy + 2];
			a[idy] += a[idy + 1];
		}
		__syncthreads();
		if(idy == 0){
			b_d[i + idx * NL] = (b_d[i + idx * NL] - a_s[0]) / D_d[i];

		}
	}//end of i loop in R x = Qtb
}


// ****************************************
// This kernel finds bins starting with kmin and stores the starting index
//
// Author: Simon Grimm
// February 2016
// *****************************************
__global__ void findCut_kernel(double *K_d, int NL, int Nxb, double kmin, int *Nxmin_d, int nbins){

	int id = threadIdx.x + blockIdx.x * blockDim.x;
	int ib = id / Nxb;

	if(id < NL - 1 && ib < nbins){
		double K = K_d[id];
		double K1 = K_d[id + 1];

		if(K <= kmin && K1 > kmin){
			int n = id - ib * Nxb + 1;
//printf("cut bin %d %d %d %d %g\n", id, ib, n, Nxb, n / ((double)(Nxb - 1))) ;
			Nxmin_d[ib] = n;

		}
		//find complete empty bins
		if(K <= kmin && id % Nxb == Nxb - 1){
//printf("empty bin %d %d\n", ib, Nxb);
			Nxmin_d[ib] = Nxb;
		}
		if(K <= kmin && id == NL - 3){
//printf("empty last bin %d %d %d\n", ib, id, Nxb);
			Nxmin_d[ib] = Nxb;
		}
	}

}

// ****************************************
// This kernel rescales bins starting with kmin to [0:1] starting from 
// the first entry bigger than kmin, stored in Nxmin 
//
// Author: Simon Grimm
// February 2015
// *****************************************
template <int nb>
__global__ void rescale_kernel(int *Nxmin_d, double *K_d, double *K2_d, int Nxb, double kmin, int f){

	int idy = threadIdx.x;
	int idx = blockIdx.x;

	int Nxmin = Nxmin_d[idx];

	if(Nxmin > 0){
		for(int k = 0; k < Nxb; k += nb){
			if(idy + k < Nxb){
				if(f == 1){
					double ii = Nxmin + (1.0 - Nxmin / ((double)(Nxb - 1))) * (k + idy); //required index position
//if(idx == 85) printf("%d %d %d %g\n", idx, Nxmin, k + idy, ii);
					if(ii >= Nxb - 1) ii = 0.999999 * (Nxb - 1);
					int il = ii / ((double)(Nxb)) * Nxb; //left index
					double Kl = K_d[idx * Nxb + il];
					double Kr = K_d[idx * Nxb + il + 1];
					double Ki = Kl + (Kr - Kl) * (ii - il);
					K2_d[idx * Nxb + k + idy] = Ki;
				}
				if(f == -1){
				//inverse transformation
					double Ki = kmin;
					if(k + idy >= Nxmin){
						double ii = (k + idy - Nxmin) / ((double)(Nxb - Nxmin - 1)) * (Nxb - 1);
						if(ii >= Nxb - 1) ii = 0.999999 * (Nxb - 1);
						int il = ii / ((double)(Nxb)) * Nxb; //left index
						double Kl = K_d[idx * Nxb + il];
						double Kr = K_d[idx * Nxb + il + 1];
//if((idy + k) % 100 == 0) printf("%d %g %d %d %g\n", idy + k, ii, k + idy - Nxmin, Nxb - 1, (k + idy - Nxmin) / ((double)(Nxb - Nxmin - 1)) * (Nxb - 1));
						Ki = Kl + (Kr - Kl) * (ii - il);
					}
					K2_d[idx * Nxb + k + idy] = Ki;
				}
//if(idx == 3) printf("K %d %d %g %g %g\n", idx * Nxb + k + idy, Nxmin, Kl, Kr, Ki);
			}
		}
	}
}

// ****************************************
// This kernel copies rescaled bins into K_d
//
// Author: Simon Grimm
// February 2015
// *****************************************
template <int nb>
__global__ void copyK2_kernel(int *Nxmin_d, double *K_d, double *K2_d, int Nxb){

	int idy = threadIdx.x;
	int idx = blockIdx.x;

	int Nxmin = Nxmin_d[idx];

	if(Nxmin > 0){

		for(int k = 0; k < Nxb; k += nb){
			if(idy + k < Nxb){
				K_d[idx * Nxb + k + idy] = K2_d[idx * Nxb + k + idy];
			}
		}
	}

}



// ****************************************
// This kernel computes the log of an array K_d
// NL is the lenght of the array
//
// Author: Simon Grimm
// January 2015
// *****************************************
__global__ void lnK_kernel(double *K_d, int NL, double unitScale){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if(id < NL){
		double K = K_d[id] * unitScale;
		K_d[id] = log(K);
//printf("%d %g %g\n", id, K_d[id], K);
	}
}

// ****************************************
// This kernel computes exp(f(x)) using the Chebyshev coefficients
// calculated previosly and stored in b_d.
// The kernel must be called at least with def_NmaxSample threads, the number of
// chebyshev coefficients
//
// The result is written again in b_d
//
// NC is the number of chebyshev coefficients
// NL the number of data points
//
// Author: Simon Grimm
// January 2015
// *****************************************
__global__ void expfx_kernel(double *b_d, int NC, int NL, double unitScale){

	int idy = threadIdx.x;
	int idx = blockIdx.x;
	__shared__ double b_s[def_NmaxSample];

	for(int k = 0; k < def_NmaxSample; k += blockDim.x){
		if(k + idy < def_NmaxSample){
			b_s[k + idy] = b_d[k + idy + idx * NL];
		}
	}
	__syncthreads();

	for(int k = 0; k < NL; k += blockDim.x){
		if(idy + k < NL){
			double x = -1.0 + (k + idy) * 2.0 / ((double)(NL));
			double d1 = 0.0;
			double d2 = 0.0;
			for(int i = NC - 1; i >= 1; --i){
				double t = d1;
				d1 = 2.0 * x * d1 - d2 + b_s[i];
				d2 = t;
			} 
			double f = x * d1 - d2 + 0.5 * b_s[0] + 0.5 * b_s[0];
			b_d[idy + k + idx * NL] = exp(f) / unitScale;
		}
	}
}

// ****************************************
// This kernel integrates a function exp(-K_d * m) using the extended Simpsons rule:
// I = 1/h (3/8 f1 + 7/6 f2 + 23/24/ f3 + f4 + ... + FN-3 + 23/23 fN-2 + 7/6 fN-1 + 3/8 fN) + O(1/N^4)
//
// The result is written in Tr_d[0]
//
// NL is the number of data points in K_d
// nb is the number of threads per block
// nTr is the number if points in the integral
// j is the index of the Integral point
// The kernel must be launched with only 1 block per bin
//
// Author: Simon Grimm
// January 2015
// *****************************************
template <int nb>
__global__ void Integrate_kernel(double *K_d, double *Tr_d, int NL, int nTr, double dTr, int *Nxmin_d, double kmin){

	int idy = threadIdx.x;
	int idx = blockIdx.x;
	__shared__ double a_s[nb];
	int Nxmin = Nxmin_d[idx];

	for(int j = 0; j < nTr; ++j){
		__syncthreads();
		a_s[idy] = 0.0;

		double m = exp((j - nTr/2) * dTr);

		__syncthreads();

		for(int k = 0; k < NL; k += nb){
			if(idy + k < NL){
				a_s[idy] += 1.0 * exp(-K_d[idy + k + idx * NL] * m);
			}
		}
		__syncthreads();

		if(nb >= 512){
			if(idy < 256){
				a_s[idy] += a_s[idy + 256];
			}
		}
		__syncthreads();

		if(nb >= 256){
			if(idy < 128){
				a_s[idy] += a_s[idy + 128];
			}
		}
		__syncthreads();

		if(nb >= 128){
			if(idy < 64){
				a_s[idy] += a_s[idy + 64];
			}
		}
		__syncthreads();
		if(idy < 3){
			//correct for Simpsons rule 
			a_s[idy] += wS_c[idy] * exp(-K_d[idy + idx * NL] * m);
			a_s[idy] += wS_c[2-idy] * exp(-K_d[NL - 1 - idy + idx * NL] * m);
		}
		__syncthreads();
		if(idy < 32){
			volatile double *a = a_s;
			if(nb >= 64) a[idy] += a[idy + 32];
			if(nb >= 32) a[idy] += a[idy + 16];
			if(nb >= 16) a[idy] += a[idy + 8];
			if(nb >= 8) a[idy] += a[idy + 4];
			if(nb >= 4) a[idy] += a[idy + 2];
			if(nb >= 2) a[idy] += a[idy + 1];
		}
		__syncthreads();

		if(idy == 0){
			Tr_d[idx * nTr + j] = a_s[0] / ((double)(NL - 1)) * (NL - Nxmin) / ((double)(NL)) + exp(-kmin * m) * Nxmin/ ((double)(NL));
	//		printf("%.20g %.20g\n", m, a_s[0] / (double)(NL));
		}
	}
}

// ****************************************
// This function set the factors for Simpsons rule and compies them
// to constant memory
//
// Author: Simon Grimm
// January 2015
// *****************************************
__host__ void SimpsonCoefficient(){
	//set correction factor for simpsons rule
	double *wS_h;
	wS_h = (double*)malloc(3 * sizeof(double));
	wS_h[0] = -5.0/8.0;
	wS_h[1] = 1.0/6.0;
	wS_h[2] = -1.0/24.0;
	cudaMemcpyToSymbol(wS_c, wS_h, 3 * sizeof(double), 0, cudaMemcpyHostToDevice);
	free(wS_h);
}

// ****************************************
// This kernel computes the terms needed for the Planck and Rosseland means. 
// It computes also the denominators. For the Planck mean it would be
// int_0^infty (2 h nu^3 / c^2 /(exp(hv/(kBT)) - 1) dnu = 2 kB^4 T^4 / ( h^3 c^2) pi^4/15
// but here we compute it also numerically to estimate the error.
//
// Author: Simon Grimm
// June 2015
// *****************************************
__global__ void Mean_kernel(double *x_d, double *Pmn_d, double *Rmn_d, double T, int Nx){

	int idy = threadIdx.x;
	int id = blockIdx.x * blockDim.x + idy;

	if(id < Nx){
		double nu = x_d[id];
		double nu3 = nu * nu * nu;

		double t1 = 2.0 * def_h * nu3 * def_c * def_c;
		double t2 = def_h * nu * def_c / (def_kB * T);
		
		double e = exp(t2);
		double e1 = e - 1.0;

		double B = t1 / e1;
		//double dB_dT = t1 * t2 * e / (T * e1 * e1);
		//for low T, e and e1 are inf, therefore dB_dt = -nan
		// e/ (T e1 e1) = e / (T (e-1)(e-1)) = e / (T (e^2 -2e + 1))
		// = 1.0 / (T (e - 2 + 1/e))
		double dB_dT = t1 * t2 / (T * (e - 2.0 + 1.0/e));


		Pmn_d[id] = B;
		Rmn_d[id] = dB_dT;
//printf("%d %g %g | %g %g %g\n", id, Pmn_d[id], Rmn_d[id], t1, t2, e);

		if(nu == 0.0){
			Pmn_d[id] = 0.0;
			Rmn_d[id] = 0.0;
		}
	}
}

template <int nb>
__global__ void IntegrateMean_kernel(double *Pmn_d, double *Rmn_d, double *x_d, double *K_d, double *means_d, int NL, int useIndividualX){

	int idy = threadIdx.x;
	int idx = blockIdx.x;
	__shared__ double a_s[nb];

	a_s[idy] = 0.0;

	__syncthreads();

	for(int k = 0; k < NL; k += nb){
		if(useIndividualX == 0){
			if(idy + k < NL){
				if(idx == 0) a_s[idy] += (Pmn_d[idy + k] * K_d[idy + k]);
				if(idx == 1){
					if(K_d[idy + k] != 0.0) a_s[idy] += (Rmn_d[idy + k] / K_d[idy + k]);
				}
				if(idx == 2) a_s[idy] += Pmn_d[idy + k];
				if(idx == 3) a_s[idy] += Rmn_d[idy + k];
			}
		}
		else{
			//Trapezoid rule for unequal spaced x
			if(idy + k + 1 < NL){
				if(idx == 0) a_s[idy] += 0.5 * ((Pmn_d[idy + k] * K_d[idy + k]) + (Pmn_d[idy + k + 1] * K_d[idy + k + 1])) * (x_d[idy + k + 1] - x_d[idy + k]);
				if(idx == 1){
					if(K_d[idy + k] != 0.0 && K_d[idy + k + 1] != 0.0) a_s[idy] += 0.5 * ((Rmn_d[idy + k] / K_d[idy + k]) + (Rmn_d[idy + k + 1] / K_d[idy + k + 1])) * (x_d[idy + k + 1] - x_d[idy + k]);
				}
				if(idx == 2) a_s[idy] += 0.5 * (Pmn_d[idy + k] + Pmn_d[idy + k + 1]) * (x_d[idy + k + 1] - x_d[idy + k]);
				if(idx == 3) a_s[idy] += 0.5 * (Rmn_d[idy + k] + Rmn_d[idy + k + 1]) * (x_d[idy + k + 1] - x_d[idy + k]);
			}
		}
	}
	__syncthreads();

	if(nb >= 512){
		if(idy < 256){
			a_s[idy] += a_s[idy + 256];
		}
	}
	__syncthreads();

	if(nb >= 256){
		if(idy < 128){
			a_s[idy] += a_s[idy + 128];
		}
	}
	__syncthreads();

	if(nb >= 128){
		if(idy < 64){
			a_s[idy] += a_s[idy + 64];
		}
	}
	__syncthreads();
	if(idy < 3 && useIndividualX == 0){
		//correct for Simpsons rule 
		if(idx == 0){
			a_s[idy] += wS_c[idy] * (Pmn_d[idy] * K_d[idy]);
			a_s[idy] += wS_c[2-idy] * (Pmn_d[NL - 1 - idy] * K_d[NL - 1 - idy]);
		}
		if(idx == 1){
			if(K_d[idy] != 0.0) a_s[idy] += wS_c[idy] * (Rmn_d[idy] / K_d[idy]);
			if(K_d[NL - 1 - idy] != 0.0) a_s[idy] += wS_c[2-idy] * (Rmn_d[NL - 1 - idy] / K_d[NL - 1 - idy]);
		}
		if(idx == 2){
			a_s[idy] += wS_c[idy] * Pmn_d[idy];
			a_s[idy] += wS_c[2-idy] * Pmn_d[NL - 1 - idy];
		}
		if(idx == 3){
			a_s[idy] += wS_c[idy] * Rmn_d[idy];
			a_s[idy] += wS_c[2-idy] * Rmn_d[NL - 1 - idy];
		}
	}
	__syncthreads();
	if(idy < 32){
		volatile double *a = a_s;
		if(nb >= 64) a[idy] += a[idy + 32];
		if(nb >= 32) a[idy] += a[idy + 16];
		if(nb >= 16) a[idy] += a[idy + 8];
		if(nb >= 8) a[idy] += a[idy + 4];
		if(nb >= 4) a[idy] += a[idy + 2];
		if(nb >= 2) a[idy] += a[idy + 1];
	}
	__syncthreads();

	if(idy == 0){
		if(idx == 0){
			means_d[0] = a_s[0];
		}
		if(idx == 1){
			means_d[1] = a_s[0];
		}
		if(idx == 2){
			means_d[2] = a_s[0];
		}
		if(idx == 3){
			means_d[3] = a_s[0];
		}
//printf("%d %.20g\n", idx, a_s[0]);
	}
}
