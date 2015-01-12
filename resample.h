
__constant__ double wS_c[3];


// ****************************************
// This kernel computes the Vandermonde like matrix for the least squares
// using Chebyshev polynomials
//
// Author: Simon Grimm
// December 2014
// *****************************************
__global__ void Vandermonde_kernel(double *A_d, double nl, int NC){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	int NL = (int)(nl);

	if(id < NL){

		double x = id * 2.0 / nl - 1.0;

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
// This kernel perform a QR decomposition if the Matrix V_d
// it overwrites the upper right part of it with R except of
// the diagonal part, which is written in D_d. C_d contains 
// the Householder scalar c. The rest of V_d is filled with the 
// Householder vectors.
//
// NL is the bigger dimension of V_d, NC the smaller
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
		scale = a_s[0];
		__syncthreads();
		a_s[idy] = 0.0;
		__syncthreads();
		for(int k = 0; k < NL; k += nb){
			if(idy + k < NL && idy + k >= i){
				double V = V_d[(idy + k) + NL * i] / scale;
				V_d[(idy + k) + NL * i] = V;
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

	__shared__ double a_s[nb];

	//Compute Qt b
	for(int j = 0; j < NC - 1; ++j){
		a_s[idy] = 0.0;
		__syncthreads();	
		for(int k = 0; k < NL; k += nb){
			if(idy + k < NL && idy + k >= j){
				a_s[idy] += V_d[(idy + k) + NL * j] * b_d[idy + k];
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
				b_d[idy + k] -= t * V_d[(idy + k) + NL * j];
			}
		}
		__syncthreads();
	}//end if j loop in Qt b

	//solve R x = Qt b

	if(idy == 0) b_d[NC - 1] /= D_d[NC - 1];
	__syncthreads();

	for(int i = NC - 2; i >= 0; --i){
		a_s[idy] = 0.0;
		__syncthreads();
		//Assume NC < nb
		if(idy < NC && idy >= i + 1){
			a_s[idy] += V_d[i + NL * idy] * b_d[idy];
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
			b_d[i] = (b_d[i] - a_s[0]) / D_d[i];

		}
	}//end of i loop in R x = Qtb
}


// ****************************************
// This kernel computes the log of an array K_d
// NL is the lenght of the array
//
// Author: Simon Grimm
// January 2015
// *****************************************
__global__ void lnK_kernel(double *K_d, int NL){

	int id = threadIdx.x + blockIdx.x * blockDim.x;

	if(id < NL){
		K_d[id] = log(K_d[id]);
	}
}

// ****************************************
// This kernel computes  exp(f(x)) using the Chebyshev coefficients
// calculated previosly and stored in b_d.
// The kernel must be called at least with NmaxSample threads, the number of
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
__global__ void expfx_kernel(double *b_d, int NC, int NL){

	int idy = threadIdx.x;
	__shared__ double b_s[NmaxSample];

	for(int k = 0; k < NmaxSample; k += blockDim.x){
		if(k + idy < NmaxSample){
			b_s[k + idy] = b_d[k + idy];
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
			b_d[idy + k] = exp(f);
		}
	}
}

// ****************************************
// This kernel integrates a Kurve exp(-K_d * m) using the extended Simpsons rule:
// I = 1/h (3/8 f1 + 7/6 f2 + 23/24/ f3 + f4 + ... + FN-3 + 23/23 fN-2 + 7/6 fN-1 + 3/8 fN) + O(1/N^4)
//
// The result is written in Tr_d[0]
//
// NL is the number of data points in K_d
// nb is the number of threads per block
// The kernel must be launched with only 1 block
//
// Author: Simon Grimm
// January 2015
// *****************************************
template <int nb>
__global__ void Integrate_kernel(double *K_d, double *Tr_d, double m, int NL){

	int idy = threadIdx.x;
	__shared__ double a_s[nb];

	a_s[idy] = 0.0;


	__syncthreads();

	for(int k = 0; k < NL; k += nb){
		if(idy + k < NL){
			a_s[idy] += 1.0 * exp(-K_d[idy + k] * m);
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
		a_s[idy] += wS_c[idy] * exp(-K_d[idy] * m);
		a_s[idy] += wS_c[2-idy] * exp(-K_d[NL - 1 - idy] * m);
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
		Tr_d[0] = a_s[0] / ((double)(NL - 1));
//		printf("%.20g %.20g\n", m, a_s[0] / (double)(NL));
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
