#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "define.h"
#include "host.h"
#include "ISO.h"
#include "voigt.h"
#include "resample.h"


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


	char qFilename[160];
	char paramFilename[160];
	sprintf(qFilename, "%s", "q.dat");
	sprintf(paramFilename, "%s", "param.dat");

	//Read prameters
	Param param;
	param.dev = 0;
	er = read_parameters(param, paramFilename, argc, argv);
	if(er == 0){
		return 0;
	}
	if(param.dev >= devCount || param.dev < 0){
		printf("Error: Devive Number is not allowed\n");
		return 0;
	}
	double time[8];

	FILE *InfoFile;
	char InfoFilename[160];
	sprintf(InfoFilename, "Info_%s.dat", param.name);
	InfoFile = fopen(InfoFilename, "w");

	int runtimeVersion;
	int driverVersion;

	cudaRuntimeGetVersion(&runtimeVersion);
	cudaDriverGetVersion(&driverVersion);

	//Determine the number of points
	int Nx = (int)((param.numax - param.numin) / param.dnu);
	//Determine the numbers of points per bin
	int Nxb = Nx / param.nbins;


	cudaSetDevice(param.dev);
	cudaDeviceProp devProp;
	for(int i = 0; i < 2; ++i){
		FILE *infofile;
		if(i == 0) infofile = InfoFile;
		if(i == 1) infofile = stdout;

		for(int j = 0; j < devCount; ++j){
			cudaGetDeviceProperties(&devProp, j);
			fprintf(infofile,"Name:%s, Major:%d, Minor:%d, Max threads per Block:%d, Max x dim:%d, #Multiprocessors:%d, Clock Rate:%d, Memory Clock Rate:%d, Global Memory:%lu, Shared memory per block: %lu",
devProp.name, devProp.major, devProp.minor, devProp.maxThreadsPerBlock, devProp.maxThreadsDim[0], devProp.multiProcessorCount,  devProp.clockRate, devProp.memoryClockRate, devProp.totalGlobalMem, devProp.sharedMemPerBlock);

		}

		fprintf(infofile, "\nVersion: %g\n", VERSION);
		fprintf(infofile, "Using device %d\n\n", param.dev);
		fprintf(infofile, "Runtime Version %d\n", runtimeVersion);
		fprintf(infofile, "Driver Version %d\n", driverVersion);

		if(Nxb < param.nC && i == 0){
			printf("Number of points per bin smaller than the number of Chebyshev coefficients: Changed nC to %d\n", Nxb);
			fprintf(infofile, "Number of points per bin smaller than the number of Chebyshev coefficients: Changed nC to %d\n", Nxb);
			param.nC = Nxb;
		}

		fprintf(infofile, "name = %s\nT = %g\nP = %g\nMolecule = %d\nnumin = %g\nnumax = %g\ndnu = %g\ncutMode = %d\ncut = %g\ndoResampling = %d\nnC = %d\ndoTransmission = %d\nnTr = %d\ndTr =  %g\ndoStoreFullK = %d\ndostoreK = %d\nnbins = %d\nkmin = %g\nqalphaL = %g\n", 
			param.name, param.T, param.P, param.nMolecule, param.numin, param.numax, param.dnu, param.cutMode, param.cut, param.doResampling, param.nC, param.doTransmission, param.nTr, param.dTr, param.doStoreFullK, param.doStoreK, param.nbins, param.kmin, param.qalphaL);
		fprintf(infofile, "Profile = %d\n", PROFILE);
	}
	fclose(InfoFile);

	
	//Compute partition function
	Partition part;
	er = ChebCoeff(qFilename, part, param.T);
	if(er == 0){
		return 0;
	}

	//Allocate Molecule properties
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

	Line L;
	
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Initial error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;

	//Allocate the memory for the Line properties
	Alloc_Line(L, m);
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Line alloc error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	//**************************
	//Read the Line list	
	//**************************
	er = readFile(m, part, L, param.qalphaL);
	if(er == 0){
		return 0;
	}
	//Copy Line data to the device
	Copy_Line(L, m);
	//************************
	printf("Number of lines: %d\n", m.NL);
	printf("Number of points: %d\n", Nx);
	printf("Number of points per bin: %d\n", Nxb);

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Line Read error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[0] = times + timems/1000000.0;
	printf("Time for input:        %g seconds\n", time[0]);

	cudaDeviceSynchronize();
	gettimeofday(&tt1, NULL);

	//****************************
	//Compute Line properties
	//***************************
	for(int k = 0; k < m.NL; k += nthmax){
		int Nk = min(nthmax, m.NL);
		S_kernel <<< (Nk + 127) / 128, 128 >>> (L.nu_d, L.S_d, L.A_d, L.EL_d, L.alphaL_d, L.alphaD_d, L.n_d, L.mass_d, L.delta_d, L.Q_d, L.ID_d, m.NL, param.T, param.P, k);
	}	

/* *************
	//print number of lines per bin
	cudaMemcpy(L.nu_h, L.nu_d, m.NL * sizeof(double), cudaMemcpyDeviceToHost);
	int nLb[param.nbins];
	for(int i = 0; i < param.nbins; ++i){
		nLb[i] = 0;
	}
	double binWidth = (param.numax - param.numin) / ((double)(param.nbins));
	printf("%g\n", binWidth);
	for(int i = 0; i < m.NL; ++i){
		int b = int(L.nu_h[i] / binWidth);
		nLb[b] += 1;
	}
	for(int i = 0; i < param.nbins; ++i){
		printf("%d, ", nLb[i]);
	}
	printf("\n");
*/ 

	//Sort the data along nu
	thrust::device_ptr<double> nu_dt = thrust::device_pointer_cast(L.nu_d);
	thrust::device_ptr<int> ID_dt = thrust::device_pointer_cast(L.ID_d);

	thrust::sort_by_key(nu_dt, nu_dt + m.NL, ID_dt);

	//Destroy Q_d to sort S_d alphaL_d and alphaD_d
        int Nk = min(nthmax, m.NL);
        for(int k = 0; k < m.NL; k += nthmax){
		Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.S_d, L.Q_d, m.NL, k);
	}
        for(int k = 0; k < m.NL; k += nthmax){
		Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.S_d, L.ID_d, m.NL, k);
	}
        for(int k = 0; k < m.NL; k += nthmax){
		Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.alphaL_d, L.Q_d, m.NL, k);
	}
        for(int k = 0; k < m.NL; k += nthmax){
		Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.alphaL_d, L.ID_d, m.NL, k);
	}
        for(int k = 0; k < m.NL; k += nthmax){
		Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.alphaD_d, L.Q_d, m.NL, k);
	}
        for(int k = 0; k < m.NL; k += nthmax){
		Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.alphaD_d, L.ID_d, m.NL, k);
	}
	//********************************


	//********************************
	//Determine which lines the block in the Line kernel has to read
	//********************************
	const int ntL = 64;	//number of threads in Line kernel
	int nLimits = (Nx + ntL - 1) / ntL;
	
	int2 *Limits_d;
	cudaMalloc((void **) &Limits_d, nLimits * sizeof(int2));

	int *MaxLimits_h, *MaxLimits_d;
	MaxLimits_h = (int*)malloc(sizeof(int));
	cudaMalloc((void **) &MaxLimits_d, sizeof(int));
	cudaMemset(MaxLimits_d, 0, sizeof(int));

	setLimits_kernel <<< (nLimits + 255) / 256, 256 >>> (Limits_d, nLimits, m.NL, param.cut);
	if(param.cut != 0.0){
		Cutoff_kernel <<< (m.NL + 255) / 256 , 256 >>> (L.nu_d, L.ID_d, Limits_d, L.alphaL_d, L.alphaD_d, ntL, param.numin, param.dnu, m.NL, nLimits, param.cut, param.cutMode);
		MaxLimits_kernel <<< (nLimits + 255) / 256, 256 >>> (Limits_d, MaxLimits_d, nLimits, m.NL);
		cudaMemcpy(MaxLimits_h, MaxLimits_d, sizeof(int), cudaMemcpyDeviceToHost);
	}
	else MaxLimits_h[0] = m.NL;

/*
//	print Limits
	int2 *Limits_h;
	Limits_h = (int2*)malloc(nLimits * sizeof(int2));
	cudaMemcpy(Limits_h, Limits_d, nLimits * sizeof(int2), cudaMemcpyDeviceToHost);
	FILE *LimitsFile;
	char LimitsFilename[160];
	sprintf(LimitsFilename, "Limits_%s_dat", param.name);
	LimitsFile = fopen(LimitsFilename, "w");

	for(int i = 0; i < nLimits; ++i){
		fprintf(LimitsFile,"%d %d %d\n", i, Limits_h[i].x, Limits_h[i].y);
	}
	fclose(LimitsFile);
	free(Limits_h);
*/
	//*********************************************



	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Line error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[1] = times + timems/1000000.0;
	printf("Time for Lines:        %g seconds\n", time[1]);

	cudaDeviceSynchronize();
	gettimeofday(&tt1, NULL);


	//***********************************
	//Compute the opacity function K(x)
	//************************************
	double *K_h, *K_d;
	int *binKey_d;
	K_h = (double*)malloc(Nx * sizeof(double));
	cudaMalloc((void **) &K_d, Nx * sizeof(double));
	cudaMalloc((void **) &binKey_d, Nx * sizeof(int));
	InitialK_kernel <<< (Nx + 511) / 512, 512 >>> (K_d, Nx, param.kmin);
	binKey_kernel <<< (Nx + 511) / 512, 512 >>> (binKey_d, Nx, Nxb);


	double cut = param.cut;
	if(cut == 0.0) cut = 1.0e30;
	for(int k = 0; k < Nx; k += nthmax){
		Nk = min(nthmax, Nx - k);
		for(int i = 0; i < MaxLimits_h[0]; i += nlmax){
			int nl = min(MaxLimits_h[0] - i, nlmax);
			//This loop reduces the running time of the kernel to a few seconds
			//A longer running time of a single kernel can cause a time out
			Line_kernel < ntL > <<< (Nk + ntL - 1) / ntL, ntL >>> (L.nu_d, L.S_d, L.alphaL_d, L.alphaD_d, K_d, param.dnu, param.numin, Nx, m.NL, Limits_d, cut, param.cutMode, nl, i, k);
		}
	}
	//*************************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("K error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[2] = times + timems/1000000.0;
	printf("Time for K(x):         %g seconds\n", time[2]);

	gettimeofday(&tt1, NULL);

	//***************************
	//Write the full line profile
	//****************************
	if(param.doStoreFullK == 1){
		cudaMemcpy(K_h, K_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
		FILE *OutFile;
		char OutFilename[160];
		sprintf(OutFilename, "Out_%s.dat", param.name);
		
		OutFile = fopen(OutFilename, "w");
		for(int j = 0; j < Nx; ++j){
			double x = param.numin + j * param.dnu;
			fprintf(OutFile, "%.20g %.20g\n", x, K_h[j]);
		}
		fclose(OutFile);
	}
	//*******************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Write error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[3] = times + timems/1000000.0;
	printf("Time for write K(x):   %g seconds\n", time[3]);

	gettimeofday(&tt1, NULL);

	//***************************************
	//Do the sorting of K for all bins
	//**************************************
	thrust::device_ptr<double> K_dt = thrust::device_pointer_cast(K_d);
	thrust::device_ptr<int> binKey_dt = thrust::device_pointer_cast(binKey_d);
	thrust::sort_by_key(K_dt, K_dt + Nx, binKey_dt);
	thrust::stable_sort_by_key(binKey_dt, binKey_dt + Nx, K_dt);
	//****************************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Sort error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[4] = times + timems/1000000.0;
	printf("Time for sort K(x):    %g seconds\n", time[4]);

	gettimeofday(&tt1, NULL);

	//********************************
	//Prepare Resampling and do QR factorization, the same for all bins
	//*********************************
	int *Nxmin_h, *Nxmin_d;		
	Nxmin_h = (int*)malloc(param.nbins * sizeof(int));
	cudaMalloc((void **) &Nxmin_d, param.nbins * sizeof(int));
	for(int i = 0; i < param.nbins; ++i){
		Nxmin_h[i] = 0;
	}
	cudaMemset(Nxmin_d, 0, param.nbins * sizeof(int));
	if(param.doResampling == 1){

		double *K2_h, *K2_d;
		K2_h = (double*)malloc(Nx * sizeof(double));
		cudaMalloc((void **) &K2_d, Nx * sizeof(double));
		cudaMemset(K2_d, 0, Nx * sizeof(double));

		findCut_kernel <<< (Nx + 511) / 512, 512 >>> (K_d, Nx, Nxb, param.kmin, Nxmin_d, param.nbins);
		rescale_kernel < 512 > <<< param.nbins, 512 >>> (Nxmin_d, K_d, K2_d, Nxb, param.kmin);
/*
cudaMemcpy(K2_h, K2_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
cudaMemcpy(K_h, K_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();

for(int i = 0; i < param.nbins; ++i){
	int il = i * Nxb;
	if(K_h[il] == param.kmin){
		for(int j = 0; j < Nxb; ++j){
			printf("%g %.20g\n", j / (double)(Nxb), K2_h[j + il]);
		}
		printf("\n\n");
	}
}
*/
		copyK2_kernel< 512 > <<< param.nbins, 512 >>> (K_d, K2_d, param.kmin, Nxb);
		cudaMemcpy(Nxmin_h, Nxmin_d, param.nbins * sizeof(int), cudaMemcpyDeviceToHost);
	

		double *V_d;			//Vandermonde like matrix for least sqaures
		double *C_d, *D_d;

		cudaMalloc((void **) &V_d, param.nC * Nxb * sizeof(double));
		cudaMalloc((void **) &C_d, param.nC * sizeof(double));
		cudaMalloc((void **) &D_d, param.nC * sizeof(double));

		Vandermonde_kernel <<< (Nxb + 511) / 512, 512 >>> (V_d, (double)(Nxb), param.nC);
		QR_kernel <512> <<< 1, 512 >>> (V_d, C_d, D_d, Nxb, param.nC);

		lnK_kernel <<< (Nx + 511) / 512, 512 >>> (K_d, Nx);
		leastSquare_kernel <512> <<< param.nbins, 512 >>> (V_d, C_d, D_d, K_d, Nxb, param.nC);

		FILE *Out3File;
		char Out3Filename[160];
		sprintf(Out3Filename, "Out_%s_cbin.dat", param.name);
		Out3File = fopen(Out3Filename, "w");
		for(int i = 0; i < param.nbins; ++i){
			int il = i * Nxb;
			cudaMemcpy(K_h + il, K_d + il, param.nC * sizeof(double), cudaMemcpyDeviceToHost);
	
			fprintf(Out3File, "%.20g %.20g ", param.kmin, Nxmin_h[i] / ((double)(Nxb)));
			fprintf(Out3File, "%.20g ", 2.0 * K_h[il + i]);
			for(int i = 1; i < param.nC; ++i){
				fprintf(Out3File, "%.20g ", K_h[il + i]);
			}
			fprintf(Out3File, "\n\n");
		}
		fclose(Out3File);
		if(param.doTransmission == 1 || param.doStoreK == 1){
			expfx_kernel <<< param.nbins, 512 >>> (K_d, param.nC, Nxb);
		}	
		cudaFree(V_d);
		cudaFree(C_d);
		cudaFree(D_d);
		cudaFree(K2_d);
		free(K2_h);
	}
	//**********************************
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Resampling error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[5] = times + timems/1000000.0;
	printf("Time for Resampling:   %g seconds\n", time[5]);

	gettimeofday(&tt1, NULL);

	//*****************************
	//Write K per bin output
	//*****************************
	if(param.doStoreK == 1){
		cudaMemcpy(K_h, K_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
		FILE *Out2File;
		char Out2Filename[160];
		sprintf(Out2Filename, "Out_%s_bin.dat", param.name);
		Out2File = fopen(Out2Filename, "w");
		for(int i = 0; i < param.nbins; ++i){
			int il = i * Nxb;
			int Nxmin = Nxmin_h[i];
			for(int j = 0; j < Nxb; ++j){
				fprintf(Out2File, "%g %.20g\n", Nxmin / ((double)(Nxb)) + j / ((double)(Nxb)) * (Nxb - Nxmin) / ((double)(Nxb)), K_h[j + il]);
			}
			fprintf(Out2File,"\n\n");
		}
		fclose(Out2File);
	}
	//******************************
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Write error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[6] = times + timems/1000000.0;
	printf("Time for write K(y):   %g seconds\n", time[6]);

	gettimeofday(&tt1, NULL);


	//*********************************
	//Calculate the Transmission function
	//*********************************
	if(param.doTransmission == 1){

		double *Tr_h, *Tr_d;
		Tr_h = (double*)malloc(param.nbins * param.nTr * sizeof(double));
		cudaMalloc((void **) &Tr_d, param.nbins * param.nTr * sizeof(double));

		//set correction factor for simpsons rule needed for resampling
		SimpsonCoefficient();
		FILE *Out3File;
		char Out3Filename[160];

		sprintf(Out3Filename, "Out_%s_tr.dat", param.name);
		Out3File = fopen(Out3Filename, "w");

		Integrate_kernel < 512 > <<< param.nbins, 512 >>> (K_d, Tr_d, Nxb, param.nTr, param.dTr, Nxmin_d, param.kmin);
		cudaMemcpy(Tr_h, Tr_d, param.nbins * param.nTr * sizeof(double), cudaMemcpyDeviceToHost);
		for(int i = 0; i < param.nbins; ++i){
			for(int j = 0; j < param.nTr; ++j){
				double m = exp((j - param.nTr/2) * param.dTr);
				fprintf(Out3File, "%.20g %.20g\n", m, Tr_h[i * param.nTr + j]);
			}
			fprintf(Out3File, "\n\n");
		}
		fclose(Out3File);
		free(Tr_h);
		cudaFree(Tr_d);
	}
	//************************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("Transmission error = %d = %s\n",error, cudaGetErrorString(error));
	if(error != 0) return 0;
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[7] = times + timems/1000000.0;
	printf("Time for Transmission: %g seconds\n", time[7]);

	InfoFile = fopen(InfoFilename, "a");
	fprintf(InfoFile,"Number of lines: %d\n", m.NL);
	fprintf(InfoFile,"Number of points: %d\n", Nx);
	fprintf(InfoFile,"Number of points per bin: %d\n", Nxb);
	fprintf(InfoFile,"Time for input:        %g seconds\n", time[0]);
	fprintf(InfoFile,"Time for Lines:        %g seconds\n", time[1]);
	fprintf(InfoFile,"Time for K(x):         %g seconds\n", time[2]);
	fprintf(InfoFile,"Time for write K(x):   %g seconds\n", time[3]);
	fprintf(InfoFile,"Time for sort K(x):    %g seconds\n", time[4]);
	fprintf(InfoFile,"Time for Resampling:   %g seconds\n", time[5]);
	fprintf(InfoFile,"Time for write K(y):   %g seconds\n", time[6]);
	fprintf(InfoFile,"Time for Transmission: %g seconds\n", time[7]);
	fclose(InfoFile);	


	free_Line(L);
	free(MaxLimits_h);
	free(K_h);
	free(Nxmin_h);

	cudaFree(Limits_d);
	cudaFree(MaxLimits_d);
	cudaFree(K_d);
	cudaFree(binKey_d);
	cudaFree(Nxmin_d);	

	error = cudaGetLastError();
	printf("Final error = %d = %s\n",error, cudaGetErrorString(error));

	return 0;
}
