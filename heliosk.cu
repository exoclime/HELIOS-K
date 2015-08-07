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
	param.useIndividualBins = 0;
	param.useOutputEdges = 0;
	param.nedges = 0;
	param.nP = 1;
	param.usePFile = 0;
	er = read_parameters(param, paramFilename, argc, argv);
	if(er == 0){
		return 0;
	}
	if(param.dev >= devCount || param.dev < 0){
		printf("Error: Devive Number is not allowed\n");
		return 0;
	}

	//If the bin file is used store the boundaries of the bins
	double *individualBins_h, *individualBins_d;
	if(param.useIndividualBins == 1){
		individualBins_h = (double*)malloc((param.nbins + 1) * sizeof(double));
		cudaMalloc((void **) &individualBins_d, (param.nbins + 1) * sizeof(double));
		er = readBinFile(param, individualBins_h);
		if(er == 0) return 0;
		param.numin = individualBins_h[0];
		param.numax = individualBins_h[param.nbins];
		cudaMemcpy(individualBins_d, individualBins_h, (param.nbins + 1) * sizeof(double), cudaMemcpyHostToDevice);

		if(param.doResampling == 1){
			printf("Error: The resampling function is not supported for the bin-file option\n");
			return 0;
		}
		if(param.doTransmission == 1){
			printf("Error: The transmission function is not supported for the bin-file option\n");
			return 0;
		}
	}
	else{
		individualBins_h = NULL;
		individualBins_d = NULL;
	}
	//If the output edges file is used store the edges
	double *outputEdges_h;
	if(param.useOutputEdges == 1){
		outputEdges_h = (double*)malloc((param.nedges + 1) * sizeof(double));
		er = readEdgesFile(param, outputEdges_h);
		if(er == 0) return 0;
	}
	else{
		outputEdges_h = NULL;
	}
	//Allocate P array 
	double *P_h;
	P_h = (double*)malloc((param.nP) * sizeof(double));
	P_h[0] = param.P;
	if(param.usePFile == 1){
		er = readPFile(param, P_h);
		if(er == 0) return 0;
	}

	double time[9];
	for(int i = 0; i < 9; ++i){
		time[i] = 0.0;
	}

	FILE *InfoFile;
	char InfoFilename[160];
	sprintf(InfoFilename, "Info_%s.dat", param.name);
	InfoFile = fopen(InfoFilename, "w");

	int runtimeVersion;
	int driverVersion;

	cudaRuntimeGetVersion(&runtimeVersion);
	cudaDriverGetVersion(&driverVersion);

	//Determine the number of points
	int Nx = (int)((param.numax - param.numin) / param.dnu) + 1;
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
			fprintf(infofile,"Name:%s, Major:%d, Minor:%d, Max threads per Block:%d, Max x dim:%d\n, #Multiprocessors:%d, Clock Rate:%d, Memory Clock Rate:%d, Global Memory:%lu, Shared memory per block: %lu",
				devProp.name, devProp.major, devProp.minor, devProp.maxThreadsPerBlock, devProp.maxThreadsDim[0],
				devProp.multiProcessorCount,  devProp.clockRate, devProp.memoryClockRate, devProp.totalGlobalMem, devProp.sharedMemPerBlock);

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
		fprintf(infofile, "name = %s\n", param.name);
		fprintf(infofile, "T = %g\n", param.T);
		if(param.usePFile == 0){
			fprintf(infofile, "P = %g\n", P_h[0]);
		}
		else{
			fprintf(infofile, "P in file: %s\n", param.PFilename);
		}
		fprintf(infofile, "useHITEMP = %d\n", param.useHITEMP);
		fprintf(infofile, "Molecule = %d\n", param.nMolecule);
		fprintf(infofile, "pathToData = %s\n", param.path);
		fprintf(infofile, "numin = %g\n", param.numin);
		fprintf(infofile, "numax = %g\n", param.numax);
		fprintf(infofile, "dnu = %g\n", param.dnu);
		fprintf(infofile, "cutMode = %d\n", param.cutMode);
		fprintf(infofile, "cut = %g\n", param.cut);
		fprintf(infofile, "doResampling = %d\n", param.doResampling);
		fprintf(infofile, "nC = %d\n", param.nC);
		fprintf(infofile, "doTransmission = %d\n", param.doTransmission);
		fprintf(infofile, "nTr = %d\n", param.nTr);
		fprintf(infofile, "dTr =  %g\n", param.dTr);
		fprintf(infofile, "doStoreFullK = %d\n", param.doStoreFullK);
		fprintf(infofile, "dostoreK = %d\n", param.doStoreK);
		fprintf(infofile, "nbins = %d\n", param.nbins);
		fprintf(infofile, "kmin = %g\n", param.kmin);
		fprintf(infofile, "qalphaL = %g\n", param.qalphaL);
		fprintf(infofile, "doMean = %d\n", param.doMean);
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
	m.NL[0] = 0;
	m.id = param.nMolecule;	//1 = H2O, 2 = CO, 5 = CO, 6 = CH4
	m.nISO = 0;

	//Initialize the Isotopologue properties for ISO.h
	Init(m, param);

	if(m.NL[0] == 0){
		printf("Molecule Id is not allowed\n");
		return 0;
	}
	
	timeval tt1;			//start time
	timeval tt2;			//end time
	long long times, timems;	//elapsed time in seconds and microseconds

	cudaDeviceSynchronize();

	Line L;

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Initial error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}

	//Allocate the memory for the Line properties
	Alloc_Line(L, m);

	double *K_h, *K_d;
	int *binKey_d;
	K_h = (double*)malloc(Nx * sizeof(double));
	cudaMalloc((void **) &K_d, param.nP * Nx * sizeof(double));
	cudaMalloc((void **) &binKey_d, Nx * sizeof(int));
	InitialK_kernel <<< (param.nP * Nx + 511) / 512, 512 >>> (K_d, param.nP * Nx, param.kmin);
	binKey_kernel <<< (Nx + 511) / 512, 512 >>> (binKey_d, Nx, Nxb, individualBins_d, param.nbins, param.numin, param.dnu);
	//int *binKey_h; 	//only needed to check the key
	//binKey_h = (int*)malloc(Nx * sizeof(int));
	//cudaMemcpy(binKey_h, binKey_d, Nx * sizeof(int), cudaMemcpyDeviceToHost);
	//for(int i = 0; i < Nx; ++i){
	//	printf("%d %g %d\n", i, param.numin + i * param.dnu, binKey_h[i]);
	//}

	const int ntL = 64;	//number of threads in Line kernel
	int nLimits = (Nx + ntL - 1) / ntL;
	int2 *Limits_d;
	cudaMalloc((void **) &Limits_d, nLimits * sizeof(int2));
	int *MaxLimits_h, *MaxLimits_d;
	MaxLimits_h = (int*)malloc(sizeof(int));
	cudaMalloc((void **) &MaxLimits_d, sizeof(int));
	//int2 *Limits_h; 	//only needed for check the Limits
	//Limits_h = (int2*)malloc(nLimits * sizeof(int2));

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Line alloc error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}

	printf("Number of points: %d\n", Nx);
	printf("Number of points per bin: %d\n", Nxb);
	InfoFile = fopen(InfoFilename, "a");
	fprintf(InfoFile,"Number of points: %d\n", Nx);
	fprintf(InfoFile,"Number of points per bin: %d\n", Nxb);
	fclose(InfoFile);
	//**************************************
	//Starting the loop around the datafiles
	//*************************************
	for(int fi = 0; fi < m.nFiles; ++fi){

		time[1] = 0;
		time[2] = 0;

		printf("Reading file %d of %d\n", fi, m.nFiles);
		gettimeofday(&tt1, NULL);
		
		printf("Number of lines: %d\n", m.NL[fi]);

		//**************************
		//Read the Line list	
		//**************************
		er = readFile(m, part, L, param.qalphaL, fi);
		if(er == 0){
			return 0;
		}
		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != 0){
			printf("Line Read error = %d = %s\n",error, cudaGetErrorString(error));
			return 0;
		}

		gettimeofday(&tt2, NULL);
		times = (tt2.tv_sec - tt1.tv_sec);
		timems = (tt2.tv_usec - tt1.tv_usec);

		time[0] = times + timems/1000000.0;
		printf("Time for input:        %g seconds\n", time[0]);

		cudaDeviceSynchronize();
		gettimeofday(&tt1, NULL);

		//start the loop around the Pressure values. only 1 iteration if no Pressure file is given
		for(int iP = 0; iP < param.nP; ++iP){

			//Copy Line data to the device
			Copy_Line(L, m, fi);
			//************************

			//****************************
			//Compute Line properties
			//***************************
			for(int k = 0; k < m.NL[fi]; k += nthmax){
				int Nk = min(nthmax, m.NL[fi]);
				S_kernel <<< (Nk + 127) / 128, 128 >>> (L.nu_d, L.S_d, L.A_d, L.EL_d, L.alphaL_d, L.alphaD_d, L.n_d, L.mass_d, L.delta_d, L.Q_d, L.ID_d, m.NL[fi], param.T, P_h[iP], k);
			}	

/* // *************
			//uncoment this only when no Presure file is given
			//print number of lines per bin
			cudaMemcpy(L.nu_h, L.nu_d, m.NL[fi] * sizeof(double), cudaMemcpyDeviceToHost);
			int nLb[param.nbins];
			for(int i = 0; i < param.nbins; ++i){
				nLb[i] = 0;
			}
			double binWidth = (param.numax - param.numin) / ((double)(param.nbins));
			printf("%g\n", binWidth);
			for(int i = 0; i < m.NL[fi]; ++i){
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

			thrust::sort_by_key(nu_dt, nu_dt + m.NL[fi], ID_dt);

			//Destroy Q_d to sort S_d alphaL_d and alphaD_d
			int Nk = min(nthmax, m.NL[fi]);
			for(int k = 0; k < m.NL[fi]; k += nthmax){
				Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.S_d, L.Q_d, m.NL[fi], k);
			}
			for(int k = 0; k < m.NL[fi]; k += nthmax){
				Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.S_d, L.ID_d, m.NL[fi], k);
			}
			for(int k = 0; k < m.NL[fi]; k += nthmax){
				Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.alphaL_d, L.Q_d, m.NL[fi], k);
			}
			for(int k = 0; k < m.NL[fi]; k += nthmax){
				Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.alphaL_d, L.ID_d, m.NL[fi], k);
			}
			for(int k = 0; k < m.NL[fi]; k += nthmax){
				Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.alphaD_d, L.Q_d, m.NL[fi], k);
			}
			for(int k = 0; k < m.NL[fi]; k += nthmax){
				Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.alphaD_d, L.ID_d, m.NL[fi], k);
			}
			//********************************


			//********************************
			//Determine which lines the block in the Line kernel has to read
			//********************************
			cudaMemset(MaxLimits_d, 0, sizeof(int));

			setLimits_kernel <<< (nLimits + 255) / 256, 256 >>> (Limits_d, nLimits, m.NL[fi], param.cut);
			if(param.cut != 0.0){
				Cutoff_kernel <<< (m.NL[fi] + 255) / 256 , 256 >>> (L.nu_d, L.ID_d, Limits_d, L.alphaL_d, L.alphaD_d, ntL, param.numin, param.dnu, m.NL[fi], nLimits, param.cut, param.cutMode);
				MaxLimits_kernel <<< (nLimits + 255) / 256, 256 >>> (Limits_d, MaxLimits_d, nLimits, m.NL[fi]);
				cudaMemcpy(MaxLimits_h, MaxLimits_d, sizeof(int), cudaMemcpyDeviceToHost);
			}
			else MaxLimits_h[0] = m.NL[fi];

/*
			//print Limits
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
			if(error != 0){
				printf("Line error = %d = %s\n",error, cudaGetErrorString(error));
				return 0;
			}

			gettimeofday(&tt2, NULL);
			times = (tt2.tv_sec - tt1.tv_sec);
			timems = (tt2.tv_usec - tt1.tv_usec);

			time[1] += times + timems/1000000.0;

			if(iP == param.nP - 1){
				printf("Time for Lines:        %g seconds\n", time[1]);
			}
			cudaDeviceSynchronize();
			gettimeofday(&tt1, NULL);


			//***********************************
			//Compute the opacity function K(x)
			//************************************
			double cut = param.cut;
			if(cut == 0.0) cut = 1.0e30;
			for(int k = 0; k < Nx; k += nthmax){
				Nk = min(nthmax, Nx - k);
				for(int i = 0; i < MaxLimits_h[0]; i += nlmax){
					int nl = min(MaxLimits_h[0] - i, nlmax);
					//This loop reduces the running time of the kernel to a few seconds
					//A longer running time of a single kernel can cause a time out
					Line_kernel < ntL > <<< (Nk + ntL - 1) / ntL, ntL >>> (L.nu_d, L.S_d, L.alphaL_d, L.alphaD_d, K_d + iP * Nx, param.dnu, param.numin, Nx, m.NL[fi], Limits_d, cut, param.cutMode, nl, i, k);
				}
			}
			//*************************************

			cudaDeviceSynchronize();
			error = cudaGetLastError();
			if(error != 0){
				printf("K error = %d = %s\n",error, cudaGetErrorString(error));
				return 0;
			}
			gettimeofday(&tt2, NULL);
			times = (tt2.tv_sec - tt1.tv_sec);
			timems = (tt2.tv_usec - tt1.tv_usec);

			time[2] += times + timems/1000000.0;
			if(iP == param.nP - 1){
				printf("Time for K(x):         %g seconds\n", time[2]);
			}
			gettimeofday(&tt1, NULL);

		} // End of pressure loop

		InfoFile = fopen(InfoFilename, "a");
		fprintf(InfoFile,"Number of lines: %d\n", m.NL[fi]);
		fprintf(InfoFile,"Time for input:        %g seconds\n", time[0]);
		fprintf(InfoFile,"Time for Lines:        %g seconds\n", time[1]);
		fprintf(InfoFile,"Time for K(x):         %g seconds\n", time[2]);
		fclose(InfoFile);
	} // End of linefile loop

	//***************************
	//Write the full line profile
	//****************************
	if(param.doStoreFullK == 1){
		FILE *OutFile;
		char OutFilename[160];
		sprintf(OutFilename, "Out_%s.dat", param.name);
			
		OutFile = fopen(OutFilename, "w");
		for(int iP = 0; iP < param.nP; ++iP){
			cudaMemcpy(K_h, K_d + iP * Nx, Nx * sizeof(double), cudaMemcpyDeviceToHost);
			for(int j = 0; j < Nx; ++j){
				double x = param.numin + j * param.dnu;
				fprintf(OutFile, "%.20g %.20g\n", x, K_h[j]);
			}
			fprintf(OutFile, "\n\n");
		}
		fclose(OutFile);
	}
	//*******************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Write error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[3] = times + timems/1000000.0;
	printf("Time for write K(x):   %g seconds\n", time[3]);

	gettimeofday(&tt1, NULL);

	//**************************************
	//compute the Planck and Rosseland means
	//**************************************
	if(param.doMean == 1){
		
		double *Pm_d;
		double *Rm_d;
		double *Pmn_d;
		double *Rmn_d;

		cudaMalloc((void **) &Pm_d, Nx * sizeof(double));
		cudaMalloc((void **) &Rm_d, Nx * sizeof(double));
		cudaMalloc((void **) &Pmn_d, Nx * sizeof(double));
		cudaMalloc((void **) &Rmn_d, Nx * sizeof(double));
	
		double *means_h;	
		means_h = (double*)malloc(4 * sizeof(double));

		FILE *Out4File;
		char Out4Filename[160];

		sprintf(Out4Filename, "Out_%s_mean.dat", param.name);
		Out4File = fopen(Out4Filename, "w");
	
		for(int iP = 0; iP < param.nP; ++iP){

			Mean_kernel <<< (Nx + 511) / 512, 512 >>> (K_d + iP * Nx, Pm_d, Rm_d, Pmn_d, Rmn_d, param.T, Nx, param.numin, param.dnu);
/*
cudaMemcpy(K_h, Pm_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
for(int i = 0; i < Nx; ++i){
	printf("%g %g\n", param.numin + i * param.dnu, K_h[i]);
}
printf("\n\n");
cudaMemcpy(K_h, Rm_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
for(int i = 0; i < Nx; ++i){
	printf("%g %g\n", param.numin + i * param.dnu, K_h[i]);
}
printf("\n\n");
cudaMemcpy(K_h, Pmn_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
for(int i = 0; i < Nx; ++i){
	printf("%g %g\n", param.numin + i * param.dnu, K_h[i]);
}
printf("\n\n");
cudaMemcpy(K_h, Rmn_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
for(int i = 0; i < Nx; ++i){
	printf("%g %g\n", param.numin + i * param.dnu, K_h[i]);
}
printf("\n\n");
*/
			IntegrateMean_kernel <512> <<< 4, 512 >>> (Pm_d, Rm_d, Pmn_d, Rmn_d, Nx);
			double sigma = 2.0 * def_kB * def_kB * def_kB * def_kB / ( def_h * def_h * def_h * def_c * def_c * 15.0) * M_PI * M_PI * M_PI * M_PI * M_PI;
			double integral1 = sigma * param.T * param.T * param.T * param.T / M_PI;
			double integral2 = M_PI / (4.0 * sigma * param.T * param.T * param.T);
		
			cudaMemcpy(means_h + 0, Pm_d, sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(means_h + 1, Rm_d, sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(means_h + 2, Pmn_d, sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(means_h + 3, Rmn_d, sizeof(double), cudaMemcpyDeviceToHost);


			fprintf(Out4File, "%.20g\n%.20g\n%.20g\n%.20g\n%.20g\n%.20g\n", means_h[0] / means_h[2], means_h[3] / means_h[1],  means_h[2], integral1, means_h[3], 1.0 / integral2);
			fprintf(Out4File, "\n\n");
		}
		
		fclose(Out4File);

		free(means_h);
		cudaFree(Pm_d);
		cudaFree(Rm_d);
		cudaFree(Pmn_d);
		cudaFree(Rmn_d);
	}
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("maen error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[4] = times + timems/1000000.0;
	printf("Time for mean K(x):    %g seconds\n", time[4]);

	gettimeofday(&tt1, NULL);



	//***************************************
	//Do the sorting of K for all bins
	//**************************************
	thrust::device_ptr<double> K_dt = thrust::device_pointer_cast(K_d);
	thrust::device_ptr<int> binKey_dt = thrust::device_pointer_cast(binKey_d);
	for(int iP = 0; iP < param.nP; ++iP){
		thrust::sort_by_key(K_dt + iP * Nx, K_dt + Nx + iP * Nx, binKey_dt);
		thrust::stable_sort_by_key(binKey_dt, binKey_dt + Nx, K_dt + iP * Nx);
	}
	//****************************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Sort error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[5] = times + timems/1000000.0;
	printf("Time for sort K(x):    %g seconds\n", time[5]);

	gettimeofday(&tt1, NULL);

	//********************************
	//Prepare Resampling and do QR factorization, the same for all bins
	// this doesn't work with individual bins
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

		double *V_d;			//Vandermonde like matrix for least sqaures
		double *C_d, *D_d;

		cudaMalloc((void **) &V_d, param.nC * Nxb * sizeof(double));
		cudaMalloc((void **) &C_d, param.nC * sizeof(double));
		cudaMalloc((void **) &D_d, param.nC * sizeof(double));

		Vandermonde_kernel <<< (Nxb + 511) / 512, 512 >>> (V_d, (double)(Nxb), param.nC);
		QR_kernel <512> <<< 1, 512 >>> (V_d, C_d, D_d, Nxb, param.nC);

		FILE *Out3File;
		char Out3Filename[160];
		sprintf(Out3Filename, "Out_%s_cbin.dat", param.name);
		Out3File = fopen(Out3Filename, "w");

		for(int iP = 0; iP < param.nP; ++iP){	
			cudaMemset(K2_d, 0, Nx * sizeof(double));
			cudaMemset(Nxmin_d, 0, param.nbins * sizeof(int));

			findCut_kernel <<< (Nx + 511) / 512, 512 >>> (K_d + iP * Nx, Nx, Nxb, param.kmin, Nxmin_d, param.nbins);
			rescale_kernel < 512 > <<< param.nbins, 512 >>> (Nxmin_d, K_d + iP * Nx, K2_d, Nxb, param.kmin);
/*
cudaMemcpy(K2_h, K2_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
cudaMemcpy(K_h, K_d + iP * Nx, Nx * sizeof(double), cudaMemcpyDeviceToHost);
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
			copyK2_kernel< 512 > <<< param.nbins, 512 >>> (K_d + iP * Nx, K2_d, param.kmin, Nxb);
			cudaMemcpy(Nxmin_h, Nxmin_d, param.nbins * sizeof(int), cudaMemcpyDeviceToHost);
	

			lnK_kernel <<< (Nx + 511) / 512, 512 >>> (K_d + iP * Nx, Nx);
			leastSquare_kernel <512> <<< param.nbins, 512 >>> (V_d, C_d, D_d, K_d + iP * Nx, Nxb, param.nC);

			for(int i = 0; i < param.nbins; ++i){
				int il = i * Nxb;
				cudaMemcpy(K_h + il, K_d + il + iP * Nx, param.nC * sizeof(double), cudaMemcpyDeviceToHost);
		
				fprintf(Out3File, "%.20g %.20g ", param.kmin, Nxmin_h[i] / ((double)(Nxb)));
				for(int i = 0; i < param.nC; ++i){
					fprintf(Out3File, "%.20g ", K_h[il + i]);
				}
				fprintf(Out3File, "\n\n");
			}
			fprintf(Out3File, "\n\n");
			if(param.doTransmission == 1 || param.doStoreK == 1){
				expfx_kernel <<< param.nbins, 512 >>> (K_d + iP * Nx, param.nC, Nxb);
			}	
		}
		fclose(Out3File);
		cudaFree(V_d);
		cudaFree(C_d);
		cudaFree(D_d);
		cudaFree(K2_d);
		free(K2_h);
	}
	//**********************************
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Resampling error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[6] = times + timems/1000000.0;
	printf("Time for Resampling:   %g seconds\n", time[6]);

	gettimeofday(&tt1, NULL);

	//*****************************
	//Write K per bin output
	//*****************************
	if(param.doStoreK == 1){
		FILE *Out2File;
		char Out2Filename[160];
		sprintf(Out2Filename, "Out_%s_bin.dat", param.name);
		Out2File = fopen(Out2Filename, "w");
		for(int iP = 0; iP < param.nP; ++iP){
			cudaMemcpy(K_h, K_d + iP * Nx, Nx * sizeof(double), cudaMemcpyDeviceToHost);
			if(param.useIndividualBins == 0){
				for(int i = 0; i < param.nbins; ++i){
					int il = i * Nxb;
					int Nxmin = Nxmin_h[i];
					int iedge = 0;
					for(int j = 0; j < Nxb; ++j){
						double y = Nxmin / ((double)(Nxb - 1)) + j / ((double)(Nxb - 1)) * (Nxb - Nxmin - 1) / ((double)(Nxb - 1));
						double y1 = Nxmin / ((double)(Nxb - 1)) + (j + 1) / ((double)(Nxb - 1)) * (Nxb - Nxmin - 1) / ((double)(Nxb - 1));
						if(param.useOutputEdges == 0){
							fprintf(Out2File, "%g %.20g\n", y, K_h[j + il]);
						}
						else{
							double edge = outputEdges_h[iedge];
							if(y <= edge && edge <= y1 && iedge < param.nedges){
								fprintf(Out2File, "%g %.20g\n", edge, (K_h[j + 1 + il] - K_h[j + il]) / ( y1 - y) * (edge - y) + K_h[j + il]);
								++iedge;
							}
						}
					}
					fprintf(Out2File,"\n\n");
				}
			}
			else{
				int ib = 0;
				int j = 0;
				int iedge = 0;
				for(int i = 0; i < Nx; ++i){
					double nu = param.numin + i * param.dnu;
					double nul = individualBins_h[ib];
					double nur = individualBins_h[ib + 1];
					
					int il = (int)((nul - param.numin) / param.dnu);
					int ir = (int)((nur - param.numin) / param.dnu);
					int Nxb = ir - il;
					if(ib == 0) ++Nxb; //take into account the first left boundary

					double y = j / ((double)(Nxb - 1));
					double y1 = (j + 1) / ((double)(Nxb - 1));

					if(param.useOutputEdges == 0){
						fprintf(Out2File, "%g %.20g\n", y, K_h[i]);
					}
					else{
						double edge = outputEdges_h[iedge];
						if(y <= edge && edge <= y1 && iedge < param.nedges){
							fprintf(Out2File, "%g %.20g\n", edge, (K_h[i + 1] - K_h[i]) / ( y1 - y) * (edge - y) + K_h[i]);
	//printf("%d %g %g %g %g %g %g\n", j, y, y1, edge, (K_h[i + 1] - K_h[i]) / ( y1 - y) * (edge - y) + K_h[i], K_h[i], K_h[i + 1]);
							++iedge;

						}
					}
					++j;

					if(nu > nur - param.dnu){
						++ib;
						j = 0;
						fprintf(Out2File,"\n\n");
						iedge = 0;
					}
				}
			}
			fprintf(Out2File,"\n\n");
		}
		fclose(Out2File);
	}
	//******************************
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Write error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[7] = times + timems/1000000.0;
	printf("Time for write K(y):   %g seconds\n", time[7]);

	gettimeofday(&tt1, NULL);

	//set correction factor for simpsons rule needed for resampling
	SimpsonCoefficient();

	//*********************************
	//Calculate the Transmission function
	//*********************************
	if(param.doTransmission == 1){

		double *Tr_h, *Tr_d;
		Tr_h = (double*)malloc(param.nbins * param.nTr * sizeof(double));
		cudaMalloc((void **) &Tr_d, param.nbins * param.nTr * sizeof(double));

		FILE *Out3File;
		char Out3Filename[160];

		sprintf(Out3Filename, "Out_%s_tr.dat", param.name);
		Out3File = fopen(Out3Filename, "w");

		for(int iP = 0; iP < param.nP; ++iP){
			Integrate_kernel < 512 > <<< param.nbins, 512 >>> (K_d + iP * Nx, Tr_d, Nxb, param.nTr, param.dTr, Nxmin_d, param.kmin);
			cudaMemcpy(Tr_h, Tr_d, param.nbins * param.nTr * sizeof(double), cudaMemcpyDeviceToHost);
			for(int i = 0; i < param.nbins; ++i){
				for(int j = 0; j < param.nTr; ++j){
					double m = exp((j - param.nTr/2) * param.dTr);
					fprintf(Out3File, "%.20g %.20g\n", m, Tr_h[i * param.nTr + j]);
				}
				fprintf(Out3File, "\n\n");
			}
			fprintf(Out3File, "\n\n");
		}
		fclose(Out3File);
		free(Tr_h);
		cudaFree(Tr_d);
	}


	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Transmission error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[8] = times + timems/1000000.0;
	printf("Time for Transmission: %g seconds\n", time[8]);

	InfoFile = fopen(InfoFilename, "a");
	fprintf(InfoFile,"Time for write K(x):   %g seconds\n", time[3]);
	fprintf(InfoFile,"Time for mean K(x):    %g seconds\n", time[4]);
	fprintf(InfoFile,"Time for sort K(x):    %g seconds\n", time[5]);
	fprintf(InfoFile,"Time for Resampling:   %g seconds\n", time[6]);
	fprintf(InfoFile,"Time for write K(y):   %g seconds\n", time[7]);
	fprintf(InfoFile,"Time for Transmission: %g seconds\n", time[8]);
	fclose(InfoFile);	


	free_Line(L);
	free(MaxLimits_h);
	free(K_h);
	free(Nxmin_h);
	free(individualBins_h);
	free(outputEdges_h);

	cudaFree(Limits_d);
	cudaFree(MaxLimits_d);
	cudaFree(K_d);
	cudaFree(binKey_d);
	cudaFree(Nxmin_d);
	cudaFree(individualBins_d);	

	error = cudaGetLastError();
	printf("Final error = %d = %s\n",error, cudaGetErrorString(error));

	return 0;
}
