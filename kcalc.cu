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
	er = read_parameters(param, paramFilename, argc, argv);
	if(er == 0){
		return 0;
	}
	if(param.dev > devCount || param.dev < 0){
		printf("Error: Devive Number is not allowed\n");
		return 0;
	}
	double time[4];

	FILE *InfoFile;
	char InfoFilename[160];
	sprintf(InfoFilename, "Info_%s.dat", param.name);
	InfoFile = fopen(InfoFilename, "w");

	cudaSetDevice(param.dev);
	for(int i = 0; i < 2; ++i){
		FILE *infofile;
		if(i == 0) infofile = InfoFile;
		if(i == 1) infofile = stdout;
		fprintf(infofile, "Version: %g\n", VERSION);
		fprintf(infofile, "name = %s\nT = %g\nP = %g\nMolecule = %d\nnumin = %g\nnumax = %g\ndnu = %g\ncutMode = %d\ncut = %g\ndoResampling = %d\nnC = %d\ndoTransmission = %d\nnTr = %d\ndTr =  %g\ndoStoreFullK = %d\ndostoreK = %d\n", 
			param.name, param.T, param.P, param.nMolecule, param.numin, param.numax, param.dnu, param.cutMode, param.cut, param.doResampling, param.nC, param.doTransmission, param.nTr, param.dTr, param.doStoreFullK, param.doStoreK);
		fprintf(infofile, "Profile = %d\n", PROFILE);
		fprintf(infofile, "Using device %d\n", param.dev);
		for(int i = 0; i < param.nbins - 1; ++i){
			fprintf(infofile, "bin %d %g - %g\n", i, param.bins[i], param.bins[i + 1]);
			if(param.bins[i + 1] - param.bins[i] < param.nC){
				fprintf(infofile, "Error: bin is smaller than nC -> Memory conflict\n");
				return 0;
			}
		}
	}
	fclose(InfoFile);

	int Nx = (int)((param.numax - param.numin) / param.dnu);

	//Compute partition function
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

	Line L;



	double *K_h, *K_d;
	double *Tr_h, *Tr_d;		//Transmission function
	double *V_d;			//Vandermonde like matrix for least sqaures
	double *C_d, *D_d;

	K_h = (double*)malloc(Nx * sizeof(double));
	cudaMalloc((void **) &K_d, Nx * sizeof(double));

	Tr_h = (double*)malloc(param.nbins * param.nTr * sizeof(double));
	cudaMalloc((void **) &Tr_d, param.nbins * param.nTr * sizeof(double));


int nl_b = 10000;	//number of lines per bin


	cudaMalloc((void **) &V_d, param.nC * nl_b * sizeof(double));
	cudaMalloc((void **) &C_d, param.nC * sizeof(double));
	cudaMalloc((void **) &D_d, param.nC * sizeof(double));

	//set correction factor for simpsons rule
	SimpsonCoefficient();
	
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	printf("error = %d = %s\n",error, cudaGetErrorString(error));


	//Allocate the memory for the Line properties
	Alloc_Line(L, m);
	//Read the Line list	
	er = readFile(m, part, L);
	if(er == 0){
		return 0;
	}
	printf("Number of lines: %d\n", m.NL);
	printf("Number of points: %d\n", Nx);

	//Copy Line data to the device
	Copy_Line(L, m);

	cudaMemset(K_d, 0, Nx * sizeof(double));

	cudaDeviceSynchronize();
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[0] = times + timems/1000000.0;
	printf("Time before S_kernel:    %g seconds\n", time[0]);

	cudaDeviceSynchronize();
	gettimeofday(&tt1, NULL);
	times = 0.0;
	timems = 0.0;

	for(int k = 0; k < m.NL; k += nthmax){
		int Nk = min(nthmax, m.NL);
		S_kernel <<< (Nk + 127) / 128, 128 >>> (L.nu_d, L.S_d, L.A_d, L.EL_d, L.alphaL_d, L.alphaD_d, L.n_d, L.mass_d, L.delta_d, L.Q_d, L.ID_d, m.NL, param.T, param.P, k);
	}	


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

	const int ntL = 64;	//number of threads in Line kernel
	
	int nLimits = (Nx + ntL - 1) / ntL;
	
	int2 *Limits_h, *Limits_d;
	Limits_h = (int2*)malloc(nLimits * sizeof(int2));
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


//cudaMemcpy(Limits_h, Limits_d, nLimits * sizeof(int2), cudaMemcpyDeviceToHost);
//for(int i = 0; i < nLimits; ++i){
//	printf("%d %d %d\n", i, Limits_h[i].x, Limits_h[i].y);
//}
	cudaDeviceSynchronize();
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[1] = times + timems/1000000.0;
	printf("Time before Line_kernel: %g seconds\n", time[1]);

	cudaDeviceSynchronize();
	gettimeofday(&tt1, NULL);
	times = 0.0;
	timems = 0.0;

//cudaDeviceSynchronize();

//return 0;

//	for(int k = 0; k < Nx; k += nthmax){
//		int Nk = min(nthmax, Nx);
//		for(int i = 0; i < m.NL; i += nlmax){
			//This loop reduces the running time of the kernel to a few seconds
			//A longer running time of a single kernel can cause a time out
//			Line_kernel < ntL, nlmax> <<< (Nk + ntL - 1) / ntL, ntL >>> (L.nu_d, L.S_d, L.alphaL_d, L.alphaD_d, K_d, param.dnu, param.numin, Nx, m.NL, param.cut, i, k);
//		}
//	}
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

	cudaDeviceSynchronize();
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[2] = times + timems/1000000.0;
	printf("Time for Line_kernel:    %g seconds\n", time[2]);

	gettimeofday(&tt1, NULL);
	times = 0.0;
	timems = 0.0;


	//Write the full line profile
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

	thrust::device_ptr<double> K_dt = thrust::device_pointer_cast(K_d);
	for(int i = 0; i < param.nbins - 1; ++i){
		//compute indexes of the bins
		int il = (int)((param.bins[i] - param.numin) / param.dnu);
		int ir = (int)((param.bins[i + 1] - param.numin) / param.dnu) - 1;
printf("%d %d %d\n", i, il, ir);

		thrust::sort(K_dt + il, K_dt + ir);

		// Do the resampling
		if(param.doResampling == 1){
			Vandermonde_kernel <<< (nl_b + 511) / 512, 512 >>> (V_d, (double)((ir - il)), param.nC);
			QR_kernel <512> <<< 1, 512 >>> (V_d, C_d, D_d, ir - il, param.nC);

			lnK_kernel <<< (ir - il + 511) / 512, 512 >>> (K_d + il, ir - il);
			leastSquare_kernel <512> <<< 1, 512 >>> (V_d, C_d, D_d, K_d + il, ir - il, param.nC);

			//Write the coefficients per bin to a file
			cudaMemcpy(K_h + il, K_d + il, param.nC * sizeof(double), cudaMemcpyDeviceToHost);
			FILE *Out3File;
			char Out3Filename[160];
			sprintf(Out3Filename, "Out_%s_cbin_%.5d.dat", param.name, i);
			Out3File = fopen(Out3Filename, "w");

			fprintf(Out3File, "%.20g ", 2.0 * K_h[il + i]);
			for(int i = 1; i < param.nC; ++i){
				fprintf(Out3File, "%.20g ", K_h[il + i]);
			}
			fclose(Out3File);

			if(param.doTransmission == 1 || param.doStoreK == 1){
				expfx_kernel <<< 1, 512 >>> (K_d +il, param.nC, ir - il);
			}	
		}
		//Write K per bin output
		if(param.doStoreK == 1){
			cudaMemcpy(K_h + il, K_d + il, (ir - il) * sizeof(double), cudaMemcpyDeviceToHost);
			FILE *Out2File;
			char Out2Filename[160];
			sprintf(Out2Filename, "Out_%s_bin_%.5d.dat", param.name, i);

			Out2File = fopen(Out2Filename, "w");
			for(int j = 0; j < (ir - il); ++j){
				fprintf(Out2File, "%g %.20g\n", j / (double)((ir - il)), K_h[j + il]);
			}
			fclose(Out2File);
		}

		//Calculate the Transmission function
		if(param.doTransmission == 1){
			FILE *Out3File;
			char Out3Filename[160];
			sprintf(Out3Filename, "Out_%s_tr_%.5d.dat", param.name, i);
			Out3File = fopen(Out3Filename, "w");

			for(int j = 0; j < param.nTr; ++j){
				double m = exp((j - param.nTr/2) * param.dTr);
				Integrate_kernel < 512 > <<< 1, 512 >>> (K_d + il, Tr_d + i * param.nTr + j, m, ir - il);
			}
			cudaMemcpy(Tr_h + i * param.nTr, Tr_d + i * param.nTr, param.nTr * sizeof(double), cudaMemcpyDeviceToHost);
			for(int j = 0; j < param.nTr; ++j){
				double m = exp((j - param.nTr/2) * param.dTr);
				fprintf(Out3File, "%.20g %.20g\n", m, Tr_h[i * param.nTr + j]);
			}
			fclose(Out3File);
		}
	}

	cudaDeviceSynchronize();
	gettimeofday(&tt2, NULL);
	times = (tt2.tv_sec - tt1.tv_sec);
	timems = (tt2.tv_usec - tt1.tv_usec);

	time[3] = times + timems/1000000.0;
	printf("Time after Line_kernel:  %g seconds\n", time[3]);

	InfoFile = fopen(InfoFilename, "a");
	fprintf(InfoFile,"Time before S_kernel:    %g seconds\n", time[0]);
	fprintf(InfoFile,"Time before Line_kernel: %g seconds\n", time[1]);
	fprintf(InfoFile,"Time for Line_kernel:    %g seconds\n", time[2]);
	fprintf(InfoFile,"Time after Line_kernel:  %g seconds\n", time[3]);
	fclose(InfoFile);	

	free_Line(L);

	error = cudaGetLastError();
	printf("error = %d = %s\n",error, cudaGetErrorString(error));

	return 0;
}
