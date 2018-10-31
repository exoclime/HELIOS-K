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


/*
// runs with biliniar interpolation
// texDescr.filterMode = cudaFilterModeLinear;
__global__ void Voigt_texture_kernel(cudaTextureObject_t K2dtex, float *K_d, int Nx, int Ny, int Nxtex, int Nytex, size_t pitch){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;


	if(idx < Nx && idy < Ny){
		float x = idx * Nxtex / float(Nx);
		float y = idy * Nytex / float(Ny);
		//float x = idx / float(Nx);
		//float y = idy / float(Ny);
	
		float K = tex2D <float> (K2dtex, x + 0.5f , y + 0.5f);
		float *row = (float *)(((char *)K_d)+(idy*pitch));
    		row[idx] = K;
//if(idy == 0) printf("%d %d %f %f %f\n", idx, idy, x * 10.0f, y * 10.0f, K);

	}
}

// runs with manual biliniar interpolation
// texDescr.filterMode = cudaFilterModePoint;
__global__ void Voigt_textureb_kernel(cudaTextureObject_t K2dtex, float *K_d, int Nx, int Ny, int Nxtex, int Nytex, size_t pitch){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;


	if(idx < Nx && idy < Ny){
		float x = idx * Nxtex / float(Nx);
		float y = idy * Nytex / float(Ny);
	
		float K00 = tex2D <float> (K2dtex, x, y);
		float K10 = tex2D <float> (K2dtex, x + 1.0f, y);
		float K01 = tex2D <float> (K2dtex, x, y + 1.0f);
		float K11 = tex2D <float> (K2dtex, x + 1.0f, y + 1.0f);

		float xx = (idx % (Nx / Nxtex)) * Nxtex / float(Nx);	
		float yy = (idy % (Ny / Nytex)) * Nytex / float(Ny);	

		float K = (1.0f - xx) * ( 1.0f - yy) * K00 + xx * (1.0f - yy) * K10 + (1.0f - xx) * yy * K01 + xx * yy * K11;

		float *row = (float *)(((char *)K_d)+(idy*pitch));
    		row[idx] = K;
//if(idy == 0) printf("%d %d %f %f | %f %f | %f %f %f %f %f\n", idx, idy, x * 10.0f / Nx, y * 10.0f / Ny, xx, yy, K00, K10, K01, K11, K);

	}
}
// runs with manual biliniar interpolation
// texDescr.filterMode = cudaFilterModePoint;
__global__ void Voigt_b_kernel(float *K2d_d, float *K_d, int Nx, int Ny, int Nxtex, int Nytex, size_t pitch){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;


	if(idx < Nx && idy < Ny){
		int x = floor(idx * Nxtex / float(Nx));
		int y = floor(idy * Nytex / float(Ny));
		
		float *row1 = (float *)(((char *)K_d)+(y*pitch)) + x;
		float K00 = *row1;
		float *row2 = (float *)(((char *)K_d)+(y*pitch)) + x + 1;
		float K10 = *row2;
		float *row3 = (float *)(((char *)K_d)+((y + 1)*pitch)) + x;
		float K01 = *row3;
		float *row4 = (float *)(((char *)K_d)+((y + 1)*pitch)) + x + 1;
		float K11 = *row4;

		float xx = (idx % (Nx / Nxtex)) * Nxtex / float(Nx);	
		float yy = (idy % (Ny / Nytex)) * Nytex / float(Ny);	

		float K = (1.0f - xx) * ( 1.0f - yy) * K00 + xx * (1.0f - yy) * K10 + (1.0f - xx) * yy * K01 + xx * yy * K11;

		float *row = (float *)(((char *)K_d)+(idy*pitch));
    		row[idx] = K;
//if(idy == 0) printf("%d %d %f %f | %f %f | %f %f %f %f %f\n", idx, idy, x * 10.0f / Nx, y * 10.0f / Ny, xx, yy, K00, K10, K01, K11, K);

	}
}

//https://stackoverflow.com/questions/34622717/bicubic-interpolation-in-c
__device__ float cubic_hermite(float A, float B, float C, float D, float t){
	float a = -A / 2.0f + (3.0f * B) / 2.0f - (3.0f * C) / 2.0f + D / 2.0f;
	float b =  A - (5.0f * B) / 2.0f + 2.0f * C - D / 2.0f;
	float c = -A / 2.0f + C / 2.0f;
	float d = B;
	float tt = t * t;

	return a * t* tt + b * tt + c * t + d;
}

// runs with manual biliniar interpolation
// texDescr.filterMode = cudaFilterModePoint;
__global__ void Voigt_bicubic_kernel(cudaTextureObject_t K2dtex, float *K_d, int Nx, int Ny, int Nxtex, int Nytex, size_t pitch){
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if(idx > 0 && idy > 0 && idx < Nx - 1&& idy < Ny - 1){
		float x = idx * Nxtex / float(Nx);
		float y = idy * Nytex / float(Ny);
	
		float K00 = tex2D <float> (K2dtex, x - 1.0f, y - 1.0f);
		float K10 = tex2D <float> (K2dtex, x       , y - 1.0f);
		float K20 = tex2D <float> (K2dtex, x + 1.0f, y - 1.0f);
		float K30 = tex2D <float> (K2dtex, x + 2.0f, y - 1.0f);

		float K01 = tex2D <float> (K2dtex, x - 1.0f, y);
		float K11 = tex2D <float> (K2dtex, x       , y);
		float K21 = tex2D <float> (K2dtex, x + 1.0f, y);
		float K31 = tex2D <float> (K2dtex, x + 2.0f, y);

		float K02 = tex2D <float> (K2dtex, x - 1.0f, y + 1.0f);
		float K12 = tex2D <float> (K2dtex, x       , y + 1.0f);
		float K22 = tex2D <float> (K2dtex, x + 1.0f, y + 1.0f);
		float K32 = tex2D <float> (K2dtex, x + 2.0f, y + 1.0f);

		float K03 = tex2D <float> (K2dtex, x - 1.0f, y + 2.0f);
		float K13 = tex2D <float> (K2dtex, x       , y + 2.0f);
		float K23 = tex2D <float> (K2dtex, x + 1.0f, y + 2.0f);
		float K33 = tex2D <float> (K2dtex, x + 2.0f, y + 2.0f);

		float xx = (idx % (Nx / Nxtex)) * Nxtex / float(Nx);	
		float yy = (idy % (Ny / Nytex)) * Nytex / float(Ny);


		float K0 = cubic_hermite(K00, K10, K20, K30, xx);
		float K1 = cubic_hermite(K01, K11, K21, K31, xx);
		float K2 = cubic_hermite(K02, K12, K22, K32, xx);
		float K3 = cubic_hermite(K03, K13, K23, K33, xx);

	
		float K = cubic_hermite(K0, K1, K2, K3, yy);
if(idx == 15 && idy == 15) printf("%d %d %g %g %g %g %g %g %g\n", idx, idy, x, y, K00, K10, K20, K30, K0, K);

		float *row = (float *)(((char *)K_d)+(idy*pitch));
    		row[idx] = K;
//if(idy == 0) printf("%d %d %f %f | %f %f | %f %f %f %f %f\n", idx, idy, x * 10.0f / Nx, y * 10.0f / Ny, xx, yy, K00, K10, K01, K11, K);

	}
}
*/

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

/*
{

double xMax = 10.0;
double yMax = 10.0;

int Nx = 10000;
int Ny = 10000;

int Nxtex = Nx + 1;
int Nytex = Ny + 1;

int Nxtexf = Nx / 10 + 1;
int Nytexf = Ny / 10 + 1;


double *K2d_h, *K2d_d;
size_t pitch;
//with pitch, the 2d memory is extendend in one dimension to set memory alignment, pitch is the new Nxtex
K2d_h = (double*)malloc( Nxtex * Nytex * sizeof(double));
cudaMallocPitch((void **) &K2d_d, &pitch, Nxtex * sizeof(double), Nytex);
//printf("%d %d %lu\n", Nxtex, Nytex, pitch);

{
	double a = (double)(M_PI * sqrt(-1.0 / log(def_TOLF * 0.5)));
	double b = (double)(1.0 / sqrt(M_PI));
	double c = (double)(2.0 * a / M_PI);
	Voigt_2d_kernel <<< dim3((Nxtex + 31) / 32, (Nytex + 31) / 32), dim3(32, 32, 1) >>> (a, b, c, K2d_d, Nxtex, Nytex, pitch, xMax, xMax);
	cudaMemcpy2D(K2d_h, Nxtex * sizeof(double), K2d_d, pitch, Nxtex * sizeof(double), Nytex, cudaMemcpyDeviceToHost);
}
/ *
for(int i = 0; i < Nxtex - 1; ++i){
	for(int j = 0; j < Nytex - 1; ++j){
		double x = i * xMax / double(Nxtex);
		double y = j * yMax / double(Nytex);
		if( x < xMax && y < yMax){
			printf("%g %g %.15g\n", x, y, K2d_h[j * Nxtex + i]);
		}
	}
}
* /

float *K2df_h, *K2df_d;
size_t pitchf;
//with pitchf, the 2d memory is extendend in one dimension to set memory alignment, pitchf is the new Nxtexf
K2df_h = (float*)malloc( Nxtexf * Nytexf * sizeof(float));
cudaMallocPitch((void **) &K2df_d, &pitchf, Nxtexf * sizeof(float), Nytexf);
//printf("%d %d %lu\n", Nxtexf, Nytexf, pitchf);

{
	float a = (float)(M_PI * sqrt(-1.0f / log(def_TOLF * 0.5f)));
	float b = (float)(1.0f / sqrt(M_PI));
	float c = (float)(2.0f * a / M_PI);
	Voigt_2df_kernel <<< dim3((Nxtexf + 31) / 32, (Nytexf + 31) / 32), dim3(32, 32, 1) >>> (a, b, c, K2df_d, Nxtexf, Nytexf, pitchf, xMax, xMax);
	cudaMemcpy2D(K2df_h, Nxtexf * sizeof(float), K2df_d, pitchf, Nxtexf * sizeof(float), Nytexf, cudaMemcpyDeviceToHost);
}
/ *
for(int i = 0; i < Nxtexf - 1; ++i){
	for(int j = 0; j < Nytexf -1; ++j){
		float x = i * xMax / float(Nxtexf - 1);
		float y = j * yMax / float(Nytexf - 1);
		if( x < xMax && y < yMax){
			printf("%g %g %.15g\n", x, y, K2df_h[j * Nxtexf + i]);
		}
	}
}

return 0;
* /
//https://stackoverflow.com/questions/41749024/edit-cuda-texture-object
cudaTextureObject_t K2dtex;

cudaResourceDesc resDescr;
memset(&resDescr, 0, sizeof(cudaResourceDesc));
resDescr.resType = cudaResourceTypePitch2D;
resDescr.res.pitch2D.desc = cudaCreateChannelDesc<float>();
resDescr.res.pitch2D.devPtr = K2df_d;
resDescr.res.pitch2D.height = Nytexf;
resDescr.res.pitch2D.pitchInBytes = pitchf;
resDescr.res.pitch2D.width = Nxtexf;


cudaTextureDesc  texDescr;
memset(&texDescr, 0, sizeof(cudaTextureDesc));
texDescr.normalizedCoords = 0;
//texDescr.filterMode = cudaFilterModeLinear;
texDescr.filterMode = cudaFilterModePoint;
texDescr.addressMode[0] = cudaAddressModeClamp;
texDescr.addressMode[1] = cudaAddressModeClamp;
texDescr.addressMode[2] = cudaAddressModeClamp;
texDescr.readMode = cudaReadModeElementType;

cudaCreateTextureObject(&K2dtex, &resDescr, &texDescr, NULL);



float *K_h, *K_d;
K_h = (float*)malloc( Nx * Ny * sizeof(float));
//with pitch, the 2d memory is extendend in one dimension to set memory alignment, pitch is the new Nx
cudaMallocPitch((void **) &K_d, &pitch, Nx * sizeof(float), Ny);


for(int t = 0; t < 1; ++t){
	//Voigt_texture_kernel <<< dim3((Nx + 31) / 32, (Ny + 31) / 32), dim3(32, 32, 1) >>> (K2dtex, K_d, Nx, Ny, Nxtexf - 1, Nytexf - 1, pitch);
	//Voigt_textureb_kernel <<< dim3((Nx + 31) / 32, (Ny + 31) / 32), dim3(32, 32, 1) >>> (K2dtex, K_d, Nx, Ny, Nxtexf -1, Nytexf - 1, pitch);
	//Voigt_b_kernel <<< dim3((Nx + 31) / 32, (Ny + 31) / 32), dim3(32, 32, 1) >>> (K2d_d, K_d, Nx, Ny, Nxtex - 1, Nytex - 1, pitch);
	Voigt_bicubic_kernel <<< dim3((Nx + 31) / 32, (Ny + 31) / 32), dim3(32, 32, 1) >>> (K2dtex, K_d, Nx, Ny, Nxtexf - 1, Nytexf - 1, pitch);
}

cudaMemcpy2D(K_h, Nx * sizeof(float), K_d, pitch, Nx * sizeof(float), Ny, cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();

for(int i = 0; i < Nx; ++i){
	for(int j = 0; j < Ny; ++j){
		double x = i * xMax / double(Nx);
		double y = j * yMax / double(Ny);
		if( x < xMax && y < yMax){
			double diff = fabs(K2d_h[j * Nxtex + i] - K_h[j * Nx + i]);
			if(diff > 5.0e-7){
				printf("%g %g %.15g %.15g %.15g\n", x, y, K2d_h[j * Nxtex + i], K_h[j * Nx + i], diff);
			}
		}
	}
}
return 0;
}
*/
	char qFilename[15][160];	//for maximal 15 isotopologues
	char paramFilename[160];
	sprintf(paramFilename, "%s", "param.dat");

	//Read prameters
	Param param;
	param.dev = 0;
	param.useIndividualBins = 0;
	param.useOutputEdges = 0;
	param.nedges = 0;
	param.nP = 1;
	param.usePFile = 0;
	param.useIndividualX = 0;
	param.useCia = 0;

	param.T = 0.0;
	param.P = 0.0;
	param.useHITEMP = 0;
	param.nMolecule = 0;
	param.numin = 0.0;
	param.numax = 0.0;
	param.dnu = 0.0;
	param.Nxb = 0;
	param.cutMode = 0;
	param.cut = 0.0;
	param.doResampling = 0;
	param.nC = 0;
	param.doTransmission = 0;
	param.nTr = 0;
	param.dTr = 0.0;
	param.doStoreFullK = 0;
	param.doStoreK = 0;
	param.nbins = 0;
	param.kmin = 0.0;
	param.qalphaL = def_qALPHA_L;
	param.doMean = 0;
	param.units = 0;	
	param.replaceFiles = 0;
	param.RLOW = 0;
	param.profile = def_PROFILE;
	
	er = read_parameters(param, paramFilename, argc, argv);
	if(er == 0){
		return 0;
	}
	if(param.dev >= devCount || param.dev < 0){
		printf("Error: Device Number is not allowed\n");
		return 0;
	}
	if(param.useIndividualX == 1 && param.RLOW == 1){
		printf("Error: bins file and RLOW  not allowed\n");
		return 0;
	}

	char filemode[16];
	if(param.replaceFiles == 0){
		sprintf(filemode, "a");
	}
	else{
		sprintf(filemode, "w");
	}

	FILE *InfoFile;
	char InfoFilename[160];
	sprintf(InfoFilename, "Info_%s.dat", param.name);
	InfoFile = fopen(InfoFilename, filemode);

	int runtimeVersion;
	int driverVersion;

	cudaRuntimeGetVersion(&runtimeVersion);
	cudaDriverGetVersion(&driverVersion);

	cudaSetDevice(param.dev);
	cudaDeviceProp devProp;
	for(int i = 0; i < 2; ++i){
		FILE *infofile;
		if(i == 0) infofile = InfoFile;
		if(i == 1) infofile = stdout;

		for(int j = 0; j < devCount; ++j){
			cudaGetDeviceProperties(&devProp, j);
			fprintf(infofile,"Name:%s, Major:%d, Minor:%d, Max threads per Block:%d, Max x dim:%d\n, #Multiprocessors:%d, Clock Rate:%d, Memory Clock Rate:%d, Global Memory:%lu, Shared memory per block: %lu\n",
				devProp.name, devProp.major, devProp.minor, devProp.maxThreadsPerBlock, devProp.maxThreadsDim[0],
				devProp.multiProcessorCount,  devProp.clockRate, devProp.memoryClockRate, devProp.totalGlobalMem, devProp.sharedMemPerBlock);

		}
	}
	if(param.Nxb != 0){
		param.useIndividualX = 1;
	}

	//If the bin file is used, store the boundaries of the bins
	double *binBoundaries_h, *binBoundaries_d;
	binBoundaries_h = (double*)malloc((param.nbins + 1) * sizeof(double));
	cudaMalloc((void **) &binBoundaries_d, (param.nbins + 1) * sizeof(double));
	if(param.useIndividualBins == 1){
		er = readBinFile(param, binBoundaries_h);
		if(er == 0) return 0;
		param.numin = binBoundaries_h[0];
		param.numax = binBoundaries_h[param.nbins];

		if(param.doResampling > 0){
			printf("Error: The resampling function is not supported for the bin-file option\n");
			return 0;
		}
		if(param.doTransmission > 0){
			printf("Error: The transmission function is not supported for the bin-file option\n");
			return 0;
		}
	}
	else{
		for(int i = 0; i < param.nbins; ++i){
			binBoundaries_h[i] = param.numin + i * (param.numax - param.numin) / ((double)(param.nbins));
		}
		binBoundaries_h[param.nbins] = param.numax;
	}
	cudaMemcpy(binBoundaries_d, binBoundaries_h, (param.nbins + 1) * sizeof(double), cudaMemcpyHostToDevice);

//for(int i = 0; i < param.nbins + 1; ++i){
//	printf("binboundaries %d %g\n", i, binBoundaries_h[i]);
//}	

	int Nx;
	if(param.useIndividualX == 0){
		Nx = (int)((param.numax - param.numin) / param.dnu + 0.5); //+ 0.5 to round correctly between double and int
		if((param.numax - param.numin) / param.dnu + 0.5 >= 2147483647){
			printf("Error: Nx too large, integer overflow. %d %g\n", Nx, (param.numax - param.numin) / param.dnu);
			return 0;
		}
printf("%g %g %g %g\n", param.numax, param.numin, param.dnu, (param.numax - param.numin) / param.dnu + 0.5);
		param.Nxb = Nx / param.nbins;
		if(Nx % param.nbins != 0){
			printf("Error: range cannot be divided evenly in bins. %d %d %g\n", Nx, param.nbins,  Nx / ((double)(param.nbins)));
			return 0;
		}
	}
	else{
		Nx = param.nbins * param.Nxb;
		if(param.nbins * param.Nxb >= 2147483647){
			printf("Error: Nx too large, integer overflow. %d %g\n", Nx, (double)(param.nbins) * (double)(param.Nxb));
			return 0;
		}
		if(param.doResampling > 0){
			printf("Error: The resampling function is not supported for unequal spacing option\n");
			return 0;
		}
		if(param.doTransmission > 0){
			printf("Error: The transmission function is not supported for unequal spacing option\n");
			return 0;
		}
	}


	int Nx1 = (Nx + 9) / 10;

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

	//Allocate Molecule properties
	Molecule m;
	m.NL[0] = 0;
	m.id = param.nMolecule;	//1 = H2O, 2 = CO, 5 = CO, 6 = CH4
	m.nISO = 0;
	m.defaultL = 0.0;
	m.defaultn = 0.0;
	//Initialize the Isotopologue properties for ISO.h
	Init(m, param, qFilename);

	for(int i = 0; i < 2; ++i){
		FILE *infofile;
		if(i == 0) infofile = InfoFile;
		if(i == 1) infofile = stdout;
		fprintf(infofile, "\nVersion: %g\n", VERSION);
		fprintf(infofile, "Using device %d\n\n", param.dev);
		fprintf(infofile, "Runtime Version %d\n", runtimeVersion);
		fprintf(infofile, "Driver Version %d\n", driverVersion);

		if(param.Nxb < param.nC && i == 0){
			printf("Number of points per bin smaller than the number of Chebyshev coefficients: Changed nC to %d\n", param.Nxb);
			fprintf(infofile, "Number of points per bin smaller than the number of Chebyshev coefficients: Changed nC to %d\n", param.Nxb);
			param.nC = param.Nxb;
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
		fprintf(infofile, "cia System = %s\n", param.ciaSystem);
		fprintf(infofile, "pathToData = %s\n", param.path);
		fprintf(infofile, "numin = %g\n", param.numin);
		fprintf(infofile, "numax = %g\n", param.numax);
		fprintf(infofile, "dnu = %g\n", param.dnu);
		fprintf(infofile, "Nnu per bin = %d\n", param.Nxb);
		fprintf(infofile, "Number of points: %d\n", Nx);
		fprintf(infofile, "cutMode = %d\n", param.cutMode);
		fprintf(infofile, "cut = %g\n", param.cut);
		fprintf(infofile, "doResampling = %d\n", param.doResampling);
		fprintf(infofile, "nC = %d\n", param.nC);
		fprintf(infofile, "doTransmission = %d\n", param.doTransmission);
		fprintf(infofile, "nTr = %d\n", param.nTr);
		fprintf(infofile, "dTr =  %g\n", param.dTr);
		fprintf(infofile, "doStoreFullK = %d\n", param.doStoreFullK);
		fprintf(infofile, "pathToK = %s\n", param.pathK);
		fprintf(infofile, "dostoreK = %d\n", param.doStoreK);
		fprintf(infofile, "nbins = %d\n", param.nbins);
		if(param.useIndividualBins == 1){
			fprintf(infofile, "use Individual bins: %s\n", param.bins);
		}
		fprintf(infofile, "kmin = %g\n", param.kmin);
		fprintf(infofile, "qalphaL = %g\n", param.qalphaL);
		fprintf(infofile, "doMean = %d\n", param.doMean);
		fprintf(infofile, "Units = %d\n", param.units);
		fprintf(infofile, "Replace files = %d\n", param.replaceFiles);
		fprintf(infofile, "default L = %g\n", m.defaultL);
		fprintf(infofile, "default n = %g\n", m.defaultn);
		fprintf(infofile, "RLOW = %d\n", param.RLOW);
		fprintf(infofile, "profile = %d\n", param.profile);
		fprintf(infofile, "def_TOL = %g\n", def_TOL);
		fprintf(infofile, "def_TOLf = %g\n", def_TOLF);
		fprintf(infofile, "def_nthmax = %d\n", def_nthmax);
		fprintf(infofile, "def_nlmax = %d\n", def_nlmax);
		fprintf(infofile, "def_maxlines = %d\n", def_maxlines);
		fprintf(infofile, "def_maxfiles = %d\n", def_maxfiles);
		fprintf(infofile, "def_NmaxSample = %d\n", def_NmaxSample);
		fprintf(infofile, "def_NXLOW = %d\n", def_NXLOW);
		if(param.useOutputEdges == 1){
			fprintf(infofile, "use output edges: %s\n", param.edges);
		}

	}
	fclose(InfoFile);
	
	//Read partition function
	Partition part;
	er = readPartition(param, param.nMolecule, qFilename, part, param.T, m);
	if(er == 0){
		return 0;
	}

	if(param.useHITEMP == 2 && m.defaultL == 0.0){
		printf("Molecule Id is not allowed for ExoMol\n");
		return 0;
	}

	//compute the mean mass
	m.meanMass = 0.0;
	for(int i = 0; i < m.nISO; ++i){
		m.meanMass += m.ISO[i].Ab * m.ISO[i].m; //mean Molar Mass (g)
	}
//printf("mean mass %g\n", m.meanMass);
	double unitScale = 1.0;

	//Set cia System properties
	ciaSystem cia;
	if(param.useCia == 1){
		er = InitCia(m, cia, param);
		if(er == 0) return 0;
	}

	if(param.useCia == 1 && param.nMolecule != 0){
		printf("Error, not allowed to use a cia system with a molecule\n");
		return 0;
	}


	if(param.units == 1){
		unitScale = 1.0 / def_NA * m.meanMass;
		if(param.useCia == 1 && param.nMolecule == 0){
			unitScale = 1.0 / def_NA * cia.mass1;
		}
		param.kmin /= unitScale;
	}	

	timeval tt1;			//start time
	timeval tt2;			//end time
	long long times, timems;	//elapsed time in seconds and microseconds
	cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
	float milliseconds;

	cudaDeviceSynchronize();

	Line L;

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Initial error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}

	//Allocate memory for Line properties
	if(param.useHITEMP < 2){
		Alloc_Line(L, m);
	}
	else{
		Alloc2_Line(L, m);
	}
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Line alloc error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	double *K_h, *K_d;
	double *K1_d;
	double *Kc_d;
	double *x_h, *x_d;
	int *binKey_d;
	int *binIndex_h, *binIndex_d;
	K_h = (double*)malloc(Nx * sizeof(double));
	x_h = (double*)malloc(Nx * sizeof(double));
	binIndex_h = (int*)malloc((param.nbins + 2) * sizeof(int));
	cudaMalloc((void **) &K_d, param.nP * Nx * sizeof(double));
	if(param.RLOW == 1){
		cudaMalloc((void **) &K1_d, Nx1 * sizeof(double));
		cudaMalloc((void **) &Kc_d, Nx * sizeof(double));
	}
	else{
		K1_d = NULL;
		Kc_d = NULL;
	}
	cudaMalloc((void **) &x_d, Nx * sizeof(double));
	cudaMalloc((void **) &binKey_d, Nx * sizeof(int));
	cudaMalloc((void **) &binIndex_d, (param.nbins + 2) * sizeof(int));

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("K alloc error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	for(int k = 0; k < param.nP * Nx; k += def_nthmax){
		int Nk = min(def_nthmax, param.nP * Nx - k);
		InitialK_kernel <<< (Nk + 511) / 512, 512 >>> (K_d, param.nP * Nx, param.kmin, k);
	}
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("bin1 error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	for(int k = 0; k < Nx; k += def_nthmax){
		int Nk = min(def_nthmax, Nx - k);
		setX_kernel <<< (Nk + 511) / 512, 512 >>> (x_d, Nx, param.numin, param.dnu, param.Nxb, param.useIndividualX, binBoundaries_d, k);
	}
	cudaMemcpy(x_h, x_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
	for(int k = 0; k < Nx; k += def_nthmax){
		int Nk = min(def_nthmax, Nx - k);
		binKey_kernel <<< (Nk + 511) / 512, 512 >>> (binKey_d, Nx, param.Nxb, binBoundaries_d, param.nbins, param.numax, x_d, param.useIndividualX, k);
	}
	for(int k = 0; k < Nx; k += def_nthmax){
		int Nk = min(def_nthmax, Nx - k);
		binIndex_kernel <<< (Nk + 511) / 512, 512 >>> (binKey_d, binIndex_d, Nx, param.nbins, k);
	}
	cudaMemcpy(binIndex_h, binIndex_d, (param.nbins + 2) * sizeof(int), cudaMemcpyDeviceToHost);


	/*			
	int *binKey_h; 	//only needed to check the key
	binKey_h = (int*)malloc(Nx * sizeof(int));
	cudaMemcpy(binKey_h, binKey_d, Nx * sizeof(int), cudaMemcpyDeviceToHost);
	for(int i = 0; i < Nx; ++i){
		int bin = binKey_h[i];
		printf("%d %.10g %d %d %d\n", i, x_h[i], bin, binIndex_h[bin], binIndex_h[bin + 1]);
	}
	*/

	const int ntL = 64;	//number of threads in Line kernel
	int nLimits = (Nx + ntL - 1) / ntL;
	int2 *Limits_d;

	cudaMalloc((void **) &Limits_d, nLimits * sizeof(int2));
	int *MaxLimits_h, *MaxLimits_d;
	MaxLimits_h = (int*)malloc(sizeof(int));
	cudaMalloc((void **) &MaxLimits_d, sizeof(int));

	if(param.useCia == 1){
		for(int iP = 0; iP < param.nP; ++iP){
			readCiaFile(param, cia, x_h, K_h, Nx, param.T, P_h[iP], m.meanMass);
			cudaMemcpy(K_d + iP * Nx, K_h, Nx * sizeof(double), cudaMemcpyHostToDevice);
		}
	}

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Alloc error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}

	if(param.nMolecule > 0 && param.doStoreFullK >= 0){
		double *nuP;
		double *ialphaDP;
		double *vyP;
		nuP = NULL;
		ialphaDP = NULL;
		vyP = NULL;
		if(Nx > def_NXLOW){
			int n = min(def_maxlines, m.NLmax);
			nuP = (double*)malloc(n * sizeof(double));
			ialphaDP = (double*)malloc(n * sizeof(double));
			vyP = (double*)malloc(n * sizeof(double));
		}

		//**************************************
		//Starting the loop around the datafiles
		//**************************************
		int fi0 = m.nFiles;
		int fi1 = 0;


		if(param.cutMode == 0 && param.cut){
			for(int fi = 0; fi < m.nFiles; ++fi){
				if(m.fileLimit[fi] - param.cut <= param.numax) fi1 = fi + 1;
				else break;
			}
			for(int fi = m.nFiles - 1; fi >= 0; --fi){
				if(m.fileLimit[fi + 1] + param.cut >= param.numin) fi0 = fi;
				else break;
			}
		}
		else{
			fi0 = 0;
			fi1 = m.nFiles;
		}

		printf("File range %d to %d\n", fi0 + 1, fi1);

		time[0] = 0.0;
		time[1] = 0.0;
		time[2] = 0.0;
	
		int fi;
		for(fi = fi0; fi < fi1; ++fi){


			FILE *dataFile;
			char dataFilename[160];
			sprintf(dataFilename, "%sbin", m.dataFilename[fi]);
			dataFile  = fopen(dataFilename, "rb");

			if(dataFile == NULL){
				printf("Error: line list file not found %s\n", dataFilename);
				return 0;
			}

			printf("Reading file %d of %d\n", fi + 1, fi1);
			printf("Number of lines: %lld\n", m.NL[fi]);

			for(long long int iL = 0LL; iL < m.NL[fi]; iL += def_maxlines){
				int NL = min(def_maxlines, m.NL[fi] - iL);
				printf("Reading Line file part %d of %d with %d lines\n", (iL + def_maxlines - 1) / def_maxlines + 1, (m.NL[fi] + def_maxlines - 1) / def_maxlines, NL);

				double timeOld = time[0];
				if(iL == 0) time[0] = 0.0;
				gettimeofday(&tt1, NULL);
				//**************************
				//Read the Line list	
				//**************************
				if(param.useHITEMP < 2){
					er = readFile(param, m, part, L, param.qalphaL, NL, dataFile);
				}
				else{
					er = readFileExomol(param, m, part, L, NL, dataFile);
				}
				if(er == 0){
					return 0;
				}
				gettimeofday(&tt2, NULL);
				times = (tt2.tv_sec - tt1.tv_sec);
				timems = (tt2.tv_usec - tt1.tv_usec);
				time[0] += times + timems/1000000.0;
				
				printf("Time for input:        %g seconds\n", time[0]);

				if(iL > 0 || fi > fi0){
					//read data before synchronization
					if(m.NL[fi - 1] > 0){
						cudaEventSynchronize(stop);
						cudaEventElapsedTime(&milliseconds, start, stop);

						time[2] += milliseconds * 0.001;
					}
					printf("Time for K(x):         %g seconds\n", time[2]);
		
					cudaDeviceSynchronize();
					error = cudaGetLastError();
					if(error != 0){
						printf("Ka error = %d = %s\n",error, cudaGetErrorString(error));
						return 0;
					}
					if(iL == 0){
						InfoFile = fopen(InfoFilename, "a");
						fprintf(InfoFile,"File %d of %d\n", fi, fi1);
						fprintf(InfoFile,"Number of lines: %lld\n", m.NL[fi - 1]);
						fprintf(InfoFile,"Time for input:        %g seconds\n", timeOld);
						fprintf(InfoFile,"Time for Lines:        %g seconds\n", time[1]);
						fprintf(InfoFile,"Time for K(x):         %g seconds\n", time[2]);
						fclose(InfoFile);
						time[1] = 0.0;
						time[2] = 0.0;
					}
				}

				cudaEventRecord(start);

				//start the loop around the Pressure values. only 1 iteration if no Pressure file is given
				for(int iP = 0; iP < param.nP; ++iP){

					//Copy Line data to the device

					if(param.RLOW == 1){
						cudaMemset(K1_d, 0, Nx1 * sizeof(double));	
						cudaMemset(Kc_d, 0, Nx * sizeof(double));	
					}

					if(param.useHITEMP < 2){
						Copy_Line(L, m, NL);
					}
					else{
						Copy2_Line(L, m, NL);
					}
					//************************

					//***************************
					//Compute Line properties
					//***************************
					if(param.useHITEMP < 2){
						for(int k = 0; k < NL; k += def_nthmax){
							int Nk = min(def_nthmax, NL);
							if(Nk > 0) S2_kernel <<< (Nk + 127) / 128, 128 >>> (L.nu_d, L.S_d, L.Sf_d, L.A_d, L.vy_d, L.vyf_d, L.ialphaD_d, L.n_d, L.delta_d, L.EL_d, L.ID_d, L.va_d, L.vb_d, L.vcut2_d, L.S1_d, L.S1f_d, NL, param.numin, param.dnu, param.cut, param.cutMode, param.profile, param.useIndividualX, param.T, P_h[iP], k);
						}	
			/* // *************
						//uncoment this only when no Pressure file is given
						//print number of lines per bin
						cudaMemcpy(L.nu_h, L.nu_d, NL * sizeof(double), cudaMemcpyDeviceToHost);
						int nLb[param.nbins];
						for(int i = 0; i < param.nbins; ++i){
							nLb[i] = 0;
						}
						double binWidth = (param.numax - param.numin) / ((double)(param.nbins));
						printf("%g\n", binWidth);
						for(int i = 0; i < NL; ++i){
							int b = int(L.nu_h[i] / binWidth);
							nLb[b] += 1;
						}
						for(int i = 0; i < param.nbins; ++i){
							printf("%d, ", nLb[i]);
						}
						printf("\n");
			 
			*/
					}
					else{
						for(int k = 0; k < NL; k += def_nthmax){
							int Nk = min(def_nthmax, NL);
							if(Nk > 0) Sf_kernel <<< (Nk + 127) / 128, 128 >>> (L.nu_d, L.S_d, L.Sf_d, L.A_d, L.vy_d, L.vyf_d, L.ialphaD_d, L.n_d, L.EL_d, L.S1_d, L.S1f_d, L.va_d, L.vb_d, L.vcut2_d, NL, param.numin, param.dnu, param.cut, param.cutMode, param.profile, param.useIndividualX, param.T, P_h[iP], k);
						}
					}
					//Sort the data along nu
					thrust::device_ptr<double> nu_dt = thrust::device_pointer_cast(L.nu_d);
					thrust::device_ptr<int> ID_dt = thrust::device_pointer_cast(L.ID_d);

					thrust::sort_by_key(nu_dt, nu_dt + NL, ID_dt);

					//Destroy Q_d to sort S_d vy_d and ialphaD_d
					int Nk = min(def_nthmax, NL);
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.S_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.S_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copyf_kernel <<< (Nk + 127) / 128, 128 >>> (L.Sf_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sortf_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.Sf_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.vy_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.vy_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copyf_kernel <<< (Nk + 127) / 128, 128 >>> (L.vyf_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sortf_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.vyf_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.ialphaD_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.ialphaD_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copyf_kernel <<< (Nk + 127) / 128, 128 >>> (L.va_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sortf_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.va_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copyf_kernel <<< (Nk + 127) / 128, 128 >>> (L.vb_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sortf_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.vb_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copyf_kernel <<< (Nk + 127) / 128, 128 >>> (L.vcut2_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sortf_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.vcut2_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.S1_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.S1_d, L.ID_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Copyf_kernel <<< (Nk + 127) / 128, 128 >>> (L.S1f_d, L.Q_d, NL, k);
					}
					for(int k = 0; k < NL; k += def_nthmax){
						if(Nk > 0) Sortf_kernel <<< (Nk + 127) / 128, 128 >>> (L.Q_d, L.S1f_d, L.ID_d, NL, k);
					}
					if(Nx > def_NXLOW){
						cudaMemcpy(nuP, L.nu_d, NL * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(ialphaDP, L.ialphaD_d, NL * sizeof(double), cudaMemcpyDeviceToHost);
						cudaMemcpy(vyP, L.vy_d, NL * sizeof(double), cudaMemcpyDeviceToHost);
					}
					//********************************

					cudaDeviceSynchronize();
					error = cudaGetLastError();
					if(error != 0){
						printf("Sort error = %d = %s\n",error, cudaGetErrorString(error));
						return 0;
					}
					if(Nx <= def_NXLOW){
						//********************************
						//Determine which lines the block in the Line kernel has to read
						//********************************
						cudaMemset(MaxLimits_d, 0, sizeof(int));

						setLimits_kernel <<< (nLimits + 255) / 256, 256 >>> (Limits_d, nLimits, NL, param.cut);
						if(param.cut != 0.0){
							if(NL > 0) Cutoff_kernel <<< (NL + 255) / 256 , 256 >>> (L.nu_d, L.ID_d, Limits_d, L.vy_d, L.ialphaD_d, ntL, param.numin, param.dnu, NL, nLimits, param.cut, param.cutMode, Nx, x_d, param.useIndividualX);
							MaxLimits_kernel <<< (nLimits + 255) / 256, 256 >>> (Limits_d, MaxLimits_d, nLimits, NL);
							cudaMemcpy(MaxLimits_h, MaxLimits_d, sizeof(int), cudaMemcpyDeviceToHost);
						}
						else MaxLimits_h[0] = NL;
						printf("Maximum number of Line Blocks %d\n", MaxLimits_h[0]);
						
	/*					//print Limits
						int2 *Limits_h;
						Limits_h = (int2*)malloc(nLimits * sizeof(int2));
						cudaMemcpy(Limits_h, Limits_d, nLimits * sizeof(int2), cudaMemcpyDeviceToHost);
						FILE *LimitsFile;
						char LimitsFilename[160];
						sprintf(LimitsFilename, "Limits_%s.dat", param.name);
						if(fi == 0){
							LimitsFile = fopen(LimitsFilename, filemode);
						}
						else{
							LimitsFile = fopen(LimitsFilename, "a");
						}

						for(int i = 0; i < nLimits; ++i){
							fprintf(LimitsFile,"%d %d %d %d\n", fi, i, Limits_h[i].x, Limits_h[i].y);
						}
						fclose(LimitsFile);
						free(Limits_h);
	*/	 				
						//*********************************************
					}

					cudaEventRecord(stop);
					cudaEventSynchronize(stop);
					error = cudaGetLastError();
					if(error != 0){
						printf("Line error = %d = %s\n",error, cudaGetErrorString(error));
						return 0;
					}
					cudaEventElapsedTime(&milliseconds, start, stop);

					time[1] += milliseconds * 0.001;

					if(iP == param.nP - 1){
						printf("Time for Lines:        %g seconds\n", time[1]);
					}

					cudaEventRecord(start);

					//***********************************
					//Compute the opacity function K(x)
					//************************************
					double cut = param.cut;
					if(cut == 0.0) cut = 1.0e30;
					if(Nx <= def_NXLOW){
						float a = (float)(M_PI * sqrt(-1.0 / log(def_TOLF * 0.5)));
						float b = (float)(1.0 / sqrt(M_PI));
						float c = (float)(2.0 * a / M_PI);
						for(int k = 0; k < Nx; k += def_nthmax){
							int Nk = min(def_nthmax, Nx - k);
							for(int i = 0; i < MaxLimits_h[0]; i += def_nlmax){
if(i % (1000 * def_nlmax) == 0){
printf("%d %d %d %d\n", MaxLimits_h[0], def_nlmax, k, ntL);
}
								int nl = min(MaxLimits_h[0] - i, def_nlmax);
								//This loop reduces the running time of the kernel to a few seconds
								//A longer running time of a single kernel can cause a time out
								Line_kernel < ntL > <<< (Nk + ntL - 1) / ntL, ntL >>> (L.Sf_d, L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, K_d + iP * Nx, x_d, Nx, NL, Limits_d, nl, i, k, param.useIndividualX, param.Nxb, binBoundaries_d, a, b, c, param.profile);
							}
						}
					}
					else{ // end Nx <= def_NXLOW
						if(param.RLOW == 0){
							const int nl = 512;
							for(int il = 0; il < NL; il += nl){ //loop over lines
								int ii11 = 0;
								int ii00 = Nx;
								if(param.useIndividualX == 0){
									for(int iil = 0; iil < nl; ++iil){
										if(il + iil < NL){
											int Inu = (int)((nuP[il + iil] - param.numin) / param.dnu);
											int ii0 = Inu - (int)(cut / param.dnu) - 1;
											int ii1 = Inu + (int)(cut / param.dnu) + 2;

	//if(iil % 10000 == 0) printf("%d %.30g %d %d %d\n", il + iil, nuP[il + iil], Inu, ii0, ii1);

											ii11 = max(ii11, ii1);
											ii00 = min(ii00, ii0);
										}
									}
								}
								else{
									double nu00 = param.numax;
									double nu11 = param.numin;
									for(int iil = 0; iil < nl; ++iil){
										if(il + iil < NL){
											double nu0 = nuP[il + iil] - cut;
											double nu1 = nuP[il + iil] + cut;

											nu11 = fmax(nu11, nu1);
											nu00 = fmin(nu00, nu0);
										}
									}
									for(int bin = 0; bin < param.nbins; ++bin){
										if(binBoundaries_h[bin + 1] > nu11){
											double dnu = (binBoundaries_h[bin + 1] - binBoundaries_h[bin]) / ((double)(param.Nxb));
											int bstart = bin * param.Nxb;
											ii11 = (nu11 - binBoundaries_h[bin]) / dnu + bstart + 2;
											break;
										}
									}
									for(int bin = 0; bin < param.nbins; ++bin){
										if(binBoundaries_h[bin + 1] > nu00){
											double dnu = (binBoundaries_h[bin + 1] - binBoundaries_h[bin]) / ((double)(param.Nxb));
											int bstart = bin * param.Nxb;
											ii00 = (nu00 - binBoundaries_h[bin]) / dnu + bstart - 1;
											break;
										}
									}
								}

								ii11 = min(Nx, ii11);
								ii00 = max(0, ii00);

								int nt = ii11 - ii00;
								int nstart = ii00;
								int nll = min(nl, NL - il);	
	if(il % 10000 == 0) printf("A %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
								for(int k = 0; k < nt; k += def_nthmax){
									int Nk = min(def_nthmax, nt - k);
									if(Nk > 0 && nll > 0){
										Line2f_kernel < nl, 0 > <<< (max(Nk, nll) + nl - 1) / nl, nl >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, K_d + iP * Nx, il, nstart, Nk, nll, param.useIndividualX, param.Nxb, binBoundaries_d, 0.0f, 0.0f, 0.0f, param.profile);
									}
									nstart += def_nthmax;
								}
							}
						} //end  param.RLOW == 0
						else{
							//lower resolution
							const int nl = 512;
							for(int il = 0; il < NL; il += nl){ //loop over lines
								int ii11 = 0;
								int ii00 = Nx1;
								for(int iil = 0; iil < nl; ++iil){
									if(il + iil < NL){
										int Inu = (int)((nuP[il + iil] - param.numin) / (param.dnu * 10));
										int ii0 = Inu - (int)(cut / (param.dnu * 10)) - 1;
										int ii1 = Inu + (int)(cut / (param.dnu * 10)) + 2;
	//if(iil % 10000 == 0) printf("%d %.30g %d %d %d\n", il + iil, nuP[il + iil], Inu, ii0, ii1);

										ii11 = max(ii11, ii1);
										ii00 = min(ii00, ii0);
									}
								}
								ii11 = min(Nx1, ii11);
								ii00 = max(0, ii00);
								int nt = ii11 - ii00;
								int nstart = ii00;
								int nll = min(nl, NL - il);	
	if(il % 10000 == 0) printf("Ac %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
								for(int k = 0; k < nt; k += def_nthmax){
									int Nk = min(def_nthmax, nt - k);
									if(Nk > 0 && nll > 0){
										Line2f_kernel < nl, -1 > <<< (max(Nk, nll) + nl - 1) / nl, nl >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, K1_d, il, nstart, Nk, nll, param.useIndividualX, (param.Nxb + 9) / 10, binBoundaries_d, 0.0f, 0.0f, 0.0f, param.profile);
									}
									nstart += def_nthmax;
								}
							}
							//lower resolution interpolation correction
							const int nlb = 512;
							for(int il = 0; il < NL; il += nlb){ //loop over lines
								int ii11 = 0;
								int ii00 = Nx;
								for(int iil = 0; iil < nlb; ++iil){
									if(il + iil < NL){

										double aD2 = 1.0 / (ialphaDP[il + iil] * ialphaDP[il + iil]);
										double aL2 = vyP[il + iil] * vyP[il + iil] * aD2;
										double Dnu2 = 1.0e6 * aD2 - aL2;
										double Dnu = 0.0;
										if(Dnu2 > 0.0){
											Dnu = sqrt(Dnu2);
											int Inu = (int)((nuP[il + iil] - param.numin) / (param.dnu));
											int ii0 = ((Inu - (int)(Dnu / (param.dnu))) / 10) * 10;
											int ii1 = ii0 + 12;

											ii11 = max(ii11, ii1);
											ii00 = min(ii00, ii0);
										}
									}
								}
								ii11 = min(Nx, ii11);
								ii00 = max(0, ii00);
								int nt = ii11 - ii00;
								int nstart = ii00;
								int nll = min(nlb, NL - il);	
	if(il % 10000 == 0) printf("Bcl %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
								for(int k = 0; k < nt; k += def_nthmax){
									int Nk = min(def_nthmax, nt - k);
									if(Nk > 0 && nll > 0){
										Line2f_kernel < nlb, 10 > <<< (max(Nk, nll) + nlb - 1) / nlb, nlb >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, Kc_d, il, nstart, Nk, nll, param.useIndividualX, param.Nxb, binBoundaries_d, 0.0f, 0.0f, 0.0f, param.profile);
									}
									nstart += def_nthmax;
								}
							}
							for(int il = 0; il < NL; il += nlb){ //loop over lines
								int ii11 = 0;
								int ii00 = Nx;
								for(int iil = 0; iil < nlb; ++iil){
									if(il + iil < NL){

										double aD2 = 1.0 / (ialphaDP[il + iil] * ialphaDP[il + iil]);
										double aL2 = vyP[il + iil] * vyP[il + iil] * aD2;
										double Dnu2 = 1.0e6 * aD2 - aL2;
										double Dnu = 0.0;
										if(Dnu2 > 0.0){
											Dnu = sqrt(Dnu2);
											int Inu = (int)((nuP[il + iil] - param.numin) / (param.dnu));
											int ii0 = ((Inu + (int)(Dnu / (param.dnu))) / 10) * 10;
											int ii1 = ii0 + 12;

											ii11 = max(ii11, ii1);
											ii00 = min(ii00, ii0);
										}
									}
								}
								ii11 = min(Nx, ii11);
								ii00 = max(0, ii00);
								int nt = ii11 - ii00;
								int nstart = ii00;
								int nll = min(nlb, NL - il);	
	if(il % 10000 == 0) printf("Bcr %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
								for(int k = 0; k < nt; k += def_nthmax){
									int Nk = min(def_nthmax, nt - k);
									if(Nk > 0 && nll > 0){
										Line2f_kernel < nlb, 11 > <<< (max(Nk, nll) + nlb - 1) / nlb, nlb >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, Kc_d, il, nstart, Nk, nll, param.useIndividualX, param.Nxb, binBoundaries_d, 0.0f, 0.0f, 0.0f, param.profile);
									}
									nstart += def_nthmax;
								}
							}
							for(int il = 0; il < NL; il += nlb){ //loop over lines
								int ii11 = 0;
								int ii00 = Nx;
								for(int iil = 0; iil < nlb; ++iil){
									if(il + iil < NL){
										int Inu = (int)((nuP[il + iil] - param.numin) / param.dnu);
										int ii0 = (Inu + (int)(cut / param.dnu)) / 10 * 10;
										int ii1 = ii0 + 12;
	//if(iil % 10000 == 0) printf("%d %.30g %d %d %d\n", il + iil, nuP[il + iil], Inu, ii0, ii1);

										ii11 = max(ii11, ii1);
										ii00 = min(ii00, ii0);
									}
								}
								ii11 = min(Nx, ii11);
								ii00 = max(0, ii00);
								int nt = ii11 - ii00;
								int nstart = ii00;
								int nll = min(nlb, NL - il);	
	if(il % 10000 == 0) printf("Acr %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
								for(int k = 0; k < nt; k += def_nthmax){
									int Nk = min(def_nthmax, nt - k);
									if(Nk > 0 && nll > 0){
										Line2f_kernel < nlb, 12 > <<< (max(Nk, nll) + nlb - 1) / nlb, nlb >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, Kc_d, il, nstart, Nk, nll, param.useIndividualX, param.Nxb, binBoundaries_d, 0.0f, 0.0f, 0.0f, param.profile);
									}
									nstart += def_nthmax;
								}
							}
							for(int il = 0; il < NL; il += nlb){ //loop over lines
								int ii11 = 0;
								int ii00 = Nx;
								for(int iil = 0; iil < nlb; ++iil){
									if(il + iil < NL){
										int Inu = (int)((nuP[il + iil] - param.numin) / param.dnu);
										int ii0 = (Inu - (int)(cut / param.dnu)) / 10 * 10;
										int ii1 = ii0 + 12;
	//if(iil % 10000 == 0) printf("%d %.30g %d %d %d\n", il + iil, nuP[il + iil], Inu, ii0, ii1);

										ii11 = max(ii11, ii1);
										ii00 = min(ii00, ii0);
									}
								}
								ii11 = min(Nx, ii11);
								ii00 = max(0, ii00);
								int nt = ii11 - ii00;
								int nstart = ii00;
								int nll = min(nlb, NL - il);	
	if(il % 10000 == 0) printf("Acl %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
								for(int k = 0; k < nt; k += def_nthmax){
									int Nk = min(def_nthmax, nt - k);
									if(Nk > 0 && nll > 0){
										Line2f_kernel < nlb, 13 > <<< (max(Nk, nll) + nlb - 1) / nlb, nlb >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, Kc_d, il, nstart, Nk, nll, param.useIndividualX, param.Nxb, binBoundaries_d, 0.0f, 0.0f, 0.0f, param.profile);
									}
									nstart += def_nthmax;
								}
							}
							for(int k = 0; k < Nx; k += def_nthmax){
								int Nk = min(def_nthmax, Nx - k);
								InterpolateX2_kernel <<< (Nk + 511) / 512, 512 >>> (K_d + iP * Nx, Kc_d, Nx, param.Nxb, param.useIndividualX, binBoundaries_d, k);
								InterpolateX1_kernel <<< (Nk + 511) / 512, 512 >>> (K_d + iP * Nx, K1_d, Nx, param.Nxb, param.useIndividualX, binBoundaries_d, k);
							}
						} // end param.RLOW == 1
						//search second order regimes of the Voigt profile
						const int nl2 = 512;
						for(int il = 0; il < NL; il += nl2){ //loop over lines
							int ii11 = 0;
							int ii00 = Nx;
							if(param.useIndividualX == 0){
								for(int iil = 0; iil < nl2; ++iil){
									if(il + iil < NL){

										double aD2 = 1.0 / (ialphaDP[il + iil] * ialphaDP[il + iil]);
										double aL2 = vyP[il + iil] * vyP[il + iil] * aD2;
										double Dnu2 = 1.0e6 * aD2 - aL2;
										double Dnu = 0.0;
										if(Dnu2 > 0.0){
											Dnu = sqrt(Dnu2);
											int Inu = (int)((nuP[il + iil] - param.numin) / param.dnu);
											int ii0 = Inu - (int)(Dnu / param.dnu) - 1;
											int ii1 = Inu + (int)(Dnu / param.dnu) + 2;
	//printf("%d %d %d\n", il + iil, ii0, ii1);

											ii11 = max(ii11, ii1);
											ii00 = min(ii00, ii0);
										}
									}
								}
							}
							else{
								double nu00 = param.numax;
								double nu11 = param.numin;
								for(int iil = 0; iil < nl2; ++iil){
									if(il + iil < NL){
										double aD2 = 1.0 / (ialphaDP[il + iil] * ialphaDP[il + iil]);
										double aL2 = vyP[il + iil] * vyP[il + iil] * aD2;
										double Dnu2 = 1.0e6 * aD2 - aL2;
										double Dnu = 0.0;
										if(Dnu2 > 0.0){
											Dnu = sqrt(Dnu2);
											double nu0 = nuP[il + iil] - Dnu;
											double nu1 = nuP[il + iil] + Dnu;

											nu11 = fmax(nu11, nu1);
											nu00 = fmin(nu00, nu0);
										}
									}
								}
								for(int bin = 0; bin < param.nbins; ++bin){
									if(binBoundaries_h[bin + 1] > nu11){
										double dnu = (binBoundaries_h[bin + 1] - binBoundaries_h[bin]) / ((double)(param.Nxb));
										int bstart = bin * param.Nxb;
										ii11 = (nu11 - binBoundaries_h[bin]) / dnu + bstart + 2;
										break;
									}
								}
								for(int bin = 0; bin < param.nbins; ++bin){
									if(binBoundaries_h[bin + 1] > nu00){
										double dnu = (binBoundaries_h[bin + 1] - binBoundaries_h[bin]) / ((double)(param.Nxb));
										int bstart = bin * param.Nxb;
										ii00 = (nu00 - binBoundaries_h[bin]) / dnu + bstart - 1;
										break;
									}
								}
							}

							ii11 = min(Nx, ii11);
							ii00 = max(0, ii00);
							int nt = ii11 - ii00;
							int nstart = ii00;
							int nll = min(nl2, NL - il);	
if(il % 10000 == 0) printf("B %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
							for(int k = 0; k < nt; k += def_nthmax){
								int Nk = min(def_nthmax, nt - k);
								if(Nk > 0 && nll > 0){
									Line2f_kernel < nl2, 1 > <<< (max(Nk, nll) + nl2 - 1) / nl2, nl2 >>> (L.S1f_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, K_d + iP * Nx, il, nstart, Nk, nll, param.useIndividualX, param.Nxb, binBoundaries_d, 0.0f, 0.0f, 0.0f, param.profile);
								}
								nstart += def_nthmax;
							}
						}

						//search higher order regimes of the Voigt profile
						const int nl3 = 512;
						float a = (float)(M_PI * sqrt(-1.0 / log(def_TOLF * 0.5)));
						float b = (float)(1.0 / sqrt(M_PI));
						float c = (float)(2.0 * a / M_PI);

						for(int il = 0; il < NL; il += nl3){ //loop over lines
							int ii11 = 0;
							int ii00 = Nx;
							if(param.useIndividualX == 0){
								for(int iil = 0; iil < nl3; ++iil){
									if(il + iil < NL){

										double aD2 = 1.0 / (ialphaDP[il + iil] * ialphaDP[il + iil]);
										double aL2 = vyP[il + iil] * vyP[il + iil] * aD2;
										double Dnu2 = 1.0e2 * aD2 - aL2;
										double Dnu = 0.0;
										if(Dnu2 > 0.0){
											Dnu = sqrt(Dnu2);
											int Inu = (int)((nuP[il + iil] - param.numin) / param.dnu);
											int ii0 = Inu - (int)(Dnu / param.dnu) - 1;
											int ii1 = Inu + (int)(Dnu / param.dnu) + 2;
	//printf("%d %d %d\n", il + iil, ii0, ii1);
											ii11 = max(ii11, ii1);
											ii00 = min(ii00, ii0);
										}
									}
								}
							}
							else{
								double nu00 = param.numax;
								double nu11 = param.numin;
								for(int iil = 0; iil < nl3; ++iil){
									if(il + iil < NL){
										double aD2 = 1.0 / (ialphaDP[il + iil] * ialphaDP[il + iil]);
										double aL2 = vyP[il + iil] * vyP[il + iil] * aD2;
										double Dnu2 = 1.0e2 * aD2 - aL2;
										double Dnu = 0.0;
										if(Dnu2 > 0.0){
											Dnu = sqrt(Dnu2);
											double nu0 = nuP[il + iil] - Dnu;
											double nu1 = nuP[il + iil] + Dnu;

											nu11 = fmax(nu11, nu1);
											nu00 = fmin(nu00, nu0);
										}
									}
								}
								for(int bin = 0; bin < param.nbins; ++bin){
									if(binBoundaries_h[bin + 1] > nu11){
										double dnu = (binBoundaries_h[bin + 1] - binBoundaries_h[bin]) / ((double)(param.Nxb));
										int bstart = bin * param.Nxb;
										ii11 = (nu11 - binBoundaries_h[bin]) / dnu + bstart + 2;
										break;
									}
								}
								for(int bin = 0; bin < param.nbins; ++bin){
									if(binBoundaries_h[bin + 1] > nu00){
										double dnu = (binBoundaries_h[bin + 1] - binBoundaries_h[bin]) / ((double)(param.Nxb));
										int bstart = bin * param.Nxb;
										ii00 = (nu00 - binBoundaries_h[bin]) / dnu + bstart - 1;
										break;
									}
								}
							}
							ii11 = min(Nx, ii11);
							ii00 = max(0, ii00);
							int nt = ii11 - ii00;
							int nstart = ii00;
							int nll = min(nl3, NL - il);	
if(il % 10000 == 0) printf("C %d %d %d %d %d\n",il, ii00, ii11, nll, nt);
							for(int k = 0; k < nt; k += def_nthmax){
								int Nk = min(def_nthmax, nt - k);
								if(Nk > 0 && nll > 0){
									Line2f_kernel < nl3, 2 > <<< (max(Nk, nll) + nl3 - 1) / nl3, nl3 >>> (L.Sf_d, L.vyf_d, L.va_d, L.vb_d, L.vcut2_d, K_d + iP * Nx, il, nstart, Nk, nll, param.useIndividualX, param.Nxb, binBoundaries_d, a, b, c, param.profile);
								}
								nstart += def_nthmax;
							}
						}

					} //end Nx > def_NXLOW
					//*************************************
					cudaEventRecord(stop);
					if(iP < param.nP - 1){
						//synchronize here only if no more data has to be read from the disk.
						//otherwise read data before synchronization
						cudaEventSynchronize(stop);
						cudaEventElapsedTime(&milliseconds, start, stop);

						time[2] += milliseconds * 0.001;
						if(iP == param.nP - 1){
							printf("Time for K(x):         %g seconds\n", time[2]);
						}
			
						cudaDeviceSynchronize();
						error = cudaGetLastError();
						if(error != 0){
							printf("Kb error = %d = %s\n",error, cudaGetErrorString(error));
							return 0;
						}
						gettimeofday(&tt1, NULL);
					}

				} // End of pressure loop

			} // End of maxLines loop
		} // End of linefile loop

		cudaEventSynchronize(stop);
		if(m.NL[fi] > 0){
			cudaEventElapsedTime(&milliseconds, start, stop);

			time[2] += milliseconds * 0.001;
		}
		printf("Time for K(x):         %g seconds\n", time[2]);

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != 0){
			printf("Kc error = %d = %s\n",error, cudaGetErrorString(error));
			return 0;
		}
			InfoFile = fopen(InfoFilename, "a");
			fprintf(InfoFile,"File %d of %d\n", fi, fi);
			fprintf(InfoFile,"Number of lines: %lld\n", m.NL[fi]);
			fprintf(InfoFile,"Time for input:        %g seconds\n", time[0]);
			fprintf(InfoFile,"Time for Lines:        %g seconds\n", time[1]);
			fprintf(InfoFile,"Time for K(x):         %g seconds\n", time[2]);
			fclose(InfoFile);
		gettimeofday(&tt1, NULL);

		free(nuP);
		free(ialphaDP);
		free(vyP);
	}
	cudaFree(Limits_d);
	cudaFree(MaxLimits_d);
	free(binBoundaries_h);
	cudaFree(binIndex_d);
	cudaFree(binBoundaries_d);	
	cudaFree(K1_d);
	cudaFree(Kc_d);

	//***************************
	//Write the full line profile
	//****************************
	if(param.doStoreFullK == 1){
		FILE *OutFile;
		char OutFilename[160];
		sprintf(OutFilename, "Out_%s.dat", param.name);
			
		OutFile = fopen(OutFilename, filemode);

		for(int iP = 0; iP < param.nP; ++iP){
			cudaMemcpy(K_h, K_d + iP * Nx, Nx * sizeof(double), cudaMemcpyDeviceToHost);
			for(int j = 0; j < Nx; ++j){

				if(param.nP == 1){
					fprintf(OutFile, "%.20g %.20g\n", x_h[j], K_h[j] * unitScale);
				}
				else{
					fprintf(OutFile, "%.20g %.20g %.20g %.20g\n", x_h[j], K_h[j] * unitScale, param.T, P_h[iP]);
				}
			}
			fprintf(OutFile, "\n\n");
		}
		fclose(OutFile);
	}
	if(param.doStoreFullK == -1){
		FILE *OutFile;
		char OutFilename[160];
		sprintf(OutFilename, "%sOut_%s.dat", param.pathK, param.name);
			
		OutFile = fopen(OutFilename, "r");
		if(OutFile == NULL){
			printf("Error: Input file not found %s\n", OutFilename);
			return 0;
		}

		for(int iP = 0; iP < param.nP; ++iP){
			for(int j = 0; j < Nx; ++j){

				if(param.nP == 1){
					double k;
					fscanf(OutFile, "%lf %lf\n", &x_h[j], &k);
					K_h[j] = k / unitScale;
				}
				else{
					double k, t, p;
					fscanf(OutFile, "%lf %lf %lf %lf\n", &x_h[j], &k, &t, &p);
					K_h[j] = k / unitScale;
				}
			}
			cudaMemcpy(K_d + iP * Nx, K_h, Nx * sizeof(double), cudaMemcpyHostToDevice);
			fscanf(OutFile, "\n\n");
		}
		fclose(OutFile);
	}
	if(param.doStoreFullK == 2){
		//write a binary file in single precision
		FILE *OutFile;
		char OutFilename[160];
		sprintf(OutFilename, "Out_%s.bin", param.name);
			
		if(param.replaceFiles == 0){
			OutFile = fopen(OutFilename, "ab");
		}
		else{
			OutFile = fopen(OutFilename, "wb");
		}

		for(int iP = 0; iP < param.nP; ++iP){
			cudaMemcpy(K_h, K_d + iP * Nx, Nx * sizeof(double), cudaMemcpyDeviceToHost);
			for(int j = 0; j < Nx; ++j){
				float Kf = (float)(K_h[j]);
				fwrite(&Kf, sizeof(float), 1, OutFile);
			}
		}
		fclose(OutFile);
	}
	if(param.doStoreFullK == -2){
		//read a binary file
		FILE *OutFile;
		char OutFilename[160];
		sprintf(OutFilename, "%sOut_%s.bin", param.pathK, param.name);
			
		OutFile = fopen(OutFilename, "rb");
		if(OutFile == NULL){
			printf("Error: Input file not found %s\n", OutFilename);
			return 0;
		}

		for(int iP = 0; iP < param.nP; ++iP){
			for(int j = 0; j < Nx; ++j){
				float Kf;
				fread(&Kf, sizeof(float), 1, OutFile);
				K_h[j] = (double)(Kf);
			}
			cudaMemcpy(K_d + iP * Nx, K_h, Nx * sizeof(double), cudaMemcpyHostToDevice);
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
	if(param.doMean > 0){
		
		double *Pmn_d;
		double *Rmn_d;

		cudaMalloc((void **) &Pmn_d, Nx * sizeof(double));
		cudaMalloc((void **) &Rmn_d, Nx * sizeof(double));
	
		double *means_h, *means_d;	
		means_h = (double*)malloc(4 * sizeof(double));
		cudaMalloc((void **) &means_d, 4 * sizeof(double));

		FILE *Out4File;
		char Out4Filename[160];

		sprintf(Out4Filename, "Out_%s_mean.dat", param.name);
		Out4File = fopen(Out4Filename, filemode);
	
		for(int iP = 0; iP < param.nP; ++iP){

			Mean_kernel <<< (Nx + 511) / 512, 512 >>> (x_d, Pmn_d, Rmn_d, param.T, Nx);
/*
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
			IntegrateMean_kernel <512> <<< 4, 512 >>> (Pmn_d, Rmn_d, x_d, K_d + iP * Nx, means_d, Nx, param.useIndividualX);
			double sigma = 2.0 * def_kB * def_kB * def_kB * def_kB / ( def_h * def_h * def_h * def_c * def_c * 15.0) * M_PI * M_PI * M_PI * M_PI * M_PI;
			double integral1 = sigma * param.T * param.T * param.T * param.T / M_PI;
			double integral2 = M_PI / (4.0 * sigma * param.T * param.T * param.T);
		
			cudaMemcpy(means_h, means_d, 4 * sizeof(double), cudaMemcpyDeviceToHost);


			if(param.nP == 1){
				fprintf(Out4File, "%.20g\n", means_h[0] / means_h[2]);
				fprintf(Out4File, "%.20g\n", means_h[3] / means_h[1]);
				fprintf(Out4File, "%.20g\n", means_h[2]);
				fprintf(Out4File, "%.20g\n", integral1);
				fprintf(Out4File, "%.20g\n", means_h[3]);
				fprintf(Out4File, "%.20g\n", 1.0 / integral2);
			}
			else{
				fprintf(Out4File, "%.20g %.20g %.20g\n", means_h[0] / means_h[2], param.T, P_h[iP]);
				fprintf(Out4File, "%.20g %.20g %.20g\n", means_h[3] / means_h[1], param.T, P_h[iP]);
				fprintf(Out4File, "%.20g %.20g %.20g\n", means_h[2], param.T, P_h[iP]);
				fprintf(Out4File, "%.20g %.20g %.20g\n", integral1, param.T, P_h[iP]);
				fprintf(Out4File, "%.20g %.20g %.20g\n", means_h[3], param.T, P_h[iP]);
				fprintf(Out4File, "%.20g %.20g %.20g\n", 1.0 / integral2, param.T, P_h[iP]);

			}
			//fprintf(Out4File, "\n\n");
		}
		
		fclose(Out4File);

		free(means_h);
		cudaFree(means_d);
		cudaFree(Pmn_d);
		cudaFree(Rmn_d);
	}
	cudaFree(x_d);
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
	//***************************************
	thrust::device_ptr<double> K_dt = thrust::device_pointer_cast(K_d);
	thrust::device_ptr<int> binKey_dt = thrust::device_pointer_cast(binKey_d);
	for(int iP = 0; iP < param.nP; ++iP){
		thrust::sort_by_key(K_dt + iP * Nx, K_dt + Nx + iP * Nx, binKey_dt);
		thrust::stable_sort_by_key(binKey_dt, binKey_dt + Nx, K_dt + iP * Nx);
	}
	cudaFree(binKey_d);
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

	//*********************************
	//Prepare Resampling and do QR factorization, the same for all bins
	// this doesn't work with individual bins
	//*********************************

//size_t free_byte;
//size_t total_byte;
//cudaMemGetInfo( &free_byte, &total_byte );
//printf("***MEMRORY %g %g %g\n", (double)(free_byte), (double)(total_byte), (double)(total_byte) - (double)(free_byte));
	int *Nxmin_h, *Nxmin_d;		
	Nxmin_h = (int*)malloc(param.nbins * sizeof(int));
	cudaMalloc((void **) &Nxmin_d, param.nbins * sizeof(int));
	for(int i = 0; i < param.nbins; ++i){
		Nxmin_h[i] = 0;
	}
	cudaMemset(Nxmin_d, 0, param.nbins * sizeof(int));
	if(param.doResampling > 0){

		double *K2_h, *K2_d;
		K2_h = (double*)malloc(Nx * sizeof(double));
		cudaMalloc((void **) &K2_d, Nx * sizeof(double));
//cudaMemGetInfo( &free_byte, &total_byte );
//printf("***MEMRORY %g %g %g\n", (double)(free_byte), (double)(total_byte), (double)(total_byte) - (double)(free_byte));

		double *V_d;			//Vandermonde like matrix for least sqaures
		double *C_d, *D_d;

		cudaMalloc((void **) &V_d, param.nC * param.Nxb * sizeof(double));
		cudaMalloc((void **) &C_d, param.nC * sizeof(double));
		cudaMalloc((void **) &D_d, param.nC * sizeof(double));

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != 0){
			printf("Resampling Allocation error = %d = %s\n",error, cudaGetErrorString(error));
			return 0;
		}

		Vandermonde_kernel <<< (param.Nxb + 511) / 512, 512 >>> (V_d, (double)(param.Nxb), param.nC);
		QR_kernel <512> <<< 1, 512 >>> (V_d, C_d, D_d, param.Nxb, param.nC);

		FILE *Out3File;
		char Out3Filename[160];
		if(param.doResampling == 1){
			sprintf(Out3Filename, "Out_%s_cbin.dat", param.name);
			Out3File = fopen(Out3Filename, filemode);
		}
		if(param.doResampling == 2){
			if(param.replaceFiles == 1){
				for(int i = 0; i < param.nbins; ++i){
					sprintf(Out3Filename, "Out_%s_cbin%.4d.dat", param.name, i);
					Out3File = fopen(Out3Filename, "w");
					fclose(Out3File);	
				}
			}
			sprintf(Out3Filename, "Out_%s_cbin%.4d.dat", param.name, 0);
			Out3File = fopen(Out3Filename, "a");
		}

		for(int iP = 0; iP < param.nP; ++iP){	
			if(param.doResampling == 2 && iP > 0){
				fclose(Out3File);
				sprintf(Out3Filename, "Out_%s_cbin%.4d.dat", param.name, 0);
				Out3File = fopen(Out3Filename, "a");
			}
			cudaMemset(K2_d, 0, Nx * sizeof(double));
			cudaMemset(Nxmin_d, 0, param.nbins * sizeof(int));

			findCut_kernel <<< (Nx + 511) / 512, 512 >>> (K_d + iP * Nx, Nx, param.Nxb, param.kmin, Nxmin_d, param.nbins);
			rescale_kernel < 512 > <<< param.nbins, 512 >>> (Nxmin_d, K_d + iP * Nx, K2_d, param.Nxb, param.kmin, 1);
/*
cudaMemcpy(K2_h, K2_d, Nx * sizeof(double), cudaMemcpyDeviceToHost);
cudaMemcpy(K_h, K_d + iP * Nx, Nx * sizeof(double), cudaMemcpyDeviceToHost);
cudaDeviceSynchronize();
//printf only cut and empty bins
for(int i = 0; i < param.nbins; ++i){
	int il = i * param.Nxb;
	if(K_h[il] == param.kmin){
		for(int j = 0; j < param.Nxb; ++j){
//			printf("%g %.20g\n", j / (double)(param.Nxb), K2_h[j + il]);
		}
//		printf("\n\n");
	}
}
//print all bins
for(int i = 0; i < Nx; ++i){
	printf("%d %.20g %.20g\n", i, K_h[i], K2_h[i]);
}
*/
			copyK2_kernel< 512 > <<< param.nbins, 512 >>> (Nxmin_d, K_d + iP * Nx, K2_d, param.Nxb);
			cudaMemcpy(Nxmin_h, Nxmin_d, param.nbins * sizeof(int), cudaMemcpyDeviceToHost);
	
			lnK_kernel <<< (Nx + 511) / 512, 512 >>> (K_d + iP * Nx, Nx);
			leastSquare_kernel <512> <<< param.nbins, 512 >>> (V_d, C_d, D_d, K_d + iP * Nx, param.Nxb, param.nC);

			for(int i = 0; i < param.nbins; ++i){
				int il = i * param.Nxb;
				cudaMemcpy(K_h + il, K_d + il + iP * Nx, param.nC * sizeof(double), cudaMemcpyDeviceToHost);
		
				fprintf(Out3File, "%.20g %.20g ", param.kmin, fmin(Nxmin_h[i] / ((double)(param.Nxb - 1)), 1.0));
				for(int ic = 0; ic < param.nC; ++ic){
					if(Nxmin_h[i] != param.Nxb) fprintf(Out3File, "%.20g ", K_h[il + ic]);
					else fprintf(Out3File, "0.0 ");
				}
				if(param.nP > 1){
					fprintf(Out3File, "%.20g %.20g ", param.T, P_h[iP]);
				}
				if(param.doResampling == 1){
					fprintf(Out3File, "\n\n");
				}
				if(param.doResampling == 2 && i < param.nbins - 1){
					fprintf(Out3File, "\n");
					fclose(Out3File);
					sprintf(Out3Filename, "Out_%s_cbin%.4d.dat", param.name, i + 1);
					Out3File = fopen(Out3Filename, "a");
				}
			}
			//fprintf(Out3File, "\n\n");
			if(param.doTransmission > 0 || param.doStoreK > 0){
				expfx_kernel <<< param.nbins, 512 >>> (K_d + iP * Nx, param.nC, param.Nxb);
				rescale_kernel < 512 > <<< param.nbins, 512 >>> (Nxmin_d, K_d + iP * Nx, K2_d, param.Nxb, param.kmin, -1);
				copyK2_kernel< 512 > <<< param.nbins, 512 >>> (Nxmin_d, K_d + iP * Nx, K2_d, param.Nxb);
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
	if(param.doStoreK > 0){
		FILE *Out2File;
		char Out2Filename[160];
		if(param.doStoreK == 1){
			sprintf(Out2Filename, "Out_%s_bin.dat", param.name);
			Out2File = fopen(Out2Filename, filemode);
		}
		if(param.doStoreK == 2){
			if(param.replaceFiles == 1){
				for(int i = 0; i < param.nbins; ++i){
					sprintf(Out2Filename, "Out_%s_bin%.4d.dat", param.name, i);
					Out2File = fopen(Out2Filename, "w");
					fclose(Out2File);	
				}
			}
			sprintf(Out2Filename, "Out_%s_bin%.4d.dat", param.name, 0);
			Out2File = fopen(Out2Filename, "a");
		}
		
		for(int iP = 0; iP < param.nP; ++iP){
			if(param.doStoreK == 2 && iP > 0){
				fclose(Out2File);
				sprintf(Out2Filename, "Out_%s_bin%.4d.dat", param.name, 0);
				Out2File = fopen(Out2Filename, "a");
			}
			cudaMemcpy(K_h, K_d + iP * Nx, Nx * sizeof(double), cudaMemcpyDeviceToHost);
			if(param.useIndividualBins == 0){
				for(int i = 0; i < param.nbins; ++i){
					int Nxb = param.Nxb;
					int il = i * Nxb;
					int iedge = 0; //index of edge
					int nedge = 0; //number of points per edge intervall
					double sedge = 0.0; //sum of points in edge intervall
					for(int j = 0; j < Nxb; ++j){
						double y = j / ((double)(Nxb - 1));
						double y1 = (j + 1) / ((double)(Nxb - 1));
						if(param.useOutputEdges == 0){
							if(param.nP == 1){
								fprintf(Out2File, "%g %.20g\n", y, K_h[j + il] * unitScale);
							}
							else{
								fprintf(Out2File, "%g %.20g %g %g %d\n", y, K_h[j + il] * unitScale, param.T, P_h[iP], j);
							}
						}
						else{
							double edge = outputEdges_h[iedge];
							++nedge;
							sedge += K_h[j + il] * unitScale;
							if(y <= edge && edge <= y1 && iedge < param.nedges){
								if(param.nP == 1){
									if(iedge > 0) fprintf(Out2File, "%g %.20g\n", 0.5 * (edge + outputEdges_h[iedge - 1]), sedge / ((double)(nedge)));
								}
								else{
									if(iedge > 0) fprintf(Out2File, "%g %.20g %g %g %d\n", 0.5 * (edge + outputEdges_h[iedge - 1]), sedge / ((double)(nedge)), param.T, P_h[iP], iedge - 1);

								}
								++iedge;
								nedge = 0;
								sedge = 0.0;
							}
						}
					}
					if(param.doStoreK == 1){
						fprintf(Out2File,"\n\n");
					}
					if(param.doStoreK == 2 && i < param.nbins - 1){
						fclose(Out2File);
						sprintf(Out2Filename, "Out_%s_bin%.4d.dat", param.name, i + 1);
						Out2File = fopen(Out2Filename, "a");
					}
				}
			}
			else{
				int ib = 0;
				int j = 0;
				int iedge = 0; //inde of edge
				int nedge = 0; //number of points per edge intervall
				double sedge = 0.0; //sum of points in edge intervall
				for(int i = 0; i < Nx; ++i){
					int il = binIndex_h[ib];
					int ir = binIndex_h[ib + 1];
					int Nxb = ir - il;

					double y = j / ((double)(Nxb - 1));
					double y1 = (j + 1) / ((double)(Nxb - 1));

					if(param.useOutputEdges == 0){
						if(param.nP == 1){
							fprintf(Out2File, "%g %.20g\n", y, K_h[i] * unitScale);
						}
						else{
							fprintf(Out2File, "%g %.20g %.20g %.20g %d\n", y, K_h[i] * unitScale, param.T, P_h[iP], j);
						}
					}
					else{
						double edge = outputEdges_h[iedge];
						++nedge;
						sedge += K_h[i] * unitScale;
						if(y <= edge && edge <= y1 && iedge < param.nedges){
							if(param.nP == 1){
								if(iedge > 0) fprintf(Out2File, "%g %.20g\n", 0.5 * (edge + outputEdges_h[iedge - 1]), sedge / ((double)(nedge)));
							}
							else{
								if(iedge > 0) fprintf(Out2File, "%g %.20g %.20g %.20g %d\n", 0.5 * (edge + outputEdges_h[iedge - 1]), sedge / ((double)(nedge)), param.T, P_h[iP], iedge - 1);
							}
							++iedge;
							nedge = 0;
							sedge = 0.0;
						}
					}
					++j;

					if(i >= ir - 1){
//printf("%d %d %d %d\n", ib, il, ir, Nxb);
						++ib;
						j = 0;
						if(param.doStoreK == 1){
							fprintf(Out2File,"\n\n");
						}
						if(param.doStoreK == 2 && ib < param.nbins){
							fclose(Out2File);
							sprintf(Out2Filename, "Out_%s_bin%.4d.dat", param.name, ib);
							Out2File = fopen(Out2Filename, "a");
						}
						iedge = 0;
					}
					if(ib >= param.nbins){
						break;
					}
				}
			}
		}//end of P loop
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
	if(param.doTransmission > 0 ){

		double *Tr_h, *Tr_d;
		Tr_h = (double*)malloc(param.nbins * param.nTr * sizeof(double));
		cudaMalloc((void **) &Tr_d, param.nbins * param.nTr * sizeof(double));

		FILE *Out3File;
		char Out3Filename[160];

		if(param.doTransmission == 1){
			sprintf(Out3Filename, "Out_%s_tr.dat", param.name);
			Out3File = fopen(Out3Filename, filemode);
		}
		if(param.doTransmission == 2){
			if(param.replaceFiles == 1){
				for(int i = 0; i < param.nbins; ++i){
					sprintf(Out3Filename, "Out_%s_tr%.4d.dat", param.name, i);
					Out3File = fopen(Out3Filename, "w");
					fclose(Out3File);	
				}
			}
			sprintf(Out3Filename, "Out_%s_tr%.4d.dat", param.name, 0);
			Out3File = fopen(Out3Filename, "a");
		}

		for(int iP = 0; iP < param.nP; ++iP){
			if(param.doTransmission == 2 && iP > 0){
				fclose(Out3File);
				sprintf(Out3Filename, "Out_%s_tr%.4d.dat", param.name, 0);
				Out3File = fopen(Out3Filename, "a");
			}
			Integrate_kernel < 512 > <<< param.nbins, 512 >>> (K_d + iP * Nx, Tr_d, param.Nxb, param.nTr, param.dTr, Nxmin_d, param.kmin);
			cudaMemcpy(Tr_h, Tr_d, param.nbins * param.nTr * sizeof(double), cudaMemcpyDeviceToHost);
			for(int i = 0; i < param.nbins; ++i){
				for(int j = 0; j < param.nTr; ++j){
					double m = exp((j - param.nTr/2) * param.dTr);
					if(param.nP == 1){
						fprintf(Out3File, "%.20g %.20g\n", m, Tr_h[i * param.nTr + j]);
					}
					else{
						fprintf(Out3File, "%.20g %.20g %.20g %.20g %d\n", m, Tr_h[i * param.nTr + j], param.T, P_h[iP], j);
					}
				}
				if(param.doTransmission == 1){
					fprintf(Out3File, "\n\n");
				}
				if(param.doTransmission == 2 && i < param.nbins - 1){
					fclose(Out3File);
					sprintf(Out3Filename, "Out_%s_tr%.4d.dat", param.name, i + 1);
					Out3File = fopen(Out3Filename, "a");
				}
			}
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


	if(param.useHITEMP < 2) free_Line(L);
	else free2_Line(L);
	free(MaxLimits_h);
	free(K_h);
	free(x_h);
	free(Nxmin_h);
	free(outputEdges_h);
	free(binIndex_h);

	cudaFree(K_d);
	cudaFree(Nxmin_d);

	error = cudaGetLastError();
	printf("Final error = %d = %s\n",error, cudaGetErrorString(error));

	return 0;
}
