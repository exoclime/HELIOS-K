#include "define.h" //must be on top for Windows compilation

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>


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

	sprintf(param.PFilename, "%s", "-");
	sprintf(param.SpeciesFilename, "%s", "-");
	sprintf(param.edges, "%s", "-");
	sprintf(param.bins, "%s", "-");
	sprintf(param.ciaSystem, "%s", "-");

	param.dev = 0;
	param.useIndividualBins = 0;
	param.useOutputEdges = 0;
	param.nedges = 0;
	param.nP = 1;
	param.usePFile = 0;
	param.useIndividualX = 0;
	param.useCia = 0;
	param.path[0] = 0;
	param.pathK[0] = 0;
	param.nSpecies = 1;
	param.useSpeciesFile = 0;
	param.useSubLorentzian = 0;

	param.T = 0.0;
	param.P = 0.0;
	param.mParamFilename[0] = 0;
	param.dataBase = 0;
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
	param.gammaF = def_gammaF;
	param.doMean = 0;
	param.units = 0;	
	param.replaceFiles = 0;
	param.profile = def_PROFILE;
	param.doTuning = def_doTuning;
	param.removePlinth = def_removePlinth;
	
	er = read_parameters(param, paramFilename, argc, argv);
	if(er == 0){
		return 0;
	}
	if(param.dev >= devCount || param.dev < 0){
		printf("Error: Device Number is not allowed\n");
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
	char InfoFilename[300];
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
				devProp.multiProcessorCount,  devProp.clockRate, devProp.memoryClockRate, (long unsigned int)(devProp.totalGlobalMem), (long unsigned int)(devProp.sharedMemPerBlock));

		}
	}
	if(param.Nxb != 0){
		param.useIndividualX = 1;
	}
	if(param.removePlinth == 1 && param.profile == 4){
		printf("Error, remove plinth is not supported for profile 4\n");
		return 0;
	}

	subLorentzianConstantCopy(param.useSubLorentzian);

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
	if(param.useSubLorentzian == 1){
		subLorentzianB(param.T);
		param.useIndividualX = 1;
		//this is needed because of the nu/nu0 factor
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
	//Allocate Species array 
	double *SpeciesA_h;	//abundance
	char **SpeciesN_h;
	SpeciesA_h = (double*)malloc(param.nSpecies * sizeof(double));
	SpeciesN_h = (char**)malloc(param.nSpecies * sizeof(char*));
	for(int i = 0; i < param.nSpecies; ++i){
		SpeciesN_h[i] = (char*)malloc(160 * sizeof(char));
	}
	if(param.useSpeciesFile == 1){
		er = readSpeciesFile(param, SpeciesN_h, SpeciesA_h);
		if(er == 0) return 0;
	}


	double time[9];
	double timeT[3];
	for(int i = 0; i < 9; ++i){
		time[i] = 0.0;
	}
	for(int i = 0; i < 3; ++i){
		timeT[i] = 0.0;
	}

	//Allocate Molecule properties
	for(int i = 0; i < 2; ++i){
		FILE *infofile;
		if(i == 0) infofile = InfoFile;
		if(i == 1) infofile = stdout;
		fprintf(infofile, "\nVersion: %g\n", VERSION);
		fprintf(infofile, "Using device %d\n\n", param.dev);
		fprintf(infofile, "Runtime Version %d\n", runtimeVersion);
		fprintf(infofile, "Driver Version %d\n", driverVersion);
		fprintf(infofile, "GIT Describe: %s\n", GIT_DESCRIBE);
		fprintf(infofile, "Build Date: %s\n", BUILD_DATE);
		fprintf(infofile, "Build Path: %s\n", BUILD_PATH);
		fprintf(infofile, "Build System: %s\n", BUILD_SYSTEM);
		fprintf(infofile, "Build Compute Capability: SM=%s\n", BUILD_SM);
		fprintf(infofile, "\n");


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
			fprintf(infofile, "Number of P values: %d\n", param.nP);
		}
		if(param.useSpeciesFile > 0){
			fprintf(infofile, "Species in file: %s\n", param.SpeciesFilename);
			fprintf(infofile, "Number of Species: %d\n", param.nSpecies);
		}
		if(param.useSubLorentzian > 0){
			fprintf(infofile, "sub-Lorentzian file: %s\n", param.subLorentzianFilename);
		}
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
		fprintf(infofile, "gammaF = %g\n", param.gammaF);
		fprintf(infofile, "doMean = %d\n", param.doMean);
		fprintf(infofile, "Units = %d\n", param.units);
		fprintf(infofile, "Replace files = %d\n", param.replaceFiles);
		fprintf(infofile, "profile = %d\n", param.profile);
		fprintf(infofile, "doTuning = %d\n", param.doTuning);
		fprintf(infofile, "def_TOL = %g\n", def_TOL);
		fprintf(infofile, "def_TOLf = %g\n", def_TOLF);
		fprintf(infofile, "def_nthmax = %d\n", def_nthmax);
		fprintf(infofile, "def_nlmax = %d\n", def_nlmax);
		fprintf(infofile, "def_maxlines = %lld\n", def_maxlines);
		fprintf(infofile, "def_maxfiles = %d\n", def_maxfiles);
		fprintf(infofile, "def_NmaxSample = %d\n", def_NmaxSample);
		if(param.useOutputEdges == 1){
			fprintf(infofile, "use output edges: %s\n", param.edges);
		}
		fprintf(infofile, "\n");

	}
	fclose(InfoFile);

	cudaEvent_t tt1;			//start time
	cudaEvent_t tt2;			//end time
	cudaEventCreate(&tt1);
	cudaEventCreate(&tt2);


	cudaEvent_t ReadStart, ReadStop;
	cudaEventCreate(&ReadStart);
	cudaEventCreate(&ReadStop);

	cudaEvent_t KStart, KStop;
	cudaEventCreate(&KStart);
	cudaEventCreate(&KStop);

	cudaEvent_t LineStart, LineStop;
	cudaEventCreate(&LineStart);
	cudaEventCreate(&LineStop);

	cudaEvent_t TuneStart, TuneStop;
	cudaEventCreate(&TuneStart);
	cudaEventCreate(&TuneStop);

	cudaEvent_t iiLimitsEvent;
	cudaEventCreate(&iiLimitsEvent);

	cudaEvent_t AEvent;
	cudaEventCreate(&AEvent);
	cudaEvent_t ALEvent;
	cudaEventCreate(&ALEvent);
	cudaEvent_t AREvent;
	cudaEventCreate(&AREvent);
	cudaEvent_t BEvent;
	cudaEventCreate(&BEvent);

	float milliseconds;


	cudaStream_t VStream[def_KSn];
	for(int i = 0; i < def_KSn; ++i){
		cudaStreamCreate(&VStream[i]);
	}
	cudaStream_t CStream[def_rBs];
	for(int i = 0; i < def_rBs; ++i){
		cudaStreamCreate(&CStream[i]);
	}
	cudaStream_t tuneStream[2];
	for(int i = 0; i < 2; ++i){
		cudaStreamCreate(&tuneStream[i]);
	}
	cudaStream_t nuLimitsStream[5];
	for(int i = 0; i < 5; ++i){
		cudaStreamCreate(&nuLimitsStream[i]);
	}

	// ************************************************************
	//calculate mean mass before starting the opacity calculation
	//needed to set kmin
	// ************************************************************
	double meanMass = 0.0;
	for(int iSpecies = 0; iSpecies < param.nSpecies; ++iSpecies){
		double Sscale = 1.0;
		if(param.nSpecies > 1){
			sprintf(param.mParamFilename, "%s", SpeciesN_h[iSpecies]);
			Sscale = SpeciesA_h[iSpecies];
		}

		Molecule m;

		if(param.useCia == 0){
			int er = Init(m, param, qFilename);
			if(er == 0) return 0;
		}
		
		//compute the mean mass
		for(int i = 0; i < m.nISO; ++i){
//incldue here mixture abundances
			meanMass += m.ISO[i].Ab * m.ISO[i].m * Sscale; //mean Molar Mass (g)
		}
	}	
	printf("mean mass %g\n", meanMass);

	//needed here already to get the cia.mass1. Initialize it again later in the main species loop
	ciaSystem cia;
	if(param.useCia == 1){
		Molecule m;
		er = InitCia(m, cia, param);
		if(er == 0) return 0;
	}

	
	double unitScale = 1.0;
	if(param.units == 1){
		unitScale = 1.0 / def_NA * meanMass;
		if(param.useCia == 1){
			unitScale = 1.0 / def_NA * cia.mass1;
		}
		param.kmin /= unitScale;
		printf("kmin %g\n", param.kmin);
	}	
	// ************************************************************


	// ****************************************************************************
	// Allocate and initialize K and x arrays
	// ****************************************************************************
	double *K_h, *K_d;
	double *KS_d;	//used in multiple y blocks
	double *x_h, *x_d;
	int *binKey_d;
	int *binIndex_h, *binIndex_d;
	K_h = (double*)malloc(Nx * sizeof(double));
	x_h = (double*)malloc(Nx * sizeof(double));
	binIndex_h = (int*)malloc((param.nbins + 2) * sizeof(int));
	cudaMalloc((void **) &K_d, param.nP * Nx * sizeof(double));
	cudaMalloc((void **) &KS_d, def_KSn * Nx * sizeof(double));

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
	for(int k = 0; k < def_KSn * Nx; k += def_nthmax){
		int Nk = min(def_nthmax, def_KSn * Nx - k);
		//kmin must be here always zero, because the different streams are added later to K_d
		InitialK_kernel <<< (Nk + 511) / 512, 512 >>> (KS_d, def_KSn * Nx, 0.0, k);
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
	// ****************************************************************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("K and x alloc error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}

	//start species loop here
	for(int iSpecies = 0; iSpecies < param.nSpecies; ++iSpecies){

		double Sscale = 1.0;	//Abundance scale for mixtures
		if(param.nSpecies > 1){
			sprintf(param.mParamFilename, "%s", SpeciesN_h[iSpecies]);
			Sscale = SpeciesA_h[iSpecies];
		}

		Molecule m;

		m.id=0;
		m.NL[0] = 0;
		m.nISO = 0;
		m.defaultL = 0.0;
		m.defaultn = 0.0;
		//Initialize the Isotopologue properties for ISO.h
		if(param.useCia == 0){
			int er = Init(m, param, qFilename);
			if(er == 0) return 0;
		}

		//print species dependent information
		InfoFile = fopen(InfoFilename, "a");
		for(int i = 0; i < 2; ++i){
			FILE *infofile;
			if(i == 0) infofile = InfoFile;
			if(i == 1) infofile = stdout;
			fprintf(infofile, "Species Name = %s\n", m.mName);
			fprintf(infofile, "dataBase = %d\n", param.dataBase);
			fprintf(infofile, "Molecule Number = %d\n", m.id);
			fprintf(infofile, "default L = %g\n", m.defaultL);
			fprintf(infofile, "default n = %g\n", m.defaultn);
			fprintf(infofile, "\n");
		}
		fclose(InfoFile);

		
		//Read partition function
		Partition part;
		er = readPartition(param, qFilename, part, param.T, m);
		if(er == 0){
			return 0;
		}

		printf("mean mass %g, Sscale %g\n", meanMass, Sscale);

		//Set cia System properties
		ciaSystem cia;
		if(param.useCia == 1){
			er = InitCia(m, cia, param);
			if(er == 0) return 0;
		}

		if(param.useCia == 1 && m.id != 0){
			printf("Error, not allowed to use a cia system with a molecule\n");
			return 0;
		}

		double *readBuffer_h, *readBuffer_d;
		int readBufferSize = 8192;
		int readBufferN;
		int readBufferCount = 0;
		int rbvs = 0;
  	
		if(param.dataBase == 2){
			//Exomol
			readBufferN = 4;
		}
		if(param.dataBase == 30){
			//Kurucz
			readBufferN = 5;
		}
		if(param.dataBase == 31){
			//NIST
			readBufferN = 5;
		}
		if(param.dataBase == 32){
			//VALD
			readBufferN = 5;
		}

		cudaHostAlloc((void **) &readBuffer_h, def_rBs * readBufferSize * readBufferN * sizeof(double), cudaHostAllocDefault);
		cudaMalloc((void **) &readBuffer_d, def_maxlines * readBufferN * sizeof(double));

//printf("Allocate read Buffer %d %d %d %lld | %d %lld\n", def_rBs, readBufferSize, readBufferN, m.NLmax,  def_rBs * readBufferSize * readBufferN, m.NLmax * readBufferN);

		Line L;


		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != 0){
			printf("Initial error = %d = %s\n",error, cudaGetErrorString(error));
			return 0;
		}

		//Allocate memory for Line properties
		if(param.dataBase < 2 || param.dataBase == 3){
			// 0 1 3
			Alloc_Line(L, m, param);
		}
		else{
			// 2 30 31 32
			Alloc2_Line(L, m, param);
		}

		if(param.useCia == 1){
			for(int iP = 0; iP < param.nP; ++iP){
				int er = readCiaFile(param, cia, x_h, K_h, Nx, param.T, P_h[iP]);
				cudaMemcpy(K_d + iP * Nx, K_h, Nx * sizeof(double), cudaMemcpyHostToDevice);
				if(er == 0){
					return 0;	
				}
			}
		}

		cudaDeviceSynchronize();
		error = cudaGetLastError();
		if(error != 0){
			printf("Line alloc error = %d = %s\n",error, cudaGetErrorString(error));
			return 0;
		}


		if(m.id > 0 && param.doStoreFullK >= 0){

			// **************************************
			// Starting the loop around the datafiles
			// **************************************
			int fi0 = m.nFiles;
			int fi1 = 0;

			if(param.cut == 0.0) param.cut = 1.0e30;

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

			printf("File range %d to %d\n", fi0, fi1 - 1);

		
			int fi;
			FILE *dataFile;
			char dataFilename[180];
			
			timeT[0] += time[0];
			time[0] = 0.0;

			//Tuning parameters for Line Kernels
			int ntAOld = 0;
			int ntA = 0;
			int ntALOld = 0;
			int ntAL = 0;
			int ntAROld = 0;
			int ntAR = 0;
			int ntBOld = 0;
			int ntB = 0;
			int ntCOld = 0;
			int ntC = 0;
			int nkA = 8;
			int nkAL = 8;
			int nkAR = 8;
			int nkB = 4;
			int nkC = 2;
			
			double c1 = def_h * def_c / (def_kB * param.T);
			double T1 = def_T0 / param.T;

			for(fi = fi0; fi < fi1; ++fi){
				timeT[1] += time[1];
				timeT[2] += time[2];

				time[1] = 0.0;
				time[2] = 0.0;
				int NL;
				int NL1;
				long long lPart;

	
				// read the first block of files outside the loop
				// the remaining reads are called at the end of the loop
				// to allow overlapping execution
				// **************************read0
				
				if(fi == fi0){
					sprintf(dataFilename, "%sbin", m.dataFilename[fi]);
					dataFile = fopen(dataFilename, "rb");

					if(dataFile == NULL){
						printf("Error: line list file not found: %s\n", dataFilename);
						return 0;
					}

					printf("Reading Line file %d of %d: %s\n", fi, fi1 - 1, dataFilename);
					printf("Number of lines: %lld\n", m.NL[fi]);


					NL = min(def_maxlines, m.NL[fi] - 0);

					lPart = (0 + def_maxlines - 1) / def_maxlines;
					cudaEventRecord(ReadStart);
					printf("Reading Line file %d of %d; part %lld of %lld with %d lines\n", fi, fi1 - 1, lPart, (m.NL[fi] + def_maxlines - 1) / def_maxlines - 1, NL);
					// **************************
					// Read the Line list	
					// **************************
					if(param.dataBase < 2 || param.dataBase == 3){
						//0 1 3
						er = readFile(param, m, part, L, param.qalphaL, NL, dataFile, Sscale, meanMass);
					}
					else {
						// 2 30 31 32
						int vs = 0;	
						for(int i = 0; i < NL; i += readBufferSize){
							er = readFileExomol(L, NL, dataFile, readBuffer_h, readBuffer_d, readBufferSize, readBufferN, readBufferCount, vs, CStream);
							readBufferCount += readBufferSize;
							++vs;
						}
					}
					if(er == 0){
						return 0;
					}
					cudaEventRecord(ReadStop);
					cudaEventSynchronize(ReadStop);
					cudaEventElapsedTime(&milliseconds, ReadStart, ReadStop);
					time[0] += milliseconds * 0.001;
				
					printf("Reading Line file %d, part %lld complete\n", fi, lPart);	
					printf("Time for input, %d %lld:        %g seconds\n", fi, lPart, time[0]);
				}
				// **************************read0
				
				for(long long int iL = 0LL; iL < m.NL[fi]; iL += def_maxlines){

					//start the loop around the Pressure values. only 1 iteration if no Pressure file is given
					for(int iP = 0; iP < param.nP; ++iP){
						
						//Copy Line data to the device
						cudaEventRecord(LineStart);

						if(param.dataBase < 2 || param.dataBase == 3){
							//0 1 3
							Copy_Line(L, m, NL);
						}
						else{
							//2 30 31 32
							double mass = m.ISO[0].m / def_NA;
							double Abundance = m.ISO[0].Ab;
							if(param.units == 0){
								Abundance *= m.ISO[0].m / meanMass;
								Sscale *= m.ISO[0].m / meanMass;
							}
							double Q = part.Q[0];

							int vs = 0;
							for(int k = 0; k < NL; k += def_nthmax / 4){
								int Nk = min(def_nthmax / 4, NL - k);
								if(Nk > 0){
									// ***************************
									// Compute Line properties 1
									// ***************************
									if(param.dataBase == 2){
										L_kernelExomol  <<< (Nk + 127) / 128, 128, 0, VStream[vs % def_KSn] >>> (readBuffer_d, L.nu_d, L.S_d, L.EL_d, L.ialphaD_d, L.A_d, L.vy_d, L.n_d, m.defaultL, m.defaultn, param.gammaF, mass, param.T, Q, Abundance, Sscale, NL, k);
									}
									if(param.dataBase == 30){
										L_kernelKurucz  <<< (Nk + 127) / 128, 128, 0, VStream[vs % def_KSn] >>> (readBuffer_d, L.nu_d, L.S_d, L.EL_d, L.ialphaD_d, L.A_d, L.vy_d, L.n_d, m.defaultL, m.defaultn, param.gammaF, mass, param.T, Q, Abundance, Sscale, NL, k);
									}
									if(param.dataBase == 31){
										L_kernelNIST  <<< (Nk + 127) / 128, 128, 0, VStream[vs % def_KSn] >>> (readBuffer_d, L.nu_d, L.S_d, L.EL_d, L.ialphaD_d, L.A_d, L.vy_d, L.n_d, m.defaultL, m.defaultn, param.gammaF, mass, param.T, Q, Abundance, Sscale, NL, k);
									}
									if(param.dataBase == 32){
										L_kernelVALD  <<< (Nk + 127) / 128, 128, 0, VStream[vs % def_KSn] >>> (readBuffer_d, L.nu_d, L.S_d, L.EL_d, L.ialphaD_d, L.A_d, L.vy_d, L.n_d, m.defaultL, m.defaultn, param.gammaF, mass, param.T, Q, Abundance, Sscale, NL, k);
									}
									// ***************************
									// Compute Line properties 2
									// ***************************
									Sf_kernel <<< (Nk + 127) / 128, 128, 0, VStream[vs % def_KSn] >>> (L.nu_d, L.S_d, L.A_d, L.vy_d, L.ialphaD_d, L.n_d, L.EL_d, L.ID_d, NL, c1, T1, P_h[iP], k);
								}
								++vs;
							}
						}
						cudaDeviceSynchronize();
						// ************************

						// ***************************
						// Compute Line properties
						// ***************************
						if(param.dataBase < 2 || param.dataBase == 3){
							// 0 1 3
							for(int k = 0; k < NL; k += def_nthmax){
								int Nk = min(def_nthmax, NL - k);
								if(Nk > 0) S2_kernel <<< (Nk + 127) / 128, 128 >>> (L.nu_d, L.S_d, L.A_d, L.vy_d, L.ialphaD_d, L.n_d, L.delta_d, L.EL_d, L.ID_d, NL, param.T, P_h[iP], k);
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
//print_kernel <<< 1, 1 >>> (L.nu_d, L.ialphaD_d, L.vy_d, L.ID_d, 500, 0);

						//Sort the data along nu
						thrust::device_ptr<double> nu_dt = thrust::device_pointer_cast(L.nu_d);
						thrust::device_ptr<int> ID_dt = thrust::device_pointer_cast(L.ID_d);

						thrust::sort_by_key(nu_dt, nu_dt + NL, ID_dt);

						//Use Sort_d and ID_d to sort S_d, vy_d and ialphaD_d
						int Nk = min(def_nthmax, NL);
						for(int k = 0; k < NL; k += def_nthmax){
							if(Nk > 0) Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.S_d, L.Sort_d, NL, k);
						}
						for(int k = 0; k < NL; k += def_nthmax){
							if(Nk > 0) Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Sort_d, L.S_d, L.ID_d, NL, k);
						}
						for(int k = 0; k < NL; k += def_nthmax){
							if(Nk > 0) Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.vy_d, L.Sort_d, NL, k);
						}
						for(int k = 0; k < NL; k += def_nthmax){
							if(Nk > 0) Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Sort_d, L.vy_d, L.ID_d, NL, k);
						}
						for(int k = 0; k < NL; k += def_nthmax){
							if(Nk > 0) Copy_kernel <<< (Nk + 127) / 128, 128 >>> (L.ialphaD_d, L.Sort_d, NL, k);
						}
						for(int k = 0; k < NL; k += def_nthmax){
							if(Nk > 0) Sort_kernel <<< (Nk + 127) / 128, 128 >>> (L.Sort_d, L.ialphaD_d, L.ID_d, NL, k);
						}
						// ********************************

						for(int k = 0; k < NL; k += def_nthmax){
							int Nk = min(def_nthmax, NL - k);
							if(Nk > 0){
								S3_kernel <<< (Nk + 127) / 128, 128 >>> (L.nu_d, L.S_d, L.S1_d, L.vy_d, L.ialphaD_d, L.Sf_d, L.S1f_d, L.vyf_d, L.vcut2_d, L.va_d, L.vb_d, param.cut, param.cutMode, param.profile, param.numin, param.dnu, param.useIndividualX, NL, k);
								if(param.removePlinth == 1 && param.cut != 0.0){
									float a = (float)(M_PI * sqrt(-1.0 / log(def_TOLF * 0.5)));
									float b = (float)(1.0 / sqrt(M_PI));
									float c = (float)(2.0 * a / M_PI);
									Plinth_kernel <<< (Nk + 127) / 128, 128 >>> (L.S1f_d, L.Sf_d, L.vyf_d, L.vcut2_d, L.plinth_d, NL, a, b, c, param.profile);
									//printPlinth_kernel <<< (Nk + 127) / 128, 128 >>> (L.plinth_d, L.nu_d, NL);
								}
							}
						}

						cudaEventRecord(LineStop);
						cudaEventSynchronize(LineStop);
						error = cudaGetLastError();
						if(error != 0){
							printf("Line error = %d = %s\n",error, cudaGetErrorString(error));
							return 0;
						}
						cudaEventElapsedTime(&milliseconds, LineStart, LineStop);

						time[1] += milliseconds * 0.001;

						printf("Time for Lines: %d %lld %d:        %g seconds\n", fi, lPart, iP, time[1]);

						cudaEventRecord(KStart);

						// ************************************
						// Compute the opacity function K(x)
						// ************************************

						int nlLimitsA = (NL + def_nlA - 1) / def_nlA;
						int nlLimitsB = (NL + def_nlB - 1) / def_nlB;
						int nlLimitsC = (NL + def_nlC - 1) / def_nlC;


						//A
						nuLimits_kernel<<< nlLimitsA, min(def_nlA, 1024), 0, nuLimitsStream[0] >>> (L.nu_d, L.ialphaD_d, L.vy_d, L.vcut2_d, L.nuLimitsA0_d, L.nuLimitsA1_d, param.numin, param.numax, def_nlA, NL, param.profile, 10);
						iiLimits_kernel <<< (nlLimitsA + 127) / 128, 128, 0, nuLimitsStream[0] >>> (L.nuLimitsA0_d, L.nuLimitsA1_d, L.iiLimitsA0_d, L.iiLimitsA1_d, binBoundaries_d, nlLimitsA, param.numin, param.dnu, Nx, param.useIndividualX, param.nbins, param.Nxb, 10);

						iiLimitsMax_kernel< 512 > <<< 1, 512 >>> (L.iiLimitsA0_d, L.iiLimitsA1_d, L.iiLimitsAT_d, Nx, nlLimitsA);

						if(param.profile == 1){	//only for voigt profiles
							//AL
							nuLimits_kernel<<< nlLimitsA, min(def_nlA, 1024), 0, nuLimitsStream[1] >>> (L.nu_d, L.ialphaD_d, L.vy_d, L.vcut2_d, L.nuLimitsAL0_d, L.nuLimitsAL1_d, param.numin, param.numax, def_nlA, NL, param.profile, 11);
							iiLimits_kernel <<< (nlLimitsA + 127) / 128, 128, 0, nuLimitsStream[1] >>> (L.nuLimitsAL0_d, L.nuLimitsAL1_d, L.iiLimitsAL0_d, L.iiLimitsAL1_d, binBoundaries_d, nlLimitsA, param.numin, param.dnu, Nx, param.useIndividualX, param.nbins, param.Nxb, 11);

							iiLimitsMax_kernel< 512 > <<< 1, 512 >>> (L.iiLimitsAL0_d, L.iiLimitsAL1_d, L.iiLimitsALT_d, Nx, nlLimitsA);

							//AR
							nuLimits_kernel<<< nlLimitsA, min(def_nlA, 1024), 0, nuLimitsStream[2] >>> (L.nu_d, L.ialphaD_d, L.vy_d, L.vcut2_d, L.nuLimitsAR0_d, L.nuLimitsAR1_d, param.numin, param.numax, def_nlA, NL, param.profile, 12);
							iiLimits_kernel <<< (nlLimitsA + 127) / 128, 128, 0, nuLimitsStream[2] >>> (L.nuLimitsAR0_d, L.nuLimitsAR1_d, L.iiLimitsAR0_d, L.iiLimitsAR1_d, binBoundaries_d, nlLimitsA, param.numin, param.dnu, Nx, param.useIndividualX, param.nbins, param.Nxb, 12);

							iiLimitsMax_kernel< 512 > <<< 1, 512 >>> (L.iiLimitsAR0_d, L.iiLimitsAR1_d, L.iiLimitsART_d, Nx, nlLimitsA);
							//B
							nuLimits_kernel<<< nlLimitsB, min(def_nlB, 1024), 0, nuLimitsStream[3] >>> (L.nu_d, L.ialphaD_d, L.vy_d, L.vcut2_d, L.nuLimitsB0_d, L.nuLimitsB1_d, param.numin, param.numax, def_nlB, NL, param.profile, 20);
							iiLimits_kernel <<< (nlLimitsB + 127) / 128, 128, 0, nuLimitsStream[3] >>> (L.nuLimitsB0_d, L.nuLimitsB1_d, L.iiLimitsB0_d, L.iiLimitsB1_d, binBoundaries_d, nlLimitsB, param.numin, param.dnu, Nx, param.useIndividualX, param.nbins, param.Nxb, 20);

							iiLimitsMax_kernel< 512 > <<< 1, 512 >>> (L.iiLimitsB0_d, L.iiLimitsB1_d, L.iiLimitsBT_d, Nx, nlLimitsB);
							//C
							nuLimits_kernel<<< nlLimitsC, min(def_nlC, 1024), 0, nuLimitsStream[4] >>> (L.nu_d, L.ialphaD_d, L.vy_d, L.vcut2_d, L.nuLimitsC0_d, L.nuLimitsC1_d, param.numin, param.numax, def_nlC, NL, param.profile, 30);
							iiLimits_kernel <<< (nlLimitsC + 127) / 128, 128, 0, nuLimitsStream[4] >>> (L.nuLimitsC0_d, L.nuLimitsC1_d, L.iiLimitsC0_d, L.iiLimitsC1_d, binBoundaries_d, nlLimitsC, param.numin, param.dnu, Nx, param.useIndividualX, param.nbins, param.Nxb, 30);

							iiLimitsMax_kernel< 512 > <<< 1, 512 >>> (L.iiLimitsC0_d, L.iiLimitsC1_d, L.iiLimitsCT_d, Nx, nlLimitsC);
						}	


						cudaEventRecord(iiLimitsEvent);
						cudaEventSynchronize(iiLimitsEvent);
						iiLimitsCheck <<< (nlLimitsA + 127) / 128, 128 >>> (L.iiLimitsA0_d,  L.iiLimitsA1_d, L.iiLimitsAL0_d,  L.iiLimitsAL1_d, L.iiLimitsAR0_d,  L.iiLimitsAR1_d, nlLimitsA);
						cudaEventRecord(iiLimitsEvent);
						cudaEventSynchronize(iiLimitsEvent);


						long long int nTA = L.iiLimitsAT_m[1] - L.iiLimitsAT_m[0];
						long long int nTAL = L.iiLimitsALT_m[1] - L.iiLimitsALT_m[0];
						long long int nTAR = L.iiLimitsART_m[1] - L.iiLimitsART_m[0];
						long long int nTB = L.iiLimitsBT_m[1] - L.iiLimitsBT_m[0];
						long long int nTC = L.iiLimitsCT_m[1] - L.iiLimitsCT_m[0];

						if(ntA < 0) ntA = 0ll;
						if(ntAL < 0) ntAL = 0ll;
						if(ntAR < 0) ntAR = 0ll;
						if(ntB < 0) ntB = 0ll;
						if(ntC < 0) ntC = 0ll;
						
						printf("A Limits %lld %lld | %lld\n", L.iiLimitsAT_m[0], L.iiLimitsAT_m[1], nTA);
						printf("AL Limits %lld %lld | %lld\n", L.iiLimitsALT_m[0], L.iiLimitsALT_m[1], nTAL);
						printf("AR Limits %lld %lld | %lld\n", L.iiLimitsART_m[0], L.iiLimitsART_m[1], nTAR);
						printf("B Limits %lld %lld | %lld\n", L.iiLimitsBT_m[0], L.iiLimitsBT_m[1], nTB);
						printf("C Limits %lld %lld | %lld\n", L.iiLimitsCT_m[0], L.iiLimitsCT_m[1], nTC);


						if(nTA > 0){
							cudaMemcpyAsync(L.iiLimitsA0_h, L.iiLimitsA0_d, nlLimitsA * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[0]);
							cudaMemcpyAsync(L.iiLimitsA1_h, L.iiLimitsA1_d, nlLimitsA * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[0]);
						}
						if(nTAL > 0){
							cudaMemcpyAsync(L.iiLimitsAL0_h, L.iiLimitsAL0_d, nlLimitsA * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[1]);
							cudaMemcpyAsync(L.iiLimitsAL1_h, L.iiLimitsAL1_d, nlLimitsA * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[1]);
						}
						if(nTAR > 0){
							cudaMemcpyAsync(L.iiLimitsAR0_h, L.iiLimitsAR0_d, nlLimitsA * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[2]);
							cudaMemcpyAsync(L.iiLimitsAR1_h, L.iiLimitsAR1_d, nlLimitsA * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[2]);
						}

						if(nTB > 0){
							cudaMemcpyAsync(L.iiLimitsB0_h, L.iiLimitsB0_d, nlLimitsB * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[3]);
							cudaMemcpyAsync(L.iiLimitsB1_h, L.iiLimitsB1_d, nlLimitsB * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[3]);
						}

						if(nTC > 0){
							cudaMemcpyAsync(L.iiLimitsC0_h, L.iiLimitsC0_d, nlLimitsC * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[4]);
							cudaMemcpyAsync(L.iiLimitsC1_h, L.iiLimitsC1_d, nlLimitsC * sizeof(long long int), cudaMemcpyDeviceToHost, nuLimitsStream[4]);
						}



						double timeOld = time[0];
						long long lPartOld = lPart;
						int fii = fi;
						if((iL < m.NL[fi] - def_maxlines || fi < fi1 - 1) && iP == param.nP - 1){
							//read the next line file while calculating the K kernels of the current file
							// **************************read iL + 1
							int iLL = iL + def_maxlines;

							if(iL >= m.NL[fi] - def_maxlines){
								iLL = 0;
								fii = fi + 1;
								timeT[0] += time[0];
								time[0] = 0.0;
						
								fclose(dataFile);	
								sprintf(dataFilename, "%sbin", m.dataFilename[fii]);
								dataFile = fopen(dataFilename, "rb");

								if(dataFile == NULL){
									printf("Error: line list file not found: %s\n", dataFilename);
									return 0;
								}

								printf("Reading Line file %d of %d: %s\n", fii, fi1 - 1, dataFilename);
								printf("Number of lines: %lld\n", m.NL[fii]);
							}
							NL1 = min(def_maxlines, m.NL[fii] - iLL);

							lPart = (iLL + def_maxlines - 1) / def_maxlines;

							cudaEventRecord(ReadStart);
							printf("Reading Line file %d of %d; part %lld of %lld with %d lines\n", fii, fi1 - 1, lPart, (m.NL[fii] + def_maxlines - 1) / def_maxlines - 1, NL);
							readBufferCount = 0;
							rbvs = 0;	
						}

						cudaDeviceSynchronize();

/*
for(int i = 0; i < nlLimitsA; ++i){
	int ni = L.iiLimitsA1_h[i] - L.iiLimitsA0_h[i];
	if(ni > 0 || i < 10) printf("iiLimitsTotal A %d %lld %lld | %d\n", i, L.iiLimitsA0_h[i], L.iiLimitsA1_h[i], ni);
}
for(int i = 0; i < nlLimitsA; ++i){
	int ni = L.iiLimitsAL1_h[i] - L.iiLimitsAL0_h[i];
	if(ni > 0 || i < 10) printf("iiLimitsTotal AL %d %lld %lld | %d\n", i, L.iiLimitsAL0_h[i], L.iiLimitsAL1_h[i], ni);
}
for(int i = 0; i < nlLimitsA; ++i){
	int ni = L.iiLimitsAR1_h[i] - L.iiLimitsAR0_h[i];
	if(ni > 0 || i < 10) printf("iiLimitsTotal AR %d %lld %lld | %d\n", i, L.iiLimitsAR0_h[i], L.iiLimitsAR1_h[i], ni);
}
for(int i = 0; i < nlLimitsB; ++i){
	int ni = L.iiLimitsB1_h[i] - L.iiLimitsB0_h[i];
	if(ni > 0 || i < 10) printf("iiLimitsTotal B %d %lld %lld | %d\n", i, L.iiLimitsB0_h[i], L.iiLimitsB1_h[i], ni);
}
for(int i = 0; i < nlLimitsC; ++i){
	int ni = L.iiLimitsC1_h[i] - L.iiLimitsC0_h[i];
	if(ni > 0 || i < 10) printf("iiLimitsTotal C %d %lld %lld | %d\n", i, L.iiLimitsC0_h[i], L.iiLimitsC1_h[i], ni);
}
*/
						if(nTA > 0){

							//A
							const int nntt = 128;
							for(int il = 0; il < NL; il += def_KSn * def_nlB){ //loop over lines
								// *************************************
								// Call A Line kernels
								ntA = Line6A_Call(L, param, KS_d, x_d, il, NL, nntt, nkA, Nx, tuneStream[0], 10, 1);
								// *************************************
								if(param.doTuning == 1){
									if(ntA > 0 && ntA < 0.6 * ntAOld || ntA > 1.6 * ntAOld){
										ntAOld = ntA;
										int nkt;
										int nktt = nkA;
										float time0;
										for(int k = 0; k < 2; ++k){
											 
											for(int j = 0; j < 8; ++j){
												if(j == 0) nkt = nkA;
												if(j > 0){
													if(k == 0) nkt = nkt * 2;
													else nkt = nkt / 2;
												}
												if(nkt > 128 || nkt < 1) break;

												cudaEventRecord(TuneStart);

												Line6A_Call(L, param, KS_d, x_d, il, NL, nntt, nkt, Nx, tuneStream[1], 10, 0);

												cudaEventRecord(TuneStop);
												cudaEventSynchronize(TuneStop);
												cudaEventElapsedTime(&milliseconds, TuneStart, TuneStop);
												printf("Selftune A %d %d %d %d %g\n", il, ntA, ntAOld, nkt, milliseconds);

												if(j == 0 && k == 0) time0 = milliseconds;
												else{
													if(milliseconds > time0) break;
													else{
														nktt = nkt;
														time0 = milliseconds;
													}
												}
											}
										}
										nkA = nktt;
										printf("Selftune A %d\n", nktt);
									}
								}
							}
							cudaEventRecord(AEvent, tuneStream[0]);

							if((iL < m.NL[fi] - def_maxlines || fi < fi1 - 1) && iP == param.nP - 1){
								// **************************
								//Read the Line list A
								// **************************
								if(param.dataBase == 2 || param.dataBase >= 30){
									// 2 30 31 32
									for(int i = readBufferCount; i < NL1; i += readBufferSize){
										//check if the A kernels have finished, otherwise use host to read more data
										int ev =  cudaEventQuery(AEvent);
										if(ev == 0) break;
//printf("read A %d %d\n", readBufferCount, ev);

										er = readFileExomol(L, NL1, dataFile, readBuffer_h, readBuffer_d, readBufferSize, readBufferN, readBufferCount, rbvs, CStream);
										readBufferCount += readBufferSize;
										++rbvs;
									}
								}

							}

						}

						if(param.profile == 1){	//only for voigt profiles
							if(nTAL > 0){

								//AL
								const int nntt = 128;
								for(int il = 0; il < NL; il += def_KSn * def_nlB){ //loop over lines
									// *************************************
									// Call AL Line kernels
									ntAL = Line6A_Call(L, param, KS_d, x_d, il, NL, nntt, nkAL, Nx, tuneStream[0], 11, 1);
									// *************************************
									if(param.doTuning == 1){
										if(ntAL > 0 && ntAL < 0.6 * ntALOld || ntAL > 1.6 * ntALOld){
											ntALOld = ntAL;
											int nkt;
											int nktt = nkAL;
											float time0;
											for(int k = 0; k < 2; ++k){
												 
												for(int j = 0; j < 8; ++j){
													if(j == 0) nkt = nkAL;
													if(j > 0){
														if(k == 0) nkt = nkt * 2;
														else nkt = nkt / 2;
													}
													if(nkt > 128 || nkt < 1) break;

													cudaEventRecord(TuneStart);

													Line6A_Call(L, param, KS_d, x_d, il, NL, nntt, nkt, Nx, tuneStream[1], 11, 0);

													cudaEventRecord(TuneStop);
													cudaEventSynchronize(TuneStop);
													cudaEventElapsedTime(&milliseconds, TuneStart, TuneStop);
													printf("Selftune AL %d %d %d %d %g\n", il, ntAL, ntALOld, nkt, milliseconds);

													if(j == 0 && k == 0) time0 = milliseconds;
													else{
														if(milliseconds > time0) break;
														else{
															nktt = nkt;
															time0 = milliseconds;
														}
													}
												}
											}
											nkAL = nktt;
											printf("Selftune AL %d\n", nktt);
										}
									}

								}
								cudaEventRecord(ALEvent, tuneStream[0]);

								if((iL < m.NL[fi] - def_maxlines || fi < fi1 - 1) && iP == param.nP - 1){
									// **************************
									//Read the Line list AL
									// **************************
									if(param.dataBase == 2 || param.dataBase >= 30){
										//2 30 31 32
										for(int i = readBufferCount; i < NL1; i += readBufferSize){
											//check if the AL kernels have finished, otherwise use host to read more data
											int ev =  cudaEventQuery(ALEvent);
											if(ev == 0) break;
	//printf("read AL %d %d\n", readBufferCount, ev);

											er = readFileExomol(L, NL1, dataFile, readBuffer_h, readBuffer_d, readBufferSize, readBufferN, readBufferCount, rbvs, CStream);
											readBufferCount += readBufferSize;
											++rbvs;
										}
									}

								}

							}
							if(nTAR > 0){

								//AR
								const int nntt = 128;
								for(int il = 0; il < NL; il += def_KSn * def_nlB){ //loop over lines
									// *************************************
									// Call AR Line kernels
									ntAR = Line6A_Call(L, param, KS_d, x_d, il, NL, nntt, nkAR, Nx, tuneStream[0], 12, 1);
									// *************************************
									if(param.doTuning == 1){
										if(ntAR > 0 && ntAR < 0.6 * ntAROld || ntAR > 1.6 * ntAROld){
											ntAROld = ntAR;
											int nkt;
											int nktt = nkAR;
											float time0;
											for(int k = 0; k < 2; ++k){
												 
												for(int j = 0; j < 8; ++j){
													if(j == 0) nkt = nkAR;
													if(j > 0){
														if(k == 0) nkt = nkt * 2;
														else nkt = nkt / 2;
													}
													if(nkt > 128 || nkt < 1) break;

													cudaEventRecord(TuneStart);

													Line6A_Call(L, param, KS_d, x_d, il, NL, nntt, nkt, Nx, tuneStream[1], 12, 0);

													cudaEventRecord(TuneStop);
													cudaEventSynchronize(TuneStop);
													cudaEventElapsedTime(&milliseconds, TuneStart, TuneStop);
													printf("Selftune AR %d %d %d %d %g\n", il, ntAR, ntAROld, nkt, milliseconds);

													if(j == 0 && k == 0) time0 = milliseconds;
													else{
														if(milliseconds > time0) break;
														else{
															nktt = nkt;
															time0 = milliseconds;
														}
													}
												}
											}
											nkAR = nktt;
											printf("Selftune AR %d\n", nktt);
										}
									}
								}
								cudaEventRecord(AREvent, tuneStream[0]);

								if((iL < m.NL[fi] - def_maxlines || fi < fi1 - 1) && iP == param.nP - 1){
									// **************************
									//Read the Line list AR
									// **************************
									if(param.dataBase == 2 || param.dataBase >= 30){
										//2 30 31 32
										for(int i = readBufferCount; i < NL1; i += readBufferSize){
											//check if the AR kernels have finished, otherwise use host to read more data
											int ev =  cudaEventQuery(AREvent);
											if(ev == 0) break;
	//printf("read AR %d %d\n", readBufferCount, ev);

											er = readFileExomol(L, NL1, dataFile, readBuffer_h, readBuffer_d, readBufferSize, readBufferN, readBufferCount, rbvs, CStream);
											readBufferCount += readBufferSize;
											++rbvs;
										}
									}

								}
							}
//cudaDeviceSynchronize();

							if(nTB > 0){
								// B
								const int nntt2 = 128;
								for(int il = 0; il < NL; il += def_KSn * def_nlB){ //loop over lines
									// *************************************
									// Call B Line kernels
									ntB = Line6B_Call(L, param, KS_d, x_d, il, NL, nntt2, nkB, Nx, tuneStream[0], 1);
									// *************************************

									if(param.doTuning == 1){
										if(ntB > 0 && ntB < 0.6 * ntBOld || ntB > 1.6 * ntBOld){
											ntBOld = ntB;
											int nkt;
											int nktt = nkB;
											float time0;
											for(int k = 0; k < 2; ++k){
												 
												for(int j = 0; j < 8; ++j){
													if(j == 0) nkt = nkB;
													if(j > 0){
														if(k == 0) nkt = nkt * 2;
														else nkt = nkt / 2;
													}
													if(nkt > 128 || nkt < 1) break;

													cudaEventRecord(TuneStart);

													Line6B_Call(L, param, KS_d, x_d, il, NL, nntt2, nkt, Nx, tuneStream[1], 0);

													cudaEventRecord(TuneStop);
													cudaEventSynchronize(TuneStop);
													cudaEventElapsedTime(&milliseconds, TuneStart, TuneStop);
													printf("Selftune B %d %d %d %d %g\n", il, ntB, ntBOld, nkt, milliseconds);

													if(j == 0 && k == 0) time0 = milliseconds;
													else{
														if(milliseconds > time0) break;
														else{
															nktt = nkt;
															time0 = milliseconds;
														}

													}

												}
											}
											nkB = nktt;
											printf("Selftune B %d\n", nktt);
										}
									}
								}
								cudaEventRecord(BEvent, tuneStream[0]);

								if((iL < m.NL[fi] - def_maxlines || fi < fi1 - 1) && iP == param.nP - 1){
									// **************************
									//Read the Line list B
									// **************************
									if(param.dataBase == 2 || param.dataBase >= 30){
										//2 30 31 32
										for(int i = readBufferCount; i < NL1; i += readBufferSize){
											//check if the B kernels have finished, otherwise use host to read more data
											int ev =  cudaEventQuery(BEvent);
											if(ev == 0) break;
//printf("read B %d %d\n", readBufferCount, ev);

											er = readFileExomol(L, NL1, dataFile, readBuffer_h, readBuffer_d, readBufferSize, readBufferN, readBufferCount, rbvs, CStream);
											readBufferCount += readBufferSize;
											++rbvs;
										}
									}

								}
							}


							//C
							if(nTC > 0){

								//search higher order regimes of the Voigt profile
								const int nntt3 = 128;


								float a = (float)(M_PI * sqrt(-1.0 / log(def_TOLF * 0.5)));
								float b = (float)(1.0 / sqrt(M_PI));
								float c = (float)(2.0 * a / M_PI);
								for(int il = 0; il < NL; il += def_KSn * def_nlC){ //loop over lines
									// *************************************
									// Call C Line kernels
									ntC = Line6C_Call(L, param, KS_d, x_d, il, NL, nntt3, nkC, Nx, a, b, c, tuneStream[0],1);
									// *************************************

									if(param.doTuning == 1){
										if(ntC > 0 && ntC < 0.6 * ntCOld || ntC > 1.6 * ntCOld){
											ntCOld = ntC;
											int nkt;
											int nktt = nkC;
											float time0;
											for(int k = 0; k < 2; ++k){
												 
												for(int j = 0; j < 8; ++j){
													if(j == 0) nkt = nkC;
													if(j > 0){
														if(k == 0) nkt = nkt * 2;
														else nkt = nkt / 2;
													}
													if(nkt > 128 || nkt < 1) break;

													cudaEventRecord(TuneStart);

													Line6C_Call(L, param, KS_d, x_d, il, NL, nntt3, nkt, Nx, a, b, c, tuneStream[1],0);

													cudaEventRecord(TuneStop);
													cudaEventSynchronize(TuneStop);
													cudaEventElapsedTime(&milliseconds, TuneStart, TuneStop);
													printf("Selftune C %d %d %d %d %g\n", il, ntC, ntCOld, nkt, milliseconds);

													if(j == 0 && k == 0) time0 = milliseconds;
													else{
														if(milliseconds > time0) break;
														else{
															nktt = nkt;
															time0 = milliseconds;
														}

													}

												}
											}
											nkC = nktt;
											printf("Selftune C %d\n", nktt);
										}
									}
								}
							}
						} //end profile 1 
						//Add now all streams together
printf("Add streams A\n");
						error = cudaGetLastError();
						if(error != 0){
							printf("K error = %d = %s\n",error, cudaGetErrorString(error));
							return 0;
						}
						cudaEventRecord(KStop);

						if((iL < m.NL[fi] - def_maxlines || fi < fi1 - 1) && iP == param.nP - 1){
							// **************************
							//Read the Line list end
							// **************************
							if(param.dataBase == 2 || param.dataBase >= 30){
								//2 30 31 32
								for(int i = readBufferCount; i < NL1; i += readBufferSize){
									er = readFileExomol(L, NL1, dataFile, readBuffer_h, readBuffer_d, readBufferSize, readBufferN, readBufferCount, rbvs, CStream);
									readBufferCount += readBufferSize;
									++rbvs;
								}
							}

						}
						if((iL < m.NL[fi] - def_maxlines || fi < fi1 - 1) && iP == param.nP - 1){
					
							if(param.dataBase < 2 || param.dataBase == 3){
								er = readFile(param, m, part, L, param.qalphaL, NL1, dataFile, Sscale, meanMass);
								if(er == 0){
									return 0;
								}
							}
						
							printf("Reading Line file %d, part %lld complete\n", fii, lPart);	
					
							cudaEventRecord(ReadStop);
							cudaEventSynchronize(ReadStop);
							cudaEventElapsedTime(&milliseconds, ReadStart, ReadStop);
							time[0] += milliseconds * 0.001;
						
							printf("Time for input, %d %lld:        %g seconds\n", fii, lPart, time[0]);
							// **************************read iL + 1
							NL = NL1;

							for(int i = 0; i < def_rBs; ++i){
								cudaStreamSynchronize(CStream[i]);
							}
						}
						

						//wait until all KS streams are complete
						cudaDeviceSynchronize();
						//collect streams and store all KS_d into K_d
						//set KS_d to zero
						AddKStreams_kernel <<< (Nx + 511) / 512, 512 >>> (K_d + iP * Nx, KS_d, def_KSn, Nx);
printf("Add streams B\n");


						// *************************************
						//synchronize here only if no more data has to be read from the disk.
						//otherwise read data before synchronization
						cudaEventSynchronize(KStop);
						cudaEventElapsedTime(&milliseconds, KStart, KStop);

						time[2] += milliseconds * 0.001;
						printf("Time for K(x):  %d %lld %d:        %g seconds\n", fi, lPartOld, iP, time[2]);
			
						cudaDeviceSynchronize();
						error = cudaGetLastError();
						if(error != 0){
							printf("Kb error = %d = %s\n",error, cudaGetErrorString(error));
							return 0;
						}
						if(iL >= m.NL[fi] - def_maxlines && iP == param.nP - 1){	
							InfoFile = fopen(InfoFilename, "a");
							fprintf(InfoFile,"File %d of %d\n", fi, fi1);
							fprintf(InfoFile,"Number of lines: %lld\n", m.NL[fi]);
							fprintf(InfoFile,"Time for input:        %g seconds\n", timeOld);
							fprintf(InfoFile,"Time for Lines:        %g seconds\n", time[1]);
							fprintf(InfoFile,"Time for K(x):         %g seconds\n", time[2]);
							fclose(InfoFile);
						}

					} // End of pressure loop

				} // End of maxLines loop
			} // End of linefile loop
			if(fi1 > fi0){
				fclose(dataFile);
			}

		}

		if(param.dataBase < 2 || param.dataBase == 3){
			free_Line(L, param);
		}
		else{
			free2_Line(L, param);
		}

	} //end species loop

	printf("\n");
 	printf("Time for input total:  %g seconds\n", timeT[0]);
 	printf("Time for Lines total:  %g seconds\n", timeT[1]);
 	printf("Time for K(x) total:   %g seconds\n", timeT[2]);

	free(binBoundaries_h);
	cudaFree(binIndex_d);
	cudaFree(binBoundaries_d);	
	cudaFree(KS_d);

	cudaEventRecord(tt1, 0);
	for(int i = 0; i < def_KSn; ++i){
		cudaStreamDestroy(VStream[i]);
	}
	for(int i = 0; i < def_rBs; ++i){
		cudaStreamDestroy(CStream[i]);
	}
	for(int i = 0; i < 2; ++i){
		cudaStreamDestroy(tuneStream[i]);
	}
	for(int i = 0; i < 5; ++i){
		cudaStreamDestroy(nuLimitsStream[i]);
	}

	// ****************************
	// Write the full line profile
	// ****************************
	if(param.doStoreFullK == 1){
		FILE *OutFile;
		char OutFilename[300];
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
		char OutFilename[500];
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
		char OutFilename[300];
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
				float Kf = (float)(K_h[j] * unitScale);
				fwrite(&Kf, sizeof(float), 1, OutFile);
			}
		}
		fclose(OutFile);
	}
	if(param.doStoreFullK == -2){
		//read a binary file
		FILE *OutFile;
		char OutFilename[500];
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
				K_h[j] = (double)(Kf) / unitScale;
			}
			cudaMemcpy(K_d + iP * Nx, K_h, Nx * sizeof(double), cudaMemcpyHostToDevice);
		}
		fclose(OutFile);
	}
	// *******************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Write error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	cudaEventRecord(tt2, 0);
	cudaEventSynchronize(tt2);
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	time[3] += milliseconds * 0.001;

	printf("Time for write K(x):   %g seconds\n", time[3]);

	cudaEventRecord(tt1, 0);

	// **************************************
	// compute the Planck and Rosseland means
	// **************************************
	if(param.doMean > 0){
		
		double *Pmn_d;
		double *Rmn_d;

		cudaMalloc((void **) &Pmn_d, Nx * sizeof(double));
		cudaMalloc((void **) &Rmn_d, Nx * sizeof(double));
	
		double *means_h, *means_d;	
		means_h = (double*)malloc(4 * sizeof(double));
		cudaMalloc((void **) &means_d, 4 * sizeof(double));

		FILE *Out4File;
		char Out4Filename[300];

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
				fprintf(Out4File, "%.20g\n", means_h[2] * param.dnu);
				fprintf(Out4File, "%.20g\n", integral1);
				fprintf(Out4File, "%.20g\n", means_h[3] * param.dnu);
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
	cudaEventRecord(tt2, 0);
	cudaEventSynchronize(tt2);
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	time[4] += milliseconds * 0.001;

	printf("Time for mean K(x):    %g seconds\n", time[4]);


	cudaEventRecord(tt1, 0);


	// ***************************************
	// Do the sorting of K for all bins
	// ***************************************
	thrust::device_ptr<double> K_dt = thrust::device_pointer_cast(K_d);
	thrust::device_ptr<int> binKey_dt = thrust::device_pointer_cast(binKey_d);
	for(int iP = 0; iP < param.nP; ++iP){
		thrust::sort_by_key(K_dt + iP * Nx, K_dt + Nx + iP * Nx, binKey_dt);
		thrust::stable_sort_by_key(binKey_dt, binKey_dt + Nx, K_dt + iP * Nx);
	}
	cudaFree(binKey_d);
	// ****************************************

	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Sort error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	cudaEventRecord(tt2, 0);
	cudaEventSynchronize(tt2);
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	time[5] += milliseconds * 0.001;

	printf("Time for sort K(x):    %g seconds\n", time[5]);

	cudaEventRecord(tt1, 0);


	// *********************************
	// Prepare Resampling and do QR factorization, the same for all bins
	// this doesn't work with individual bins
	// *********************************

//size_t free_byte;
//size_t total_byte;
//cudaMemGetInfo( &free_byte, &total_byte );
//printf("***MEMORY %g %g %g\n", (double)(free_byte), (double)(total_byte), (double)(total_byte) - (double)(free_byte));
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
//printf("***MEMORY %g %g %g\n", (double)(free_byte), (double)(total_byte), (double)(total_byte) - (double)(free_byte));

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
		char Out3Filename[300];
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
	// **********************************
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Resampling error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	cudaEventRecord(tt2, 0);
	cudaEventSynchronize(tt2);
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	time[6] += milliseconds * 0.001;

	printf("Time for Resampling:   %g seconds\n", time[6]);

	cudaEventRecord(tt1, 0);


	// *****************************
	// Write K per bin output
	// *****************************
	if(param.doStoreK > 0){
		FILE *Out2File;
		char Out2Filename[300];
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
	// ******************************
	cudaDeviceSynchronize();
	error = cudaGetLastError();
	if(error != 0){
		printf("Write error = %d = %s\n",error, cudaGetErrorString(error));
		return 0;
	}
	cudaEventRecord(tt2, 0);
	cudaEventSynchronize(tt2);
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	time[7] += milliseconds * 0.001;

	printf("Time for write K(y):   %g seconds\n", time[7]);

	cudaEventRecord(tt1, 0);


	//set correction factor for simpsons rule needed for resampling
	SimpsonCoefficient();

	// *********************************
	// Calculate the Transmission function
	// *********************************
	if(param.doTransmission > 0 ){

		double *Tr_h, *Tr_d;
		Tr_h = (double*)malloc(param.nbins * param.nTr * sizeof(double));
		cudaMalloc((void **) &Tr_d, param.nbins * param.nTr * sizeof(double));

		FILE *Out3File;
		char Out3Filename[300];

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
	cudaEventRecord(tt2, 0);
	cudaEventSynchronize(tt2);
	cudaEventElapsedTime(&milliseconds, tt1, tt2);
	time[8] += milliseconds * 0.001;

	printf("Time for Transmission: %g seconds\n", time[8]);

	InfoFile = fopen(InfoFilename, "a");
	fprintf(InfoFile,"\n");
 	fprintf(InfoFile,"Time for input total:	 %g seconds\n", timeT[0]);
 	fprintf(InfoFile,"Time for Lines total:  %g seconds\n", timeT[1]);
 	fprintf(InfoFile,"Time for K(x) total: 	 %g seconds\n", timeT[2]);
	fprintf(InfoFile,"Time for write K(x):   %g seconds\n", time[3]);
	fprintf(InfoFile,"Time for mean K(x):    %g seconds\n", time[4]);
	fprintf(InfoFile,"Time for sort K(x):    %g seconds\n", time[5]);
	fprintf(InfoFile,"Time for Resampling:   %g seconds\n", time[6]);
	fprintf(InfoFile,"Time for write K(y):   %g seconds\n", time[7]);
	fprintf(InfoFile,"Time for Transmission: %g seconds\n", time[8]);
	fclose(InfoFile);	


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
