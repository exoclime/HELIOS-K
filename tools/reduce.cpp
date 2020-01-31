#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


int main(int argc, char*argv[]){

	char X[160];
	char expP[6];
	char outdir[160];
	char indir[160];
	int m = 1;

	double dnu = 1e-5;	//resolution in wavenumbers
	int Dnu = 1000;		//points in wavenumbers per file
	int reduce = 1;		//amount of reducing the data files

	double numin = 0;
	double numax = 0;

	int binary = 1;		//write binary files or not
	FILE *binFile;
	FILE *outFile;
	char binFilename[160];
	char outFilename[160];

	sprintf(X, "1500_p000");
	sprintf(outdir, ".");
	sprintf(indir, ".");

	int setname = 0;
	double Tmin = 0.0;
	double Tmax = 2900.0;
        //Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-b") == 0){
			binary = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-reduce") == 0){
			reduce = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dnu") == 0){
			dnu = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-Dnu") == 0){
			Dnu = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-numin") == 0){
			numin = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-numax") == 0){
			numax = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outdir") == 0){
			sprintf(outdir, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-indir") == 0){
			sprintf(indir, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-outname") == 0){
			sprintf(X, "%s", argv[i + 1]);
			setname = 1;
		}
		else if(strcmp(argv[i], "-Tmax") == 0){
			Tmax = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-Tmin") == 0){
			Tmin = atof(argv[i + 1]);
		}
		else{
			printf("error: console argument not known\n");
			return 0;
		}

	}

	if(-log10(dnu) < 0.0){
		printf("error, resolution too low\n");
		return 0;
	}

	int Nx = Dnu * (int)((1.0 + 0.5 * dnu) / dnu);
	//int Nx = (int)(1.0 / dnu * Dnu);
//printf("%d %d %.40g %.40g\n", Nx, Dnu, dnu, (1.0 + 0.5 * dnu) / dnu);
	float *K;
	K = (float*)malloc(Nx * sizeof(float));

	int extract = 1;

printf("extract %d, setname %d, Nx %d\n", extract, setname, Nx);

	for(int ip = 0; ip < 34; ++ip){

		if(ip == 0) sprintf(expP, "n800");
		if(ip == 1) sprintf(expP, "n766");
		if(ip == 2) sprintf(expP, "n733");
		if(ip == 3) sprintf(expP, "n700");
		if(ip == 4) sprintf(expP, "n666");
		if(ip == 5) sprintf(expP, "n633");
		if(ip == 6) sprintf(expP, "n600");
		if(ip == 7) sprintf(expP, "n566");
		if(ip == 8) sprintf(expP, "n533");
		if(ip == 9) sprintf(expP, "n500");
		if(ip == 10) sprintf(expP, "n466");
		if(ip == 11) sprintf(expP, "n433");
		if(ip == 12) sprintf(expP, "n400");
		if(ip == 13) sprintf(expP, "n366");
		if(ip == 14) sprintf(expP, "n333");
		if(ip == 15) sprintf(expP, "n300");
		if(ip == 16) sprintf(expP, "n266");
		if(ip == 17) sprintf(expP, "n233");
		if(ip == 18) sprintf(expP, "n200");
		if(ip == 19) sprintf(expP, "n166");
		if(ip == 20) sprintf(expP, "n133");
		if(ip == 21) sprintf(expP, "n100");
		if(ip == 22) sprintf(expP, "n066");
		if(ip == 23) sprintf(expP, "n033");
		if(ip == 24) sprintf(expP, "p000");
		if(ip == 25) sprintf(expP, "p033");
		if(ip == 26) sprintf(expP, "p066");
		if(ip == 27) sprintf(expP, "p100");
		if(ip == 28) sprintf(expP, "p133");
		if(ip == 29) sprintf(expP, "p166");
		if(ip == 30) sprintf(expP, "p200");
		if(ip == 31) sprintf(expP, "p233");
		if(ip == 32) sprintf(expP, "p266");
		if(ip == 33) sprintf(expP, "p300");


		for(int iT = 0; iT < 29; ++iT){ 

			int T;
			if(iT < 14) T = 50 + iT * 50;
			else if(iT < 14 + 8) T = 800 + (iT - 14) * 100;
			else T = 1700 + (iT - 14 - 8) * 200;

			if(T > Tmax) continue;
			if(T < Tmin) continue;
	
			if(extract == 1){
				if(binary == 1){
					if (setname == 0) sprintf(outFilename, "%s/Out_%05d_%05d_%05d_%s.bin", outdir, (int)(numin), (int)(numax), T, expP);
					else sprintf(outFilename, "%s/Out_%s.bin", outdir, X);
					outFile = fopen(outFilename, "wb");
				}
				else{
					if (setname == 0)sprintf(outFilename, "%s/Out_%05d_%05d_%05d_%s.dat", outdir, (int)(numin), (int)(numax), T, expP);
					else sprintf(outFilename, "%s/Out_%s.dat", outdir, X);
					outFile = fopen(outFilename, "w");
				}
			}
			
			for(int nu = (int)(numin) / Dnu; nu < (int)((numax + Dnu - 1)) / Dnu; ++nu){
				sprintf(binFilename, "%s/Out_%05d_%05d_%05d_%s.bin", indir, nu * Dnu, (nu + 1) * Dnu, T, expP);
				printf("%s    ", binFilename);

				binFile = fopen(binFilename, "rb");
				if(extract == 0){
					if(binary == 1){
						if (setname == 0) sprintf(outFilename, "%s/Out_%05d_%05d_%05d_%s.bin", outdir, nu * Dnu, (nu + 1) * Dnu, T, expP);
						else sprintf(outFilename, "%s/Out_%s.bin", outdir, X);
						outFile = fopen(outFilename, "wb");
					}
					else{
						if (setname == 0) sprintf(outFilename, "%s/Out_%05d_%05d_%05d_%s.dat", outdir, nu * Dnu, (nu + 1) * Dnu, T, expP);
						else sprintf(outFilename, "%s/Out_%s.dat", outdir, X);
						outFile = fopen(outFilename, "w");
					}
				}
				printf("%s\n", outFilename);
				fread(K, sizeof(float), Nx, binFile);
				for(int i = 0; i < Nx; ++i){
					if(i % reduce == 0 && nu * Dnu + i * dnu >= numin && nu * Dnu + i * dnu< numax){
						if(binary == 1){
							fwrite(&K[i], sizeof(float), 1, outFile);	
						}
						else{
							fprintf(outFile,"%.20g %.20g\n", nu * Dnu + i * dnu, K[i]);
						}
					}
				}
				fclose(binFile);
				if(extract == 0) fclose(outFile);
			}
			if(extract == 1) fclose(outFile);
		}
	}
	return 0;
}
