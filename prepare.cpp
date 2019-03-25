#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "define.h"
#include "ISO.h"


// ******************************************************************
//This Function reads the Hitran or Hitemp data files
//Author: Simon Grimm
//September 2016
// *******************************************************************
int readFile(Molecule &m, int fi, int dataBase, char *iso, int *nLines){
	FILE *dataFile, *OutFile;
	char dataFilename[160], OutFilename[160];

	sprintf(dataFilename, "%spar", m.dataFilename[fi]);
	sprintf(OutFilename, "%sbin", m.dataFilename[fi]);

	if(dataBase == 1){
		if(strcmp(iso, " 11") == 0){
			sprintf(OutFilename, "01_%05d-%05d_1H2-16O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 12") == 0){
			sprintf(OutFilename, "01_%05d-%05d_1H2-18O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 13") == 0){
			sprintf(OutFilename, "01_%05d-%05d_1H2-17O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 14") == 0){
			sprintf(OutFilename, "01_%05d-%05d_1H-2H-16O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 15") == 0){
			sprintf(OutFilename, "01_%05d-%05d_1H-2H-18O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 16") == 0){
			sprintf(OutFilename, "01_%05d-%05d_1H-2H-17O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}


		if(strcmp(iso, " 21") == 0){
			sprintf(OutFilename, "02_%05d-%05d_12C-16O2__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 22") == 0){
			sprintf(OutFilename, "02_%05d-%05d_13C-16O2__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 23") == 0){
			sprintf(OutFilename, "02_%05d-%05d_16O-12C-18O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 24") == 0){
			sprintf(OutFilename, "02_%05d-%05d_16O-12C-17O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 25") == 0){
			sprintf(OutFilename, "02_%05d-%05d_16O-13C-18O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 26") == 0){
			sprintf(OutFilename, "02_%05d-%05d_16O-13C-17O__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}
		if(strcmp(iso, " 27") == 0){
			sprintf(OutFilename, "02_%05d-%05d_12C-18O2__HITEMP2010.bin", m.fileLimit[fi], m.fileLimit[fi + 1]);
		}

		if(strcmp(iso, " 51") == 0){
			sprintf(OutFilename, "05_12C-16O__HITEMP2010.bin");
		}
		if(strcmp(iso, " 52") == 0){
			sprintf(OutFilename, "05_13C-16O__HITEMP2010.bin");
		}
		if(strcmp(iso, " 53") == 0){
			sprintf(OutFilename, "05_12C-18O__HITEMP2010.bin");
		}
		if(strcmp(iso, " 54") == 0){
			sprintf(OutFilename, "05_12C-17O__HITEMP2010.bin");
		}
		if(strcmp(iso, " 55") == 0){
			sprintf(OutFilename, "05_13C-18O__HITEMP2010.bin");
		}
		if(strcmp(iso, " 56") == 0){
			sprintf(OutFilename, "05_13C-17O__HITEMP2010.bin");
		}

		if(strcmp(iso, " 81") == 0){
			sprintf(OutFilename, "08_14N-16O__HITEMP2010.bin");
		}
		if(strcmp(iso, " 82") == 0){
			sprintf(OutFilename, "08_15N-16O__HITEMP2010.bin");
		}
		if(strcmp(iso, " 83") == 0){
			sprintf(OutFilename, "08_14N-18O__HITEMP2010.bin");
		}

		if(strcmp(iso, "131") == 0){
			sprintf(OutFilename, "13_16O-1H__HITEMP2010.bin");
		}
		if(strcmp(iso, "132") == 0){
			sprintf(OutFilename, "13_18O-1H__HITEMP2010.bin");
		}
		if(strcmp(iso, "133") == 0){
			sprintf(OutFilename, "13_16O-2H__HITEMP2010.bin");
		}
	}

	printf("reading file %s\n", dataFilename);
	printf("Output file %s\n", OutFilename);
	dataFile  = fopen(dataFilename, "r");
	OutFile = fopen(OutFilename, "wb");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", dataFilename);
		return 0;
	}
	//read line list file           

	char c1[3];
	char c2[2];
	char c3[13];
	char c4[11];
	char c5[11];
	char c6[6];
	char c7[6];
	char c8[11];
	char c9[5];
	char c10[9];
	char c11[16];
	char c12[16];
	char c13[16];
	char c14[16];
	char c15[7];
	char c16[13];
	char c17[2];
	char c18[8];
	char c19[8];

	char skip[6];
	int count = 0;	
	double numax = 0.0;

	for(int i = 0; i < m.NL[fi]; ++i){
		//fgets(skip, 1, dataFile);
		//fgets(cid, 4, dataFile);         //Use combined notation for Id (AFGL and molecule + abundance number
		fgets(c1, 3, dataFile);
		fgets(c2, 2, dataFile);
		fgets(c3, 13, dataFile);
		fgets(c4, 11, dataFile);
		fgets(c5, 11, dataFile);
		fgets(c6, 6, dataFile);
		fgets(c7, 6, dataFile);
		fgets(c8, 11, dataFile);
		fgets(c9, 5, dataFile);
		fgets(c10, 9, dataFile);

		fgets(c11, 16, dataFile);
		fgets(c12, 16, dataFile);
		fgets(c13, 16, dataFile);
		fgets(c14, 16, dataFile);
		fgets(c15, 7, dataFile);
		fgets(c16, 13, dataFile);
		fgets(c17, 2, dataFile);
		fgets(c18, 8, dataFile);
		fgets(c19, 8, dataFile);
		fgets(skip, 6, dataFile);

		int id = atoi(c1);
		char cid[4];
		sprintf(cid, "%2d%s", id, c2); 

		double nu = strtod(c3, NULL);
		double S = strtod(c4, NULL);
		double A = strtod(c5, NULL);
		double delta = strtod(c10, NULL);
		double EL = strtod(c8, NULL);
		double gammaAir = strtod(c6, NULL);
		double gammaSelf = strtod(c7, NULL);
		double n = strtod(c9, NULL);

		double mass;
		double Q0 = -1000.0;


		//Assign the Isotopologue properties
		for(int j = 0; j < m.nISO; ++j){
//if(i < 10) printf("%d|%s|%s|\n", id, cid, m.ISO[j].cid);
			if(strcmp(cid, m.ISO[j].cid) == 0){
			//if(cid == m.ISO[j].cid){
				mass = m.ISO[j].m / def_NA;
				Q0 = m.ISO[j].Q;
			}
		}
		if(Q0 < -800){
			if(iso == 0){
				printf("Error in assigning isotopologue indices\n");
				return 0;
			}
		}

		S /= mass;

		S = S * Q0;

		if(strcmp(iso, cid) == 0 || strcmp(iso, "") == 0){
			//fwrite(&id, sizeof(int), 1, OutFile);
			fwrite(&cid, 4*sizeof(char), 1, OutFile);
			fwrite(&nu, sizeof(double), 1, OutFile);
			fwrite(&S, sizeof(double), 1, OutFile);
			fwrite(&EL, sizeof(double), 1, OutFile);
			fwrite(&A, sizeof(double), 1, OutFile);
			fwrite(&delta, sizeof(double), 1, OutFile);
			fwrite(&gammaAir, sizeof(double), 1, OutFile);
			fwrite(&gammaSelf, sizeof(double), 1, OutFile);
			fwrite(&n, sizeof(double), 1, OutFile);
if(i < 10 || i > m.NL[fi] - 10) printf("%s | %s %.20g %.20g %.20g %.20g %g\n", iso, cid, nu, S, EL, A, mass);
			++count;
			numax = fmax(numax, nu);
		}
	}
	nLines[fi] = count;
printf("File %d Number of lines %d, numax %g\n", fi, count, numax);
	fclose(dataFile);
	fclose(OutFile);

	return 1;
}


int main(int argc, char*argv[]){

        Param param;
        param.dev = 0;
	param.dataBase = 0;
	param.path[0] = 0;
	param.mParamFilename[0] = 0;

	Molecule m;
        m.NL[0] = 0;
        m.id = 0;
        m.nISO = 0;
	char iso[4];
	iso[0] = 0;

	//Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-M") == 0){
			sprintf(param.mParamFilename, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-ISO") == 0){
			sprintf(iso, "%3s", argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}
	}


printf("ISO |%s|\n", iso);

	char qFilename[15][160];
	Init(m, param, qFilename);

	int *nLines;
	nLines = (int*)malloc(m.nFiles * sizeof(int));
	for(int i = 0; i < m.nFiles; ++i){
		readFile(m, i, param.dataBase, iso, nLines);
	}
	for(int i = 0; i < m.nFiles; ++i){
		printf("%d\n", nLines[i]);
	}

	return 0;
}
