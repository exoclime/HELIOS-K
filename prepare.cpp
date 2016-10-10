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
int readFile(Molecule &m, int fi){
	FILE *dataFile, *OutFile;
	char dataFilename[160], OutFilename[160];

	sprintf(dataFilename, "%spar", m.dataFilename[fi]);
	sprintf(OutFilename, "%sbin", m.dataFilename[fi]);
	printf("reading file %s\n", dataFilename);
	dataFile  = fopen(dataFilename, "r");
	OutFile = fopen(OutFilename, "wb");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", dataFilename);
		return 0;
	}
	//read line list file           

	char c1[4];
	//char c2[2];
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

	for(int i = 0; i < m.NL[fi]; ++i){
		fgets(skip, 1, dataFile);
		fgets(c1, 4, dataFile);         //Use combined notation for Id (AFGL and molecule + abundance number
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

		double nu = strtod(c3, NULL);
		double S = strtod(c4, NULL);
		double A = strtod(c5, NULL);
		double delta = strtod(c10, NULL);
		double EL = strtod(c8, NULL);
		double gammaAir = strtod(c6, NULL);
		double gammaSelf = strtod(c7, NULL);
		double n = strtod(c9, NULL);

		double mass;
		double Q0;

		//Assign the Isotopologue properties
		for(int j = 0; j < m.nISO; ++j){
			if(id == m.ISO[j].id){
				mass = m.ISO[j].m / def_NA;
				Q0 = m.ISO[j].Q;
			}
		}

		S /= mass;

		S = S * Q0;

		fwrite(&id, sizeof(int), 1, OutFile);
		fwrite(&nu, sizeof(double), 1, OutFile);
		fwrite(&S, sizeof(double), 1, OutFile);
		fwrite(&EL, sizeof(double), 1, OutFile);
		fwrite(&A, sizeof(double), 1, OutFile);
		fwrite(&delta, sizeof(double), 1, OutFile);
		fwrite(&gammaAir, sizeof(double), 1, OutFile);
		fwrite(&gammaSelf, sizeof(double), 1, OutFile);
		fwrite(&n, sizeof(double), 1, OutFile);
//if(i < 1000) printf("%d %.20g %.20g %.20g %.20g %g\n", id, nu, S, EL, A, mass);
	}
	fclose(dataFile);
	fclose(OutFile);

	return 1;
}


int main(int argc, char*argv[]){

        Param param;
        param.dev = 0;
	param.useHITEMP = 1;
	sprintf(param.path, "");

	Molecule m;
        m.NL[0] = 0;
        m.id = 5; //1 = H2O, 2 = CO, 5 = CO, 6 = CH4
        m.nISO = 0;


	//Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-HITEMP") == 0){
			param.useHITEMP = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-M") == 0){
			m.id = atoi(argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}
	}


	Init(m, param);

	for(int i = 0; i < m.nFiles; ++i){
		readFile(m, i);
	}

	return 0;
}
