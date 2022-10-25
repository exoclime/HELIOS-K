#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "ISO.h"
#include <algorithm>



int readTransitions(Molecule &m, long long int nT, double mass){

	FILE *transFile, *OutFile;
	char transFilename[300], OutFilename[300];

	sprintf(transFilename, "%ssuper", m.dataFilename[0]);
	sprintf(OutFilename, "%sbin", m.dataFilename[0]);

	transFile = fopen(transFilename, "r");
	OutFile = fopen(OutFilename, "wb");
	if(transFile == NULL){
		printf("Error: line list file not found %s\n", transFilename);
		return 0;
	}

	double nu;
	double S;

	double numax = 0.0;

	for(long long int i = 0LL; i < nT + 1; ++i){

		fscanf(transFile, "%lf", &nu);
		fscanf(transFile, "%lf", &S);
		
if(i < 100 || i % 100000 == 0 || i > nT - 10) printf("%.8g %g\n", nu, S); 


		//S = gU * A /(8.0 * M_PI * def_c * nu * nu * mass);
		S /= mass;

		if(nu == 0.0) S = 0.0;

		fwrite(&nu, sizeof(double), 1, OutFile);
		fwrite(&S, sizeof(double), 1, OutFile);
		numax = fmax(nu, numax);

		if(feof(transFile)){
printf("\n numax %lld %lld %g\n", nT, i, numax);
			break;
		}

	}
	fclose(transFile);
	fclose(OutFile);


	return 1;

}

int main(int argc, char*argv[]){

        Param param;
        param.dev = 0;
	param.path[0] = 0;
	param.mParamFilename[0] = 0;

	Molecule m;

	//Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-M") == 0){
			sprintf(param.mParamFilename, "%s", argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}
	}

	char qFilename[15][160];
	Init(m, param, qFilename);

	double mass = m.ISO[0].m;	//Molar Mass (g/mol)
	mass /= def_NA;			//mass in g

	long long int nT = m.NLmax;//20000000000LL;
	printf("Molecule %s, mass:%g\n", param.mParamFilename, mass);
	readTransitions(m, nT, mass);

	return 0;

}
