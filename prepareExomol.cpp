#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "ISO.h"
#include <algorithm>



// ******************************************************************
//This Function reads the Exomol data files
//Author: Simon Grimm
//September 2016
// *******************************************************************
int readStates(Molecule &m, int *id, double *E, int *g){

	FILE *dataFile;
	char statesFilename[180];
	sprintf(statesFilename, "%s.states", m.mName);
	dataFile = fopen(statesFilename, "r");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", statesFilename);
		return 0;
	}
	char c1[30];
	char c2[30];
	char c2b[30];
	char c3[30];
	char c4[151];
	char skip[151];

	for(int i = 0; i < m.nStates; ++i){
		//fgets(c1, 13, dataFile);
		//fgets(c2, 14, dataFile);
		//fgets(c3, 8, dataFile);
		
		
		fscanf(dataFile, "%s", c1);
		fscanf(dataFile, "%s", c2);
		fscanf(dataFile, "%s", c3);
	

		/*		
		//Use this for other energy level columns (TIO DUO)
		fscanf(dataFile, "%s", c1);
		fscanf(dataFile, "%s", c2b);
		fscanf(dataFile, "%s", c3);

		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", skip);
		fscanf(dataFile, "%s", c2);

		double Eb = strtod(c2b, NULL);
		*/

		fgets(c4, 250, dataFile);
	
		id[i] = atoi(c1);
		E[i] = strtod(c2, NULL);
		g[i] = atoi(c3);

//if(Eb != E[i]) printf("EDuo %d %d %.20g %.20g\n", i, id[i], E[i], Eb);

if(i < 10 || i > m.nStates - 10 || i % 1000000 == 0) printf("s %d %.20g %d\n", id[i], E[i], g[i]);
	}
	printf("states file complete\n");
	return 1;
}


int readTransitions(Molecule &m, int *id, double *E, int *g, long long int nT, double mass, int fi){
	FILE *transFile, *OutFile;
	char transFilename[160], OutFilename[160];

	printf("reading file %d\n", fi);
	sprintf(transFilename, "%strans", m.dataFilename[fi]);
	sprintf(OutFilename, "%sbin", m.dataFilename[fi]);
	transFile = fopen(transFilename, "r");
	OutFile = fopen(OutFilename, "wb");
	if(transFile == NULL){
		printf("Error: line list file not found %s\n", transFilename);
		return 0;
	}

	char c1[15];
	char c2[15];
	char c3[15];
	char c4[25];

	char skip[100];
	
	int state0;
	int state1;
	double A;

	double nu, nu1;
	double EL;
	int gU;
	double S;

	double numax = 0.0;
	int nuError = 0;

	for(long long int i = 0LL; i < nT + 1; ++i){

		if(m.ntcol == 3){
			//fgets(c1, 13, transFile);
			//fgets(c2, 14, transFile);
			//fgets(c3, 12, transFile);
			//fgets(skip, 2, transFile);
			fscanf(transFile, "%s", c1);
			fscanf(transFile, "%s", c2);
			fscanf(transFile, "%s", c3);
			
if(i < 100 || i % 100000 == 0) printf("||%s|%s|%s||\n", c1, c2, c3); 

			state1 = atoi(c1);
			state0 = atoi(c2);
			A = strtod(c3, NULL);
			
			EL = E[state0 - 1];
			nu = E[state1 - 1] - E[state0 - 1];
			gU = g[state1 - 1];
		}
		if(m.ntcol == 4){
			fscanf(transFile, "%s", c1);
			fscanf(transFile, "%s", c2);
			fscanf(transFile, "%s", c3);
			fscanf(transFile, "%s", c4);

			state1 = atoi(c1);
			state0 = atoi(c2);
			A = strtod(c3, NULL);
			nu1 = strtod(c4, NULL);	

			EL = E[state0 - 1];
			gU = g[state1 - 1];

			//use always the energy levels
			nu = E[state1 - 1] - E[state0 - 1];
			if(fabs(nu - nu1) > 1.0e-6){
if(i < 100) printf("nu %lld %.20g %.20g %.20g %.20g %.20g %d %d %d %.20g %.20g\n", i, nu, nu1, S, EL, A, gU, state0, state1, E[state0 - 1], E[state1 - 1]); 
				nuError = 1;
			}

if(i < 100 || i % 100000 == 0) printf("%s|%s|%s|%s\n", c1, c2, c3, c4); 
		
			EL = E[state0 - 1];
			gU = g[state1 - 1];
		}

		S = gU * A /(8.0 * M_PI * def_c * nu * nu * mass);
		if(nu == 0.0) S = 0.0;
if(i < 100 || i % 100000 == 0) printf("%lld %.20g %.20g %.20g %.20g %d %d %d %.20g %.20g\n", i, nu, S, EL, A, gU, state0, state1, E[state0 - 1], E[state1 - 1]); 
		fwrite(&nu, sizeof(double), 1, OutFile);
		fwrite(&S, sizeof(double), 1, OutFile);
		fwrite(&EL, sizeof(double), 1, OutFile);
		fwrite(&A, sizeof(double), 1, OutFile);
		numax = fmax(nu, numax);
		if(feof(transFile)){
printf("%lld %.20g %.20g %.20g %.20g %d %d %d %.20g %.20g\n", i, nu, S, EL, A, gU, state0, state1, E[state0 - 1], E[state1 - 1]); 
printf("\n %g\n", numax);
			break;
		}

	}
	fclose(transFile);
	fclose(OutFile);

	if(nuError == 1){
		printf("************* nu in 4th column not correct ******************\n");
	}

	return 1;

}

int main(int argc, char*argv[]){

        Param param;
        param.dev = 0;
	param.path[0] = 0;
	param.mParamFilename[0] = 0;

	Molecule m;
        m.NL[0] = 0;
        m.id = 0;
        m.nISO = 0;

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

	double mass = m.ISO[0].m;	//Molar Mass (g)
	mass /= def_NA;
	int *id, *g;
	double *E;
	
	id = (int*)malloc(m.nStates * sizeof(int));
	E = (double*)malloc(m.nStates * sizeof(double));
	g = (int*)malloc(m.nStates * sizeof(int));

	readStates(m, id, E, g);

	long long int nT = 20000000000LL;
	for(int i = 0; i < m.nFiles; ++i){
		printf("id %d, file: %d, mass:%g\n", m.id, i, mass);
		readTransitions(m, id, E, g, nT, mass, i);
	}

	free(id);
	free(E);
	free(g);
	return 0;

}
