#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "ISO.h"
#include <algorithm>

// ******************************************************************
//This Function reads the maximal states index from the Exomol states files
//Author: Simon Grimm
//October 2022
// *******************************************************************
int readStatesN(Molecule &m, int &nStates){

	nStates = 0;

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
	char c3[30];
	char c4[251];

	int id;

	for(int i = 0; i < m.nStates; ++i){
		//fgets(c1, 13, dataFile);
		//fgets(c2, 14, dataFile);
		//fgets(c3, 8, dataFile);
		
		fscanf(dataFile, "%s", c1);
		fscanf(dataFile, "%s", c2);
		fscanf(dataFile, "%s", c3);

		fgets(c4, 250, dataFile);
	
		id = atoi(c1);

		nStates = max(nStates, id);

	}
	++nStates;
	printf("states file max index: %d %d\n", m.nStates, nStates);
	return 1;
}


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
	char c4[251];
	char skip[151];

	int ii;
	for(int i = 0; i < m.nStates; ++i){
		//fgets(c1, 13, dataFile);
		//fgets(c2, 14, dataFile);
		//fgets(c3, 8, dataFile);
		
		
		fscanf(dataFile, "%s", c1);
		fscanf(dataFile, "%s", c2);
		fscanf(dataFile, "%s", c3);
	

		/*		
		//Use this for other energy level columns (TiO DUO)
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
	
		ii = atoi(c1);
		id[ii] = atoi(c1);
		E[ii] = strtod(c2, NULL);
		g[ii] = atoi(c3);

//if(Eb != E[i]) printf("EDuo %d %d %.20g %.20g\n", i, id[i], E[i], Eb);

if(i < 10 || i > m.nStates - 10 || i % 1000000 == 0) printf("s %d %d %.20g %d\n", i, id[i], E[i], g[i]);
	}
	printf("states file complete\n");
	return 1;
}


int readTransitions(Molecule &m, int *id, double *E, int *g, long long int nT, double mass, int fi){
	FILE *transFile, *OutFile;
	char transFilename[300], OutFilename[300];

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

			if(id[state0] == -1 || id[state1] == -1){
				printf("Error, state not valid %d %d\n", state0, state1);
				return 0;
			}
			
			EL = E[state0];
			nu = E[state1] - E[state0];
			gU = g[state1];
			nu1 = 0.0;
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

			if(id[state0] == -1 || id[state1] == -1){
				printf("Error, state not valid %d %d\n", state0, state1);
				return 0;
			}

			EL = E[state0];
			gU = g[state1];

//printf("%lld %d %d %g %g | %g %d\n", i, state0, state1, A, nu1, EL, gU);

			//use always the energy levels
			nu = E[state1] - E[state0];
			if(fabs(nu - nu1) > 1.0e-6){
if(i < 100) printf("nu %lld %.20g %.20g %.20g %.20g %.20g %d %d %d %.20g %.20g\n", i, nu, nu1, S, EL, A, gU, state0, state1, E[state0], E[state1]); 
				nuError = 1;
			}

if(i < 100 || i % 100000 == 0) printf("%s|%s|%s|%s\n", c1, c2, c3, c4); 
		
		}

		S = gU * A /(8.0 * M_PI * def_c * nu * nu * mass);
		if(nu == 0.0) S = 0.0;
if(i < 100 || i % 100000 == 0) printf("%lld %.20g %.20g %.20g %.20g %.20g %d %d %d %.20g %.20g\n", i, nu, nu1, S, EL, A, gU, state0, state1, E[state0], E[state1]); 
		fwrite(&nu, sizeof(double), 1, OutFile);
		fwrite(&S, sizeof(double), 1, OutFile);
		fwrite(&EL, sizeof(double), 1, OutFile);
		fwrite(&A, sizeof(double), 1, OutFile);
		numax = fmax(nu, numax);
		
		if(i > nT - 10 || feof(transFile)){
printf("%lld %.20g %.20g %.20g %.20g %d %d %d %.20g %.20g\n", i, nu, S, EL, A, gU, state0, state1, E[state0], E[state1]); 

		}

		if(feof(transFile)){
printf("\n %lld numax %g\n", i, numax);
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
	int *id, *g;
	double *E;

	int nStates;
	readStatesN(m, nStates);
	
	id = (int*)malloc(nStates * sizeof(int));
	E = (double*)malloc(nStates * sizeof(double));
	g = (int*)malloc(nStates * sizeof(int));

	for(int i = 0; i < nStates; ++i){
		id[i] = -1;
		E[i] = 0.0;
		g[i] = 0;
	}

	readStates(m, id, E, g);

	long long int nT = 20000000000LL;
	int er = 0;
	for(int i = 0; i < m.nFiles; ++i){
		printf("Molecule %s, file: %d, mass:%g\n", param.mParamFilename, i, mass);
		er = readTransitions(m, id, E, g, nT, mass, i);
		if(er <= 0){
			return 0;
		}
	}

	free(id);
	free(E);
	free(g);
	return 0;

}
