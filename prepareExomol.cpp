#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "define.h"
#include "ISO.h"
#include <algorithm>

//Works only for Methane and H2O for now


// ******************************************************************
//This Function reads the Exomol data files
//Author: Simon Grimm
//September 2016
// *******************************************************************
int readStates(Molecule &m, int *id, double *E, int *g){

	FILE *dataFile;
	char statesFilename[160];
	if(m.id == 1) {
		sprintf(statesFilename, "1H2-16O__BT2.states");
	}
	if(m.id == 6) {
		sprintf(statesFilename, "12C-1H4__YT10to10.states");
	}
	if(m.id == 11) {
		sprintf(statesFilename, "14N-1H3__BYTe.states");
	}
	if(m.id == 23) {
		sprintf(statesFilename, "1H-12C-14N__Harris.states");
	}
	if(m.id == 31) {
		sprintf(statesFilename, "1H2-32S__AYT2.states");
	}
	dataFile = fopen(statesFilename, "r");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", statesFilename);
		return 0;
	}
	if(m.id == 1){
		char c1[13];
		char c2[15];
		char c3[9];
		char c4[58];
	
		for(int i = 0; i < m.nStates; ++i){
			fgets(c1, 13, dataFile);
			fgets(c2, 15, dataFile);
			fgets(c3, 9, dataFile);
			fgets(c4, 58, dataFile);

			id[i] = atoi(c1);
			E[i] = strtod(c2, NULL);
			g[i] = atoi(c3);
if(i < 10) printf("%d %d %.40g %d\n", i, id[i], E[i], g[i]);
			if(i % 1000000 == 0) printf("read states line %d\n", i);
		}
	}
	if(m.id == 6){
		char c1[13];
		char c2[15];
		char c3[9];
		char c4[135];
	
		for(int i = 0; i < m.nStates; ++i){
			fgets(c1, 13, dataFile);
			fgets(c2, 15, dataFile);
			fgets(c3, 9, dataFile);
			fgets(c4, 135, dataFile);

			id[i] = atoi(c1);
			E[i] = strtod(c2, NULL);
			g[i] = atoi(c3);
if(i < 10) printf("s %d %.20g %d\n", id[i], E[i], g[i]);
			if(i % 1000000 == 0) printf("read states line %d\n", i);
		}
	}
	if(m.id == 11){
		char c1[14];
		char c2[15];
		char c3[9];
		char c4[110];
	
		for(int i = 0; i < m.nStates; ++i){
			fgets(c1, 13, dataFile);
			fgets(c2, 14, dataFile);
			fgets(c3, 8, dataFile);
			fgets(c4, 109, dataFile);

			id[i] = atoi(c1);
			E[i] = strtod(c2, NULL);
			g[i] = atoi(c3);
if(i < 10) printf("s %d %.20g %d\n", id[i], E[i], g[i]);
			if(i % 1000000 == 0) printf("read states line %d\n", i);
		}
	}
	if(m.id == 23){
		char c1[14];
		char c2[15];
		char c3[9];
		char c4[110];
	
		for(int i = 0; i < m.nStates; ++i){
			fgets(c1, 13, dataFile);
			fgets(c2, 14, dataFile);
			fgets(c3, 8, dataFile);
			fgets(c4, 109, dataFile);

			id[i] = atoi(c1);
			E[i] = strtod(c2, NULL);
			g[i] = atoi(c3);
if(i < 10) printf("s %d %.20g %d\n", id[i], E[i], g[i]);
			if(i % 1000000 == 0) printf("read states line %d\n", i);
		}
	}
	if(m.id == 31){
		char c1[14];
		char c2[15];
		char c3[9];
		char c4[59];
	
		for(int i = 0; i < m.nStates; ++i){
			fgets(c1, 13, dataFile);
			fgets(c2, 14, dataFile);
			fgets(c3, 8, dataFile);
			fgets(c4, 58, dataFile);

			id[i] = atoi(c1);
			E[i] = strtod(c2, NULL);
			g[i] = atoi(c3);
if(i < 10) printf("s %d %.20g %d\n", id[i], E[i], g[i]);
			if(i % 1000000 == 0) printf("read states line %d\n", i);
		}
	}
	printf("states file complete\n");
	return 1;
}


int readTransitions(Molecule &m, int *id, double *E, int *g, int nT, double mass, int fi){
	FILE *transFile, *OutFile;
	char transFilename[160], OutFilename[160];

	printf("reading file %d ", fi);
	sprintf(transFilename, "%strans", m.dataFilename[fi]);
	sprintf(OutFilename, "%sbin", m.dataFilename[fi]);
	transFile = fopen(transFilename, "r");
	OutFile = fopen(OutFilename, "wb");
	if(transFile == NULL){
		printf("Error: line list file not found %s\n", transFilename);
		return 0;
	}

	char c1[13];
	char c2[14];
	char c3[13];
	char c4[19];

	char skip[3];
	
	int state0;
	int state1;
	double A;

	double nu;
	double EL;
	int gU;
	double S;

	for(int i = 0; i < nT + 1; ++i){

		if(m.ntcol == 3){
			fgets(c1, 13, transFile);
			fgets(skip, 1, transFile);
			fgets(c2, 14, transFile);
			fgets(skip, 1, transFile);
			fgets(c3, 13, transFile);
			fgets(skip, 2, transFile);
			

			state1 = atoi(c1);
			state0 = atoi(c2);
			A = strtof(c3, NULL);
			
			EL = E[state0 - 1];
			nu = E[state1 - 1] - E[state0 - 1];
			gU = g[state1 - 1];
		}
		if(m.ntcol == 4){
			fgets(c1, 13, transFile);
			fgets(skip, 1, transFile);
			fgets(c2, 14, transFile);
			fgets(skip, 1, transFile);
			fgets(c3, 12, transFile);
			fgets(skip, 6, transFile);
			fgets(c4, 13, transFile);
			fgets(skip, 2, transFile);


			

			state1 = atoi(c1);
			state0 = atoi(c2);
			A = strtof(c3, NULL);
			nu = strtof(c4, NULL);	
//printf("%s|%s|%s|%s\n", c1, c2, c3, c4); 
		
			EL = E[state0 - 1];
			gU = g[state1 - 1];
		}

		S = gU * A /(8.0 * M_PI * def_c * nu * nu * mass);
		if(nu == 0.0) S = 0.0;
if(i < 10 || i % 10000 == 0) printf("%d %.20g %.20g %.20g %.20g %d %d %.20g %.20g\n", i, nu, S, EL, A, state0, state1, E[state0 - 1], E[state1 - 1]); 
		fwrite(&nu, sizeof(double), 1, OutFile);
		fwrite(&S, sizeof(double), 1, OutFile);
		fwrite(&EL, sizeof(double), 1, OutFile);
		fwrite(&A, sizeof(double), 1, OutFile);
		if(feof(transFile)){
printf("%d %.10g %d %d %s\n", i, nu, state0, state1, c4);
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
	param.useHITEMP = 2;
	sprintf(param.path, "");

	Molecule m;
        m.NL[0] = 0;
        m.id = 23; //1 = H2O, 2 = CO, 5 = CO, 6 = CH4
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

	char qFilename[160];
	Init(m, param, qFilename);

	double mass = m.ISO[0].m;	//Molar Mass (g)
	mass /= def_NA;
	int *id, *g;
	double *E;
	
	id = (int*)malloc(m.nStates * sizeof(int));
	E = (double*)malloc(m.nStates * sizeof(double));
	g = (int*)malloc(m.nStates * sizeof(int));

	readStates(m, id, E, g);

	int nT = 200000000;
	for(int i = 0; i < m.nFiles; ++i){
		readTransitions(m, id, E, g, nT, mass, i);
	}

	free(id);
	free(E);
	free(g);
	return 0;

}
