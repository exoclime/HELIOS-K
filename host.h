// ****************************************
// This function computes the Chebyshev polynomials as a function of T
// n must by greater than 1.
// *****************************************
__host__ void Chebyshev(double T, double *Cheb, int n){
	Cheb[0] = 1.0;
	Cheb[1] = T;
	for(int i = 2; i < n; ++i){
		Cheb[i] = 2.0 * T * Cheb[i - 1] - Cheb[i - 2];
	}
}

// ****************************************
// This function reads the file q.dat and computes for each Isotopologue
// the corresponding Partition function Q(T)

//Author: Simon Grimm
//November 2014
// *****************************************
__host__ int ChebCoeff(char *qFilename, Partition &part, double T){

	//Calculate Chebychev polynomial
	double Cheb[NCheb];
	Chebyshev(T, Cheb, NCheb);

	//Read Chebychev Coefficients from q file	
	FILE *qFile;
	//Check size of q.dat file
	qFile = fopen(qFilename, "r");
	if(qFile == NULL){
		printf("Error: q.dat file not found\n");
		return 0;
	}
	int j;
	for(j = 0; j < 100; ++j){
		int id;
		double coeff;
		int er = fscanf (qFile, "%d", &id);
		if (er <= 0) break;
		for(int i = 0; i < NCheb; ++i){
			fscanf (qFile, "%lf", &coeff);
		}
	}
	fclose(qFile);
	part.n = j;
	
	part.id = (int*)malloc(j * sizeof(int));
	part.Q = (double*)malloc(j * sizeof(double));
	
	qFile = fopen(qFilename, "r");
	for(j = 0; j < 100; ++j){
		int id;
		double coeff;
		double Q = 0.0;
		int er = fscanf (qFile, "%d", &id);
		if (er <= 0) break;
		for(int i = 0; i < NCheb; ++i){
			fscanf (qFile, "%lf", &coeff);
			Q += coeff * Cheb[i];
		}
		part.id[j] = id;
		part.Q[j] = Q;
	}
	fclose(qFile);
	return 1;
}


__host__ int read_parameters(Param &param, char *paramFilename, int argc, char*argv[]){
	//Read parameters from param.dat file
	FILE *paramFile;
	param.dev = 0;
	paramFile = fopen(paramFilename, "r");
		char skip[160];
		char skip2[160];
		//read name
		fgets(skip, 7, paramFile);
		fgets(param.name, 160, paramFile);
		//fscanf (paramFile, "%s", &param.name);
		fgets(skip2, 3, paramFile);
		//read T
		fgets(skip, 4, paramFile);
		fscanf (paramFile, "%lf", &param.T);
		fgets(skip2, 3, paramFile);
		//read P
		fgets(skip, 4, paramFile);
		fscanf (paramFile, "%lf", &param.P);
		fgets(skip2, 3, paramFile);
		//read Molecule
		fgets(skip, 11, paramFile);
		fscanf (paramFile, "%d", &param.nMolecule);
		fgets(skip2, 3, paramFile);
		//read numin
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%lf", &param.numin);
		fgets(skip2, 3, paramFile);
		//read numax
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%lf", &param.numax);
		fgets(skip2, 3, paramFile);
		//read dnu
		fgets(skip, 6, paramFile);
		fscanf (paramFile, "%lf", &param.dnu);
		fgets(skip2, 3, paramFile);
		//read cutMode
		fgets(skip, 10, paramFile);
		fscanf (paramFile, "%d", &param.cutMode);
		fgets(skip2, 3, paramFile);
		//read cut
		fgets(skip, 6, paramFile);
		fscanf (paramFile, "%lf", &param.cut);
		fgets(skip2, 3, paramFile);
		//read doResampling
		fgets(skip, 15, paramFile);
		fscanf (paramFile, "%d", &param.doResampling);
		fgets(skip2, 3, paramFile);
		//read nC
		fgets(skip, 5, paramFile);
		fscanf (paramFile, "%d", &param.nC);
		fgets(skip2, 3, paramFile);
		//read doTransmission
		fgets(skip, 17, paramFile);
		fscanf (paramFile, "%d", &param.doTransmission);
		fgets(skip2, 3, paramFile);
		//read nTr
		fgets(skip, 6, paramFile);
		fscanf (paramFile, "%d", &param.nTr);
		fgets(skip2, 3, paramFile);
		//read dTr
		fgets(skip, 6, paramFile);
		fscanf (paramFile, "%lf", &param.dTr);
		fgets(skip2, 3, paramFile);
		if(param.nC > NmaxSample){
			printf("nC larger than NmaxSample, reduced to %d\n", NmaxSample);
			param.nC = NmaxSample;
		}
		//read doStoreFullK
		fgets(skip, 15, paramFile);
		fscanf (paramFile, "%d", &param.doStoreFullK);
		fgets(skip2, 3, paramFile);
		//read doStoreK
		fgets(skip, 12, paramFile);
		fscanf (paramFile, "%d", &param.doStoreK);
		fgets(skip2, 3, paramFile);
		//read nbins
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%d", &param.nbins);
		fgets(skip2, 3, paramFile);
		//read kmin
		fgets(skip, 7, paramFile);
		fscanf (paramFile, "%lf", &param.kmin);
		fgets(skip2, 3, paramFile);
		//read qalphaL
		fgets(skip, 10, paramFile);
		fscanf (paramFile, "%lf", &param.qalphaL);
		fgets(skip2, 3, paramFile);

	fclose(paramFile);

	//Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-name") == 0){
			sprintf(param.name, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-T") == 0){
			param.T = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-P") == 0){
			param.P = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-M") == 0){
			param.nMolecule = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-numin") == 0){
			param.numin = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-numax") == 0){
			param.numax = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dnu") == 0){
			param.dnu = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-cutM") == 0){
			param.cutMode = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-cut") == 0){
			param.cut = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dR") == 0){
			param.doResampling = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-nC") == 0){
			param.nC = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dT") == 0){
			param.doTransmission = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-nTr") == 0){
			param.nTr = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dTr") == 0){
			param.dTr = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dSF") == 0){
			param.doStoreFullK = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dSS") == 0){
			param.doStoreK = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-nbins") == 0){
			param.nbins = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-kmin") == 0){
			param.kmin = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-q") == 0){
			param.qalphaL = atof(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dev") == 0){
			param.dev = atoi(argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}

	}
	return 1;
}

__host__ int readFile(Molecule &m, Partition &part, Line &L, double qalphaL){
	FILE *dataFile;
	dataFile  = fopen(m.dataFilename, "r");
	if(dataFile == NULL){
		printf("Error: line list file not found\n");
		return 0;
	}
	//read line list file		

	char c1[3];
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
	
	char skip[5];

//int count[40];
//for(int cc = 0; cc < 40; ++cc){
//count[cc] = 0;	
//}
	for(int i = 0; i < m.NL; ++i){
	
		fgets(skip, 1, dataFile);
		//fgets(c1, 3, dataFile);
		//fgets(c2, 2, dataFile);
		fgets(c1, 4, dataFile);		//Use combined notation for Id (AFGL and molecule + abundance number
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
		
		L.nu_h[i] = strtod(c3, NULL);		
		L.S_h[i] = strtod(c4, NULL);		
		L.A_h[i] = strtod(c5, NULL);		
		L.delta_h[i] = strtod(c10, NULL);
		L.EL_h[i] = strtod(c8, NULL);		
		
		double gammaAir = strtod(c6, NULL);
		double gammaSelf = strtod(c7, NULL);
		L.alphaL_h[i] = (1.0 - qalphaL) * gammaAir + qalphaL * gammaSelf;
		L.n_h[i] = strtod(c9, NULL);
		L.alphaD_h[i] = L.n_h[i];
		
		int id= std::atoi(c1);
		int idAFGL;
//count[0] += 1;
		//Assign the Isotopologue properties
		for(int j = 0; j < m.nISO; ++j){
			if(id == m.ISO[j].id){
				L.mass_h[i] = m.ISO[j].m;
				L.Q_h[i] = m.ISO[j].Q;
				idAFGL = m.ISO[j].AFGL;
			}
		}
		double Q;
		//Assign the Partition function
		for(int j = 0; j < part.n; ++j){
			if(idAFGL == part.id[j]){
				Q = part.Q[j];
			}
		}
		L.Q_h[i] /= exp(Q);
		
//if(i < 10) printf("%d %d %d %g %g\n", i, id, idAFGL, exp(Q), Q_h[i]);
	}
//for(int cc = 0; cc < 40; ++cc){
//printf("%d %d\n", cc, count[cc]);	
//}
	fclose(dataFile);
	return 1;
}


__host__ void Alloc_Line(Line &L, Molecule &m){
	L.nu_h = (double*)malloc(m.NL * sizeof(double));
	L.S_h = (double*)malloc(m.NL * sizeof(double));
	L.A_h = (double*)malloc(m.NL * sizeof(double));
	L.delta_h = (double*)malloc(m.NL * sizeof(double));
	L.EL_h = (double*)malloc(m.NL * sizeof(double));
	L.alphaL_h = (double*)malloc(m.NL * sizeof(double));
	L.alphaD_h = (double*)malloc(m.NL * sizeof(double));
	L.n_h = (double*)malloc(m.NL * sizeof(double));
	L.mass_h = (double*)malloc(m.NL * sizeof(double));
	L.Q_h = (double*)malloc(m.NL * sizeof(double));
	L.ID_h = (int*)malloc(m.NL * sizeof(int));

	cudaMalloc((void **) &L.nu_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.S_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.A_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.delta_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.EL_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.alphaL_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.alphaD_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.n_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.mass_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.Q_d, m.NL * sizeof(double));
	cudaMalloc((void **) &L.ID_d, m.NL * sizeof(int));
}

__host__ void Copy_Line(Line &L, Molecule &m){
	cudaMemcpy(L.nu_d, L.nu_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.S_d, L.S_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.A_d, L.A_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.delta_d, L.delta_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.EL_d, L.EL_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.alphaL_d, L.alphaL_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.alphaD_d, L.alphaD_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.n_d, L.n_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.mass_d, L.mass_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.Q_d, L.Q_h, m.NL * sizeof(double), cudaMemcpyHostToDevice);
}

__host__ void free_Line(Line &L){
	free(L.nu_h);
	free(L.S_h);
	free(L.A_h);
	free(L.delta_h);
	free(L.EL_h);
	free(L.alphaL_h);
	free(L.alphaD_h);
	free(L.n_h);
	free(L.mass_h);
	free(L.Q_h);
	free(L.ID_h);

	cudaFree(L.nu_d);
	cudaFree(L.S_d);
	cudaFree(L.A_d);
	cudaFree(L.delta_d);
	cudaFree(L.EL_d);
	cudaFree(L.alphaL_d);
	cudaFree(L.alphaD_d);
	cudaFree(L.n_d);
	cudaFree(L.mass_d);
	cudaFree(L.Q_d);
	cudaFree(L.ID_d);


}
