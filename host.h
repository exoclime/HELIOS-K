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
	for(j = 0; j < 1000; ++j){
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
	for(j = 0; j < 1000; ++j){
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

// ****************************************
// This function reads the partition function file from ExoMol
// and interpolates Q

//Author: Simon Grimm
//August 2016
// *****************************************
__host__ int readPartitionExomol(char *qFilename, Partition &part, double T){

	//Read Chebychev Coefficients from q file	
	FILE *qFile;
	qFile = fopen(qFilename, "r");
	if(qFile == NULL){
		printf("Error: partition file not found %s\n", qFilename);
		return 0;
	}
	double T0, T1;
	double q0, q1, q;
	T1 = 0;
	q1 = 1.0;
	q = 1.0;
	int er = 0;
	for(int j = 0; j < 150; ++j){
		T0 = T1;
		q0 = q1;
		er = fscanf (qFile, "%lf", &T1);
		er = fscanf (qFile, "%lf", &q1);
		if (er <= 0) break;
		if(T0 < T && T1 >= T){
			double tt = (T - T0) / (T1 - T0);
			q = (q1 - q0) * tt + q0;
			break;
		}
	}
	fclose(qFile);
	part.n = 1;
	
	part.id = (int*)malloc(sizeof(int));
	part.Q = (double*)malloc(sizeof(double));


	part.id = 0;
	part.Q[0] = q;
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
		fscanf (paramFile, "%s", &param.name);
		fgets(skip2, 3, paramFile);
		//read T
		fgets(skip, 4, paramFile);
		fscanf (paramFile, "%lf", &param.T);
		fgets(skip2, 3, paramFile);
		//read P
		fgets(skip, 4, paramFile);
		fscanf (paramFile, "%lf", &param.P);
		fgets(skip2, 3, paramFile);
		//read PFile
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%s", &param.PFilename);
		fgets(skip2, 3, paramFile);
		//read HITEMP
		fgets(skip, 12, paramFile);
		fscanf (paramFile, "%d", &param.useHITEMP);
		fgets(skip2, 3, paramFile);
		//read Molecule
		fgets(skip, 11, paramFile);
		fscanf (paramFile, "%d", &param.nMolecule);
		fgets(skip2, 3, paramFile);
		//read ciaSystem
		fgets(skip, 12, paramFile);
		fscanf (paramFile, "%s", &param.ciaSystem);
		fgets(skip2, 3, paramFile);
		//read path
		fgets(skip, 13, paramFile);
		fscanf (paramFile, "%s", &param.path);
		fgets(skip2, 3, paramFile);
		if(strcmp(param.path, "numin") == 0){
			sprintf(param.path, "");
		
			//fgets(skip, 3, paramFile);
			fscanf (paramFile, "%lf", &param.numin);
			fgets(skip2, 3, paramFile);
		}
		else{
		//read numin
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%lf", &param.numin);
		fgets(skip2, 3, paramFile);
		}
		//read numax
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%lf", &param.numax);
		fgets(skip2, 3, paramFile);
		//read dnu
		fgets(skip, 6, paramFile);
		fscanf (paramFile, "%lf", &param.dnu);
		fgets(skip2, 3, paramFile);
		//read Nxb
		fgets(skip, 14, paramFile);
		fscanf (paramFile, "%d", &param.Nxb);
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
		//read pathK
		fgets(skip, 10, paramFile);
		fscanf (paramFile, "%s", &param.pathK);
		fgets(skip2, 3, paramFile);
		if(strcmp(param.pathK, "doStoreSK") == 0){
			sprintf(param.pathK, "");
		
			//fgets(skip, 3, paramFile);
			fscanf (paramFile, "%d", &param.doStoreK);
			fgets(skip2, 3, paramFile);
		}
		else{
		//read doStoreK
		fgets(skip, 12, paramFile);
		fscanf (paramFile, "%d", &param.doStoreK);
		fgets(skip2, 3, paramFile);
		}
		//read nbins
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%d", &param.nbins);
		fgets(skip2, 3, paramFile);
		//read binsfile
		fgets(skip, 11, paramFile);
		fscanf (paramFile, "%s", &param.bins);
		fgets(skip2, 3, paramFile);
		//read putputEdges
		fgets(skip, 18, paramFile);
		fscanf (paramFile, "%s", &param.edges);
		fgets(skip2, 3, paramFile);
		//read kmin
		fgets(skip, 7, paramFile);
		fscanf (paramFile, "%lf", &param.kmin);
		fgets(skip2, 3, paramFile);
		//read qalphaL
		fgets(skip, 10, paramFile);
		fscanf (paramFile, "%lf", &param.qalphaL);
		fgets(skip2, 3, paramFile);
		//read doMean
		fgets(skip, 9, paramFile);
		fscanf (paramFile, "%d", &param.doMean);
		fgets(skip2, 3, paramFile);
		//read Units
		fgets(skip, 8, paramFile);
		fscanf (paramFile, "%d", &param.units);
		fgets(skip2, 3, paramFile);
		//read ReplaceFiles
		fgets(skip, 15, paramFile);
		fscanf (paramFile, "%d", &param.replaceFiles);
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
		else if(strcmp(argv[i], "-HITEMP") == 0){
			param.useHITEMP = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-M") == 0){
			param.nMolecule = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-path") == 0){
			sprintf(param.path, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-pathK") == 0){
			sprintf(param.pathK, "%s", argv[i + 1]);
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
		else if(strcmp(argv[i], "-Mean") == 0){
			param.doMean = atoi(argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}

	}

	if(strcmp(param.bins, "-") != 0){
		printf("Use bins in file %s\n", param.bins);
		param.useIndividualBins = 1;
		FILE *binsfile;
		binsfile = fopen(param.bins, "r");
		if(binsfile == NULL){
			printf("Error: bins file not found: %s\n", param.bins);
			return 0;
		}

		double b;
		int er;
		for(int i = 0; i < 1000000; ++i){
			er = fscanf(binsfile, "%lf", &b);
			if(er <= 0){
				param.nbins = i - 1;
				break;
			}
		}		
		fclose(binsfile);	
	}
	if(strcmp(param.edges, "-") != 0){
		printf("Use output edges in file %s\n", param.edges);
		param.useOutputEdges = 1;
		FILE *edgesfile;
		edgesfile = fopen(param.edges, "r");
		if(edgesfile == NULL){
			printf("Error: output edges file not found: %s\n", param.edges);
			return 0;
		}

		double b;
		int er;
		for(int i = 0; i < 1000000; ++i){
			er = fscanf(edgesfile, "%lf", &b);
			if(er <= 0){
				param.nedges = i;
				break;
			}
		}		
		fclose(edgesfile);	
	}
	if(strcmp(param.PFilename, "-") != 0){
		printf("Use Pressure in file %s\n", param.PFilename);
		param.usePFile = 1;
		FILE *Pfile;
		Pfile = fopen(param.PFilename, "r");
		if(Pfile == NULL){
			printf("Error: Pressure file not found: %s\n", param.PFilename);
			return 0;
		}

		double b;
		int er;
		for(int i = 0; i < 1000000; ++i){
			er = fscanf(Pfile, "%lf", &b);
			if(er <= 0){
				param.nP = i;
				break;
			}
		}		
		fclose(Pfile);	
	}
	if(strcmp(param.ciaSystem, "-") != 0){
		param.useCia = 1;
	}

	return 1;
}

__host__ int readBinFile(Param &param, double *binBoundaries_h){
	FILE *binsfile;
	binsfile = fopen(param.bins, "r");
	int er;
	double binsOld = -100000.0;
	for(int i = 0; i < param.nbins + 1; ++i){
		er = fscanf(binsfile, "%lf", &binBoundaries_h[i]);
		if(er <= 0) return 0;
		if(binBoundaries_h[i] <= binsOld){
			printf("Error; bin boundaries not monotonic growing\n");
			return 0;
		}
		binsOld = binBoundaries_h[i];

		//printf("%g\n", binBoundaries_h[i]);
	}
	fclose(binsfile);	
	return 1;
}

__host__ int readEdgesFile(Param &param, double *outputEdges_h){
	FILE *edgesfile;
	edgesfile = fopen(param.edges, "r");
	int er;
	double edgesOld = -100000.0;
	for(int i = 0; i < param.nedges; ++i){
		er = fscanf(edgesfile, "%lf", &outputEdges_h[i]);
		if(er <= 0) return 0;
		if(outputEdges_h[i] <= edgesOld){
			printf("Error; output edges not monotonic growing\n");
			return 0;
		}
		edgesOld = outputEdges_h[i];

		//printf("%g\n", outputEdges_h[i]);
	}
	fclose(edgesfile);	
	return 1;
}

__host__ int readPFile(Param &param, double *P_h){
	FILE *Pfile;
	Pfile = fopen(param.PFilename, "r");
	int er;
	for(int i = 0; i < param.nP; ++i){
		er = fscanf(Pfile, "%lf", &P_h[i]);
		if(er <= 0) return 0;
	}
	fclose(Pfile);	
	return 1;
}

// ******************************************************************
//This Function reads the Hitran or Hitemp data files
//Author Simon Grimm
//January 2015
// *******************************************************************
__host__ int readFile(Molecule &m, Partition &part, Line &L, double qalphaL, int fi){
	FILE *dataFile;
	dataFile  = fopen(m.dataFilename[fi], "r");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", m.dataFilename[fi]);
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

//int count[40];
//for(int cc = 0; cc < 40; ++cc){
//count[cc] = 0;	
//}
	for(int i = 0; i < m.NL[fi]; ++i){
		fgets(skip, 1, dataFile);
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
	
//		fgets(skip, 36, dataFile);
	
		L.nu_h[i] = strtod(c3, NULL);		
		L.S_h[i] = strtod(c4, NULL);		
		L.A_h[i] = strtod(c5, NULL);		
		L.delta_h[i] = strtod(c10, NULL);
		L.EL_h[i] = strtod(c8, NULL);		
		double gammaAir = strtod(c6, NULL);
		double gammaSelf = strtod(c7, NULL);
		L.alphaL_h[i] = (1.0 - qalphaL) * gammaAir + qalphaL * gammaSelf;
		L.n_h[i] = strtod(c9, NULL);
		L.alphaD_h[i] = 0.0;

//printf("%g %g %g %g %g %g %g %g\n", L.nu_h[i], L.S_h[i], L.A_h[i], gammaAir, gammaSelf, L.EL_h[i], L.n_h[i], L.delta_h[i]);		
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
		double Q = 0.0;
		int Qcheck = 0;
		//Assign the Partition function
		for(int j = 0; j < part.n; ++j){
			if(idAFGL == part.id[j]){
				Q = part.Q[j];
				Qcheck = 1;
			}
		}
		if(Qcheck == 0){
			printf("Error: partition function for AFGL %d not found\n", idAFGL);
			return 0;

		}
		L.Q_h[i] /= exp(Q);
		
//printf("%d %d %d %g %g %d %g %d\n", i, id, idAFGL, exp(Q), L.Q_h[i], part.n, Q, Qcheck);
	}
//for(int cc = 0; cc < 40; ++cc){
//printf("%d %d\n", cc, count[cc]);	
//}
	fclose(dataFile);
	return 1;
}
// ******************************************************************
//This Function reads the prepared Exomol files (user prepareExomol.cpp)
//Author Simon Grimm
//August 2016
// *******************************************************************
__host__ int readFileExomol(Molecule &m, Partition &part, Line &L, int fi, double T){
	FILE *dataFile;
	dataFile  = fopen(m.dataFilename[fi], "rb");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", m.dataFilename[fi]);
		return 0;
	}
	//read line list file		

	double mass = m.ISO[0].m / def_NA;
	double Q = part.Q[0];
	double c = def_h * def_c / (def_kB * T);
	double A, EL;
	for(int i = 0; i < m.NL[fi]; ++i){
	
		fread(&L.nu_h[i], sizeof(double), 1, dataFile);		
		fread(&L.S_h[i], sizeof(double), 1, dataFile);		
		fread(&EL, sizeof(double), 1, dataFile);		
		fread(&A, sizeof(double), 1, dataFile);		


		L.S_h[i] *= exp(-c * EL) * (1.0 - exp(-c * L.nu_h[i])) / Q;
		L.alphaD_h[i] = def_c / L.nu_h[i] * sqrt( mass / (2.0 * def_kB * T));      //inverse Doppler halfwdith
		L.alphaL_h[i] = A / (4.0 * M_PI * def_c);
		L.ID_h[i] = i % maxlines;

		if(L.nu_h[i] == 0.0){
			L.S_h[i] = 0.0;
			L.alphaD_h[i] = 0.0;
			L.alphaL_h[i] = 0.0;
		}

//if(i < 10000) printf("%d %g %g %g\n", i, L.nu_h[i], L.S_h[i], L.alphaD_h[i]);		
		
	}
	fclose(dataFile);
	return 1;
}
// ******************************************************************
//This Function computes the Lorentz halfwidths from the ExoMol default values
//Author Simon Grimm
//August 2016
// *******************************************************************
__host__ int alphaLExomol(Molecule &m, Line &L, int fi, double T, double P){
	for(int i = 0; i < m.NL[fi]; ++i){
//read this numbers for ExoMol define Files
		L.alphaL_h[i] += 0.0700 * pow(296.0 / T, 0.5) * (P / 0.986923); 
//if(i < 10000) printf("%d %g %g %g %g\n", i, L.nu_h[i], L.S_h[i], L.alphaD_h[i], L.alphaL_h[i]);
	}
	return 1;
}


__host__ void readCiaFile(Param param, ciaSystem cia, double *x_h, double *K_h, int Nx, double T, double P, double meanMass){

	FILE *ciaFile;
	ciaFile = fopen(cia.dataFilename, "r");

	char skip[160];
	char c1[12];
	char c2[12];
	char c3[9];
	char c4[9];
	char c5[12];
	double numin;
	double numax;
	int Ncia;
	double Tcia0, Tcia1;
	double max;

	double nu0, nu1;
	double cia0, cia1;

	Tcia1 = 0.0;

	for(int i = 0; i < Nx; ++i){
		K_h[i] = 0.0;
	}

	for(int it = 0; it < cia.Nsets; ++it){
		
		Tcia0 = Tcia1;

		fgets(skip, 21, ciaFile);
		fgets(c1, 11, ciaFile);
		fgets(c2, 11, ciaFile);
		fgets(c3, 8, ciaFile);
		fgets(c4, 8, ciaFile);
		fgets(c5, 11, ciaFile);
		fgets(skip, 37, ciaFile);
	 
		numin = strtod(c1, NULL);
		numax = strtod(c2, NULL);
		Ncia = atoi(c3);
		Tcia1 = strtod(c4, NULL);
		max = strtod(c5, NULL);

//printf("%g %g %d %g %g %g\n", numin, numax, Ncia, Tcia1, Tcia0, max);

		fscanf(ciaFile, "%lf", &nu0);
		fscanf(ciaFile, "%lf", &cia0);
		fscanf(ciaFile, "%lf", &nu1);
		fscanf(ciaFile, "%lf", &cia1);

		int nc = 2;
		
		for(int i = 0; i < Nx; ++i){

			while(x_h[i] >= nu1){
				nu0 = nu1;
				cia0 = cia1;
				fscanf(ciaFile, "%lf", &nu1);
				fscanf(ciaFile, "%lf", &cia1);
//printf("C %g %g %d %d\n", nu1, cia1, i, nc);
				++nc;
				if(nc == Ncia) break;
			}

			if(nu0 <= x_h[i] && x_h[i] < nu1){
				if(Tcia1 <= T){
					K_h[i] = cia0 + (cia1 - cia0) * (x_h[i] - nu0) / (nu1 - nu0);
//if(Tcia1 == 400) printf("A %g %g %g %g %g %g %g\n", nu0, x_h[i], nu1, cia0, cia1, K_h[i], Tcia1);
				}
				if(T < Tcia1 && Tcia0 > 0.0){
					double K0 = K_h[i];
					double K1 = cia0 + (cia1 - cia0) * (x_h[i] - nu0) / (nu1 - nu0);
					K_h[i] = K0 + (K1 - K0) * (T - Tcia0) / (Tcia1 - Tcia0);
//printf("B %g %g %g %g %g %g %g %g %g\n", nu0, x_h[i], nu1, cia0, cia1, K_h[i], T, K0, K1);
				}
			}
		
			if(nu1 == numax) break;

			if(i == Nx - 1){
				while(nc < Ncia){
				fscanf(ciaFile, "%lf", &nu1);
				fscanf(ciaFile, "%lf", &cia1);
//printf("D %g %g %d %d\n", nu1, cia1, i, nc);
				++nc;
				}
			}


		}
		fgets(skip, 3, ciaFile);
		if(T <= Tcia1) break;
	}	

	//See HITRAN cia Paper, Richard et al 2012
	//P0 is set to 1 atm
	//T0 is set to 273.15 K 
	double rho1 = P * 273.15 / T * def_amagat; // numerical density in molecules cm^-3


	for(int i = 0; i < Nx; ++i){
		if(Tcia1 < T){
			K_h[i] = param.kmin;
		}	
		else{
			K_h[i] *= rho1; // K in cm^2 / molecule
			K_h[i] *=def_NA / cia.mass1; //K in cm^2 / g  //should be masss of he in h2he
			K_h[i] += param.kmin;
		}
	}

	fclose(ciaFile);
	for(int i = 0; i < Nx; ++i){
//printf("%g %g %g\n", x_h[i], K_h[i], T);

	}
}


__host__ void Alloc_Line(Line &L, Molecule &m){
	L.nu_h = (double*)malloc(m.NLmax * sizeof(double));
	L.S_h = (double*)malloc(m.NLmax * sizeof(double));
	L.A_h = (double*)malloc(m.NLmax * sizeof(double));
	L.delta_h = (double*)malloc(m.NLmax * sizeof(double));
	L.EL_h = (double*)malloc(m.NLmax * sizeof(double));
	L.alphaL_h = (double*)malloc(m.NLmax * sizeof(double));
	L.alphaD_h = (double*)malloc(m.NLmax * sizeof(double));
	L.n_h = (double*)malloc(m.NLmax * sizeof(double));
	L.mass_h = (double*)malloc(m.NLmax * sizeof(double));
	L.Q_h = (double*)malloc(m.NLmax * sizeof(double));
	L.ID_h = (int*)malloc(m.NLmax * sizeof(int));

	int n = min(maxlines, m.NLmax);

	cudaMalloc((void **) &L.nu_d, n * sizeof(double));
	cudaMalloc((void **) &L.S_d, n * sizeof(double));
	cudaMalloc((void **) &L.A_d, n * sizeof(double));
	cudaMalloc((void **) &L.delta_d, n * sizeof(double));
	cudaMalloc((void **) &L.EL_d, n * sizeof(double));
	cudaMalloc((void **) &L.alphaL_d, n * sizeof(double));
	cudaMalloc((void **) &L.alphaD_d, n * sizeof(double));
	cudaMalloc((void **) &L.n_d, n * sizeof(double));
	cudaMalloc((void **) &L.mass_d, n * sizeof(double));
	cudaMalloc((void **) &L.Q_d, n * sizeof(double));
	cudaMalloc((void **) &L.ID_d, n * sizeof(int));
}
__host__ void Alloc2_Line(Line &L, Molecule &m){
	L.nu_h = (double*)malloc(m.NLmax * sizeof(double));
	L.S_h = (double*)malloc(m.NLmax * sizeof(double));
	L.alphaL_h = (double*)malloc(m.NLmax * sizeof(double));
	L.alphaD_h = (double*)malloc(m.NLmax * sizeof(double));
	L.Q_h = (double*)malloc(m.NLmax * sizeof(double));
	L.ID_h = (int*)malloc(m.NLmax * sizeof(int));

	int n = min(maxlines, m.NLmax);

	cudaMalloc((void **) &L.nu_d, n * sizeof(double));
	cudaMalloc((void **) &L.S_d, n * sizeof(double));
	cudaMalloc((void **) &L.alphaL_d, n * sizeof(double));
	cudaMalloc((void **) &L.alphaD_d, n * sizeof(double));
	cudaMalloc((void **) &L.Q_d, n * sizeof(double));
	cudaMalloc((void **) &L.ID_d, n * sizeof(int));
}

__host__ void Copy_Line(Line &L, Molecule &m, int i){

	int n = min(maxlines, m.NL[i]);

	cudaMemcpy(L.nu_d, L.nu_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.S_d, L.S_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.A_d, L.A_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.delta_d, L.delta_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.EL_d, L.EL_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.alphaL_d, L.alphaL_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.alphaD_d, L.alphaD_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.n_d, L.n_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.mass_d, L.mass_h, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.Q_d, L.Q_h, n * sizeof(double), cudaMemcpyHostToDevice);
}
__host__ void Copy2_Line(Line &L, Molecule &m, int i, int iL){

	int n = min(maxlines, m.NL[i]);

	cudaMemcpy(L.nu_d, L.nu_h + iL, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.S_d, L.S_h + iL, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.alphaL_d, L.alphaL_h + iL, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.alphaD_d, L.alphaD_h + iL, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.Q_d, L.Q_h + iL, n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.ID_d, L.ID_h + iL, n * sizeof(int), cudaMemcpyHostToDevice);
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
__host__ void free2_Line(Line &L){
	free(L.nu_h);
	free(L.S_h);
	free(L.alphaL_h);
	free(L.alphaD_h);
	free(L.Q_h);
	free(L.ID_h);

	cudaFree(L.nu_d);
	cudaFree(L.S_d);
	cudaFree(L.alphaL_d);
	cudaFree(L.alphaD_d);
	cudaFree(L.Q_d);
	cudaFree(L.ID_d);
}

