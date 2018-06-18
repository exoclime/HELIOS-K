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
__host__ int ChebCoeff(Param &param, char *qFilename, Partition &part, double T){

	//Calculate Chebychev polynomial
	double Cheb[NCheb];
	Chebyshev(T, Cheb, NCheb);

	//Read Chebychev Coefficients from q file	
	FILE *qFile;
	//Check size of q.dat file
	qFile = fopen(qFilename, "r");
	if(qFile == NULL){
		printf("Error: q.dat file not found. Path:%s\n", param.path);
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
//March 2018
// *****************************************
__host__ int readPartition(Param &param, int nMolecule, char (*qFilename)[160], Partition &part, double T,  Molecule &m){
	part.n = 1;
	
	part.id = (int*)malloc(sizeof(int));
	part.Q = (double*)malloc( m.nISO * sizeof(double));


	part.id = 0;

	for(int i = 0; i < m.nISO; ++i){

		//Read Chebychev Coefficients from q file	
		FILE *qFile;
		qFile = fopen(qFilename[i], "r");
		printf("Read partition function: %s\n", qFilename[i]);
		if(qFile == NULL){
			printf("Error: partition file not found %s. Path: %s\n", qFilename[i], param.path);
			return 0;
		}
		double T0, T1;
		double q0, q1, q;
		double skip;
		T1 = 0;
		q1 = 1.0;
		q = 1.0;
		int er = 0;
		for(int j = 0; j < 100000; ++j){
			T0 = T1;
			q0 = q1;
			if(m.npfcol == 3){
				er = fscanf (qFile, "%lf", &T1);
				er = fscanf (qFile, "%lf", &q1);
				er = fscanf (qFile, "%lf", &skip);
			}
			else if(m.npfcol == 2){
				er = fscanf (qFile, "%lf", &T1);
				er = fscanf (qFile, "%lf", &q1);
			}
			else{
				printf("Error: partition file not specified\n");
				return 0;
			}
			if(T0 == T1 && T0 < T && j > 0){
				printf("Error: partition function not valid for given temperature %g %g\n", T0, T);
				return 0;
			}
			if (er <= 0) break;
			if(T0 < T && T1 >= T){
				double tt = (T - T0) / (T1 - T0);
				q = (q1 - q0) * tt + q0;
				if(j == 0) q = q1;
				break;
			}
			if(j == 100000 - 1){
				printf("Error: partition function not complete\n");
				return 0;
			}
		}
		fclose(qFile);
//printf("Q %d %g\n", i, q);
		part.Q[i] = q;
	}
	return 1;
}

__host__ int read_parameters(Param &param, char *paramFilename, int argc, char*argv[]){
	//Read parameters from param.dat file
	FILE *paramFile;
	paramFile = fopen(paramFilename, "r");

	param.profile = PROFILE;

	char sp[160];

	for(int j = 0; j < 70; ++j){ //loop around all lines in the param.dat file
		int c;
		for(int i = 0; i < 50; ++i){
			c = fgetc(paramFile);
			if(c == EOF) break;
			sp[i] = char(c);
			if(c == '=' || c == ':'){
				sp[i + 1] = '\0';
				break;
			}
		}
		if(c == EOF) break;
		//read name
		if(strcmp(sp, "name =") == 0){
			fscanf (paramFile, "%s", param.name);
			fgets(sp, 3, paramFile);
		}
		//read T
		else if(strcmp(sp, "T =") == 0){
			fscanf (paramFile, "%lf", &param.T);
			fgets(sp, 3, paramFile);
		}
		//read P
		else if(strcmp(sp, "P =") == 0){
			fscanf (paramFile, "%lf", &param.P);
			fgets(sp, 3, paramFile);
		}
		//read PFile
		else if(strcmp(sp, "PFile =") == 0){
			fscanf (paramFile, "%s", param.PFilename);
			fgets(sp, 3, paramFile);
		}
		//read HITEMP
		else if(strcmp(sp, "useHITEMP =") == 0){
			fscanf (paramFile, "%d", &param.useHITEMP);
			fgets(sp, 3, paramFile);
		}
		//read Molecule
		else if(strcmp(sp, "Molecule =") == 0){
			fscanf (paramFile, "%d", &param.nMolecule);
			fgets(sp, 3, paramFile);
		}
		//read ciaSystem
		else if(strcmp(sp, "ciaSystem =") == 0){
			fscanf (paramFile, "%s", param.ciaSystem);
			fgets(sp, 3, paramFile);
		}
		//read path
		else if(strcmp(sp, "pathToData =") == 0){
			fscanf (paramFile, "%s", param.path);
			fgets(sp, 3, paramFile);
	
			if(strcmp(param.path, "numin") == 0){
				param.path[0] = 0;
		
				//fgets(skip, 3, paramFile);
				fscanf (paramFile, "%lf", &param.numin);
				fgets(sp, 3, paramFile);
			}
		}
		//read numin
		else if(strcmp(sp, "numin =") == 0){
			fscanf (paramFile, "%lf", &param.numin);
			fgets(sp, 3, paramFile);
		}
		//read numax
		else if(strcmp(sp, "numax =") == 0){
			fscanf (paramFile, "%lf", &param.numax);
			fgets(sp, 3, paramFile);
		}
		//read dnu
		else if(strcmp(sp, "dnu =") == 0){
			fscanf (paramFile, "%lf", &param.dnu);
			fgets(sp, 3, paramFile);
		}
		//read Nxb
		else if(strcmp(sp, "Nnu per bin =") == 0){
			fscanf (paramFile, "%d", &param.Nxb);
			fgets(sp, 3, paramFile);
		}
		//read cutMode
		else if(strcmp(sp, "cutMode =") == 0){
			fscanf (paramFile, "%d", &param.cutMode);
			fgets(sp, 3, paramFile);
		}
		//read cut
		else if(strcmp(sp, "cut =") == 0){
			fscanf (paramFile, "%lf", &param.cut);
			fgets(sp, 3, paramFile);
		}
		//read doResampling
		else if(strcmp(sp, "doResampling =") == 0){
			fscanf (paramFile, "%d", &param.doResampling);
			fgets(sp, 3, paramFile);
		}
		//read nC
		else if(strcmp(sp, "nC =") == 0){
			fscanf (paramFile, "%d", &param.nC);
			fgets(sp, 3, paramFile);
			if(param.nC > NmaxSample){
				printf("nC larger than NmaxSample, reduced to %d\n", NmaxSample);
				param.nC = NmaxSample;
			}
		}
		//read doTransmission
		else if(strcmp(sp, "doTransmission =") == 0){
			fscanf (paramFile, "%d", &param.doTransmission);
			fgets(sp, 3, paramFile);
		}
		//read nTr
		else if(strcmp(sp, "nTr =") == 0){
			fscanf (paramFile, "%d", &param.nTr);
			fgets(sp, 3, paramFile);
		}
		//read dTr
		else if(strcmp(sp, "dTr =") == 0){
			fscanf (paramFile, "%lf", &param.dTr);
			fgets(sp, 3, paramFile);
		}
		//read doStoreFullK
		else if(strcmp(sp, "doStoreFullK =") == 0){
			fscanf (paramFile, "%d", &param.doStoreFullK);
			fgets(sp, 3, paramFile);
		}
		//read pathK
		else if(strcmp(sp, "pathToK =") == 0){
			fscanf (paramFile, "%s", param.pathK);
			fgets(sp, 3, paramFile);
		
			if(strcmp(param.pathK, "doStoreSK") == 0){
				param.pathK[0] = 0;
		
				//fgets(skip, 3, paramFile);
				fscanf (paramFile, "%d", &param.doStoreK);
				fgets(sp, 3, paramFile);
			}
			else{
				//read doStoreK
				if(strcmp(sp, "doStoreSK =") == 0){
					fscanf (paramFile, "%d", &param.doStoreK);
					fgets(sp, 3, paramFile);
				}
			}
		}
		//read nbins
		else if(strcmp(sp, "nbins =") == 0){
			fscanf (paramFile, "%d", &param.nbins);
			fgets(sp, 3, paramFile);
		}
		//read binsfile
		else if(strcmp(sp, "binsFile =") == 0){
			fscanf (paramFile, "%s", param.bins);
			fgets(sp, 3, paramFile);
		}
		//read outputEdges
		else if(strcmp(sp, "OutputEdgesFile =") == 0){
			fscanf (paramFile, "%s", param.edges);
			fgets(sp, 3, paramFile);
		}
		//read kmin
		else if(strcmp(sp, "kmin =") == 0){
			fscanf (paramFile, "%lf", &param.kmin);
			fgets(sp, 3, paramFile);
		}
		//read qalphaL
		else if(strcmp(sp, "qalphaL =") == 0){
			fscanf (paramFile, "%lf", &param.qalphaL);
			fgets(sp, 3, paramFile);
		}
		//read doMean
		else if(strcmp(sp, "doMean =") == 0){
			fscanf (paramFile, "%d", &param.doMean);
			fgets(sp, 3, paramFile);
		}
		//read Units
		else if(strcmp(sp, "Units =") == 0){
			fscanf (paramFile, "%d", &param.units);
			fgets(sp, 3, paramFile);
		}
		//read ReplaceFiles
		else if(strcmp(sp, "ReplaceFiles =") == 0){
			fscanf (paramFile, "%d", &param.replaceFiles);
			fgets(sp, 3, paramFile);
		}
		//read RLOW
		else if(strcmp(sp, "RLOW =") == 0){
			fscanf (paramFile, "%d", &param.RLOW);
			fgets(sp, 3, paramFile);
		}
		//read profile
		else if(strcmp(sp, "profile =") == 0){
			fscanf (paramFile, "%d", &param.profile);
			fgets(sp, 3, paramFile);
		}
		else{
			printf("Undefined line in param.dat file: line %d\n", j);
			return 0;
		}
	}

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
__host__ int readFile(Param param, Molecule &m, Partition &part, Line &L, double qalphaL, int fi){
	FILE *dataFile;
	char dataFilename[160];
	sprintf(dataFilename, "%sbin", m.dataFilename[fi]);
	dataFile = fopen(dataFilename, "rb");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", dataFilename);
		return 0;
	}
	//read line list file		
	double gammaAir, gammaSelf;
	double mass;
	int id, idAFGL;
	double S;

	for(int i = 0; i < m.NL[fi]; ++i){
		fread(&id, sizeof(int), 1, dataFile);
		fread(&L.nu_h[i], sizeof(double), 1, dataFile);
		fread(&S, sizeof(double), 1, dataFile);
		fread(&L.EL_h[i], sizeof(double), 1, dataFile);
		fread(&L.A_h[i], sizeof(double), 1, dataFile);
		fread(&L.delta_h[i], sizeof(double), 1, dataFile);
		fread(&gammaAir, sizeof(double), 1, dataFile);
		fread(&gammaSelf, sizeof(double), 1, dataFile);
		fread(&L.n_h[i], sizeof(double), 1, dataFile);

		double Q = 0.0;
		for(int j = 0; j < m.nISO; ++j){
//if(i < 10) printf("%d %d\n", id, m.ISO[j].id);
			if(id == m.ISO[j].id){
				mass = m.ISO[j].m / def_NA;
				idAFGL = m.ISO[j].AFGL;
				if(m.npfcol != 0){
					Q = part.Q[j];
				}
			}
		}

		if(m.npfcol == 0){
			//use q.dat file
			int Qcheck = 0;
			//Assign the Partition function
			for(int j = 0; j < part.n; ++j){
				if(idAFGL == part.id[j]){
					Q = exp(part.Q[j]);
					Qcheck = 1;
				}
			}
			if(Qcheck == 0){
				printf("Error: partition function for AFGL %d not found. %d %d\n", idAFGL, i, id);
				return 0;

			}
		}

		L.vy_h[i] = (1.0 - qalphaL) * gammaAir + qalphaL * gammaSelf;
		S /= Q;
		L.S_h[i] = S;
		L.ialphaD_h[i] = def_c * sqrt( mass / (2.0 * def_kB * param.T));      //inverse Doppler halfwdith, 1.0/nu is missing here and inserted later
                L.ID_h[i] = i % maxlines;
//if(i < 10) printf("%d %g %g %g %g %g %g\n", i, L.nu_h[i], L.S_h[i], L.ialphaD_h[i], L.EL_h[i], 0.0, Q);
		
	}

	fclose(dataFile);
	return 1;
}
// ******************************************************************
//This Function reads the prepared Exomol files (use prepareExomol.cpp)
//Author Simon Grimm
//August 2016
// *******************************************************************
__host__ int readFileExomol(Param param, Molecule &m, Partition &part, Line &L, int fi){
	FILE *dataFile;
	char dataFilename[160];
	sprintf(dataFilename, "%sbin", m.dataFilename[fi]);
	dataFile  = fopen(dataFilename, "rb");

	if(dataFile == NULL){
		printf("Error: line list file not found %s\n", dataFilename);
		return 0;
	}
	//read line list file		

	double mass = m.ISO[0].m / def_NA;
	double Q = part.Q[0];
	double c = def_h * def_c / (def_kB * param.T);
	double A, EL;
	double S;
	double GammaN = 0.0; 	//natural broadening parameter
	for(int i = 0; i < m.NL[fi]; ++i){
	
		fread(&L.nu_h[i], sizeof(double), 1, dataFile);		
		fread(&S, sizeof(double), 1, dataFile);		
		fread(&EL, sizeof(double), 1, dataFile);		
		fread(&A, sizeof(double), 1, dataFile);		
		//include the following for Kurucz
		if(param.useHITEMP == 30){
			fread(&GammaN, sizeof(double), 1, dataFile);		
		}
//if(i < 100) printf("%d %g %g %g %g\n", i, L.nu_h[i], S, EL, A);
		L.ialphaD_h[i] = def_c / L.nu_h[i] * sqrt( mass / (2.0 * def_kB * param.T));      //inverse Doppler halfwdith
		S *= exp(-c * EL) * (1.0 - exp(-c * L.nu_h[i])) / Q * L.ialphaD_h[i];
		L.vy_h[i] = A / (4.0 * M_PI * def_c) + GammaN / (2.0 * M_PI);						//alphaL
		if(param.useIndividualX == 0){
			L.va_h[i] = (float)((param.numin - L.nu_h[i]) * L.ialphaD_h[i]);
			L.vb_h[i] = (float)(param.dnu * L.ialphaD_h[i]);
		}
		else{
			L.va_h[i] = (float)(-L.nu_h[i] * L.ialphaD_h[i]);
			L.vb_h[i] = (float)(L.ialphaD_h[i]);
		}

		L.vcut2_h[i] = (float)(param.cut * param.cut * L.ialphaD_h[i] * L.ialphaD_h[i]); //square of modified cut lenght
		if(param.cutMode == 2){
			L.vcut2_h[i] = (float)(param.cut * param.cut);
		}
		L.ID_h[i] = i % maxlines;
		L.S_h[i] = S;
		L.S1_h[i] = 0.0;
// if(i < 100) printf("%d %g %g %g %g %g %g %g %g %g\n", i, L.nu_h[i], L.S_h[i], L.ialphaD_h[i], EL, exp(-c * L.nu_h[i]), Q, GammaN, L.vy_h[i], exp(-c * EL));

		if(L.nu_h[i] == 0.0){
			L.S1_h[i] = 0.0;
			L.ialphaD_h[i] = 0.0;
			L.vy_h[i] = 0.0;
			L.va_h[i] = 0.0f;
			L.vb_h[i] = 0.0f;
			L.vcut2_h[i] = 0.0f;
		}

//if(i < 10000) printf("%d %g %g %g\n", i, L.nu_h[i], L.S_h[i], L.ialphaD_h[i]);		
		
	}
	fclose(dataFile);
	return 1;
}
// ******************************************************************
//This Function computes the Lorentz halfwidths from the ExoMol default values
//Author Simon Grimm
//August 2016
// *******************************************************************
__host__ int alphaLExomol(Param param, Molecule &m, Line &L, int fi, double T, double P){
	for(int i = 0; i < m.NL[fi]; ++i){
		L.vy_h[i] += (m.defaultL * pow(296.0 / T, m.defaultn) * (P / 0.986923));
		L.vy_h[i] *= L.ialphaD_h[i]; 
		if(param.profile < 4){
			L.S1_h[i] = L.S_h[i] * L.vy_h[i] / M_PI;
		}
		else{
			L.S1_h[i] = L.S_h[i];
		}
		if(param.cutMode == 1){
			L.vcut2_h[i] = (float)(param.cut * param.cut * L.vy_h[i] * L.vy_h[i]);
		}
//if(i < 100) printf("%d %g %g %g %g %g %g %g\n", i, L.nu_h[i], L.S_h[i], (float)(L.S_h[i]), L.ialphaD_h[i], L.vy_h[i], (m.defaultL * pow(296.0 / T, m.defaultn) * (P / 0.986923)));
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
	int Ncia;
	double Tcia0, Tcia1;

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
	 
		Ncia = atoi(c3);
		Tcia1 = strtod(c4, NULL);
		double numax = strtod(c2, NULL);
//double max = strtod(c5, NULL);
//double numin = strtod(c1, NULL);
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
	L.S1_h = (double*)malloc(m.NLmax * sizeof(double));
	L.A_h = (double*)malloc(m.NLmax * sizeof(double));
	L.delta_h = (double*)malloc(m.NLmax * sizeof(double));
	L.EL_h = (double*)malloc(m.NLmax * sizeof(double));
	L.vy_h = (double*)malloc(m.NLmax * sizeof(double));
	L.va_h = (float*)malloc(m.NLmax * sizeof(float));
	L.vb_h = (float*)malloc(m.NLmax * sizeof(float));
	L.vcut2_h = (float*)malloc(m.NLmax * sizeof(float));
	L.ialphaD_h = (double*)malloc(m.NLmax * sizeof(double));
	L.n_h = (double*)malloc(m.NLmax * sizeof(double));
	L.Q_h = (double*)malloc(m.NLmax * sizeof(double));
	L.ID_h = (int*)malloc(m.NLmax * sizeof(int));

	int n = min(maxlines, m.NLmax);

	cudaMalloc((void **) &L.nu_d, n * sizeof(double));
	cudaMalloc((void **) &L.S_d, n * sizeof(double));
	cudaMalloc((void **) &L.Sf_d, n * sizeof(float));
	cudaMalloc((void **) &L.S1_d, n * sizeof(double));
	cudaMalloc((void **) &L.S1f_d, n * sizeof(float));
	cudaMalloc((void **) &L.A_d, n * sizeof(double));
	cudaMalloc((void **) &L.delta_d, n * sizeof(double));
	cudaMalloc((void **) &L.EL_d, n * sizeof(double));
	cudaMalloc((void **) &L.vy_d, n * sizeof(double));
	cudaMalloc((void **) &L.vyf_d, n * sizeof(float));
	cudaMalloc((void **) &L.va_d, n * sizeof(float));
	cudaMalloc((void **) &L.vb_d, n * sizeof(float));
	cudaMalloc((void **) &L.vcut2_d, n * sizeof(float));
	cudaMalloc((void **) &L.ialphaD_d, n * sizeof(double));
	cudaMalloc((void **) &L.n_d, n * sizeof(double));
	cudaMalloc((void **) &L.Q_d, n * sizeof(double));
	cudaMalloc((void **) &L.ID_d, n * sizeof(int));
}
__host__ void Alloc2_Line(Line &L, Molecule &m){
	L.nu_h = (double*)malloc(m.NLmax * sizeof(double));
	L.S_h = (double*)malloc(m.NLmax * sizeof(double));
	L.S1_h = (double*)malloc(m.NLmax * sizeof(double));
	L.vy_h = (double*)malloc(m.NLmax * sizeof(double));
	L.va_h = (float*)malloc(m.NLmax * sizeof(float));
	L.vb_h = (float*)malloc(m.NLmax * sizeof(float));
	L.vcut2_h = (float*)malloc(m.NLmax * sizeof(float));
	L.ialphaD_h = (double*)malloc(m.NLmax * sizeof(double));
	L.Q_h = (double*)malloc(m.NLmax * sizeof(double));
	L.ID_h = (int*)malloc(m.NLmax * sizeof(int));

	int n = min(maxlines, m.NLmax);

	cudaMalloc((void **) &L.nu_d, n * sizeof(double));
	cudaMalloc((void **) &L.S_d, n * sizeof(double));
	cudaMalloc((void **) &L.Sf_d, n * sizeof(float));
	cudaMalloc((void **) &L.S1_d, n * sizeof(double));
	cudaMalloc((void **) &L.S1f_d, n * sizeof(float));
	cudaMalloc((void **) &L.vy_d, n * sizeof(double));
	cudaMalloc((void **) &L.vyf_d, n * sizeof(float));
	cudaMalloc((void **) &L.va_d, n * sizeof(float));
	cudaMalloc((void **) &L.vb_d, n * sizeof(float));
	cudaMalloc((void **) &L.vcut2_d, n * sizeof(float));
	cudaMalloc((void **) &L.ialphaD_d, n * sizeof(double));
	cudaMalloc((void **) &L.Q_d, n * sizeof(double));
	cudaMalloc((void **) &L.ID_d, n * sizeof(int));
}

__host__ void Copy_Line(Line &L, Molecule &m, int iL, int NL){

	cudaMemcpy(L.nu_d, L.nu_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.S_d, L.S_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.S1_d, L.S1_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.vy_d, L.vy_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.va_d, L.va_h + iL, NL * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L.vb_d, L.vb_h + iL, NL * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L.vcut2_d, L.vcut2_h + iL, NL * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L.A_d, L.A_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.delta_d, L.delta_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.EL_d, L.EL_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.ialphaD_d, L.ialphaD_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.n_d, L.n_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.ID_d, L.ID_h + iL, NL * sizeof(int), cudaMemcpyHostToDevice);
}

__host__ void Copy2_Line(Line &L, Molecule &m, int iL, int NL){

	cudaMemcpy(L.nu_d, L.nu_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.S_d, L.S_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.S1_d, L.S1_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.vy_d, L.vy_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.va_d, L.va_h + iL, NL * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L.vb_d, L.vb_h + iL, NL * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L.vcut2_d, L.vcut2_h + iL, NL * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L.ialphaD_d, L.ialphaD_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.Q_d, L.Q_h + iL, NL * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(L.ID_d, L.ID_h + iL, NL * sizeof(int), cudaMemcpyHostToDevice);
}

__host__ void free_Line(Line &L){
	free(L.nu_h);
	free(L.S_h);
	free(L.S1_h);
	free(L.A_h);
	free(L.delta_h);
	free(L.EL_h);
	free(L.vy_h);
	free(L.va_h);
	free(L.vb_h);
	free(L.vcut2_h);
	free(L.ialphaD_h);
	free(L.n_h);
	free(L.Q_h);
	free(L.ID_h);

	cudaFree(L.nu_d);
	cudaFree(L.S_d);
	cudaFree(L.Sf_d);
	cudaFree(L.S1_d);
	cudaFree(L.S1f_d);
	cudaFree(L.A_d);
	cudaFree(L.delta_d);
	cudaFree(L.EL_d);
	cudaFree(L.vy_d);
	cudaFree(L.vyf_d);
	cudaFree(L.va_d);
	cudaFree(L.vb_d);
	cudaFree(L.vcut2_d);
	cudaFree(L.ialphaD_d);
	cudaFree(L.n_d);
	cudaFree(L.Q_d);
	cudaFree(L.ID_d);
}
__host__ void free2_Line(Line &L){
	free(L.nu_h);
	free(L.S_h);
	free(L.S1_h);
	free(L.vy_h);
	free(L.va_h);
	free(L.vb_h);
	free(L.vcut2_h);
	free(L.ialphaD_h);
	free(L.Q_h);
	free(L.ID_h);

	cudaFree(L.nu_d);
	cudaFree(L.S_d);
	cudaFree(L.Sf_d);
	cudaFree(L.S1_d);
	cudaFree(L.S1f_d);
	cudaFree(L.vy_d);
	cudaFree(L.vyf_d);
	cudaFree(L.va_d);
	cudaFree(L.vb_d);
	cudaFree(L.vcut2_d);
	cudaFree(L.ialphaD_d);
	cudaFree(L.Q_d);
	cudaFree(L.ID_d);
}

