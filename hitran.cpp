#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <math.h>

#define def_NA  6.0221412927e23  //Avogadro Constant  1/mol

// ./hitran -M 1 -ISO 1 -in hit16

// ******************************************************************
//This Function reads the Hitran or Hitemp data files
//and creates the <species name>.param and *.bin files
//Author: Simon Grimm
//May 2019
// *******************************************************************

int main(int argc, char *argv[]){

	char iso[4];
	sprintf(iso, "%s", "");
	int M = -1;

	char name[160];
	sprintf(name, "%s", "");
	char inname[160];
	sprintf(inname, "%s", "");

	//Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-M") == 0){
			M = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-ISO") == 0){
			sprintf(iso, "%1s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-in") == 0){
			sprintf(inname, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-out") == 0){		//optional
			sprintf(name, "%s", argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}
	}
	
	//M = 2;
	//sprintf(inname, "cdsd_");
	
	//sprintf(iso, "%1s", "1");
	//sprintf(name, "12C-16O2__CDSD_4000");

	//sprintf(iso, "%1s", "2");
	//sprintf(name, "13C-16O2__CDSD_4000");

	//sprintf(iso, "%1s", "3");
	//sprintf(name, "16O-12C-18O__CDSD_4000");

	//sprintf(iso, "%1s", "4");
	//sprintf(name, "16O-12C-17O__CDSD_4000");


	if(M <= 0){
		printf("Error, no -M agrument given\n");
		return 0;
	}
	if(strcmp(inname, "") == 0){
		printf("Error, no -in agrument given\n");
		return 0;
	}
	if(strcmp(name, "") == 0){
		sprintf(name, "%s", inname);
	}


	// ********************************************************
	//Read isotopologue parameters
	FILE *IsoFile;
	IsoFile = fopen("Hitran_species.dat", "r");
	if(IsoFile == NULL){
		printf("Error, Hitran_species.dat file not found\n");
		return 0;
	}

	char index[20][4];
	char qfile[20][32];
	char c[4];
	double abundance[20];
	double m[20];
	double Q0[20];
	double g[20];
	char formula[20][32];
	char skip[4];


	int nISO = 0;
	for(int i = 0; i < 200; ++i){

		fscanf(IsoFile, "%s %lf %lf %lf %s %lf %s\n", index[nISO], &abundance[nISO], &m[nISO], &Q0[nISO], qfile[nISO], &g[nISO], formula[nISO]);

		//ensure that leading space exists
		char iii[4];
		sprintf(iii, "%3s", index[nISO]);
	
		sprintf(index[nISO], "%3s", iii);
		
		//extract molecule index
		char ii[3];
		sprintf(ii, "%.2s", iii);
		int molecule = atoi(ii);

		char iiso[1];
		sprintf(iiso, "%s", &iii[2]);

		//printf("%d |%s|%s| %g | %d %d |%s|\n", nISO, index[nISO], iii, abundance[nISO], molecule, M, iiso);

		if(iso[0] == 0){
			if(molecule == M){
				printf("ISOm %d %s %g %g %g %s %g\n", molecule, index[nISO], abundance[nISO], m[nISO], Q0[nISO], qfile[nISO], g[nISO]);
				++nISO;
			}
		}
		else{
			if(molecule == M && (strcmp(iiso, iso) == 0)){
				printf("ISO1 %d %s %g %g %g %s %g\n", molecule, index[nISO], abundance[nISO], m[nISO], Q0[nISO], qfile[nISO], g[nISO]);
				abundance[nISO] = 1.0;
				++nISO;
			}
		}
	}

	fclose(IsoFile);
	// **************************************************


	int count = 0;
	double numax = 0.0;
	double nuOld = 0.0;
	int nFiles = 0;


	int filesCount[1000];
	int filesRange[1000];

	//read all files in the current directory
	struct dirent **namelist;
	int n;
	n = scandir(".", &namelist, NULL, alphasort);
	printf("total number of files in directory %d\n", n);

	int nn = 0;
	int range0, range1;
	int len = 0;
	int slen = 0;
	//scan and count desired files 
	for(int j = 0; j < n; ++j){
		std::string str = std::string(namelist[j]->d_name);
		//printf("--- %d %s %s %s\n", j, namelist[j]->d_name, str.c_str(), inname);

		//check if file name contains name and starts with Molecule id
		char s2[32];
		sprintf(s2,"%02d_", M);
		char s3[32];
		sprintf(s3,".par");

		//extract file range from string
		int ppos = int(str.find(inname));	//position of name in string
		std::string sbstr;
		std::string sbstr0;
		std::string sbstr1;
		if(ppos > 3){
			sbstr = str.substr(3, ppos - 4);
			slen = sbstr.length();
			sbstr0 = sbstr.substr(0, slen/2);
			sbstr1 = sbstr.substr(slen/2+1, slen/2);
			range0 = std::atoi(sbstr0.c_str());
			range1 = std::atoi(sbstr1.c_str());
		}
		else{
			sbstr = "";
			slen = 0;
			sbstr0 = "";
			sbstr1 = "";
			range0 = 0;
			range1 = 0;
		}
		len = str.length();

		// find inname in string and find molecule in string and find ".par" in file name
		if(str.find(inname) != std::string::npos && str.rfind(s2, 0) == 0 && str.compare(len - 4, 4,s3) == 0){
			++nn;
			printf("----- %d %s %d %d %d |%s| %d %d \n", nn, str.c_str(), ppos, len, slen, sbstr.c_str(), range0, range1);
		}
	}
	printf("number of data files in directory %d\n", nn);
	
	FILE *binFile;
	char binFileName[160];

	FILE *paramFile;
	char paramFileName[160];

	if(iso[0] == 0){
		sprintf(paramFileName, "%02d_%s.param", M, name);
	}
	else{
		sprintf(paramFileName, "%02d_%s_%s.param", M, iso, name);
	}
	paramFile = fopen(paramFileName, "w");

	fprintf(paramFile, "Database = 0\n");
	fprintf(paramFile, "Molecule number = %d\n", M);
	if(iso[0] == 0){
		fprintf(paramFile, "Name = %02d_%s\n", M, name);
	}
	else{
		fprintf(paramFile, "Name = %02d_%s_%s\n", M, iso, name);

	}
	fprintf(paramFile, "Number of Isotopologues = %d\n", nISO);
	fprintf(paramFile, "#ID Abundance      Q(296K)   g     Molar Mass(g)  partition file :\n");
	for(int i = 0; i < nISO; ++i){
		fprintf(paramFile, "%s %12.10g %8.6g %5.5g %14.14g %s\n", index[i], abundance[i], Q0[i], g[i], m[i], qfile[i]);

	}

	double nu = 0.0;
	for(int j = 0; j < n; ++j){
		std::string str = std::string(namelist[j]->d_name);

		//printf("-- %s %s %s\n", namelist[j]->d_name, str.c_str(), inname);
		char s2[32];
		sprintf(s2,"%02d_", M);
		char s3[32];
		sprintf(s3, ".par");
		
		//extract file range from string
		int ppos = int(str.find(inname));	//position of name in string
		std::string sbstr;
		std::string sbstr0;
		std::string sbstr1;
		if(ppos > 3){
			sbstr = str.substr(3, ppos - 4);
			slen = sbstr.length();
			sbstr0 = sbstr.substr(0, slen/2);
			sbstr1 = sbstr.substr(slen/2+1, slen/2);
			range0 = std::atoi(sbstr0.c_str());
			range1 = std::atoi(sbstr1.c_str());
		}
		else{
			sbstr = "";
			slen = 0;
			sbstr0 = "";
			sbstr1 = "";
			range0 = 0;
			range1 = 0;
		}
		int len = str.length();

		if(str.find(inname) != std::string::npos && str.rfind(s2, 0) == 0 && str.compare(len - 4, 4,s3) == 0){

			printf("---- %d %s\n", n, str.c_str());
			FILE *dataFile;
			dataFile = fopen(str.c_str(), "r");

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

			char skip[161];

			for(int i = 0; i < 100000000; ++i){

				if(fgets(c1, 3, dataFile) == NULL){
					int nnu;
					if(slen == 0){
						nnu = ceil(nu);
					}
					else{
						nnu = range1;
					}
					printf("end %d %g %g %d %d\n", i, nu, nuOld, count, nnu);
					filesCount[nFiles] = count;
					filesRange[nFiles] = nnu;
					++nFiles;
					count = 0;
					break;
				}
				fgets(c2, 2, dataFile);
				fgets(c3, 13, dataFile);
				fgets(c4, 11, dataFile);
				fgets(c5, 11, dataFile);
				fgets(c6, 6, dataFile);
				fgets(c7, 6, dataFile);
				fgets(c8, 11, dataFile);
				fgets(c9, 5, dataFile);
				fgets(c10, 9, dataFile);

				fgets(skip, 160, dataFile);

				int id = atoi(c1);
				char cid[4];
				sprintf(cid, "%2d%s", id, c2);

				nu = strtod(c3, NULL);
				double S = strtod(c4, NULL);
				double A = strtod(c5, NULL);
				double delta = strtod(c10, NULL);
				double EL = strtod(c8, NULL);
				double gammaAir = strtod(c6, NULL);
				double gammaSelf = strtod(c7, NULL);
				double n = strtod(c9, NULL);

				//if(i < 1000) printf("**** %d %g %g %g %g %g %g %g %g | %g |%s|%s|\n", i, nu, S, A, delta, EL, gammaAir, gammaSelf, n, nuOld, c2, iso);


				double mass = 0.0;
				double q0 = -1000.0;

				/*
				//open next binary file
				if(int(nu / 100) != int(nuOld / 100)){
						
					int nnu = int(nu / 100) * 100;
					printf("********** %g %g %d %d %d\n", nu, nuOld, count, int(nuOld / 100) * 100, nnu);
					filesCount[nFiles] = count;
					filesRange[nFiles] = nnu;
					count = 0;
					++nFiles;
					fclose(binFile);
					sprintf(binFileName, "%s__%05d_%05d.bin", name, nnu, nnu + 100);
					binFile = fopen(binFileName, "wb");
				}
		
				*/
				if(strcmp(c2, iso) == 0 || strcmp(iso, "") == 0){
					//printf("**** %d %s %g %g %g %g %g %g %g %g | %g |%s|%s|\n", i, cid, nu, S, A, delta, EL, gammaAir, gammaSelf, n, nuOld, c2, iso);
					//first entry
					if(nFiles == 0){
						int nnu;
						if(slen == 0){
							nnu = int(nu);
						}
						else{
							nnu = range0;
						}
						printf("********** %g %g %d %d\n", nu, nuOld, count, nnu);
						filesCount[nFiles] = count;
						filesRange[nFiles] = nnu;

						++nFiles;
					}
					if(count == 0){
						if(nn == 1){
							if(strcmp(iso, "") == 0){
								sprintf(binFileName, "%02d_%s.bin", M, name);
							}
							else{
								sprintf(binFileName, "%02d_%s_%s.bin", M, iso, name);
							}
						}
						else{
							if(strcmp(iso, "") == 0){
								sprintf(binFileName, "%02d_%s_%s.bin", M, name, sbstr.c_str());
							}
							else{
								sprintf(binFileName, "%02d_%s_%s_%s.bin", M, iso, name, sbstr.c_str());
							}
						}
						printf("******** open bin file %s\n", binFileName);

						binFile = fopen(binFileName, "wb");
					}


					//Assign the Isotopologue properties
					for(int k = 0; k < nISO; ++k){
						//printf("k %d |%s|%s| \n", k, cid, index[k]);
						if(strcmp(cid, index[k]) == 0){
							mass = m[k] / def_NA;
							q0 = Q0[k];
							break;
						}
					}
					if(q0 < -800){
						printf("Error in assigning isotopologue properties\n");
					}

					S /= mass;
					S = S * q0;

					//if(count < 1000 || count % 100000 == 0) printf("%s | %s %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n", iso, cid, nu, S, A, gammaAir, gammaSelf, EL, n, delta, mass, q0);
	
					fwrite(&cid, 4*sizeof(char), 1, binFile);
					fwrite(&nu, sizeof(double), 1, binFile);
					fwrite(&S, sizeof(double), 1, binFile);
					fwrite(&EL, sizeof(double), 1, binFile);
					fwrite(&A, sizeof(double), 1, binFile);
					fwrite(&delta, sizeof(double), 1, binFile);
					fwrite(&gammaAir, sizeof(double), 1, binFile);
					fwrite(&gammaSelf, sizeof(double), 1, binFile);
					fwrite(&n, sizeof(double), 1, binFile);
	
					++count;
					numax = fmax(numax, nu);
				}
				nuOld = nu;

			}
			printf("nFiles %d,  Number of lines %d, numax %g\n", nFiles, count, numax);

			fclose(dataFile);
		}
	}

	for(int i = 0; i < nFiles; ++i){
		printf("F %d %d\n", filesCount[i], filesRange[i]);
	}

	fprintf(paramFile, "Number of columns in partition File = 2\n");
	fprintf(paramFile, "Number of line/transition files = %d\n", nFiles - 1);
	fprintf(paramFile, "Number of lines per file :\n");
	for(int i = 0; i < nFiles - 1; ++i){
        	fprintf(paramFile, "%d\n", filesCount[i + 1]);
	}
	fprintf(paramFile, "Line file limits :\n");
	for(int i = 0; i < nFiles; ++i){
        	fprintf(paramFile, "%d\n", filesRange[i]);

	}	

	fprintf(paramFile, "#ExoMol :\n");
	fprintf(paramFile, "Number of states = 0\n");
	fprintf(paramFile, "Number of columns in transition files = 0\n");
	fprintf(paramFile, "Default value of Lorentzian half-width for all lines = 0\n");
	fprintf(paramFile, "Default value of temperature exponent for all lines = 0\n");
	fprintf(paramFile, "Version = 0\n");

	fclose(paramFile);
	return 0;
}

