#include <algorithm> //max
using namespace std;
//**************************************************
// This function initializes the Properties for the Isotopologues
// The values are taken from HITRAN molparam.txt
// The ids must be adapted in order to be consistent with the line list files in HITRAN
// The parameters are:
// m.id: Id of molecule
// m.NL: Number of Lines in the Molecule file
// m.nISO: Number of Isotopologues in the Molecule file
// The last line specifies the name of the Molecule file.

//Author: Simon Grimm
//November 2014
// ****************************************************
int Init(Molecule &m, Param &param, char (*qFilename)[160]){
	FILE *pFile;
	char pFileName[800];

	/*
	//Atoms and ions from Kurucz
	if(param.useHITEMP == 30){
		sprintf(pFileName, "%sgfnew%04d.param", param.path, m.id);
	}

	//Atoms and ions from NIST
	if(param.useHITEMP == 31){
		sprintf(pFileName, "%sNIST%04d.param", param.path, m.id);
	}
	*/

 	sprintf(pFileName, "%s%s.param", param.path, param.mParamFilename);



	pFile = fopen(pFileName, "r");
	if(pFile == NULL){
		printf("Error, no molecule.param file %s\n", pFileName);
		return 0;
	}
	char sp[160];
	fgets(sp, 11, pFile);
	if(strcmp(sp, "Database =") != 0){
		printf("Error in molecule.param file, Database\n");
		return 0;
	}
	fscanf (pFile, "%d", &param.dataBase);
	fgets(sp, 4, pFile);

	fgets(sp, 18, pFile);
	if(strcmp(sp, "Molecule number =") != 0){
		printf("Error in molecule.param file, Molecule number\n");
		return 0;
	}
	fscanf (pFile, "%d", &m.id);
	fgets(sp, 4, pFile);

	printf("Read Molecule param file: |%s| %d %d\n", param.mParamFilename, param.dataBase, m.id);

	fgets(sp, 7, pFile);
	if(strcmp(sp, "Name =") != 0){
		printf("Error in molecule.param file, Name\n");
		return 0;
	}
	fscanf (pFile, "%s", m.mName);
	fgets(sp, 4, pFile);
	//can be "Number of Isotopes =" or "Number of Isotopes ="
	//scan until "="
	fscanf(pFile, "%[^=]s", sp);
	if(strcmp(sp, "Number of Isotopes ") != 0 && strcmp(sp, "Number of Isotopologues ") != 0){
		printf("Error in molecule.param file, Number of Isotopes, |%s|\n", sp);
		return 0;
	}
	//scan "="
	fgets(sp, 2, pFile);

	fscanf (pFile, "%d", &m.nISO);
	//printf("nISO %d\n", m.nISO);
	fgets(sp, 4, pFile);
	m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));

	fgets(sp, 100, pFile);

	for(int i = 0; i < m.nISO; ++i){
		fscanf (pFile, "%s", sp);
		sprintf(m.ISO[i].cid, "%3s", sp);
		fscanf (pFile, "%lf", &m.ISO[i].Ab);
		fscanf (pFile, "%lf", &m.ISO[i].Q);
		fscanf (pFile, "%d",  &m.ISO[i].g);
		fscanf (pFile, "%lf", &m.ISO[i].m);
		fscanf (pFile, "%s", sp);
		sprintf(qFilename[i], "%s%s", param.path, sp);
//printf("%s %g %s\n", m.ISO[i].cid, m.ISO[i].Ab, qFilename[i]);
	}
	fgets(sp, 4, pFile);
	fgets(sp, 38, pFile);
	if(strcmp(sp, "Number of columns in partition File =") != 0){
		printf("Error in molecule.param file, Number of columns in partition File\n");
		return 0;
	}
	fscanf (pFile, "%d", &m.npfcol);
	fgets(sp, 4, pFile);

	fgets(sp, 34, pFile);
	if(strcmp(sp, "Number of line/transition files =") != 0){
		printf("Error in molecule.param file, Number of line/transition files\n");
		return 0;
	}
	fscanf (pFile, "%d", &m.nFiles);
	fgets(sp, 4, pFile);

	fgets(sp, 100, pFile);

	m.NLmax = 0;
	for(int i = 0; i < m.nFiles; ++i){
		fscanf (pFile, "%lld", &m.NL[i]);
		m.NLmax = max(m.NLmax, m.NL[i]);
	}
	fgets(sp, 4, pFile);

	fgets(sp, 19, pFile);
	for(int i = 0; i < m.nFiles + 1; ++i){
		fscanf (pFile, "%d", &m.fileLimit[i]);
//printf("%d\n", m.fileLimit[i]);
	}
	for(int i = 0; i < m.nFiles; ++i){
		if(param.dataBase == 0){
			//Hitran or HITEMP
			if(m.nFiles > 1){
				sprintf(m.dataFilename[i], "%s%s_%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			else{
				sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
			}
		}
		if(param.dataBase == 2){
			//ExoMol
			if(m.nFiles > 1){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			else{
				sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
			}
		}
		if(param.dataBase == 3){
			//IAO (CDSD)
			if(m.nFiles > 1){
				sprintf(m.dataFilename[i], "%s%s__%05d_%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			else{
				sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
			}
		}
		if(param.dataBase == 30){
			//Kurucz
			sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
		}
		if(param.dataBase == 31){
			//NIST
			sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
		}
		if(param.dataBase == 32){
			//VALD
			sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
		}
		
	}
	fgets(sp, 4, pFile);
	fgets(sp, 100, pFile);
	fgets(sp, 19, pFile);
	if(strcmp(sp, "Number of states =") != 0){
		printf("Error in species.param file, Number of states\n");
		return 0;
	}
	fscanf (pFile, "%d", &m.nStates);
	fgets(sp, 4, pFile);

	fgets(sp, 40, pFile);
	if(strcmp(sp, "Number of columns in transition files =") != 0){
		printf("Error in species.param file, Number of columns in transition files\n");
		return 0;
	}
	fscanf (pFile, "%d", &m.ntcol);
	fgets(sp, 4, pFile);

	fgets(sp, 55, pFile);
	if(strcmp(sp, "Default value of Lorentzian half-width for all lines =") != 0){
		printf("Error in species.param file, Default value of Lorentzian half-width for all lines\n");
		return 0;
	}
	fscanf (pFile, "%lf", &m.defaultL);
	fgets(sp, 4, pFile);

	fgets(sp, 54, pFile);
	if(strcmp(sp, "Default value of temperature exponent for all lines =") != 0){
		printf("Error in species.param file, Default value of temperature exponent for all lines\n");
		return 0;
	}
	fscanf (pFile, "%lf", &m.defaultn);
	fgets(sp, 4, pFile);

	fgets(sp, 10, pFile);
	if(strcmp(sp, "Version =") != 0){
		printf("Error in species.param file, Version\n");
		return 0;
	}
	fscanf (pFile, "%d", &m.version);
	fgets(sp, 4, pFile);

	fclose(pFile);

	return 1;
}


int InitCia(Molecule &m, ciaSystem &cia, Param param){
	m.NLmax = 0;
	cia.Nsets = 0;
	cia.mass1 = 1.0;
	if(strcmp(param.ciaSystem, "H2-H2") == 0){
		cia.Nsets = 113;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H2_2011.cia");
		cia.mass1 = 2.0 * 1.00794; //mass of H2 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-H2_eq") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H2_eq_2011.cia");
		cia.mass1 = 2.0 * 1.00794; //mass of H2 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-H2_norm") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H2_norm_2011.cia");
		cia.mass1 = 2.0 * 1.00794; //mass of H2 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-He") == 0){
		cia.Nsets = 339;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-He_2011.cia");
		cia.mass1 = 4.002602; //mass of He in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-He_eq") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-He_eq_2011.cia");
		cia.mass1 = 4.002602; //mass of He in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-He_norm") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-He_norm_2011.cia");
		cia.mass1 = 4.002602; //mass of He in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-CH4_eq") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-CH4_eq_2011.cia");
		cia.mass1 = 16.04246; //mass of CH4 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-CH4_norm") == 0){
		cia.Nsets = 10;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-CH4_norm_2011.cia");
		cia.mass1 = 16.04246; //mass of CH4 in g / mol
	}
	else if(strcmp(param.ciaSystem, "H2-H") == 0){
		cia.Nsets = 4;
		sprintf(cia.dataFilename, "%s%s", param.path, "H2-H_2011.cia");
		cia.mass1 = 1.00794; //mass of H in g / mol
	}
	else{
		printf("Error: cia System not found %s\n", param.ciaSystem);
		return 0;
	}
	return 1;
}

