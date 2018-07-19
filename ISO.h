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
void Init(Molecule &m, Param &param, char (*qFilename)[160]){
	FILE *pFile;
	char pFileName[160];

	if(m.id == 1){//H2O
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "01_hit16.param");
		}
		if(param.useHITEMP == 1){
	        	sprintf(pFileName, "%s%s", param.path, "01_HITEMP2010.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-16O__BT2.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-16O__POKAZATEL.param");
		}
                if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-18O__HotWat78.param");
		}
                if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-17O__HotWat78.param");
		}
                if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "1H-2H-16O__VTT.param");
		}
	}
	if(m.id == 2){//CO2
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "02_hit16.param");
		}
		if(param.useHITEMP == 1){
	        	sprintf(pFileName, "%s%s", param.path, "02_HITEMP2010.param");
		}
	}
	if(m.id == 3){//O3
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "03_hit16.param");
		}
	}
	if(m.id == 4){//N2O
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "04_hit16.param");
		}
	}
	if(m.id == 5){//CO
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "05_hit16.param");
		}
		if(param.useHITEMP == 1){
	        	sprintf(pFileName, "%s%s", param.path, "05_HITEMP2010.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-16O__Li2015.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczCO.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczCOax.param");
		}
	}
	if(m.id == 6){//CH4
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "06_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-1H4__YT10to10.param");
		}
	}
	if(m.id == 7){//O2
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "07_hit16.param");
		}
	}
	if(m.id == 8){//NO
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "08_hit16.param");
		}
		if(param.useHITEMP == 1){
	        	sprintf(pFileName, "%s%s", param.path, "08_HITEMP2010.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "14N-16O__NOname.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "15N-16O__NOname.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "14N-18O__NOname.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "15N-18O__NOname.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "14N-17O__NOname.param");
		}
		if(param.useHITEMP == 7){
	        	sprintf(pFileName, "%s%s", param.path, "15N-17O__NOname.param");
		}
	}
	if(m.id == 9){//SO2
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "09_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "32S-16O2__ExoAmes.param");
		}
	}
	if(m.id == 11){//NH3
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "11_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "14N-1H3__BYTe.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "15N-1H3__BYTe-15.param");
		}
	}
	if(m.id == 12){//HNO3
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "12_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H-14N-16O3__AIJS.param");
		}
	}
	if(m.id == 13){//OH
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "13_hit16.param");
		}
		if(param.useHITEMP == 1){
	        	sprintf(pFileName, "%s%s", param.path, "13_HITEMP2010.param");
		}
	}
	if(m.id == 15){//HCl
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "15_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H-35Cl__Yueqi.param");
		}
	}
	if(m.id == 20){//H2CO
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "20_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-12C-16O__AYTY.param");
		}
	}
	if(m.id == 23){//HCN
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "23_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H-12C-14N__Harris.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "1H-13C-14N__Larner.param");
		}
	}
	if(m.id == 25){//H2O2
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "25_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-16O2__APTY.param");
		}
	}
	if(m.id == 26){//C2H2
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "26_hit16.param");
		}
	}
	if(m.id == 28){//PH3
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "28_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "31P-1H3__SAlTY.param");
		}
	}
	if(m.id == 31){//H2S
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "31_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-32S__AYT2.param");
		}
	}
	if(m.id == 47){//SO3
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "47_hit16.param");
		}
	}
	if(m.id == 78){//CH3F
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-1H3-19F__OYKYT.param");
		}
	}
	if(m.id == 79){//SiH4
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-1H4__OY2T.param");
		}
	}
	if(m.id == 80){//VO
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "51V-16O__VOMYT.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczVO.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczVOax.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczVObx.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczVOcx.param");
		}
		if(param.useHITEMP == 7){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczVOmyt.param");
		}
	}
	if(m.id == 81){//TiO
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "curbestplusdiag_TiO_iota_48_Refined.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "Plez2012.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "Plez2012-norlander.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "Plez2012-philips.param");
		}
		if(param.useHITEMP == 7){
	        	sprintf(pFileName, "%s%s", param.path, "Plez2012-polfits.param");
		}
		if(param.useHITEMP == 8){
	        	sprintf(pFileName, "%s%s", param.path, "Schwenke1998.param");
		}
	}
	if(m.id == 82){//FeH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "56Fe-1H__Yueqi.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "KuruczFeHfx.param");
		}
	}
	if(m.id == 83){//AlO
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "27Al-16O__ATP.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "26Al-16O__ATP.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "27Al-17O__ATP.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "27Al-18O__ATP.param");
		}
	}
	if(m.id == 84){//SiO
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-16O__EBJT.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-17O__EBJT.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-18O__EBJT.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "29Si-16O__EBJT.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "30Si-16O__EBJT.param");
		}
	}
	if(m.id == 85){//CaO
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "40Ca-16O__VBATHY.param");
		}
	}
	if(m.id == 86){//SiH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-1H__SiGHTLY.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "29Si-1H__SiGHTLY.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "30Si-1H__SiGHTLY.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-2H__SiGHTLY.param");
		}
	}
	if(m.id == 87){//CaH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "40Ca-1H__Yadin.param");
		}
	}
	if(m.id == 88){//H3+
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "1H3_p__MiZATeP.param");
		}
	}
	if(m.id == 89){//PO
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "31P-16O__POPS.param");
		}
	}
	if(m.id == 90){//MgH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "24Mg-1H__Yadin.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "25Mg-1H__Yadin.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "26Mg-1H__Yadin.param");
		}
	}
	if(m.id == 91){//NaH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "23Na-1H__Rivlin.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "23Na-2H__Rivlin.param");
		}
	}
	if(m.id == 92){//AlH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "27Al-1H__AlHambra.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "27Al-2H__AlHambra.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "26Al-1H__AlHambra.param");
		}
	}
	if(m.id == 93){//CrH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "52Cr-1H__Yueqi.param");
		}
	}
	if(m.id == 94){//BeH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "9Be-1H__Yadin.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "9Be-1H__Darby-Lewis.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "9Be-2H__Darby-Lewis.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "9Be-3H__Darby-Lewis.param");
		}
	}
	if(m.id == 95){//TiH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "48Ti-1H__Yueqi.param");
		}
	}
	if(m.id == 97){//ScH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "45Sc-1H__LYT.param");
		}
	}
	if(m.id == 99){//NH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "14N-1H__Yueqi.param");
		}
	}
	if(m.id == 100){//CH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-1H__Yueqi.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "13C-1H__Yueqi.param");
		}
	}
	if(m.id == 101){//SH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "32S-1H__SNaSH.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "33S-1H__SNaSH.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "34S-1H__SNaSH.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "36S-1H__SNaSH.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "32S-2H__SNaSH.param");
		}
	}
	if(m.id == 102){//PN
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "31P-14N__YYLT.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "31P-15N__YYLT.param");
		}
	}
	if(m.id == 103){//KCl
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "39K-35Cl__Barton.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "39K-37Cl__Barton.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "41K-35Cl__Barton.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "41K-37Cl__Barton.param");
		}
	}
	if(m.id == 104){//Nal
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "23Na-35Cl__Barton.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "23Na-37Cl__Barton.param");
		}
	}
	if(m.id == 105){//NC
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-14N__Yueqi.param");
		}
	}
	if(m.id == 106){//C2
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C2__8states.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "12C-13C__8states.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "13C2__8states.param");
		}
	}
	if(m.id == 107){//CS
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-32S__JnK.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "12C-34S__JnK.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "13C-32S__JnK.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "12C-33S__JnK.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "12C-36S__JnK.param");
		}
		if(param.useHITEMP == 7){
	        	sprintf(pFileName, "%s%s", param.path, "13C-33S__JnK.param");
		}
		if(param.useHITEMP == 8){
	        	sprintf(pFileName, "%s%s", param.path, "13C-34S__JnK.param");
		}
		if(param.useHITEMP == 9){
	        	sprintf(pFileName, "%s%s", param.path, "13C-36S__JnK.param");
		}
	}
	if(m.id == 108){//CP
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-31P__Yueqi.param");
		}
	}
	if(m.id == 109){//PO
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "31P-32S__POPS.param");
		}
	}
	if(m.id == 110){//NS
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "14N-32S__SNaSH.param");
		}
	}
	if(m.id == 111){//SiS
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-32S__UCTY.param");
		}
	}
	if(m.id == 112){//HeH+
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "4He-1H_p__Engel.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "3He-1H_p__Engel.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "4He-2H_p__Engel.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "3He-2H_p__Engel.param");
		}
	}
	pFile = fopen(pFileName, "r");
	if(pFile == NULL){
		printf("Error, no molecule.param file %s\n", pFileName);
	}
	char sp[160];
	fgets(sp, 11, pFile);
	if(strcmp(sp, "Database =") != 0){
		printf("Error in molecule.param file, Database\n");
	}
	fscanf (pFile, "%d", &param.useHITEMP);
	fgets(sp, 4, pFile);

	fgets(sp, 18, pFile);
	if(strcmp(sp, "Molecule number =") != 0){
		printf("Error in molecule.param file, Molecule number\n");
	}
	fscanf (pFile, "%d", &m.id);
	fgets(sp, 4, pFile);

	fgets(sp, 7, pFile);
	if(strcmp(sp, "Name =") != 0){
		printf("Error in molecule.param file, Name\n");
	}
	fscanf (pFile, "%s", m.mName);
	fgets(sp, 4, pFile);
	fgets(sp, 21, pFile);
	if(strcmp(sp, "Number of Isotopes =") != 0){
		printf("Error in molecule.param file, Number of Isotopes\n");
	}
	fscanf (pFile, "%d", &m.nISO);
	fgets(sp, 4, pFile);
	m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));

	fgets(sp, 100, pFile);

	for(int i = 0; i < m.nISO; ++i){
		fscanf (pFile, "%s",  sp);
		sprintf(m.ISO[i].cid, "%3s", sp);
		fscanf (pFile, "%lf", &m.ISO[i].Ab);
		fscanf (pFile, "%lf", &m.ISO[i].Q);
		fscanf (pFile, "%d",  &m.ISO[i].g);
		fscanf (pFile, "%lf", &m.ISO[i].m);
		fscanf (pFile, "%s", sp);
		sprintf(qFilename[i], "%s%s", param.path, sp);
printf("%s %g %s\n", m.ISO[i].cid, m.ISO[i].Ab, qFilename[i]);
	}
	fgets(sp, 4, pFile);
	fgets(sp, 38, pFile);
	if(strcmp(sp, "Number of columns in partition File =") != 0){
		printf("Error in molecule.param file, Number of columns in partition File\n");
	}
	fscanf (pFile, "%d", &m.npfcol);
	fgets(sp, 4, pFile);

	fgets(sp, 34, pFile);
	if(strcmp(sp, "Number of line/transition files =") != 0){
		printf("Error in molecule.param file, Number of line/transition files\n");
	}
	fscanf (pFile, "%d", &m.nFiles);
	fgets(sp, 4, pFile);

	fgets(sp, 100, pFile);

	m.NLmax = 0;
	for(int i = 0; i < m.nFiles; ++i){
		fscanf (pFile, "%d", &m.NL[i]);
printf("%d\n", m.NL[i]);
		m.NLmax = max(m.NLmax, m.NL[i]);
	}
	fgets(sp, 4, pFile);

	fgets(sp, 19, pFile);
	for(int i = 0; i < m.nFiles + 1; ++i){
		fscanf (pFile, "%d", &m.fileLimit[i]);
printf("%d\n", m.fileLimit[i]);
	}
	for(int i = 0; i < m.nFiles; ++i){
		if(param.useHITEMP == 0){
			sprintf(m.dataFilename[i], "%s%02d_%s.", param.path, m.id, m.mName);
		}
		if(param.useHITEMP == 1){
			if(m.nFiles > 1){
				sprintf(m.dataFilename[i], "%s%02d_%05d-%05d_%s.", param.path, m.id, m.fileLimit[i], m.fileLimit[i + 1], m.mName);
			}
			else{
				sprintf(m.dataFilename[i], "%s%02d_%s.", param.path, m.id, m.mName);
			}
		}
		if(param.useHITEMP == 2){
			if(m.nFiles > 1){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			else{
				sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
			}
		}
		
	}
	fgets(sp, 4, pFile);
	fgets(sp, 100, pFile);
	fgets(sp, 19, pFile);
	if(strcmp(sp, "Number of states =") != 0){
		printf("Error in molecule.param file, Number of states\n");
	}
	fscanf (pFile, "%d", &m.nStates);
	fgets(sp, 4, pFile);

	fgets(sp, 40, pFile);
	if(strcmp(sp, "Number of columns in transition files =") != 0){
		printf("Error in molecule.param file, Number of columns in transition files\n");
	}
	fscanf (pFile, "%d", &m.ntcol);
	fgets(sp, 4, pFile);

	fgets(sp, 55, pFile);
	if(strcmp(sp, "Default value of Lorentzian half-width for all lines =") != 0){
		printf("Error in molecule.param file, Default value of Lorentzian half-width for all lines\n");
	}
	fscanf (pFile, "%lf", &m.defaultL);
	fgets(sp, 4, pFile);

	fgets(sp, 54, pFile);
	if(strcmp(sp, "Default value of temperature exponent for all lines =") != 0){
		printf("Error in molecule.param file, Default value of temperature exponent for all lines\n");
	}
	fscanf (pFile, "%lf", &m.defaultn);
	fgets(sp, 4, pFile);

	fgets(sp, 10, pFile);
	if(strcmp(sp, "Version =") != 0){
		printf("Error in molecule.param file, Version\n");
	}
	fscanf (pFile, "%d", &m.version);
	fgets(sp, 4, pFile);

	fclose(pFile);
	if(m.id == 300){//Li 7
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){300,  7,  1.0,    0.0,    0,     6.941};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0300";
			sprintf(m.mName, "%s", "gfnew0300");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2863;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 42720;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){300,  7,  1.0,    0.0,    0,     6.941};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 400){//Be 9
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){400,  9,  1.0,    0.0,    0,     9.01218};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0400";
			sprintf(m.mName, "%s", "gfnew0400");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 3832;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 106385;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){400,  9,  1.0,    0.0,    0,     9.01218};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 401){//Be+ 9
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){401,  9,  1.0,    0.0,    0,     9.01218};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0401";
			sprintf(m.mName, "%s", "gfnew0401");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 897;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 142450;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){401,  9,  1.0,    0.0,    0,     9.01218};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 500){//B 11
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){500,  11,  1.0,    0.0,    0,     10.811};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0500";
			sprintf(m.mName, "%s", "gfnew0500");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2753;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 100682;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){500,  11,  1.0,    0.0,    0,     10.811};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 501){//B+ 11
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){501,  11,  1.0,    0.0,    0,     10.811};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0501";
			sprintf(m.mName, "%s", "gfnew0501");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 5155;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 246859;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){501,  11,  1.0,    0.0,    0,     10.811};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 502){//B+2 11
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){502,  11,  1.0,    0.0,    0,     10.811};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0502";
			sprintf(m.mName, "%s", "gfnew0502");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 859;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 297766;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){502,  11,  1.0,    0.0,    0,     10.811};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 601){//C+ 12
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){601,  12,  1.0,    0.0,    0,     12.011};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0601";
			sprintf(m.mName, "%s", "gfnew0601");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 13887;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 236750;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){601,  12,  1.0,    0.0,    0,     12.011};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 602){//C+2 12
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){602,  12,  1.0,    0.0,    0,     12.011};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0602";
			sprintf(m.mName, "%s", "gfnew0602");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8673;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 423110;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){602,  12,  1.0,    0.0,    0,     12.011};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 700){//N 14
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){700,  14,  1.0,    0.0,    0,     14.0067};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0700";
			sprintf(m.mName, "%s", "gfnew0700");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 14522;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 163400;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){700,  14,  1.0,    0.0,    0,     14.0067};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 701){//N+ 14
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){701,  14,  1.0,    0.0,    0,     14.0067};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0701";
			sprintf(m.mName, "%s", "gfnew0701");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 4142;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 277952;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){701,  14,  1.0,    0.0,    0,     14.0067};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 702){//N+2 14
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){702,  14,  1.0,    0.0,    0,     14.0067};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0702";
			sprintf(m.mName, "%s", "gfnew0702");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 12772;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 479735;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){702,  14,  1.0,    0.0,    0,     14.0067};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 800){//O 16
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){800,  16,  1.0,    0.0,    0,     15.9994};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0800";
			sprintf(m.mName, "%s", "gfnew0800");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 13496;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 150032;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){800,  16,  1.0,    0.0,    0,     15.9994};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 801){//O+ 16
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){801,  16,  1.0,    0.0,    0,     15.9994};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0801";
			sprintf(m.mName, "%s", "gfnew0801");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9617;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 287048;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){801,  16,  1.0,    0.0,    0,     15.9994};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 802){//O+2 16
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){802,  16,  1.0,    0.0,    0,     15.9994};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0802";
			sprintf(m.mName, "%s", "gfnew0802");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8678;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 474689;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){802,  16,  1.0,    0.0,    0,     15.9994};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 900){//F 19
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){900,  19,  1.0,    0.0,    0,     18.9984};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0900";
			sprintf(m.mName, "%s", "gfnew0900");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 5463;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 171001;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){900,  19,  1.0,    0.0,    0,     18.9984};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 901){//F+ 19
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){901,  19,  1.0,    0.0,    0,     18.9984};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0901";
			sprintf(m.mName, "%s", "gfnew0901");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9367;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 290143;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){901,  19,  1.0,    0.0,    0,     18.9984};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 902){//F+2 19
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){902,  19,  1.0,    0.0,    0,     18.9984};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew0902";
			sprintf(m.mName, "%s", "gfnew0902");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9852;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 491593;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){902,  19,  1.0,    0.0,    0,     18.9984};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1000){//Ne 20
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1000,  20,  1.0,    0.0,    0,     20.1797};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1000";
			sprintf(m.mName, "%s", "gfnew1000");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 15362;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 390091;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1000,  20,  1.0,    0.0,    0,     20.1797};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1001){//Ne+ 20
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1001,  20,  1.0,    0.0,    0,     20.1797};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1001";
			sprintf(m.mName, "%s", "gfnew1001");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 18467;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 349275;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1001,  20,  1.0,    0.0,    0,     20.1797};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1002){//Ne+2 20
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1002,  20,  1.0,    0.0,    0,     20.1797};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1002";
			sprintf(m.mName, "%s", "gfnew1002");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9476;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 534501;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1002,  20,  1.0,    0.0,    0,     20.1797};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1100){//Na 23
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1100,  23,  1.0,    0.0,    0,     22.9898};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1100";
			sprintf(m.mName, "%s", "gfnew1100");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8677;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 312357;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1100,  23,  1.0,    0.0,    0,     22.9898};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1101){//Na+ 23
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1101,  23,  1.0,    0.0,    0,     22.9898};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1101";
			sprintf(m.mName, "%s", "gfnew1101");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 4337;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 637401;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1101,  23,  1.0,    0.0,    0,     22.9898};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1102){//Na+2 23
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1102,  23,  1.0,    0.0,    0,     22.9898};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1102";
			sprintf(m.mName, "%s", "gfnew1102");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1990;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 552417;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1102,  23,  1.0,    0.0,    0,     22.9898};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1200){//Mg 24
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1200,  24,  1.0,    0.0,    0,     24.305};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1200";
			sprintf(m.mName, "%s", "gfnew1200");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 11319;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 96437;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1200,  24,  1.0,    0.0,    0,     24.305};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1201){//Mg+ 24
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1201,  24,  1.0,    0.0,    0,     24.305};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1201";
			sprintf(m.mName, "%s", "gfnew1201");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2832;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 505781;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1201,  24,  1.0,    0.0,    0,     24.305};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1202){//Mg+2 24
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1202,  24,  1.0,    0.0,    0,     24.305};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1202";
			sprintf(m.mName, "%s", "gfnew1202");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1643;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 940701;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1202,  24,  1.0,    0.0,    0,     24.305};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1300){//Al 27
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1300,  27,  1.0,    0.0,    0,     26.9815};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1300";
			sprintf(m.mName, "%s", "gfnew1300");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 4296;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 84080;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1300,  27,  1.0,    0.0,    0,     26.9815};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1301){//Al+ 27
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1301,  27,  1.0,    0.0,    0,     26.9815};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1301";
			sprintf(m.mName, "%s", "gfnew1301");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 5799;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 199001;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1301,  27,  1.0,    0.0,    0,     26.9815};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1302){//Al+2 27
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1302,  27,  1.0,    0.0,    0,     26.9815};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1302";
			sprintf(m.mName, "%s", "gfnew1302");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2839;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 727501;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1302,  27,  1.0,    0.0,    0,     26.9815};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1400){//Si 28
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1400,  28,  1.0,    0.0,    0,     28.0855};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1400";
			sprintf(m.mName, "%s", "gfnew1400");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 10635;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 100431;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1400,  28,  1.0,    0.0,    0,     28.0855};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1401){//Si+ 28
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1401,  28,  1.0,    0.0,    0,     28.0855};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1401";
			sprintf(m.mName, "%s", "gfnew1401");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 3056;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 157484;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1401,  28,  1.0,    0.0,    0,     28.0855};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1402){//Si+2 28
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1402,  28,  1.0,    0.0,    0,     28.0855};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1402";
			sprintf(m.mName, "%s", "gfnew1402");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2974;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 244934;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1402,  28,  1.0,    0.0,    0,     28.0855};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1500){//P 31
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1500,  31,  1.0,    0.0,    0,     30.9738};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1500";
			sprintf(m.mName, "%s", "gfnew1500");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 12291;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 128339;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1500,  31,  1.0,    0.0,    0,     30.9738};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1501){//P+ 31
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1501,  31,  1.0,    0.0,    0,     30.9738};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1501";
			sprintf(m.mName, "%s", "gfnew1501");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2969;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 150889;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1501,  31,  1.0,    0.0,    0,     30.9738};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1502){//P+2 31
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1502,  31,  1.0,    0.0,    0,     30.9738};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1502";
			sprintf(m.mName, "%s", "gfnew1502");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2225;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 254723;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1502,  31,  1.0,    0.0,    0,     30.9738};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1600){//S 32
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1600,  32,  1.0,    0.0,    0,     32.066};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1600";
			sprintf(m.mName, "%s", "gfnew1600");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 24734;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 107646;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1600,  32,  1.0,    0.0,    0,     32.066};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1601){//S+ 32
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1601,  32,  1.0,    0.0,    0,     32.066};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1601";
			sprintf(m.mName, "%s", "gfnew1601");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 7297;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 184643;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1601,  32,  1.0,    0.0,    0,     32.066};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1602){//S+2 32
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1602,  32,  1.0,    0.0,    0,     32.066};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1602";
			sprintf(m.mName, "%s", "gfnew1602");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 445;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 238196;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1602,  32,  1.0,    0.0,    0,     32.066};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1700){//Cl 35
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1700,  35,  1.0,    0.0,    0,     35.4527};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1700";
			sprintf(m.mName, "%s", "gfnew1700");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 17530;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 131792;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1700,  35,  1.0,    0.0,    0,     35.4527};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1701){//Cl+ 35
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1701,  35,  1.0,    0.0,    0,     35.4527};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1701";
			sprintf(m.mName, "%s", "gfnew1701");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9593;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 188248;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1701,  35,  1.0,    0.0,    0,     35.4527};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1702){//Cl+2 35
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1702,  35,  1.0,    0.0,    0,     35.4527};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1702";
			sprintf(m.mName, "%s", "gfnew1702");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 968;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 258891;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1702,  35,  1.0,    0.0,    0,     35.4527};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1800){//Ar 40
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1800,  40,  1.0,    0.0,    0,     39.948};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1800";
			sprintf(m.mName, "%s", "gfnew1800");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 16650;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 234976;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1800,  40,  1.0,    0.0,    0,     39.948};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1801){//Ar+ 40
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1801,  40,  1.0,    0.0,    0,     39.948};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1801";
			sprintf(m.mName, "%s", "gfnew1801");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 17190;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 219247;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1801,  40,  1.0,    0.0,    0,     39.948};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1802){//Ar+2 40
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1802,  40,  1.0,    0.0,    0,     39.948};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1802";
			sprintf(m.mName, "%s", "gfnew1802");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1923;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 286010;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1802,  40,  1.0,    0.0,    0,     39.948};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1900){//K 39
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1900,  39,  1.0,    0.0,    0,     39.0983};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1900";
			sprintf(m.mName, "%s", "gfnew1900");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 3123;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 179888;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1900,  39,  1.0,    0.0,    0,     39.0983};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1901){//K+ 39
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1901,  39,  1.0,    0.0,    0,     39.0983};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1901";
			sprintf(m.mName, "%s", "gfnew1901");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1125;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 364601;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1901,  39,  1.0,    0.0,    0,     39.0983};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 1902){//K+2 39
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){1902,  39,  1.0,    0.0,    0,     39.0983};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew1902";
			sprintf(m.mName, "%s", "gfnew1902");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 213;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 250859;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){1902,  39,  1.0,    0.0,    0,     39.0983};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2000){//Ca 40
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2000,  40,  1.0,    0.0,    0,     40.078};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2000";
			sprintf(m.mName, "%s", "gfnew2000");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 25410;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 72290;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2000,  40,  1.0,    0.0,    0,     40.078};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2001){//Ca+ 40
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2001,  40,  1.0,    0.0,    0,     40.078};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2001";
			sprintf(m.mName, "%s", "gfnew2001");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 3467;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 319401;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2001,  40,  1.0,    0.0,    0,     40.078};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2002){//Ca+2 40
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2002,  40,  1.0,    0.0,    0,     40.078};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2002";
			sprintf(m.mName, "%s", "gfnew2002");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2973;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 492851;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2002,  40,  1.0,    0.0,    0,     40.078};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2100){//Sc 45
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2100,  45,  1.0,    0.0,    0,     44.9559};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2100";
			sprintf(m.mName, "%s", "gfnew2100");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 16252;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 63730;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2100,  45,  1.0,    0.0,    0,     44.9559};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2101){//Sc+ 45
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2101,  45,  1.0,    0.0,    0,     44.9559};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2101";
			sprintf(m.mName, "%s", "gfnew2101");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 5402;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 91166;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2101,  45,  1.0,    0.0,    0,     44.9559};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2102){//Sc+2 45
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2102,  45,  1.0,    0.0,    0,     44.9559};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2102";
			sprintf(m.mName, "%s", "gfnew2102");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1313;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 194595;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2102,  45,  1.0,    0.0,    0,     44.9559};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2200){//Ti 48
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2200,  48,  1.0,    0.0,    0,     47.88};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2200";
			sprintf(m.mName, "%s", "gfnew2200");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 36050;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 49359;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2200,  48,  1.0,    0.0,    0,     47.88};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2201){//Ti+ 48
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2201,  48,  1.0,    0.0,    0,     47.88};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2201";
			sprintf(m.mName, "%s", "gfnew2201");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9318;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 82371;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2201,  48,  1.0,    0.0,    0,     47.88};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2202){//Ti+2 48
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2202,  48,  1.0,    0.0,    0,     47.88};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2202";
			sprintf(m.mName, "%s", "gfnew2202");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 4179;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 201982;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2202,  48,  1.0,    0.0,    0,     47.88};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2300){//V 51
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2300,  51,  1.0,    0.0,    0,     50.9415};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2300";
			sprintf(m.mName, "%s", "gfnew2300");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 211129;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 52637;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2300,  51,  1.0,    0.0,    0,     50.9415};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2301){//V+ 51
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2301,  51,  1.0,    0.0,    0,     50.9415};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2301";
			sprintf(m.mName, "%s", "gfnew2301");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 21482;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 93274;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2301,  51,  1.0,    0.0,    0,     50.9415};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2302){//V+2 51
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2302,  51,  1.0,    0.0,    0,     50.9415};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2302";
			sprintf(m.mName, "%s", "gfnew2302");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 10318;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 191598;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2302,  51,  1.0,    0.0,    0,     50.9415};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2400){//Cr 52
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2400,  52,  1.0,    0.0,    0,     51.9961};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2400";
			sprintf(m.mName, "%s", "gfnew2400");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 38788;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 66094;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2400,  52,  1.0,    0.0,    0,     51.9961};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2401){//Cr+ 52
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2401,  52,  1.0,    0.0,    0,     51.9961};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2401";
			sprintf(m.mName, "%s", "gfnew2401");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 95312;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 121333;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2401,  52,  1.0,    0.0,    0,     51.9961};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2402){//Cr+2 52
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2402,  52,  1.0,    0.0,    0,     51.9961};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2402";
			sprintf(m.mName, "%s", "gfnew2402");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 16013;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 179978;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2402,  52,  1.0,    0.0,    0,     51.9961};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2500){//Mn 55
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2500,  55,  1.0,    0.0,    0,     54.938};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2500";
			sprintf(m.mName, "%s", "gfnew2500");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 42891;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 68339;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2500,  55,  1.0,    0.0,    0,     54.938};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2501){//Mn+ 55
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2501,  55,  1.0,    0.0,    0,     54.938};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2501";
			sprintf(m.mName, "%s", "gfnew2501");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 61276;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 117399;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2501,  55,  1.0,    0.0,    0,     54.938};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2502){//Mn+2 55
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2502,  55,  1.0,    0.0,    0,     54.938};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2502";
			sprintf(m.mName, "%s", "gfnew2502");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 18145;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 232265;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2502,  55,  1.0,    0.0,    0,     54.938};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2600){//Fe 56
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2600,  56,  1.0,    0.0,    0,     55.847};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2600";
			sprintf(m.mName, "%s", "gfnew2600");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 127897;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 65591;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2600,  56,  1.0,    0.0,    0,     55.847};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2601){//Fe+ 56
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2601,  56,  1.0,    0.0,    0,     55.847};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2601";
			sprintf(m.mName, "%s", "gfnew2601");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 127757;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 132155;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2601,  56,  1.0,    0.0,    0,     55.847};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2602){//Fe+2 56
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2602,  56,  1.0,    0.0,    0,     55.847};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2602";
			sprintf(m.mName, "%s", "gfnew2602");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 37795;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 219781;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2602,  56,  1.0,    0.0,    0,     55.847};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2700){//Co 59
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2700,  59,  1.0,    0.0,    0,     58.9332};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2700";
			sprintf(m.mName, "%s", "gfnew2700");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 249130;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 59948;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2700,  59,  1.0,    0.0,    0,     58.9332};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2701){//Co+ 59
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2701,  59,  1.0,    0.0,    0,     58.9332};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2701";
			sprintf(m.mName, "%s", "gfnew2701");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 24515;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 112008;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2701,  59,  1.0,    0.0,    0,     58.9332};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2702){//Co+2 59
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2702,  59,  1.0,    0.0,    0,     58.9332};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2702";
			sprintf(m.mName, "%s", "gfnew2702");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9706;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 195704;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2702,  59,  1.0,    0.0,    0,     58.9332};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2800){//Ni 59
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2800,  59,  1.0,    0.0,    0,     58.6934};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2800";
			sprintf(m.mName, "%s", "gfnew2800");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 16804;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 58897;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2800,  59,  1.0,    0.0,    0,     58.6934};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2801){//Ni+ 59
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2801,  59,  1.0,    0.0,    0,     58.6934};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2801";
			sprintf(m.mName, "%s", "gfnew2801");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 56546;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 138842;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2801,  59,  1.0,    0.0,    0,     58.6934};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2802){//Ni+2 59
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2802,  59,  1.0,    0.0,    0,     58.6934};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2802";
			sprintf(m.mName, "%s", "gfnew2802");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 21415;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 229781;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2802,  59,  1.0,    0.0,    0,     58.6934};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2900){//Cu 64
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2900,  64,  1.0,    0.0,    0,     63.546};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2900";
			sprintf(m.mName, "%s", "gfnew2900");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 18087;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 88019;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2900,  64,  1.0,    0.0,    0,     63.546};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2901){//Cu+ 64
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2901,  64,  1.0,    0.0,    0,     63.546};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2901";
			sprintf(m.mName, "%s", "gfnew2901");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 15077;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 153458;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2901,  64,  1.0,    0.0,    0,     63.546};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2902){//Cu+2 64
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2902,  64,  1.0,    0.0,    0,     63.546};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew2902";
			sprintf(m.mName, "%s", "gfnew2902");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 17590;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 261763;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){2902,  64,  1.0,    0.0,    0,     63.546};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3000){//Zn 65
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3000,  65,  1.0,    0.0,    0,     65.39};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew3000";
			sprintf(m.mName, "%s", "gfnew3000");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 6282;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 166419;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){3000,  65,  1.0,    0.0,    0,     65.39};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3001){//Zn+ 65
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3001,  65,  1.0,    0.0,    0,     65.39};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew3001";
			sprintf(m.mName, "%s", "gfnew3001");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 968;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 143243;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){3001,  65,  1.0,    0.0,    0,     65.39};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3002){//Zn+2 65
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3002,  65,  1.0,    0.0,    0,     65.39};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew3002";
			sprintf(m.mName, "%s", "gfnew3002");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 12681;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 300009;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){3002,  65,  1.0,    0.0,    0,     65.39};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3800){//Sr 88
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3800,  88,  1.0,    0.0,    0,     87.62};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew3800";
			sprintf(m.mName, "%s", "gfnew3800");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 22776;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 61872;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){3800,  88,  1.0,    0.0,    0,     87.62};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3801){//Sr+ 88
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3801,  88,  1.0,    0.0,    0,     87.62};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew3801";
			sprintf(m.mName, "%s", "gfnew3801");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 674;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 75312;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){3801,  88,  1.0,    0.0,    0,     87.62};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3900){//Y 89
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3900,  89,  1.0,    0.0,    0,     88.9059};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew3900";
			sprintf(m.mName, "%s", "gfnew3900");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 5654;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 49043;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){3900,  89,  1.0,    0.0,    0,     88.9059};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3901){//Y+ 89
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3901,  89,  1.0,    0.0,    0,     88.9059};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew3901";
			sprintf(m.mName, "%s", "gfnew3901");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 7588;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 90213;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){3901,  89,  1.0,    0.0,    0,     88.9059};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4000){//Zr 91
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4000,  91,  1.0,    0.0,    0,     91.224};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4000";
			sprintf(m.mName, "%s", "gfnew4000");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 6200;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 51900;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4000,  91,  1.0,    0.0,    0,     91.224};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4001){//Zr+ 91
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4001,  91,  1.0,    0.0,    0,     91.224};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4001";
			sprintf(m.mName, "%s", "gfnew4001");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1834;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 61862;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4001,  91,  1.0,    0.0,    0,     91.224};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4002){//Zr+2 91
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4002,  91,  1.0,    0.0,    0,     91.224};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4002";
			sprintf(m.mName, "%s", "gfnew4002");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2360;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 158488;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4002,  91,  1.0,    0.0,    0,     91.224};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4100){//Nb 93
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4100,  93,  1.0,    0.0,    0,     92.9064};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4100";
			sprintf(m.mName, "%s", "gfnew4100");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 117714;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 51094;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4100,  93,  1.0,    0.0,    0,     92.9064};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4101){//Nb+ 93
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4101,  93,  1.0,    0.0,    0,     92.9064};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4101";
			sprintf(m.mName, "%s", "gfnew4101");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 28653;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 78371;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4101,  93,  1.0,    0.0,    0,     92.9064};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4102){//Nb+2 93
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4102,  93,  1.0,    0.0,    0,     92.9064};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4102";
			sprintf(m.mName, "%s", "gfnew4102");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 4009;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 128012;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4102,  93,  1.0,    0.0,    0,     92.9064};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4200){//Mo 96
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4200,  96,  1.0,    0.0,    0,     95.94};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4200";
			sprintf(m.mName, "%s", "gfnew4200");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 13862;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 59237;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4200,  96,  1.0,    0.0,    0,     95.94};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4201){//Mo+ 96
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4201,  96,  1.0,    0.0,    0,     95.94};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4201";
			sprintf(m.mName, "%s", "gfnew4201");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 13272;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 87034;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4201,  96,  1.0,    0.0,    0,     95.94};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4202){//Mo+2 96
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4202,  96,  1.0,    0.0,    0,     95.94};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4202";
			sprintf(m.mName, "%s", "gfnew4202");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 12387;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 163672;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4202,  96,  1.0,    0.0,    0,     95.94};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4300){//Tc 98
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4300,  98,  1.0,    0.0,    0,     97.9072};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4300";
			sprintf(m.mName, "%s", "gfnew4300");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8815;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 52376;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4300,  98,  1.0,    0.0,    0,     97.9072};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4301){//Tc+ 98
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4301,  98,  1.0,    0.0,    0,     97.9072};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4301";
			sprintf(m.mName, "%s", "gfnew4301");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 119;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 61420;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4301,  98,  1.0,    0.0,    0,     97.9072};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4400){//Ru 101
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4400,  101,  1.0,    0.0,    0,     101.07};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4400";
			sprintf(m.mName, "%s", "gfnew4400");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8383;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 53718;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4400,  101,  1.0,    0.0,    0,     101.07};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4401){//Ru+ 101
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4401,  101,  1.0,    0.0,    0,     101.07};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4401";
			sprintf(m.mName, "%s", "gfnew4401");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 5340;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 90166;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4401,  101,  1.0,    0.0,    0,     101.07};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4402){//Ru+2 101
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4402,  101,  1.0,    0.0,    0,     101.07};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4402";
			sprintf(m.mName, "%s", "gfnew4402");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 70;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 83998;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4402,  101,  1.0,    0.0,    0,     101.07};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4500){//Rh 103
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4500,  103,  1.0,    0.0,    0,     102.906};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4500";
			sprintf(m.mName, "%s", "gfnew4500");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2332;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 52066;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4500,  103,  1.0,    0.0,    0,     102.906};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4501){//Rh+ 103
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4501,  103,  1.0,    0.0,    0,     102.906};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4501";
			sprintf(m.mName, "%s", "gfnew4501");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1527;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 87395;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4501,  103,  1.0,    0.0,    0,     102.906};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4502){//Rh+2 103
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4502,  103,  1.0,    0.0,    0,     102.906};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4502";
			sprintf(m.mName, "%s", "gfnew4502");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 3969;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 135855;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4502,  103,  1.0,    0.0,    0,     102.906};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4600){//Pd 106
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4600,  106,  1.0,    0.0,    0,     106.42};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4600";
			sprintf(m.mName, "%s", "gfnew4600");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 2996;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 63849;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4600,  106,  1.0,    0.0,    0,     106.42};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4601){//Pd+ 106
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4601,  106,  1.0,    0.0,    0,     106.42};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew4601";
			sprintf(m.mName, "%s", "gfnew4601");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 4558;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 128423;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){4601,  106,  1.0,    0.0,    0,     106.42};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 5600){//Ba 137
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){5600,  137,  1.0,    0.0,    0,     137.327};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew5600";
			sprintf(m.mName, "%s", "gfnew5600");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 9218;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 41184;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){5600,  137,  1.0,    0.0,    0,     137.327};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 5601){//Ba+ 137
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){5601,  137,  1.0,    0.0,    0,     137.327};
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		m.npfcol = 0;
		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 0;
		sprintf(m.dataFilename[0], "%s%s", param.path, ".");
		if(param.useHITEMP == 0){
			printf("Error: no Hitran data for this molecule\n");
		}
		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
		if(param.useHITEMP == 30){
			char name[] = "gfnew5601";
			sprintf(m.mName, "%s", "gfnew5601");
			m.defaultL = 0.0;
			m.defaultn = 0.0;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1956;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 75946;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){5601,  137,  1.0,    0.0,    0,     137.327};
			//version =  gfallwn08oct17.dat
		}
	}

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

