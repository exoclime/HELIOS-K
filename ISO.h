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
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "02_12C-16O2__HITEMP2010.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "02_13C-16O2__HITEMP2010.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "02_16O-12C-18O__HITEMP2010.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "02_16O-12C-17O__HITEMP2010.param");
		}
		if(param.useHITEMP == 7){
	        	sprintf(pFileName, "%s%s", param.path, "02_16O-13C-18O__HITEMP2010.param");
		}
		if(param.useHITEMP == 8){
	        	sprintf(pFileName, "%s%s", param.path, "02_16O-13C-17O__HITEMP2010.param");
		}
		if(param.useHITEMP == 9){
	        	sprintf(pFileName, "%s%s", param.path, "02_12C-18O2__HITEMP2010.param");
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
		if(param.useHITEMP == 5){
			sprintf(pFileName, "%s%s", param.path, "05_12C-16O__HITEMP2010.param");
		}
		if(param.useHITEMP == 6){
			sprintf(pFileName, "%s%s", param.path, "05_13C-16O__HITEMP2010.param");
		}
		if(param.useHITEMP == 7){
			sprintf(pFileName, "%s%s", param.path, "05_12C-18O__HITEMP2010.param");
		}
		if(param.useHITEMP == 8){
			sprintf(pFileName, "%s%s", param.path, "05_12C-17O__HITEMP2010.param");
		}
		if(param.useHITEMP == 9){
			sprintf(pFileName, "%s%s", param.path, "05_13C-18O__HITEMP2010.param");
		}
		if(param.useHITEMP == 10){
			sprintf(pFileName, "%s%s", param.path, "05_13C-17O__HITEMP2010.param");
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
		if(param.useHITEMP == 8){
	        	sprintf(pFileName, "%s%s", param.path, "08_14N-16O__HITEMP2010.param");
		}
		if(param.useHITEMP == 9){
	        	sprintf(pFileName, "%s%s", param.path, "08_15N-16O__HITEMP2010.param");
		}
		if(param.useHITEMP == 10){
	        	sprintf(pFileName, "%s%s", param.path, "08_14N-18O__HITEMP2010.param");
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
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "13_16O-1H__HITEMP2010.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "13_18O-1H__HITEMP2010.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "13_16O-2H__HITEMP2010.param");
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
	if(m.id == 24){//CH3Cl
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C-1H3-35Cl__OYT.param");
		}
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "12C-1H3-37Cl__OYT.param");
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
	if(m.id == 38){//C2H4
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "12C2-1H4__MaYTY.param");
		}
	}
	if(m.id == 46){//CS
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
	if(m.id == 47){//SO3
		if(param.useHITEMP == 0){
	        	sprintf(pFileName, "%s%s", param.path, "47_hit16.param");
		}
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "32S-16O3__UYT2.param");
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
		if(param.useHITEMP == 9){
	        	sprintf(pFileName, "%s%s", param.path, "curbestplusdiag_TiO_iota_48_Refined_PGopher.param");
		}
		if(param.useHITEMP == 10){
	        	sprintf(pFileName, "%s%s", param.path, "curbestplusdiag_TiO_iota_48_Refined_PGopher_scaled.param");
		}
		if(param.useHITEMP == 11){
	        	sprintf(pFileName, "%s%s", param.path, "curbestplusdiag_TiO_ELv1.0D_48.param");
		}
		if(param.useHITEMP == 12){
	        	sprintf(pFileName, "%s%s", param.path, "TiO_ELv1.0D_48.param");
		}
		if(param.useHITEMP == 13){
	        	sprintf(pFileName, "%s%s", param.path, "TiO_ELv1.0M_48.param");
		}
		if(param.useHITEMP == 14){
	        	sprintf(pFileName, "%s%s", param.path, "TiO_ELv1.0P_48.param");
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
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "1H2-2H_p__ST.param");
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
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "24Mg-1H__Yueqi.param");
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
	if(m.id == 96){//LiH
		if(param.useHITEMP == 2){
	        	sprintf(pFileName, "%s%s", param.path, "7Li-1H__CLT.param");
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
		if(param.useHITEMP == 3){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-33S__UCTY.param");
		}
		if(param.useHITEMP == 4){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-34S__UCTY.param");
		}
		if(param.useHITEMP == 5){
	        	sprintf(pFileName, "%s%s", param.path, "28Si-36S__UCTY.param");
		}
		if(param.useHITEMP == 6){
	        	sprintf(pFileName, "%s%s", param.path, "29Si-32S__UCTY.param");
		}
		if(param.useHITEMP == 7){
	        	sprintf(pFileName, "%s%s", param.path, "29Si-34S__UCTY.param");
		}
		if(param.useHITEMP == 8){
	        	sprintf(pFileName, "%s%s", param.path, "29Si-36S__UCTY.param");
		}
		if(param.useHITEMP == 9){
	        	sprintf(pFileName, "%s%s", param.path, "30Si-32S__UCTY.param");
		}
		if(param.useHITEMP == 10){
	        	sprintf(pFileName, "%s%s", param.path, "30Si-33S__UCTY.param");
		}
		if(param.useHITEMP == 11){
	        	sprintf(pFileName, "%s%s", param.path, "30Si-34S__UCTY.param");
		}
		if(param.useHITEMP == 12){
	        	sprintf(pFileName, "%s%s", param.path, "30Si-36S__UCTY.param");
		}
		if(param.useHITEMP == 13){
	        	sprintf(pFileName, "%s%s", param.path, "29Si-33S__UCTY.param");
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

	//Atoms and ions from Kurucz
	if(param.useHITEMP == 30){
        	sprintf(pFileName, "%sgfnew%04d.param", param.path, m.id);
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
//printf("%s %g %s\n", m.ISO[i].cid, m.ISO[i].Ab, qFilename[i]);
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
		m.NLmax = max(m.NLmax, m.NL[i]);
	}
	fgets(sp, 4, pFile);

	fgets(sp, 19, pFile);
	for(int i = 0; i < m.nFiles + 1; ++i){
		fscanf (pFile, "%d", &m.fileLimit[i]);
//printf("%d\n", m.fileLimit[i]);
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
		if(param.useHITEMP == 30){
			sprintf(m.dataFilename[i], "%s%s.", param.path, m.mName);
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

