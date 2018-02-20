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
void Init(Molecule &m, Param param, char *qFilename){
	if(m.id == 1){//H2O
		m.nISO = 6;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){11,  161,	.997317E+00,    1.7464E+02,    1,     18.010565};
		m.ISO[1] = (Isotopologue){12,  181,	1.99983E-03,    1.7511E+02,    1,     20.014811};
		m.ISO[2] = (Isotopologue){13,  171,	3.71884E-04,    1.0479E+03,    6,     19.014780};
		m.ISO[3] = (Isotopologue){14,  162,	3.10693E-04,    8.5901E+02,    6,     19.016740};
		m.ISO[4] = (Isotopologue){15,  182,	6.23003E-07,    8.7519E+02,    6,     21.020985};
		m.ISO[5] = (Isotopologue){16,  172,	1.15853E-07,    5.2204E+03,   36,     20.020956};
		m.nFiles = 1;		//number of data files
		//HITRAN2012
		m.NL[0] = 224515;	//number of lines
		m.NLmax = 224515;	//same as the number of lines

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 30000;
		sprintf(m.dataFilename[0], "%s%s", param.path, "01_hit12.");
		//HITRAN2008
		//m.NL[0] = 69201;	//number of lines
		//m.NLmax = 69201;	//same as the number of lines
		//sprintf(m.dataFilename[0], "%s%s", param.path, "01_hit08.");
		
		if(param.useHITEMP == 1){
			m.nFiles = 34;			//number of data files
			m.NL[ 0] = 2048387;		//number of lines per data file
			m.NL[ 1] = 3470074;
			m.NL[ 2] = 3625575;
			m.NL[ 3] = 3500494;
			m.NL[ 4] = 4972481;
			m.NL[ 5] = 3561630;
			m.NL[ 6] = 3795745;
			m.NL[ 7] = 3907835;
			m.NL[ 8] = 3969887;
			m.NL[ 9] = 3820310;
			m.NL[10] = 4960367;
			m.NL[11] = 4163192;
			m.NL[12] = 4363563;
			m.NL[13] = 4336535;
			m.NL[14] = 3758715;
			m.NL[15] = 3283806;
			m.NL[16] = 2971465;
			m.NL[17] = 2992331;
			m.NL[18] = 3076556;
			m.NL[19] = 2708657;
			m.NL[20] = 2299459;
			m.NL[21] = 4583911;
			m.NL[22] = 2421691;
			m.NL[23] = 2834046;
			m.NL[24] = 4079576;
			m.NL[25] = 3784869;
			m.NL[26] = 3452571;
			m.NL[27] = 2680328;
			m.NL[28] = 3154700;
			m.NL[29] = 2671165;
			m.NL[30] = 2220333;
			m.NL[31] = 1735773;
			m.NL[32] = 3573701;
			m.NL[33] = 1461436;

			m.NLmax = 0;
			for(int i = 0; i < m.nFiles; ++i){
				m.NLmax = max(m.NLmax, m.NL[i]);
			}

			m.fileLimit[ 0] = 0;
			m.fileLimit[ 1] = 50;
			m.fileLimit[ 2] = 150;
			m.fileLimit[ 3] = 250;
			m.fileLimit[ 4] = 350;
			m.fileLimit[ 5] = 500;
			m.fileLimit[ 6] = 600;
			m.fileLimit[ 7] = 700;
			m.fileLimit[ 8] = 800;
			m.fileLimit[ 9] = 900;
			m.fileLimit[10] = 1000;
			m.fileLimit[11] = 1150;
			m.fileLimit[12] = 1300;
			m.fileLimit[13] = 1500;
			m.fileLimit[14] = 1750;
			m.fileLimit[15] = 2000;
			m.fileLimit[16] = 2250;
			m.fileLimit[17] = 2500;
			m.fileLimit[18] = 2750;
			m.fileLimit[19] = 3000;
			m.fileLimit[20] = 3250;
			m.fileLimit[21] = 3500;
			m.fileLimit[22] = 4150;
			m.fileLimit[23] = 4500;
			m.fileLimit[24] = 5000;
			m.fileLimit[25] = 5500;
			m.fileLimit[26] = 6000;
			m.fileLimit[27] = 6500;
			m.fileLimit[28] = 7000;
			m.fileLimit[29] = 7500;
			m.fileLimit[30] = 8000;
			m.fileLimit[31] = 8500;
			m.fileLimit[32] = 9000;
			m.fileLimit[33] = 11000;
			m.fileLimit[34] = 30000;

			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s01_%05d-%05d_HITEMP2010.", param.path, m.fileLimit[i], m.fileLimit[i + 1]);
			}
		}
		if(param.useHITEMP == 2){
			char name[] = "1H2-16O__BT2";
			m.nStates = 221097;
			m.nFiles = 16;
			m.ntcol = 3;
			m.defaultL = 0.07;
			m.defaultn = 0.5;

			m.NL[0] = 17490214;
			m.NL[1] = 17022667;
			m.NL[2] = 16530697;
			m.NL[3] = 16098481;
			m.NL[4] = 30866787;
			m.NL[5] = 29161189;
			m.NL[6] = 13954798;
			m.NL[7] = 26727622;
			m.NL[8] = 37249656;
			m.NL[9] = 44635822;
			m.NL[10] = 39325124;
			m.NL[11] = 50083781;
			m.NL[12] = 52289428;
			m.NL[13] = 76679377;
			m.NL[14] = 31640191;
			m.NL[15] = 6050421;

			m.fileLimit[ 0] = 0;
			m.fileLimit[ 1] = 250; 
			m.fileLimit[ 2] = 500;
			m.fileLimit[ 3] = 750;
			m.fileLimit[ 4] = 1000;
			m.fileLimit[ 5] = 1500;
			m.fileLimit[ 6] = 2000;
			m.fileLimit[ 7] = 2250;
			m.fileLimit[ 8] = 2750;
			m.fileLimit[ 9] = 3500;
			m.fileLimit[10] = 4500;
			m.fileLimit[11] = 5500;
			m.fileLimit[12] = 7000;
			m.fileLimit[13] = 9000;
			m.fileLimit[14] = 14000;
			m.fileLimit[15] = 20000; 
			m.fileLimit[16] = 30000; 
			
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename, "%s%s%s", param.path, name, ".pf");

			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){11,  161,	1.0,    0.0,    0,     18.010565};
		}


	}
	if(m.id == 2){//CO2
		m.nFiles = 1;		//number of data files
		m.NL[0] = 471847;
		m.NLmax = 471847;
		m.nISO = 10;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){21,  626,	.984204E+00,    2.8694E+02,    1,     43.989830};
		m.ISO[1] = (Isotopologue){22,  636,	1.10574E-02,    5.7841E+02,    2,     44.993185};
		m.ISO[2] = (Isotopologue){23,  628,	3.94707E-03,    6.0948E+02,    1,     45.994076};
		m.ISO[3] = (Isotopologue){24,  627,	7.33989E-04,    3.5527E+03,    6,     44.994045};
		m.ISO[4] = (Isotopologue){25,  638,	4.43446E-05,    1.2291E+03,    2,     46.997431};
		m.ISO[5] = (Isotopologue){26,  637,	8.24623E-06,    7.1629E+03,   12,     45.997400};
		m.ISO[6] = (Isotopologue){27,  828,	3.95734E-06,    3.2421E+02,    1,     47.998322};
		m.ISO[7] = (Isotopologue){28,  728,	1.47180E-06,    3.7764E+03,    6,     46.998291};
		m.ISO[8] = (Isotopologue){29,  727,	1.36847E-07,    1.1002E+04,    1,     45.998262};
		m.ISO[9] = (Isotopologue){20,  838,	4.44600E-08,    6.5350E+02,    2,     49.001675};
		//m.ISO[0] = (Isotopologue){31,  837	1.65354E-08,    7.6152E+03,   12,     48.001646};

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 12785;
		sprintf(m.dataFilename[0], "%s%s", param.path, "02_hit12.");

		if(param.useHITEMP == 1){
			m.nFiles = 20;			//number of data files
			m.NL[ 0] = 259836;		//number of lines per data file
			m.NL[ 1] = 730397;
			m.NL[ 2] = 1140808;
			m.NL[ 3] = 943281;
			m.NL[ 4] = 378616;
			m.NL[ 5] = 456531;
			m.NL[ 6] = 576241;
			m.NL[ 7] = 1103020;
			m.NL[ 8] = 507509;
			m.NL[ 9] = 212373;
			m.NL[10] = 214133;
			m.NL[11] = 939432;
			m.NL[12] = 1046856;
			m.NL[13] = 213769;
			m.NL[14] = 138258;
			m.NL[15] = 866153;
			m.NL[16] = 450170;
			m.NL[17] = 177444;
			m.NL[18] = 344265;
			m.NL[19] = 494516;
			m.NLmax = 0;

			for(int i = 0; i < m.nFiles; ++i){
				m.NLmax = max(m.NLmax, m.NL[i]);
			}


			m.fileLimit[ 0] = 0;
			m.fileLimit[ 1] = 500;
			m.fileLimit[ 2] = 625;
			m.fileLimit[ 3] = 750;
			m.fileLimit[ 4] = 1000;
			m.fileLimit[ 5] = 1500;
			m.fileLimit[ 6] = 2000;
			m.fileLimit[ 7] = 2125;
			m.fileLimit[ 8] = 2250;
			m.fileLimit[ 9] = 2500;
			m.fileLimit[10] = 3000;
			m.fileLimit[11] = 3250;
			m.fileLimit[12] = 3500;
			m.fileLimit[13] = 3750;
			m.fileLimit[14] = 4000;
			m.fileLimit[15] = 4500;
			m.fileLimit[16] = 5000;
			m.fileLimit[17] = 5500;
			m.fileLimit[18] = 6000;
			m.fileLimit[19] = 6500;
			m.fileLimit[20] = 12785;
	
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s02_%05d-%05d_HITEMP2010.", param.path, m.fileLimit[i], m.fileLimit[i + 1]);
			}
		}

	}
	if(m.id == 3){//O3
		m.nFiles = 1;		//number of data files
		m.NL[0] = 422116;
		m.NLmax = 422116;
		m.nISO = 5;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){31,  666,  .992901E+00,    3.4838E+03,    1,     47.984745};
		m.ISO[1] = (Isotopologue){32,  668,  3.98194E-03,    7.4657E+03,    1,     49.988991};
		m.ISO[2] = (Isotopologue){33,  686,  1.99097E-03,    3.6471E+03,    1,     49.988991};
		m.ISO[3] = (Isotopologue){34,  667,  7.40475E-04,    4.3331E+04,    6,     48.988960};
		m.ISO[4] = (Isotopologue){35,  676,  3.70237E-04,    2.1405E+04,    6,     48.988960};

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 6997;

		sprintf(m.dataFilename[0], "%s%s", param.path, "03_hit12.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
	}
	if(m.id == 4){//N20
		m.nFiles = 1;		//number of data files
		m.NL[0] = 47843;
		m.NLmax = 47843;
		m.nISO = 5;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){41,  446,  .990333E+00,    5.0018E+03,    9,     44.001062};
		m.ISO[1] = (Isotopologue){42,  456,  3.64093E-03,    3.3619E+03,    6,     44.998096};
		m.ISO[2] = (Isotopologue){43,  546,  3.64093E-03,    3.4586E+03,    6,     44.998096};
		m.ISO[3] = (Isotopologue){44,  448,  1.98582E-03,    5.3147E+03,    9,     46.005308};
		m.ISO[4] = (Isotopologue){45,  447,  3.69280E-04,    3.0971E+04,   54,     45.005278};

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 7797;

		sprintf(m.dataFilename[0], "%s%s", param.path, "04_hit08.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
	}
	if(m.id == 5){//CO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 4606;
		m.NLmax = 4606;
		m.nISO = 6;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){51,  26,	.986544E+00,    1.0712E+02,    1,     27.994915};
		m.ISO[1] = (Isotopologue){52,  36,	1.10836E-02,    2.2408E+02,    2,     28.998270};
		m.ISO[2] = (Isotopologue){53,  28,	1.97822E-03,    1.1247E+02,    1,     29.999161};
		m.ISO[3] = (Isotopologue){54,  27,	3.67867E-04,    6.5934E+02,    6,     28.999130};
		m.ISO[4] = (Isotopologue){55,  38,	2.22250E-05,    2.3582E+02,    2,     31.002516};
		m.ISO[5] = (Isotopologue){56,  37,	4.13292E-06,    1.3809E+03,   12,     30.002485};

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 8465;

		sprintf(m.dataFilename[0], "%s%s", param.path, "05_hit12.");

		if(param.useHITEMP == 1){
			m.nFiles = 1;			//number of data files
			m.NL[ 0] = 113631;		//number of lines per data file

			m.NLmax = 113631;              //The naximum of number of lines per file

			m.fileLimit[ 0] = 0;
			m.fileLimit[ 1] = 8465;

			sprintf(m.dataFilename[0], "%s%s", param.path, "05_HITEMP2010.");
		}
	}
	if(m.id == 6){//CH4
		m.nFiles = 1;		//number of data files
		m.NL[0] = 468013;
		m.NLmax = 468013;
		m.nISO = 4;
		if(param.useHITEMP < 2){
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
			m.ISO[0] = (Isotopologue){61,  211,	.988274E+00,    5.9052E+02,    1,     16.031300};
			m.ISO[1] = (Isotopologue){62,  311,	1.11031E-02,    1.1808E+03,    2,     17.034655};
			m.ISO[2] = (Isotopologue){63,  212,	6.15751E-04,    4.7954E+03,    3,     17.037475};
			m.ISO[3] = (Isotopologue){64,  312,	6.91785E-06,    9.5990E+03,    6,     18.040830};

			m.fileLimit[ 0] = 0;
			m.fileLimit[ 1] = 11510;

			sprintf(m.dataFilename[0], "%s%s", param.path, "06_hit12.");
		}

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			char name[] = "12C-1H4__YT10to10";
			m.nStates = 8194057;
			m.defaultL = 0.0488;
			m.defaultn = 0.4;
			m.nFiles = 121;
			m.ntcol = 3;
			m.NL[0] = 7312353;
			m.NL[1] = 7417002;
			m.NL[2] = 7376496;
			m.NL[3] = 7142235;
			m.NL[4] = 6958489;
			m.NL[5] = 6826055;
			m.NL[6] = 6927182;
			m.NL[7] = 7355982;
			m.NL[8] = 8099195;
			m.NL[9] = 9170704;
			m.NL[10] = 10534275;
			m.NL[11] = 11965717;
			m.NL[12] = 13427471;
			m.NL[13] = 14635993;
			m.NL[14] = 15452960;
			m.NL[15] = 16000801;
			m.NL[16] = 16074877;
			m.NL[17] = 15895827;
			m.NL[18] = 15539903;
			m.NL[19] = 15342439;
			m.NL[20] = 15417227;
			m.NL[21] = 15995449;
			m.NL[22] = 17154818;
			m.NL[23] = 18909097;
			m.NL[24] = 21143787;
			m.NL[25] = 23723810;
			m.NL[26] = 26258491;
			m.NL[27] = 28596066;
			m.NL[28] = 30473271;
			m.NL[29] = 31670520;
			m.NL[30] = 32123630;
			m.NL[31] = 31962252;
			m.NL[32] = 31518336;
			m.NL[33] = 30888318;
			m.NL[34] = 30641730;
			m.NL[35] = 31016121;
			m.NL[36] = 32241309;
			m.NL[37] = 34523575;
			m.NL[38] = 37683972;
			m.NL[39] = 41542421;
			m.NL[40] = 45836433;
			m.NL[41] = 50010056;
			m.NL[42] = 53812186;
			m.NL[43] = 56771278;
			m.NL[44] = 58675019;
			m.NL[45] = 59526692;
			m.NL[46] = 59382179;
			m.NL[47] = 58748933;
			m.NL[48] = 58016946;
			m.NL[49] = 57846624;
			m.NL[50] = 58766039;
			m.NL[51] = 61141061;
			m.NL[52] = 65035751;
			m.NL[53] = 70477647;
			m.NL[54] = 77030152;
			m.NL[55] = 84264072;
			m.NL[56] = 91513066;
			m.NL[57] = 98038320;
			m.NL[58] = 103322250;
			m.NL[59] = 106982708;
			m.NL[60] = 108811943;
			m.NL[61] = 109184115;
			m.NL[62] = 108619315;
			m.NL[63] = 107960176;
			m.NL[64] = 108143731;
			m.NL[65] = 109804609;
			m.NL[66] = 113608962;
			m.NL[67] = 119753622;
			m.NL[68] = 128004234;
			m.NL[69] = 137999858;
			m.NL[70] = 148907218;
			m.NL[71] = 159745634;
			m.NL[72] = 169484099;
			m.NL[73] = 177202988;
			m.NL[74] = 182227325;
			m.NL[75] = 184374532;
			m.NL[76] = 183731473;
			m.NL[77] = 181025138;
			m.NL[78] = 177116059;
			m.NL[79] = 173078686;
			m.NL[80] = 169903575;
			m.NL[81] = 168394558;
			m.NL[82] = 168861373;
			m.NL[83] = 171512176;
			m.NL[84] = 175702628;
			m.NL[85] = 180878242;
			m.NL[86] = 186026797;
			m.NL[87] = 190151054;
			m.NL[88] = 192273315;
			m.NL[89] = 191886168;
			m.NL[90] = 188754046;
			m.NL[91] = 183142759;
			m.NL[92] = 175565071;
			m.NL[93] = 167048551;
			m.NL[94] = 158457440;
			m.NL[95] = 150713640;
			m.NL[96] = 144343634;
			m.NL[97] = 139647279;
			m.NL[98] = 136660014;
			m.NL[99] = 135119314;
			m.NL[100] = 129857500;
			m.NL[101] = 124001960;
			m.NL[102] = 118694873;
			m.NL[103] = 112466840;
			m.NL[104] = 105812803;
			m.NL[105] = 98096160;
			m.NL[106] = 89250592;
			m.NL[107] = 80746261;
			m.NL[108] = 71850011;
			m.NL[109] = 63609263;
			m.NL[110] = 56595055;
			m.NL[111] = 50561731;
			m.NL[112] = 45803978;
			m.NL[113] = 42396202;
			m.NL[114] = 39569085;
			m.NL[115] = 37584924;
			m.NL[116] = 35758387;
			m.NL[117] = 33822634;
			m.NL[118] = 31956905;
			m.NL[119] = 29207069;
			m.NL[120] = 11;

			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 100;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){61,  211,     1.0,    0.0,    0,     16.031300};
		
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			sprintf(qFilename, "%s%s%s", param.path, name, ".pf");
		}
	}
	if(m.id == 11){//NH3
		m.nFiles = 1;		//number of data files
		m.NL[0] = 46392;
		m.NLmax = 46392;
		m.nISO = 2;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){111,  4111,  .995872E+00,    1.7252E+03,    3,     17.026549};
		m.ISO[1] = (Isotopologue){112,  5111,  3.66129E-03,    1.1527E+03,    2,     18.023583};

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 7000;

		sprintf(m.dataFilename[0], "%s%s", param.path, "11_hit12.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}

		if(param.useHITEMP == 2){
			char name[] = "14N-1H3__BYTe";
			m.nStates = 4167360;
			m.nFiles = 120;
			m.ntcol = 3;
			m.defaultL = 0.053;
			m.defaultn = 0.5;
			m.NL[0] = 982907;
			m.NL[1] = 1037545;
			m.NL[2] = 1062896;
			m.NL[3] = 1082795;
			m.NL[4] = 1113869;
			m.NL[5] = 1178862;
			m.NL[6] = 1259075;
			m.NL[7] = 1330374;
			m.NL[8] = 1392503;
			m.NL[9] = 1448135;
			m.NL[10] = 1473957;
			m.NL[11] = 1525679;
			m.NL[12] = 1609146;
			m.NL[13] = 1721833;
			m.NL[14] = 1844544;
			m.NL[15] = 1967673;
			m.NL[16] = 2058265;
			m.NL[17] = 2102558;
			m.NL[18] = 2129362;
			m.NL[19] = 2144368;
			m.NL[20] = 2186309;
			m.NL[21] = 2271167;
			m.NL[22] = 2373033;
			m.NL[23] = 2501654;
			m.NL[24] = 2622434;
			m.NL[25] = 2723233;
			m.NL[26] = 2802387;
			m.NL[27] = 2898849;
			m.NL[28] = 3023510;
			m.NL[29] = 3195291;
			m.NL[30] = 3405551;
			m.NL[31] = 3617557;
			m.NL[32] = 3785845;
			m.NL[33] = 3913479;
			m.NL[34] = 3991611;
			m.NL[35] = 4039425;
			m.NL[36] = 4118227;
			m.NL[37] = 4238610;
			m.NL[38] = 4397129;
			m.NL[39] = 4594598;
			m.NL[40] = 4801158;
			m.NL[41] = 4972457;
			m.NL[42] = 5120515;
			m.NL[43] = 5289932;
			m.NL[44] = 5467013;
			m.NL[45] = 5693057;
			m.NL[46] = 5974425;
			m.NL[47] = 6256099;
			m.NL[48] = 6503065;
			m.NL[49] = 6706294;
			m.NL[50] = 6860855;
			m.NL[51] = 6960417;
			m.NL[52] = 7079819;
			m.NL[53] = 7234325;
			m.NL[54] = 7442462;
			m.NL[55] = 7711422;
			m.NL[56] = 8010525;
			m.NL[57] = 8276442;
			m.NL[58] = 8531943;
			m.NL[59] = 8789752;
			m.NL[60] = 9055326;
			m.NL[61] = 9368728;
			m.NL[62] = 9745838;
			m.NL[63] = 10130895;
			m.NL[64] = 10481180;
			m.NL[65] = 10780525;
			m.NL[66] = 11006291;
			m.NL[67] = 11163868;
			m.NL[68] = 11333508;
			m.NL[69] = 11540059;
			m.NL[70] = 11809529;
			m.NL[71] = 12166183;
			m.NL[72] = 12557951;
			m.NL[73] = 12937073;
			m.NL[74] = 13293545;
			m.NL[75] = 13639572;
			m.NL[76] = 13986280;
			m.NL[77] = 14379711;
			m.NL[78] = 14829655;
			m.NL[79] = 15288793;
			m.NL[80] = 15725255;
			m.NL[81] = 16113244;
			m.NL[82] = 16426095;
			m.NL[83] = 16667848;
			m.NL[84] = 16904219;
			m.NL[85] = 17157959;
			m.NL[86] = 17476169;
			m.NL[87] = 17898755;
			m.NL[88] = 18383515;
			m.NL[89] = 18882490;
			m.NL[90] = 19373298;
			m.NL[91] = 19837694;
			m.NL[92] = 20239006;
			m.NL[93] = 20646768;
			m.NL[94] = 21096539;
			m.NL[95] = 21547715;
			m.NL[96] = 21982458;
			m.NL[97] = 22373923;
			m.NL[98] = 22702270;
			m.NL[99] = 22934619;
			m.NL[100] = 22339283;
			m.NL[101] = 21108601;
			m.NL[102] = 19952989;
			m.NL[103] = 18926591;
			m.NL[104] = 18011268;
			m.NL[105] = 17192700;
			m.NL[106] = 16452057;
			m.NL[107] = 15692461;
			m.NL[108] = 14961423;
			m.NL[109] = 14165788;
			m.NL[110] = 13409218;
			m.NL[111] = 12705679;
			m.NL[112] = 12056386;
			m.NL[113] = 11407255;
			m.NL[114] = 10751832;
			m.NL[115] = 10124085;
			m.NL[116] = 9474612;
			m.NL[117] = 8854492;
			m.NL[118] = 8273327;
			m.NL[119] = 7750733;

			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 100;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename, "%s%s%s", param.path, name, ".pf");
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){111,  4111,  1.0,    0.0,    0,     17.026549};
		}

	}
	if(m.id == 23){//HCN
		m.nFiles = 1;		//number of data files
		m.NL[0] = 4253;
		m.NLmax = 4253;
		m.nISO = 3;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){231,  124,  .985114E+00,    8.9529E+02,    6,     27.010899};
		m.ISO[1] = (Isotopologue){232,  134,  1.10676E-02,    1.8403E+03,   12,     28.014254};
		m.ISO[2] = (Isotopologue){233,  125,  3.62174E-03,    6.2141E+02,    4,     28.007933};

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 3424;

		sprintf(m.dataFilename[0], "%s%s", param.path, "23_hit08.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			char name[] = "1H-12C-14N__Harris";
			m.defaultL = 0.084;
			m.defaultn = 0.5;
			m.nStates = 168110;
			m.nFiles = 1;
			m.ntcol = 4;
			m.NL[0] = 34418408;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 17586;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename, "%s%s%s", param.path, name, ".pf");
			sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
			m.ISO[0] = (Isotopologue){231,  124,  1.0,    0.0,    0,     27.010899};
		}
	}
	if(m.id == 26){//C2H2
		m.nFiles = 1;		//number of data files
		m.NL[0] = 20410;
		m.NLmax = 20410;
		m.nISO = 3;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){261,  1221,  .977599E+00,    4.1403E+02,    1,     26.015650};
		m.ISO[1] = (Isotopologue){262,  1231,  2.19663E-02,    1.6562E+03,    8,     27.019005};
		m.ISO[2] = (Isotopologue){263,  1222,  3.04550E-04,    1.5818E+03,    6,     27.021825};

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 9890;
		sprintf(m.dataFilename[0], "%s%s", param.path, "26_hit12.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
	}
	if(m.id == 31){//H2S
		m.nFiles = 1;		//number of data files
		m.NL[0] = 54235;
		m.NLmax = 54235;
		m.nISO = 3;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){311,  121,  .949884E+00,    5.0307E+02,    1,     33.987721};
		m.ISO[1] = (Isotopologue){312,  141,  4.21369E-02,    5.0435E+02,    1,     35.983515};
		m.ISO[1] = (Isotopologue){313,  131,  7.49766E-03,    2.0149E+03,    4,     34.987105};


		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 11330;
		sprintf(m.dataFilename[0], "%s%s", param.path, "31_hit12.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			char name[] = "1H2-32S__AYT2";
			m.nStates = 220618;
			m.defaultL = 0.07;
			m.defaultn = 0.5;

			m.nFiles = 35;
			m.ntcol = 3;
			m.NL[0] = 7226596;
			m.NL[1] = 7089219;
			m.NL[2] = 7227349;
			m.NL[3] = 7407306;
			m.NL[4] = 7808607;
			m.NL[5] = 8393450;
			m.NL[6] = 9177141;
			m.NL[7] = 10224670;
			m.NL[8] = 11359454;
			m.NL[9] = 11926204;
			m.NL[10] = 9383613;
			m.NL[11] = 6651647;
			m.NL[12] = 4578893;
			m.NL[13] = 3007385;
			m.NL[14] = 1891135;
			m.NL[15] = 1107584;
			m.NL[16] = 590239;
			m.NL[17] = 295196;
			m.NL[18] = 132366;
			m.NL[19] = 52319;
			m.NL[20] = 25093;
			m.NL[21] = 16301;
			m.NL[22] = 12454;
			m.NL[23] = 9865;
			m.NL[24] = 7923;
			m.NL[25] = 6177;
			m.NL[26] = 4697;
			m.NL[27] = 3448;
			m.NL[28] = 2525;
			m.NL[29] = 1807;
			m.NL[30] = 1192;
			m.NL[31] = 719;
			m.NL[32] = 352;
			m.NL[33] = 185;
			m.NL[34] = 69;

			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 1000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}

			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){311,  121,  1.0,    0.0,    0,     33.987721};
		
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			sprintf(qFilename, "%s%s%s", param.path, name, ".pf");
		}
	}
	if(m.id == 80){//VO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};


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
			char name[] = "51V-16O__VOMYT";
			m.defaultL = 0.07;
			m.defaultn = 0.5;
			m.nStates = 638958;
			m.nFiles = 7;
			m.ntcol = 4;
			m.NL[0] = 22585666;
			m.NL[1] = 31085102;
			m.NL[2] = 37095100;
			m.NL[3] = 42355006;
			m.NL[4] = 49301568;
			m.NL[5] = 56596408;
			m.NL[6] = 38112774;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 5000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename, "%s%s%s", param.path, name, ".pf");
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			//insert ISO information here
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};
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

