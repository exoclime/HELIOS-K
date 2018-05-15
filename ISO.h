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
void Init(Molecule &m, Param param, char (*qFilename)[160]){
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

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		sprintf(qFilename[3], "%s%s", param.path, "q.dat");
		sprintf(qFilename[4], "%s%s", param.path, "q.dat");
		sprintf(qFilename[5], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

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
			sprintf(m.mName, "%s", "1H2-16O__BT2");
			m.nStates = 221097;
			m.nFiles = 16;
			m.ntcol = 3;
			m.npfcol = 3;
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
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");

			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){11,  161,	1.0,    0.0,    0,     18.010565};
		}
		if(param.useHITEMP == 3){
			char name[] = "1H2-16O__POKAZATEL";
			sprintf(m.mName, "%s", "1H2-16O__POKAZATEL");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nStates = 810269;
			m.nFiles = 412;
			m.ntcol = 3;
			m.npfcol = 2;
			m.NL[0] = 63545284;
			m.NL[1] = 62953585;
			m.NL[2] = 62406329;
			m.NL[3] = 61806710;
			m.NL[4] = 61211769;
			m.NL[5] = 60645422;
			m.NL[6] = 60081707;
			m.NL[7] = 59527108;
			m.NL[8] = 58997522;
			m.NL[9] = 58454032;
			m.NL[10] = 57898326;
			m.NL[11] = 57354828;
			m.NL[12] = 56822761;
			m.NL[13] = 56299929;
			m.NL[14] = 55785749;
			m.NL[15] = 55257484;
			m.NL[16] = 54737905;
			m.NL[17] = 54234508;
			m.NL[18] = 53719800;
			m.NL[19] = 53209511;
			m.NL[20] = 52710632;
			m.NL[21] = 52216798;
			m.NL[22] = 51728862;
			m.NL[23] = 51242841;
			m.NL[24] = 50767793;
			m.NL[25] = 50283335;
			m.NL[26] = 49806042;
			m.NL[27] = 49334764;
			m.NL[28] = 48870521;
			m.NL[29] = 48408485;
			m.NL[30] = 47944583;
			m.NL[31] = 47483363;
			m.NL[32] = 47032973;
			m.NL[33] = 46586524;
			m.NL[34] = 46139607;
			m.NL[35] = 45686142;
			m.NL[36] = 45247489;
			m.NL[37] = 44819995;
			m.NL[38] = 44366512;
			m.NL[39] = 43949647;
			m.NL[40] = 43525354;
			m.NL[41] = 43092528;
			m.NL[42] = 42680142;
			m.NL[43] = 42267798;
			m.NL[44] = 41843019;
			m.NL[45] = 41430274;
			m.NL[46] = 41035539;
			m.NL[47] = 40625567;
			m.NL[48] = 40219144;
			m.NL[49] = 39832882;
			m.NL[50] = 39433757;
			m.NL[51] = 39038150;
			m.NL[52] = 38659097;
			m.NL[53] = 38276595;
			m.NL[54] = 37896041;
			m.NL[55] = 37514544;
			m.NL[56] = 37127988;
			m.NL[57] = 36765424;
			m.NL[58] = 36392834;
			m.NL[59] = 36024995;
			m.NL[60] = 35663751;
			m.NL[61] = 35298617;
			m.NL[62] = 34944150;
			m.NL[63] = 34593236;
			m.NL[64] = 34240714;
			m.NL[65] = 33880586;
			m.NL[66] = 33525576;
			m.NL[67] = 33194493;
			m.NL[68] = 32847204;
			m.NL[69] = 32512662;
			m.NL[70] = 32174920;
			m.NL[71] = 31837631;
			m.NL[72] = 31526218;
			m.NL[73] = 31173580;
			m.NL[74] = 30848655;
			m.NL[75] = 30516886;
			m.NL[76] = 30207120;
			m.NL[77] = 29886457;
			m.NL[78] = 29576533;
			m.NL[79] = 29256380;
			m.NL[80] = 28941993;
			m.NL[81] = 28648968;
			m.NL[82] = 28344141;
			m.NL[83] = 28034160;
			m.NL[84] = 27734214;
			m.NL[85] = 27438885;
			m.NL[86] = 27144461;
			m.NL[87] = 26852091;
			m.NL[88] = 26562807;
			m.NL[89] = 26268028;
			m.NL[90] = 25992322;
			m.NL[91] = 25705538;
			m.NL[92] = 25423068;
			m.NL[93] = 25137620;
			m.NL[94] = 24867033;
			m.NL[95] = 24594433;
			m.NL[96] = 24311872;
			m.NL[97] = 24044102;
			m.NL[98] = 23786062;
			m.NL[99] = 23507410;
			m.NL[100] = 23263918;
			m.NL[101] = 23000316;
			m.NL[102] = 22735511;
			m.NL[103] = 22477652;
			m.NL[104] = 22232337;
			m.NL[105] = 21972714;
			m.NL[106] = 21717726;
			m.NL[107] = 21471849;
			m.NL[108] = 21222524;
			m.NL[109] = 20981716;
			m.NL[110] = 20734222;
			m.NL[111] = 20499995;
			m.NL[112] = 20271442;
			m.NL[113] = 20029427;
			m.NL[114] = 19804057;
			m.NL[115] = 19571228;
			m.NL[116] = 19335526;
			m.NL[117] = 19117904;
			m.NL[118] = 18886005;
			m.NL[119] = 18666584;
			m.NL[120] = 18441986;
			m.NL[121] = 18225164;
			m.NL[122] = 18015981;
			m.NL[123] = 17791131;
			m.NL[124] = 17574629;
			m.NL[125] = 17368280;
			m.NL[126] = 17158791;
			m.NL[127] = 16952362;
			m.NL[128] = 16743232;
			m.NL[129] = 16540839;
			m.NL[130] = 16334996;
			m.NL[131] = 16139076;
			m.NL[132] = 15940267;
			m.NL[133] = 15736472;
			m.NL[134] = 15550474;
			m.NL[135] = 15358235;
			m.NL[136] = 15162945;
			m.NL[137] = 14980680;
			m.NL[138] = 14785642;
			m.NL[139] = 14593689;
			m.NL[140] = 14411900;
			m.NL[141] = 14230937;
			m.NL[142] = 14049115;
			m.NL[143] = 13869788;
			m.NL[144] = 13693722;
			m.NL[145] = 13524016;
			m.NL[146] = 13346961;
			m.NL[147] = 13173937;
			m.NL[148] = 13001843;
			m.NL[149] = 12836207;
			m.NL[150] = 12677489;
			m.NL[151] = 12508092;
			m.NL[152] = 12339773;
			m.NL[153] = 12176175;
			m.NL[154] = 12011252;
			m.NL[155] = 11849977;
			m.NL[156] = 11689452;
			m.NL[157] = 11539565;
			m.NL[158] = 11386333;
			m.NL[159] = 11233555;
			m.NL[160] = 11076811;
			m.NL[161] = 10925692;
			m.NL[162] = 10783635;
			m.NL[163] = 10641509;
			m.NL[164] = 10486323;
			m.NL[165] = 10346100;
			m.NL[166] = 10204824;
			m.NL[167] = 10057861;
			m.NL[168] = 9921926;
			m.NL[169] = 9773592;
			m.NL[170] = 9633923;
			m.NL[171] = 9501177;
			m.NL[172] = 9372751;
			m.NL[173] = 9239376;
			m.NL[174] = 9108502;
			m.NL[175] = 8977426;
			m.NL[176] = 8850163;
			m.NL[177] = 8731626;
			m.NL[178] = 8600205;
			m.NL[179] = 8469901;
			m.NL[180] = 8354958;
			m.NL[181] = 8231521;
			m.NL[182] = 8109897;
			m.NL[183] = 7980537;
			m.NL[184] = 7866505;
			m.NL[185] = 7751710;
			m.NL[186] = 7633609;
			m.NL[187] = 7529115;
			m.NL[188] = 7404183;
			m.NL[189] = 7297966;
			m.NL[190] = 7189256;
			m.NL[191] = 7080320;
			m.NL[192] = 6974491;
			m.NL[193] = 6865014;
			m.NL[194] = 6758330;
			m.NL[195] = 6653612;
			m.NL[196] = 6549313;
			m.NL[197] = 6448230;
			m.NL[198] = 6341090;
			m.NL[199] = 6245922;
			m.NL[200] = 6144754;
			m.NL[201] = 6047808;
			m.NL[202] = 5950372;
			m.NL[203] = 5856896;
			m.NL[204] = 5763730;
			m.NL[205] = 5670113;
			m.NL[206] = 5583587;
			m.NL[207] = 5487976;
			m.NL[208] = 5399430;
			m.NL[209] = 5311892;
			m.NL[210] = 5223913;
			m.NL[211] = 5141011;
			m.NL[212] = 5054055;
			m.NL[213] = 4965864;
			m.NL[214] = 4883922;
			m.NL[215] = 4805252;
			m.NL[216] = 4720428;
			m.NL[217] = 4645826;
			m.NL[218] = 4562340;
			m.NL[219] = 4483201;
			m.NL[220] = 4408052;
			m.NL[221] = 4329692;
			m.NL[222] = 4255398;
			m.NL[223] = 4178932;
			m.NL[224] = 4108751;
			m.NL[225] = 4034007;
			m.NL[226] = 3961668;
			m.NL[227] = 3891961;
			m.NL[228] = 3823382;
			m.NL[229] = 3750682;
			m.NL[230] = 3682246;
			m.NL[231] = 3619996;
			m.NL[232] = 3555110;
			m.NL[233] = 3482747;
			m.NL[234] = 3423532;
			m.NL[235] = 3359484;
			m.NL[236] = 3293612;
			m.NL[237] = 3231720;
			m.NL[238] = 3174049;
			m.NL[239] = 3113984;
			m.NL[240] = 3051475;
			m.NL[241] = 2997391;
			m.NL[242] = 2941919;
			m.NL[243] = 2883283;
			m.NL[244] = 2829925;
			m.NL[245] = 2773571;
			m.NL[246] = 2722646;
			m.NL[247] = 2668326;
			m.NL[248] = 2615504;
			m.NL[249] = 2565382;
			m.NL[250] = 2512457;
			m.NL[251] = 2461276;
			m.NL[252] = 2411607;
			m.NL[253] = 2357275;
			m.NL[254] = 2312711;
			m.NL[255] = 2265136;
			m.NL[256] = 2215405;
			m.NL[257] = 2172406;
			m.NL[258] = 2126703;
			m.NL[259] = 2079734;
			m.NL[260] = 2037120;
			m.NL[261] = 1996408;
			m.NL[262] = 1951206;
			m.NL[263] = 1909330;
			m.NL[264] = 1870015;
			m.NL[265] = 1829489;
			m.NL[266] = 1787365;
			m.NL[267] = 1749418;
			m.NL[268] = 1706917;
			m.NL[269] = 1668210;
			m.NL[270] = 1634206;
			m.NL[271] = 1594081;
			m.NL[272] = 1561266;
			m.NL[273] = 1524338;
			m.NL[274] = 1486528;
			m.NL[275] = 1459049;
			m.NL[276] = 1426572;
			m.NL[277] = 1391368;
			m.NL[278] = 1360892;
			m.NL[279] = 1328212;
			m.NL[280] = 1299511;
			m.NL[281] = 1269161;
			m.NL[282] = 1237303;
			m.NL[283] = 1207013;
			m.NL[284] = 1176729;
			m.NL[285] = 1146250;
			m.NL[286] = 1115908;
			m.NL[287] = 1091822;
			m.NL[288] = 1061334;
			m.NL[289] = 1031940;
			m.NL[290] = 1008240;
			m.NL[291] = 984076;
			m.NL[292] = 958478;
			m.NL[293] = 932601;
			m.NL[294] = 910169;
			m.NL[295] = 885270;
			m.NL[296] = 862238;
			m.NL[297] = 840305;
			m.NL[298] = 816332;
			m.NL[299] = 795097;
			m.NL[300] = 771329;
			m.NL[301] = 750466;
			m.NL[302] = 730919;
			m.NL[303] = 710066;
			m.NL[304] = 686190;
			m.NL[305] = 670611;
			m.NL[306] = 651176;
			m.NL[307] = 632699;
			m.NL[308] = 614830;
			m.NL[309] = 596672;
			m.NL[310] = 582422;
			m.NL[311] = 564097;
			m.NL[312] = 550924;
			m.NL[313] = 534325;
			m.NL[314] = 517564;
			m.NL[315] = 502065;
			m.NL[316] = 487280;
			m.NL[317] = 473291;
			m.NL[318] = 458215;
			m.NL[319] = 442091;
			m.NL[320] = 425589;
			m.NL[321] = 413527;
			m.NL[322] = 399363;
			m.NL[323] = 385721;
			m.NL[324] = 372955;
			m.NL[325] = 360682;
			m.NL[326] = 349944;
			m.NL[327] = 337237;
			m.NL[328] = 327431;
			m.NL[329] = 316290;
			m.NL[330] = 305728;
			m.NL[331] = 293821;
			m.NL[332] = 282571;
			m.NL[333] = 274358;
			m.NL[334] = 263242;
			m.NL[335] = 252795;
			m.NL[336] = 243813;
			m.NL[337] = 234264;
			m.NL[338] = 226171;
			m.NL[339] = 217573;
			m.NL[340] = 208468;
			m.NL[341] = 201875;
			m.NL[342] = 194596;
			m.NL[343] = 186633;
			m.NL[344] = 179917;
			m.NL[345] = 173142;
			m.NL[346] = 167686;
			m.NL[347] = 160652;
			m.NL[348] = 154675;
			m.NL[349] = 148613;
			m.NL[350] = 142142;
			m.NL[351] = 136632;
			m.NL[352] = 130219;
			m.NL[353] = 125425;
			m.NL[354] = 119123;
			m.NL[355] = 113675;
			m.NL[356] = 107743;
			m.NL[357] = 102563;
			m.NL[358] = 98131;
			m.NL[359] = 91965;
			m.NL[360] = 88649;
			m.NL[361] = 84326;
			m.NL[362] = 80712;
			m.NL[363] = 75947;
			m.NL[364] = 72654;
			m.NL[365] = 69178;
			m.NL[366] = 65585;
			m.NL[367] = 61641;
			m.NL[368] = 58630;
			m.NL[369] = 54629;
			m.NL[370] = 52365;
			m.NL[371] = 48769;
			m.NL[372] = 45718;
			m.NL[373] = 43112;
			m.NL[374] = 40754;
			m.NL[375] = 37906;
			m.NL[376] = 36021;
			m.NL[377] = 34287;
			m.NL[378] = 32013;
			m.NL[379] = 30205;
			m.NL[380] = 28424;
			m.NL[381] = 27158;
			m.NL[382] = 25391;
			m.NL[383] = 24482;
			m.NL[384] = 23041;
			m.NL[385] = 21616;
			m.NL[386] = 20305;
			m.NL[387] = 18771;
			m.NL[388] = 17431;
			m.NL[389] = 15980;
			m.NL[390] = 14676;
			m.NL[391] = 13241;
			m.NL[392] = 12105;
			m.NL[393] = 10701;
			m.NL[394] = 9523;
			m.NL[395] = 8717;
			m.NL[396] = 7739;
			m.NL[397] = 7120;
			m.NL[398] = 6561;
			m.NL[399] = 5908;
			m.NL[400] = 5498;
			m.NL[401] = 4816;
			m.NL[402] = 4306;
			m.NL[403] = 3623;
			m.NL[404] = 2841;
			m.NL[405] = 2232;
			m.NL[406] = 1760;
			m.NL[407] = 1305;
			m.NL[408] = 845;
			m.NL[409] = 519;
			m.NL[410] = 187;
			m.NL[411] = 20;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 100;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){11,  161,	1.0,    0.0,    0,     18.010565};
			//version = 20180501
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

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		sprintf(qFilename[3], "%s%s", param.path, "q.dat");
		sprintf(qFilename[4], "%s%s", param.path, "q.dat");
		sprintf(qFilename[5], "%s%s", param.path, "q.dat");
		sprintf(qFilename[6], "%s%s", param.path, "q.dat");
		sprintf(qFilename[7], "%s%s", param.path, "q.dat");
		sprintf(qFilename[8], "%s%s", param.path, "q.dat");
		sprintf(qFilename[9], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

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
	
		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		sprintf(qFilename[3], "%s%s", param.path, "q.dat");
		sprintf(qFilename[4], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

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

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		sprintf(qFilename[3], "%s%s", param.path, "q.dat");
		sprintf(qFilename[4], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

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

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		sprintf(qFilename[3], "%s%s", param.path, "q.dat");
		sprintf(qFilename[4], "%s%s", param.path, "q.dat");
		sprintf(qFilename[5], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

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

		if(param.useHITEMP == 2){

			sprintf(m.mName, "%s", "12C-16O__Li2015");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nStates = 6383;
			m.nFiles = 1;
			m.ntcol = 3;
			m.npfcol = 2;
			m.NL[0] = 145636;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 22000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){51,  26,  1.0,    0.0,    0,      28.0101};
			//version = 20170101

		}
		if(param.useHITEMP == 3){
			//from Vald database
			sprintf(m.mName, "%s", "KuruczCO");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 161139;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 33332;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){51,  26,  1.0,    0.0,    0,      28.0101};

		}
		if(param.useHITEMP == 4){

			sprintf(m.mName, "%s", "KuruczCOax");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 396946;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 89819;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){51,  26,  1.0,    0.0,    0,      28.0101};

		}
	}
	if(m.id == 6){//CH4
		m.nFiles = 1;		//number of data files
		m.NL[0] = 468013;
		m.NLmax = 468013;
		m.nISO = 4;

		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){61,  211,	.988274E+00,    5.9052E+02,    1,     16.031300};
		m.ISO[1] = (Isotopologue){62,  311,	1.11031E-02,    1.1808E+03,    2,     17.034655};
		m.ISO[2] = (Isotopologue){63,  212,	6.15751E-04,    4.7954E+03,    3,     17.037475};
		m.ISO[3] = (Isotopologue){64,  312,	6.91785E-06,    9.5990E+03,    6,     18.040830};

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		sprintf(qFilename[3], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 11510;

		sprintf(m.dataFilename[0], "%s%s", param.path, "06_hit12.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			sprintf(m.mName, "%s", "12C-1H4__YT10to10");
			m.nStates = 8194057;
			m.defaultL = 0.0488;
			m.defaultn = 0.4;
			m.nFiles = 121;
			m.ntcol = 3;
			m.npfcol = 2;
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
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
		}
	}
	if(m.id == 7){//O2
		m.nFiles = 1;		//number of data files
		m.NL[0] = 14085;
		m.NLmax = 14085;
		m.nISO = 3;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){71,  66,  0.995262,    215.73,    1,     31.98983};
		m.ISO[1] = (Isotopologue){72,  68,  0.003991,    455.23,    1,     33.994076};
		m.ISO[2] = (Isotopologue){73,  67,  7.422350E-4, 2658.12,   6,     32.994045};

		sprintf(qFilename[0], "%s%s", param.path, "q36.txt");
		sprintf(qFilename[1], "%s%s", param.path, "q37.txt");
		sprintf(qFilename[2], "%s%s", param.path, "q38.txt");
		m.npfcol = 2;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 17273;

		sprintf(m.dataFilename[0], "%s%s", param.path, "07_hit16.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
		}
	}
	if(m.id == 8){//NO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 105079;
		m.NLmax = 105079;
		m.nISO = 3;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){81,  46,   	0.993974,    1142.13,   3,     29.997989};
		m.ISO[1] = (Isotopologue){82,  56,   	0.003654,    789.26,    2,     30.995023};
		m.ISO[2] = (Isotopologue){83,  48,   	0.001993,    1204.44,   3,     32.002234};

		sprintf(qFilename[0], "%s%s", param.path, "q39.txt");
		sprintf(qFilename[1], "%s%s", param.path, "q40.txt");
		sprintf(qFilename[2], "%s%s", param.path, "q41.txt");
		m.npfcol = 2;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 9274;

		sprintf(m.dataFilename[0], "%s%s", param.path, "08_hit16.");

		if(param.useHITEMP == 1){
			m.nFiles = 1;			//number of data files
			m.NL[ 0] = 115610;		//number of lines per data file
			m.NLmax =  115610;              //The naximum of number of lines per file

			m.fileLimit[ 0] = 0;
			m.fileLimit[ 1] = 9274;

			sprintf(m.dataFilename[0], "%s%s", param.path, "08_HITEMP2010.");
		}
		if(param.useHITEMP == 2){
			char name[] = "14N-16O__NOname";
			sprintf(m.mName, "%s", "14N-16O__NOname");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nStates = 21688;
			m.nFiles = 1 ;
			m.ntcol = 4;
			m.npfcol = 2;
			m.NL[0] = 2280366;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 40000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  46,  1.0,    0.0,    0,      29.997988};
			//version = 20170517
		}
	}
	if(m.id == 9){//SO2
		m.nFiles = 1;		//number of data files
		m.NL[0] = 95121;
		m.NLmax = 95121;
		m.nISO = 2;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){91,  626,  0.945678,    6340.30,    1,     63.961901};
		m.ISO[1] = (Isotopologue){92,  646,  0.041950,    6368.98,    1,     65.957695};
		sprintf(qFilename[0], "%s%s", param.path, "q42.txt");
		sprintf(qFilename[1], "%s%s", param.path, "q43.txt");
		m.npfcol = 2;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 4093;
		sprintf(m.dataFilename[0], "%s%s", param.path, "09_hit16.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			char name[] = "32S-16O2__ExoAmes";
			sprintf(m.mName, "%s", "32S-16O2__ExoAmes");
			m.defaultL = 0.1063 ;
			m.defaultn = 0.695 ;
			m.nStates = 3270270;
			m.nFiles = 80 ;
			m.ntcol = 3;
			m.npfcol = 2;
			m.NL[0] = 10844718;
			m.NL[1] = 12016717;
			m.NL[2] = 12071854;
			m.NL[3] = 13058966;
			m.NL[4] = 14836227;
			m.NL[5] = 14694847;
			m.NL[6] = 14431756;
			m.NL[7] = 13953275;
			m.NL[8] = 14300070;
			m.NL[9] = 16133463;
			m.NL[10] = 17720307;
			m.NL[11] = 17723715;
			m.NL[12] = 19042028;
			m.NL[13] = 19557779;
			m.NL[14] = 18543254;
			m.NL[15] = 18488857;
			m.NL[16] = 17322481;
			m.NL[17] = 17324579;
			m.NL[18] = 17756889;
			m.NL[19] = 17962764;
			m.NL[20] = 19513981;
			m.NL[21] = 21030137;
			m.NL[22] = 20426471;
			m.NL[23] = 22135743;
			m.NL[24] = 22929268;
			m.NL[25] = 21532143;
			m.NL[26] = 21966693;
			m.NL[27] = 21443358;
			m.NL[28] = 20316907;
			m.NL[29] = 20563461;
			m.NL[30] = 20268467;
			m.NL[31] = 20231123;
			m.NL[32] = 21375929;
			m.NL[33] = 20963690;
			m.NL[34] = 22729805;
			m.NL[35] = 23219958;
			m.NL[36] = 21840928;
			m.NL[37] = 22165918;
			m.NL[38] = 22522086;
			m.NL[39] = 21787138;
			m.NL[40] = 21683169;
			m.NL[41] = 21234543;
			m.NL[42] = 20784096;
			m.NL[43] = 20524452;
			m.NL[44] = 19777914;
			m.NL[45] = 21380948;
			m.NL[46] = 21527563;
			m.NL[47] = 20782894;
			m.NL[48] = 21488130;
			m.NL[49] = 22202666;
			m.NL[50] = 20926756;
			m.NL[51] = 20361370;
			m.NL[52] = 20005221;
			m.NL[53] = 19426510;
			m.NL[54] = 18503939;
			m.NL[55] = 18166461;
			m.NL[56] = 18330680;
			m.NL[57] = 17948571;
			m.NL[58] = 18029424;
			m.NL[59] = 18834413;
			m.NL[60] = 18504440;
			m.NL[61] = 17381887;
			m.NL[62] = 16550617;
			m.NL[63] = 16032378;
			m.NL[64] = 15602179;
			m.NL[65] = 15076992;
			m.NL[66] = 14823808;
			m.NL[67] = 14295603;
			m.NL[68] = 13638078;
			m.NL[69] = 13368442;
			m.NL[70] = 12927122;
			m.NL[71] = 11809307;
			m.NL[72] = 11048886;
			m.NL[73] = 11048556;
			m.NL[74] = 10619552;
			m.NL[75] = 9584555;
			m.NL[76] = 8828811;
			m.NL[77] = 8319879;
			m.NL[78] = 7520527;
			m.NL[79] = 6608600;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 100;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){91,  626,  1.0,    0.0,    0,      63.961900};
			//version = 20170131

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

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 7000;

		sprintf(m.dataFilename[0], "%s%s", param.path, "11_hit12.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}

		if(param.useHITEMP == 2){
			sprintf(m.mName, "%s", "14N-1H3__BYTe");
			m.nStates = 4167360;
			m.nFiles = 120;
			m.ntcol = 3;
			m.npfcol = 2;
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
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){111,  4111,  1.0,    0.0,    0,     17.026549};
		}

	}
	if(m.id == 13){//OH
		m.nFiles = 1;		//number of data files
		m.NL[0] = 33058;
		m.NLmax = 33058;
		m.nISO = 3;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){131,  61,  0.997473,    80.35,    2,     17.00274};
		m.ISO[1] = (Isotopologue){132,  81,  0.002000,    80.88,    2,     19.006986};
		m.ISO[2] = (Isotopologue){133,  62,  1.553710E-4, 209.32,   3,     18.008915};

		sprintf(qFilename[0], "%s%s", param.path, "q48.txt");
		sprintf(qFilename[1], "%s%s", param.path, "q49.txt");
		sprintf(qFilename[2], "%s%s", param.path, "q50.txt");
		m.npfcol = 2;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 35875;

		sprintf(m.dataFilename[0], "%s%s", param.path, "13_hit16.");

		if(param.useHITEMP == 1){
			m.nFiles = 1;			//number of data files
			m.NL[ 0] = 41557;		//number of lines per data file
			m.NLmax =  41557;              //The naximum of number of lines per file

			m.fileLimit[ 0] = 0;
			m.fileLimit[ 1] = 19268;

			sprintf(m.dataFilename[0], "%s%s", param.path, "13_HITEMP2010.");
		}
		if(param.useHITEMP == 2){
			printf("Error: no EXOMOL data for this molecule\n");
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

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 3424;

		sprintf(m.dataFilename[0], "%s%s", param.path, "23_hit08.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			sprintf(m.mName, "%s", "1H-12C-14N__Harris");
			m.defaultL = 0.084;
			m.defaultn = 0.5;
			m.nStates = 168110;
			m.nFiles = 1;
			m.ntcol = 4;
			m.npfcol = 2;
			m.NL[0] = 34418408;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 17586;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
			sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
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

		sprintf(qFilename[0], "%s%s", param.path, "q.dat");
		sprintf(qFilename[1], "%s%s", param.path, "q.dat");
		sprintf(qFilename[2], "%s%s", param.path, "q.dat");
		m.npfcol = 0;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 9890;
		sprintf(m.dataFilename[0], "%s%s", param.path, "26_hit12.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
	}
	if(m.id == 28){//PH3
		m.nFiles = 1;		//number of data files
		m.NL[0] = 22189;
		m.NLmax = 22189;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){281,  1111,  0.999533,    3249.44,    2,     33.997238};

		sprintf(qFilename[0], "%s%s", param.path, "q79.txt");
		m.npfcol = 2;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 3602;

		sprintf(m.dataFilename[0], "%s%s", param.path, "28_hit16.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			char name[] = "31P-1H3__SAlTY";
			sprintf(m.mName, "%s", "31P-1H3__SAlTY");
			m.defaultL = 0.0750 ;
			m.defaultn = 0.530 ;
			m.nStates = 9787832;
			m.nFiles = 100 ;
			m.ntcol = 3;
			m.npfcol = 2;
			m.NL[0] = 18996939;
			m.NL[1] = 18543301;
			m.NL[2] = 17331258;
			m.NL[3] = 15855561;
			m.NL[4] = 15178350;
			m.NL[5] = 15189419;
			m.NL[6] = 16618207;
			m.NL[7] = 19049749;
			m.NL[8] = 22612913;
			m.NL[9] = 26696382;
			m.NL[10] = 29957530;
			m.NL[11] = 31762865;
			m.NL[12] = 31076117;
			m.NL[13] = 29411843;
			m.NL[14] = 27259513;
			m.NL[15] = 26294909;
			m.NL[16] = 26477699;
			m.NL[17] = 28542757;
			m.NL[18] = 32350239;
			m.NL[19] = 37466098;
			m.NL[20] = 43403119;
			m.NL[21] = 47914839;
			m.NL[22] = 50475617;
			m.NL[23] = 49698193;
			m.NL[24] = 47448367;
			m.NL[25] = 44610078;
			m.NL[26] = 43271540;
			m.NL[27] = 43891634;
			m.NL[28] = 47015336;
			m.NL[29] = 52690941;
			m.NL[30] = 59803089;
			m.NL[31] = 67755745;
			m.NL[32] = 73764049;
			m.NL[33] = 77010826;
			m.NL[34] = 76063255;
			m.NL[35] = 73022174;
			m.NL[36] = 69610036;
			m.NL[37] = 68179081;
			m.NL[38] = 69753487;
			m.NL[39] = 74480890;
			m.NL[40] = 82458905;
			m.NL[41] = 92039102;
			m.NL[42] = 102271933;
			m.NL[43] = 109837013;
			m.NL[44] = 113535295;
			m.NL[45] = 112347913;
			m.NL[46] = 108512167;
			m.NL[47] = 104885818;
			m.NL[48] = 103740267;
			m.NL[49] = 106780337;
			m.NL[50] = 113721492;
			m.NL[51] = 124518419;
			m.NL[52] = 136994069;
			m.NL[53] = 149579396;
			m.NL[54] = 158709160;
			m.NL[55] = 162448161;
			m.NL[56] = 161000720;
			m.NL[57] = 156346617;
			m.NL[58] = 152896864;
			m.NL[59] = 152520061;
			m.NL[60] = 157494207;
			m.NL[61] = 167162708;
			m.NL[62] = 181208067;
			m.NL[63] = 196998289;
			m.NL[64] = 211718242;
			m.NL[65] = 222271512;
			m.NL[66] = 225683107;
			m.NL[67] = 223816943;
			m.NL[68] = 218488235;
			m.NL[69] = 215198980;
			m.NL[70] = 215981874;
			m.NL[71] = 223437189;
			m.NL[72] = 236772204;
			m.NL[73] = 254689972;
			m.NL[74] = 274264076;
			m.NL[75] = 291435430;
			m.NL[76] = 303245945;
			m.NL[77] = 306404436;
			m.NL[78] = 303349306;
			m.NL[79] = 297385325;
			m.NL[80] = 294169715;
			m.NL[81] = 297058048;
			m.NL[82] = 307750262;
			m.NL[83] = 325806201;
			m.NL[84] = 348611531;
			m.NL[85] = 373049061;
			m.NL[86] = 393807899;
			m.NL[87] = 406712172;
			m.NL[88] = 409482824;
			m.NL[89] = 404253838;
			m.NL[90] = 397079438;
			m.NL[91] = 393091618;
			m.NL[92] = 397441390;
			m.NL[93] = 410615106;
			m.NL[94] = 432522811;
			m.NL[95] = 459569517;
			m.NL[96] = 487363916;
			m.NL[97] = 509550607;
			m.NL[98] = 519490420;
			m.NL[99] = 520409264;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 100;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, name, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){281,  1111,  1.0,    0.0,    0,      33.997237};
			//version = 20170131

		}
	}
	if(m.id == 31){//H2S
		m.nFiles = 1;		//number of data files
		m.NL[0] = 54228;
		m.NLmax = 54228;
		m.nISO = 3;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){311,  121,  .949884E+00,    505.79,    1,     33.987721};
		m.ISO[1] = (Isotopologue){312,  141,  4.21369E-02,    504.35,    1,     35.983515};
		m.ISO[2] = (Isotopologue){313,  131,  7.49766E-03,    2014.94,    4,     34.987105};

		sprintf(qFilename[0], "%s%s", param.path, "q81.txt");
		sprintf(qFilename[1], "%s%s", param.path, "q82.txt");
		sprintf(qFilename[2], "%s%s", param.path, "q83.txt");
		m.npfcol = 2;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 11330;
		sprintf(m.dataFilename[0], "%s%s", param.path, "31_hit16.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
			sprintf(m.mName, "%s", "1H2-32S__AYT2");
			m.nStates = 220618;
			m.defaultL = 0.07;
			m.defaultn = 0.5;

			m.nFiles = 35;
			m.ntcol = 3;
			m.npfcol = 2;
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
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
		}
	}
	if(m.id == 47){//SO3
		m.nFiles = 1;		//number of data files
		m.NL[0] = 14295;
		m.NLmax = 14295;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){471,  26,  0.943400,    7783.30,    1,     79.95682};

		sprintf(qFilename[0], "%s%s", param.path, "q114.txt");
		m.npfcol = 2;

		m.fileLimit[ 0] = 0;
		m.fileLimit[ 1] = 2825;
		sprintf(m.dataFilename[0], "%s%s", param.path, "47_hit16.");

		if(param.useHITEMP == 1){
			printf("Error: no HITEMP data for this molecule\n");
		}
		if(param.useHITEMP == 2){
		}
	}
	if(m.id == 80){//VO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){801,  5116,  1.0,    0.0,    0,     66.938871};
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
			sprintf(m.mName, "%s", "51V-16O__VOMYT");
			m.defaultL = 0.07;
			m.defaultn = 0.5;
			m.nStates = 638958;
			m.nFiles = 7;
			m.ntcol = 4;
			m.npfcol = 2;
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
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
			for(int i = 0; i < m.nFiles; ++i){
				sprintf(m.dataFilename[i], "%s%s__%05d-%05d.", param.path, m.mName, m.fileLimit[i], m.fileLimit[i + 1]);
			}
			//insert ISO information here
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};
		}
		if(param.useHITEMP == 3){
			sprintf(m.mName, "%s", "KuruczVO");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 4509519;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 28096;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};

		}
		if(param.useHITEMP == 4){
			sprintf(m.mName, "%s", "KuruczVOax");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1611872;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 17182;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};

		}
		if(param.useHITEMP == 5){
			sprintf(m.mName, "%s", "KuruczVObx");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1841400;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 21473;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};

		}
		if(param.useHITEMP == 6){
			sprintf(m.mName, "%s", "KuruczVOcx");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 1056247;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 28097;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};

		}
		if(param.useHITEMP == 7){
			sprintf(m.mName, "%s", "KuruczVOmyt");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 49274820;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 35011;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){81,  16,  1.0,    0.0,    0,     66.938871};

		}
	}
	if(m.id == 81){//TI0
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){811,  4816,  1.0,    0.0,    0,     63.942862};
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
			printf("Error: no Exomol data for this molecule\n");
		}
		if(param.useHITEMP == 3){

			sprintf(m.mName, "%s", "Plez2012");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8262872 ;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 31587.9 ;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}

			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){811,  4816,  1.0,    0.0,    0,     63.942862};
		}
		if(param.useHITEMP == 4){

			sprintf(m.mName, "%s", "Plez2012-norlander");
			m.defaultL = 0.0700;                                                                                                                                                                       
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8259589 ;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 32467.8 ;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}

			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){811,  4816,  1.0,    0.0,    0,     63.942862};
		}
		if(param.useHITEMP == 5){

			sprintf(m.mName, "%s", "Plez2012-philips");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 39093 ;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 22610.1 ;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){811,  4816,  1.0,    0.0,    0,     63.942862};
		}
		if(param.useHITEMP == 6){

			sprintf(m.mName, "%s", "Plez2012-polfits");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8339450 ;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 31587.9 ;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}

			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){811,  4816,  1.0,    0.0,    0,     63.942862};
		}
		if(param.useHITEMP == 7){

			sprintf(m.mName, "%s", "Schwenke1998");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nFiles = 1;
			m.npfcol = 2;
			m.NL[0] = 8964328 ;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 29402.1 ;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}

			sprintf(qFilename[0], "%s%s%s", param.path, m.mName, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, m.mName);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){811,  4816,  1.0,    0.0,    0,     63.942862};
		}
	}
	if(m.id == 82){//FeH
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){821,  561,  1.0,    0.0,    0,     56.942762};
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
			char name[] = "56Fe-1H__Yueqi"; 
			sprintf(m.mName, "%s", "56Fe-1H__Yueqi");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nStates = 3563;
			m.nFiles = 1 ;
			m.ntcol = 3;
			m.npfcol = 2;
			m.NL[0] = 93040;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 7476;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){821,  561,  1.0,    0.0,    0,      56.942762};
			//version = 20160726
		}
		if(param.useHITEMP == 3){
			char name[] = "KuruczFeHfx"; 
			sprintf(m.mName, "%s", "KuruczFeHfx");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nFiles = 1 ;
			m.npfcol = 2;
			m.NL[0] = 111404;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 16139;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){821,  561,  1.0,    0.0,    0,      56.942762};
			//version = 20160726
		}
	}
	if(m.id == 83){//AlO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){831,  2716,  1.0,    0.0,    0,      42.976454};
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
			char name[] = "27Al-16O__ATP";
			sprintf(m.mName, "%s", "27Al-16O__ATP");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nStates = 94209;
			m.nFiles = 1 ;
			m.ntcol = 3;
			m.npfcol = 2;
			m.NL[0] = 4945580;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 35000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){831,  2716,  1.0,    0.0,    0,      42.976454};
			//version = 20160726
		}
	}
	if(m.id == 84){//SiO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){841,  2816,  1.0,    0.0,    0,      43.971842};
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
			char name[] = "28Si-16O__EBJT";
			sprintf(m.mName, "%s", "28Si-16O__EBJT");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nStates = 24306;
			m.nFiles = 1 ;
			m.ntcol = 4;
			m.npfcol = 2;
			m.NL[0] = 254675;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 6050;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){841,  2816,  1.0,    0.0,    0,      43.971842};
			//version = 20160726
		}
	}
	if(m.id == 85){//CaO
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){851,  4016,  1.0,    0.0,    0,      55.957506};
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
			char name[] = "40Ca-16O__VBATHY";
			sprintf(m.mName, "%s", "40Ca-16O__VBATHY");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nStates = 130660;
			m.nFiles = 1 ;
			m.ntcol = 4;
			m.npfcol = 2;
			m.NL[0] = 28417920;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 25000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){851,  4016,  1.0,    0.0,    0,      55.957506};
			//version = 20160726
		}
	}
	if(m.id == 86){//SiH
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){861,  281,  1.0,    0.0,    0,      28.98475156};
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
			char name[] = "28Si-1H__SiGHTLY";
			sprintf(m.mName, "%s", "28Si-1H__SiGHTLY");
			m.defaultL = 0.0700;
			m.defaultn = 0.500;
			m.nStates = 11785;
			m.nFiles = 1;
			m.ntcol = 4;
			m.npfcol = 2;
			m.NL[0] = 1724841;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 32000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){861,  281,  1.0,    0.0,    0,      28.98475156};
			//version = 20171101

		}
	}
	if(m.id == 87){//CaH
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){871,  401,  1.0,    0.0,    0,      40.970416};
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
			char name[] = "40Ca-1H__Yadin";
			sprintf(m.mName, "%s", "40Ca-1H__Yadin");
			m.defaultL = 0.0700 ;
			m.defaultn = 0.500 ;
			m.nStates = 1892;
			m.nFiles = 1 ;
			m.ntcol = 4;
			m.npfcol = 2;
			m.NL[0] = 26980;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 15278;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){871,  401,  1.0,    0.0,    0,      40.970416};
			//version = 20160726
		}
	}
	if(m.id == 88){//H3+
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//			 id	AFGL	Abundance	Q(296K)		gj	Molar Mass(g)
		m.ISO[0] = (Isotopologue){881,  111,  1.0,    0.0,    0,      3.023475};
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
			char name[] = "1H3_p__MiZATeP";
			sprintf(m.mName, "%s", "1H3_p__MiZATeP");
			m.defaultL = 0.07 ;
			m.defaultn = 0.500 ;
			m.nStates = 158721;
			m.nFiles = 1 ;
			m.ntcol = 3;
			m.npfcol = 3;
			m.NL[0] = 127542657;
			m.NLmax = 0;
			for(int i = 0; i < m.nFiles + 1; ++i){
				m.fileLimit[i] = i * 25000;
				m.NLmax = max(m.NLmax, m.NL[i]);
			}
			sprintf(qFilename[0], "%s%s%s", param.path, name, ".pf");
				sprintf(m.dataFilename[0], "%s%s.", param.path, name);
			m.nISO = 1;
			m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
			m.ISO[0] = (Isotopologue){881,  11,  1.0,    0.0,    0,      3.023475};
			//version = 20170330
		}
	}
       	if(m.id == 300){//Li 7
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){300,  7,  1.0,    0.0,    0,     7.016};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){300,  7,  1.0,    0.0,    0,     7.016};
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
		m.ISO[0] = (Isotopologue){400,  9,  1.0,    0.0,    0,     9.0122};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){400,  9,  1.0,    0.0,    0,     9.0122};
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
		m.ISO[0] = (Isotopologue){401,  9,  1.0,    0.0,    0,     9.0122};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){401,  9,  1.0,    0.0,    0,     9.0122};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		m.ISO[0] = (Isotopologue){601,  12,  1.0,    0.0,    0,     12};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){601,  12,  1.0,    0.0,    0,     12};
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
		m.ISO[0] = (Isotopologue){602,  12,  1.0,    0.0,    0,     12};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){602,  12,  1.0,    0.0,    0,     12};
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
		m.ISO[0] = (Isotopologue){700,  14,  1.0,    0.0,    0,     14.0031};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){700,  14,  1.0,    0.0,    0,     14.0031};
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
		m.ISO[0] = (Isotopologue){701,  14,  1.0,    0.0,    0,     14.0031};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){701,  14,  1.0,    0.0,    0,     14.0031};
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
		m.ISO[0] = (Isotopologue){702,  14,  1.0,    0.0,    0,     14.0031};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){702,  14,  1.0,    0.0,    0,     14.0031};
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
		m.ISO[0] = (Isotopologue){800,  16,  1.0,    0.0,    0,     15.9949};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){800,  16,  1.0,    0.0,    0,     15.9949};
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
		m.ISO[0] = (Isotopologue){801,  16,  1.0,    0.0,    0,     15.9949};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){801,  16,  1.0,    0.0,    0,     15.9949};
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
		m.ISO[0] = (Isotopologue){802,  16,  1.0,    0.0,    0,     15.9949};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){802,  16,  1.0,    0.0,    0,     15.9949};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		m.ISO[0] = (Isotopologue){1000,  20,  1.0,    0.0,    0,     19.9924};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1000,  20,  1.0,    0.0,    0,     19.9924};
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
		m.ISO[0] = (Isotopologue){1001,  20,  1.0,    0.0,    0,     19.9924};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1001,  20,  1.0,    0.0,    0,     19.9924};
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
		m.ISO[0] = (Isotopologue){1002,  20,  1.0,    0.0,    0,     19.9924};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1002,  20,  1.0,    0.0,    0,     19.9924};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		m.ISO[0] = (Isotopologue){1200,  24,  1.0,    0.0,    0,     23.985};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1200,  24,  1.0,    0.0,    0,     23.985};
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
		m.ISO[0] = (Isotopologue){1201,  24,  1.0,    0.0,    0,     23.985};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1201,  24,  1.0,    0.0,    0,     23.985};
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
		m.ISO[0] = (Isotopologue){1202,  24,  1.0,    0.0,    0,     23.985};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1202,  24,  1.0,    0.0,    0,     23.985};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
                m.nFiles = 1;           //number of data files
                m.NL[0] = 0;
                m.NLmax = 0;
                m.nISO = 1;
                m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
                //                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
                m.ISO[0] = (Isotopologue){1400,  28,  1.0,    0.0,    0,     27.9769};
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
                if(param.useHITEMP == 3){
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
                        m.ISO[0] = (Isotopologue){1400,  28,  1.0,    0.0,    0,     27.9769};
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
		m.ISO[0] = (Isotopologue){1401,  28,  1.0,    0.0,    0,     27.9769};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1401,  28,  1.0,    0.0,    0,     27.9769};
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
		m.ISO[0] = (Isotopologue){1402,  28,  1.0,    0.0,    0,     27.9769};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1402,  28,  1.0,    0.0,    0,     27.9769};
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
		m.ISO[0] = (Isotopologue){1500,  31,  1.0,    0.0,    0,     30.9737};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1500,  31,  1.0,    0.0,    0,     30.9737};
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
		m.ISO[0] = (Isotopologue){1501,  31,  1.0,    0.0,    0,     30.9737};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1501,  31,  1.0,    0.0,    0,     30.9737};
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
		m.ISO[0] = (Isotopologue){1502,  31,  1.0,    0.0,    0,     30.9737};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1502,  31,  1.0,    0.0,    0,     30.9737};
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
		m.ISO[0] = (Isotopologue){1600,  32,  1.0,    0.0,    0,     31.9721};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1600,  32,  1.0,    0.0,    0,     31.9721};
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
		m.ISO[0] = (Isotopologue){1601,  32,  1.0,    0.0,    0,     31.9721};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1601,  32,  1.0,    0.0,    0,     31.9721};
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
		m.ISO[0] = (Isotopologue){1602,  32,  1.0,    0.0,    0,     31.9721};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1602,  32,  1.0,    0.0,    0,     31.9721};
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
		m.ISO[0] = (Isotopologue){1700,  35,  1.0,    0.0,    0,     34.9689};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1700,  35,  1.0,    0.0,    0,     34.9689};
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
		m.ISO[0] = (Isotopologue){1701,  35,  1.0,    0.0,    0,     34.9689};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1701,  35,  1.0,    0.0,    0,     34.9689};
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
		m.ISO[0] = (Isotopologue){1702,  35,  1.0,    0.0,    0,     34.9689};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1702,  35,  1.0,    0.0,    0,     34.9689};
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
		m.ISO[0] = (Isotopologue){1800,  40,  1.0,    0.0,    0,     39.9624};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1800,  40,  1.0,    0.0,    0,     39.9624};
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
		m.ISO[0] = (Isotopologue){1801,  40,  1.0,    0.0,    0,     39.9624};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1801,  40,  1.0,    0.0,    0,     39.9624};
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
		m.ISO[0] = (Isotopologue){1802,  40,  1.0,    0.0,    0,     39.9624};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1802,  40,  1.0,    0.0,    0,     39.9624};
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
		m.ISO[0] = (Isotopologue){1900,  39,  1.0,    0.0,    0,     38.9637};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1900,  39,  1.0,    0.0,    0,     38.9637};
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
		m.ISO[0] = (Isotopologue){1901,  39,  1.0,    0.0,    0,     38.9637};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1901,  39,  1.0,    0.0,    0,     38.9637};
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
		m.ISO[0] = (Isotopologue){1902,  39,  1.0,    0.0,    0,     38.9637};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){1902,  39,  1.0,    0.0,    0,     38.9637};
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
		m.ISO[0] = (Isotopologue){2000,  40,  1.0,    0.0,    0,     39.9626};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2000,  40,  1.0,    0.0,    0,     39.9626};
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
		m.ISO[0] = (Isotopologue){2001,  40,  1.0,    0.0,    0,     39.9626};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2001,  40,  1.0,    0.0,    0,     39.9626};
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
		m.ISO[0] = (Isotopologue){2002,  40,  1.0,    0.0,    0,     39.9626};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2002,  40,  1.0,    0.0,    0,     39.9626};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		m.ISO[0] = (Isotopologue){2200,  48,  1.0,    0.0,    0,     47.9479};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2200,  48,  1.0,    0.0,    0,     47.9479};
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
		m.ISO[0] = (Isotopologue){2201,  48,  1.0,    0.0,    0,     47.9479};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2201,  48,  1.0,    0.0,    0,     47.9479};
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
		m.ISO[0] = (Isotopologue){2202,  48,  1.0,    0.0,    0,     47.9479};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2202,  48,  1.0,    0.0,    0,     47.9479};
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
		m.ISO[0] = (Isotopologue){2300,  51,  1.0,    0.0,    0,     50.944};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2300,  51,  1.0,    0.0,    0,     50.944};
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
		m.ISO[0] = (Isotopologue){2301,  51,  1.0,    0.0,    0,     50.944};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2301,  51,  1.0,    0.0,    0,     50.944};
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
		m.ISO[0] = (Isotopologue){2302,  51,  1.0,    0.0,    0,     50.944};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2302,  51,  1.0,    0.0,    0,     50.944};
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
		m.ISO[0] = (Isotopologue){2400,  52,  1.0,    0.0,    0,     51.9405};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2400,  52,  1.0,    0.0,    0,     51.9405};
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
		m.ISO[0] = (Isotopologue){2401,  52,  1.0,    0.0,    0,     51.9405};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2401,  52,  1.0,    0.0,    0,     51.9405};
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
		m.ISO[0] = (Isotopologue){2402,  52,  1.0,    0.0,    0,     51.9405};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2402,  52,  1.0,    0.0,    0,     51.9405};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		m.ISO[0] = (Isotopologue){2600,  56,  1.0,    0.0,    0,     55.9349};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2600,  56,  1.0,    0.0,    0,     55.9349};
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
		m.ISO[0] = (Isotopologue){2601,  56,  1.0,    0.0,    0,     55.9349};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2601,  56,  1.0,    0.0,    0,     55.9349};
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
		m.ISO[0] = (Isotopologue){2602,  56,  1.0,    0.0,    0,     55.9349};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2602,  56,  1.0,    0.0,    0,     55.9349};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
	if(m.id == 2800){//Ni 58
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2800,  58,  1.0,    0.0,    0,     57.9353};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2800,  58,  1.0,    0.0,    0,     57.9353};
			//version =  gfallwn08oct17.dat
		}
	}
        if(m.id == 2801){//Ni+ 58
                m.nFiles = 1;           //number of data files
                m.NL[0] = 0;
                m.NLmax = 0;
                m.nISO = 1;
                m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
                //                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
                m.ISO[0] = (Isotopologue){2801,  58,  1.0,    0.0,    0,     57.9353};
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
                if(param.useHITEMP == 3){
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
                        m.ISO[0] = (Isotopologue){2801,  58,  1.0,    0.0,    0,     57.9353};
                        //version =  gfallwn08oct17.dat
                }
        }
	if(m.id == 2802){//Ni+2 58
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2802,  58,  1.0,    0.0,    0,     57.9353};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2802,  58,  1.0,    0.0,    0,     57.9353};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2900){//Cu 63
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2900,  63,  1.0,    0.0,    0,     62.9296};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2900,  63,  1.0,    0.0,    0,     62.9296};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2901){//Cu+ 63
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2901,  63,  1.0,    0.0,    0,     62.9296};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2901,  63,  1.0,    0.0,    0,     62.9296};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 2902){//Cu+2 63
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){2902,  63,  1.0,    0.0,    0,     62.9296};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){2902,  63,  1.0,    0.0,    0,     62.9296};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3000){//Zn 64
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3000,  64,  1.0,    0.0,    0,     63.9291};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){3000,  64,  1.0,    0.0,    0,     63.9291};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3001){//Zn+ 64
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3001,  64,  1.0,    0.0,    0,     63.9291};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){3001,  64,  1.0,    0.0,    0,     63.9291};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 3002){//Zn+2 64
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){3002,  64,  1.0,    0.0,    0,     63.9291};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){3002,  64,  1.0,    0.0,    0,     63.9291};
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
		m.ISO[0] = (Isotopologue){3800,  88,  1.0,    0.0,    0,     87.9056};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){3800,  88,  1.0,    0.0,    0,     87.9056};
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
		m.ISO[0] = (Isotopologue){3801,  88,  1.0,    0.0,    0,     87.9056};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){3801,  88,  1.0,    0.0,    0,     87.9056};
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
		m.ISO[0] = (Isotopologue){3900,  89,  1.0,    0.0,    0,     88.9058};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){3900,  89,  1.0,    0.0,    0,     88.9058};
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
		m.ISO[0] = (Isotopologue){3901,  89,  1.0,    0.0,    0,     88.9058};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){3901,  89,  1.0,    0.0,    0,     88.9058};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4000){//Zr 90
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4000,  90,  1.0,    0.0,    0,     89.9043};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4000,  90,  1.0,    0.0,    0,     89.9043};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4001){//Zr+ 90
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4001,  90,  1.0,    0.0,    0,     89.9043};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4001,  90,  1.0,    0.0,    0,     89.9043};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4002){//Zr+2 90
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4002,  90,  1.0,    0.0,    0,     89.9043};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4002,  90,  1.0,    0.0,    0,     89.9043};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4100){//Nb 83
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4100,  83,  1.0,    0.0,    0,     82.9064};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4100,  83,  1.0,    0.0,    0,     82.9064};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4101){//Nb+ 83
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4101,  83,  1.0,    0.0,    0,     82.9064};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4101,  83,  1.0,    0.0,    0,     82.9064};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4102){//Nb+2 83
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4102,  83,  1.0,    0.0,    0,     82.9064};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4102,  83,  1.0,    0.0,    0,     82.9064};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4200){//Mo 98
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4200,  98,  1.0,    0.0,    0,     97.9054};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4200,  98,  1.0,    0.0,    0,     97.9054};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4201){//Mo+ 98
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4201,  98,  1.0,    0.0,    0,     97.9054};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4201,  98,  1.0,    0.0,    0,     97.9054};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4202){//Mo+2 98
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4202,  98,  1.0,    0.0,    0,     97.9054};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4202,  98,  1.0,    0.0,    0,     97.9054};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
	if(m.id == 4400){//Ru 104
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4400,  104,  1.0,    0.0,    0,     103.905};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4400,  104,  1.0,    0.0,    0,     103.905};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4401){//Ru+ 104
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4401,  104,  1.0,    0.0,    0,     103.905};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4401,  104,  1.0,    0.0,    0,     103.905};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 4402){//Ru+2 104
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){4402,  104,  1.0,    0.0,    0,     103.905};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4402,  104,  1.0,    0.0,    0,     103.905};
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		if(param.useHITEMP == 3){
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
		m.ISO[0] = (Isotopologue){4600,  106,  1.0,    0.0,    0,     105.904};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4600,  106,  1.0,    0.0,    0,     105.904};
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
		m.ISO[0] = (Isotopologue){4601,  106,  1.0,    0.0,    0,     105.904};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){4601,  106,  1.0,    0.0,    0,     105.904};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 5600){//Ba 138
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){5600,  138,  1.0,    0.0,    0,     137.905};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){5600,  138,  1.0,    0.0,    0,     137.905};
			//version =  gfallwn08oct17.dat
		}
	}
	if(m.id == 5601){//Ba+ 138
		m.nFiles = 1;		//number of data files
		m.NL[0] = 0;
		m.NLmax = 0;
		m.nISO = 1;
		m.ISO = (Isotopologue*)malloc(m.nISO * sizeof(Isotopologue));
		//                       id     AFGL    Abundance       Q(296K)         gj      Molar Mass(g)
		m.ISO[0] = (Isotopologue){5601,  138,  1.0,    0.0,    0,     137.905};
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
		if(param.useHITEMP == 3){
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
			m.ISO[0] = (Isotopologue){5601,  138,  1.0,    0.0,    0,     137.905};
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

