#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[]){


	char M[160];
	M[0] = 0;

	int nnMax = 0;
	int TMax = 63;
	const int PMax = 33;

	int dev = 0;

	int S = 1;

	//Read console input arguments
	for(int i = 1; i < argc; i += 2){
		if(strcmp(argv[i], "-M") == 0){
			sprintf(M, "%s", argv[i + 1]);
		}
		else if(strcmp(argv[i], "-nMax") == 0){
			nnMax = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-TMax") == 0){
			TMax = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-dev") == 0){
			dev = atoi(argv[i + 1]);
		}
		else if(strcmp(argv[i], "-tp") == 0){
			S = atoi(argv[i + 1]);
		}
		else{
			printf("Error: Console arguments not valid!\n");
			return 0;
		}
	}


	printf("%s %d\n", M, nnMax);

	if(strcmp(M, "") == 0){
		printf("Error, no -M value\n");
		return 0;
	}

	if(nnMax <= 0){
		printf("Error, no -nMax value\n");
		return 0;
	}

	FILE* scriptfile;
	char filename[160];


	//const int dtt = 1;
	const int dtt = TMax + 1;
	const int dpp = 34;
	const int dnn = nnMax + 1;

	//cut hybrid
	//const int dpp = 1;

	double P[PMax];

	for(int nn = 0; nn < nnMax; nn += dnn){
printf("%d %d\n", nn, nnMax);
	//for(int nn = 0; nn <= 5; ++nn){

		for(int tt = 0; tt <= TMax; tt += dtt){
			
			for(int pp = 0; pp <= PMax; pp += dpp){


				sprintf(filename, "script%02d_%02d_%02d", nn, tt, pp);
				scriptfile = fopen(filename, "w");


				//needed on slurm cluster
				fprintf(scriptfile, "#!/bin/bash\n");	
/*
				fprintf(scriptfile, "#SBATCH -J %02d_%02d_%02d\n", nn, tt, pp);	
				//fprintf(scriptfile, "#SBATCH --output ./slurm%s_%02d_%02d_%02d.out\n", M, nn, tt, pp);	
				//fprintf(scriptfile, "#SBATCH -n 1 --gres=gpu:gtx1080ti:1 -p gpu\n");	
				fprintf(scriptfile, "#SBATCH -n 1 --gres gpu:1 -p production\n");
				fprintf(scriptfile, "#SBATCH --time 1-0\n");

				//fprintf(scriptfile, "#SBATCH --qos=job_gpu\n");
				fprintf(scriptfile, "id=$SLURM_ARRAY_TASK_ID\n");
*/

				fprintf(scriptfile, "\n");	

				//temperature
				for(int i = 0; i < TMax + 1; ++i){
					if(i < 14)          fprintf(scriptfile, "T[%d]=%05g\n", i, 50.0 + i * 50.0); 
					else if(i < 14 + 8) fprintf(scriptfile, "T[%d]=%05g\n", i, 800.0 + (i - 14) * 100.0);
					else                fprintf(scriptfile, "T[%d]=%05g\n", i, 1700.0 + (i - 14 - 8) * 200.0);
				}
				fprintf(scriptfile, "\n");	



				int ni = 12; //number of magnitudes
				int nj = 3; //dex
				double pstart = -8.0;	//lowest exponent

				for(int i = 0; i < ni; ++i){
					for(int j = 0; j < nj; ++j){
						if(i == ni - 1 && j == 1) break;
						double p = pstart + 1.0 * i + 1.0/((double)(nj)) * j;
						P[i * nj + j] = pow(10.0, p) * 0.986923267;

						//printf("%d %g %g %g\n", i * nj + j, p, P, P * 0.986923267); //convert bar in atm
						fprintf(scriptfile, "P[%d]=%.20g\n", i * nj + j, P[i * nj + j]);
					}
				}
				fprintf(scriptfile, "\n");	

				for(int i = 0; i < ni; ++i){
					for(int j = 0; j < nj; ++j){
						if(i == ni - 1 && j == 1) break;
						int expp = (int)((pstart + 1 * i + 1.0/((double)(nj)) * j) * 100);

						if(expp < 0) fprintf(scriptfile, "expP[%d]=\"n%03d\"\n", i * nj + j, abs(expp)); //convert bar in atm
						else fprintf(scriptfile, "expP[%d]=\"p%03d\"\n", i * nj + j, abs(expp)); //convert bar in atm
					}
				}

				fprintf(scriptfile, "\n\n");


				//nu
				for(int i = 0; i < 55; ++i){
					fprintf(scriptfile, "NU[%d]=%05g\n", i, i * 1000.0);
				}

				// ************************************************************
				if(S == 1){		
					fprintf(scriptfile, "\n");	
		
//cut hybrid
//					fprintf(scriptfile, "cut=%g\n", 25.0 * P[pp]);

					fprintf(scriptfile, "for nu in {%d..%d..1}\n", nn, nn + dnn - 1);
					fprintf(scriptfile, "do\n");	
					fprintf(scriptfile, "  for t in {%d..%d..1}\n", tt, tt + dtt - 1);	
					fprintf(scriptfile, "  do\n");	
					fprintf(scriptfile, "    for p in {%d..%d..1}\n", pp, pp + dpp - 1);	
					fprintf(scriptfile, "    do\n");	
					fprintf(scriptfile, "      echo ${P[$p]}\n");	
					//fprintf(scriptfile, "      srun heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} -cut $cut} \n", dev, M);	
	//cut hybrid
					//fprintf(scriptfile, "      srun heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} \n", dev, M);	
					fprintf(scriptfile, "      ./heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} \n", dev, M);	
					fprintf(scriptfile, "    done\n");	
					fprintf(scriptfile, "  done\n");	
					fprintf(scriptfile, "done\n");	

					fprintf(scriptfile, "\n");	
					fprintf(scriptfile, "echo \"DONE\" \n");	
					fclose(scriptfile);
				}
				// ************************************************************


				// ************************************************************
				if(S == 2){
					fprintf(scriptfile, "t=$(($id %% %d))\n", dtt);
		
					fprintf(scriptfile, "\n");	
		
//cut hybrid
//					fprintf(scriptfile, "cut=%g\n", 25.0 * P[pp]);

					fprintf(scriptfile, "for nu in {%d..%d..1}\n", nn, nn + dnn - 1);
					fprintf(scriptfile, "do\n");	
					fprintf(scriptfile, "  for p in {%d..%d..1}\n", pp, pp + dpp - 1);	
					fprintf(scriptfile, "  do\n");	
					fprintf(scriptfile, "    echo ${P[$p]}\n");	
					//fprintf(scriptfile, "      srun heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} -cut $cut} \n", dev, M);	
	//cut hybrid
					//fprintf(scriptfile, "    srun heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} \n", dev, M);	
					fprintf(scriptfile, "    ./heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} \n", dev, M);	
					fprintf(scriptfile, "  done\n");	
					fprintf(scriptfile, "done\n");	

					fprintf(scriptfile, "\n");	
					fprintf(scriptfile, "echo \"DONE\" \n");	
					fclose(scriptfile);
				}
				// ************************************************************


				// ************************************************************
				if(S == 3){
					fprintf(scriptfile, "t=$(($id / %d))\n", dpp);
					fprintf(scriptfile, "p=$(($id %% %d))\n", dpp);
		
					fprintf(scriptfile, "\n");	
		
//cut hybrid
//					fprintf(scriptfile, "cut=%g\n", 25.0 * P[pp]);

					fprintf(scriptfile, "for nu in {%d..%d..1}\n", nn, nn + dnn - 1);
					fprintf(scriptfile, "do\n");	
					fprintf(scriptfile, "  echo ${P[$p]}\n");	
					//fprintf(scriptfile, "      srun heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} -cut $cut} \n", dev, M);	
	//cut hybrid
					//fprintf(scriptfile, "  srun heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} \n", dev, M);	
					fprintf(scriptfile, "  ./heliosk -dev %d -M %s -T ${T[$t]} -P ${P[$p]} -numin ${NU[$nu]} -numax ${NU[$(($nu+1))]} -name ${NU[$nu]}\"_\"${NU[$(($nu+1))]}\"_\"${T[$t]}\"_\"${expP[$p]} \n", dev, M);	
					fprintf(scriptfile, "done\n");	

					fprintf(scriptfile, "\n");	
					fprintf(scriptfile, "echo \"DONE\" \n");	
					fclose(scriptfile);
				}
				// ************************************************************
			}
		}
	}
	return 0;
}

