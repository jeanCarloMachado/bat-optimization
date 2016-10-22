//gcc -c mersenne.c -lm
//gcc mersenne.o -o alg algorithm.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include"mersenne.h"

#ifndef TRUE
#define TRUE    1
#define FALSE   0
#endif

#define FAIL    0

#define SEQ_MAX_LENGTH 100  //Protein max length
#define PI 3.14159265

#define lb -PI //-180.0 /*lower bound of the parameters. */
#define ub PI //180.0 /*upper bound of the parameters. lb and ub can be defined as arrays for the problems of which parameters have different bounds*/

#define DIM_MAX 105  //105 for 55 aminoacids

/* Control Parameters of the search algorithm*/
//#define POP_SIZE 40 /* The number of candidate solutions*/
//#define MAX_ITER 700000 /*The number of iterations*/
//#define RUN 100 /*Algorithm can run many times in order to see its robustness*/
//Problem definitions
//#define DIM 2 					//number of problem variables
#define TOL 0.00001 //Tolerancia 10^-5

char sequence[SEQ_MAX_LENGTH]; //protein sequence (max length = SEQ_MAX_LENGTH)
int ProteinSize;
int NP;
int maxCycle;
int maxrun;
int D;

typedef struct 
{
	int Tipo; 		//tipo A (1) ou B (-1)
	double x, y, z; //posição no plano cartesiano
}amino;

amino *Aminoacido;


typedef struct estruturaMorcego 
{
        double x[DIM_MAX]; //Posicao
        double v[DIM_MAX]; //Velocidade
        double f; //Frequencia
        double fitness;
        double A;
        double r;
       double r0;
} estruturaMorcego; 

estruturaMorcego *morcego; //morcego[NP]

double *solucaoTemp; //solucaoTemp[DIM];

time_t ti; 
double tp;


/* aux functions */
int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer);
void GetParameters(char **argv);
void showParameters();
void freeArrays();
void AllocArrays();
void report(int run, double bestf, int iter, double tp);
void report_variables(int r);
void objfgen(int run, int iter, double bestf);

double AB3D_Function(double V[]);		//3D AB function
void Rotation(double x,double y,double z,double a,double b,double c,double u,double v,double w,double t,double *rx,double *ry,double *rz);

#define FUNCAO 0
/* 0=Rastrigin                  DIM=2   lb=-5.12   ub=5.12   Result:502920 +- 134880(100%)  Otimo:0.00
   1=Schaffer                   DIM=?   lb=-100.00 ub=100.00 Result:?                       Otimo:0.00
   2=Griewank                   DIM=2   lb=-600.00 ub=600.00 Result:391680  +- 189280(100%) Otimo:0.00
   3=Ackley                     DIM=128 lb=-30.00  ub=30.00  Result:277320  +- 92680(100%)  Otimo:0.00
   4=Rosenbrock                 DIM=16  lb=-2.048  ub=2.048  Result:316920  +- 131720(100%) Otimo:1.00*/

//Encontrar novo Loudness. Para Yang melhor 0.70 - 0.90. Quanto maior ALFA, mais devagar loudness cai

//Metaheuristic Bat Algorithms p7  ... in our implementation 0.90
//A New Metaheuristic Bat-Inspired Algorithm p5 ... we have used 0.90 in our simulations
#define ALFA 0.90

//Controlar a exponencial negativa, para controlar a busca local. Para Yang melhor 0.70 - 0.90

//Metaheuristic Bat Algorithms p7  ... in our implementation 0.50
//A New Metaheuristic Bat-Inspired Algorithm p5 ... we have used 0.90 in our simulations
#define L 0.90

//Metaheuristic Bat Algorithms p5 ... In our implementation, we will use 0 and 100,
//                                    depending the domain size of the problem of interest
//A New Metaheuristic Bat-Inspired Algorithm p5 ... In our implementation, we will use 0 and 100,
#define F_MIN 0.00
#define F_MAX 100.00

//Global variables
double **pop;
double *fo;
double *best;
double *gbest;
double *U;


double bestfo;         //best fo value
double gBestfo = 0.00; //Fitness da melhor solucao entre todas as execucoes
double gQtdFoRun; //Quantidade de avaliacoes por execucao
double gR; //Execucao da melhor solucao encontrada
double gt; //Tempo total da execucao


    //A New Metaheuristic Bat-Inspired Algorithm p3 ... we can also use A0 = 1 and Amin = 0
    //                                           p5 ... A0 can typically be [1, 2]
    const double LOUDNESS_MIN = 0.00;
    const double LOUDNESS_MAX = 1.00;

    //Metaheuristic Bat Algorithms p3 ... r [0, 1]
    //A New Metaheuristic Bat-Inspired Algorithm p3 ... r [0, 1]
    //                                           p5 ... can be around zero
const double PULSE_RATE_MIN = 0.00;
const double PULSE_RATE_MAX = 1.00;

    //Ajustar a frequencia
    //Metaheuristic Bat Algorithms p5 ... [0, 1] is a random vector
    //A New Metaheuristic Bat-Inspired Algorithm p4 ... [0, 1] is a random vector drawn from a uniform distribution
    const double BETA_MIN = -1.00;
    const double BETA_MAX = 1.00;

    //Random walk
    const double E_MIN = -0.10;
    const double E_MAX = 0.10;

//Functions declarations
double randon( double inferior, double superior);
void batAlgorithmMatlab();

double randon( double inferior, double superior)
{
  double aux = (double)inferior + ((superior - inferior)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));
   return aux;
}


/*Main program of the search algorithm*/
int main(int argc, char **argv)
{
	int i, j, k, r, it;
	int iterbest;

	srand(time(NULL));
	MT_seed();

	GetParameters(argv);	//read input file
	showParameters();

	AllocArrays();			//alloc arrays

	for(i=0;i<ProteinSize;i++)
   	{
		if (sequence[i] == 'A' || sequence[i] == 'a')
			Aminoacido[i].Tipo = 1;
		if (sequence[i] == 'B' || sequence[i] == 'b')
			Aminoacido[i].Tipo = -1;
   	}

 	for (r=0;r<maxrun;r++)
	{
		ti = time(NULL);

		//Init population
		for (j=0;j<NP;j++)//each individual
		{
			fo[j]  = 0.0;
			for (k=0; k<D;k++) //each dimension
			{
				best[k] = 0.0;
				pop[j][k] = randon(lb,ub);
			}
		}
		gQtdFoRun = 0.00;

		for (it=0;it<maxCycle;it++)
	   	{
			//Objective function calculation
			for (i = 0;i<NP;i++)
			{
				for (j = 0;j<D;j++)
					U[j] = pop[i][j];
	
				fo[i] = AB3D_Function(U);
			}
	
			//Best current solution identification.
			bestfo = fo[0];

			for(i=0;i<NP;i++)
			{
				if (fo[i]<bestfo) //for minimization
				{
    	    		bestfo=fo[i];
    	    		for(j=0;j<D;j++)
    	    		   best[j]=pop[i][j];

   		        }
			}

    	    batAlgorithmMatlab(); //BAT ALGORITHM


	        //Atualizar melhor entre todas as execuoes
	        if (it == 0)
	        {
	            gBestfo = bestfo;
	            gR = r;
	            for(j=0;j<D;j++)
	            {
	                gbest[j] = best[j];
	            }
	
	        }

			bestfo = fo[0];
			for(i=1;i<NP;i++)
			{
				if (fo[i] < bestfo)
				{
					for(j=0;j<D;j++)
					{
				        best[j]=pop[i][j];
					}
					bestfo = fo[i];
				}
			}	


			if (bestfo < gBestfo)  //minimization.
			{	
				gBestfo=bestfo;
	        	for(j=0;j<D;j++)
	        	   gbest[j]=best[j];

				iterbest = it; //iteration best fitness
			}	

			//report obj function value
			objfgen(r, it, gBestfo); //output --> fitness plot

		} //end for iter

		//report
		tp = (double) (time(NULL) - ti);
		report(r, gBestfo, iterbest, tp);
		report_variables(r);

	}//end FOR RUN

	freeArrays();
	return 1;
}//end MAIN



void batAlgorithmMatlab()
{
    int n=0, k=0, d=0;
    double fitnessTotal=0.00;
    int tMelhorMorcego=0;
    double t;
    int contTempo, indiceAleatorio;
	double solucaoTempFitness;
    double totalLoudness = 0.00, mediaLoudness = 0.00;

    //Criar morcego
    for (n=0; n < NP; n++) 
	{
        for (d=0; d < D; d++) 
		{
            morcego[n].x[d] = pop[n][d]; //Atribuir valores gerados pelo framework
            morcego[n].v[d] = 0.00;
            morcego[n].f = 0.00;
        }
        morcego[n].fitness = fo[n];
        morcego[n].A = 1.00;
        morcego[n].r = 0.00;
        morcego[n].r0 = randon(PULSE_RATE_MIN, PULSE_RATE_MAX);

        //Atualiza indice melhor morcego da execucao
        if (morcego[n].fitness < bestfo) 
		{
            bestfo = morcego[n].fitness;
            for(d=0; d < D; d++) 
			{
                best[d] = morcego[n].x[d];
            }

            tMelhorMorcego = 0;
        }
    }

    t = 0.00;
    contTempo = 0;



    //Iniciar busca
    //for ( t = 1; t <=maxCycle ; t++ ) {
//    while (bestfo > TOL) 
//	{
//        if (t > (double)maxCycle) 
//		{
//            break;
 //       }

        //Para acompanhar a execucao
//       if (contTempo == 10000) {
//            contTempo = 0;
//            printf("Execucao;%.0f;Tempo;%.0f;Melhor fitness;%8.10f;Quantidade Avaliacao;%.0f;\n", r, t,  bestfo, gQtdFoRun);
//        }

        fitnessTotal = 0.00;
        for (n=0; n < NP; n++) 
		{
            for (d=0; d < D; d++) 
			{
                /*Gerar novas solucoes com os ajustes abaixo.
                Obs: as solucoes sao temporarias ate satisfazer*/

                morcego[n].f = F_MIN + (F_MAX - F_MIN) * randon(BETA_MIN, BETA_MAX);

                //Definir velocidade
                morcego[n].v[d] = morcego[n].v[d] + (morcego[n].x[d] - best[d]) * morcego[n].f;
                if (morcego[n].v[d] > ub)
                    morcego[n].v[d] = randon(lb, ub);
                else if (morcego[n].v[d] < lb)
                    morcego[n].v[d] = randon(lb, ub);

                //Definir posicao
                solucaoTemp[d] = morcego[n].x[d] + morcego[n].v[d];
            }

            //Pulse rate
            if (randon(PULSE_RATE_MIN, PULSE_RATE_MAX) > morcego[n].r) 
			{
                totalLoudness = 0;
                for(k=0; k < NP; k++) 
				{
                    totalLoudness += morcego[k].A;
                }
                mediaLoudness = totalLoudness / (double) NP;

                for (d=0; d < D; d++) 
				{
                    solucaoTemp[d] = best[d] + randon(E_MIN, E_MAX) * mediaLoudness;
                }
            }

            //Generate a new solution by flying randomly
            indiceAleatorio = MT_randInt(D-1);
            solucaoTemp[indiceAleatorio] = randon(lb,ub);

            gQtdFoRun += 1.00; //Numero de avaliacoes
            solucaoTempFitness = 0.00;
            solucaoTempFitness = AB3D_Function(solucaoTemp);
            if ((solucaoTempFitness <= morcego[n].fitness) && (randon(LOUDNESS_MIN, LOUDNESS_MAX) < morcego[n].A)) 
			{
                for (d=0; d < D; d++) 
				{
                    morcego[n].x[d] = solucaoTemp[d];
                }
                morcego[n].fitness = solucaoTempFitness;

                //Aumentar pulse rate
                morcego[n].r = morcego[n].r0 * (1.00 - exp(-L * t));

                //Diminuir loudness
                morcego[n].A = ALFA * morcego[n].A;
            }

            //Atualiza indice melhor morcego
            if (solucaoTempFitness <= bestfo) 
			{
                bestfo = solucaoTempFitness;
                for(d=0; d < D; d++) 
				{
                    best[d] = solucaoTemp[d];
                }

                tMelhorMorcego = t;
            }
            fitnessTotal += morcego[n].fitness;
        } //end FOR NP


        t += 1.00;
        gt = t;
        contTempo += 1;
   // } //end FORmaxCycle

    //Atribuir valores dos morcegos para o framework
    for (n=0; n < NP; n++) 
	{
        for (d=0; d < D; d++) 
		{
            pop[n][d] = morcego[n].x[d];
        }
        fo[n] = morcego[n].fitness;
    }
}



/* Dynamic array allocation */
/* mol and sites */
void AllocArrays()	
{
	int i;

	Aminoacido = (amino *) malloc ((ProteinSize) * sizeof (amino));

	pop = malloc (NP * sizeof(double *));
	for (i = 0; i < NP; i++)
		pop[i] = malloc (D * sizeof(double));

	fo = (double *) malloc (NP * sizeof (double));
	best = (double *) malloc (D * sizeof (double));
	gbest = (double *) malloc (D * sizeof (double));
	U = (double *) malloc (D * sizeof (double));

	morcego = (estruturaMorcego *) malloc (NP * sizeof (estruturaMorcego));
	solucaoTemp = (double *) malloc (D * sizeof (double));

}

/* Free arrays */
void freeArrays()
{
	int i;

	free(Aminoacido);

	for (i = 0; i < NP; i++)
		free(pop[i]);	
	free(pop);	

	free(fo);		
	free(best);
	free(gbest);		
	free(U);
	free(morcego);
	free(solucaoTemp);
}

/*Input file reading*/
void GetParameters(char **argv)
{
	FILE *file = fopen( argv[1], "r" );
    
	if (file == 0)
	{
    	printf( "Could not open file!\n" );
	}
    else 
    {
		ffscanf("sequence", file, "%s", &sequence);
		ProteinSize = strlen(sequence);
		D = 2 * ProteinSize - 5; 				//i.e.: 13 ProteinSize --> D = 21
		ffscanf("NP", file, "%d", &NP);
		ffscanf("maxCycle", file, "%d", &maxCycle);
		ffscanf("runtime", file, "%d", &maxrun);
    }
    fclose(file);
}


void showParameters()
{
	printf("***PARAMETERS***\n");
	printf("sequence = %s\n", sequence);
	printf("ProteinSize = %d\n", ProteinSize);
	printf("NP = %d\n", NP);
	printf("maxCycle = %d\n", maxCycle);
	printf("maxrun = %d\n", maxrun);
	printf("****************\n");
}



int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer)
{
	char buffer[100];
	int len;
    int commentflag = 0;
	char *pch;
	char *pch2;
	
	do
	{
		if(fgets(buffer, 99, fp) == 0) return FAIL;
		buffer[99] = '\0';
		len = strlen(buffer);
	    if (buffer[len - 1] == '\n') buffer[len - 1] = '\0';

	    switch (commentflag) 
		{
		    case 0:
				if (strstr(buffer, "/*") != 0) 
				{
			    	commentflag = 1;
			    	if (strstr(buffer, "*/") != 0)
						commentflag = 2;
				}
				break;

		    case 1:
				if (strstr(buffer, "*/") != 0)
				    commentflag = 2;
				break;

			    case 2:
				if (strstr(buffer, "/*") != 0) 
				{
				    commentflag = 1;
				    if (strstr(buffer, "*/") != 0)
					commentflag = 2;
				}
				else
				    commentflag = 0;
				break;
			    }	
	}while(commentflag != 0);	

	//separate field name: token = "="
	if (strstr (buffer, fieldname) != 0)
	{
		pch = strtok (buffer,"=");

		while (pch != NULL)
		{
			pch2 = pch;
    		pch = strtok (NULL, "= ");
		}
		sscanf(pch2, format, inbuffer);
		return 1;//ok
	}
	else return 0; 
}

/******** Output files ********/

//Report : overall results
//obs.: besteverf = best ever obj function value
void report(int run, double bestf, int iter, double tp)
{
	FILE *file;	
	char *file_name;					

	file_name = (char *)malloc(100*sizeof(char));
	sprintf(file_name,"report_%s.txt", sequence); 
	file = fopen(file_name, "a+");
	free(file_name);  	

	if(run == 0) fprintf(file, "Run\tBest\tIteration(best)\tTp(s)\n");
	
	fprintf(file, "%d\t%lf\t%d\t%lf\n", run, bestf, iter, tp);

	fclose(file);
}

//report variables (rotation and torsion angles)
void report_variables(int r)
{
	FILE *file;	
	char *file_name;					
	int j;

	file_name = (char *)malloc(100*sizeof(char));
	sprintf(file_name,"variables_run%d_%s.txt", r, sequence);
	file = fopen(file_name, "a+");

	fprintf(file, "Variables (Best result): \n");

	for(j=0;j<D;j++)
	{
		fprintf(file, "gbest[%d]: %lf\n",j+1,gbest[j]);
	}	

	fclose(file);

}

//Obj function x iteration
//obs.: bestf = best obj function value
void objfgen(int r, int iter, double bestf)
{
	FILE *file;	
	char *file_name;					

	file_name = (char *)malloc(100*sizeof(char));
	sprintf(file_name,"fitnessgen_run%d_%s.txt", r, sequence);
	file = fopen(file_name, "a+");
	free(file_name);  	

	if(iter == 0) fprintf(file, "Iteration\tFitness(best)\n");
	
	fprintf(file, "%d\t%lf\n", iter, bestf);

	fclose(file);
}

double AB3D_Function(double V[])
{

    double w1, w2, LJ, LJ_aux;
	double r, rx, ry, rz, u1,u2,u3,v1,v2,v3;
    int i,j;
	double nx,ny,nz,vx,vy,vz,tx,ty,tz,ux,uy,uz; 

    int dim = ProteinSize - 2;


//avaliacao debug
/*
V[0] = -2.00123; //-1.99143;//-2.00123;// * 180/PI;
V[1] = 2.00765;//-2.005;//2.00765;// * 180/PI;
V[2] = -1.98636;//-1.97603;//-1.98636;// * 180/PI;
V[3] = -1.9938;//-1.9807;//-1.9938;// * 180/PI;
V[4] =  2.01299;//2.00377;//2.01299;// * 180/PI;
V[5] = -1.97374;//-1.97429;//-1.97374;// * 180/PI;
V[6] = 1.98693;//-1.9888;//1.98693;// * 180/PI;
V[7] = 1.98266;//-1.9702;//1.98266;// * 180/PI;
V[8] =  2.00189;//-1.98978;//2.00189;//* 180/PI;
V[9] =  2.01195;//2.00129;//2.01195;// * 180/PI;
V[10] =  -1.97379;//-1.9735;//-1.97379;// * 180/PI;
V[11] =   1.36669;//1.78119;//1.36669;// * 180/PI;
V[12] =  -0.282944;//0.298257;//-0.282944;// * 180/PI;
V[13] = -0.773825;// 0.776777;//-0.773825;// * 180/PI;
V[14] = -1.3366;//1.30329;//-1.3366;// * 180/PI;
V[15] =  0.310387;//-0.301256;//0.310387;// * 180/PI;
V[16] =  -2.50039;//-0.624863;//-2.50039;// * 180/PI;
V[17] =  2.44781;//  0.684238;//2.44781;// * 180/PI;
V[18] =  -2.30847;//-0.839701;//-2.30847;// * 180/PI;
V[19] = 1.34701;//-1.31482;//1.34701;// * 180/PI;
V[20] =  -0.351089;//0.33129;//-0.351089;// * 180/PI;
*/
//



	//polar --> cartesian
    Aminoacido[0].x = Aminoacido[0].y = Aminoacido[0].z = 0;
    Aminoacido[1].x = Aminoacido[1].z = 0;  Aminoacido[1].y = 1;

	nx=0;
	ny=0; 
	nz=1;

    for (i = 1; i < ProteinSize-1; i++)
    {  
		if(V[i-1]>PI)
			V[i-1] = randon(0.0,PI);
		if(V[i-1]<-PI)
			V[i-1] = randon(-PI,0.0);
		if(V[i+dim-1]>PI)
			V[i+dim-1] = randon(0.0,PI);
		if(V[i+dim-1]<-PI)
			V[i+dim-1] = randon(-PI,0.0);	


		tx=0;ty=0;tz=0;
		vx = Aminoacido[i].x - Aminoacido[i-1].x;
		vy = Aminoacido[i].y - Aminoacido[i-1].y;
		vz = Aminoacido[i].z - Aminoacido[i-1].z;


     Rotation(Aminoacido[i].x+vx, Aminoacido[i].y+vy, Aminoacido[i].z+vz, Aminoacido[i].x, Aminoacido[i].y, Aminoacido[i].z,nx,ny,nz,V[i-1],&tx,&ty,&tz);

	 if (i == 1)
	 	 Rotation(tx,ty,tz,Aminoacido[i].x,Aminoacido[i].y,Aminoacido[i].z,vx,vy,vz,0,&Aminoacido[i+1].x,&Aminoacido[i+1].y,&Aminoacido[i+1].z);		
	 else
	 	Rotation(tx,ty,tz,Aminoacido[i].x,Aminoacido[i].y,Aminoacido[i].z,vx,vy,vz,V[(i-2)+dim],&Aminoacido[i+1].x,&Aminoacido[i+1].y,&Aminoacido[i+1].z);


		ux = Aminoacido[i+1].x - Aminoacido[i].x;
	   	uy = Aminoacido[i+1].y - Aminoacido[i].y;
   		uz = Aminoacido[i+1].z - Aminoacido[i].z;

  		nx = vy*uz-vz*uy;
  		ny = vz*ux-vx*uz;
  		nz = vx*uy-vy*ux;

    }

	//bond angles
    w1 = 0;
    for (i = 1; i < ProteinSize-1; i++) 
	{
     	u1 = Aminoacido[i].x - Aminoacido[i-1].x;
     	u2 = Aminoacido[i].y - Aminoacido[i-1].y;
     	u3 = Aminoacido[i].z - Aminoacido[i-1].z;
     	v1 = Aminoacido[i+1].x - Aminoacido[i].x;
    	v2 = Aminoacido[i+1].y - Aminoacido[i].y;
     	v3 = Aminoacido[i+1].z - Aminoacido[i].z;
     	w1 += (u1*v1 + u2*v2 + u3*v3)/(sqrt(pow(u1,2) + pow(u2,2) + pow(u3,2)) * sqrt(pow(v1,2) + pow(v2,2) + pow(v3,2)));
	}


	//torsion angles
    w2 = 0;
    for (i = 1; i < ProteinSize-2; i++)

    {   
     	u1 = Aminoacido[i].x - Aminoacido[i-1].x;
    	u2 = Aminoacido[i].y - Aminoacido[i-1].y;
     	u3 = Aminoacido[i].z - Aminoacido[i-1].z;
     	v1 = Aminoacido[i+2].x - Aminoacido[i+1].x;
     	v2 = Aminoacido[i+2].y - Aminoacido[i+1].y;
     	v3 = Aminoacido[i+2].z - Aminoacido[i+1].z;
     	w2 = w2 - ((u1*v1 + u2*v2 + u3*v3)/(sqrt(pow(u1,2) + pow(u2,2) + pow(u3,2))*sqrt(pow(v1,2) + pow(v2,2)+pow(v3,2))))/2; 
    }

	//Lennard-Jones
	LJ = 0;
    for (i = 0; i < ProteinSize-2; i++)
    {   for (j = i + 2; j < ProteinSize; j++)
        {   
			rx = Aminoacido[i].x - Aminoacido[j].x;
			ry = Aminoacido[i].y - Aminoacido[j].y;
			rz = Aminoacido[i].z - Aminoacido[j].z;
			r = sqrt(rx * rx + ry * ry + rz * rz);

			 LJ_aux = 4 * (1/pow(r, 12) - 1/pow(r, 6));

			if ((Aminoacido[i].Tipo != Aminoacido[j].Tipo) || (Aminoacido[i].Tipo == -1 && Aminoacido[j].Tipo == -1))
			{	
				 LJ_aux *= 0.5;
			}
			
			LJ +=  LJ_aux;
        }
    }

//debug
//	printf("w1=%lf\nw2=%lf\nLJ=%lf\n", w1, w2, LJ); 	printf("Fitness = %lf\n", w1 + w2 + LJ);
//

    return(w1 + w2 + LJ);

}

void Rotation(double x,double y,double z,double a,double b,double c,double u,double v,double w,double t,double *rx,double *ry,double *rz)
{
   double u2 = u*u;
   double v2 = v*v;
   double w2 = w*w;
   double cost = cos(t);
   double sint = sin(t);
   double l2 = u2+v2+w2;
   double l = sqrt(l2);

   *rx = (a*(v2+w2)+u*(-b*v-c*w+u*x+v*y+w*z)+(-a*(v2+w2)+u*(b*v+c*w-v*y-w*z)+(v2+w2)*x)*cost+l*(-c*v+b*w-w*y+v*z)*sint)/l2;
   *ry = (b*(u2+w2)+v*(-a*u-c*w+u*x+v*y+w*z)+(-b*(u2+w2)+v*(a*u+c*w-u*x-w*z)+(u2+w2)*y)*cost+l*(c*u-a*w+w*x-u*z)*sint)/l2;
   *rz = (c*(u2+v2)+w*(-a*u-b*v+u*x+v*y+w*z)+(-c*(u2+v2)+w*(a*u+b*v-u*x-v*y)+(u2+v2)*z)*cost+l*(-b*u+a*v-v*x+u*y)*sint)/l2;
}

