//gcc -c mersenne.c -lm
//gcc mersenne.o -o alg algorithm.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include"mersenne.h"

#define FAIL    0

/* Control Parameters of the search algorithm*/
int POP_SIZE;  /* The number of candidate solutions*/
int MAX_ITER; /*The number of iterations*/
int FUNCTION;
//Problem definitions
int DIM;    //number of problem variables
int lb; //lower bound of the variables
int ub; //upper bound of the variables
int RUN;  /*Algorithm can run many times in order to see its robustness*/

//Global variables
	double **pop; //[POP_SIZE][DIM];  //population of candidate solutions.
	double *fo; //[POP_SIZE];      //objective function value.
	double *best; //[DIM];           //best solution found
	double bestfo;         //best fo value
	int best_index;     	  //index for the best solution

//Functions declarations
double randon( double inferior, double superior);


//aux functions
int GetParameters(char **argv);
void showParameters();
void freeArrays();
void AllocArrays();

void prepararObjFunc(void);
double objfunc(double sol[]);

double randon( double inferior, double superior)
{
//  double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));
  double aux = (double)inferior + ((superior - inferior)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));
 
  return aux;
}


/*Main program of the search algorithm*/
int main(int argc, char **argv)
{
	int i, j, k, r;
	double *U; //[DIM]; // vetor para cálculo da func obj.

	srand(time(NULL));
	MT_seed();

	if (GetParameters(argv) == -1)	//read input file
		return 0;
	showParameters();

	AllocArrays();			//alloc arrays
	U = (double*)malloc(DIM * sizeof(double));

	prepararObjFunc();

	for (r=0;r<RUN;r++)
	{	
		printf("RUN: %d\n",r);
		//Init population
		for (j=0;j<POP_SIZE;j++)//each individual
		{
			fo[j]  = 0.0;
			for (k=0; k<DIM;k++) //each dimension
			{
				best[k] = 0.0;
				pop[j][k] = randon(lb,ub);
			}
		}
		bestfo = 0.0;
		best_index = 0;
		//Objective function calculation
		for (i = 0;i<POP_SIZE;i++)
		{
			for (j = 0;j<DIM;j++)
				U[j] = pop[i][j];
	
			fo[i] = objfunc(U);
		}
	
		//Best current solution identification.
		bestfo = fo[0];
		best_index = 0;
		for(i=0;i<POP_SIZE;i++)
		{
			if (fo[i]<=bestfo) //for minimization
			{
	        		bestfo=fo[i];
	        		for(j=0;j<DIM;j++)
	        		   best[j]=pop[i][j];
	
				best_index = i;
	        	}
		}

///***** Implement the search algorithm <<<<<<<<<<<<<<<<<<<<<<
//Loop de Iterações.
	
		//just to test.
		for (j=0;j<POP_SIZE;j++)
		{
			printf("Variables: ");		
			for (k=0; k<DIM;k++) //variables
			{
				printf("%7.4f ",pop[j][k]);
			}
			printf(" Fo:");
			printf("%8.4f \n",fo[j]);
		}
		printf("\nBest current solution: ");
		for (k=0; k<DIM;k++) //variables
		{
			printf("%7.4f ",best[k]);
		}
		printf(" Fo:");
		printf("%8.4f \n",bestfo);

	}//end FOR RUN
	freeArrays();
	free(U);
}//end MAIN

/* Dynamic array allocation */
/* mol and sites */
void AllocArrays()	
{
	int i,j;

	pop = malloc (POP_SIZE*sizeof(double*));
        for (j = 0;j < POP_SIZE; j++)
		pop[j] = (double*)malloc (DIM * sizeof(double));

	fo = (double*)malloc (POP_SIZE*sizeof(double));

	best = (double*)malloc (DIM*sizeof(double));
}


/* Free arrays */
void freeArrays()
{
	int i,j;
	free(fo);
	free(best);
	for (i = 0; i < POP_SIZE; i++)
		free(pop[i]);
	free(pop);
}


//
/*Input file reading*/
int GetParameters(char **argv)
{
	FILE *file = fopen( argv[1], "r" );
    
	if (file == 0)
	{
    	printf( "Could not open ini file! Usage ./<exec> <file.in>\n" );
	return -1;
	}
    else 
    {
		ffscanf("RUN", file, "%d", &RUN);
		ffscanf("MAX_ITER", file, "%d", &MAX_ITER);
		ffscanf("POP_SIZE", file, "%d", &POP_SIZE);
		ffscanf("DIM", file, "%d", &DIM); 
		ffscanf("FUNCTION",file, "%d", &FUNCTION);
	return 1;
    }
    fclose(file);
}


void showParameters()
{
	printf("***PARAMETERS***\n");
	printf("RUNS = %d\n", RUN);
	printf("MAX_ITER = %d\n", MAX_ITER);
	printf("POP_SIZE = %d\n", POP_SIZE);
	printf("DIM = %d\n", DIM);
	switch (FUNCTION)
	{
		case 0:
			printf("FUNCTION = %s\n","Rastrigin");
			break;
		case 1:
			printf("FUNCTION = %s\n","Schaffer");
			break;
		case 2:
			printf("FUNCTION = %s\n","Griewank");
			break;
		case 3:
			printf("FUNCTION = %s\n","Ackley");
			break;
		case 4:
			printf("FUNCTION = %s\n","Rosenbrock");
			break;
		case 5:
			printf("FUNCTION = %s\n","Sphere");
			break;
	}
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

void prepararObjFunc(void)
{
    switch (FUNCTION) {
    case 0: //Rastrigin
        lb = -5.12;
        ub = 5.12;
        break;
    case 1: //Schaffer
        lb = -100.00;
        ub = 100.00;
        break;
    case 2: //Griewank
        lb = -600.00;
        ub = 600.00;
        break;
    case 3: //Ackley
/*
-	Dimension: n arbitrary
-       Domain:   -32 <= | x_i | <= 32.0
-       Minimum 0 at x_i = 0.0
*/
        lb = -32.00;
        ub = 32.00;
        break;
    case 4: //Rosenbrock
        lb = -30.00;
        ub = 30.00;
        break;
    case 5: //Sphere
	lb = -100.00;
	ub = 100.00;
	break;
    default:
        printf("Info: Invalid function.\n") ;
        exit(0);
    }
}

double objfunc(double sol[])
{
    int j, i;
    double top = 0.00 , top1 = 0.00, top2 = 0.00;
    double aux = 0.0;
    double aux1 = 0.0;

    switch (FUNCTION) {
    case 0: //Rastrigin

        for(j=0;j<DIM;j++)
        {
            top=top+(pow(sol[j],(double)2)-10*cos(2*M_PI*sol[j])+10);
        }
        return top;

    case 1: //Schaffer

        top1 = 0;
        for(j=0;j<DIM;j++)
        {
        top=top+(pow(sol[j],(double)2));
        }
        top = pow(top,(double)0.25);
        for(j=0;j<DIM;j++)
        {
        top1=top1+(pow(sol[j],(double)2));
        }
        top1=pow(top1,(double)0.1);
        top1 = pow(sin(50*top1),(double)2) +1.0;

        return top*top1;

    case 2: //Griewank

        top=0;
        top1=0;
        top2=1;
        for(j=0;j<DIM;j++)
        {
        top1=top1+pow((sol[j]),(double)2);
        top2=top2*cos((((sol[j])/sqrt((double)(j+1)))*M_PI)/180);
        }
        top=(1/(double)4000)*top1-top2+1;

        return top;

    case 3: //Ackley

        for (i = 0; i < DIM; i++)
        {
        aux += sol[i]*sol[i];
        }
        for (i = 0; i < DIM; i++)
        {
        aux1 += cos(2.0*M_PI*sol[i]);
        }

        return (-20.0*(exp(-0.2*sqrt(1.0/(float)DIM*aux)))-exp(1.0/(float)DIM*aux1)+20.0+exp(1));

    case 4: //Rosenbrock

        for (i = 0; i < DIM-1; i++)
        {
            top=top+100.*pow((sol[i+1] - pow(sol[i],2.)),2) + pow((1. - sol[i]),2);
        }

       return top;

    case 5: //Sphere
	for(j=0;j<DIM;j++)
	{
		top=top+sol[j]*sol[j];
	}

	return top;

    default:
        printf("Info: Invalid function.\n") ;
        exit(0);
    }
}

