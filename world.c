/* Copyright Dassault Systemes, 1999, 2010 */

/*=========================================================================*\
 *
 * This template file can be used to integrate an optimization algorithm
 * into Isight via an external executable program. 
 *
 * To integrate your optimization algorithm, find the "main" routine at 
 * the bottom of this file, and insert a call to your algorithm function 
 * in the indicated location "INSERT A CALL TO YOUR OPTIMIZATION ROUTINE HERE"
 *
 * At the top of the file there are several utility routines for allocating
 * memory, reading a technique options file, reading a problem setup file,
 * executing a single design point, executing a set of design points.
 *
 * All communication between Isight and this executable program is done
 * via writing and reading files.
 *
\*=========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <math.h>


#define NUM_COMMAND_LINE_ARGS 6
#define DV_TYPE_REAL 0
#define DV_TYPE_INTEGER 1
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define LOG_FILE_NAME "woa-test-log-file.txt"
#define pi 3.14

static int     numDesignVariables        = 0;
static int     numConstraints            = 0;
static int     numEqConstraints          = 0;
static char**  designVariableNames       = NULL;
static double* designVariableValues      = NULL;
static double* designVariableLowerBounds = NULL;
static double* designVariableUpperBounds = NULL;
static int*    designVariableTypes       = NULL;

static char*   goRunSignalFileName       = NULL;
static char*   inputValsFileName         = NULL;
static char*   resultsFileName           = NULL;
static FILE*   logFile                   = NULL;

/* 
 * NOTE! Change the following value to match the number of 
 * technique options in your optimization plug-in 
 */
#define MAX_NUM_TECHNIQUE_OPTIONS 15
static double  techniqueOptions[MAX_NUM_TECHNIQUE_OPTIONS];
static char techniqueOptionNames[MAX_NUM_TECHNIQUE_OPTIONS][128];

/*=========================================================================*\
 * Utility function that calls calloc() and in case of memory
 * allocation failure prints a message to stderr and exits the program.
\*=========================================================================*/
static void * allocateArray(size_t nElem, size_t elemSize, char* arrayName) {
	void * memPointer = calloc(nElem, elemSize);
	if(memPointer == NULL) {
		fprintf(stderr, "Failed to allocate memory for %s\n", arrayName);
		exit(1);
	}
	return memPointer;
}

/*=========================================================================*\
 * This function reads technique options from a file.
 * The format of the technique options file is as follows:
 *
 *     -------------------
 *    |Val_1 OptionName_1 |
 *    |Val_2 OptionName_2 |
 *    |...                |
 *     -------------------
 *
 * NOTE! the order of technique options in the file is not guaranteed to
 * match the displayed order in Isight. When using the values of the
 * technique options, always check the name of the option.
 *
\*=========================================================================*/
static void readTechniqueOptionsFile(char* fName) {

	FILE* f;
	int i;

	/* open the file for reading */
	if((f = fopen(fName, "r")) == NULL) {
		fprintf(stderr,"Failed to open file \"%s\" for reading\n", fName);
		exit(1);
	}

	fprintf(logFile,"Reading TECHNIQUE OPTIONS from \"%s\"\n",fName);

	/* read the values and names, store in the global arrays */
	for(i=0; i<MAX_NUM_TECHNIQUE_OPTIONS; i++) {
		if(feof(f)) {
			break;
		}
		if(fscanf(f,"%lf %[^\n]\n", &techniqueOptions[i],techniqueOptionNames[i]) != 2) {
			fprintf(stderr,"Invalid file format, \"%s\", line=%d\n",fName,i+1);
			exit(1);
		}
		fprintf(logFile,"  %d. %s = %lg\n", i+1, techniqueOptionNames[i], techniqueOptions[i]);
	}
}

/*=========================================================================*\
 * This function reads problem setup from a file.
 * The format of the problem setup file is as follows. The first 3 lines
 * are the sizing parameters (# of design variables, # of constraints (all),
 * and # of equality constraints. The rest of the file - design variable
 * definitions, one line for each design variable:
 *
 *     -------------------------------------------- 
 *    |n1 Number of design variables               |
 *    |n2 Number of constraints (total)            |
 *    |n3 Number of equality (target) constraints  |
 *    |Val_1 LB_1 UB_1 Type_1 VariableName_1       |
 *    |Val_2 LB_2 UB_2 Type_2 VariableName_2       |
 *    |...                                         |
 *     -------------------------------------------- 
 *
 *     where
 *       Val_1 - starting value of design variable #1 (DV#1)
 *       LB_1 and UB_1 - lower/upper bounds of DV#1
 *       Type_1 - value type of DV#1 (0 - real, 1 - integer)
 *       VariableName_1 - name of DV#1 (need not be used)
 *
\*=========================================================================*/
static void readProblemSetupFile(char* fName) {

	FILE* f;
	char strBuffer[256];
	int i;
	int numVals;

	/* open the file for reading */
	if((f = fopen(fName, "r")) == NULL) {
		fprintf(stderr,"Failed to open file \"%s\" for reading\n", fName);
		exit(1);
	}

	/* read sizing parameters first, lines 1,2,3 */
	fprintf(logFile,"Reading PROBLEM SETUP data from \"%s\"\n",fName);
	numVals = fscanf(f, "%d %*[^\n]\n %d %*[^\n]\n %d %*[^\n]\n", 
			&numDesignVariables, &numConstraints, &numEqConstraints);
	if(numVals != 3) {
		fprintf(stderr,"Invalid file format, \"%s\", " 
			"lines 1-3 must have problem sizing parameters (integer values)\n",
			fName);
		exit(1);
	}

	/* allocate global arrays for storing design variable and constraint data */
	designVariableNames = (char**)allocateArray(numDesignVariables, sizeof(char*), 
			"design variable names array");
	for(i=0; i<numDesignVariables; i++) {
		designVariableNames[i] = (char*)allocateArray(128, sizeof(char), 
				"design variable names array");
	}
	designVariableValues = (double*)allocateArray(numDesignVariables, sizeof(double), 
			"design variable values array");
	designVariableLowerBounds = (double*)allocateArray(numDesignVariables, sizeof(double), 
			"design variable lower bounds array");
	designVariableUpperBounds = (double*)allocateArray(numDesignVariables, sizeof(double), 
			"design variable upper bounds array");
	designVariableTypes = (int*)allocateArray(numDesignVariables, sizeof(int), 
			"design variable types array");

	/* read data into global arrays */
	for(i=0; i<numDesignVariables; i++) {
		numVals = fscanf(f, "%lf %lf %lf %d %[^\n]\n",
			&designVariableValues[i],&designVariableLowerBounds[i],
			&designVariableUpperBounds[i],&designVariableTypes[i],
			designVariableNames[i]);
		if(numVals < 5) {
			fprintf(stderr,"Invalid file format, \"%s\":%d\n",fName,i+1);
			exit(1);
		}
		fprintf(logFile,"%d. %s = %lg, [%lg...%lg] type=%d\n", i+1, designVariableNames[i],
				designVariableValues[i], designVariableLowerBounds[i], 
				designVariableUpperBounds[i],designVariableTypes[i]);
	}

	/* close the file */
	fclose(f);
}

/*=========================================================================*\
 *
 * This routine runs a single design point analysis; constraint and 
 * objective values are stored in the supplied arrays/pointers
 *
 * args:
 * inputValues - array of input values for subflow run
 * constraintValues - array for storing constraint values; if this
 *         value is NULL, the constraint values are not stored
 * objectiveValue - pointer to a variable for storing objective val
 * penaltyValue - pointer to a variable for storing penalty value
 *
\*=========================================================================*/
static void executeDesignPoint(double* inputValues, double* constraintValues, 
		double* objectiveValue, double* penaltyValue) {

	int i;
	FILE* signalFile;
	FILE* inputValsFile;
	FILE* resultsFile;
	double tempVal;

	/* write input values into the specified file */
	fprintf(logFile,"Opening input values file \"%s\"\n", inputValsFileName);
	remove(inputValsFileName);
	inputValsFile = fopen(inputValsFileName, "w+");
	if(inputValsFile == NULL) {
		fprintf(logFile,"Failed to open file \"%s\" for writing\n", 
				inputValsFileName);
		fprintf(stderr, "Failed to open file \"%s\" for writing\n", 
				inputValsFileName);
		exit(1);
	}
	fprintf(logFile,"Requesting execution of 1 design point, input values:\n");
	for(i=0; i<numDesignVariables; i++) {
		tempVal = inputValues[i];
		/* for integer variable types, round the value */
		if(designVariableTypes[i] == DV_TYPE_INTEGER) {
			tempVal = round(tempVal);
		}
		fprintf(inputValsFile, "%.15g ", tempVal);
		fprintf(logFile,"  DV#%d = %.15g\n",i+1,tempVal);
	}
	fprintf(inputValsFile, "\n");
	fclose(inputValsFile);

	/* 
	 * create the signal file to tell Isight that the input
	 * values file is ready and a sublfow evaluation must be done 
	 */
	signalFile = fopen(goRunSignalFileName,"w");
	fclose(signalFile);

	/* wait for Isight to delete the signal file */
	fprintf(logFile,"Waiting for Isight to delete the signal file...\n");
	while(1) {
		signalFile = fopen(goRunSignalFileName,"r");
		if(signalFile == NULL) {
			/* file not found, Isight is done executing subflow runs */
			break;
		}
		fclose(signalFile);
		Sleep(1000);	
	
	}

	/* try to open results file and read constraints and objective values */
	resultsFile = fopen(resultsFileName,"r");
	if(resultsFile == NULL) {
		fprintf(stderr, "Failed to open file \"%s\" for readin\n", 
				resultsFileName);
		exit(1);
	}

	/* read constraints and obj values from the results file */
	fprintf(logFile,"Constraint/Obj values:\n");
	for(i=0; i<numConstraints; i++) {
		if(fscanf(resultsFile, "%lf", &tempVal) != 1) {
			fprintf(stderr, 
				"Failed to read constraint value #%d from file \"%s\"\n", 
				i+1, resultsFileName);
			exit(1);
		}
		if(constraintValues != NULL) {
			constraintValues[i] = tempVal;
		}
		fprintf(logFile,"  CON#%d = %lg\n", i+1, tempVal);
	}
	if(fscanf(resultsFile, "%lf", objectiveValue) != 1) {
		fprintf(stderr, "Failed to read objective value from file \"%s\"\n", 
				resultsFileName);
		exit(1);
	}
	fprintf(logFile,"  OBJ = %lg\n", *objectiveValue);
	if(fscanf(resultsFile, "%lf", penaltyValue) != 1) {
		fprintf(stderr, "Failed to read penalty value from file \"%s\"\n", 
				resultsFileName);
		exit(1);
	}
	fprintf(logFile,"  PEN = %lg\n", *penaltyValue);

	/* remove the results file for the next iteration */
	remove(resultsFileName);
}

/*=========================================================================*\
 * This function runs multiple design points analysis; constraint and 
 * objective values are stored in the supplied arrays
 *
 * args:
 * numDesignPoints - number of design points to be executed
 * inputValues - matrix of input values for subflow runs
 * constraintValues - matrix for storing constraint values; if this
 *         value is NULL, the constraint values are not stored
 * objectiveValues - array for storing objective values
 * penaltyValues - array for storing penalty values; if this value
 *         is NULL, penalty values are not stored
 *
\*=========================================================================*/
static void executeDesignPointSet(int numDesignPoints, double** inputValues, 
        double** constraintValues, double* objectiveValues, 
        double* penaltyValues) {

    int i,row;
    FILE* signalFile;
    FILE* inputValsFile;
    FILE* resultsFile;
    double tempVal;

    /* write all input values into the file */
    remove(inputValsFileName);
    fprintf(logFile,"Opening input values file \"%s\"\n", inputValsFileName); 
    fflush(logFile);
    inputValsFile = fopen(inputValsFileName, "w");
    if(inputValsFile == NULL) {
        fprintf(stderr, "Failed to open file \"%s\" for writing\n", 
                inputValsFileName);
        exit(1);
    }
    fprintf(logFile,"Requesting execution of %d design points, input values:\n", numDesignPoints); 
    fflush(logFile);
    for(row=0; row<numDesignPoints; row++) {
        fprintf(logFile,"-%d-\n",row+1); 
        fflush(logFile);
        for(i=0; i<numDesignVariables; i++) {
            tempVal = inputValues[row][i];
            /* for integer variable types, round the value */
            if(designVariableTypes[i] == DV_TYPE_INTEGER) {
                tempVal = round(tempVal);
            }
            fprintf(logFile,"  %d. \"%s\" = %.15g\n",i+1,designVariableNames[i],tempVal);
            fflush(logFile);
            fprintf(inputValsFile, "%.15g ", tempVal);
        }
        fprintf(inputValsFile, "\n");
    }
    fclose(inputValsFile);

    /* 
     * create the signal file to tell Isight that a 
     * sublfow evaluation must be done 
     */
    signalFile = fopen(goRunSignalFileName,"w");
    fclose(signalFile);

    /* wait for Isight to delete the signal file */
    fprintf(logFile,"Waiting for Isight to delete the signal file...\n"); 
    fflush(logFile);
    while(1) {
        signalFile = fopen(goRunSignalFileName,"r");
        if(signalFile == NULL) {
            /* file not found, Isight is done executing subflow runs */
            break;
        }
        fclose(signalFile);

        Sleep(1000);    /* windows only - sleep for 1000 milliseconds */

    }

    /* try to open results file and read constraints and objective values */
    resultsFile = fopen(resultsFileName,"r");
    if(resultsFile == NULL) {
        fprintf(stderr, "Failed to open file \"%s\" for readin\n", 
                resultsFileName);
        exit(1);
    }

    /* read constraints and obj values from the results file */
    fprintf(logFile,"Constraint/Obj values:\n"); 
    fflush(logFile);
    for(row=0; row<numDesignPoints; row++) {
        fprintf(logFile,"-%d-\n",row+1); 
        fflush(logFile);
        for(i=0; i<numConstraints; i++) {
            if(fscanf(resultsFile, "%lf", &tempVal) != 1) {
                fprintf(stderr, 
                    "Failed to read constraint value #%d from file \"%s\"\n", 
                    i+1, resultsFileName);
                exit(1);
            }
            if(constraintValues != NULL) {
                constraintValues[row][i] = tempVal;
            }
            fprintf(logFile,"  CON#%d = %lg\n", i+1, tempVal); 
            fflush(logFile);
        }
        if(fscanf(resultsFile, "%lf", &tempVal) != 1) {
            fprintf(stderr, "Failed to read objective value from file \"%s\"\n", 
                    resultsFileName);
            exit(1);
        }
        if(objectiveValues != NULL) {
            objectiveValues[row] = tempVal;
        }
        fprintf(logFile,"  OBJ = %lg\n", tempVal); 
        fflush(logFile);
        if(fscanf(resultsFile, "%lf", &tempVal) != 1) {
            fprintf(stderr, "Failed to read penalty value from file \"%s\"\n", 
                    resultsFileName);
            exit(1);
        }
        if(penaltyValues != NULL) {
            penaltyValues[row] = tempVal;
        }
        fprintf(logFile,"  PEN = %lg\n", tempVal); 
        fflush(logFile);
    }

    /* remove the results file for the next iteration */
    remove(resultsFileName);
}


/*=========================================================================*\
 * InRange - check that the coordinate is in range and correct
 * value - design variable value to be checked
 * minimum - lower bound
 * maximum - upper bound
 * returns double - adjusted design variable value
\*=========================================================================*/
static double InRange(double value, double minimum, double maximum) {

    if (value < minimum) {
		if((2*minimum-value)<maximum){
			return 2*minimum-value;
		}else{
			return maximum;
		}
        
    }
    if (value > maximum) {		
		if((2*maximum-value)>minimum){
			return 2*maximum-value;
		}else{
			return minimum;
		}       
    }
    return value;
}

/*=========================================================================*\
 * evaluateDesign - evaluate design point
 * x - double[] array with design variable values
 * returns double - objective+penalty for the design point
\*=========================================================================*/
static double evaluateDesign(double* x) {

    double objectiveValue, penaltyValue;
    executeDesignPoint(x, NULL, &objectiveValue, &penaltyValue);
    return (penaltyValue + objectiveValue);
}



/*=========================================================================*\
 * woa test main
\*=========================================================================*/

static double* world(int nvars, double* startpt, double* bl, double* bu,int popsize, int evalmax, double umut, int* parmType) {
				
	double* xbest = NULL;
    int i, row,j;	
	double fbefore;		
    double** xMatrix = NULL;
	double* objectiveValuesVector = NULL;
	double* penaltyValuesVector = NULL;
	double* popfit = NULL;
	double a,A,C,l,p,po;
	int nn,n1,n2,n3;
	
	int nmut;
	int* idx;
	int temp;
	int numIterations = evalmax/popsize-1;
    double CR,F;
	
	 xbest = (double*)allocateArray(nvars, sizeof(double), "bestx array");
     popfit = (double*)allocateArray(popsize, sizeof(double), "best popfit array");

    xMatrix = (double**)allocateArray(popsize, sizeof(double*), 
			"input values matrix");
	for(i=0; i<popsize; i++) {
		xMatrix[i] = (double*)allocateArray(nvars,
				sizeof(double), "input values matrix");
	}
	objectiveValuesVector = (double*)allocateArray(popsize, sizeof(double), 
			"objective values array");
	penaltyValuesVector = (double*)allocateArray(popsize, sizeof(double), 
			"penalty values array");

     /* Starting point analysis */
    fprintf(logFile, "Starting point analysis...\n");
    fflush(logFile);
	fbefore = evaluateDesign(startpt);
		
    for(i=0;i<popsize;i++){
		for(row=0;row<nvars;row++){
			xMatrix[i][row]= ((double)rand())/((double)RAND_MAX)*(bu[row]-bl[row])+bl[row];
			xbest[row] = startpt[row];
		}
	}
	
	executeDesignPointSet(popsize, xMatrix, NULL, objectiveValuesVector, penaltyValuesVector);
	fprintf(logFile, "intial population evaluated end...\n");
    fflush(logFile);
	
	for(i=0;i<popsize;i++){
		popfit[i]= objectiveValuesVector[i]+ penaltyValuesVector[i];
		if(popfit[i] < fbefore){
			fbefore = popfit[i];
			for(row=0;row<nvars;row++){
				xbest[row] = xMatrix[i][row];
			}
		}	
	}
	
	fprintf(logFile, "start main loop \n");
    fflush(logFile);	

	for(j=0; j<numIterations; j++) {

		fprintf(logFile,"======Iteration #%d==========\n", j+1);
		 fflush(logFile);
	         a = 2.0*(numIterations-j)/numIterations;
			 F = 0.5 + 0.5*j/numIterations;
			 CR = 0.9-0.5*j/numIterations;
			 
         for(i=0;i<popsize;i++){		
			 A = 2*a*rand()/RAND_MAX-a;
			 C = 2.0*rand()/RAND_MAX;
			 l = 2.0*rand()/RAND_MAX-1;
			 p = 1.0*rand()/RAND_MAX;
			 po = 1.0*rand()/RAND_MAX;
			 if(p<0.5){
				 if(fabs(A)<1){					 
						 if( po < CR){
							 for(row=0;row<nvars;row++){
								xMatrix[i][row] = InRange((xbest[row]-A*fabs(C*xbest[row]-xMatrix[i][row])),bl[row],bu[row]);
						      }
						 }
						 else{    
						       for(row=0;row<nvars;row++){
									n1 = rand()%popsize;	
									n2 = rand()%popsize;								
									n3 = rand()%popsize;	

									xMatrix[i][row] = InRange((xMatrix[n1][row]+(xMatrix[n2][row]-xMatrix[n3][row])*F),bl[row],bu[row]);
								}												 
					     }		
				 }
				 else{
					 nn = rand()%popsize;
					 for(row=0;row<nvars;row++){
						 xMatrix[i][row] = InRange(xMatrix[i][row]-A*fabs(C*xMatrix[nn][row]-xMatrix[i][row]),bl[row],bu[row]);
					 }
				 }
			 }
			 else{
				 for(row=0;row<nvars;row++){
					 xMatrix[i][row] = InRange(fabs(xbest[row]-xMatrix[i][row])*exp(l)*cos(2*pi*l)+xbest[row],bl[row],bu[row]);
				 }
			 }
		 }
		 
		 
		  nmut = round(umut*popsize/3);	
		  fprintf(logFile," nmut= %d\n", nmut); 
		   fflush(logFile);
	      idx = (int*)allocateArray(nmut, sizeof(int), "index array");
	   
	            for(i=1;i<nmut;i++){
			         idx[i] = round(popsize/3.0)+ rand()%(round(popsize/3.0))+1;
			         temp =  idx[i];					
					 
					for(row=0;row<nvars;row++){
			          xMatrix[temp][row] = rand()/RAND_MAX*(bu[row]-bl[row])+bl[row];
					}	
				}			
	            	  	  
 	  
		 executeDesignPointSet(popsize, xMatrix, NULL, objectiveValuesVector, penaltyValuesVector);
		 
		 	for(i=0; i<popsize; i++){
		         popfit[i]= objectiveValuesVector[i]+ penaltyValuesVector[i];
		         if(popfit[i]<fbefore){
			           fbefore = popfit[i];
					for(row=0;row<nvars;row++){
				       xbest[row] = xMatrix[i][row];
					 
			        }
					 
		        }	
	        }
			
	    fprintf(logFile," bestf= %lg\n", fbefore); 
		fflush(logFile);
	}
		
	return xbest;
}




/*=========================================================================*\
 *
 * Entry point for the program. Isight will create a technique options file 
 * and a problem setup file, then call this program to execute your
 * optimization algorithm.
 *
 * The command syntax is as follows:
 *    "this.exe options problemSetup inputVals goRunSignal runsDoneSignal runResults"
 *
 * where:
 *    <0> this.exe - the name of this executable program
 *    <1> options - name of the file with technique options
 *    <2> problemSetup - name of the file with problem size parameters and design 
 *        variable values, bounds, types, and names
 *    <3> inputsVals - name of the file where this program must write input 
 *        vals for next subflow run(s), Isight will then try to read this
 *        file and execute the subflows for each row of values
 *    <4> goRunSignal - name of the file that this program must create to 
 *        tell Isight that the input values file is ready and to run subflow(s);
 *        Isight will delete this file when the output results are ready!
 *    <5> runResults - name of the file where Isight will put values of 
 *        constraints, objective, and penalty for all executed subflow runs
 *
\*=========================================================================*/
int main(int argc, char** argv) {

	int i;
	double* x = NULL;	
	char* techniqueOptionsFileName = NULL;
	char* problemSetupFileName = NULL;

	/* create local variables for your technique options here */

	int popsize = 25;
	int evalmax = 500;
	double umut = 0.2;

	/* validate command line arguments */
	if(argc < NUM_COMMAND_LINE_ARGS) {
		fprintf(stderr,"Invalid command, argc = %d, expected argc = %d\n", 
				argc, NUM_COMMAND_LINE_ARGS);
		exit(1);
	}

	/* set file names from the command line arguments */
	techniqueOptionsFileName  = argv[1];
	problemSetupFileName      = argv[2];
	inputValsFileName         = argv[3];
	goRunSignalFileName       = argv[4];
	resultsFileName           = argv[5];

	/* 
	 * open technique log file for all log messages;
	 * this file name must be specified in SDK-Generator, on the 
	 * "Technique Options" tab
	 * */
	logFile = stdout;
	logFile = fopen(LOG_FILE_NAME,"w");

	/* 
	 * debugging printout; 
	 * all stdOut prints can be seen as DEBUG 
	 * messages in Isight job log 
	 */
	fprintf(stdout,"**** Easy Native Optimization Technique ****\n");
	fprintf(stdout,"**** args = \n");
	for(i=0; i<argc; i++) {
		fprintf(stdout,"    %s\n", argv[i]);
	}
	fprintf(stdout,"\n");

	/* 
	 * read technique options from the file;
	 * the values of technique options and their names will
	 * be stored in the global arrays
	 */
	readTechniqueOptionsFile(techniqueOptionsFileName);

	/* 
	 * check option names in the global array, and assign 
	 * technique option values to local variables 
	 * based on their names
	  */ 
	for(i=0; i<MAX_NUM_TECHNIQUE_OPTIONS; i++) {

		if(strcmp(techniqueOptionNames[i],"popsize")==0) {
            popsize = (int)techniqueOptions[i];
        }

		if(strcmp(techniqueOptionNames[i],"max evaluations")==0) {
            evalmax = (int)techniqueOptions[i];
        }	
		
		if(strcmp(techniqueOptionNames[i],"percentage of mutation")==0) {
            umut = techniqueOptions[i];
        }	
			
		
	}


	/* 
	 * read problem setup data from the file;
	 * the values of all problem setup parameters will
	 * be stored in global arrays
	 */
	readProblemSetupFile(problemSetupFileName);

	/*
	 * INSERT A CALL TO YOUR OPTIMIZATION ROUTINE HERE
	 * and remove the rest of the code from this "main" routine
	 * except the last statement "exit(0)".
	 */

	/* allocate arrays for inputs, constraints, etc. */
	x = (double*)allocateArray(numDesignVariables, sizeof(double), 
			"input values array");
	    for(i=0; i<numDesignVariables; i++) {
        x[i] = designVariableValues[i];
    }		
			

 /* call woa algorithm */
    world(numDesignVariables, x, designVariableLowerBounds, designVariableUpperBounds,
        popsize,evalmax,umut, designVariableTypes);

	/* 
	 * exit with return code = 0; 
	 * any non-zero return code means an error
	 */
	exit(0);
}

