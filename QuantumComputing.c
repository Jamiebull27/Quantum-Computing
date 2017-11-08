/* Computational Physics and Modelling- Project 2, Option D **
** A program to simulate operations of a quantum computer **

** Programming projects referred to in the code can be found in: **
** D. Candela, **
** Undergraduate computational physics projects on quantum computing **
** Am. J. Phys. 83, 688 (2015) **
** http://dx.doi.org/10.1119/1.4922296 **

** 28/01/2017 **
** Jamie Bull */

#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define F_NAME_MAX 256    //Max number of characters allowed for an output file name, need to change check_file, project_1 if this is changed




/* Stores values required to perform operations on a set of qubits */
typedef struct qubits{
	
	int n;    //Number of qubits.
	int m;    //Number of elements (2^n), stored here for simplicity in functions.
	int **matrix;    //Stores values corresponding to energy level of each qubit, 1 being excited state, 0 being ground state
					 //Needs  m x n  elements.
	double complex *element;    //Amplitude for each element, needs m elements.
	
}qubits;


/* Union to store int and double complex in sparse matrix, will always know the type being stored */
typedef union sparseMatrix {
    int intVal ;
	double complex cDoubleVal;
} sparseMatrix;


/* Checks if a file allready exists, gives user option to overwrite or create new file */
static void check_file( char file[F_NAME_MAX] ){
	
	if( access(file, F_OK) == -1 ){
		return;    //File does not allready exist
	}
	else{
		int validChoice = 0;
		char choice;
		while( validChoice == 0 ) {
			
		    printf("\n%s all-ready exists would you like to:\n"
		           "1. Overwrite %s\n"
		    	   "2. Create a new file\n", file, file);
		    scanf(" %c", &choice);
		    if( choice == '1' ){
				validChoice = 1;
		    	FILE *f=fopen(file, "w");    /*Clearing file*/
		    	if( f == NULL ){
		    		fprintf(stderr, "Error opening %s\n", file);
		    		exit(errno);
		    	}
		    	else{
					fclose(f);
		    		return; 
		    	}
		    }
		    else if( choice == '2' ){
				validChoice = 1;
		        printf("Enter a new name for your file. (With extension, max %i characters)\n", F_NAME_MAX);
				scanf(" %256s", file);    //Needs to be changed if F_NAME_MAX is changed.
				if( access(file, F_OK)!=-1 ){
					validChoice = 0;
				}
	        }
			else{
				printf("Invalid choice, please try again.\n");
			}
	    }
	}
	return;
}


/* Displays a qubit using Dirac bra-ket notation */
static void display_wave_function( qubits *wf ){
	
	char** displayElement = malloc( wf->m * sizeof(char*) );
	int i;
	for( i=0; i<wf->m; i++ ){
		displayElement[i] = malloc( (wf->n + 2.0) * sizeof(char) + 1.0 );
	}
	
	int j;
	for( i=0; i<wf->m; i++ ){
		/* Creating "Ket" vectors to display matrix */
		displayElement[i][0] = '|';
		for( j=0; j<wf->n; j++ ){
		    if( wf->matrix[i][j] == 0 ){
			    displayElement[i][j+1] = '0';
		    }
			else{
				displayElement[i][j+1] = '1';
			}
		}
		displayElement[i][wf->n+1] = '>';
	}
	/* Displaying "Ket" vectors */
	int nDisplayed = 0;
	printf("|Psi> = ");
	for( i=0; i<wf->m; i++ ){
		if( wf->element[i] != 0 ){
			if( nDisplayed != 0 ){
				printf("\n      + ");
			}
			if( wf->element[i] != 1 ){
				if( cimag( wf->element[i] ) == 0.0 ){
					printf("%g", creal(wf->element[i]));
				}
				else{
					printf("%g + %gi", creal(wf->element[i]), cimag(wf->element[i]));
				}
			}
			for ( j=0; j<wf->n+2; j++ ){
				printf("%c", displayElement[i][j]);
			}
			nDisplayed++;
		}
	}
	printf("\n");
	
	for( i=0; i<wf->m; i++ ){
		free(displayElement[i]);
	}
	free(displayElement);
}


/* Sets the all qubits in the wavefunction to the ground state */
static void set_wf_gs( qubits *wf ){
	
	int i;
	for( i=1; i<wf->m; i++ ){
		wf->element[i] = 0.0 + 0.0*I;
	}
	wf->element[0] = 1.0;
}


/* Initialises all ellements of the qubits structure, WARNING: contains malloc without free */
static void initialise_qubits( qubits *wf ){
	
	int i, j;
	
	printf("Enter the number of qubits to use in the calculations.\n");
	scanf(" %i", &wf->n);
	wf->m = pow( 2, wf->n );
	wf->matrix = malloc( wf->m * sizeof(int*) );
	if( wf->matrix == NULL ){
		fprintf( stderr, "Error allocating memory %i.\n", errno );
		exit(errno);
	}
	
	for( i=0; i<wf->m; i++ ){
		wf->matrix[i] = malloc( wf->n * sizeof(int) );
		if( wf->matrix[i] == NULL ){
		fprintf( stderr, "Error allocating memory %i.\n", errno );
		exit(errno);
		}
	}
	
	wf->element = malloc( wf->m * sizeof(double complex) );
	if( wf->element == NULL ){
		fprintf( stderr, "Error allocating memory %i.\n", errno );
		exit(errno);
	}
	
	/* Creating binary elements for wf->matrix */
	for( i=0; i<wf->m; i++ ){
		int decimal = i;
		for( j=wf->n-1; j>=0; j--){    //Counting down to ensure first qubit is on the left.
			wf->matrix[i][j] = decimal % 2;
			decimal = decimal / 2;
		}
	}
	
	set_wf_gs( wf );    //Initialise to ground state.
	printf("Qubits will start in state:\n");
	display_wave_function( wf );
}


/* Takes user input to scan a qubit from user input */
static void scan_elements( qubits *wf ){

	double buffer;
	
	printf("\nEnter the %i real components on the amplitude for the states of the qubit\n", wf->m);
	       
	int i;
	for( i=0; i<wf->m; i++ ){
		scanf(" %lg", &buffer);
		wf->element[i] = buffer + 0*I;
	}
	printf("You have entered:\n");
	display_wave_function( wf );
}


/* Converst a probability amplitude to a probability */
static double calc_probability( double complex amplitude ){

	return pow( creal(amplitude), 2.0 ) + pow( cimag(amplitude), 2.0 );
}


/* Rounds a double to the nearest integer */
static int round_nearest_int( double x ){
	
	if( 0.0 <= (x - (int)x)    &&    (x - (int)x) < 0.5 ){
		return (int)x;
	}
	return (int)x + 1;
}


/* Random number generator, produces U[low, high] */
static double rng(double low, double high){
	
	return low+(rand()/(double)RAND_MAX)*(high-low);
}


/* Randomly evolves a set of qubits through time, simulating a measurement of the system */
static void evolve_qubits( qubits *wf, double *evolvedQubit ){
	
	int probStart = 0, i, j;
	double qubitBuffer[wf->m], probSpacing[(wf->m)+1];
	
	probSpacing[0] = 0.0;
	for( i=0; i<wf->m; i++ ){
		/* Split up range from 0 to probStart (total probability) into sections scaled by the probability of each state */
		qubitBuffer[i] = 0.0;
		probSpacing[i+1] = calc_probability( wf->element[i] );
		probStart += calc_probability( wf->element[i] );    //Calculate the total probability at the start to ensure that it stays the same throughout the calculation.
	}
	probStart = round_nearest_int( probStart );    //Must be integer to prevent problems int the loop.
	for( i=0; i<probStart; i++ ){
		double rand = rng(0.0, probStart);
		if( rand == 0.0 ){
			qubitBuffer[0]++;
		}
		else{
			double prevProbSpace = probSpacing[0], nextProbSpace = prevProbSpace + probSpacing[1];
		    for( j=0; j<wf->m; j++){
			    if( prevProbSpace < rand  &&  rand <= nextProbSpace ){
				    qubitBuffer[j]++;
			    }
				/* Move to next spacing */
				prevProbSpace = nextProbSpace;
				nextProbSpace += probSpacing[j+1];
		    }
		}
	}
	/* Convert from probability to amplitude (ignoring phase) */
	for( i=0; i<wf->m; i++ ){
		evolvedQubit[i] = sqrt(qubitBuffer[i]);
	}
	
	int probEnd = 0;
	for( i=0; i<wf->m; i++ ){
		probEnd += calc_probability( wf->element[i] ); 
	}
	probEnd = round_nearest_int( probEnd );
	
	if( probEnd < probStart    ||    probStart < probEnd ){
		errno = -2;
		fprintf( stderr, "An error has occured (-2: Particle number not conserved, probStart = %i, probEnd = %i).\n", probStart, probEnd);
		exit(errno);
	}
}


/* Ammends a file adding the probability of the states to a new line */
static void save_wf( FILE *f, int m, double waveFunction[m] ){
	
	int i;
	
	fprintf( f, "\n");
	for( i=0; i<m; i++ ){
		fprintf( f, "%lg ", calc_probability( waveFunction[i] ) );
	}
}


/* Allocates memory for a sparse matrix with nMax non zero elements */
static sparseMatrix** sparse_matrix_alloc( int nMax ){
	
	sparseMatrix **matrix = malloc( nMax * (sizeof(sparseMatrix*)) );
	if( matrix == NULL ){
		fprintf( stderr, "Error allocating memory %i.\n", errno );
		exit(errno);
	}
	
	int i;
	for( i=0; i<nMax; i++ ){
		matrix[i] = malloc( 3 * sizeof(double complex) );    //safest to allocate memory for largest element of union
		if( matrix[i] == NULL ){
			fprintf( stderr, "Error allocating memory %i.\n", errno );
			exit(errno);
		}
		matrix[i][0].intVal = 0;
		matrix[i][1].intVal = 0;
		matrix[i][2].cDoubleVal = 0.0;
	}
	return matrix;
}


/* Multiplies the wavefunction by the sparse matrix operator */
static void matrix_operate( sparseMatrix **operator, qubits *wf, int nMax ){
	
	int i, j, k;
	double originalWf[ wf->m ];    //Need to make a copy of the original wave-function
	
	for( i=0; i<wf->m; i++ ){
		originalWf[i] = wf->element[i];
		wf->element[i] = 0.0 + 0.0*I;
	}
	
	for( i=0; i<wf->m; i++ ){
		for( j=0; j<wf->m; j++ ){
			for( k=0; k<nMax; k++ ){
				if( operator[k][0].intVal == i    &&    operator[k][1].intVal == j ){
					wf->element[i] += operator[k][2].cDoubleVal * creal( originalWf[j] )    +    operator[k][2].cDoubleVal * cimag( originalWf[j] );
				}
			}
		}
	}
}


/* Frees a sparse matrix which has been created by sparse_matrix_alloc */
static void sparse_matrix_free( sparseMatrix **matrix, int nMax ){
	
	int i;
	for( i=0; i<nMax; i++ ){
		free(matrix[i]);
	}
	free(matrix);
}


/* Scales a 2x2 gate up to operate on a single qubit in a system of n qubits, NOTE: qubit = 0 operates on first qubit in set */
static void scale_operator_2x2( qubits *wf, int qubit, sparseMatrix **scaledOperator, double complex operator[2][2] ){
	
	int i, j, k, delta, element = 0;
	
	for( i=0; i<wf->m; i++ ){
		for( j=0; j<wf->m; j++ ){
			delta = 0;
			for( k=0; k<wf->n; k++ ){
				/* Kronecker delta for every qubit but the qubit which is being operated on */
				if( wf->matrix[i][k] == wf->matrix[j][k]    &&    k != qubit ){    
					delta++;
				}
			}
			/* All Kronecker deltas are 1 */
			if( delta == wf->n - 1 ){
				scaledOperator[element][0].intVal = i;
				scaledOperator[element][1].intVal = j;
				scaledOperator[element][2].cDoubleVal = operator[ wf->matrix[i][qubit] ][ wf->matrix[j][qubit] ];
				element++;
			}
		}
	}
}


/* Scales a 4x4 gate up to operate on a 2 qubits in a system of wf->n qubits, (qubit 1 controlles qubit2)*/
static void scale_operator_4x4( qubits *wf, int qubit1, int qubit2, sparseMatrix **scaledOperator, double complex operator[4][4] ){
	
	int i, j, k, element = 0, check;
	
	/* Constructing Hadamard gate for wf->n qubits */
	for( i=0; i<wf->m; i++ ){
		for( j=0; j<wf->m; j++ ){
			check = 0;
			for( k=0; k<wf->n; k++ ){
				if( wf->matrix[i][k] == wf->matrix[j][k]    &&    k != qubit1    &&    k != qubit2 ){
					check++;
				}
			}
			if( check == wf->n-2 ){
				int element1Decimal=0, element2Decimal=0;
				/* Converting elements from binary to decimal ( minus any qubits which are not being operated on )*/
				int m = 1;
				element1Decimal += wf->matrix[i][qubit1] * pow( 2, m );    //qubit which controlls other qubit must come first
				element2Decimal += wf->matrix[j][qubit1] * pow( 2, m );
				m--;
				element1Decimal += wf->matrix[i][qubit2] * pow( 2, m );
				element2Decimal += wf->matrix[j][qubit2] * pow( 2, m );
				
				if( operator[element1Decimal][element2Decimal] != 0 ){
					scaledOperator[element][0].intVal = i;
					scaledOperator[element][1].intVal = j;
					scaledOperator[element][2].cDoubleVal = operator[element1Decimal][element2Decimal];
					element++;
				}
			}
		}
	}
}


/* Fills a given array with elements of the 2x2 Hadamard gate */
static void fill_array_hadamard( double complex hadamard[2][2] ){
	
	hadamard[0][0] = 1.0 / sqrt(2.0);
	hadamard[0][1] = 1.0 / sqrt(2.0);
	hadamard[1][0] = 1.0 / sqrt(2.0);
	hadamard[1][1] = -1.0 / sqrt(2.0);
}


/* Fills a given array with elements of the 2x2 phase shift gate, changing the phase by theta */
static void fill_array_phase( double complex phaseShift[2][2], double theta ){
	
	phaseShift[0][0] = 1.0;
	phaseShift[0][1] = 0;
	phaseShift[1][0] = 0;
	phaseShift[1][1] = cos(theta) + sin(theta)*I;
}


/* Creates a Grover diffusion operator for n qubits */
static void fill_array_grover( int m, sparseMatrix **grover ){
	
	int i, j, element = 0;
	
	for( i=0; i<m; i++ ){
		for( j=0; j<m; j++ ){
			if( i == j ){
				if( i == 0 ){
					grover[element][0].intVal = i;
					grover[element][1].intVal = j;
					grover[element][2].cDoubleVal = -1.0;
					element++;
				}
				else{
					grover[element][0].intVal = i;
					grover[element][1].intVal = j;
					grover[element][2].cDoubleVal = 1.0;
					element++;
				}
			}
		}
	} 
}


/* Creates a quantum oracle matrix for n qubits */ 
static void fill_array_oracle( int m, sparseMatrix **oracle, int question ){
	
	int i, j, element = 0;
	
	for( i=0; i<m; i++ ){
		for( j=0; j<m; j++ ){
			if( i == j ){
				if( i == question ){
					oracle[element][0].intVal = i;
					oracle[element][1].intVal = j;
					oracle[element][2].cDoubleVal = -1.0;
					element++;
				}
				else{
					oracle[element][0].intVal = i;
					oracle[element][1].intVal = j;
					oracle[element][2].cDoubleVal = 1.0;
					element++;
				}
			}
		}
	} 
}


/* Fills a given array with elements of the 4x4 CNOT gate */
static void fill_array_CNOT( double complex CNOT[4][4] ){
	
	int i, j;
	for( i=0; i<4; i++ ){
		for( j=0; j<4; j++ ){
			CNOT[i][j] = 0.0;
		}
	}
	CNOT[0][0] = 1.0;
	CNOT[1][1] = 1.0;
	CNOT[2][3] = 1.0;
	CNOT[3][2] = 1.0;
}


/* Fills a given array with elements of the 4x4 controlled phase shift gate, changing the phase by theta */
static void fill_array_cPhase( double complex cPhase[4][4], double theta ){
	
	int i, j;
	for( i=0; i<4; i++ ){
		for( j=0; j<4; j++ ){
			cPhase[i][j] = 0.0;
		}
	}
	cPhase[2][2] = 1.0;
	cPhase[2][3] = 0;
	cPhase[3][2] = 0;
	cPhase[3][3] = cos(theta) + sin(theta)*I;
}


/* Performs an IQFT using Hadamard gates and cPhase gates */
static void inverse_quantum_fourier_transform( qubits *wf ){
	
	int nMaxLarge = 2 * wf->m, nMaxSmall = wf->m;
	double complex hadamard[2][2], cPhase1[4][4], cPhase2[4][4];
	
	fill_array_hadamard( hadamard );
	fill_array_cPhase( cPhase1, M_PI/2.0 );
	fill_array_cPhase( cPhase2, M_PI/4.0 );
	
	sparseMatrix **hOperator1 = sparse_matrix_alloc( nMaxLarge );
	scale_operator_2x2( wf, 0, hOperator1, hadamard );
	matrix_operate( hOperator1, wf, nMaxLarge );
	sparse_matrix_free( hOperator1, nMaxLarge );
	
	sparseMatrix **pOperator1 = sparse_matrix_alloc( nMaxSmall );
	scale_operator_4x4( wf, 0, 1, pOperator1, cPhase1 );
	matrix_operate( pOperator1, wf, nMaxSmall );
	sparse_matrix_free( pOperator1, nMaxSmall );
	
	sparseMatrix **pOperator2 = sparse_matrix_alloc( nMaxSmall );
	scale_operator_4x4( wf, 0, 2, pOperator2, cPhase2 );
	matrix_operate( pOperator2, wf, nMaxSmall );
	sparse_matrix_free( pOperator2, nMaxSmall );
	
	sparseMatrix **hOperator2 = sparse_matrix_alloc( nMaxLarge );
	scale_operator_2x2( wf, 1, hOperator2, hadamard );
	matrix_operate( hOperator2, wf, nMaxLarge );
	sparse_matrix_free( hOperator2, nMaxLarge );
	
	sparseMatrix **pOperator3 = sparse_matrix_alloc( nMaxSmall );
	scale_operator_4x4( wf, 1, 2, pOperator3, cPhase1 );
	matrix_operate( pOperator3, wf, nMaxSmall );
	sparse_matrix_free( pOperator3, nMaxSmall );
	
	sparseMatrix **hOperator3 = sparse_matrix_alloc( nMaxLarge );
	scale_operator_2x2( wf, 2, hOperator3, hadamard );
	matrix_operate( hOperator3, wf, nMaxLarge );
	sparse_matrix_free( hOperator3, nMaxLarge );
}


/* Function to implement programming project 1 from the reference given at the top of the program */
static void project_1( qubits *wf ){
	
	int evolveN;
	char file[F_NAME_MAX];
	
	scan_elements( wf );
	printf("How many times would you like to evolve the qubits?\n");
	scanf(" %i", &evolveN);
	printf("Enter a name for the output file. (With extension, max %i characters)\n", F_NAME_MAX);
	scanf(" %256s", file);		
	check_file( file );
	FILE *f = fopen(file, "a");
	
	srand((unsigned)time(NULL));
	printf("\nEvolving qubits...\n");
	int i;
	double evolvedQubit[wf->m];    //evolve_qubits only returns real values
	
	fprintf( f, "# Output file created by QuantumComputing #\n"
				"# data contains probabilities of individual states in a set of qubits, being evolved and measured through time #");
	for( i=0; i<evolveN; i++ ){
		evolve_qubits( wf, evolvedQubit );
		save_wf( f, wf->m, evolvedQubit );
	}
	fclose(f);
}


/* Function to implement programming project 2 from the reference given at the top of the program */
static void project_2( qubits *wf ){
	
	int nMax = 2 * wf->m;    //Max number of non zero values in sparse matrices.
	
	double complex hadamard[2][2], phaseShift[2][2];
	fill_array_hadamard( hadamard ); 
	fill_array_phase( phaseShift, M_PI ); 
	
	/* Part a */
	printf("Applying a Hadamard gate to qubit 2:\n");
	sparseMatrix **hOperatorA = sparse_matrix_alloc( nMax );
	scale_operator_2x2( wf, 1, hOperatorA, hadamard );
	matrix_operate( hOperatorA, wf, nMax );
	sparse_matrix_free( hOperatorA, nMax );
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
	
	/* Part b */
	printf("Applying a Hadamard gates to all qubits:\n");
	
	int i;
	for( i=0; i<wf->n; i++ ){
		sparseMatrix **hOperatorB = sparse_matrix_alloc( nMax );
		scale_operator_2x2( wf, i, hOperatorB, hadamard );
		matrix_operate( hOperatorB, wf, nMax );
		sparse_matrix_free( hOperatorB, nMax );
	}
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
	
	/* Part c */
	printf("Applying two Hadamard gates to qubit 3:\n");
	sparseMatrix **hOperatorC = sparse_matrix_alloc( nMax );
	scale_operator_2x2( wf, 2, hOperatorC, hadamard );
	matrix_operate( hOperatorC, wf, nMax );
	matrix_operate( hOperatorC, wf, nMax );
	sparse_matrix_free( hOperatorC, nMax );
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
	
	/* Part d */
	printf("Applying a Hadamard gate to qubit 3, followed by a phase shift of pi, followed by another Hadamard gate:\n");
	
	sparseMatrix **hOperatorD1 = sparse_matrix_alloc( nMax );
	scale_operator_2x2( wf, 2, hOperatorD1, hadamard );
	matrix_operate( hOperatorD1, wf, nMax );
	sparse_matrix_free( hOperatorD1, nMax );
	
	sparseMatrix **pOperatorD = sparse_matrix_alloc( nMax );
	scale_operator_2x2( wf, 2, pOperatorD, phaseShift );
	matrix_operate( pOperatorD, wf, nMax );
	sparse_matrix_free( pOperatorD, nMax );
	
	sparseMatrix **hOperatorD2 = sparse_matrix_alloc( nMax );
	scale_operator_2x2( wf, 2, hOperatorD2, hadamard );
	matrix_operate( hOperatorD2, wf, nMax );
	sparse_matrix_free( hOperatorD2, nMax );
	
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
}


/* Performs Grover's quantum search */
static void grover_quantum_search( qubits *wf ){
	
	int nMax = 2 * wf->m;    //Max number of non zero values in sparse matrices.
	
	/* Creating oracle and Grover matrices */
	sparseMatrix **oracle = sparse_matrix_alloc( wf->m ), **grover = sparse_matrix_alloc( wf->m );    //Do not need to use nMax for these operators.
			
	fill_array_grover( wf->m, grover );
	printf("Enter a question for the oracle (0 <= question < %i).\n", wf->m );
	int question;
	scanf(" %i", &question);
		
	/* Allocating operators and putting qubits into superposition of all states */
	printf("Creating operators...\n");    //This may take a long time so tell user what program is doing.
	double complex hadamard[2][2];
	sparseMatrix ***hOperator = malloc( wf->n * sizeof(sparseMatrix**) );    //Need one operator for each qubit.
	
	fill_array_oracle( wf->m, oracle, question );
	fill_array_hadamard( hadamard ); 
	
	int i;
	for( i=0; i<wf->n; i++ ){
		hOperator[i] = sparse_matrix_alloc( nMax );
		scale_operator_2x2( wf, i, hOperator[i], hadamard );
		matrix_operate( hOperator[i], wf, nMax );
	}
			
	int optRep = round_nearest_int( (M_PI/4.0) * sqrt( pow(2.0, wf->n) ) );    //Optimum number of repetitions
	
	printf("Program will perform %i itterations.\n", optRep);
			
	int j;
	printf("Performing calculation...\n");
	for( i=0; i<optRep; i++ ){
		matrix_operate( oracle, wf, wf->m );
		for( j=0; j<wf->n; j++ ){
			matrix_operate( hOperator[j], wf, nMax );
		}
		matrix_operate( grover, wf, wf->m );
		for( j=0; j<wf->n; j++ ){
			matrix_operate( hOperator[j], wf, nMax );
		}
	}
	/* Free operators */
	
	for( i=0; i<wf->n; i++ ){
		sparse_matrix_free( hOperator[i], nMax );
	}
	free(hOperator);
	sparse_matrix_free( oracle, wf->m );
	sparse_matrix_free( grover, wf->m );
	
	/* Finding element with greatest probability aplitude */
	int elementMax = 0;
	for( i=1; i<wf->m; i++ ){
		if( cabs( wf->element[i] ) > cabs( wf->element[elementMax] ) ){
			elementMax = i;
		}
	}
			
	printf("The final wavefunction is:\n");
	display_wave_function( wf );
	printf("The question was %i with a probability of %g%%.\n", question, calc_probability(wf->element[question])*100 );
}


/* Function to implement programming project 4 from the reference given at the top of the program */
static void project_4( qubits *wf ){
	
	double complex CNOT[4][4], hadamard[2][2];
	fill_array_hadamard( hadamard );
	fill_array_CNOT( CNOT );
	int nMaxLarge = 2 * wf->m, nMaxSmall = wf->m;
	
	/* Part a */
	printf("\nApplying a Hadamard gate to qubit 2, followed by a CNOT gate with qubit 2 controlling qubit 3:\n");
	
	sparseMatrix **CNOTOperator23A = sparse_matrix_alloc( nMaxSmall ), 
				 **hOperatorA = sparse_matrix_alloc( nMaxLarge );
	
	scale_operator_2x2( wf, 1, hOperatorA, hadamard );
	matrix_operate( hOperatorA, wf, nMaxLarge);
	
	scale_operator_4x4( wf, 1, 2, CNOTOperator23A, CNOT );
	matrix_operate( CNOTOperator23A, wf, nMaxSmall);
	
	sparse_matrix_free( hOperatorA, nMaxLarge );
	sparse_matrix_free( CNOTOperator23A, nMaxSmall );
	
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
	
	/* Part b */
	printf("\nApplying a Hadamard gate to qubit 2, followed by a CNOT gate with qubit 2 controlling qubit 3, "
		   "followed by a CNOT gate with qubit 2 controlling qubit 1 (producing a cat state if number of qubits = 3):\n");
		   
	sparseMatrix **CNOTOperator23B = sparse_matrix_alloc( nMaxSmall ), 
				 **CNOTOperator21B = sparse_matrix_alloc( nMaxSmall ),
				 **hOperatorB = sparse_matrix_alloc( nMaxLarge );
				 
	scale_operator_2x2( wf, 1, hOperatorB, hadamard );
	matrix_operate( hOperatorB, wf, nMaxLarge);
	
	scale_operator_4x4( wf, 1, 2, CNOTOperator23B, CNOT );
	matrix_operate( CNOTOperator23B, wf, nMaxSmall);
	
	scale_operator_4x4( wf, 1, 0, CNOTOperator21B, CNOT );
	matrix_operate( CNOTOperator21B, wf, nMaxSmall);
	
	sparse_matrix_free( hOperatorB, nMaxLarge );
	sparse_matrix_free( CNOTOperator23B, nMaxSmall );
	sparse_matrix_free( CNOTOperator21B, nMaxSmall );
	
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
	
	/* Part c */
	printf("\nApplying 2 Hadamard gates to qubit 2:\n");
	
	sparseMatrix **hOperatorC = sparse_matrix_alloc( nMaxLarge );
	
	scale_operator_2x2( wf, 1, hOperatorC, hadamard );
	matrix_operate( hOperatorC, wf, nMaxLarge);
	matrix_operate( hOperatorC, wf, nMaxLarge);
	
	sparse_matrix_free( hOperatorC, nMaxLarge );
	
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
	
	/* Part d */
	printf("\nApplying a Hadamard gate to qubit 2, followed by a CNOT gate with qubit 2 controlling qubit 3, followed by a Hadamard gate:\n");
	
	sparseMatrix **hOperatorD = sparse_matrix_alloc( nMaxLarge ),
				 **CNOTOperator23D = sparse_matrix_alloc( nMaxSmall );
	
	scale_operator_2x2( wf, 1, hOperatorD, hadamard );
	matrix_operate( hOperatorD, wf, nMaxLarge);
	
	scale_operator_4x4( wf, 1, 2, CNOTOperator23D, CNOT );
	matrix_operate( CNOTOperator23D, wf, nMaxSmall);
	
	/* Can unobserve qubit 2 by applying another CNOT gate, meaning that the wf does not collapse **
	** scale_operator_4x4( wf, 1, 2, CNOTOperator23D, CNOT ); **
	** matrix_operate( CNOTOperator23D, wf, nMaxSmall); */
	
	matrix_operate( hOperatorD, wf, nMaxLarge);
	
	sparse_matrix_free( hOperatorD, nMaxLarge );
	sparse_matrix_free( CNOTOperator23D, nMaxSmall );
	
	printf("The wavefunction is now:\n");
	display_wave_function( wf );
	set_wf_gs( wf );
}


/* Performs Shor's quantum factoring algorithm for 7 qubits INCOMPLETE */
static void shors_algorithm( qubits *wf ){
	
	int L = 3, nMax = 2 * wf->m;
	double complex hadamard[2][2];
	fill_array_hadamard( hadamard );
	
	int i;
	for( i=0; i<L; i++ ){
		sparseMatrix **hOperator = sparse_matrix_alloc( nMax );
		scale_operator_2x2( wf, i, hOperator, hadamard );
		matrix_operate( hOperator, wf, nMax );
		sparse_matrix_free( hOperator, nMax );
	}
	//Insert additional gates here
	inverse_quantum_fourier_transform( wf );
	display_wave_function( wf );
	
}


/* Menu to allow user to select operations to perform */
static void quantum_computing_menu( qubits *wf ){
	
	int choice = 0, again = 1;
	while( again == 1){
		while( choice == 0 ){
			printf("\nChoose the task which you would like the program to perform:\n"
				   "1. Programming project 1: Simulate measurement of the N-qubit register.\n"
				   "2. Programming project 2: First full quantum computations.\n"
				   "3. Programming project 3: Perform Grover's quantum search.\n"
				   "4. Programming project 4: Computations using CNOT gates.\n"
				   "5. Programming project 7 : Shor's algorithm (INCOMPLETE).\n"
				   "6. Exit.\n");
			scanf(" %i", &choice);
			if( choice != 1    &&    choice != 2    &&    choice != 3    &&    choice != 4
				&& choice != 5    &&    choice != 6 ){
					
				printf("Invalid selection, please select from the options given.\n");
			}
		}
		if( choice == 1 ){
			project_1( wf );
		}
		else if( choice == 2 ){
			project_2( wf );
		}
		else if( choice == 3 ){
			grover_quantum_search( wf );
		}
		else if( choice == 4 ){
			project_4( wf );
		}
		else if( choice == 5 ){
			printf("WARNING Shor's algorithm is incomplete.\n");
			shors_algorithm( wf );
		}
		else if ( choice == 6 ){
			again = 0;
			break;
		}
		printf("Would you like to select another operation? (y/n)\n");
		
		/* Taking yes/no input */
		int valid = 0;	//Loop always starts
		while( valid == 0 ){
			char ans;
			scanf(" %c", &ans);
			if( ans == 'y' ){
				again = 1;
				break;
			}
			else if( ans == 'n' ){
				again = 0;
				break;
			}
			printf("Invalid selection, please use y/n.\n");
		}
		printf("\n");
		choice = 0;
	}
}



int main(){
	
	/* Allocating qubits structure */
	qubits *wf = malloc( sizeof(qubits) );
	if( wf == NULL ){
		fprintf(stderr, "Error allocating memory (%i).\n", errno);
		exit(errno);
	}
	initialise_qubits( wf );
	quantum_computing_menu(  wf );
	
	/* Freeing qubits structure */
	int i;
	for( i=0; i<wf->n; i++ ){
		free(wf->matrix[i]);
	}
	free(wf->matrix);
	free(wf->element);
	free(wf);
}