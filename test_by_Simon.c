//test.c by Xiuge Chen 

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

//Globals
#define STRAINS 4
#define CHROMS 3
double PRB_RECOM = 0.01;
double PRB_WRONG = 0.001;
double PRB_SAME_BASE = 0.99;
int LENGTH = 10;

char * SEQ[STRAINS][CHROMS];
char * TEST[CHROMS];



void ComputeMatrix ();
void Initiation ();

struct ScoreMatrix {
	int * strain;
	int * prestrain;
};

int main (void) {
	
	
	Initiation ();
	
	//test if the initiation is correct
	
	/*
	int i = 0, j = 0;
	for (i = 0; i < STRAINS; i++) {
		for (j = 0; j < CHROMS; j++) {
			printf("%s\n", SEQ[i][j]);
		}
	}
	
	for (j = 0; j < CHROMS; j++) {
			printf("%s\n", TEST[j]);
	}
	*/
	
	ComputeMatrix ();
	
	return 0;
}

//function for artificial initiation
void Initiation () {
	
	int i = 0, j = 0;
	
	for (i = 0; i < STRAINS; i++) {
		for (j = 0; j < CHROMS; j++) {
			SEQ[i][j] = (char*)malloc(LENGTH * sizeof(char));
		}
	}
	
	SEQ[0][0] = "TTTTTTTTTT";
	SEQ[0][1] = "GCGCGCGCGC";
	SEQ[0][2] = "ATCGATCGAT";
	
	SEQ[1][0] = "AAAAAAAAAA";
	SEQ[1][1] = "ATATATATAT";
	SEQ[1][2] = "GCTTGCTTGC";
	
	SEQ[2][0] = "GGGGGGGGGG";
	SEQ[2][1] = "AGAGAGAGAG";
	SEQ[2][2] = "TTAATTAAT";
	
	SEQ[3][0] = "CCCCCCCCCC";
	SEQ[3][1] = "TCTCTCTCTC";
	SEQ[3][2] = "CGTTCGTTCG";
	
	
	for (i = 0; i < CHROMS; i++) {
		TEST[i] = (char*)malloc(LENGTH * sizeof(char));
	}
	
	TEST[0] = "AATTAAGGCC";
	TEST[1] = "GCGCAGATAT";
	TEST[2] = "ATCGCGTAAT";
	
	return;
}

void ComputeMatrix () {

	int numChr = 0, position = 0, numStrain = 0, numMaxScore = 0, i = 0;
	double score = 0.0, maxScore = 0.0;
	int * strainMaxScore[CHROMS];
	int numBases = 0;
	struct ScoreMatrix * matrix[CHROMS];
	
	for (numChr = 0; numChr < CHROMS; numChr++) {
	    
	    printf("Chromosome %d:\n", numChr + 1);
		
		numBases = strlen(SEQ[0][numChr]);
		
		printf("numBases: %d\n",numBases);
		
		matrix[numChr] = (struct ScoreMatrix *)malloc(numBases * (sizeof(int*) * 2));
		
		for (position = 0; position < numBases; position++) {
			
			matrix[numChr][position].strain = NULL;
			maxScore = 0.0;
			
			//count the number of strains have the same max score
			for (numStrain = 0; numStrain < STRAINS; numStrain++) {
				
				//probability of the same base
				if (TEST[numChr][position] == SEQ[numStrain][numChr][position]) {
					score = PRB_SAME_BASE;
				}
				else {
					score = 1 - PRB_SAME_BASE;
				}
				
				//probability of recombination
				if (position != 0) {
				    for (i = 0; i < numMaxScore; i++) {
				    	if (numStrain == matrix[numChr][position].prestrain[i]) {
					    	score *= 1 - PRB_RECOM;
				    	}
				    	else {
					    	score *= PRB_RECOM;
				    	}
			    	}
				}
				
				//probability of wrong base
				score *= 1 - PRB_WRONG;
	
				//count the number of strains have the same max score
				//if have more than one, build a array to store it
				if ((score - maxScore) > 1e-10) {
				
					maxScore = score;
					numMaxScore = 1;
					
					matrix[numChr][position].strain = (int*)malloc(sizeof(int));
					*matrix[numChr][position].strain = numStrain;
				}
				else if (fabs(score - maxScore) < 1e-10) {
					numMaxScore++;
					
					strainMaxScore[numChr] = (int*)malloc(numMaxScore * sizeof(int));
					
					for (i = 0; i < numMaxScore - 1; i++) {
						strainMaxScore[numChr][i] = matrix[numChr][position].strain[i];
					}
					strainMaxScore[numChr][i] = numStrain;
					
					free(matrix[numChr][position].strain);
					matrix[numChr][position].strain = (int*)malloc(numMaxScore * sizeof(int));
					
					for (i = 0; i < numMaxScore; i++) {
						matrix[numChr][position].strain[i] = strainMaxScore[numChr][i];
					}
				}	
			}
			
			//print out the order of the strain that has the highest probability
			for (i = 0; i < numMaxScore;) {
				printf("%d", matrix[numChr][position].strain[i] + 1);
				i++;
				if (i != numMaxScore) {
				    printf(" or ");
				}
			}
			
			printf(", ");
			
			//assign the order of the strain that has the highest probability to the order of the prestrain variable in next position.
			for (i = 0; i < numMaxScore && position < numBases; i++) {
			    matrix[numChr][position + 1].prestrain = (int*)malloc(numMaxScore * sizeof(int));
		        matrix[numChr][position + 1].prestrain[i] = matrix[numChr][position].strain[i];
	    	}
		}
		
		printf("\n");
	}
	
	return;
}