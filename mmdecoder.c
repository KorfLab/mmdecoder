/* mmdecoder.c */
//3rd version
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include "toolbox.h"

/* Globals */
#define STRAINS 55
#define TEST_STRAINS 1
#define CHROMS 22
double PR_TRANS = 0.02;
double PR_EMIT = 0.0001;
char * DATA_FILE = "mmdecoder.db";
char * TEST_FILE = "TEST.db";
//for database cases
ik_tvec STRvec;               /* array of 55 strain names */
ik_map  STRmap;               /* map of 55 names to numbers */
ik_tvec CHRvec;               /* array of 22 chrom names */
ik_map  CHRmap;               /* map of 22 names to numbers */
char * SEQ[STRAINS][CHROMS];  /* [strain][chrom] = SNP sequence */
//for TEST cases
ik_tvec TEST_STRvec;               /* array of 55 strain names */
ik_map  TEST_STRmap;               /* map of 55 names to numbers */
ik_tvec TEST_CHRvec;               /* array of 22 chrom names */
ik_map  TEST_CHRmap;               /* map of 22 names to numbers */
char * TEST_SEQ[STRAINS][CHROMS];  /* [strain][chrom] = SNP sequence */

//read database file
static void read_database (const char * database_filename) {
	
	ik_pipe databasefile;
	char line[20000], strain[64], chrom[64], seq[20000];
	void *v;
	int *ir, c, s, chr_count = 0, str_count = 0;
	
	// check letters, only {A, C, G, T, H, N}
	
	//initiation
	STRvec = ik_tvec_new();
	STRmap = ik_map_new();
	CHRvec = ik_tvec_new();
	CHRmap = ik_map_new();
	databasefile = ik_pipe_open(database_filename, "r");
	fprintf(stderr, "reading data file %s\n", database_filename);
	
	while (fgets(line, sizeof(line), databasefile->stream) != NULL) {
		if (sscanf(line, "%s %s %s", strain, chrom, seq) != 3) {
			ik_exit(1, "parsing error\n");
		}
		
		v = ik_map_get(CHRmap, chrom);
		if (v == NULL) {
			ir = ik_malloc(sizeof(int));
			*ir = chr_count;
			chr_count++;
			ik_map_set(CHRmap, chrom, ir);
			ik_tvec_push(CHRvec, chrom);
			c = *ir;
		} else {
			ir = v;
			c = *ir;
		}
		
		v = ik_map_get(STRmap, strain);
		if (v == NULL) {
			ir = ik_malloc(sizeof(int));
			*ir = str_count;
			str_count++;
			ik_map_set(STRmap, strain, ir);
			ik_tvec_push(STRvec, strain);
			s = *ir;
		} else {
			ir = v;
			s = *ir;
		}
		
		SEQ[s][c] = ik_malloc(strlen(seq)+1);
		strcpy(SEQ[s][c], seq);
	}
}

//read TEST file
static void read_TESTdata (const char * TEST_filename) {
	
	ik_pipe TESTfile;
	char line[20000], strain[64], chrom[64], seq[20000];
	void *v;
	int *ir, c, s, chr_count = 0, str_count = 0;
	
	// check letters, only {A, C, G, T, H, N}
	
	//initiation
	TEST_STRvec = ik_tvec_new();
	TEST_STRmap = ik_map_new();
	TEST_CHRvec = ik_tvec_new();
	TEST_CHRmap = ik_map_new();
	TESTfile = ik_pipe_open(TEST_filename, "r");
	fprintf(stderr, "reading data file %s\n", TEST_filename);
	
	while (fgets(line, sizeof(line), TESTfile->stream) != NULL) {
		if (sscanf(line, "%s %s %s", strain, chrom, seq) != 3) {
			ik_exit(1, "parsing error\n");
		}
		
		v = ik_map_get(TEST_CHRmap, chrom);
		if (v == NULL) {
			ir = ik_malloc(sizeof(int));
			*ir = chr_count;
			chr_count++;
			ik_map_set(TEST_CHRmap, chrom, ir);
			ik_tvec_push(TEST_CHRvec, chrom);
			c = *ir;
		} else {
			ir = v;
			c = *ir;
		}
		
		v = ik_map_get(TEST_STRmap, strain);
		if (v == NULL) {
			ir = ik_malloc(sizeof(int));
			*ir = str_count;
			str_count++;
			ik_map_set(TEST_STRmap, strain, ir);
			ik_tvec_push(TEST_STRvec, strain);
			s = *ir;
		} else {
			ir = v;
			s = *ir;
		}
		
		TEST_SEQ[s][c] = ik_malloc(strlen(seq)+1);
		strcpy(TEST_SEQ[s][c], seq);
	}
}

void fasta (void) {

	int s, c;
	
	printf("size %d\n", STRvec->size);
	
	for (s = 0; s < STRAINS; s++) {
		for (c = 0; c < CHROMS; c++) {
			printf(">%s ", STRvec->elem[s]);
			printf("-Chr%s\n", CHRvec->elem[c]);
			printf("%s\n", SEQ[s][c]);
		}
	}
}

void dist (void) {
	
	//s1 is the current strain of database, s2 is the current strain of test
	//If want to calculate dist for specific chromosome, need to initiate c, 0 for chrom M, 1 for chrom X, 2 for chrom Y, 3 for chrom 1 and so on.
	int s1, s2, c, base;
	double score = 0.0;
	
	/*uncomment if need to calculate dist for specific chromosome
	printf ("for chr %s \n", CHRvec->elem[c]); */
	
	printf("\t%s\t", STRvec->elem[0]);
	
	for (s1 = 0; s1 < STRAINS; s1++) {
		for (s2 = s1; s2 < STRAINS; s2++) {
			for (c = 0; c < CHROMS; c++) {
				for (base = 0; SEQ[s1][c][base] != '\0'; base++) {
			
					if (SEQ[s1][c][base] == 'H' || SEQ[s2][c][base] == 'H') score += 0.5;
				
					else if (SEQ[s1][c][base] == 'N' || SEQ[s2][c][base] == 'N') score += 0.75;
				
					else if (SEQ[s1][c][base] != SEQ[s2][c][base]) score += 1;
					
				}
			}
			
			printf("%.2lf ", score);
			score = 0;
		}
		
		printf("\n");
		printf("\t%s\t", STRvec->elem[s1 + 1]);
	}
}

void viterbi (int c, int s2) {
	//s2 is the current strain of test, c is the current chromosome
	
	//t indicates the base are reaching, s1 is the current strain of database
	//size is the size of current chromosome
	//pb is the value of highest probability reaching the current base.
	int i, j, base, s1, size, path_size, trace_size;
	double score, max_score, pb = 1.0;
	int path[11000][56];
	int trace[11000][56];
	
	///* 2, compute and get path way matrix *///
			
	for (base = 0; TEST_SEQ[s2][c][base] != '\0'; base++) {
		
		max_score = 0.0;
		path_size = 0;
					
		for (s1 = 0; s1 < STRAINS; s1++) {
					
			//when we have the same base matching
			if (TEST_SEQ[s2][c][base] == SEQ[s1][c][base]) {
				
				/*add the probability of same base matching to the score
				there is a probability of having correct base when there might be an error existing.
				since we need consider both the test and database, so I multiply it twice */
				score = pb * (1 - PR_EMIT) * (1 - PR_EMIT);
				
				if (base != 0) {
					//suppose there is recombination happened
					score = score * (PR_TRANS / 55);
					
					for (i = 0; path[base - 1][i] != -1; i++) {
						/*when the base of current database comes from the same strains as the previous one.
						i.e. there might be no recombination happen
						offset the pervious calculation and do new calculation based on no recombination*/
						if (s1 == path[base - 1][i]) {
							score = score / (PR_TRANS / 55) * (1 - PR_TRANS);
							break;
						}
					}
				}	
			}
			//when we have the different base matching
			else {
				score = pb * (PR_EMIT) * (PR_EMIT);
				
				if (base != 0) {
					//suppose there is recombination happened
					score = score * (PR_TRANS / 55);
					
					for (i = 0; path[base - 1][i] != -1; i++) {
						/*when the base of current database comes from the same strains as the previous one.
						i.e. there might be no recombination happen
						offset the pervious calculation and do new calculation based on no recombination*/
						if (s1 == path[base - 1][i]) {
							score = score / (PR_TRANS / 55) * (1 - PR_TRANS);
							break;
						}
					}
				}
			}
					
			//if we have a new max score, i.e. a new possible origin strain
			if (score >= max_score) {
				
				//store the potential into the list
				path[base][path_size] = s1;
				//printf("store %d in base %d as %d\n", s1, base, path_size);
						
				max_score = score;
				path_size++;
			}
		}

		path[base][path_size] = -1;		//-1 indicates the ending of each row
		pb = max_score;
		
	}
	size = base;	
	
	///* 3, trace back *///
	/*starting with the last but one base, excluding those strains only show in the front 
	since we consider that recombination are less likely to happen and a certain length of bases are more likely to origin from the same strain*/
	for (i = 0; path[base - 1][i] != -1; i++) {
		trace[base - 1][i] = path[base - 1][i];
	}
	trace[base - 1][i] = -1;
	
	for (base = base - 2; base >= 0; base--) {
		trace_size = 0;
		
		for (i = 0; path[base][i] != -1; i++) {
			for (j = 0; trace[base + 1][j] != -1; j++) {
				
				if (path[base][i] == trace[base + 1][j]) {
					trace[base][trace_size] = path[base][i];
					trace_size++;
					break;
				}
			}
		}
		trace[base][trace_size] = -1;
		
		//recombination happened
		if (trace_size == 0) {
			for (i = 0; path[base][i] != -1; i++) {
				trace[base][i] = path[base][i];
			}
			trace[base][i] = -1;
		}
	}
	
	///* 4, report *///
	
	for (i = 0; i < 10; i++) {
		for (j = 0; trace[i][j] != -1; j++) {
			printf("%d, ", trace[i][j]);
		}
		printf("\n");
	}
	
	printf("size: %d\n", size);
	printf("The probability of having this outcome is: %lf\n\n", pb);
	
	///* 5, termination *///
	
}

void origin (void) {

	//s2 is the current strain of test, c is the current chromosome
	int s2, c;
	
	for (s2 = 0; s2 < TEST_STRAINS; s2++) {
		for (c = 0; c < CHROMS; c++) {
			printf("for strain %s -> chromosome %s: \n", TEST_STRvec->elem[s2], TEST_CHRvec->elem[c]);
			
			viterbi (c, s2);
		}
	}
}

static char * usage = "\
usage: mmdecoder [paramters]\n\
mandatory:\n\
  -d <file>  database file, mmdecoder.db if not specified\n\
commands:\n\
  -decode      run viterbi and posterior decoding algorithms\n\
    -f <file>  file of variant data\n\
    -x <float> exit probability, 0.0001 if not specified]\n\
    -y <float> error probability, 0.0001 if not specified]\n\
  -tree        create a distance tree in NEXUS format\n\
    -f <file>  optional additional variant file\n\
  -fasta       export sequences as fasta format\n\
  -dist	       export the distances matrix among strings for different chromosomes\n\
  -origin	   calculate the origin of TEST/new strings\n\
";

int main (int argc, char ** argv) {

	if (argc == 1) {
		fprintf(stderr, "%s", usage);
		exit(1);
	}
	
	/* command line options */
	ik_set_program_name(argv[0]);
	ik_register_option("-d", 1);
	ik_register_option("-f", 1);
	ik_register_option("-x", 1);
	ik_register_option("-y", 1);
	ik_register_option("-decode", 0);
	ik_register_option("-fasta", 0);
	ik_register_option("-tree", 0);
	ik_register_option("-dist", 0);
	ik_register_option("-origin", 0);
	ik_parse_options(&argc, argv);
	
	if (ik_option("-x")) PR_TRANS = atof(ik_option("-x"));
	if (ik_option("-y")) PR_EMIT = atof(ik_option("-y"));
	if (ik_option("-d")) DATA_FILE = ik_option("-d");
	
	read_database(DATA_FILE);
	read_TESTdata(TEST_FILE);
	
	if (ik_option("-fasta")) {
		fasta();
	}
	
	if (ik_option("-dist")) {
		dist();
	}
	
	if (ik_option("-origin")) {
		origin();
	}
	
	printf("exit: %f\n", PR_TRANS);
	printf("emit: %f\n", PR_EMIT);
	
	return 0;
}
