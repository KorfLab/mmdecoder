/* mmdecoder.c 2nd edition*/
//add the -fasta and -origin function, but the -origin function is still need to be fixed
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#include "toolbox.h"

/* Globals */
#define STRAINS 55
#define CHROMS 22
double EXIT_PRB = 0.0001;
double MIN_EMIT = 0.0001;
double PRB_RECOM = 0.01;
double PRB_WRONG = 0.001;
double PRB_SAME_BASE = 0.999;
char * DATA_FILE = "mmdecoder.db";
//for database cases
ik_tvec STRvec;               /* array of 55 strain names */
ik_map  STRmap;               /* map of 55 names to numbers */
ik_tvec CHRvec;               /* array of 22 chrom names */
ik_map  CHRmap;               /* map of 22 names to numbers */
char * SEQ[STRAINS][CHROMS];  /* [strain][chrom] = SNP sequence */
//for test cases
ik_tvec TEST_STRvec;              /* array of test strain names */
ik_map  TEST_STRmap;              /* map of test names to numbers */
ik_tvec TEST_CHRvec;              /* array of test chrom names */
ik_map  TEST_CHRmap;              /* map of test names to numbers */
char * TEST[STRAINS][CHROMS];	  /* [strain][chrom] = TEST SNP sequence */

typedef struct OriginMatrix {
	int *** num_strain;
	int **** strain;
	int *** num_prestrain;
	int **** prestrain;
} OriginMatrix;

typedef struct TraceBackMatrix {
	int ** num_strain;
	int **** strain;
} TraceBackMatrix;

static void read_database (const char *filename) {
	ik_pipe datafile;
	char line[20000], strain[64], chrom[64], seq[20000];
	void *v;
	int *ir, c, s, chr_count = 0, str_count = 0;
	
	// check letters, only {A, C, G, T, H, N}
	
	STRvec = ik_tvec_new();
	STRmap = ik_map_new();
	CHRvec = ik_tvec_new();
	CHRmap = ik_map_new();
	
	datafile = ik_pipe_open(filename, "r");
	fprintf(stderr, "reading data file %s\n", filename);
	while (fgets(line, sizeof(line), datafile->stream) != NULL) {
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

static void read_test (char *filename, int *num_strings) {
	ik_pipe datafile;
	char line[20000], strain[64], chrom[64], seq[20000];
	void *v;
	int *ir, c, s, chr_count = 0, str_count = 0;
	
	// check letters, only {A, C, G, T, H, N}
	
	TEST_STRvec = ik_tvec_new();
	TEST_STRmap = ik_map_new();
	TEST_CHRvec = ik_tvec_new();
	TEST_CHRmap = ik_map_new();
	
	datafile = ik_pipe_open(filename, "r");
	fprintf(stderr, "reading data file %s\n", filename);
	while (fgets(line, sizeof(line), datafile->stream) != NULL) {
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
		
		TEST[s][c] = ik_malloc(strlen(seq)+1);
		strcpy(TEST[s][c], seq);
	}
	
	*num_strings = s + 1;
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

void Create_Matrix(OriginMatrix * matrix, int num_strings) {

	int test_string, chr, num_bases, base;

	matrix->num_strain = (int ***) malloc(num_strings * sizeof(int **));
	matrix->num_prestrain = (int ***) malloc(num_strings * sizeof(int **));
	matrix->strain = (int ****) malloc(num_strings * sizeof(int ***));
	matrix->prestrain = (int ****) malloc(num_strings * sizeof(int ***));
	
	for (test_string = 0; test_string < num_strings; test_string++) {
	
		matrix->num_strain[test_string] = (int **) malloc(CHROMS * sizeof(int *));
		matrix->num_prestrain[test_string] = (int **) malloc(CHROMS * sizeof(int *));
		matrix->strain[test_string] = (int ***) malloc(CHROMS * sizeof(int **));
		matrix->prestrain[test_string] = (int ***) malloc(CHROMS * sizeof(int **));
		
		for (chr = 0; chr < CHROMS; chr++) {
		
			num_bases = strlen(SEQ[0][chr]);
			
			matrix->num_strain[test_string][chr] = (int *) malloc(num_bases * sizeof(int));
			matrix->num_prestrain[test_string][chr] = (int *) malloc(num_bases * sizeof(int));
			matrix->strain[test_string][chr] = (int **) malloc(num_bases * sizeof(int *));
			matrix->prestrain[test_string][chr] = (int **) malloc(num_bases * sizeof(int *));
			
			for (base = 0; base < num_bases; base++) {
				
				matrix->strain[test_string][chr][base] = (int *) malloc(sizeof(int));
				matrix->prestrain[test_string][chr][base] = (int *) malloc(sizeof(int));
			}
		}
	}
}

void Destroy_Matrix(OriginMatrix * matrix, int num_strings) {

	int test_string, chr, num_bases, base;
	
	for (test_string = 0; test_string < num_strings; test_string++) {
		for (chr = 0; chr < CHROMS; chr++) {
		
			num_bases = strlen(SEQ[0][chr]);
			
			for (base = 0; base < num_bases; base++) {
				
				if (base) {
					free(matrix->prestrain[test_string][chr][base]);
				}
				free(matrix->strain[test_string][chr][base]);
			}
		}
	}
	
	for (test_string = 0; test_string < num_strings; test_string++) {
		for (chr = 0; chr < CHROMS; chr++) {
			free(matrix->num_strain[test_string][chr]);
			free(matrix->num_prestrain[test_string][chr]);
			free(matrix->strain[test_string][chr]);
			free(matrix->prestrain[test_string][chr]);
		}
	}

	for (test_string = 0; test_string < num_strings; test_string++) {
		free(matrix->num_strain[test_string]);
		free(matrix->num_prestrain[test_string]);
		free(matrix->strain[test_string]);
		free(matrix->prestrain[test_string]);
	} 

	free(matrix->num_strain);
	free(matrix->num_prestrain);
	free(matrix->strain);
	free(matrix->prestrain);
}

void Print_Results(TraceBackMatrix trace_back_matrix, int num_strings) {
	
	int test_string, chr, num_bases, base, path;
	
	for (test_string = 0; test_string < num_strings; test_string++) {
	
		for (chr = 0; chr < CHROMS; chr++) {
				
			for (path = 0; path < trace_back_matrix.num_strain[test_string][chr]; path++) {
			
				printf("%s %s ", TEST_STRvec->elem[test_string], CHRvec->elem[chr]);
					
				num_bases = strlen(SEQ[0][chr]);
			
				//printf("number of bases: %d\n", num_bases);
			
				for (base = num_bases - 1; base >= 0; base--) {
					
					if (base == num_bases - 1) {
						printf("[1-");
					}
					else if (base == 0) {
						printf("%d](%s)", num_bases, STRvec->elem[(trace_back_matrix.strain[test_string][chr][path][base]) - 1]);
					}
					else {
						if (trace_back_matrix.strain[test_string][chr][path][base - 1] != trace_back_matrix.strain[test_string][chr][path][base]) {
							printf("%d](%s)[%d-", num_bases - base, STRvec->elem[(trace_back_matrix.strain[test_string][chr][path][base]) - 1], num_bases - base + 1);
						}
					}
				}
			}
			
			printf("\n");
		}
	}
}

void Trace_back (OriginMatrix matrix, int num_strings) {

	int test_string, chr, num_bases, base;
	int i, j, k;
	int num_paths;
	int num_max_scores;
	double score, max_score;
	TraceBackMatrix trace_back_matrix;

	//allocate memory for trace_back_matrix
	trace_back_matrix.num_strain = (int **) malloc(num_strings * sizeof(int *));
	trace_back_matrix.strain = (int ****) malloc(num_strings * sizeof(int ***));
	
	for (test_string = 0; test_string < num_strings; test_string++) {
	
		trace_back_matrix.num_strain[test_string] = (int *) malloc(CHROMS * sizeof(int));
		trace_back_matrix.strain[test_string] = (int ***) malloc(CHROMS * sizeof(int **));
		
		for (chr = 0; chr < CHROMS; chr++) {
		
			num_bases = strlen(SEQ[0][chr]);
	
			trace_back_matrix.strain[test_string][chr] = (int **) malloc(sizeof(int *));
			
			trace_back_matrix.strain[test_string][chr][0] = (int *) malloc(num_bases * sizeof(int));
		}
	}
	
	//trace back
	for (test_string = 0; test_string < num_strings; test_string++) {
		for (chr = 0; chr < CHROMS; chr++) {
		
			num_bases = strlen(SEQ[0][chr]);
			num_paths = 1;
			
			for (base = num_bases - 1; base >= 0; base--) {
			
				max_score = 0.0;
				num_max_scores = 1;
				
				for (i = 0; i < matrix.num_strain[test_string][chr][base]; i++) {
					for (j = 0; j < matrix.num_prestrain[test_string][chr][base]; j++) {
					
						//calculate the score, only consider the probability of the recombination
						if (matrix.strain[test_string][chr][base][i] == matrix.prestrain[test_string][chr][base][j]) {
							score = 1 - PRB_RECOM;
						}
						else {
							score = PRB_RECOM;
						}
						
						if (score > max_score) {
						
							max_score = score;
							trace_back_matrix.strain[test_string][chr][num_max_scores - 1][base] = matrix.num_strain[test_string][chr][base];
						}
						else if (fabs(score - max_score) < 1e-10) {
						
							num_max_scores++;
							
							trace_back_matrix.strain[test_string][chr] = (int **)realloc(trace_back_matrix.strain[test_string][chr], num_max_scores * sizeof(int *));
							trace_back_matrix.strain[test_string][chr][num_max_scores - 1] = (int *)malloc(num_bases * sizeof(int));
							
							for (k = 0; k < base; k++) {
								trace_back_matrix.strain[test_string][chr][num_max_scores - 1][k] = trace_back_matrix.strain[test_string][chr][0][k];
							}
							trace_back_matrix.strain[test_string][chr][num_max_scores - 1][base] = matrix.num_strain[test_string][chr][base];
						}
					}
				}
			}
			
			trace_back_matrix.num_strain[test_string][chr] = num_paths;
		}
	}
	
	//print out the results
	Print_Results(trace_back_matrix, num_strings);
	/*for (test_string = 0; test_string < num_strings; test_string++) {
		for (chr = 0; chr < CHROMS; chr++) {
			for (i = 0; i < trace_back_matrix.num_strain[test_string][chr]; i++) {
			
				printf("%s %s ", TEST_STRvec->elem[test_string], CHRvec->elem[chr]);
				
				num_bases = strlen(SEQ[0][chr]);
			
				for (base = num_bases - 1; base >= 0; base--) {
					printf("%d", trace_back_matrix.strain[test_string][chr][i][base]);
				}
				
				printf("\n");
			}
		}
	}*/
	
	//free all the memory
	for (test_string = 0; test_string < num_strings; test_string++) {
		for (chr = 0; chr < CHROMS; chr++) {
			for (i = 0; i < trace_back_matrix.num_strain[test_string][chr]; i++) {
				free(trace_back_matrix.strain[test_string][chr][i]);
			}
		}
	}
	
	for (test_string = 0; test_string < num_strings; test_string++) {
		for (chr = 0; chr < CHROMS; chr++) {
			free(trace_back_matrix.strain[test_string][chr]);
		}
	}
	for (test_string = 0; test_string < num_strings; test_string++) {
		free(trace_back_matrix.num_strain[test_string]);
		free(trace_back_matrix.strain[test_string]);
	}
	free(trace_back_matrix.num_strain);
	free(trace_back_matrix.strain);
}

void ComputeOrigin (int num_strings) {

	int test_string, chr, num_bases, base, data_string;
	int i;
	double score, max_score;
	int num_max_scores;
	bool same_base;
	OriginMatrix matrix;
	
	Create_Matrix(&matrix, num_strings);

	for (test_string = 0; test_string < num_strings; test_string++) {
	
		//printf("test string: %d\n", test_string + 1);
		
		for (chr = 0; chr < CHROMS; chr++) {
		
			//printf("chr: %d\n", chr + 1);
		
			num_bases = strlen(SEQ[0][chr]);
			//printf("number of bases: %d\n", num_bases);
			
			for (base = 0; base < num_bases; base++) {
			
				//printf("base: %d\n", base + 1);
				max_score = 0;
				num_max_scores = 1;
			
				for (data_string = 0; data_string < STRAINS; data_string++) {
				
					//printf("comparing to data string: %d\n", data_string + 1);
					
					//calculate the score
					//probability of the same base
					if (TEST[test_string][chr][base] == SEQ[data_string][chr][base]) {
						score = PRB_SAME_BASE;
					}
					else {
						score = 1 - PRB_SAME_BASE;
					}
				
					//probability of recombination
					if (base) {
					
						same_base = false;
					
				    	for (i = 0; i < matrix.num_prestrain[test_string][chr][base]; i++) {
				    		if (data_string == matrix.prestrain[test_string][chr][base][i]) {
					    		same_base = true;
				    		}
			    		}
			    		
			    		if (same_base) {
			    			score *= 1 - PRB_RECOM;
			    		}
			    		else {
			    			score *= PRB_RECOM;
			    		}
					}
					else {
						matrix.num_prestrain[test_string][chr][base] = 0;
						free(matrix.prestrain[test_string][chr][base]);
						matrix.prestrain[test_string][chr][base] = NULL;
					}
				
					//probability of wrong base
					score *= 1 - PRB_WRONG;
					
					//printf("score: %lf \n", score);
					
					if (score > max_score) {
					
						max_score = score;
						matrix.strain[test_string][chr][base][num_max_scores - 1] = data_string;
					}
					else if (fabs(score - max_score) < 1e-10) {
					
						num_max_scores++;
						
						matrix.strain[test_string][chr][base] = (int *) realloc(matrix.strain[test_string][chr][base], num_max_scores * sizeof(int));
						
						matrix.strain[test_string][chr][base][num_max_scores - 1] = data_string;
					}
				}
				
				matrix.num_strain[test_string][chr][base] = num_max_scores;
				
				if (base != num_bases - 1) {
				
					matrix.prestrain[test_string][chr][base + 1] = (int *) realloc(matrix.prestrain[test_string][chr][base + 1], num_max_scores * sizeof(int));
				
					for (i = 0; i < num_max_scores; i++) { 
						matrix.prestrain[test_string][chr][base + 1][i] = matrix.strain[test_string][chr][base][i];
					}
					
					matrix.num_prestrain[test_string][chr][base + 1] = num_max_scores;
				}
				
			}
		}
	}
	
	/*for (test_string = 0; test_string < num_strings; test_string++) {
		printf("test string: %d\n", test_string + 1);
	
		for (chr = 0; chr < CHROMS; chr++) {
			printf("chr: %d\n", chr + 1);
			num_bases = strlen(SEQ[0][chr]);
			printf("number of bases: %d\n", num_bases);
			
			for (base = 0; base < num_bases; base++) {
				printf("base: %d\n", base + 1);
				printf("%d %d\n", matrix.num_strain[test_string][chr][base], matrix.num_prestrain[test_string][chr][base]);
			}
		}
	}*/
	
	//trace back and print out the results
	Trace_back (matrix, num_strings);
	
	Destroy_Matrix(&matrix, num_strings);
}

void origin (void) {

	char file_name[100];
	int name_length = 0, num_strings = 0;
	
	//get the name of the file
	printf("please input the name of the file containing the test string: ");
	for (name_length = 0; file_name[name_length - 1] != '\n'; name_length++) {
		scanf("%c", &file_name[name_length]);
	}
	file_name[name_length - 1] = '\0';
	
	//open the file and read test
	read_test(file_name, &num_strings);
	
	ComputeOrigin(num_strings);
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
  -origin	   annotate the genetic origins of a mouse straing in fasta format\n\
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
	ik_register_option("-origin", 0);
	ik_parse_options(&argc, argv);
	
	if (ik_option("-x")) EXIT_PRB = atof(ik_option("-x"));
	if (ik_option("-y")) MIN_EMIT = atof(ik_option("-y"));
	if (ik_option("-d")) DATA_FILE = ik_option("-d");
	
	read_database(DATA_FILE);
	
	if (ik_option("-fasta")) {
		fasta();
	}
	
	if (ik_option("-origin")) {
		origin();
	}
	
	printf("exit: %f\n", EXIT_PRB);
	printf("emit: %f\n", MIN_EMIT);
	
	return 0;
}
