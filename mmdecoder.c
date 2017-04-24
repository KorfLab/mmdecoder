/* mmdecoder.c */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include "toolbox.h"

/* Globals */
#define STRAINS 55
#define CHROMS 22
double EXIT_PRB = 0.0001;
double MIN_EMIT = 0.0001;
char * DATA_FILE = "mmdecoder.db.gz";
ik_tvec STRvec;               /* array of 55 strain names */
ik_map  STRmap;               /* map of 55 names to numbers */
ik_tvec CHRvec;               /* array of 22 chrom names */
ik_map  CHRmap;               /* map of 22 names to numbers */
char * SEQ[STRAINS][CHROMS];  /* [strain][chrom] = SNP sequence */

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

void fasta (void) {
	int s, c;
	
	printf("size %d\n", STRvec->size);
	
	for (s = 0; s < STRAINS; s++) {
		for (c = 0; c < CHROMS; c++) {
			printf(">%s\n", STRvec->elem[s]);
		}
	}
}

static char * usage = "\
usage: mmdecoder [paramters]\n\
mandatory:\n\
  -d <file>  database file, mmdecoder.dat.gz if not specified\n\
commands:\n\
  -decode      run viterbi and posterior decoding algorithms\n\
    -f <file>  file of variant data\n\
    -x <float> exit probability, 0.0001 if not specified]\n\
    -y <float> error probability, 0.0001 if not specified]\n\
  -tree        create a distance tree in NEXUS format\n\
    -f <file>  optional additional variant file\n\
  -fasta       export sequences as fasta format\n\
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
	ik_parse_options(&argc, argv);
	
	if (ik_option("-x")) EXIT_PRB = atof(ik_option("-x"));
	if (ik_option("-y")) MIN_EMIT = atof(ik_option("-y"));
	if (ik_option("-d")) DATA_FILE = ik_option("-d");
	
	read_database(DATA_FILE);
	
	if (ik_option("-fasta")) {
		fasta();
	}
	
	printf("exit: %f\n", EXIT_PRB);
	printf("emit: %f\n", MIN_EMIT);
	
	return 0;
}

