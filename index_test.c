// Tries aligning reads against an index read from file

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"
#include "fileio.h"
#include "rdtscll.h"
#include "time.h"

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s seqfile indexfile\n", argv[0]);
    exit(-1);
  }
  char *seq, *seqfile, *indexfile, *buf, c;
  fm_index *fmi;
  int len;
  int i, j, k, jj;
  FILE *sfp, *ifp;
  seqfile = argv[1];
  indexfile = argv[2];
  sfp = fopen(seqfile, "rb");
  if (sfp == 0) {
    fprintf(stderr, "Could not open sequence\n");
    exit(-1);
  }
  fseek(sfp, 0L, SEEK_END);
  len = ftell(sfp);
  rewind(sfp);
  seq = malloc(len/4+1);
  for (i = 0; i < len/4 + 1; ++i) {
    switch(fgetc(sfp)) {
    case 'C': c = 64; break;
    case 'G': c = 128; break;
    case 'T': c = 192; break;
    default: c = 0;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 16; break;
    case 'G': c ^= 32; break;
    case 'T': c ^= 48;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 4; break;
    case 'G': c ^= 8; break;
    case 'T': c ^= 12;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 1; break;
    case 'G': c ^= 2; break;
    case 'T': c ^= 3;
    }
    seq[i] = c;
  }
  // Handle the last character (which is at seq[len/4]
  c = 0;
  for (i = 0; i < len&3; ++i) {
    switch(fgetc(sfp)) {
    case 'C': c ^= 64 >> (2 * i); break;
    case 'G': c ^= 128 >> (2 * i); break;
    case 'T': c ^= 192 >> (2 * i);
    }
    seq[len/4] = c;
  }
  fclose(sfp);
  
  // Open index file
  ifp = fopen(indexfile, "rb");
  if (ifp == 0) {
    fprintf(stderr, "Could not open index file");
    exit(-1);
  }
  fmi = read_index(seq, ifp);
  fclose(ifp);
  // And we're done! Well, okay, we might want to align some sequences.
  int seqlen = 30;
  buf = malloc(seqlen); // The C/C++ standard guarantees that sizeof(char) == 1
  srand(time(0));
  long long a, b;
  rdtscll(a);
  
  for (i = 0; i < 1000000; ++i) {
    // Pick some randomish location to start from (i.e. anywhere from 0
    // to len-16)
    j = rand() % (len-seqlen);
    for (k = 0; k < seqlen; ++k) {
      buf[k] = getbase(seq, j+k);
    }
    jj = locate(fmi, buf, seqlen);
    if (j != jj && j != -1) {
      printf("Ruh roh ");
      printf("%d %d\n", j, jj); }
  }
  rdtscll(b);
  fprintf(stderr, "Took %lld cycles to search 1000000 %dbp sequences\n",
	  b-a, seqlen);
  fprintf(stderr, "(%f seconds), over a genome of length %d\n", 
	 ((double)(b-a)) / 2500000000, len);
  destroy_fmi(fmi);
  free(seq);
  free(buf);
  return 0;
}
