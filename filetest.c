#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "rdtscll.h"
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"
#include "fileio.h"

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// Regression test for the file I/O functionality

// Writes an index to file, then reads it back and tries aligning reads
// against it

int main(int argc, char **argv) {
  // We take our input filename from argv
  int len, i, j, k, jj;
  char *seq, *buf;
  unsigned char c;
  long long a, b;
  fm_index *fmi;
  FILE *fp;
  if (argc == 1) {
    printf("Usage: searchtest seq_file");
    exit(-1);
  }
  fp = fopen(argv[1], "rb");
  fseek(fp, 0L, SEEK_END);
  len = ftell(fp);
  fseek(fp, 0L, SEEK_SET); // This is *technically* not portable
  seq = malloc(len/4+1);
  for (i = 0; i < len/4 + 1; ++i) {
    switch(fgetc(fp)) {
    case 'C': c = 64; break;
    case 'G': c = 128; break;
    case 'T': c = 192; break;
    default: c = 0;
    }
    switch(fgetc(fp)) {
    case 'C': c ^= 16; break;
    case 'G': c ^= 32; break;
    case 'T': c ^= 48;
    }
    switch(fgetc(fp)) {
    case 'C': c ^= 4; break;
    case 'G': c ^= 8; break;
    case 'T': c ^= 12;
    }
    switch(fgetc(fp)) {
    case 'C': c ^= 1; break;
    case 'G': c ^= 2; break;
    case 'T': c ^= 3;
    }
    seq[i] = c;
  }
  // Handle the last character (which is at seq[len/4]
  c = 0;
  for (i = 0; i < len&3; ++i) {
    switch(fgetc(fp)) {
    case 'C': c ^= 64 >> (2 * i); break;
    case 'G': c ^= 128 >> (2 * i); break;
    case 'T': c ^= 192 >> (2 * i);
    }
    seq[len/4] = c;
  }
  fclose(fp);
  // Now that we've loaded the sequence (ish) we can build an fm-index on it
  fmi = make_fmi(seq, len);

  // Write the index to a (temporary) file
  FILE *f = tmpfile();
  write_index(fmi, f);
  rewind(f);
  destroy_fmi(fmi);
  fmi = read_index(seq, f);
  fclose(f);
  if (fmi == NULL) {
    fprintf(stderr, "Error reading from file\n");
    exit(-1);
  }

  // Do some fun tests (load up a length 30 sequence (starting from anywhere
  // on the "genome") and backwards search for it on the fm-index (and we're
  // going to fix locate() now too)
  buf = malloc(16); // The C/C++ standard guarantees that sizeof(char) == 1
  srand(time(0));
  rdtscll(a);
  for (i = 0; i < 1000000; ++i) {
    // Pick some randomish location to start from (i.e. anywhere from 0
    // to len-16)
    j = rand() % (len-16);
    for (k = 0; k < 16; ++k) {
      buf[k] = getbase(seq, j+k);
    }
    jj = locate(fmi, buf, 16);
    if (j != jj && j != -1) {
      printf("Ruh roh ");
      printf("%d %d\n", j, jj); }
  }
  rdtscll(b);
  fprintf(stderr, "Took %lld cycles to search 1000000 16bp sequences\n",
	 b-a);
  fprintf(stderr, "(%f seconds), over a genome of length %d\n", 
	 ((double)(b-a)) / 2500000000, len);
  destroy_fmi(fmi);
  free(seq);
  free(buf);
  return 0;
}
