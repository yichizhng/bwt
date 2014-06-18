#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "rdtscll.h"
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// A variation on the search test; a preliminary implementation of gapped
// alignment

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

  // Testing!
  // Pull two 15-nt sequences from the genome and try aligning

  buf = malloc(30); // The C/C++ standard guarantees that sizeof(char) == 1
  srand(time(0));
  //rdtscll(a);
  for (i = 0; i < 10; ++i) {
    // Pick two random spots on the genome
    j = rand() % (len-15);
    for (k = 0; k < 15; ++k) {
      buf[k] = getbase(seq, j+k);
    }
    int jj = rand() % (len-15);
    while (abs(j-jj) < 15)
      jj = rand() % (len-15);
    for (k = 0; k < 15; ++k) {
      buf[k+15] = getbase(seq, jj+k);
    }
    // Now try searching (using loc_search(), not locate(), since we
    // are no longer looking for an exact match on one part of the genome)
    // for the string
    //jj = locate(fmi, buf, 30);
    int start, end, start2, end2;
    int nmatched = mms(fmi, buf, 30, &start, &end);
    int nmatched2 = mms(fmi, buf, 30-nmatched, &start2, &end2);
    for (int kk = start; kk < end; ++kk) {
      printf("%d %d\n", unc_sa(fmi, kk), jj);
      printf("%d bases matched\n", nmatched);
      for (int iii = (30 - nmatched); iii < 30; ++iii)
	putchar("ACGT"[buf[iii]]);
      putchar('\n');
      printseq(seq, unc_sa(fmi, kk), nmatched);
    }
    for (int kk = start2; kk < end2; ++kk) {
      printf("%d %d\n", unc_sa(fmi, kk), j);
      printf("%d bases matched\n", nmatched2);
      for (int iii = 0; iii < nmatched2; ++iii)
	putchar("ACGT"[buf[iii]]);
      putchar('\n');
      printseq(seq, unc_sa(fmi, kk), nmatched2);
    }
    putchar('\n');
  }
  //rdtscll(b);
  destroy_fmi(fmi);
  free(seq);
  free(buf);
  return 0;
}
