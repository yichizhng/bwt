// Build and write an index to file

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"
#include "fileio.h"

// Command line switches:
// -h specifies histogram sort, -c specifies SACA-K
// If not specified, -h will be used below 500 million nt, -c otherwise.
// First argument will be used as the filename, second as the output file

// TODO: Implement that switch (should probably be passed as a parameter
// to make_fmi())

int main(int argc, char **argv) {
  int mode = 0, len, i;
  char *seqfile, *indexfile;
  char *seq;
  fm_index *fmi;
  char c;
  
  if (argc < 3) {
    fprintf(stderr, "Usage: %s seqfile indexfile [-c|h]\n", argv[0]);
  }
  seqfile = argv[1];
  indexfile = argv[2];
  /*
  // Parse last argument if present
  if (argc > 3)
    if (argv[3][0] == '-') {
      if (argv[3][1] == 'h')
	mode = 1;
      else if (argv[3][0] == 'c')
	mode = -1;
      else
	printf("Invalid switch\n");
	}*/
  FILE *ifp, *ofp;
  ifp = fopen(seqfile, "rb");
  if (ifp == 0) {
    fprintf(stderr, "Couldn't open index file\n");
  }
  ofp = fopen(indexfile, "w"); // wx may be better, but that's a C2011 thing
  if (ofp == 0) {
    fprintf(stderr, "Couldn't write to output file\n");
  }
  fseek(ifp, 0L, SEEK_END);
  len = ftell(ifp);
  rewind(ifp);
  seq = malloc(len/4+1);
  for (i = 0; i < len/4 + 1; ++i) {
    switch(fgetc(ifp)) {
    case 'C': c = 64; break;
    case 'G': c = 128; break;
    case 'T': c = 192; break;
    default: c = 0;
    }
    switch(fgetc(ifp)) {
    case 'C': c ^= 16; break;
    case 'G': c ^= 32; break;
    case 'T': c ^= 48;
    }
    switch(fgetc(ifp)) {
    case 'C': c ^= 4; break;
    case 'G': c ^= 8; break;
    case 'T': c ^= 12;
    }
    switch(fgetc(ifp)) {
    case 'C': c ^= 1; break;
    case 'G': c ^= 2; break;
    case 'T': c ^= 3;
    }
    seq[i] = c;
  }
  // Handle the last character (which is at seq[len/4]
  c = 0;
  for (i = 0; i < len&3; ++i) {
    switch(fgetc(ifp)) {
    case 'C': c ^= 64 >> (2 * i); break;
    case 'G': c ^= 128 >> (2 * i); break;
    case 'T': c ^= 192 >> (2 * i);
    }
    seq[len/4] = c;
  }
  fclose(ifp);

  // Make the fmi
  fmi = make_fmi(seq, len);
  write_index(fmi, ofp);
  fclose(ofp);
  destroy_fmi(fmi);
  free(seq);
  return 0;
}
