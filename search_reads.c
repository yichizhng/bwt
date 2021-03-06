// Tries aligning reads from a file against an index and sequence read from
// file

// usage: search_reads seqfile indexfile readfile

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

// Reminder to self: buf length (i.e. maximum read length) is currently
// hardcoded; change to a larger value (to align longer reads) or make it
// dynamic

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s seqfile indexfile readfile\n", argv[0]);
    exit(-1);
  }
  char *seq, *seqfile, *indexfile, *readfile, *buf = malloc(256*256), *revbuf = malloc(256*256), c;
  fm_index *fmi;
  int len;
  int i, j, k, jj;
  FILE *sfp, *ifp, *rfp;
  seqfile = argv[1];
  indexfile = argv[2];
  readfile = argv[3];
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

  // And now we go read the index file
  rfp = fopen(readfile, "r");
  if (rfp == 0) {
    fprintf(stderr, "Could not open reads file");
    exit(-1);
  }
  // Read one line ("read") and try aligning it
  
  printf("Beginning alignment\n");
  int nread = 0;
  while (!feof(rfp)) {
    if (! fgets(buf, 256*256-1, rfp))
      continue;
    int forward_pos, backward_pos;
    int forward_match = 0, backward_match = 0;
    // fgets() writes the ending newline if present, so we need to remove
    // that
    if (buf[strlen(buf)-1] == '\n')
      buf[strlen(buf)-1] = 0;
    int len = strlen(buf);
    for (int k = 0; k < strlen(buf); ++k)
      revbuf[k] = buf[strlen(buf)-k-1];
    revbuf[strlen(buf)] = 0;
    while (len > 20 /* Replace with user-specified constant? */) {
      // Try aligning against the end of the read (MMS)
      int start, end;
      int matched = mms(fmi, buf, len, &start, &end);
      if (matched >= 20) {
	// Got an anchor length of >20
	// Print out the matches
	//printf("\n%d anchor(s) found with length %d for read %d\n", end - start, matched, nread);
	//for (int j = start; j < end; ++j)
	//printf("Starting at position %d\n", unc_sa(fmi, j));
	forward_match++;
	len -= matched;
	forward_pos = unc_sa(fmi, start);
      }
      else {
	len -= 1; // this constant should probably be bigger than 1 for performance
		  // reasons
      }
    }
    len = strlen(revbuf);
     while (len > 20 /* Replace with user-specified constant? */) {
      // Try aligning against the end of the read (MMS)
      int start, end;
      int matched = mms(fmi, revbuf, len, &start, &end);
      if (matched >= 20) {
	// Got an anchor length of >20
	// Print out the matches
	//printf("\n%d anchor(s) found with length %d for read %d\n", end - start, matched, nread);
	//for (int j = start; j < end; ++j)
	//printf("Starting at position %d\n", unc_sa(fmi, j));
	backward_match++;
	len -= matched;
	backward_pos = unc_sa(fmi, start);
      }
      else {
	len -= 1; // this constant should probably be bigger than 1 for performance
		  // reasons
      }
    }
    if (forward_match && backward_match && (abs(forward_pos - backward_pos) < 10000)) {
      printf("\nRead %d: Aligned both forward (%d) and backward (%d)\n",
             nread, forward_match, backward_match);
      printf("At locations %d and %d respectively\n", forward_pos, backward_pos);
      printf("%s\n", buf);
    }

    nread++;
  }
  fclose(rfp);
  
  free(buf);
  free(revbuf);
  destroy_fmi(fmi);
  free(seq);
  return 0;
}
