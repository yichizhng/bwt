#include <stdio.h>
#include <stdlib.h>
#include "seqindex.h"
#include "fmitest.h"

// Implements reading FM-indexes from file and writing them to file
// Something that I really should have done a while back -.-

typedef struct _fmi {
  char *bwt;
  int *idxs;
  int **rank_index;
  unsigned char* lookup;
  int endloc;
  int C[5];
  int len;
} fm_index;

// We really don't need to know how to construct a FM-index to do this

// Return 0 on success, the file error code otherwise
int write_fmi(FILE *out, fm_index *fmi) {
  // Not really going to bother checking whether out is valid, that'll be
  // the caller's problem
  // The first thing we should really write out is the length
  fwrite(&(fmi->len), sizeof(int), 1, out);
  
  // Write out the Burrows-Wheeler transform
  fwrite(fmi->bwt, sizeof(char), (fmi->len + 3)/4, out);

  // Write out idxs (the compressed suffix array)
  fwrite(fmi->idxs, sizeof(int), (fmi->len + 1) / 32, out);

  // Don't write out the rank index, it's easy enough to reconstruct
  // and is kind of a waste of disk space
  // Don't bother writing the lookup table, we can construct it
  // Write out endloc
  fwrite(&(fmi->endloc), sizeof(int), 1, out);
  
  // Write out C
  fwrite(fmi->C, sizeof(int), 5, out);

  // Return the error code
  return ferror(out);
}

// Returns 0 on success, the file error code otherwise
int read_fmi(FILE *in, fm_index *fmi) {
  // fmi is assumed to be properly allocated already, although its fields
  // are not.

  // Keep in mind that we have to reconstruct the lookup table (trivial) and
  // the rank index (trivial yet annoying)
  int len;
  
  // Read in everything that we actually stored
  fread(&(fmi->len), sizeof(int), 1, in);
  len = fmi->len;
  fmi->bwt = malloc((len + 3)/4);
  // Read bwt from file
  fread(fmi->bwt, sizeof(char), (len+3)/4, in);
  fmi->idxs = malloc((fmi->len+1)/32);
  fread(fmi->idxs, sizeof(int), (len+1)/32, in);
  fread(&(fmi->endloc), sizeof(int), 1, in);
  fread(fmi->C, sizeof(int), 5, in);
  // Construct the lookup table and rank index

  fmi->lookup = lookup_table();

  // Compression ratio on rank_index is hardcoded at 32.
  // TODO: change this code if that ever changes, or maybe code a macro for it
  // If we allow it to be user-changable we could write it to the file as well
  // This seems rather inadvisable though; never give the user too many
  // options ;)
  fmi->rank_index = seq_index(fmi->bwt, len, 32, fmi->lookup);
  return ferror(in);
}

// This thingy will build an FM-index and then write it to disk, because
// that is a thing that should be done
int main(int argc, char **argv) {
  // Usage: ./fmibuild genome_name index_name
  FILE *gen, *idx;
  fm_index *fmi;
  char *seq;
  int len;
  if (argc != 3) {
    exit(-1); // Me, give meaningful error messages? Never!
  }
  // Open up the genome; I suppose for now we may as well assume it's not
  // compressed, which means we must now compress it -.-;; (if it were
  // we would be able to use fread directly)

  // TODO: actually do that :) steal the code from searchtest.c probably, that
  // part is working as far as I can tell
  gen = fopen(argv[1], "rb");
  fseek(seq, 0L, SEEK_END); // Technically not portable
  len = ftell(seq);
  fseek(seq, 0L, SEEK_SET);
  
}
