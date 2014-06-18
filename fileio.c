// Functions to write an index to file and read it back
// Does not actually store the original sequence; that seems pointless

#include "seqindex.h"
#include <stdio.h>
#include <stdlib.h>

void write_index(const fm_index *fmi, FILE *f) {
  // Writes the FM-index to file... well, the parts that take
  // time to actually generate.
  fwrite(&fmi->len, sizeof(int), 1, f);
  fwrite(fmi->C, sizeof(int), 5, f);
  fwrite(&fmi->endloc, sizeof(int), 1, f);
  fwrite(fmi->idxs, sizeof(int), ((fmi->len+1)/32), f);
  fwrite(fmi->bwt, 1, (fmi->len+3)/4, f);
  // C standard guarantees sizeof(char) to be 1
  return;
}

fm_index *read_index(const char *seq, FILE *f) {
  // Reads the BWT from file and reconstructs the FM-index
  // Returns a newly allocated FM-index.
  // Doesn't check for running out of memory, run at your own risk
  fm_index *fmi = malloc(sizeof(fm_index));
  fread(&fmi->len, sizeof(int), 1, f);
  fread(fmi->C, 5, sizeof(int), f);
  fread(&fmi->endloc, sizeof(int), 1, f);
  fmi->idxs = malloc((fmi->len+1)/32 * sizeof(int));
  fread(fmi->idxs, sizeof(int), (fmi->len+1)/32, f);
  fmi->bwt = malloc((fmi->len+3)/4);
  fread(fmi->bwt, 1, (fmi->len+3)/4, f);
  
  fmi->lookup = lookup_table();
  fmi->rank_index = seq_index(fmi->bwt, fmi->len, 16, fmi->lookup);
  return fmi;
}

int main(void) {
  printf("Hello, world!");
  return 0;
}
