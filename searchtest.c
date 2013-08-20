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

typedef struct _fmi {
	char *bwt;
	int *idxs;
	int **rank_index;
	unsigned char* lookup;
	int endloc;
	int C[5];
	int len;
} fm_index;

void destroy_fmi (fm_index *fmi) {
	int i;
	free(fmi->bwt);
	free(fmi->idxs);
	for (i = 0; i <= (fmi->len+15)/16; ++i)
		free(fmi->rank_index[i]);
	free(fmi->rank_index);
	free(fmi->lookup);
	free(fmi);
}

fm_index *make_fmi(const char *str, int len) {
  int *idxs, i;
  fm_index *fmi;
  idxs = histsort(str, len); // i.e. SA
  // idxs = csuff_arr(str, len);
  fmi = malloc(sizeof(fm_index));
  fmi->idxs = malloc((len+1) / 32 * sizeof(int));
  // idxs is probably more properly referred to as "CSA"
  for (i = 0; i < (len+1)/32; ++i)
    fmi->idxs[i] = idxs[32 * i];
  // TODO: Is this right?
  fmi->bwt = malloc((len+3)/4);
  fmi->len = len;
  fmi->endloc = sprintcbwt(str, idxs, len, fmi->bwt);
  free(idxs);
  fmi->lookup = lookup_table();
  fmi->rank_index = seq_index(fmi->bwt, len, 16, fmi->lookup);
  fmi->C[0] = 1;
  fmi->C[1] = 1         + fmi->rank_index[(len+15)/16][0];
  fmi->C[2] = fmi->C[1] + fmi->rank_index[(len+15)/16][1];
  fmi->C[3] = fmi->C[2] + fmi->rank_index[(len+15)/16][2];
  fmi->C[4] = fmi->C[3] + fmi->rank_index[(len+15)/16][3];
  return fmi;
}

int rank(const fm_index *, char, int);

static inline int lf(const fm_index *fmi, int idx) {
	if (idx == fmi->endloc)
		return 0;
	//if (idx > fmi->endloc)
	//	idx--;
	return fmi->C[getbase(fmi->bwt,idx - (idx > fmi->endloc))] +
		 rank(fmi, getbase(fmi->bwt,idx - (idx > fmi->endloc)), idx);
}

int rank(const fm_index *fmi, char c, int idx) {
	if (idx > fmi->endloc)
		idx--;
	return seq_rank(fmi->bwt, fmi->rank_index, 16, idx, c, fmi->lookup);
}

int reverse_search(const fm_index *fmi, const char *pattern, int len) {
  int start, end, i;
  start = fmi->C[pattern[len-1]]+1;
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      return 0;
    }
    start = fmi->C[pattern[i]] + 
      rank(fmi, pattern[i], start-1)+1;
    end = fmi->C[pattern[i]] +
      rank(fmi, pattern[i], end);
  }
  return end - start+1;
}

int unc_sa(const fm_index *fmi, int idx) {
  // Calculates SA[idx] given an enhancedish fm-index (a partial suffix array
  // and a BWT and some other stuff)
  int i;
  for (i = 0; idx & 31; ++i) {
    // Use the LF-mapping to find the rotation previous to idx
    idx = lf(fmi, idx);
  }
  return fmi->idxs[idx/32] + i;
}

int locate(const fm_index *fmi, const char *pattern, int len) {
  // Find the (first) instance of a given sequence in a given fm-index
  // Returns -1 if none are found
  int start, end, i;
  start = fmi->C[pattern[len-1]]+1;
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      return -1;
    }
    start = fmi->C[pattern[i]] + 
      rank(fmi, pattern[i], start-1);
    end = fmi->C[pattern[i]] +
      rank(fmi, pattern[i], end);
  }
  if (end - start != 1)
    printf("Warning: multiple matches found (returned first)\n");
  return unc_sa(fmi, start);
}

void loc_search(const fm_index *fmi, const char *pattern, int len,
	int *sp, int *ep) {
	// Searches for a pattern in fmi and returns the start and
	// end indices. This is to be used for seed searches (as such
	// it would be called with len=14 instead of, say, 100 (the latter
	// puts unrealistic demands on read accuracy)), which gives us some
	// seeds which we can extend (e.g. via one of the dynamic algorithms)
	// The advantage of this over backtracking search (e.g. Bowtie) is
	// that backtracking takes O(|p|^(1+e)) which is very painful
	// (if at least somewhat feasible for small |p|? we can use it for
	// micro-exons if that comes up)

	// This function is implemented essentially identically to the
	// previous, it just stores start and end into pointers (this
	// being better than, say, returning a struct).
	int start, end, i;
	start = fmi->C[pattern[len-1]]+1;
	end = fmi->C[pattern[len-1]+1];
	for (i = len-2; i >= 0; --i) {
		if (end <= start) {
			break;
		}
		start = fmi->C[pattern[i]] + 
			rank(fmi, pattern[i], start-1)+1;
		end = fmi->C[pattern[i]] +
			rank(fmi, pattern[i], end);
	}
	*sp = start;
	*ep = end;
}

// A sort of real-life test to see how effective we are at actually searching
// patterns over a sequence

int main(int argc, char **argv) {
  // We take our input filename from argv
  int len, i, j, k, jj;
  char *seq, *buf, c;
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
    case 'C':
      c = 64;
      break;
    case 'G':
      c = 128;
      break;
    case 'T':
      c = 192;
      break;
    default:
      c = 0;
      break;
    }
    switch(fgetc(fp)) {
    case 'C':
      c ^= 16;
      break;
    case 'G':
      c ^= 32;
      break;
    case 'T':
      c ^= 48;
    }
    switch(fgetc(fp)) {
    case 'C':
      c ^= 4;
      break;
    case 'G':
      c ^= 8;
      break;
    case 'T':
      c ^= 12;
      break;
    }
    switch(fgetc(fp)) {
    case 'C':
      c ^= 1;
      break;
    case 'G':
      c ^= 2;
      break;
    case 'T':
      c ^= 3;
      break;
    }
    seq[i] = c;
  }
  // Handle the last character (which is at seq[len/4]
  c = 0;
  for (i = 0; i < len&3; ++i) {
    switch(fgetc(fp)) {
    case 'C':
      c ^= 64 >> (2 * i);
      break;
    case 'G':
      c ^= 128 >> (2 * i);
      break;
    case 'T':
      c ^= 192 >> (2 * i);
    }
    seq[len/4] = c;
  }
  // Now that we've loaded the sequence (ish) we can build an fm-index on it
  fmi = make_fmi(seq, len);
  // Do some fun tests (load up a length 20 sequence (starting from anywhere
  // on the "genome") and backwards search for it on the fm-index (and we're
  // going to fix locate() now too)
  buf = malloc(30); // The C/C++ standard guarantees that sizeof(char) == 1
  srand(time(0));
  for (i = 0; i < 1000; ++i) {
    // Pick some randomish location to start from (i.e. anywhere from 0
    // to len-21)
    j = rand() % (len-30);
    for (k = 0; k < 30; ++k) {
      buf[k] = getbase(seq, j+k);
    }
    jj = locate(fmi, buf, 30);
    printf("%d %d\n", j, jj);
  }
  return 0;
}
