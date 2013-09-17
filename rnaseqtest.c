#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "rdtscll.h"
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"
#include "smw.h"

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
  // TODO: Write a mechanism (command line argument perhaps) to switch
  // between the algorithms
  // csuff_arr() uses less memory but is slower for all but the most
  // extreme cases
  // histsort(), on the other hand, is cache-friendly and multithreaded
  // idxs = csuff_arr(str, len);
  fmi = malloc(sizeof(fm_index));
  fmi->idxs = malloc((len+1) / 32 * sizeof(int));
  // idxs is probably more properly referred to as "CSA"
  for (i = 0; i < (len+1)/32; ++i)
    fmi->idxs[i] = idxs[32 * i];
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
  return fmi->C[getbase(fmi->bwt,idx - (idx > fmi->endloc))] +
    rank(fmi, getbase(fmi->bwt,idx - (idx > fmi->endloc)), idx);
}

int rank(const fm_index *fmi, char c, int idx) {
	if (idx > fmi->endloc)
		idx--;
	return seq_rank(fmi->bwt, fmi->rank_index, 16, idx, c, fmi->lookup);
}

// Runs in O(m) time
int reverse_search(const fm_index *fmi, const char *pattern, int len) {
  int start, end, i;
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      return 0;
    }
    start = fmi->C[pattern[i]] + 
      rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] +
      rank(fmi, pattern[i], end);
  }
  return end - start+1;
}

int unc_sa(const fm_index *fmi, int idx) {
  // Calculates SA[idx] given an fm-index ("enhancedish partial suffix array"?)
  int i;
  for (i = 0; idx & 31; ++i) {
    // Use the LF-mapping to find the rotation previous to idx
    idx = lf(fmi, idx);
  }
  return fmi->idxs[idx/32] + i;
}

// Runs in O(log(n) + m) time
int locate(const fm_index *fmi, const char *pattern, int len) {
  // Find the (first[0]) instance of a given sequence in a given fm-index
  // Returns -1 if none are found
  // [0] "first" in terms of location in the suffix array; i.e. the match
  // whose corresponding rotation (or equivalently, suffix) comes first
  // lexicographically; this is largely irrelevant in any real usage
  int start, end, i;
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      return -1;
    }
    start = fmi->C[pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] + rank(fmi, pattern[i], end);
  }
  //if (end - start != 1)
  //  printf("Warning: multiple matches found (returned first)\n");
  return unc_sa(fmi, start);
}

// Runs in O(m) time
void loc_search(const fm_index *fmi, const char *pattern, int len,
		int *sp, int *ep) {
  int start, end, i;
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      break;
    }
    start = fmi->C[pattern[i]] + 
      rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] +
      rank(fmi, pattern[i], end);
  }
  *sp = start;
  *ep = end;
}

// GCC suggests that static inline functions are as fast as macros
static inline int fnw_acc(int *values, int i, int j, int width) {
  // Width is len2 in the construction of values
  return values[i * (width2 + 1) + j];
}

// Performs a maximum mappable suffix search; returns the position of the
// suffix on the genome if the length found is at least equal to the cutoff,
// -1 otherwise; stores the length of the suffix in len_p. Returns an arbitrary
// (in fact, the first in the FM-index) match if multiple ones are found
// (TODO: fix that if necessary)

// One trick we can do is to drop the cutoff by some (significant) amount
// to do a search with around two nts skipped (on the pattern); this allows
// us to do a lot less stitching if we find the pattern to be immediately
// extensible (question: how fast is unc_sa()? And how fast can we calculate
// the inverse SA to "extend" searches?). Alternatively we could N-W it,
// but that's pretty slow(ish)
int mms_search(const fm_index *fmi, const char *pattern, int len,
	       int *len_p, int cutoff) {
  // Structurally almost identical to locate(), except for the way we use
  // len...
  int start, end, i;
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      // Then the longest mappable suffix starts at the next nucleotide
      *len_p = len - i + 1;
      if (len - i > cutoff) {
	return unc_sa(fmi, start) + 1;
      }
      return -1;
    }
    start = fmi->C[pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] + rank(fmi, pattern[i], end);
  }
  // Matched the entire pattern; return the first match (note that this
  // is not likely to be useful; we need some notice of where the previous
  // match was to limit our possibilities, or maybe calculate SA^-1 to
  // allow us to continue a search from elsewhere)
  *len_p = len;
  return unc_sa(fmi, start);
}

void rna_seq(const fm_index *fmi, const char *pattern, int len) {
  // Use maximum mappable *suffix* (generated via backward search, as opposed
  // to binary search, which is slow) combined with the Needleman-Wunsch
  // algorithm (for stitching purposes mostly) to align RNA sequences

  // Motivation: STARS uses maximum mappable prefix generated via binary
  // search on the FM-index; this is obviously silly (and O(m log(n))), and
  // so we're better off doing a backwards search on the suffix instead
  // to get O(m + log(n)) time
}



// A sort of real-life test to see how effective we are at actually searching
// patterns over a sequence

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
  // Do some fun tests (load up a length 30 sequence (starting from anywhere
  // on the "genome") and backwards search for it on the fm-index (and we're
  // going to fix locate() now too)
  buf = malloc(12); // The C/C++ standard guarantees that sizeof(char) == 1
  srand(time(0));
  rdtscll(a);
  for (i = 0; i < 1000000; ++i) {
    // Pick some randomish location to start from (i.e. anywhere from 0
    // to len-21)
    j = rand() % (len-12);
    for (k = 0; k < 12; ++k) {
      buf[k] = getbase(seq, j+k);
    }
    jj = locate(fmi, buf, 12);
    //if (j != jj && j != -1)
    //        printf("Ruh roh ");
    //    printf("%d %d\n", j, jj);
  }
  rdtscll(b);
  fprintf(stderr, "Took %lld cycles to search 1000000 12bp sequences\n",
	 b-a);
  fprintf(stderr, "(%f seconds), over a genome of length %d\n", 
	 ((double)(b-a)) / 2500000000, len);
  destroy_fmi(fmi);
  free(seq);
  free(buf);
  return 0;
}
