// Calculates a rank index for a given compressed BWT sequence
// We have no particular need for the select() operation, although we can
// implement it in nearly constant time if we allow certain constraints
// regarding our nucleotide sequences (which are roughly true); nucleotide
// sequences are also essentially incompressible, so there's no point
// using a wavelet tree or RRR

// TODO (maybe): We can optimize some code if blocksize is required to be
// a power of 2 rather than a multiple of 4

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqindex.h"
#include "histsortcomp.h"

static inline unsigned char getbase(const char *str, int idx) {
  // Gets the base at the appropriate index
  return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

int **seq_index(char *bwt, int len, int blocksize,const unsigned char *tbl) {
  // len is, as usual, the length of the bwt. bwt is in compressed form.
  // blocksize is assumed to be a multiple of 4 and is the number of
  // base pairs per block.
  // index, unlike lookup_table, is a 2D array, largely because I feel
  // like it.
  // tbl is the return value of lookup_table() (because who wants to call
  // it more than once?).
  int i, j, **index;
  unsigned char c; // Apparently C does something odd when casting
  // unsigned char to unsigned int... by sign extending
  // which is certainly not what I want
  
  // Compiling with full warnings will give you tons of rubbish about
  // array indexing with a char; I know perfectly well what I'm doing.
  index = malloc((1+((len+blocksize-1)/blocksize)) * sizeof(int *));
  
  index[0] = calloc(4, sizeof(int));
  // Do the first loop separately to avoid complicated logic
  for (j = 1; j < 1 + (len/blocksize); ++j) {
    index[j] = malloc(4 * sizeof(int));
    index[j][0] = index[j-1][0];
    index[j][1] = index[j-1][1];
    index[j][2] = index[j-1][2];
    index[j][3] = index[j-1][3];
    for (i = -blocksize/4; i; ++i) {
      // Beautifully obfuscated loop unroll
      c = bwt[j*(blocksize/4) + i];
      index[j][0] += tbl[4 * c];
      index[j][1] += tbl[4 * c + 1];
      index[j][2] += tbl[4 * c + 2];
      index[j][3] += tbl[4 * c + 3];
    }
  }
  // Handle some edge cases; notice how j went to len/blocksize rather
  // than (len + (blocksize-1)) / blocksize
  if (len % blocksize) {
    index[j] = malloc(4 * sizeof(int));
    index[j][0] = index[j-1][0];
    index[j][1] = index[j-1][1];
    index[j][2] = index[j-1][2];
    index[j][3] = index[j-1][3];
    // Count up the base pairs in the last block; this loop
    // will also be unrolled, just for old time's sake
    for (i = -blocksize/4; i < ((len%blocksize) - blocksize)/4 - !!(len%4); ++i) {
      c = bwt[j*(blocksize/4) + i];
      index[j][0] += tbl[4 * c];
      index[j][1] += tbl[4 * c + 1];
      index[j][2] += tbl[4 * c + 2];
      index[j][3] += tbl[4 * c + 3];
    }
    // Now we handle the last few base pairs (there are between
    // 0 and 3 of them)
    for (i = 0; i < (len&3); ++i) {
      // Bitwise manipulations to make indexing easier
      c = getbase(bwt, (len & 0xFFFFFFC) ^ i);
      index[j][c]++;
    }
  }
  //for (i = 0; i < 1+((len+blocksize-1)/blocksize); ++i) {
  //	printf("%d %d %d %d\n", index[i][0], index[i][1],
  //		index[i][2], index[i][3]);
  //}
  return index;
}

// Calculates the index of character c (between 0 and 3) at index idx
// in the given sequence index in the given Burrows-Wheeler transform
int seq_rank(unsigned char *bwt, int **sidx, int blocksize, int idx, char c, const unsigned char * lookup) {
  int x, i;
  // First we look up the appropriate block prefix sum
  x = sidx[idx/blocksize][c];
  // Now we iterate through the block
  /* unoptimized version */
  //for (i = 0; i < idx % blocksize; i ++) {
  //	if (getbase(bwt,(idx/blocksize)*blocksize + i) == c) {
  //		++y;
  //	}
  //}
  for (i = (idx/blocksize)*blocksize/4; i < idx/4; i++) {
    x += lookup[4*bwt[i] + c];
  }
  for (i = 0 ; i < idx % 4; ++i) {
    if (getbase(bwt, (idx & 0xFFFFFFFC) ^ i) == c)
      ++x;
  }
  // printf("Rank of %c is %d at index %d\n", '0'+c, x, idx);
  return x;
}

// TODO: Maybe hardcode this table or something, stop it from being called
// multiple times (e.g. we could just write a function to print it out)
unsigned char * lookup_table() {
  // Calculates the lookup table for one byte of the sequence (i.e.
  // 4 base pairs). 256 possible combinations * 4 entries per byte
  // = 1024 bytes for our table. Used to speed up some iterations/searches
  unsigned char *tbl = malloc(1024);
  memset(tbl, 0, 1024);
  int i, j, k, l;
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      for (k = 0; k < 4; ++k)
	for (l = 0; l < 4; ++l) {
	  // Brilliantly obfsucated loop
	  // In essence, tbl[4x + i] (for i < 4, x < 256) is
	  // the count of base pair i in the byte represented by x
	  tbl[257*i + 64*j + 16*k + 4*l]++;
	  tbl[256*i + 65*j + 16*k + 4*l]++;
	  tbl[256*i + 64*j + 17*k + 4*l]++;
	  tbl[256*i + 64*j + 16*k + 5*l]++;
	}
  // There are possibly more intuitive ways to calculate this, but this
  // makes for the simplest code, and it's enough to know what the lookup table
  // is without knowing how to construct one
  return tbl;
}

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

// Comment: rather memory intensive
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
  fmi->idxs = malloc((1 + (len / 32)) * sizeof(int));
  // idxs is probably more properly referred to as "CSA"
  for (i = 0; i < (1+(len / 32)); ++i)
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

// The same, but using SACA-K instead 
fm_index *make_fmi_sacak(const char *str, int len) {
  int *idxs, i;
  fm_index *fmi;
   idxs = csuff_arr(str, len);
  fmi = malloc(sizeof(fm_index));
  fmi->idxs = malloc((1 + (len / 32)) * sizeof(int));
  for (i = 0; i < (1+(len / 32)); ++i)
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

int lf(const fm_index *fmi, int idx) {
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
  int i, x;
  for (i = 0; idx & 31; ++i) {
    // Use the LF-mapping to find the rotation previous to idx
    idx = lf(fmi, idx);
  }
  x = fmi->idxs[idx/32] + i;
  if (x > fmi->len)
    x -= fmi->len + 1;
  return x;
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
  if (end - start != 1)
    printf("Warning: multiple matches found (returned first)\n");
  return unc_sa(fmi, start);
}

// Runs in O(m) time. Finds the locations of all matches to the pattern.
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

// Finds the maximum mappable suffix of the pattern; returns the length
// matched (starting at the end of the pattern) and stores the range
// of matches in sp and ep.
int mms(const fm_index *fmi, const char *pattern, int len, int *sp, int *ep) {
  int start, end, i;
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      break;
    }
    *sp = start;
    *ep = end;
    start = fmi->C[pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] + rank(fmi, pattern[i], end);
  }
  if (end <= start) // Didn't finish matching
    return len - i - 2;
  else { // Finished matching
    *sp = start;
    *ep = end;
    return len - i - 1;
  }
}

// Prints part of a compressed sequence in more human readable format
void printseq(const char *seq, int startidx, int len) {
  const char *nts = "ACGT";
  for (int i = 0; i < len; ++i) {
    putchar(nts[getbase(seq, startidx+i)]);
  }
  putchar('\n');
}
