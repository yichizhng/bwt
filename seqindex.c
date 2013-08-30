// Calculates a rank index for a given compressed BWT sequence
// We have no particular need for the select() operation, although we can
// implement it in nearly constant time if we allow certain constraints
// regarding our nucleotide sequences (which are roughly true); nucleotide
// sequences are also essentially incompressible, so there's no point
// using a wavelet tree or RRR

// TODO (maybe): I can optimize some code if blocksize is required to be
// a power of 2 rather than a multiple of 4

#include <stdio.h> /*TODO: remove, once I get debugging done */
#include <stdlib.h>
#include <string.h>
#include "seqindex.h"
#include "histsortcomp.h"

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// TODO: Fix code logic, this is a jumbled mess. It does work, though.
// In particular a lot of my indexing no longer makes sense to me
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

unsigned char * lookup_table() {
  // Calculates the lookup table for one byte of the sequence (i.e.
  // 4 base pairs). 256 possible combinations * 4 entries per byte
  // = 1024 bytes for our table.
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
  return tbl;
}
