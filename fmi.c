// Implements a FM-index

#include "hwt.h"

/* Goal: implement approximate backwards search. Backwards search on
 * a FM index is linear in pattern length (i.e. not dependent on index
 * length), but the "approximate" part makes life harder for us (maybe we
 * can use more memory to compensate?). We need to somehow avoid exponential
 * time behavior, but allowing more errors is what causes that... 
 * (I believe the naive algorithm gives us O(p^(e+1)) where e is the number 
 * of errors allowed; this is bad
 */

typedef struct _fi {
	// For our purposes, the alphabet has four characters in it
	// Which shall be known as A ((char) 0), C ((char) 1), G ((char) 2),
	// and T ((char) 3).

	struct wt_block *blocks_before;
	struct wt_block *blocks_after;
	int block_size;
	// I suspect this should be hard-coded, after doing some testing
	// to see what size works best

	int count[4]; // Cumulative counts over the string of
	// (lexicographically) smaller characters

	// In theory we could store the Occ table here, but that can
	// be calculated in constant time, so there's not all that much	
	// point in using 17x the memory to speed that up.

	int before_len; // The number of blocks in before_len
	int before_last; // Number of base pairs in the last block
	int after_len;
	int after_last
} fm_index;

void fm_index_init(fm_index *id, const char *bwt) {
	// TODO: bwt is currently considered to be an uncompressed string
	// representation of the Burrows-Wheeler transform of the DNA
	// sequence; therefore its characters are 0, 1, 2, 3, and $.
	int i, j;
	// id is assumed to already have sufficient space allocated
	// e.g. id = malloc(sizeof(fm_index));

	id->bwt = bwt;
	// TODO: We should actually index bwt in some way and put
	// that (self-)index in our struct. See hwt.c
	id->count[0] = 0;
	id->count[1] = 0;
	id->count[2] = 0;
	id->count[3] = 0;
	// Iterate over the thingy, adding up the thingies
	for (i = strlen(bwt); i; --i)
		if !(bwt[i] & ~3) // That is, if bwt[i] is 0-3
			id->count[bwt[i]]++;
	id->len = len;
}

int occ(fm_index *id, char c, int k) {
	// c indicates the character we're interested in, k the position
	// in the string
	// Returns the count of character c within the first k characters
	// of the Burrows-Wheeler transformed string
	// TODO: implement (this is largely dependent on the way we
	// compress the BWT: I think fixed-block wavelet tree is the
	// way we're going)
	return -1;
}

int count(fm_index *id, char *pattern) {
	// Counts the occurences of the pattern as a substring of the
	// indexed file; this runs in O(P) (that is, it does NOT depend
	// on index length)
	// TODO: implement (it's some trickery with occ and suffixes)
	return -1;
}

int locate(fm_index *id, char *pattern) {
	// Finds a single occurence of the pattern in O(log(n)) time
	// TODO: implement
	return -1;
}

char *extract(fm_index *id, int pos, int len) {
	// Extracts the substring from pos to len of the original string
	// Returns a null-terminated string of length len+1
	// TODO: implement
	return NULL;
}
