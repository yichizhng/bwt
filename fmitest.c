#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "rdtscll.h"
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"

/* Some calculations
 * Space needed for uncompressed SA representation of human genome
 * 3.6 * 10^9 nts * 8.25 bytes/nt = 2.97*10^10 bytes = 28324.127MB = 27.66GB
 * This is an amount which is quite possible to store in main memory on a
 * modern server (or a very powerful consumer desktop)
 */

// Now with pthreads inside (and no concurrency primitives, because we have
// no data races!) for speed purposes

// Note that backwards search is O(n log m) if we assume m << 4^n; this is
// largely irrelevant (and a product of statistics)

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

typedef struct _fmi {
	char *bwt; // Burrows-Wheeler transformed string
	int *idxs; // Partial index (taken from the suffix array) for
	// reverse BWT
	int **rank_index;
	unsigned char* lookup;
	int endloc; // In particular, idxs[endloc] == 0
	int C[5]; // Prefix sum of counts of symbols
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
  // Take the bwt of str. Note that str should be given in compressed
  // form.
  int *idxs, i;
  fm_index *fmi;
  // histsort() builds the suffix array for the string
  idxs = histsort(str, len);
  // Replacing histogram sort with saca-k
  // idxs = csuff_arr(str, len);
	fmi = malloc(sizeof(fm_index));
	//fmi->idxs = idxs;
	fmi->idxs = malloc((len+1) / 32 * sizeof(int));
	for (i = 0; i < (len+1)/32; ++i)
		fmi->idxs[i] = idxs[32 * i];
	fmi->bwt = malloc((len+3)/4);
	fmi->len = len;
	// Uncompressed version (for the stuff commented out in main)
	// sprintbwt(fmi->bwt, str, idxs, len+1);
	// Compressed version (for actual usage)
	// Note that sprintcbwt is a rather slow function, due mostly to the insane
	// number of cache misses.
	fmi->endloc = sprintcbwt(str, idxs, len, fmi->bwt);
	free(idxs); //Goodbye memory usage
	// Calculate rank index (see seqindex.c)
	fmi->lookup = lookup_table();
	// Block size currently hardcoded at 16 base pairs (4 bytes)
	// TODO: Should probably be increased in the future; this is doable
	// without recoding the lookup table (we only need to require that
	// blocksize is a multiple of 16), but I'm not awake enough now
	fmi->rank_index = seq_index(fmi->bwt, len, 16, fmi->lookup);
	// Experiment with block sizes?
	// Calculate C (based on fmi->rank_index)
	fmi->C[0] = 1;
	fmi->C[1] = 1         + fmi->rank_index[(len+15)/16][0];
	fmi->C[2] = fmi->C[1] + fmi->rank_index[(len+15)/16][1];
	fmi->C[3] = fmi->C[2] + fmi->rank_index[(len+15)/16][2];
	fmi->C[4] = fmi->C[3] + fmi->rank_index[(len+15)/16][3];
	return fmi;
}

int rank(const fm_index *, char, int);

// Implements the last-first mapping on the FM-index
// LF(i) = C[L[i]] + Occ(L[i], i)
static inline int lf(const fm_index *fmi, int idx) {
	if (idx == fmi->endloc)
		return 0;
	//if (idx > fmi->endloc)
	//	idx--;
	return fmi->C[getbase(fmi->bwt,idx - (idx > fmi->endloc))] +
		 rank(fmi, getbase(fmi->bwt,idx - (idx > fmi->endloc)), idx);
}

int rank(const fm_index *fmi, char c, int idx) {
	// Calculates the number of character c which occur before index idx
	// Also known as Occ() in relation to the FM-index

	// Implementation using the rank index
	// If idx is after the $ symbol, _decrease_ idx by one (obviously
	// we will never ask for the rank of the $ symbol)
	if (idx > fmi->endloc)
		idx--;
	return seq_rank(fmi->bwt, fmi->rank_index, 16, idx, c, fmi->lookup);
}

// TODO: SSE implementation of one of the dynamic algorithms (for seed+extend
// or stitching)?

// This function is also known as Count(p) in the context of FM-indexes
// It runs in O(|p|), but note that it's actually O(|p| log m) if
// m << 4^p and some statistical assumptions hold.
int reverse_search(const fm_index *fmi, const char *pattern, int len) {
	// Searches for a pattern in fmi and returns the number
	// of matches (for now it also prints the indices of said matches)

	// pattern is taken in uncompressed form; because the patterns
	// are relatively short (on the order of 100 base pairs) this doesn't
	// really affect anything.

	// len is the length of the pattern.
	int start, end, i;
//	int j, x;
	//printf("Searching for pattern ");
	//for (i = 0; i < len; ++i)
	//	putchar('0' + pattern[i]);
	//printf(":\n");
	start = fmi->C[pattern[len-1]];
	end = fmi->C[pattern[len-1]+1];

	// Iterate backwards through the string (because the fmi implements
	// a sort of suffix array, sorted in an order which makes it very
	// efficient for backwards search)
	for (i = len-2; i >= 0; --i) {
		if (end <= start) {
	//		puts("Pattern not found");
			return 0;
		}
		start = fmi->C[pattern[i]] + 
			rank(fmi, pattern[i], start);
		end = fmi->C[pattern[i]] +
			rank(fmi, pattern[i], end);
	}
	
	// Output indices (remove later if this ends up in actual code,
	// it's slow and has no purpose here (In particular, printing
	// a particular match is a different function, called locate())
/*	printf("Matches at indices: ");
	for (i = start-1; i < end; ++i) {
//		printf("%d ", fmi->idxs[i]);
		x = i;
		for (j = 0; x & 31; ++j) {
			x = lf(fmi, x);
		}
		// TODO: Is that right?
		printf("%d ", j + fmi->idxs[x] - 1);
	}
	putchar('\n');
*/	return end - start+1;
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

// Note that this is single-threaded; we parallelize by calling it
// on different threads. Since it doesn't modify the index the only
// thing we have to worry about is cache locality (which is a real concern,
// but not really something we can fix)
void rna_seq(const fm_index *fmi, const char *pattern, int len) {
	// Copy pattern into a different aray for cyclic search
	char *pat = malloc(2*len);
	int i, j, start, end;
	memcpy(pat, pattern, len);
	memcpy(pat+len, pattern, len);
	for (i = len; i; --i) {
		loc_search(fmi, pattern+i, 14, &start, &end);
	 	if (start < end) {
		}
	}
	// Identify non-mapped segments and deal with them
	free(pat);
}

struct thread_args {
	fm_index *fmi;
	int num;
	char *pats;
	int start;
};

void *worker(void *arg) {
	int i, j, max;
	max = (((struct thread_args*)arg)->start) + 
		(((struct thread_args*)arg)->num);
	//printf("Starting thread to process %d patterns\n",
	//	max/4);
	for (i = ((int) ((struct thread_args*)arg)->start); i < max; i ++) {
		reverse_search(((struct thread_args*)arg)->fmi, 
			       ((struct thread_args*)arg)->pats + i, 12);
	}
	return NULL;
}

void startWorkers(fm_index *fmi, int num, char *pats) {
	// Starts 4 threads, each of which will process num/4 patterns
	// There are no synchronization primitives because they're unnecessary.
	pthread_t threads[4];
	int i;
	void *status;
	struct thread_args **ta = malloc(4 * sizeof (struct thread_args*));
	for (i = 0; i < 4; ++i) {
		ta[i] = malloc(sizeof(struct thread_args));
		ta[i]->fmi = fmi;
		ta[i]->num = num/4;
		ta[i]->pats = pats;
		ta[i]->start = i * (num/4);
	}
	for (i = 0; i < 4; ++i) {
		pthread_create(&threads[i], NULL, worker, (void *)ta[i]);
	}
	for (i = 0; i < 4; ++i) {
		pthread_join(threads[i], &status);
		// Ignoring return code, we didn't generate one anyway
		free(ta[i]);
	}
	free(ta);
	return;
}

// Misfeature: Index construction is O(n log n) on average; this is fast enough
// to dominate SACA-K below a billion base pairs or so, but uses too much
// memory

// For our use case (RNA sequencing) we can reverse search without backtracking
// on some fixed size seed (e.g. 14); STAR uses maximum mappable seed,
// which is sort of interesting as well
int main(int argc, char **argv) {
	
	// Performance testing code
	// To be honest, this should have the same performance as
	// histsortcomptest.c's testing code.

	// Single-threaded it takes about 4.6 seconds to construct the index
	// for 10^7 elements ; this suggests ~90 seconds for 10^8 elements
	// and 1800 seconds (30 minutes) for 10^9. Fortunately index-building
	// is a one-time thing, so runtime is not all that significant.

	// For some reason multithreaded index building is slower (overhead
	// and cache conflicts no doubt) below about 10 million base pairs
	// OTOH 4-threaded histogram sort is faster than SACA-K below
	// a billion nucleotides, so memory usage aside this should be
	// fine (and they make servers with enough memory to run it anyway)
	// Performance testing code
	int len, i, j, k;
	char *str, *pats;
	fm_index *fmi;
	unsigned long long a, b;
	len = atoi(argv[1]);
	str = malloc(len/4 + 1);
	for (i = 0; i < len/4 + 1; ++i)
		str[i] = rand();
	switch (len % 4) {
	case 1:
	  str[len/4] &= 0xC0;
	  break;
	case 2:
	  str[len/4] &= 0xF0;
	  break;
	case 3:
	  str[len/4] &= 0xFC;
	  break;
	case 0:
	  str[len/4] = 0;
	}
	rdtscll(a);
	fmi = make_fmi(str, len);
	rdtscll(b);
	printf("Built index with %d base pairs in %lld cycles (%f s)\n",
		len, b-a, ((double)(b-a)) / 2500000000.);
	printf("(%f cycles per base pair (%e seconds))\n", ((double)(b-a)) \
		/ len, ((double)(b-a)) / (len * 2500000000.));
	pats = malloc(sizeof(char) * 10000011);
	for (i = 0; i < 10000011; ++i) {
	  pats[i] = rand() & 3;
	}
	// Comment: This is huge in 64-bit mode and thus makes valgrind
	// runs somewhat untenable (at least with that number of bps to
	// search; we can, of course, always scale things down but this
	// reveals more of the performance overhead from threading and less
	// asymptotic performance (which is what we actually care about)
	rdtscll(a);
	startWorkers(fmi, 10000000, pats);
	rdtscll(b);
	printf("Searched 10000000 12bp sequences in %lld cycles (%f s)\n",
		b-a, ((double)(b-a)) / 2500000000.);
	printf("(%f cycles per base pair (%e seconds))\n", ((double)(b-a)) \
		/ 120000000., ((double)(b-a)) / 300000000000000000.);
	
	destroy_fmi(fmi);
	// With current settings the FMI with 1000000 base pairs takes
	// 1.5 million bytes to store (50% auxiliary data)
	// which is kind of bad but yeah
	free(str);
	free(pats);
	/* Correctness verification code
	int x;
	char *str, *pat;
	fm_index *fmi;	
	str = malloc(14 * sizeof(char));
	pat = malloc(4 * sizeof(char));
	pat[0] = 1;
	pat[1] = 0;
	pat[2] = 1;
	strcpy(str, ";;;;;I@7456;;");
	// String is 0323102110000313031003110312
	// BWT is 2101013010$133123331031020000
	// 8 18 4 23
	fmi = make_fmi(str, 52);
	x = reverse_search(fmi, pat, 2); // Search for "03" (AT) over str
	printf("%d matches\n", x);
	destroy_fmi(fmi);
	free(str);
	free(pat);
	return 0; */
}
