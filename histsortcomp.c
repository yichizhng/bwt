// Implements a histogram sort of the suffix array of an array of width 2
// unsigned ints using O(n) auxiliary space. This one works with the
// compressed form, which allows us to pack 4 base pairs into each char
// at the cost of using some bit shifting (which actually
// speeds up the program at larger array size - I blame caches)
// In any case, this uses a bit less memory than the silly version

// On the other hand, the index only needs to be built once for a given
// genome, and that can be done on a fast machine with enough memory
// (which had better be 64-bit, if we're going to use genomes with billions
// of base pairs)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "histsortcomp.h"

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

static int *saux;

// base is a pointer to the array of prefixes
// len is one more than the number of base pairs (it's the length of the bwt)
// arr is a pointer to the array of indices being sorted (i.e. we're
// repesenting each string using the prefix)
// aux is the auxiliary array being sorted into, although some trickery
// ensures that the original array will end up sorted
// depth is the current "string index" being considered.
void histhelper(const char *base, int len, int *arr, int *aux, int start, int end, int depth) {
	// We are responsible for sorting the array from index arr[start] to
	// arr[end-1], and are currently considering the base at
	// position depth (that is, getbase(base,arr[i]+depth))
	
	// There are four buckets, plus one for a string which has just
	// ended. Each of those four buckets is then sorted recursively.
	int lens[4] = {0}, i, t;
	int ptrs[4];
	// Base cases, where the array is already sorted
	if ((end - start) == 1) {
		if (arr == saux)
			aux[start] = arr[start];
			// We want the original array to be sorted
		return;
	}
	if (end == start)
		return;
	for (i = start; i != end; ++i) {
		if(arr[i] + depth != len) { 
			// Not <, the others were handled already
			lens[getbase(base,arr[i]+depth)]++;
		}
		else {
			aux[start] = arr[i];
			arr[i] = arr[start];
			arr[start] = aux[start];
			++start; // This element is known to be in the right
			// position now, so don't bother processing it again
		}
	}
	// lens will now contain the numbers of elements which should
	// go into each bucket, which allows us to assign appropriate
	// values for ptrs
	// In essence, this is a 3-pivot quicksort with no comparisons,
	// which makes it slightly faster
	ptrs[0] = start;
	ptrs[1] = ptrs[0] + lens[0];
	ptrs[2] = ptrs[1] + lens[1];
	ptrs[3] = ptrs[2] + lens[2];
	// I don't know why I kept lens[3] around either, maybe a sanity
	// check (or to avoid a branch prediction error). It lets me
	// find the end of the array, but I know where that is anyway
	for (i = start; i != end; ++i) {
		aux[ptrs[getbase(base,arr[i]+depth)]] = arr[i];
		++ptrs[getbase(base,arr[i]+depth)];
	}
	histhelper(base, len, aux, arr,
		start, start + lens[0], depth+1);
	histhelper(base, len, aux, arr,
		start+lens[0], start+lens[0]+lens[1], depth+1);
	histhelper(base, len, aux, arr,
		start+lens[0]+lens[1], end - lens[3],depth+1);
	histhelper(base, len, aux, arr,
		end - lens[3], end, depth+1);
}

struct hh_args {
	char *base;
	int len;
	int *arr;
	int *aux;
	int start;
	int end;
};

void *hh_wrapper(void *arg) {
	struct hh_args *args = (struct hh_args*)arg;
	//printf("%d %d\n", args->start, args->end);
	histhelper(args->base, args->len, args->arr, args->aux,
		args->start, args->end, 1);
	return NULL;
}

// The depth = 0 call to histhelper. Making this a separate function allows
// us to parallelize properly (we need 4 threads, and conveniently enough
// we have four buckets, and absolutely no data races between them!)
void histzero(const char *base, int len, int *arr, int *aux, int start, int end) {
	int lens[4] = {0}, i, t;
	int ptrs[4];
	struct hh_args *args[4];
	pthread_t threads[4];
	void *status;
	// Assumptions about our input allow us to do this iteration
	// somewhat more efficiently (as if it matters all that much)
	for (i = start; i != end; ++i) {
		lens[getbase(base,arr[i])]++;
	}
	ptrs[0] = start;
	ptrs[1] = ptrs[0] + lens[0];
	ptrs[2] = ptrs[1] + lens[1];
	ptrs[3] = ptrs[2] + lens[2];
	aux[0] = arr[0];
	for (i = start; i != end; ++i) {
		aux[ptrs[getbase(base,arr[i])]] = arr[i];
		++ptrs[getbase(base,arr[i])];
	}
	for (i = 0; i < 4; ++i) {
		args[i] = malloc(sizeof(struct hh_args));
		args[i]->base = base;
		args[i]->len = len;
		args[i]->arr = aux;
		args[i]->aux = arr;
	}
	args[0]->start = start;
	args[1]->start = start + lens[0];
	args[2]->start = start + lens[0] + lens[1];
	args[3]->start = end - lens[3];
	args[0]->end = args[1]->start;
	args[1]->end = args[2]->start;
	args[2]->end = args[3]->start;
	args[3]->end = end;
	for (i = 0; i < 4; ++i) {
		pthread_create(&threads[i], NULL, hh_wrapper, (void *) args[i]);
	}
	for (i = 0; i < 4; ++i) {
		pthread_join(threads[i], &status);
		free(args[i]);
	}
}

// len is the number of base pairs, but we're going to increment that to
// avoid some odd indexing problems (in particular, len becomes the length
// of the bwt'd string)
// Returns something quite like a suffix array
int * histsort(const char *str, int len) {
	int *arr = malloc((len+1) * sizeof(int));
	int i;
	int *aux = malloc((len+1) * sizeof(int));
	arr[0] = len; // Note that the last rotation leaves the $ in
	// front, so it will certainly be sorted here
	for (i = 1; i <= len; ++i)
		arr[i] = i-1;
	saux = aux;
	if (len < 10000000) { // This figure was established experimentally
		histhelper(str, len, arr, aux, 0, len+1, 0);
	}
	else {
		histzero(str, len, arr, aux, 1, len+1);
	}
	// By the power of magic and handwaving, arr will end up sorted
	free(aux);
	return arr;
	// arr is the suffix array
}

// A wrapper function for histsort that simply prints the BWT as a string
char * makebwt(const char *str, int len) {
	int *idxs;
	char *buf = malloc(len+2); // +2, because we'd like it to be
	// null-terminated (to play nice with printf and so on)
	idxs = histsort(str, len);
	sprintbwt(buf, str, idxs, len);
	free(idxs);
	return buf;
}

// Make compressed bwt; this is for testing purposes, I need the indices
// for the FM index building
int makecbwt(const char *str, int len, char *out) {
	int *idxs, i;
	idxs = histsort(str, len);
	i = sprintcbwt(str, idxs, len, out);
	free(idxs);
	return i;
}

// Same as above, but uses saca-k to calculate the SA
int saca_makecbwt(const char *str, int len, char *out) {
	int *idxs, i;
	idxs = csuff_arr(str, len); // WILL TOTALLY WORK WINK WINK NUDGE NUDGE
	i = sprintcbwt(str, idxs, len, out);
	free(idxs);
	return i;
}

// Now, that's not very useful for the real usage, so we can write a
// compressed version of the bwt instead
// Returns the position of the '$' character (which will not be printed
// in out). out is assumed to have at least len/4 bytes of space
// allocated. str is assumed to be in compressed form.
// len is the length of str, not idxs
int sprintcbwt(const char *str, int *idxs, int len, char *out) {
	int i, d;
	char c = 0, u=3;
	for (i=0; i<=len; ++i) {
		if (idxs[i]) {
			c ^= getbase(str, idxs[i]-1)<<(2*u);
			if (u-- == 0) {
				out[i/4] = c;
				c = 0;
				u = 3;
			}
		}
		else {
			d = i;
			break; // So we don't have to pollute our loop with
			// conditionals
		}
	}
	for (++i; i<=len; ++i) {
		c ^= getbase(str, idxs[i]-1)<<(2*u);
		if (u-- == 0) {
			out[(i-1)/4] = c;
			c = 0;
			u = 3;
		}
	}
	// "flush" the rest of the buffer
	if (u != 3)
		out[((len+3)/4)-1] = c;
	 /*
	for (i = 0; i < len; ++i) {
		if (i==d)
			putchar('$');
		putchar('0'+getbase(out, i));
	}
	putchar('\n'); */
	return d;
}

void putsg(const char *str, int idx, int len) {
	// Prints a rotation of str corresponding to idx
	int i;
	// len is the length of str
	printf("%3d: ", idx);
	for (i = idx; i != len; ++i) {
		putchar('0'+getbase(str, i));
	}
	putchar('$');
	for (i = 0; i != idx; ++i) {
		putchar('0'+getbase(str, i));
	}
	putchar('\n');
}

void sprintbwt(char *out, const char *str, const int *bwt, int len) {
	// Prints the bwt into a string; bwt is assumed to be the return
	// value of histsort, and out should be at least len+1 bytes long.
	// Note that len is the length of the bwt, not str.
	int i;
	for (i = 0; i != len; ++i)
		out[i] = bwt[i] ? ('0' + getbase(str, bwt[i]-1)) : '$';
	out[i] = 0;
}

void putbwt(const char *str, const int *bwt, int len) {
	// Prints the bwt; bwt is assumed to be the return value
	// of histsort
	int i;
	for (i = 0; i <= len; ++i)
		if(bwt[i])
			putchar('0' + getbase(str,bwt[i]-1));
		else
			putchar('$');
	putchar('\n');
}
