// Implements a histogram sort (i.e. a bucket sort in O(n) auxiliary space)
// (It can be done in O(1) auxiliary space, but at the cost of speed; this
// seems unnecessary in the modern world with even "modest" servers having
// 128GB of memory)
// which is applied to the suffixes of a string, which is assumed to only
// contain (char)0, (char)1, (char)2, and (char)3.

// The idea here is that we don't want to use quicksort here, because it's
// asymptotically slow; there are theoretically O(n) algorithms (e.g.
// DC3) but they're 1) complicated 2) bad cache locality 3) hard to parallelize

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "histsort.h"

static char **saux;

void histhelper(char **arr, char **aux, char *e, int start, int end, int depth) {
	// We are responsible for sorting the array from arr[start] to
	// arr[end-1], and are currently considering the character at
	// position depth.
	
	// There are four buckets, plus one for a string which has just
	// ended. Each of those four buckets is then sorted recursively.
	int lens[4] = {0}, i, t;
	char **ptrs[4];
	// Base cases, wherein the array is already sorted
	if ((end - start) == 1) {
		if (arr == saux)
			aux[start] = arr[start];
			// We need to make sure the original array is sorted
		return;
	}
	if (end == start)
		return;
	for (i = start; i != end; ++i) {
		if(arr[i] + depth != e) {
			lens[arr[i][depth]]++;
		}
		else {
			aux[start] = arr[i];
			arr[i] = arr[start];
			arr[start] = aux[start];
			++start; // This element is known to be in the right
			// position now, so skip on further iterations
		}
	}
	// lens will now contain the numbers of elements which should
	// go into each bucket
	ptrs[0] = aux + start;
	ptrs[1] = ptrs[0] + lens[0];
	ptrs[2] = ptrs[1] + lens[1];
	ptrs[3] = ptrs[2] + lens[2];
	// I don't know why I kept lens[3] around either, maybe a sanity
	// check (or to avoid a branch prediction error). It lets me
	// find the end of the array, but I know where that is anyway
	for (i = start; i != end; ++i) {
		// None of the strings have ended (because we've taken
		// care of that already) so we can just use arr[i][depth]
		// to index
		*ptrs[arr[i][depth]] = arr[i];
		++ptrs[arr[i][depth]];
	}
	// And now recursively call histhelper on each of the four buckets
	histhelper(aux, arr, e, start, start + lens[0], depth+1);
	histhelper(aux, arr, e, start+lens[0], start+lens[0]+lens[1],depth+1);
	histhelper(aux, arr, e, start+lens[0]+lens[1], end - lens[3],depth+1);
	histhelper(aux, arr, e, end - lens[3], end, depth+1);
}

void ippass(char **arr, int start, int end, int depth) {
	// Does a bucketing of arr from start to end in place
	// 2-pass, because a 1-pass version of this would be pretty complicated
	// and have quite a few swaps (and that would be bad-ish (not really))
}

void histsort(char **arr, int len) {
	// We assume that the array is currently in descending order of
	// length (that is, the elements are str, str+1, str+2, ..., str+len-1)
	char *end = (*arr) + len - 1; // Because I don't want to waste memory
	// adding the $ symbol, we simply keep track of where the array
	// end should be
	char **aux = malloc(len * sizeof(char *));
	saux = aux; // So we know which array we want the data to be in
	// Allows us to flip-flop the data between two arrays
	histhelper(arr, aux, end, 0, len, 0);
	// By the power of magic and handwaving, arr will end up sorted
	free(aux);
}

void putsg(char *base, char *str, int len) {
	// Does some trickery to print things correctly
	// This makes some pretty specific assumptions about the string
	// namely that it has a copy of itself appended to the end
	// The value of len passed here should be 1 less than elsewhere
	// which I guess is an indexing error on my part :)
	int i, m;
	m = len - strlen(str);
	printf("%s", str);
	putchar('$');
	for(i = 0; i < m; ++i)
		putchar(base[i]);
	putchar('\n');	
}

void putbwt(char *base, char **arr, int len) {
	// arr is assumed to be a sorted prefix list (which points into base)
	// len should actually be strlen(base) + 1.
	// This function is fully general (i.e. does not assume that the 
	// elements of base consist of 0, 1, 2, and 3) and so can be called
	// after the list is processed to put it back into, say, ACGT form.
	int i;
	for (i = 0; i < len; ++i) {
		if (arr[i] == base)
			putchar('$');
		else
			putchar(arr[i][-1]);
	}
	putchar('\n');
}
