#include "histsort.h"
#include "sacak.h"
#include "rdtscll.h"
#include <stdlib.h>
#include <stdio.h>

// A more efficient sort... ish. It's O(n log n) because the string length
// also increases as we go, but quicksort is O(n log^2 n) for this, because
// the comparison function isn't constant time.
// TODO: optimize for space (Done, results in histsortcomp.c)

int main(int argc, char **argv) {
	/* Testing code (to verify correctness)
	char *str, **arr;
	int i;
	str = malloc(10 * sizeof(char));
	arr = malloc(10 * sizeof(char *));
	strcpy(str, "032310211");
	for (i = 0; i < 10; ++i) {
		arr[i] = str + i;
		str[i] -= '0'; // To get them in the right range
	}
	histsort(arr, 10);
	for (i = 0; i < 10; ++i) {
		str[i] += '0';
	}
	//for (i = 0; i < 10; ++i) {
	//	putsg(str, arr[i], 9);
	//}
	putbwt(str, arr, 10);
	return 0;
	*/
	// Performance testing code
	int len;
	char *str, **arr;
	int i;
	unsigned long long a, b;
	unsigned int *SA;
	len = atoi(argv[1]);
	str = malloc((len+1) * sizeof(char));
	arr = malloc((len+1) * sizeof(char *));
	for (i = 0; i <= len; ++i) {
		arr[i] = str + i;
		str[i] = rand() & 3;
	}
	str[len] = 0;
	rdtscll(a);
	SA = suff_arr(str, len);
	rdtscll(b);
	puts("Using SACA-K");
	printf("Processed %d nucleotides in %lld cycles (%f seconds)\n",
		len, (b-a), ((double)(b-a)) / 2500000000.);
	free(SA);
	rdtscll(a);
	histsort(arr, len+1);
	rdtscll(b);
	//for (i = 1; i < len; ++i) {
	//	if (SA[i] + str != arr[i]) {
	//		printf("Mismatch at index %d (%X %X)\n", i,
	//			SA[i]+str, arr[i]);
	//		break;
	//	}
	//}
	puts("Using histogram sort");
	printf("Processed %d nucleotides in %lld cycles (%f seconds)\n",
		len, (b-a), ((double)(b-a)) / 2500000000.);
	//if (i == len)
	//	printf("Sequences match :)\n");
	free(arr);
	free(str);
	return 0;
	// Doesn't print anything to screen for performance reasons
}
