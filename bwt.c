#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int sort_len;

static int compar(const void * a, const void * b) {
	// In essence, strcmp()
	/*
	int i, j, c;
	for (i = 0, j = 0; i < sort_len; ++i, ++j) {
		c = (*(char **)a)[i] - (*(char **)b)[j];
		if (c)
			return c;
	}
	return 0; */
	return strcmp(*(char **)a, *(char **)b);
	// Note that *a and *b are overlapping strings; however, this
	// does not really matter for our purposes.
}

char * bwt (char *in, int len) {
	// Performs the Burrows-Wheeler transform on the string in
	// which is len bytes long. Returns a newly allocated string
 	// containing the result.

	// in[len-1] should be (char)0; furthermore, (char)0 should not
	// appear anywhere else in in, unless you like undefined results.

	// TODO: Reimplement in linear (or at worst linearithmetic) time; most
	// likely by calculating the suffix array (SA-IS, SACA-K); we can try some
	// sort of modified bucket sort though, with fiveish buckets? (four
	// if we're clever about indexing) for a easily parallelized
	// linearithmetic version
	int i;
	char **c = malloc(len * sizeof(char *)), *tmp = malloc(len*2+1), *out;
	sort_len = len;
	// The index-building will be atomic and single-threaded in any case
	// So we don't have to worry about things like race conditions
	
	// Copy the string and a copy of it so that we can index into this
	// new string to get the rotations
	memcpy(tmp, in, len);
	tmp[len] = 0;
	memcpy(tmp+len+1, in, len);
	for (i = 0; i < len; ++i)
		c[i] = tmp+i;
	qsort(c, len, sizeof(char *), compar);
	out = malloc(len+1);
	for (i = 0; i < len; ++i) {
		printf("%d\n", c[i]-tmp);
		out[i] = c[i][len-1];
		}
	free(c);
	free(tmp);
	return out;
}

int main(int argc, char **argv) {	
	char *str, *out;
	int i, len;
	len = strlen(argv[1]);
	str = malloc(len+2);
	strcpy(str+1, argv[1]);
	//len = atoi(argv[1]);
	//str = malloc(len+2);
	//for(i = 1; i < len+1; ++i)
	//	str[i] = '0'+(rand()&3);
	*str = 0; // Prepend a EOF. This will be sorted before any other
	// character
	out = bwt(str, len+1);
	// Replace the EOF with a $ for output purposes
	for (i = 0; i < len+1; ++i)
		if (!out[i]) {
			out[i] = '$';
			break;
		}
	puts(out);
	free(str);
	free(out);
	return 0;
}
