#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int sort_len;

static int compar(const void * a, const void * b) {
	return strcmp(*(char **)a, *(char **)b);
	// Note that *a and *b are overlapping strings; however, this
	// does not really matter for our purposes.
}

char * bwt (const char *in, int len) {
  // Performs the Burrows-Wheeler transform on the string in
  // which is len bytes long. Returns a newly allocated string
  // containing the result.
  // in containing a zero char will cause undefined results (well, to be
  // exact, they will be defined, just not useful)
  int i;
  char **c = malloc(len * sizeof(char *)), *tmp = malloc(len*2+1), *out;
  sort_len = len;
  // Copy the string and a copy of it so that we can index into this
  // new string to get the rotations
  memcpy(tmp, in, len);
  //tmp[len] = 0;
  memcpy(tmp+len, in, len);
  for (i = 0; i < len; ++i)
    c[i] = tmp+i;
  // Use library quicksort to compute suffix array
  qsort(c, len, sizeof(char *), compar);
  out = malloc(len+1);
  for (i = 0; i < len; ++i) {
    out[i] = c[i][len-1];
  }
  free(c);
  free(tmp);
  return out;
}

int main(int argc, char **argv) {	
	char *str, *out;
	int i, len;
	if (argc == 1)
	  return 0; // Not putting up with that crap
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
