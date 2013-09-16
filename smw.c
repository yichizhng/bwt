// Slow implementation of Smith-Waterman algorithm for local alignments
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
static inline int max(int a, int b, int c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

// Note that this implementation takes the full O(m*n) memory; it is possible
// to do with less, but much more annoying
int** smw(const char *str1, int len1, const char *str2, int len2) {
  // Allocate a 2-D array for values
  int **values, i, j;
  values = malloc((len1+1) * sizeof(void *));
  for (i = 0; i <= len1; ++i) {
    values[i] = malloc((len2+1) * sizeof(int));
    values[i][0] = -i;
  }
  for (j = 1; j <= len2; ++j) {
    values[0][j] = -j;
  }
  // The first row and column of the matrix are 0's, and have already
  // been assigned in the previous loop
  
  // For now we use +2 for a match and -1 for indel or mismatch
  
  for (i = 1; i <= len1; ++i) {
    for (j = 1; j <= len2; ++j) {
      if (str1[i-1] == str2[j-1]) {
	values[i][j] = 2 + values[i-1][j-1];
	continue;
	// Taking a match is always better than taking an indel if
	// we use linear gap penalty (instead of affine)
      }
      // Find the maximum of the three possibilities and subtract 1
      values[i][j] = max(values[i][j-1], values[i-1][j], values[i-1][j-1]) - 1;
      //if (values[i][j] < 0) // uncomment to use smith-waterman
      //values[i][j] = 0; // keep commented to use needleman-wunsch
      // smith-waterman is harder to backtrack, if you're wondering
      // what the problem is, so it's better to just pass appropriate length
      // strings in to get local alignments that way
    }
  }
  // In order to retrieve the correct path we need to iterate backward from
  // (i, j), but we can leave that for a different function
  return values;
}

int main(int argc, char **argv) {
  // Take two inputs from arguments
  int **val, i, j, k, ii;
  char *buf1, *buf2;
  if (argc < 3)
    return -1;
  val = smw(argv[1], strlen(argv[1]),argv[2], strlen(argv[2]));
  // Print the optimal local alignment; we can do this by recursing backwards
  // and storing (or better, by reversing the string before aligning!)
  i = strlen(argv[1]);
  j = strlen(argv[2]);
  buf1 = malloc(i + j + 1);
  buf2 = malloc(i + j + 1);
  buf1[i+j] = 0;
  buf2[i+j] = 0;
  k = 0;
  // This is as long as we can possibly need; it is possible to calculate
  // the number of indels first, but that hardly seems necessary ;)
  
  // Complicated loop logic to print the thing properly -.-
  // There are a lot of pitfalls in string indexing to be avoided here
  while (i && j) {
    // We can check whether we've taken a match (or *can* take a match) in order
    // to get to an optimal alignment (note that there can be multiple)
    // Since we would take a mismatch the same way, we can deal with that too
    if (((argv[1][i-1] == argv[2][j-1]) && (val[i][j] - val[i-1][j-1] == 2)) ||
	val[i][j] - val[i-1][j-1] == -1) {
      // We take a (mis)match here; increment k and print the character to both
      // strings
      --i;
      --j;
      // The point of incrementing k is to not break string indexing
      ++k;
      buf1[i+j+k] = argv[1][i];
      buf2[i+j+k] = argv[2][j];
    }
    else if (val[i][j] - val[i][j-1] == -1) {
      // Then we only decrement j (and don't increment k)
      --j;
      buf1[i+j+k] = '-';
      buf2[i+j+k] = argv[2][j];
    }
    else if (val[i][j] - val[i-1][j] == -1) {
      // Decrement i
      --i;
      buf1[i+j+k] = argv[1][i];
      buf2[i+j+k] = '-';
    }
    else {
      // crash and burn
      exit(-1);
    }
  }

  // If there's an indel at the beginning we get rid of that now
  // Do it in two loops just to simplify logic
  ii = i;
  while (i) {
    --i;
    buf1[i+j+k] = argv[1][i];
    buf2[i+j+k] = '-';
  }
  i = ii;
  while (j) {
    --j;
    buf1[i+j+k] = '-';
    buf2[i+j+k] = argv[2][j];
  }
  // Finally, print both buffers, using k to figure out where to start
  printf("%s\n%s\n", buf1+k, buf2+k);
  free(buf1);
  free(buf2);

  for (i = 0; i <= strlen(argv[1]); ++i) {
    for (j = 0; j <= strlen(argv[2]); ++j) {
      printf("%3d", val[i][j]);
    }
    putchar('\n');
    free(val[i]);
  }
  free(val);
}
