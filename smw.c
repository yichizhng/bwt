// Slowish implementation of Smith-Waterman algorithm for local alignments
// There are optimizations involving SIMD instructions but they don't really
// seem that worthwhile
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "rdtscll.h"

static inline int max(int a, int b, int c) {
  return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

// TODO: more typedefs? I don't personally find C's default types confusing,
// but a lot of "production" code I've seen is full of it. Maybe they're too
// dependent on IDEs?

// TODO: move this struct into another file, most likely. I have a lot of
// reused / reusable code which I should really reorganize
// A struct for returning information relating to a gap stitching
struct gap_stitch {
  int len1; // Number of nucleotides matched onto the first part
  // Note that this includes indels on both sides!
  int len2; // Analogous
  // Hint: strlen(pattern) == strlen(gen_pat) == len1 + len2 + 1
  char *pattern; // Printable representation of the alignment of the pattern
  // '/' is used to denote the gap
  char *gen_pat; // Parts of the genome that the pattern was mapped onto
};

struct gap_stitch *nw_stitch(const char *pat, int pat_len, const char *gen_l,
			     const char *gen_r, int gen_len) {
  int *valuesl, *valuesr, i;
  char *patrev, *gen_rrev;
  // We need to reverse pat and gen_r to use them as input for nw_gap_l
  patrev = malloc(pat_len);
  gen_rrev = malloc(gen_len);
  for (i = 0; i < pat_len; ++i) {
    patrev[i] = pat[pat_len - i];
  }
  for (i = 0; i < gen_len; ++i) {
    gen_rrev[i] = gen_r[gen_len - i];
  }
}

// Aligns from the left; obviously we can reuse this to align from the right
// just by reversing everything
// To be used as a component in aligning pretty much everything.
int *nw_gap_l(const char *pat, int pat_len, const char *gen, int gen_len) {
  // gen_len should probably only be a bit larger than pat_len; we can't
  // have *that* many indels
  
  // The main difference is that the return array now has two extra columns
  // which denote the highest score found in each row and the location of
  // that score (for backtracking purposes)
  
  int *values, i, j, m, m_pos;
  values = malloc((pat_len + 1) * (gen_len + 3) * sizeof(int));
  // Initialize first row
  for (j = 0; j <= gen_len; ++j) {
    values[j] = -j;
  }
  values[gen_len + 1] = 0;
  values[gen_len + 2] = 0;

  // Deal with rest of matrix
  for (i = 1; i <= pat_len; ++i) {
    // Initialize first column
    values[i * (gen_len + 3)] = -i;
    m = -i;
    m_pos = 0;
    for (j = 1; j <= gen_len; ++j) {
      // Update cell appropriately
      values[i * (gen_len + 3) + j] = 
	max(values[(i-1) * (gen_len + 3) + j - 1] +
	    ((pat[i-1] == gen[j-1])?2:-1),
	    values[i * (gen_len + 3) + j - 1] - 1,
	    values[(i-1) * (gen_len + 3) + j] - 1);
      if (values[i * (gen_len + 3) + j] > m) {
	m = values[i * (gen_len + 3) + j];
	m_pos = j;
      }
    }
    // Write the appropriate m and m_pos
    values[i * (gen_len + 3) + gen_len + 1] = m;
    values[i * (gen_len + 3) + gen_len + 2] = m_pos;
  }
  return values;
}

// Optimization of needleman-wunsch by using a 1d output array; is faster by a
// pretty hilarious factor (>40x) due to the lack of another level of
// indirection and/or cache optimizations and/or lack of malloc() calls
int nw_fast(const char *str1, int len1, const char *str2, int len2) {
  int *values, i, j;
  int mx = -1000000, maxloc = 0;
  values = malloc((len1 + 1) * (len2 + 1) * sizeof(int));
  // "Zero" first row
  values[0] = 0;
  for (j = 1; j <= len2; ++j) {
    values[j] = -5-3*j;
  }
  for (i = 1; i <= len1; ++i) {
    // Zero first column
    values[i * (len2 + 1)] = -5-3*i;
    for (j = 1; j <= len2; ++j) { /*
      int skip1 = (((values[(i-1) * (len2 + 1) + j - 1] - 
		     values[(i-2) * (len2 + 1) + j - 1]) == -3) ||
		   ((values[(i-1) * (len2 + 1) + j - 1] - 
		     values[(i-2) * (len2 + 1) + j - 1]) == -8)) ? 0 : -5;
      int skip2 = (((values[(i-1) * (len2 + 1) + j - 1] - 
		     values[(i-1) * (len2 + 1) + j - 2]) == -3) ||
		   ((values[(i-1) * (len2 + 1) + j - 1] - 
		   values[(i-1) * (len2 + 1) + j - 2]) == -8)) ? 0 : -5; */
      int skip1 = 0, skip2 = 0;
      if (i > 1 && j>1) {
	skip1 = ((values[(i-1) * (len2 + 1) + j - 1] - 2 ==
		  values[(i-2) * (len2 + 1) + j - 2]) ||
		 (values[(i-1) * (len2 + 1) + j - 1] + 3 ==
		  values[(i-1) * (len2 + 1) + j - 2]) ||
		 (values[(i-1) * (len2 + 1) + j - 1] + 8 ==
		  values[(i-1) * (len2 + 1) + j - 2])) ? -5 : 0;
	skip2 = ((values[(i-1) * (len2 + 1) + j - 1] - 2 ==
		  values[(i-2) * (len2 + 1) + j - 2]) ||
		 (values[(i-1) * (len2 + 1) + j - 1] + 3 ==
		  values[(i-2) * (len2 + 1) + j - 1]) ||
		 (values[(i-1) * (len2 + 1) + j - 1] + 8 ==
		  values[(i-2) * (len2 + 1) + j - 1])) ? -5 : 0;
      }
      // Update cell appropriately
      values[i * (len2 + 1) + j] =
	max(values[(i-1) * (len2 + 1) + j - 1] + ((str1[i-1] == str2[j-1])?2:-6),
	    values[i * (len2 + 1) + j - 1] - 3 + skip1,
	    values[(i-1) * (len2 + 1) + j] - 3 + skip2);
      if ((j == len2) && values[i * (len2 + 1) + j] > mx) {
	mx = values[i * (len2 + 1) + j];
	maxloc = i;
      }
    }
  }
  free(values);
  return maxloc - 1;
}

// Note that this implementation takes the full O(m*n) memory; it is possible
// to do with less (especially if we only want the optimal
// alignment), but much more annoying
// The slow and obnoxious implementation, returns a 2D array with results
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

/*

// Testing routine; should probably be removed or moved to a different file
int main(int argc, char **argv) {
  // Take two inputs from arguments
  long long int a, b;
  int **val, i, j, k, ii, *fval;
  char *buf1, *buf2;
  if (argc < 3)
    return -1;
  rdtscll(a);
  val = smw(argv[1], strlen(argv[1]),argv[2], strlen(argv[2]));
  rdtscll(b);
  // Do something with a and b :)
  printf("%f\n", (double)(b-a) / ((double)(strlen(argv[1]) * strlen(argv[2]))));
  // Print the optimal global alignment; we can do this by recursing backwards
  // and storing (or better, by reversing the string before aligning!)

  rdtscll(a);
  fval = nw_fast(argv[1], strlen(argv[1]), argv[2], strlen(argv[2]));
  rdtscll(b);
  printf("%f\n", (double)(b-a) / ((double)(strlen(argv[1]) * strlen(argv[2]))));


  // Note that backtracking on a Needleman-Wunsch matrix is easier than
  // on a Smith-Waterman, because we don't have to deal with all the zeroes
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

  for (i = 0; i <= strlen(argv[1]); ++i) {
    for (j = 0; j <= strlen(argv[2]); ++j) {
      printf("%3d", fval[i * (1+strlen(argv[2])) + j]);
    }
    putchar('\n');
  }
  free(fval);
}

*/
