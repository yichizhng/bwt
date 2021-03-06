#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "rdtscll.h"
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"
#include "smw.h"

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// Calculates floor(log_4(x)) + 1; as far as I'm concerned, this is the same
// as a ceiling (and will lead to one of the oddest edge cases if I leave it
// like this :])
int log_4(int x) {
  union {
    int i;
    float f;
  } blah;
  blah.f = (float)(2*x);
  return ((blah.i >> 24) - 63);
}

// GCC suggests that static inline functions are as fast as macros
// 
static inline int fnw_acc(int *values, int i, int j, int width) {
  // Width is len2 in the construction of values
  return values[i * (width + 1) + j];
}

// Performs a maximum mappable suffix search; returns the position of the
// suffix on the genome if the length found is at least equal to the cutoff,
// -1 otherwise; stores the length of the suffix in len_p. Returns an arbitrary
// (in fact, the first in the FM-index) match if multiple ones are found
// (TODO: fix that if necessary)

// TODO: write a more specialized version for "continuing" searches (e.g. take
// a gap of 3 on either side); stitching is trivial by comparison
int mms_search(const fm_index *fmi, const char *pattern, int len,
	       int *len_p, int cutoff) {
  // Structurally almost identical to locate(), except for the way we use
  // len...
  int start, end, i;
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      // Then the longest mappable suffix starts at the next nucleotide
      *len_p = len - i + 1;
      if (len - i > cutoff) {
	return unc_sa(fmi, start) + 1;
      }
      return -1;
    }
    start = fmi->C[pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] + rank(fmi, pattern[i], end);
  }
  // Matched the entire pattern; return the first match (note that this
  // is not likely to be useful; we need some notice of where the previous
  // match was to limit our possibilities, or maybe calculate SA^-1 to
  // allow us to continue a search from elsewhere)
  *len_p = len;
  return unc_sa(fmi, start);
}

// Another variant of locate() (actually more a variant of mms_search());
// this one takes into account the position (on the genome) of the previous
// MMS found and tries to continue it if possible; to compensate for this
// we allow a much smaller value for cutoff (noting that the relation of
// cutoff to statistical significance is exponential!)
int mms_continue(const fm_index *fmi, const char *pattern, int len,
		 int *len_p, int cutoff, int lastpos) {
  // Tries to continue a MMS search with the knowledge of the previous position;
  // this allows us to (usually) continue the search if we have some mismatch
  // or indel rather than a splice site
  int start, end, i, j, pos;
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      *len_p = len - i + 1;
      if (len - i > cutoff) {
	// Check every current match; return one which is before lastpos if
	// possible, no match otherwise (this may seem harsh, but meh to you
	// too. This is really unlikely to matter unless we have highly
	// repetitive/similar sequences in the genome)
	for (j = start; j < end; ++j) {
	  pos = unc_sa(fmi, j) + 1;
	  if (lastpos - (pos + (len - i)) > 3) {
	    // We do allow for a small overlap, since we can fit that into
	    // the stitching routine
	    return pos;
	    // This has the potential to be wrong, but it's fairly unlikely.
	    // (i.e. we might skip over the correct match to get the wrong one)
	  }
	}
      }
      return -1;
    }
    if (len - i == cutoff) {
      // Check all currently valid matches. If one is "near" lastpos (which
      // I'll hardcode as being within 6 nucleotides, or 2 codons, allowing a
      // fairly large tolerance in terms of mismatches), use it (and ignore all
      // other matches); otherwise continue the search as normal.
      
      for (j = start; j < end; ++j) {
	// Check the position of j
	pos = unc_sa(fmi, j);
	fprintf(stderr, "%d\n", lastpos - (pos+cutoff));
	if ((pos < lastpos) && (lastpos - (pos + cutoff) <= 6)) {
	  start = j;
	  end = j+1;
	  printf("Should be right?\n"); // Hmm.
	  break;
	}
      }
    }
    start = fmi->C[pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] + rank(fmi, pattern[i], end);
  }
  *len_p = len;
  return unc_sa(fmi, start);
}

// Another variant of locate(); this one does not try to continue the previous
// search, so this has a loop removed from mms_continue.
int mms_gap(const fm_index *fmi, const char *pattern, int len,
		 int *len_p, int cutoff, int lastpos) {
  // Tries to continue a MMS search with the knowledge of the previous position;
  // We abandon explicitly isolating the match that continues the previous
  // search (one interesting "optimization" to try would be to isolate the
  // match closest to the previous search and see whether that improves things)
  int start, end, i, j, pos;
  int score; // Keep a running tally of the score of the alignment so far
  // TODO: borrow scoring ideas from bowtie or something
  start = fmi->C[pattern[len-1]];
  end = fmi->C[pattern[len-1]+1];
  for (i = len-2; i >= 0; --i) {
    if (end <= start) {
      *len_p = len - i + 1;
      if (len - i > cutoff) {
	for (j = start; j < end; ++j) {
	  pos = unc_sa(fmi, j) + 1;
	  if (lastpos - (pos + (len - i)) > 3) {
	    return pos;
	  }
	}
      }
      return -1;
    }
    start = fmi->C[pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] + rank(fmi, pattern[i], end);
  }
  *len_p = len;
  return unc_sa(fmi, start);
}

// "Traditional" gapped alignment search
// Unidirectional search with no real tricks; meant to be fast, but doesn't
// guarantee best alignment or finding an alignment if one exists
// Also supports ungapped pretty much by default, it's probably better at that
// than you would think

// TODO: helper functions to stitch on the "head" and "tail" of the search
// And some stuff to help print all the crap this is going to produce
// It seems unlikely that a given strand of RNA will span more than, say,
// 10 exons, but who knows? (Hardcoded limits are for systems people)
void rna_seq(const fm_index *fmi, const char *pattern, int len) {
  // Use maximum mappable *suffix* (generated via backward search, as opposed
  // to binary search, which is slow) combined with the Needleman-Wunsch
  // algorithm (for stitching purposes mostly) to align RNA sequences

  // Motivation: STARS uses maximum mappable prefixes generated via binary
  // search on the FM-index; this is obviously silly (and O(m log(n))), and
  // so we're better off doing a backwards search to get (functionally
  // equivalent) maximum mappable suffixes instead to get O(m + log(n)) time
  // We could just reverse the genome to get equivalent results anyway
  
  // TODO: Adding a reverse FM-index may allow us to double anchor the search
  // for better accuracy? (Particularly if we have paired reads or something)

  int i, mmslen, mmspos, genpos, nextpos;
  // We begin indexing from the end of the pattern. Reverse search to find
  // a "statistically significant" match (len - log_4(fmi->len) > 2, for
  // example, makes a decent cutoff). Use needleman-wunsch or a specialized
  // variation (hint: it's possible to generate statistics on the highest
  // score per row and where it's located, which lets us do a stitch alignment)
  i = len;
  // TODO: replace 14 with some appropriate expression. I've written the log4
  // function (by abusing unions in a somewhat non-portable way). log2 is
  // apparently only in C11 (and not implemented in gcc!); check your local
  // version of tgmath.h to see if you have a library implementation
  mmspos = mms_search(fmi, pattern, i, &mmslen, 14);
  while ((mmspos == -1) && i > 14) {
    --i;
    mmspos = mms_search(fmi, pattern, i, &mmslen, 14);
    // We need *somewhere* to start...
  }
  // Now that we have a starting position, "stitch" the stuff behind it if
  // necessary
  // TODO: write the helper function to do that
  i -= mmslen; // LOL forgot that.
  //fprintf(stderr, "%d\n", i);
  while (i > 18) { // There's kind of a point where we should just give up
    genpos = mmspos;
    // Skip ahead 3 nucleotides (i.e. 1 codon; this deals with deletions of
    // entire codons, which is pretty common since it's less likely to cause
    // nonsense errors); it also allows us to catch shorter runs of
    // mismatches should that be necessary
    // We could probably skip by a smaller amount, but I foresee problems with
    // either too much or too little, and this is something that could probably
    // be tuned, or perhaps given as a command line option
    i -= 3;
    // Try continuing the search from before
    nextpos = mms_continue(fmi, pattern, i, &mmslen, 10, genpos);
    if (nextpos != -1) {
      // TODO: Stitch the matches as appropriate
      i -= mmslen;
    }
    else {
      // Assume that there's a gap (there might not be (multiple indels/
      // mismatches near each other) but who knows?)
      while (i > 14) {
	--i;
	// Try to align starting here
	nextpos = mms_gap(fmi, pattern, i, &mmslen, 14, genpos);
	if (nextpos != -1) {
	  // TODO: do something with this
	  i -= mmslen;
	  break;
	}
      }
    }
    // Main loop of function
  }
  // TODO: Do something at the end of the thingy
  printf("%d ", mmspos);
  // TODO: something more useful
}



// A sort of real-life test to see how effective we are at actually searching
// patterns over a sequence

int main(int argc, char **argv) {
  // We take our input filename from argv
  int len, i, j, k, jj;
  char *seq, *buf;
  unsigned char c;
  long long a, b;
  fm_index *fmi;
  FILE *fp;
  if (argc == 1) {
    printf("Usage: searchtest seq_file");
    exit(-1);
  }
  fp = fopen(argv[1], "rb");
  fseek(fp, 0L, SEEK_END);
  len = ftell(fp);
  fseek(fp, 0L, SEEK_SET); // This is *technically* not portable
  seq = malloc(len/4+1);
  // Complicated crap to take input and put it into compressed form
  // I should probably refactor this code at some point, it seems useful
  // The inverse function is much easier to write thankfully
  for (i = 0; i < len/4 + 1; ++i) {
    switch(fgetc(fp)) {
    case 'C': c = 64; break;
    case 'G': c = 128; break;
    case 'T': c = 192; break;
    default: c = 0;
    }
    switch(fgetc(fp)) {
    case 'C': c ^= 16; break;
    case 'G': c ^= 32; break;
    case 'T': c ^= 48;
    }
    switch(fgetc(fp)) {
    case 'C': c ^= 4; break;
    case 'G': c ^= 8; break;
    case 'T': c ^= 12;
    }
    switch(fgetc(fp)) {
    case 'C': c ^= 1; break;
    case 'G': c ^= 2; break;
    case 'T': c ^= 3;
    }
    seq[i] = c;
  }
  // Handle the last character (which is at seq[len/4]
  c = 0;
  for (i = 0; i < len&3; ++i) {
    switch(fgetc(fp)) {
    case 'C': c ^= 64 >> (2 * i); break;
    case 'G': c ^= 128 >> (2 * i); break;
    case 'T': c ^= 192 >> (2 * i);
    }
    seq[len/4] = c;
  }
  fclose(fp);
  // Now that we've loaded the sequence (ish) we can build an fm-index on it
  fmi = make_fmi(seq, len);
  // Do some fun tests (load up a length 30 sequence (starting from anywhere
  // on the "genome") and backwards search for it on the fm-index (and we're
  // going to fix locate() now too)
  buf = malloc(50); // The C/C++ standard guarantees that sizeof(char) == 1
  srand(time(0));
  rdtscll(a);
  for (i = 0; i < 1000000; ++i) {
    // Pick some randomish location to start from (i.e. anywhere from 0
    // to len-21)
    j = rand() % (len-50);
    for (k = 0; k < 50; ++k) {
      buf[k] = getbase(seq, j+k);
    }
    // Now we raaandomly screw up one nucleotide in the middle
    k = 20 + (rand() % 10);
    buf[k]^=(char)3;
    rna_seq(fmi, buf, 50);
    //jj = locate(fmi, buf, 50);
    printf("%d\n", j);
  }
  rdtscll(b);
  fprintf(stderr, "Took %lld cycles to search 1000000 50bp sequences (twice)\n",
	 b-a);
  fprintf(stderr, "(%f seconds), over a genome of length %d\n", 
	 ((double)(b-a)) / 2500000000, len);
  destroy_fmi(fmi);
  free(seq);
  free(buf);
  return 0;
}
