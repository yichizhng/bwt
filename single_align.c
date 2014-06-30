// Tries aligning reads from a file against an index and sequence read from
// file, assuming that they are not spliced reads
// This, of course, requires that we put another function together.

// usage: single_align seqfile indexfile readfile

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"
#include "fileio.h"
#include "rdtscll.h"
#include "time.h"

static inline unsigned char getbase(const char *str, int idx) {
	// Gets the base at the appropriate index
	return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// Continues a MMS search
int mms_continue(const fm_index *fmi, const char *pattern, int len, int *sp, int *ep) {
  int start, end, i;
  start = *sp;
  end = *ep;
  for (i = len-1; i >= 0; --i) {
    if (end <= start) {
      break;
    }
    *sp = start;
    *ep = end;
    start = fmi->C[pattern[i]] + rank(fmi, pattern[i], start);
    end = fmi->C[pattern[i]] + rank(fmi, pattern[i], end);
  }
  if (end <= start) // Didn't finish matching
    return len - i - 2;
  else { // Finished matching
    *sp = start;
    *ep = end;
    return len - i - 1;
  }
}

// Tries continuing a mms search with mismatch; returns upon finding any continuation with at least 6 matching nts
int mms_mismatch(const fm_index *fmi, const char *seq, const char *pattern, int len, int *sp, int *ep, int *penalty) {
  // If there are too many matches, don't even bother
  //  if (*ep - *sp > 10)
  //    return -1;
  if (len < 2) { // nothing to do, really
    *penalty = -6;
    return 1;
  }
  int best_align = 0;
  int best_pos = -1;
  for (int i = *sp; i < *ep; ++i) {
    // Reads the start and end from sp and ep instead of using the last
    // character of the sequence. It assumes that we have a mismatch at that
    // point (mms returns if that happens or it finished)
    // and tries the following things to try aligning it
    
    // 1) Assume that there was a substitution at that point. Use LF() to skip
    // to the next nt and decrement len, then try aligning
    {
      int loc = unc_sa(fmi, i);
      char sub_c = getbase(seq, loc-1);
      int sub_idx = fmi->C[sub_c] + rank(fmi, sub_c, i), ins_idx = sub_idx;
      int sub_end = sub_idx + 1, sub_align;
      sub_align = mms_continue(fmi, pattern, len-1, &sub_idx, &sub_end) + 1;
      best_align = sub_align;
      best_pos = sub_idx;
      if (sub_align > 6 || sub_align == len) {
	*penalty = -6;
	break;
      }

      // 1.5) Assume that there was an insertion (on the genome) at that point of up to three nts
      // Use LF() to skip one, two, and three nts and _don't_ decrement len, then try aligning for each of those
      int bleh = ins_idx;

      int ins_end = ins_idx + 1, ins_align;
      ins_align = mms_continue(fmi, pattern, len, &ins_idx, &ins_end);
      if (ins_align > 5 || ins_align == len) {
	best_align = sub_align;
	best_pos = sub_idx;
	*penalty = -8;
	break;
      }

      // two!
      sub_c = getbase(seq, loc-2);
      ins_idx = fmi->C[sub_c] + rank(fmi, sub_c, bleh);
      int blah = ins_idx;
      ins_align = mms_continue(fmi, pattern, len, &ins_idx, &ins_end);
      if (ins_align > 5 || ins_align == len) {
	best_align = sub_align;
	best_pos = sub_idx;
	*penalty = -11;
	break;
      }

      // three!
      sub_c = getbase(seq, loc-3);
      ins_idx = fmi->C[sub_c] + rank(fmi, sub_c, blah);
      ins_align = mms_continue(fmi, pattern, len, &ins_idx, &ins_end);
      if (ins_align > 5 || ins_align == len) {
	best_align = sub_align;
	best_pos = sub_idx;
	*penalty = -14;
	break;
      }
    }
    
    // 2) Assume that there was a deletion (on the genome) at that point.
    // Ignore up to three nts and start aligning again
    {
      // This one is a lot simpler because we don't actually need to
      // figure out the character
      int del_idx = i, del_end = del_idx + 1, del_align;
      del_align = mms_continue(fmi, pattern, len-1, &del_idx, &del_end) + 1;
      if (del_align > 6 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*penalty = -8;
	break;
      }

      del_idx = i;
      del_end = del_idx + 1;
      del_align = mms_continue(fmi, pattern, len-2, &del_idx, &del_end) + 2;
      if (del_align > 7 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*penalty = -11;
	break;
      }

      del_idx = i;
      del_end = del_idx + 1;
      del_align = mms_continue(fmi, pattern, len-3, &del_idx, &del_end) + 3;
      if (del_align > 8 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*penalty = -14;
	break;
      }
    }
  }
  *sp = best_pos;
  *ep = best_pos + 1;
  return best_align;
}

// Reminder to self: buf length (i.e. maximum read length) is currently
// hardcoded; change to a larger value (to align longer reads) or make it
// dynamic

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s seqfile indexfile readfile\n", argv[0]);
    exit(-1);
  }
  char *seq, *seqfile, *indexfile, *readfile, *buf = malloc(256*256), *revbuf = malloc(256*256), c;
  fm_index *fmi;
  int len;
  int i, j, k, jj;
  FILE *sfp, *ifp, *rfp;
  seqfile = argv[1];
  indexfile = argv[2];
  readfile = argv[3];
  sfp = fopen(seqfile, "rb");
  if (sfp == 0) {
    fprintf(stderr, "Could not open sequence\n");
    exit(-1);
  }
  fseek(sfp, 0L, SEEK_END);
  len = ftell(sfp);
  rewind(sfp);
  seq = malloc(len/4+1);
  for (i = 0; i < len/4 + 1; ++i) {
    switch(fgetc(sfp)) {
    case 'C': c = 64; break;
    case 'G': c = 128; break;
    case 'T': c = 192; break;
    default: c = 0;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 16; break;
    case 'G': c ^= 32; break;
    case 'T': c ^= 48;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 4; break;
    case 'G': c ^= 8; break;
    case 'T': c ^= 12;
    }
    switch(fgetc(sfp)) {
    case 'C': c ^= 1; break;
    case 'G': c ^= 2; break;
    case 'T': c ^= 3;
    }
    seq[i] = c;
  }
  // Handle the last character (which is at seq[len/4]
  c = 0;
  for (i = 0; i < len&3; ++i) {
    switch(fgetc(sfp)) {
    case 'C': c ^= 64 >> (2 * i); break;
    case 'G': c ^= 128 >> (2 * i); break;
    case 'T': c ^= 192 >> (2 * i);
    }
    seq[len/4] = c;
  }
  fclose(sfp);
  
  // Open index file
  ifp = fopen(indexfile, "rb");
  if (ifp == 0) {
    fprintf(stderr, "Could not open index file");
    exit(-1);
  }
  fmi = read_index(seq, ifp);
  fclose(ifp);

  // And now we go read the index file
  rfp = fopen(readfile, "r");
  if (rfp == 0) {
    fprintf(stderr, "Could not open reads file");
    exit(-1);
  }
  // Read one line ("read") and try aligning it
  
  int naligned = 0;
  int nread = 0;
  while (!feof(rfp)) {
    // Align the read using mms and mms_mismatch (which is a sort of wrapper
    // for the correct calls to mms_continue)
    if (! fgets(buf, 256*256-1, rfp))
      continue;
    nread++;
    if (buf[strlen(buf)-1] == '\n')
      buf[strlen(buf)-1] = 0;
    int len = strlen(buf);
    for (int i = 0; i < len; ++i) {
      // Replace with "compressed" characters
      switch(buf[i]) {
      case 'A':
	buf[i] = 0;
	revbuf[len-i-1] = 3;
	break;
      case 'C':
	buf[i] = 1;
	revbuf[len-i-1] = 2;
	break;
      case 'T':
	buf[i] = 3;
	revbuf[len-i-1] = 0;
	break;
      case 'G':
	buf[i] = 2;
	revbuf[len-i-1] = 1;
	break;
      default:
	buf[i] = 2;
	revbuf[len-i-1] = 2;
	break;
	// TODO: handle 'N' correctly
      }
    }

    int aligned = 0;

    int score = 0;
    int thresh = (int) (-0.6 * (1+len)); // the same one bowtie2 uses
    while(len) {
      if (score <= thresh) {
	break;
      }
      int start, end;
      int matched = mms(fmi, buf, len, &start, &end);
      if (end - start > 10) {
      	len -= 1;
      	score -= 6;
      	continue;
      }
      // Try continuing from these results
      int res_len = len - matched;
      int tscore = score;
      while(res_len && tscore > thresh) {
	int penalty;
	int matched_cont = mms_mismatch(fmi, seq, buf, res_len, &start, &end, &penalty);
	//printf("%d\n", matched_cont);
	if (matched_cont == -1) {
	  tscore = thresh;
	  break; // too many matches
	}
	tscore += penalty;
	res_len -= matched_cont;
      }
      if (tscore <= thresh) {
	len -= 1;
	score -= 6;
	continue;
      }
      else {
	// we're good
	printf("%d %d\n", nread, unc_sa(fmi, start) + 1);
	aligned = 1;
	naligned++;
	break;
      }
    }

    if (!aligned) {
      // Try aligning as a reversed antisense strand
      int score = 0;
      int thresh = (int) (-0.6 * (1+len)); // the same one bowtie2 uses
      while(len) {
	if (score <= thresh) {
	  printf("%d 0\n", nread);
	  break;
	}
	int start, end;
	int matched = mms(fmi, revbuf, len, &start, &end);
	//	printf("Matched %d\n", matched);
	if (end - start > 10) {
	  len -= 1;
	  score -= 6;
	  continue;
	}
	// Try continuing from these results
	int res_len = len - matched;
	int tscore = score;
	while(res_len && tscore > thresh) {
	  int penalty;
	  int matched_cont = mms_mismatch(fmi, seq, revbuf, res_len, &start, &end, &penalty);
	  //printf("%d\n", matched_cont);
	  if (matched_cont == -1) {
	    tscore = thresh;
	    break; // too many matches
	  }
	  tscore += penalty;
	  res_len -= matched_cont;
	}
	if (tscore <= thresh) {
	  len -= 1;
	  score -= 6;
	  continue;
	}
	else {
	  // we're good
	  printf("%d %d\n", nread, unc_sa(fmi, start) + 1);
	  naligned++;
	  break;
	}
      }
    }
  }
  fclose(rfp);
  fprintf(stderr, "%d of %d reads aligned\n", naligned, nread);
  
  free(buf);
  free(revbuf);
  destroy_fmi(fmi);
  free(seq);
  return 0;
}
