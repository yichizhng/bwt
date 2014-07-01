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
#include "smw.h"

static inline unsigned char getbase(const char *str, int idx) {
  if (idx<0) idx=0;
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
// Last argument is the difference between the return value and the number of nts on the genome matched (from -3 to 3).
int mms_mismatch(const fm_index *fmi, const char *seq, const char *pattern, int len, int *sp, int *ep, int *genomeskips) {
  // If there are too many matches, don't even bother
  //  if (*ep - *sp > 10)
  //    return -1;
  if (len < 2) { // nothing to do, really
    int loc = unc_sa(fmi, *sp);
    char sub_c = getbase(seq, loc-1);
    *sp = fmi->C[sub_c] + rank(fmi, sub_c, *sp);
    *ep = *sp + 1;
    *genomeskips = 0;
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
	*genomeskips = 0;
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
	*genomeskips = 1;
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
	*genomeskips = 2;
	break;
      }

      // three!
      sub_c = getbase(seq, loc-3);
      ins_idx = fmi->C[sub_c] + rank(fmi, sub_c, blah);
      ins_align = mms_continue(fmi, pattern, len, &ins_idx, &ins_end);
      if (ins_align > 5 || ins_align == len) {
	best_align = sub_align;
	best_pos = sub_idx;
	*genomeskips = 3;
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
	*genomeskips = -1;
	break;
      }

      del_idx = i;
      del_end = del_idx + 1;
      del_align = mms_continue(fmi, pattern, len-2, &del_idx, &del_end) + 2;
      if (del_align > 7 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*genomeskips = -2;
	break;
      }

      del_idx = i;
      del_end = del_idx + 1;
      del_align = mms_continue(fmi, pattern, len-3, &del_idx, &del_end) + 3;
      if (del_align > 8 || del_align == len) {
	best_align = del_align;
	best_pos = del_idx;
	*genomeskips = -3;
	break;
      }
    }
  }
  *sp = best_pos;
  *ep = best_pos + 1;
  return best_align;
}

// Pass in the required anchor length. No mismatch will be allowed.
int align_read_anchored(const fm_index *fmi, const char *seq, const char *pattern, int len, int anchor_len) {
  int olen = len;
  int anchmisses = len/5, nmisses;
  // Here we require an anchor to start in the last 20% of the read
  int curgap = 0;
  int curpos = -1;
  int endpos;
  int anchlen;
  // Look for an anchor of length at least anchor_len (try 20 or so, or maybe
  // log_4(fmi->len)+1)
  while (len > anchor_len && anchmisses > 0) {
    nmisses = 0;
    while ((len > anchor_len) && (anchmisses > 0)) {
      int seglen = mms(fmi, pattern, len, &curpos, &endpos);
      if (seglen < anchor_len || endpos - curpos > 1) {
	anchmisses--;
	len -= 2;
	continue;
      }
      else {
	len -= seglen;
	anchlen = seglen;
	nmisses = anchmisses + len/5;
	curpos = unc_sa(fmi, curpos);
	break;
      }
    }
    
    if (nmisses < 1)
      continue;

    // In the second loop we try to extend our anchor backwards
    while ((len > nmisses) && (len > 4) && (nmisses > 0)) {
      for (curgap = 1; curgap < 4; ++curgap) {
	int start, end;
	int seglen = mms(fmi, pattern, len-curgap, &start, &end);
	int matched = 0;
	for (int i = start; i < end; ++i) {
	  if (abs(unc_sa(fmi, i) + seglen - curpos) <= curgap) {
	    // TODO: write proper scoring function, the number of misses
	    // is not going to be curgap
	    nmisses -= curgap;
	    matched = 1;
	    curpos = unc_sa(fmi, i);
	    len -= seglen + curgap;
	    curgap = 0;
	    break;
	  }
	}
	if (matched)
	  break;
	else
	  continue;
      }
      if (curgap)
	nmisses = 0;
    }
    if (nmisses > 0) {
      if (len < 3)
	return curpos - len;
      // Set up matrix for N-W alignment
      char *buf = malloc(len + 10);
      for (int i = 0; i < len + 10; ++i)
	buf[i] = getbase(seq, curpos - i);
      char *buf2 = malloc(len);
      for (int i = 0; i < len; ++i)
	buf2[i] = pattern[len-1-i];
      int x = nw_fast(buf, len+10, buf2, len);
      free(buf);
      free(buf2);
      //printf("%d %d\t", x, len);
      return curpos - x;
    }

    len -= anchlen;
    anchmisses -= anchlen / 10;
  }
  if (len > nmisses || nmisses < 1)
    return 0;

  switch(len) { /*
  case 0:
  case 1:
  case 2:
    return curpos - len;
    break; */
    // Taking the mismatch is always cheaper than opening a gap (8+5 > 6+6)
  case 3:
    if ( (pattern[0] != getbase(seq, curpos-3)) &&
	 (pattern[1] != getbase(seq, curpos-2)) &&
	 (pattern[2] != getbase(seq, curpos-1))) {
      if (pattern[0] == getbase(seq, curpos-4) &&
	  pattern[1] == getbase(seq, curpos-3) &&
	  pattern[2] == getbase(seq, curpos-2)) {
	return curpos - len - 1;
      }
      else if ((pattern[0] == getbase(seq, curpos-2)) &&
	       (pattern[1] == getbase(seq, curpos-1)))
	return curpos - len - 1;
    }
    break;

  default:
    return curpos - len;
    break;
  }
}

int align_read(const fm_index *fmi, const char *seq, const char *pattern, int len, int thresh) {
  int starts[10], lens[10], nsegments;
  int penalty;
  int nmisses = len/10;
  int olen = len;
  for (nsegments = 0; nsegments < 10; nsegments++) {
    if (len < 10)
      break;
    int start, end;
    int seglen = mms(fmi, pattern, len, &start, &end);
    if (seglen < thresh) {
      int mlen = mms_mismatch(fmi, seq, pattern, len - seglen, &start, &end, &penalty);
      if (mlen + seglen > 2 * thresh) {
	len -= seglen + mlen + 3;
	starts[nsegments] = start;
	lens[nsegments] = seglen + mlen;
	continue;
      }
      if (!nmisses--)
	return 0;
      len -= 3;
      nsegments--;
      if (nsegments > -1) {
	starts[nsegments] -= 3;
	lens[nsegments] += 3;
      }
      continue;
    }
    if ((len - seglen == 0) || ((len - seglen > 10) && end - start == 1)) {
      starts[nsegments] = start;
      lens[nsegments] = seglen;
      len -= seglen + 3;
      continue;
    }
    // Otherwise try continuing the search
    int mlen = mms_mismatch(fmi, seq, pattern, len - seglen, &start, &end, &penalty);
    len -= seglen + mlen + 3;
    starts[nsegments] = start;
    lens[nsegments] = seglen + mlen;
  }
  int totlen = lens[0];
  if (nsegments == 10)
    return 0; // Too many segments
  else {
    // For each segment check whether it's within 6 nts of the next
    
    for (int i = 0; i < nsegments - 1; ++i) {
      if (abs(unc_sa(fmi, starts[i+1]) + lens[i+1] - unc_sa(fmi, starts[i])) < 7) {
	totlen += lens[i+1];
	continue;
      } 
      else
	return 0; // Gapped
    }
  }
  if (3 * totlen > 2 * olen)
    return unc_sa(fmi, starts[nsegments-1]) - len;
  return 0;
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
    if (!fgets(buf, 256*256-1, rfp))
      break;
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
    //    int thresh = (int) (-1.2 * (1+len));

    //    int pos = align_read(fmi, seq, buf, len, 10);
    int pos = align_read_anchored(fmi, seq, buf, len, 12);
    if (pos) {
      naligned++;
      printf("%d\n", pos + 1);
    }
    else {
      //      pos = align_read(fmi, seq, revbuf, len, 10);
      pos = align_read_anchored(fmi, seq, revbuf, len, 12);
      if (pos) {
	naligned++;
	printf("%d\n", pos + 1);
      }
      else
	printf("0\n");
    }

    /*
    while(len) {
      if (score <= thresh) {
	break;
      }
      int start, end;
      int matched = mms(fmi, buf, len, &start, &end);
      if (matched < 10) {
      	len -= 1;
      	score -= 3;
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
	score -= 3;
	continue;
      }
      else {
	// we're good
	printf("%d\n", unc_sa(fmi, start) + 1);
	aligned = 1;
	naligned++;
	break;
      }
    }

    if (!aligned) {
      // Try aligning as a reverse complement
      int score = 0;
      while(len) {
	if (score <= thresh) {
	  printf("0\n");
	  break;
	}
	int start, end;
	int matched = mms(fmi, revbuf, len, &start, &end);
	//	printf("Matched %d\n", matched);
	if (matched < 10) {
	  len -= 1;
	  score -= 3;
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
	  score -= 3;
	  continue;
	}
	else {
	  // we're good
	  printf("%d\n", unc_sa(fmi, start) + 1);
	  naligned++;
	  break;
	}
      }
    }
    */
  }
  fclose(rfp);
  fprintf(stderr, "%d of %d reads aligned\n", naligned, nread);
  
  free(buf);
  free(revbuf);
  destroy_fmi(fmi);
  free(seq);
  return 0;
}
