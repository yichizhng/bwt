#ifndef _SEQINDEX_H
#define _SEQINDEX_H

// The function to build the sequence index are here, as are the functions
// relating to the actual FM-index, as well as the struct definition thereof


int **seq_index(char *, int, int, const unsigned char *);

int seq_rank(unsigned char *, int **, int, int, char, const unsigned char *);

unsigned char *lookup_table();

typedef struct _fmi {
	char *bwt;
	int *idxs;
	int **rank_index;
	unsigned char* lookup;
	int endloc;
	int C[5];
	int len;
} fm_index;

// As the name suggests; deallocates all memory allocated for fmi, including
// fmi itself
void destroy_fmi(fm_index *fmi);

// Creates a FM-index from a given sequence using multithreaded histogram
// sort (allocating memory dynamically)
fm_index *make_fmi(const char *str, int len);

// Creates a FM-inde from a give sequence using SACA-K (allocating memory
// dynamically)
fm_index *make_fmi_sacak(const char *str, int len);

// Calculates the rank of a given symbol at a given index (i.e. the number
// of times the symbol has appeared up to that point) using the FM-index
// (Roughly constant time; this depends on implementation)
int rank(const fm_index *fmi, char c, int idx);

// Calculates the LF column mapping using the FM-index (constant time)
int lf(const fm_index *fmi, int idx);

// Searches for an exact pattern over the indexed sequence; returns
// the "first" (in this context, this means the rotation which appears
// first lexicographically) match, printing a warning message if there
// are multiple matches.
// pattern should be given uncompressed but in 0-3 form.
// Linear time in len * complexity of rank()
int reverse_search(const fm_index *fmi, const char *pattern, int len);

// Calculates SA[idx] from the FM-index
int unc_sa(const fm_index *fmi, int idx);

// Same as reverse_search
int locate(const fm_index *fmi, const char *pattern, int len);

// Same as locate, but returns the indices all matches (in the BWT; this means
// that you need to retrieve their indices using unc_sa) via sp and ep
void loc_search(const fm_index *fmi, const char *pattern, int len, int *sp, int *ep);

// Finds the maximum mappable suffix of a given pattern; returns the number
// of bases matched, storing matches in sp and ep as per loc_search
int mms(const fm_index *fmi, const char *pattern, int len, int *sp, int *ep);

// Prints part of a compressed sequence in human-readable form
void printseq(const char *seq, int startidx, int len);

#endif /* _SEQINDEX_H */
