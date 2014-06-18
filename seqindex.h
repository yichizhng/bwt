#ifndef _SEQINDEX_H
#define _SEQINDEX_H

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

void destroy_fmi(fm_index *fmi);

fm_index *make_fmi(const char *str, int len);

int rank(const fm_index *fmi, char c, int idx);

int lf(const fm_index *fmi, int idx);

int reverse_search(const fm_index *fmi, const char *pattern, int len);

int unc_sa(const fm_index *fmi, int idx);

int locate(const fm_index *fmi, const char *pattern, int len);

void loc_search(const fm_index *fmi, const char *pattern, int len, int *sp, int *ep);

int mms(const fm_index *fmi, const char *pattern, int len, int *sp, int *ep);

#endif /* _SEQINDXE_H */
