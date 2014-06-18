#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "rdtscll.h"
#include "histsortcomp.h"
#include "seqindex.h"
#include "csacak.h"

/* Some calculations
 * Space needed for uncompressed SA representation of human genome
 * 3.6 * 10^9 nts * 8.25 bytes/nt = 2.97*10^10 bytes = 28324.127MB = 27.66GB
 * This is an amount which is quite possible to store in main memory on a
 * modern server (or a very powerful consumer desktop)
 */

// Now with pthreads inside (and no concurrency primitives, because we have
// no data races!) for speed purposes

// Note that backwards search is O(n log m) if we assume m << 4^n; this is
// largely irrelevant (and a product of statistics)

static inline unsigned char getbase(const char *str, int idx) {
  // Gets the base at the appropriate index
  return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

// Note that this is single-threaded; we parallelize by calling it
// on different threads. Since it doesn't modify the index the only
// thing we have to worry about is cache locality (which is a real concern,
// but not really something we can fix)
void rna_seq(const fm_index *fmi, const char *pattern, int len) {
  // Copy pattern into a different aray for cyclic search
  char *pat = malloc(2*len);
  int i, j, start, end;
  memcpy(pat, pattern, len);
  memcpy(pat+len, pattern, len);
  for (i = len; i; --i) {
    loc_search(fmi, pattern+i, 14, &start, &end);
    if (start < end) {
    }
  }
  // Identify non-mapped segments and deal with them
  free(pat);
}

struct thread_args {
  fm_index *fmi;
  int num;
  char *pats;
  int start;
};

void *worker(void *arg) {
  int i, j, max;
  max = (((struct thread_args*)arg)->start) + 
    (((struct thread_args*)arg)->num);
  //printf("Starting thread to process %d patterns\n",
  //	max/4);
  for (i = ((int) ((struct thread_args*)arg)->start); i < max; i ++) {
    reverse_search(((struct thread_args*)arg)->fmi, 
		   ((struct thread_args*)arg)->pats + i, 12);
  }
  return NULL;
}

void startWorkers(fm_index *fmi, int num, char *pats) {
  // Starts 4 threads, each of which will process num/4 patterns
  // There are no synchronization primitives because they're unnecessary.
  pthread_t threads[4];
  int i;
  void *status;
  struct thread_args **ta = malloc(4 * sizeof (struct thread_args*));
  for (i = 0; i < 4; ++i) {
    ta[i] = malloc(sizeof(struct thread_args));
    ta[i]->fmi = fmi;
    ta[i]->num = num/4;
    ta[i]->pats = pats;
    ta[i]->start = i * (num/4);
  }
  for (i = 0; i < 4; ++i) {
    pthread_create(&threads[i], NULL, worker, (void *)ta[i]);
  }
  for (i = 0; i < 4; ++i) {
    pthread_join(threads[i], &status);
    // Ignoring return code, we didn't generate one anyway
    free(ta[i]);
  }
  free(ta);
  return;
}

// Misfeature: Index construction is O(n log n) on average; this is fast enough
// to dominate SACA-K below a billion base pairs or so, but uses too much
// memory

// For our use case (RNA sequencing) we can reverse search without backtracking
// on some fixed size seed (e.g. 14); STAR uses maximum mappable seed,
// which is sort of interesting as well
int main(int argc, char **argv) {
	
  // Performance testing code
  // To be honest, this should have the same performance as
  // histsortcomptest.c's testing code.
  
  // Single-threaded it takes about 4.6 seconds to construct the index
  // for 10^7 elements ; this suggests ~90 seconds for 10^8 elements
  // and 1800 seconds (30 minutes) for 10^9. Fortunately index-building
  // is a one-time thing, so runtime is not all that significant.
  
  // For some reason multithreaded index building is slower (overhead
  // and cache conflicts no doubt) below about 10 million base pairs
  // OTOH 4-threaded histogram sort is faster than SACA-K below
  // a billion nucleotides, so memory usage aside this should be
  // fine (and they make servers with enough memory to run it anyway)
  // Performance testing code
  int len, i, j, k;
  char *str, *pats;
  fm_index *fmi;
  unsigned long long a, b;
  if (argc != 2) {
    fprintf(stderr, "Usage: fmitest len\n");
    exit(-1);
  }

  len = atoi(argv[1]);
  str = malloc(len/4 + 1);
  for (i = 0; i < len/4 + 1; ++i)
    str[i] = rand();
  switch (len % 4) {
  case 1:
    str[len/4] &= 0xC0;
    break;
  case 2:
    str[len/4] &= 0xF0;
    break;
  case 3:
    str[len/4] &= 0xFC;
    break;
  case 0:
    str[len/4] = 0;
  }
  rdtscll(a);
  fmi = make_fmi(str, len);
  rdtscll(b);
  printf("Built index with %d base pairs in %lld cycles (%f s)\n",
	 len, b-a, ((double)(b-a)) / 2500000000.);
  printf("(%f cycles per base pair (%e seconds))\n", ((double)(b-a))	\
	 / len, ((double)(b-a)) / (len * 2500000000.));
  pats = malloc(sizeof(char) * 10000011);
  for (i = 0; i < 10000011; ++i) {
    pats[i] = rand() & 3;
  }
  // Comment: This is huge in 64-bit mode and thus makes valgrind
  // runs somewhat untenable (at least with that number of bps to
  // search; we can, of course, always scale things down but this
  // reveals more of the performance overhead from threading and less
  // asymptotic performance (which is what we actually care about)
  rdtscll(a);
  startWorkers(fmi, 10000000, pats);
  rdtscll(b);
  printf("Searched 10000000 12bp sequences in %lld cycles (%f s)\n",
	 b-a, ((double)(b-a)) / 2500000000.);
  printf("(%f cycles per base pair (%e seconds))\n", ((double)(b-a))	\
	 / 120000000., ((double)(b-a)) / 300000000000000000.);
  
  destroy_fmi(fmi);
  // With current settings the FMI with 1000000 base pairs takes
  // 1.5 million bytes to store (50% auxiliary data)
  // which is kind of bad but yeah
  free(str);
  free(pats);
}
