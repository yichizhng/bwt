#include "histsortcomp.h"
//#include "csacak.h"
#include "rdtscll.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static inline unsigned char getbase(const char *str, int idx) {
  // Gets the base at the appropriate index
  return ((str[idx>>2])>>(2*(3-(idx&3)))) & 3;
}

int main(int argc, char **argv) {
  /*
  // Testing code to verify correctness
  char *str, *bwt;
  int i;
  str = malloc(8 * sizeof(char));
  strcpy(str, ";I@7456"); // 0x3b494037343536
  // 03 23 10 21 10 00 03 13 03 10 03 11 03 12
  //for (i = 0; i < 7; ++i) {
  //	printf("%x", str[i]);
  //}
  //putchar('\n');
  bwt = makebwt(str, 28);
  puts(bwt);
  free(bwt);
  free(str); 
  */
  
  // Performance testing code
  unsigned int len;
  char *str, *bwt, *saca_bwt;
  unsigned int i;
  unsigned long long a, b;
  len = atoi(argv[1]);
  str = malloc(1+len/4 * sizeof(char));
  bwt = malloc(1+len/4 * sizeof(char));
  saca_bwt = malloc(1+len/4 * sizeof(char));
  for (i = 0; i < (len+1)/4; ++i) {
    str[i] = (char)rand();
  }
  str[len/4] &= (0xFF << 2*(4 - (len&3)));
  printf("Sentinel character is %d\n", getbase(str, len));
  // (Obfuscatedly) appends the sentinel character ;)
  //bwt = makebwt(str, len);
  rdtscll(a);
  i = makecbwt(str, len, bwt);
  rdtscll(b);
  puts("Using parallelized histogram sort");
  printf("%u nucleotides processed in %lld cycles (%f seconds)\n",
	 len, (b-a), ((double)(b-a)) / 2500000000.);
  rdtscll(a);
  i = saca_makecbwt(str, len, saca_bwt);
  rdtscll(b);
  puts("Using SACA-K");
  printf("%u nucleotides processed in %lld cycles (%f seconds)\n",
	 len, (b-a), ((double)(b-a)) / 2500000000.);
  // saca_makecbwt has the major differences of 1) not using nearly
  // as much memory and 2) not being as fast (due to saca-k being
  // damn near impossible to parallelize nicely)
  
  /*
    for (i = 0; i < len; ++i) {
    if (getbase(bwt, i) != getbase(saca_bwt, i)) {
    printf("Mismatch at base %u, %u %u", 
    i, getbase(bwt, i), getbase(saca_bwt,i));
    }
    }
    if (i == len)
    printf("Sequences match\n"); */
  free(str);
  free(bwt);
  free(saca_bwt);
  /* Saves ~7.5% memory, some performance gain at large array sizes
     Peak memory usage is about 9.25 bytes * len */
  return 0;
}
