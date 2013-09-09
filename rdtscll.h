#ifndef _RDTSCLL_H
#define _RDTSCLL_H

// Everyone's favorite implementation of rdtscll. This one only works on 32-bit
// platforms, you'll have to figure out the 64-bit one yourself. (Hint: trying
// this on a 64-bit platform will reveal that you can't quite use %edx:%eax,
// because %rax is the extension of %eax now, and the upper 32 bits still end
// up in %edx...
/*

  #define rdtscll(val)				\
  __asm__ __volatile__ ("rdtsc" : "=A" (val));

/*/

#define rdtscll(val) \
  __asm__ __volatile__ ("rdtsc\n\tshl $32, %%rdx\n\tor %%rdx, %%rax" :	\
			"=A" (val) : : "rdx");


#endif /* _RDTSCLL_H */
