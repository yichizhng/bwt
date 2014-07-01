#include <stdio.h>
#include <stdlib.h>
#include "stack.h"

stack *stack_make() {
  stack *s = malloc(sizeof(stack));
  if (!s)
    return 0;
  s->size = 0;
  s->cap = 10;
  s->counts = malloc(s->size * sizeof(int));
  if (!s->counts) {
    free(s);
    return 0;
  }
  s->chars = malloc(s->size);
  if (!s->chars) {
    free(s->counts);
    free(s);
    return 0;
  }
  return s;
}

// Destroys the stack and prints its contents (in CIGAR format, which is
// more or less a RLE)
void stack_destroy(stack *s) {
  while(s->size) {
    s->size--;
    printf("%d%c", s->counts[s->size], s->chars[s->size]);
  }
  free(s->counts);
  free(s->chars);
  free(s);
}

// Pushes some number of somethings onto the stack
void stack_push(stack *s, char c, int count) {
  if (s->chars[s->size-1] == c) {
    s->counts[s->size-1] += count;
    return;
  }
  else {
    if (s->size == s->cap) {
      s->cap += 10;
      int *newcounts = realloc(s->counts, s->cap * sizeof(int));
      char *newchars = realloc(s->chars, s->cap);
      if (newcounts)
	s->counts = newcounts;
      else
	return;
      if (newchars)
	s->chars = newchars;
      else
	return;
    }
    s->counts[s->size] = count;
    s->chars[s->size] = c;
    s->size++;
    return;
  }
}
