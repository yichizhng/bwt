#include <stdio.h>
#include <stdlib.h>
#include "stack.h"

stack *stack_make() {
  stack *s = malloc(sizeof(stack));
  if (!s)
    return 0;
  s->size = 0;
  s->cap = 10;
  s->counts = malloc(s->cap * sizeof(int));
  if (!s->counts) {
    free(s);
    return 0;
  }
  s->chars = malloc(s->cap);
  if (!s->chars) {
    free(s->counts);
    free(s);
    return 0;
  }
  return s;
}

// Destroys the stack and prints its contents (in CIGAR format, which is
// more or less a RLE)
void stack_print_destroy(stack *s) {
  putchar(' ');
  while(s->size) {
    s->size--;
    printf("%d%c", s->counts[s->size], s->chars[s->size]);
  }
  printf("\n");
  free(s->counts);
  free(s->chars);
  free(s);
}

// destroy the stack silently
void stack_destroy(stack *s) {
  while(s->size) {
    s->size--;
  }
  free(s->counts);
  free(s->chars);
  free(s);
}

// "Flips" stack 1 onto stack 2 (so that the bottom of stack 1 is now on
// the top of stack 2. Destroys stack 1.
void stack_flip(stack *s1, stack *s2) {
  while (s1->size) {
    s1->size--;
    stack_push(s2, s1->chars[s1->size], s1->counts[s1->size]);
  }
  free(s1->counts);
  free(s1->chars);
  free(s1);
}

// Pushes some number of somethings onto the stack
void stack_push(stack *s, char c, int count) {
  if (s->size  && (s->chars[s->size-1] == c)) {
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
