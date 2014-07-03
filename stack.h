#ifndef STACK_H_
#define STACK_H_

typedef struct stack_ {
  int size;
  int cap;
  int *counts;
  char *chars;
} stack;

stack *stack_make();

void stack_print_destroy(stack *s);

void stack_destroy(stack *s);

void stack_flip(stack *s1, stack *s2);

void stack_push(stack *s, char c, int count);

#endif /* STACK_H_ */
